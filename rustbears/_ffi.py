"""Ownership and NumPy marshalling for the private Rust C ABI."""
import ctypes as ct
from functools import lru_cache
from pathlib import Path
import sys
import threading

import numpy as np

_P = ct.POINTER(ct.c_double)
_ARGS = [ct.POINTER(_P), ct.POINTER(ct.c_size_t), ct.c_size_t]
_PRECISION = {"legacy_f32": 0, "strict_f64": 1, "legacy_fused": 2}
_EMPTY = np.empty(0, dtype=np.float64)


def _resolve_native_library_path() -> Path:
    if sys.platform == "win32":
        filename = "bgb_ffi.dll"
    elif sys.platform == "darwin":
        filename = "libbgb_ffi.dylib"
    else:
        filename = "libbgb_ffi.so"
    packaged = Path(__file__).with_name(filename)
    if packaged.exists():
        return packaged
    dev_candidate = Path(__file__).resolve().parents[1] / "rust" / "target" / "release" / filename
    if dev_candidate.exists():
        return dev_candidate
    return packaged


@lru_cache(maxsize=1)
def library():
    dll = ct.CDLL(str(_resolve_native_library_path()))
    ptr = ct.c_void_p
    size = ct.c_size_t
    for name, args, returns in (
        ("bgb_m0_abi_version", [], ct.c_uint),
        ("bgb_kernel", [ct.c_uint, *_ARGS, ct.c_uint], ptr),
        ("bgb_result_error", [ptr], ct.c_char_p),
        ("bgb_result_count", [ptr], size),
        ("bgb_result_len", [ptr, size], size),
        ("bgb_result_data", [ptr, size], _P),
        ("bgb_result_free", [ptr], None),
        ("bgb_context_create", [*_ARGS, size], ptr),
        ("bgb_strat_context_create", [*_ARGS, size], ptr),
        ("bgb_context_error", [ptr], ct.c_char_p),
        ("bgb_context_free", [ptr], None),
        ("bgb_context_evaluate", [ptr, *_ARGS, ct.c_uint], ptr),
        ("bgb_context_stats", [ct.POINTER(size)], None),
        ("bgb_register_dgemm", [ptr], None),
        ("bgb_register_dgesv", [ptr], None),
    ):
        f = getattr(dll, name)
        f.argtypes, f.restype = args, returns
    if dll.bgb_m0_abi_version() != 1:
        raise ImportError("unsupported RustBEARS native ABI")
    _register_blas(dll)
    return dll


def _register_blas(dll):
    # SciPy owns the native capsule for the duration of the Python process.
    # Validate its pointer signature before crossing the ABI boundary.
    from scipy.linalg import cython_blas, cython_lapack

    signatures = (
        ("dgemm", cython_blas,
         b"void (char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *)",
         b"__pyx_t_5scipy_6linalg_11cython_blas_d"),
        ("dgesv", cython_lapack, b"void (int *, int *, double *, int *, int *, double *, int *, int *)",
         b"__pyx_t_5scipy_6linalg_13cython_lapack_d"),
    )
    name_of = ct.pythonapi.PyCapsule_GetName
    name_of.argtypes = [ct.py_object]
    name_of.restype = ct.c_char_p
    address_of = ct.pythonapi.PyCapsule_GetPointer
    address_of.argtypes = [ct.py_object, ct.c_char_p]
    address_of.restype = ct.c_void_p
    for symbol,module,signature,cython_double in signatures:
        capsule=module.__pyx_capi__.get(symbol)
        if capsule is None:
            continue
        name=name_of(capsule)
        if name and name.replace(cython_double,b"double")==signature:
            getattr(dll,f"bgb_register_{symbol}")(address_of(capsule,name))


def _buffers(values):
    arrays = [np.require(x, dtype=np.float64, requirements=("F", "A")).ravel(order="F")
              for x in values]
    return (arrays, (_P * len(arrays))(*(a.ctypes.data_as(_P) for a in arrays)),
            (ct.c_size_t * len(arrays))(*(a.size for a in arrays)))


def _reply(dll, reply):
    if not reply:
        raise MemoryError("Rust returned no result")
    try:
        error = dll.bgb_result_error(reply)
        if error:
            raise ValueError(error.decode("utf-8", errors="replace"))
        out = []
        for i in range(dll.bgb_result_count(reply)):
            n = dll.bgb_result_len(reply, i)
            out.append(np.ctypeslib.as_array(dll.bgb_result_data(reply, i), (n,)).copy()
                       if n else _EMPTY.copy())
        return out
    finally:
        dll.bgb_result_free(reply)


def kernel(op, *buffers, precision="strict_f64"):
    dll = library()
    arrays, ptrs, lengths = _buffers(buffers)
    return _reply(dll, dll.bgb_kernel(op, ptrs, lengths, len(arrays), _PRECISION[precision]))


class CompiledContext:
    """A uniquely owned Rust context. Evaluation and close are serialized."""
    def __init__(self, config, *, stratified=False, num_threads=1):
        self._ptr = None
        self._lock = threading.RLock()
        self._lib = library()
        if (
            isinstance(num_threads, (bool, np.bool_))
            or not isinstance(num_threads, (int, np.integer))
            or not (1 <= int(num_threads) <= 65536)
        ):
            raise ValueError("num_threads must be an integer between 1 and 65536")
        num_threads = int(num_threads)
        arrays, ptrs, lengths = _buffers(config)
        create = self._lib.bgb_strat_context_create if stratified else self._lib.bgb_context_create
        self._ptr = create(ptrs, lengths, len(arrays), num_threads)
        if not self._ptr:
            raise MemoryError("Rust returned no context")
        error = self._lib.bgb_context_error(self._ptr)
        if error:
            self.close()
            raise ValueError(error.decode("utf-8", errors="replace"))

    def evaluate(self, buffers, precision):
        arrays, ptrs, lengths = _buffers(buffers)
        with self._lock:
            if self._ptr is None:
                raise RuntimeError("compiled context is closed")
            reply = self._lib.bgb_context_evaluate(
                self._ptr, ptrs, lengths, len(arrays), _PRECISION[precision])
            return _reply(self._lib, reply)

    def close(self):
        with self._lock:
            if self._ptr is not None:
                self._lib.bgb_context_free(self._ptr)
                self._ptr = None

    def __enter__(self):
        if self._ptr is None:
            raise RuntimeError("compiled context is closed")
        return self

    def __exit__(self, *_):
        self.close()

    def __copy__(self):
        raise TypeError("CompiledContext owns a native Rust handle and cannot be copied")

    def __deepcopy__(self, memo):
        raise TypeError("CompiledContext owns a native Rust handle and cannot be deep-copied")

    def __reduce__(self):
        raise TypeError("CompiledContext owns a native Rust handle and cannot be pickled")

    def __del__(self):
        if getattr(self, "_ptr", None) is not None:
            self.close()


def context_stats():
    values = (ct.c_size_t * 4)()
    library().bgb_context_stats(values)
    return dict(created=values[0] + values[2], dropped=values[1] + values[3],
                single_created=values[0], stratified_created=values[2])
