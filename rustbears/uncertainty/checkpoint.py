"""Content-addressed manifests and atomic, lossless checkpoint records.

The JSONL index commits an already-written NPZ shard. An interrupted shard
write is harmless; an interrupted final index line is discarded on resume.
Only one writer may own a store (workers return records to their coordinator).
"""
from dataclasses import fields, is_dataclass, replace
import base64
import hashlib
import json
import os
from pathlib import Path
import threading
import uuid

import numpy as np


class CheckpointMismatchError(ValueError):
    """A checkpoint describes a different task or has damaged content."""


def canonical_manifest(value):
    """Encode types, IEEE floating bytes and complete array contents exactly."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, np.generic):
        return canonical_manifest(value.item())
    if isinstance(value, float):
        return {"__float64__": np.float64(value).tobytes().hex()}
    if isinstance(value, Path):
        return {"__path__": str(value)}
    if isinstance(value, np.ndarray):
        a = np.ascontiguousarray(value)
        if a.dtype.hasobject:
            raise TypeError("object arrays cannot enter checkpoint manifests")
        return {"__array__": a.dtype.str, "shape": list(a.shape),
                "bytes": base64.b64encode(a.tobytes()).decode("ascii")}
    if is_dataclass(value):
        return {"__dataclass__": type(value).__module__ + "." + type(value).__qualname__,
                "fields": {f.name: canonical_manifest(getattr(value, f.name)) for f in fields(value)}}
    if isinstance(value, dict) or hasattr(value, "items"):
        pairs = [[canonical_manifest(k), canonical_manifest(v)] for k, v in value.items()]
        return {"__mapping__": sorted(pairs, key=lambda kv: _json(kv[0]))}
    if isinstance(value, (tuple, list)):
        return {"__sequence__": type(value).__name__, "items": [canonical_manifest(x) for x in value]}
    raise TypeError(f"unsupported checkpoint manifest type: {type(value).__name__}")


def _json(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def task_fingerprint(manifest):
    return hashlib.sha256(_json(canonical_manifest(manifest)).encode()).hexdigest()


def model_manifest(fit_or_model):
    """Freeze resolved inputs, including matrices loaded from external files."""
    from .. import __version__
    from .._ffi import library, _resolve_native_library_path
    from ..model import FitResult
    fit = fit_or_model if isinstance(fit_or_model, FitResult) else None
    m = fit.model if fit is not None else fit_or_model
    result = dict(software_version=__version__, ffi_abi=int(library().bgb_m0_abi_version()),
                  native_sha256=hashlib.sha256(_resolve_native_library_path().read_bytes()).hexdigest(),
                  rng_layout="rustbears-task-identity-v1", tree=m.tree, geography=m.geog,
                  config=replace(m.config, tree=m.tree, geog=m.geog,
                                 distances=m.dist, envdistances=m.env,
                                 dispersal_multipliers=m.manual, area_of_areas=m.area_size),
                  model=m.model, states=m.states, areas=m.areas,
                  tips=m.tips, root_likelihood_weights=m.root_likelihood_weights,
                  root_sampling_probs=m.root_sampling_probs)
    for name in ("times", "masks", "dist", "env", "manual", "area_size", "fixes", "fixlikes"):
        if hasattr(m, name):
            result[name] = getattr(m, name)
    if fit is not None:
        result.update(fitted_params=fit.params, fitted_lnL=fit.lnL,
                      converged=fit.converged, precision=fit.precision_mode)
    return result


def _record_encode(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return {"__nonfinite_float__": value.hex()}
    if isinstance(value, dict):
        if not all(isinstance(k, str) for k in value):
            raise TypeError("record keys must be strings")
        return {k: _record_encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_record_encode(x) for x in value]
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError("records must contain JSON values; use arrays for numerical arrays")


def _record_decode(value):
    if isinstance(value, dict):
        if set(value) == {"__nonfinite_float__"}:
            return float.fromhex(value["__nonfinite_float__"])
        return {k: _record_decode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_record_decode(x) for x in value]
    return value


class CheckpointStore:
    def __init__(self, path, manifest, resume=False):
        self.path = Path(path)
        self.shards = self.path.with_name(self.path.name + ".shards")
        self.task_fingerprint = task_fingerprint(manifest)
        self._lock = threading.Lock()
        self._entries = {}
        header = {"schema": "rustbears.checkpoint.v1", "task_fingerprint": self.task_fingerprint,
                  "manifest": canonical_manifest(manifest)}
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.shards.mkdir(parents=True, exist_ok=True)
        if self.path.exists() and self.path.stat().st_size:
            if not resume:
                raise FileExistsError("checkpoint exists; explicitly use resume=True")
            raw = self.path.read_bytes()
            lines = raw.splitlines(keepends=True)
            try:
                observed = json.loads(lines[0])
            except (ValueError, IndexError) as error:
                raise CheckpointMismatchError("invalid checkpoint header") from error
            if observed != header:
                raise CheckpointMismatchError("checkpoint task fingerprint/manifest mismatch")
            committed = len(lines[0])
            for index, line in enumerate(lines[1:], 1):
                if not line.endswith(b"\n") and index == len(lines) - 1:
                    break
                try:
                    entry = json.loads(line)
                    rid = self._id(entry["record"]["id"])
                    if rid in self._entries:
                        raise ValueError("duplicate record")
                    self._verify_shard(entry)
                except (ValueError, KeyError, OSError) as error:
                    raise CheckpointMismatchError("corrupt committed checkpoint record") from error
                self._entries[rid] = entry
                committed += len(line)
            if committed != len(raw):
                with self.path.open("r+b") as handle:
                    handle.truncate(committed)
        else:
            tmp = self.path.with_name(self.path.name + "." + uuid.uuid4().hex + ".tmp")
            with tmp.open("w") as handle:
                handle.write(_json(header) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(tmp, self.path)

    @staticmethod
    def _id(value):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError("checkpoint record id must be a positive integer")
        return int(value)

    def _verify_shard(self, entry):
        if entry.get("shard") is None:
            return
        name = entry["shard"]
        if Path(name).name != name:
            raise ValueError("invalid shard path")
        if hashlib.sha256((self.shards / name).read_bytes()).hexdigest() != entry["shard_sha256"]:
            raise CheckpointMismatchError("checkpoint array shard checksum mismatch")

    @property
    def completed_ids(self):
        return tuple(sorted(self._entries))

    @property
    def records(self):
        return tuple(_record_decode(self._entries[k]["record"]) for k in self.completed_ids)

    def load_arrays(self, record_id):
        entry = self._entries[self._id(record_id)]
        self._verify_shard(entry)
        if entry["shard"] is None:
            return {}
        with np.load(self.shards / entry["shard"], allow_pickle=False) as z:
            return {key: z[key].copy() for key in z.files}

    def append(self, record, arrays=None):
        record = _record_encode(record)
        rid = self._id(record["id"])
        with self._lock:
            if rid in self._entries:
                raise ValueError(f"checkpoint record {rid} already committed")
            entry = {"record": record, "shard": None, "shard_sha256": None}
            if arrays:
                values = {str(k): np.asarray(v) for k, v in arrays.items()}
                if any(a.dtype.hasobject for a in values.values()):
                    raise TypeError("checkpoint arrays cannot have object dtype")
                name = f"{rid:012d}-{uuid.uuid4().hex}.npz"
                tmp = self.shards / (name + ".tmp")
                with tmp.open("wb") as handle:
                    np.savez_compressed(handle, **values)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(tmp, self.shards / name)
                entry.update(shard=name, shard_sha256=hashlib.sha256((self.shards / name).read_bytes()).hexdigest())
            with self.path.open("a") as handle:
                handle.write(_json(entry) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            self._entries[rid] = entry


def frame_arrays(frame, prefix):
    """Losslessly persist event-table columns without pickle/object arrays."""
    import pandas as pd
    result = {}
    dtypes = {}
    for index, column in enumerate(frame.columns):
        series = frame[column]
        key = f"{prefix}_{index}"
        if isinstance(series.dtype, pd.StringDtype):
            # pandas 3 infers a StringDtype for ordinary text columns. Keep
            # its storage and missing-value convention, rather than silently
            # restoring it as an object column (or changing pandas.NA to NaN).
            result[key] = np.asarray([
                _json(None if pd.isna(v) else str(v)) for v in series
            ], dtype=str)
            dtypes[column] = {
                "kind": "pandas_string", "storage": series.dtype.storage,
                "na_value": "pd.NA" if series.dtype.na_value is pd.NA else "nan",
            }
        elif series.dtype.kind == "O":
            encoded = [_json(_record_encode(v)) for v in series]
            result[key] = np.asarray(encoded, dtype=str)
            dtypes[column] = "json"
        else:
            result[key] = series.to_numpy(copy=True)
            dtypes[column] = str(series.dtype)
    return {"columns": list(frame.columns), "dtypes": dtypes, "prefix": prefix}, result


def restore_frame(metadata, arrays):
    import pandas as pd
    values = {}
    for index, column in enumerate(metadata["columns"]):
        a = arrays[f"{metadata['prefix']}_{index}"]
        dtype = metadata["dtypes"][column]
        if isinstance(dtype, dict) and dtype.get("kind") == "pandas_string":
            kwargs = {"storage": dtype["storage"]}
            if dtype["na_value"] == "nan":
                kwargs["na_value"] = np.nan
            values[column] = pd.Series([json.loads(v) for v in a],
                                       dtype=pd.StringDtype(**kwargs))
        elif dtype == "json":
            values[column] = pd.Series([_record_decode(json.loads(v)) for v in a],
                                       dtype=object)
        else:
            values[column] = pd.Series(a, dtype=dtype)
    return pd.DataFrame(values)
