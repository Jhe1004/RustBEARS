"""One CPU budget for outer tasks and their native likelihood evaluations."""
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass

import numpy as np
from threadpoolctl import threadpool_limits


def positive_int(value, name, maximum=None):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    if value < 1 or maximum is not None and value > maximum:
        raise ValueError(f"{name} must be positive" + (f" and <= {maximum}" if maximum else ""))
    return int(value)


@dataclass(frozen=True)
class ParallelConfig:
    """Outer workers own independent models; native workers never multiply them.

    BLAS is limited once for the whole executor lifetime. A ThreadPoolExecutor
    is used because the native kernel releases the GIL; task inputs must not
    contain compiled contexts or sampler handles.
    """
    workers: int = 1
    native_threads: int = 1

    def __post_init__(self):
        positive_int(self.workers, "workers", 256)
        positive_int(self.native_threads, "native_threads", 256)
        if self.workers > 1 and self.native_threads != 1:
            raise ValueError("outer parallelism requires native_threads=1")

    @property
    def inner_threads(self):
        return self.native_threads


def resolve_parallel(parallel, num_threads):
    if parallel is not None:
        if not isinstance(parallel, ParallelConfig):
            raise TypeError("parallel must be a ParallelConfig")
        return parallel
    return ParallelConfig(workers=positive_int(num_threads, "num_threads", 256))


@contextmanager
def parallel_scope(config):
    """Yield an ordered map function; no per-worker thread-limit restoration."""
    with threadpool_limits(limits=1):
        if config.workers == 1:
            yield map
        else:
            with ThreadPoolExecutor(max_workers=config.workers) as executor:
                yield executor.map
