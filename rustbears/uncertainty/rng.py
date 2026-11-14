"""Stable random-task identities for new uncertainty workflows.

Legacy single-tree BSM passes its seed directly to Rust and does not use this
module. New workflows derive a seed from the complete logical task identity,
independently of scheduling, worker count, batch size, or retry order.

SHA-256 domain separation does not make a 64-bit seed collision-free. Declare
the complete batch through ``derive_task_seeds`` before dispatch: it rejects
duplicate task identities and observed truncated-seed collisions. Checkpoints
should retain the identity and full ``task_fingerprint`` as well as the seed.
"""
from dataclasses import asdict, dataclass
import hashlib
import json
from numbers import Integral

_DOMAIN = b"RustBEARS/uncertainty/task-seed/v1\x00"
_MAX_U64 = (1 << 64) - 1


def _uint64(value, name):
    if isinstance(value, bool) or not isinstance(value, Integral) or not 0 <= value <= _MAX_U64:
        raise ValueError(f"{name} must be an integer between 0 and 2**64-1")
    return int(value)


@dataclass(frozen=True)
class TaskIdentity:
    """A logical task key; IDs must be stable labels, never worker positions.

    ``purpose`` distinguishes simulation, optimization, mapping, or diagnostics;
    ``role`` distinguishes null/alternative fits or other tasks within a phase.
    Zero is a valid index. Callers should use a single documented index scheme.
    """
    purpose: str
    model_id: str
    tree_id: str
    outer_rep: int = 0
    inner_rep: int = 0
    role: str = ""
    layer_id: int = 0
    step_id: int = 0

    def __post_init__(self):
        for name in ("purpose", "model_id", "tree_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a nonempty stable string identifier")
        if not isinstance(self.role, str):
            raise ValueError("role must be a string")
        for name in ("outer_rep", "inner_rep", "layer_id", "step_id"):
            object.__setattr__(self, name, _uint64(getattr(self, name), name))

    def to_dict(self):
        """Return the exact serializable fields used to derive this task seed."""
        return asdict(self)


class SeedCollisionError(ValueError):
    """Two different task identities produced the same 64-bit seed."""


def canonical_task_payload(master_seed, task):
    """Canonical JSON text shared with the R wrapper (domain prefix excluded)."""
    master_seed = _uint64(master_seed, "master_seed")
    if not isinstance(task, TaskIdentity):
        raise TypeError("task must be a TaskIdentity")
    payload = {"master_seed": master_seed, "task": task.to_dict()}
    return json.dumps(payload, sort_keys=True, ensure_ascii=True,
                      allow_nan=False, separators=(",", ":"))


def _task_digest(master_seed, task):
    canonical = canonical_task_payload(master_seed, task).encode("ascii")
    return hashlib.sha256(_DOMAIN + canonical).digest()


def task_fingerprint(master_seed, task):
    """Full SHA-256 hex digest for audit/checkpoint identity, not a native seed."""
    return _task_digest(master_seed, task).hex()


def derive_seed(master_seed, task):
    """Derive one u64 seed; batch callers must use ``derive_task_seeds``."""
    return int.from_bytes(_task_digest(master_seed, task)[:8], "big")


def derive_task_seeds(master_seed, tasks):
    """Return ``{TaskIdentity: seed}``, rejecting repeated keys or seed collisions.

    Keys include every semantic dimension relevant to the workflow. A retry
    reuses its previously assigned seed; it is not declared as a second task.
    Collision errors are explicit and are never repaired by silently changing
    seeds, since that would make results depend on dispatch order.
    """
    _uint64(master_seed, "master_seed")
    result, owners = {}, {}
    for task in tasks:
        if not isinstance(task, TaskIdentity):
            raise TypeError("every task must be a TaskIdentity")
        if task in result:
            raise ValueError(f"duplicate random task identity: {task}")
        seed = derive_seed(master_seed, task)
        if seed in owners:
            raise SeedCollisionError(f"64-bit task seed collision between {owners[seed]} and {task}")
        result[task] = seed
        owners[seed] = task
    return result
