"""Shared prospective-budget validation; this module never executes an experiment.

The statistical design may be frozen separately before Gate 2/3. The final
budget can then reference that design and the completed gate without a circular
hash dependency. Experiment-specific settings remain the caller's responsibility.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def load_locked_budget(path, *, source_script, execute=False):
    path = Path(path).resolve()
    raw = path.read_bytes()
    budget = json.loads(raw)
    if budget.get("schema") != "rustbears.formal_budget.v1" or budget.get("status") != "locked":
        raise ValueError("budget must use schema rustbears.formal_budget.v1 and status locked")
    locked = datetime.fromisoformat(budget["locked_utc"].replace("Z", "+00:00"))
    if locked.tzinfo is None or locked > datetime.now(timezone.utc):
        raise ValueError("budget lock must have a past timezone-aware timestamp")
    if not isinstance(budget.get("jobs"), list) or not budget["jobs"]:
        raise ValueError("locked budget must contain a nonempty jobs list")
    script = Path(source_script).resolve()
    relative = str(script.relative_to(REPO))
    for name, source in ((relative, script), ("scripts/v12_formal_protocol.py", Path(__file__))):
        if budget.get("source_sha256", {}).get(name) != hashlib.sha256(source.read_bytes()).hexdigest():
            raise ValueError("locked budget does not match source script " + name)
    # The optional upstream design is immutable and must predate the final lock.
    design = budget.get("design_evidence")
    if design is not None:
        design_path = Path(design["path"])
        if not design_path.is_absolute():
            design_path = path.parent / design_path
        if hashlib.sha256(design_path.read_bytes()).hexdigest() != design.get("sha256"):
            raise ValueError("locked budget does not match its frozen statistical design")
        design_payload = json.loads(design_path.read_text())
        if design_payload.get("jobs") != budget["jobs"]:
            raise ValueError("final budget settings differ from the frozen statistical design")
    evidence = budget.get("gate_evidence", {})
    gate_path = Path(evidence.get("path", ""))
    if not gate_path.is_absolute():
        gate_path = path.parent / gate_path
    if execute:
        if not gate_path.is_file() or hashlib.sha256(gate_path.read_bytes()).hexdigest() != evidence.get("sha256"):
            raise ValueError("formal execution requires matching hashed Gate 2/3 evidence")
        gate = json.loads(gate_path.read_text())
        if gate.get("gate2_passed") is not True or gate.get("gate3_passed") is not True:
            raise ValueError("Gate 2 and Gate 3 must both pass before formal execution")
    return budget, dict(path=str(path), sha256=hashlib.sha256(raw).hexdigest(),
                       locked_utc=locked.isoformat(), payload=budget,
                       gate_verification_required_before_execution=True)


def select_budget_job(budget, **criteria):
    """Require one fully specified job; no fallback to a smaller or smoke budget."""
    matches = [job for job in budget["jobs"] if all(job.get(key) == value for key, value in criteria.items())]
    if len(matches) != 1:
        raise ValueError("budget requires exactly one matching job for " + repr(criteria))
    return matches[0]
