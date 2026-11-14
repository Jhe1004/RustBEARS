"""Lagrange PHYLIP and BioGeoBEARS layer files."""
from dataclasses import dataclass, field
from pathlib import Path
from collections.abc import Mapping
import re

import numpy as np


@dataclass(frozen=True)
class TipRanges:
    areas: tuple[str, ...]
    ranges: dict[str, tuple[int, ...]]
    ambiguous: dict[str, tuple[int, ...]] = field(default_factory=dict)

    def likelihoods(self, labels, states, allow_null_tips=False, use_ambiguities=False):
        unknown = set(labels).symmetric_difference(self.ranges)
        if unknown:
            raise ValueError(f"tree and geography taxon sets differ: {sorted(unknown)}")
        rows = np.zeros((len(labels), len(states)), order="F")
        columns = {tuple(state): i for i, state in enumerate(states)}
        all_areas = set(range(len(self.areas)))
        for i, tip in enumerate(labels):
            ambig = self.ambiguous.get(tip, ())
            if ambig:
                if not use_ambiguities:
                    raise ValueError(
                        f"ambiguous '?' geography for {tip}; set use_ambiguities=True to expand tip ambiguities"
                    )
                if len(ambig) == len(self.areas):
                    rows[i, :] = 1.0
                else:
                    req = set(self.ranges[tip])
                    disallowed = all_areas - req - set(ambig)
                    for j, st in enumerate(states):
                        if st and req.issubset(st) and disallowed.isdisjoint(st):
                            rows[i, j] = 1.0
                if not np.any(rows[i, :] > 0):
                    raise ValueError(f"tip range for {tip} is outside the state space")
                continue
            state = self.ranges[tip]
            if not state and not allow_null_tips:
                raise ValueError(f"null tip range for {tip}")
            try:
                rows[i, columns[state]] = 1
            except KeyError as error:
                raise ValueError(f"tip range for {tip} is outside the state space") from error
        return rows


def _cell_to_bit(value, taxon):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"invalid boolean geography entry for {taxon}")
    if isinstance(value, str):
        token = value.strip()
        if token in ("0", "1", "?"):
            return token
        raise ValueError(f"invalid geography entry {value!r} for {taxon}")
    if isinstance(value, (int, float, np.number)):
        if not np.isfinite(value) or value not in (0, 1):
            raise ValueError(f"non-binary geography entry {value!r} for {taxon}")
        return "1" if value == 1 else "0"
    raise ValueError(f"invalid geography entry {value!r} for {taxon}")


def read_geog(source):
    if isinstance(source, TipRanges):
        return source
    if isinstance(source, Mapping):
        rows = list(source.items())
        if not rows:
            raise ValueError("empty tip geography")
        n = len(rows[0][1])
        areas = [chr(65 + i) if i < 26 else str(i + 1) for i in range(n)]
    elif hasattr(source, "columns") and hasattr(source, "iterrows"):
        areas = list(map(str, source.columns))
        rows = [(str(t), "".join(_cell_to_bit(x, t) for x in row)) for t, row in source.iterrows()]
        n = len(areas)
    else:
        lines = [l.strip() for l in Path(source).read_text().splitlines() if l.strip()]
        if not lines:
            raise ValueError("empty PHYLIP geography")
        m = re.fullmatch(r"(\d+)\s+(\d+)(?:\s*\(([^)]+)\))?", lines[0])
        if not m:
            raise ValueError("invalid Lagrange PHYLIP header")
        count, n = int(m[1]), int(m[2])
        areas = m[3].split() if m[3] else [chr(65+i) for i in range(n)]
        rows = [l.rsplit(None, 1) for l in lines[1:]]
        if len(rows) != count:
            raise ValueError("incorrect taxon count in PHYLIP geography")
    if n < 1 or len(areas) != n or len(set(areas)) != n:
        raise ValueError("invalid area names")
    found = {}
    ambiguous = {}
    for label, binary in rows:
        if label in found or not isinstance(binary, str) or len(binary) != n or set(binary) - {"0", "1", "?"}:
            raise ValueError(f"invalid or repeated geography record for {label}")
        found[label] = tuple(i for i, x in enumerate(binary) if x == "1")
        q_idx = tuple(i for i, x in enumerate(binary) if x == "?")
        if q_idx:
            ambiguous[label] = q_idx
    return TipRanges(tuple(areas), found, ambiguous)


def read_times(source):
    if source is None:
        return ()
    if isinstance(source, (str, Path)):
        source = [x for l in Path(source).read_text().splitlines() for x in l.split("#", 1)[0].split()]
    result = tuple(float(x) for x in source)
    if not np.all(np.isfinite(result)) or any(a >= b for a,b in zip((0, *result), result)):
        raise ValueError("timeperiods must be finite, positive and strictly increasing")
    return result


def read_layer_matrices(source, n, layers, default=1.0, *, missing_diagonal=False, areas=None, raw_layers=None):
    if source is None:
        return (np.full((n,n), default, order="F"),) * layers
    if isinstance(source, (str, Path)):
        tokens = []
        row = 0
        def numeric(value):
            return float("nan") if value.upper()=="NA" else float(value)

        for line in Path(source).read_text().splitlines():
            fields = line.split("#", 1)[0].split()
            if not fields or fields[0].upper() == "END":
                continue
            try:
                values = [numeric(v) for v in fields]
            except ValueError as error:
                if len(fields)==n and row%n==0 and len(set(fields))==n and (
                    areas is None or tuple(fields)==tuple(areas)
                ):
                    continue  # R's files may repeat the column header for each layer.
                if len(fields)!=n+1 or areas is not None and fields[0]!=areas[row%n]:
                    raise ValueError(f"invalid matrix row: {line!r}") from error
                values=[numeric(v) for v in fields[1:]]
            if len(values)!=n:
                raise ValueError(f"invalid matrix width in row {row+1}")
            tokens.extend(values)
            row+=1
        a = np.array(tokens).reshape(-1,n,n)
    else:
        a = np.asarray(source, dtype=float)
    if a.shape == (n,n):
        a = a[None,...]
    if raw_layers is not None and raw_layers > layers and a.shape == (raw_layers,n,n):
        a = a[:layers]
    if a.shape == (1,n,n):
        a = np.broadcast_to(a, (layers,n,n))
    if missing_diagonal and a.shape == (layers,n,n):
        a = np.array(a)
        diagonal = np.arange(n)
        diag = a[:,diagonal,diagonal]
        a[:,diagonal,diagonal] = np.where(np.isnan(diag) | (diag==0),1.0,diag)
    if a.shape != (layers,n,n) or not np.isfinite(a).all() or (a<0).any():
        raise ValueError(f"expected {layers} nonnegative finite {n}×{n} matrices")
    return tuple(np.asfortranarray(m) for m in a)
