"""Rooted binary Newick with ape-compatible node IDs and pruningwise edges."""
from dataclasses import dataclass
import math
from pathlib import Path
import re

import numpy as np


@dataclass
class _Node:
    label: str = ""
    length: float = 0.0
    children: tuple = ()
    index: int = -1


@dataclass
class PhyloTree:
    labels: tuple[str, ...]
    children: tuple[tuple[int, ...], ...]
    lengths: np.ndarray
    ages: np.ndarray
    tip_count: int
    root_edge: float = 0.0
    fossils_older_than: float = 0.001

    @property
    def root(self):
        return self.tip_count

    @property
    def tip_labels(self):
        return self.labels[:self.tip_count]

    @property
    def parent(self):
        out = np.full(len(self.labels), -1, dtype=int)
        for node, children in enumerate(self.children):
            out[list(children)] = node
        return out

    @property
    def edges(self):
        return pruningwise_edges(self.children, self.root)

    def node_table(self):
        import pandas as pd

        cladewise = []
        edge_stack = list(reversed(self.children[self.root]))
        while edge_stack:
            child = edge_stack.pop()
            cladewise.append(child)
            for grand in reversed(self.children[child]):
                edge_stack.append(grand)
        edge_numbers = {c: i + 1 for i, c in enumerate(cladewise)}
        descendants = {}
        post_stack = [(self.root, False)]
        while post_stack:
            node, expanded = post_stack.pop()
            ch = self.children[node]
            if not ch:
                descendants[node] = (node + 1,)
            elif expanded:
                descendants[node] = tuple(t for c in ch for t in descendants[c])
            else:
                post_stack.append((node, True))
                for c in reversed(ch):
                    post_stack.append((c, False))
        rows = []
        parent = self.parent
        for i, label in enumerate(self.labels):
            children = self.children[i]
            rows.append({
                "node": i + 1,
                "node.type": "root" if i == self.root else "internal" if children else "tip",
                "parent_br": edge_numbers.get(i),
                "edge.length": self.root_edge if i == self.root else self.lengths[i],
                "ancestor": int(parent[i] + 1) if parent[i] >= 0 else None,
                "daughter_nds": tuple(c + 1 for c in children),
                "time_bp": round(float(self.ages[i]), 7),
                "label": label or (f"inNode{i+1}" if children else label),
                "fossils": bool(self.ages[i] > self.fossils_older_than) if not children else None,
                "left_desc_nodes": descendants[children[0]] if children else (),
                "right_desc_nodes": descendants[children[1]] if children else (),
            })
        return pd.DataFrame(rows)


_TOKEN = re.compile(r"\s*(\[[^\]]*\]|'(?:[^']|'')*'|[^\s()\[\]:;,']+|[():;,])")


def read_newick(source, *, fossils_older_than=0.001):
    if isinstance(source, PhyloTree):
        return source
    if isinstance(source, Path):
        text = source.read_text()
    else:
        raw_str = str(source)
        if not any(c in raw_str for c in "();"):
            text = Path(raw_str).read_text()
        else:
            is_existing_file = False
            if "\n" not in raw_str and "\r" not in raw_str and len(raw_str) <= 4096:
                try:
                    is_existing_file = Path(raw_str).is_file()
                except OSError:
                    is_existing_file = False
            text = Path(raw_str).read_text() if is_existing_file else raw_str
    tokens = [m.group(1) for m in _TOKEN.finditer(text) if not m.group(1).startswith("[")]
    if _TOKEN.sub("", text).strip():
        raise ValueError("unrecognized Newick input")
    i = 0

    def take():
        nonlocal i
        if i == len(tokens):
            raise ValueError("unexpected end of Newick tree")
        value = tokens[i]
        i += 1
        return value

    def unquote(raw):
        return raw[1:-1].replace("''", "'") if raw.startswith("'") else raw

    def read_length():
        nonlocal i
        if i < len(tokens) and tokens[i] == ":":
            take()
            length = float(take())
            if not np.isfinite(length) or length < 0:
                raise ValueError("Newick branch lengths must be finite and nonnegative")
            return length
        return 0.0

    frames = []
    current = None
    while True:
        if current is None:
            value = take()
            if value == "(":
                frames.append([])
                continue
            if value in (")", ",", ":", ";"):
                raise ValueError("missing Newick tip label")
            current = _Node(unquote(value), read_length(), ())
        if not frames:
            root = current
            break
        frames[-1].append(current)
        current = None
        sep = take()
        if sep == ",":
            continue
        if sep != ")" or len(frames[-1]) != 2:
            raise ValueError("the likelihood engine requires a rooted binary tree")
        children = tuple(frames.pop())
        label = ""
        if i < len(tokens) and tokens[i] not in (":", ",", ")", ";"):
            label = unquote(take())
        current = _Node(label, read_length(), children)

    if take() != ";" or i != len(tokens) or not root.children:
        raise ValueError("expected one rooted Newick tree ending in ';'")
    tips, internal, preorder = [], [], []
    order_stack = [root]
    while order_stack:
        node = order_stack.pop()
        preorder.append(node)
        (internal if node.children else tips).append(node)
        for child in reversed(node.children):
            order_stack.append(child)
    nodes = tips + internal
    if len({x.label for x in tips}) != len(tips):
        raise ValueError("tree tip labels must be unique")
    for idx, node in enumerate(nodes):
        node.index = idx
    heights_list = [0.0] * len(nodes)
    for node in preorder:
        base_h = heights_list[node.index]
        for c in node.children:
            child_h = base_h + float(c.length)
            if not math.isfinite(child_h):
                raise ValueError("Newick root-to-node path lengths overflowed finite float64 range")
            heights_list[c.index] = child_h
    heights = np.asarray(heights_list, dtype=float)
    ages = float(heights[:len(tips)].max()) - heights
    if not np.isfinite(ages).all():
        raise ValueError("Newick node ages must be finite")
    ages[np.abs(ages) < 1e-12] = 0
    return PhyloTree(tuple(n.label for n in nodes), tuple(tuple(c.index for c in n.children) for n in nodes),
                     np.array([n.length if n is not root else 0.0 for n in nodes]), ages, len(tips),
                     root.length, fossils_older_than)

def pruningwise_edges(children, root):
    """Collect ready sibling edge groups in cladewise order, in waves."""
    stack = [(root, False)]
    depths = [0] * len(children)
    ready = []
    while stack:
        node, expanded = stack.pop()
        ch = children[node]
        if not ch:
            continue
        if expanded:
            d = 1 + max(depths[c] for c in ch)
            depths[node] = d
            ready.append((d, node))
        else:
            stack.append((node, True))
            for c in reversed(ch):
                stack.append((c, False))
    return np.array([(p,c) for _,p in sorted(ready) for c in children[p]], dtype=int)
