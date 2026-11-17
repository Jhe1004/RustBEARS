"""Slice a rooted binary tree and compile layer-local pruning buffers."""
from dataclasses import dataclass

import numpy as np

from ._ffi import _EMPTY
from .io import read_newick, read_times
from .io.newick import pruningwise_edges
from .operators import _states


@dataclass
class Piece:
    nodes: tuple[int, ...]
    sources: tuple[int, ...]
    edges: np.ndarray
    lengths: np.ndarray
    duration: float
    rows: tuple[int, ...] = ()

    @property
    def branch(self):
        return len(self.nodes) == 1

    @property
    def master_root(self):
        return self.nodes[0] if self.branch else self.nodes[len(self.sources)]


@dataclass
class TreeSections:
    tree: object
    timeperiods: tuple
    pieces: tuple
    rows: list

    @property
    def master_table(self):
        import pandas as pd

        return pd.DataFrame(self.rows)


def effective_times(times, root_age):
    """Keep every occupied layer and make its final bound cover the root."""
    if not times:
        return ()
    oldest = int(np.searchsorted(times, root_age, side="left"))
    keep = min(oldest + 1, len(times))
    if times[keep - 1] <= root_age:
        return (*times[:keep - 1], root_age + 1.0)
    return times[:keep]


def section_the_tree(tree, timeperiods, *, min_dist_between_node_and_stratum_line=1e-5, exact_ages=False):
    tree = read_newick(tree)
    times = effective_times(read_times(timeperiods),float(tree.ages[tree.root]))
    ages = tree.ages
    root = tree.root
    nt = tree.tip_count
    parents = tree.parent
    if not times:
        return TreeSections(tree,(),(),[])
    if any(np.min(np.abs(ages - age)) < min_dist_between_node_and_stratum_line for age in times[:-1]):
        raise ValueError("a node is too close to a stratum boundary")

    def add_row(layer, piece, node, kind, idx, length, node_type, lower, upper):
        node_age = float(ages[node]) if exact_ages else round(float(ages[node]),7)
        rows.append({"node": node+1, "node.type": "root" if node==root else "tip" if node<nt else "internal",
                     "edge.length": tree.lengths[node], "time_bp": node_age, "label": tree.labels[node],
                     "fossils": node<nt and ages[node]>tree.fossils_older_than,
                     "stratum":layer,"piecenum":piece,"piececlass":kind,"time_top":lower,"time_bot":upper,
                     "SUBnode":idx+1,"SUBedge.length":length,"SUBnode.type":node_type})

    rows = []
    for node in range(nt):
        add_row(0,0,node,"orig_tip",node,tree.lengths[node],"orig_tip",0,0)
    current = tuple(range(nt))
    active_lengths = tree.lengths.copy()
    if np.any(ages[:nt] > tree.fossils_older_than):
        active_lengths[:nt] += ages[:nt]
    original_ages = np.asarray(ages, dtype=float) if exact_ages else np.round(ages,7)
    layers = []
    relative = np.asarray(times).copy()

    for i, bot in enumerate(times, 1):
        slab = float(relative[i-1])
        relative -= slab
        top = bot-slab
        slots = {node:j for j,node in enumerate(current)}
        depths = {}
        measure_stack = [(root, 0.0)]
        while measure_stack:
            node, depth = measure_stack.pop()
            depths[node] = depth
            if node not in slots:
                for child in reversed(tree.children[node]):
                    measure_stack.append((child, depth + active_lengths[child]))
        height=max(depths[t] for t in current)

        if i == len(times):
            tops = [root]
        else:
            tops = []
            frontier_stack = [root]
            while frontier_stack:
                node = frontier_stack.pop()
                if node in slots or height - depths[node] <= slab:
                    tops.append(node)
                else:
                    for child in reversed(tree.children[node]):
                        frontier_stack.append(child)
        layer = []
        retained = []
        for j, base in enumerate(tops, 1):
            if base in slots:
                seg_len = (min(bot, float(ages[parents[base]])) - max(top, float(ages[base]))) if exact_ages else slab
                add_row(i,j,base,"subbranch",slots[base],active_lengths[base],"tip",top,bot)
                layer.append(Piece((base,),(slots[base],),np.empty((0,2),int),_EMPTY,seg_len,(len(rows)-1,)))
                retained.append((base,active_lengths[base],depths[base]))
                continue

            tips, inner, heights = [], [], {}
            visit_stack = [(base, 0.0)]
            while visit_stack:
                node, depth = visit_stack.pop()
                if node in slots:
                    tips.append(node)
                    heights[node] = depth
                else:
                    inner.append(node)
                    for child in reversed(tree.children[node]):
                        visit_stack.append((child, depth + active_lengths[child]))
            tips.sort(key=slots.__getitem__)
            nodes = tuple(tips+inner)
            idx = {node:k for k,node in enumerate(nodes)}
            children = tuple(tuple(idx[c] for c in tree.children[n]) if n in inner else () for n in nodes)
            edge = pruningwise_edges(children,len(tips))
            if exact_ages:
                lengths = np.asarray([float(ages[parents[nodes[c]]]) - max(top, float(ages[nodes[c]])) for _,c in edge], dtype=float)
                duration = 0.0 if base==root else float(bot - ages[base])
            else:
                lengths = np.asarray([active_lengths[nodes[c]] for _,c in edge])
                for k,(_,child) in enumerate(edge):
                    master=nodes[child]
                    if master in slots and original_ages[master]-top>1e-10:
                        lengths[k]-=original_ages[master]-top
                duration = 0.0 if base==root else slab-max(heights.values())
            local_lengths = {int(c):length for (_,c),length in zip(edge,lengths)}
            first = len(rows)
            for k,node in enumerate(nodes):
                role = "root" if node==base else "tip" if node in slots else "internal"
                add_row(i,j,node,"subtree",k,local_lengths.get(k,duration),role,top,bot)
            layer.append(Piece(nodes,tuple(slots[t] for t in tips),edge,lengths,duration,
                               tuple(range(first,len(rows)))))
            last=tips[-1]
            retained.append((base,active_lengths[base]+heights[last],depths[last]))
        current = tuple(piece.master_root for piece in layer)
        layers.append(tuple(layer))
        if i<len(times):
            # ape::drop.tip retains the last tip of each cut clade. The
            # historical chainsaw shortens those collapsed edges using prt's
            # ages rounded to seven decimal places at every cut.
            retained_height=max(d for _,_,d in retained)
            for node,length,depth in retained:
                age=(retained_height-depth) if exact_ages else round(retained_height-depth,7)
                active_lengths[node]=length+(age-slab)
                if active_lengths[node]<=0:
                    raise ValueError("a cut produced a nonpositive remaining branch")
    return TreeSections(tree,times,tuple(layers),rows)


def config_single(nt, n, edges, lengths, tips, null, root_edge, qstates, cstates, na,
                  fixes=(), fixlikes=_EMPTY, minimum=1e-6, ntraits=0, prior=_EMPTY):
    qbuf = _states(qstates) if qstates else (_EMPTY,_EMPTY)
    cbuf = _states(cstates) if cstates else (_EMPTY,_EMPTY)
    return [[nt,n,2*nt-1,minimum,null,True,root_edge,1,3,ntraits], edges, lengths, tips,
            fixes,fixlikes,_EMPTY,prior,*qbuf,*cbuf,[na if qstates else 0,na if cstates else 0]]


def compile_sections(sections, states, masks, null, areas, fixes, fixlikes, minimum, ntraits=0, root_prior=_EMPTY, exact_ages=False):
    tree, layers, rows = sections.tree, sections.pieces, sections.rows
    nt, nn, n = tree.tip_count, len(tree.labels), len(states)*max(1,ntraits)
    # Master rows for original tips precede the layer pieces, as in M3.
    ttl = sorted((i for i,r in enumerate(rows) if r["SUBnode.type"] in ("orig_tip","root","internal")),
                 key=lambda i: rows[i]["node"])
    down = sorted((i for i,r in enumerate(rows) if r["piececlass"]=="subtree" and r["SUBnode.type"]!="root"),
                  key=lambda i: rows[i]["node"])
    routes = [{node:(j,k) for j,p in enumerate(layer) for k,node in enumerate(p.nodes)} for layer in layers]
    if any(t not in routes[0] for t in range(nt)):
        raise ValueError("stratified compilation requires all tips to reach the present stratum")
    tips = [(0,*routes[0][t]) for t in range(nt)]
    config = [[n,nt,nn,len(rows),tree.root,len(layers),null],ttl,down,np.asarray(tips).ravel()]
    widths = []
    for i,layer in enumerate(layers):
        mask = masks[i]
        trait_mask = np.tile(mask,max(1,ntraits))
        columns = np.flatnonzero(trait_mask)
        local = tuple(s for s,keep in zip(states,mask) if keep)
        cstates = local[1:] if null else local
        if not cstates:
            raise ValueError(f"stratum {i+1} has no occupied geographic state")
        widths.append(max(map(len,cstates)))
        state_buffers = (_EMPTY,)*4 if ntraits else (*_states(local),*_states(cstates))
        config.extend(([len(layer),areas],columns,mask,*state_buffers))
        for p in layer:
            older = routes[i+1][p.master_root] if i+1<len(layers) else _EMPTY
            if p.branch:
                r = rows[p.rows[0]]
                age,top,bot = r["time_bp"],r["time_top"],r["time_bot"]
                fossil = r["fossils"]
                if exact_ages:
                    ok = bool(np.isfinite(p.duration) and p.duration > 0 and r["edge.length"] >= minimum and not fossil)
                    config.extend(([0,p.sources[0],p.rows[0],p.duration,ok,p.duration,ok,p.duration],older))
                else:
                    downflag,upflag = not (fossil and age>top), age>=top
                    length = p.duration
                    if fossil and top<=age<bot:
                        downflag = upflag = True
                        length -= age-top
                    if r["edge.length"]<minimum:
                        downflag=upflag=False
                    config.extend(([0,p.sources[0],p.rows[0],p.duration,downflag,length,upflag,length],older))
            else:
                ntips=len(p.sources)
                positions = sorted((p.nodes.index(node), k) for k,node in enumerate(fixes)
                                   if node in p.nodes[ntips:])
                fn=[j for j,k in positions]
                fl=fixlikes[[k for j,k in positions],:][:,trait_mask] if fn else _EMPTY
                piece_prior = (root_prior[columns]
                               if i==len(layers)-1 and p.master_root==tree.root and len(root_prior)
                               else _EMPTY)
                config.extend(([1],p.sources,p.rows,p.nodes,
                               [j==ntips for j in range(len(p.nodes))],
                               older if p.duration>0 else _EMPTY,fn,fl))
                config.extend(config_single(ntips,len(columns),p.edges,p.lengths,_EMPTY,null,p.duration,
                                            (),(),areas,fn,fl,minimum,prior=piece_prior))
    return config,widths
