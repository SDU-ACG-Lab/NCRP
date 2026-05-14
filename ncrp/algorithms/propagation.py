"""
algorithms/propagation.py – Layered BFS label propagation.

Starting from all labeled seeds, BFS expands layer by layer.  At each
layer, an unlabeled node inherits the weighted-majority label from its
already-labeled neighbors in the previous layer.

The weight of each neighbor vote is the normalized overlap length stored
in the graph edge.
"""

from __future__ import annotations

import logging
from collections import defaultdict, deque
from typing import Dict

from ncrp.algorithms.graph import Graph

logger = logging.getLogger("ncrp")

Labels = Dict[int, int]


def label_propagation(graph: Graph, init_labels: Labels) -> Labels:
    """
    Spread labels across the graph via weighted majority BFS.

    Parameters
    ----------
    graph:       Overlap graph.
    init_labels: Seed labels after refinement.

    Returns
    -------
    Extended label mapping that includes propagated assignments.
    """
    logger.info("Label Propagation...")

    if not init_labels:
        logger.warning("  No seed labels – propagation skipped.")
        return dict(init_labels)

    max_id = max(
        max(graph.keys(), default=0),
        max(init_labels.keys(), default=0),
    ) + 1

    labels = dict(init_labels)
    dist = [10 ** 15] * max_id

    # Seed all labeled nodes at distance 0
    q: deque = deque()
    for s in labels:
        if s < max_id:
            dist[s] = 0
            q.append(s)

    # BFS to compute distances
    while q:
        u = q.popleft()
        d = dist[u]
        neis, _ = graph.get(u, ((), ()))
        for v in neis:
            if v < max_id and dist[v] == 10 ** 15:
                dist[v] = d + 1
                q.append(v)

    # Bucket unlabeled nodes by BFS layer
    buckets: Dict[int, list] = defaultdict(list)
    max_dist = 0
    for v in graph:
        dv = dist[v] if v < max_id else 10 ** 15
        if 0 < dv < 10 ** 15:
            buckets[dv].append(v)
            if dv > max_dist:
                max_dist = dv

    # Layer-by-layer weighted majority vote
    added = 0
    for d in range(1, max_dist + 1):
        for v in buckets[d]:
            neis, wts = graph.get(v, ((), ()))
            acc: Dict[int, float] = {}
            for i, u in enumerate(neis):
                if u < max_id and dist[u] == d - 1:
                    lab = labels.get(u)
                    if lab is not None:
                        acc[lab] = acc.get(lab, 0.0) + wts[i]
            if acc:
                # Tie-break: higher weight wins; on tie, lower taxon ID wins
                best = max(acc.items(), key=lambda x: (x[1], -x[0]))[0]
                labels[v] = best
                added += 1

    logger.info(f"  Propagated: {added} reads")
    return labels
