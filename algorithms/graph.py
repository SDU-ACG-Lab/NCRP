"""
algorithms/graph.py – Graph construction and finalization.

The core graph representation uses parallel arrays (array.array) for
cache-friendly neighbor iteration, stored as::

    graph[node_id] = (neighbors: array('I'), weights: array('f'))

where weights are overlap lengths normalized by the longest read.
"""

from __future__ import annotations

import logging
from array import array
from typing import Dict, Optional, Tuple

from ncrp.io.readers import AdjDict, Lengths

logger = logging.getLogger("ncrp")

# Graph type alias: node → (neighbor_ids, edge_weights)
Graph = Dict[int, Tuple[array, array]]


def build_graph(
    adj: AdjDict,
    lengths: Lengths,
    default_len: int = 1000,
) -> Tuple[Graph, float]:
    """
    Convert an adjacency dict into the compact array-based graph.

    Edge weights are normalized to [0, 1] by the maximum read length
    observed in *lengths* (or *default_len* if *lengths* is empty).

    Parameters
    ----------
    adj:         Raw adjacency from a reader {u: {v: overlap}}.
    lengths:     Read lengths {read_id: length}.
    default_len: Fallback max length when *lengths* is empty.

    Returns
    -------
    graph:  Compact graph representation.
    l_max:  The normalization constant (max read length).
    """
    l_max: float = float(max(lengths.values())) if lengths else float(default_len)
    graph: Graph = {}

    for u, nbrs in adj.items():
        if not nbrs:
            continue
        neis: array = array("I")
        wts: array = array("f")
        for v, ov in nbrs.items():
            neis.append(v)
            wts.append(ov / l_max)
        graph[u] = (neis, wts)

    node_count = len(graph)
    edge_count = sum(len(neis) for neis, _ in graph.values()) // 2
    logger.info(f"Graph: {node_count} nodes, {edge_count} edges, l_max={l_max:.0f}")
    return graph, l_max
