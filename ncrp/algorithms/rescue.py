"""
algorithms/rescue.py – Rescue step for isolated / unlabeled reads.

Reads that are unlabeled after propagation AND have no edges in the main
graph (isolated nodes) are candidates for rescue.  The rescue step re-scans
the raw overlap file and assigns a label if there is a dominant labeled
neighbor with sufficient overlap evidence.

Acceptance criteria
-------------------
1. Best candidate score >= rescue_min_score.
2. Best score >= second-best score * rescue_sig_ratio   (significance gap).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, Tuple

from ncrp.algorithms.graph import Graph
from ncrp.core.id_map import IdMap

logger = logging.getLogger("ncrp")

Labels = Dict[int, int]


def _collect_candidates(
    file_path: str,
    file_type: str,
    rid_map: IdMap,
    targets: set,
    final_labels: Labels,
) -> Dict[int, Dict[int, int]]:
    """
    Scan *file_path* and accumulate overlap evidence for *targets*.

    Returns candidates[target_id][label] = total_overlap_score.
    """
    candidates: Dict[int, Dict[int, int]] = defaultdict(lambda: defaultdict(int))

    def _add(u_s: str, v_s: str, ov: int) -> None:
        u = rid_map.get_int(u_s)
        v = rid_map.get_int(v_s)
        if u == v:
            return
        if u in targets and v in final_labels:
            candidates[u][final_labels[v]] += ov
        elif v in targets and u in final_labels:
            candidates[v][final_labels[u]] += ov

    with open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            if file_type == "paf":
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                try:
                    ov = int(fields[9])
                except ValueError:
                    continue
                _add(fields[0], fields[5], ov)
            elif file_type == "edges":
                parts = line.split()
                if len(parts) != 3:
                    continue
                try:
                    ov = int(parts[2])
                except ValueError:
                    continue
                _add(parts[0], parts[1], ov)

    return candidates


def rescue_isolated(
    file_path: str,
    file_type: str,
    rid_map: IdMap,
    final_labels: Labels,
    graph: Graph,
    sig_ratio: float = 1.2,
    min_score: int = 50,
) -> Tuple[Labels, int]:
    """
    Attempt to rescue unlabeled, graph-isolated reads.

    Parameters
    ----------
    file_path:     Path to overlap file (PAF or edges).
    file_type:     ``'paf'`` or ``'edges'``.
    rid_map:       IdMap for read IDs.
    final_labels:  Current label assignments (after propagation).
    graph:         Overlap graph (used to identify isolated nodes).
    sig_ratio:     Minimum ratio of best/second-best score for acceptance.
    min_score:     Minimum raw score for the best candidate.

    Returns
    -------
    updated_labels: New dict with rescued reads added.
    rescued_count:  Number of reads newly classified.
    """
    logger.info("Rescue Step...")

    targets: set = set()
    for r in range(len(rid_map)):
        if r not in final_labels:
            if r not in graph or not graph[r][0]:
                targets.add(r)

    if not targets:
        logger.info("  No isolated unclassified reads – rescue skipped.")
        return final_labels, 0

    logger.info(f"  Candidates for rescue: {len(targets)}")
    candidates = _collect_candidates(file_path, file_type, rid_map, targets, final_labels)

    updated = dict(final_labels)
    rescued = 0

    for rid, scores in candidates.items():
        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        best_label, best_score = ranked[0]

        if best_score < min_score:
            continue
        if len(ranked) > 1 and best_score < ranked[1][1] * sig_ratio:
            continue

        updated[rid] = best_label
        rescued += 1

    logger.info(f"  Rescued: {rescued}")
    return updated, rescued
