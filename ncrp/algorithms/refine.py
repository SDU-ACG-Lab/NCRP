"""
algorithms/refine.py – Two-stage label refinement.

Phase A – Correction (Local Consensus)
    For each labeled read, if ALL neighbors share a single consensus label
    that differs from the read's own label, and neighbor count >= min_support,
    correct the label.

Phase B – Removal (Global BFS)
    Run a multi-source BFS from all labeled seeds. For each seed, collect the
    labels of the nearest seeds from other components (the "BFS frontier").
    If any conflict is detected, the seed's label is removed, unless an
    ablation strategy says otherwise.

    Ablation options
    ----------------
    keep_single   (--keep-single)
        Do NOT remove if the nearest frontier has only one label type
        (even if it differs from the seed's own label).

    keep_majority (--keep-majority)
        Do NOT remove if the seed's own label accounts for > 50% of
        frontier votes.
"""

from __future__ import annotations

import logging
from collections import defaultdict, deque
from typing import Dict

from ncrp.algorithms.graph import Graph

logger = logging.getLogger("ncrp")

Labels = Dict[int, int]


# ---------------------------------------------------------------------------
# Phase A – Correction
# ---------------------------------------------------------------------------

def phase_a_correction(
    graph: Graph,
    labels: Labels,
    min_support: int = 2,
) -> Labels:
    """
    Correct labels whose neighborhood is unanimously a different label.

    Returns a *new* dict; the input is not mutated.
    """
    logger.info("Phase A: Correction (Local Consensus)...")
    current = dict(labels)
    corrections = 0

    for s, self_lab in list(current.items()):
        neis, _ = graph.get(s, ((), ()))
        neighbor_labs = [current[u] for u in neis if u in current]

        if not neighbor_labs:
            continue

        unique = set(neighbor_labs)
        if len(unique) == 1:
            consensus = next(iter(unique))
            if len(neighbor_labs) >= min_support and consensus != self_lab:
                current[s] = consensus
                corrections += 1

    logger.info(f"  Corrections: {corrections}")
    return current


# ---------------------------------------------------------------------------
# Phase B – Removal
# ---------------------------------------------------------------------------

def phase_b_removal(
    graph: Graph,
    labels: Labels,
    keep_single: bool = False,
    keep_majority: bool = False,
) -> Labels:
    """
    Remove ambiguous labels via global BFS frontier analysis.

    Each labeled seed expands outward; when two seeds' BFS frontiers meet,
    the opposing label is recorded. Ambiguous seeds (those with conflicting
    frontier evidence) are removed unless a keep strategy overrides.

    Parameters
    ----------
    graph:        The overlap graph.
    labels:       Current label assignments.
    keep_single:  Strategy 1 – retain if frontier is single-label.
    keep_majority: Strategy 2 – retain if self-label has > 50% frontier votes.

    Returns
    -------
    A pruned copy of *labels*.
    """
    logger.info("Phase B: Removal (Global BFS)...")
    if keep_single:
        logger.info("  Strategy: keep_single enabled")
    if keep_majority:
        logger.info("  Strategy: keep_majority enabled")

    active_seeds = list(labels.keys())
    if not active_seeds:
        return dict(labels)

    max_id = max(
        max(graph.keys(), default=0),
        max(active_seeds),
    ) + 1

    # BFS state
    owner = [-1] * max_id   # which seed owns this node
    dist = [-1] * max_id    # BFS distance from seed

    # Per-seed: label counts at the closest BFS frontier distance
    best_d: Dict[int, int] = {s: 10 ** 15 for s in active_seeds}
    best_counts: Dict[int, Dict[int, int]] = {
        s: defaultdict(int) for s in active_seeds
    }

    q: deque = deque()
    for s in active_seeds:
        if s < max_id:
            owner[s] = s
            dist[s] = 0
            q.append(s)

    while q:
        u = q.popleft()
        own_u = owner[u]
        du = dist[u]
        neis, _ = graph.get(u, ((), ()))

        for v in neis:
            if v >= max_id:
                continue
            if owner[v] == -1:
                owner[v] = own_u
                dist[v] = du + 1
                q.append(v)
            else:
                own_v = owner[v]
                if own_v == own_u:
                    continue
                # Boundary edge: two different seeds meet
                cand_dist = du + dist[v] + 1
                for seed, opp_seed in ((own_u, own_v), (own_v, own_u)):
                    if seed not in best_d:
                        continue
                    opp_label = labels[opp_seed]
                    if cand_dist < best_d[seed]:
                        best_d[seed] = cand_dist
                        best_counts[seed] = defaultdict(int)
                        best_counts[seed][opp_label] = 1
                    elif cand_dist == best_d[seed]:
                        best_counts[seed][opp_label] += 1

    refined = dict(labels)
    removed = 0

    for s, self_lab in labels.items():
        counts = best_counts.get(s)
        if not counts:
            continue  # isolated or no frontier → safe

        unique_count = len(counts)
        total_votes = sum(counts.values())
        self_votes = counts.get(self_lab, 0)

        # Detect conflict
        has_conflict = (unique_count > 1) or (self_lab not in counts)
        if not has_conflict:
            continue  # unanimous self-supporting frontier

        # --- Ablation keep strategies ---
        should_keep = False

        if keep_single and unique_count == 1:
            should_keep = True  # frontier is uniform (even if different)

        if keep_majority and total_votes > 0:
            ratio = self_votes / total_votes
            if ratio > 0.5:
                should_keep = True

        if not should_keep:
            del refined[s]
            removed += 1

    logger.info(f"  Removed: {removed} ambiguous labels")
    return refined


# ---------------------------------------------------------------------------
# Unified entry point
# ---------------------------------------------------------------------------

def refine_labels(
    graph: Graph,
    init_labels: Labels,
    min_support: int = 2,
    skip_correction: bool = False,
    skip_removal: bool = False,
    keep_single: bool = False,
    keep_majority: bool = False,
) -> Labels:
    """
    Run Phase A and/or Phase B refinement.

    Parameters
    ----------
    graph, init_labels: Input graph and seed labels.
    min_support:        Phase A minimum neighbor count.
    skip_correction:    Bypass Phase A.
    skip_removal:       Bypass Phase B.
    keep_single:        Phase B ablation: keep_single strategy.
    keep_majority:      Phase B ablation: keep_majority strategy.

    Returns
    -------
    Refined label dict.
    """
    current = init_labels

    if not skip_correction:
        current = phase_a_correction(graph, current, min_support)
    else:
        logger.info("Phase A: Skipped.")

    if not skip_removal:
        current = phase_b_removal(graph, current, keep_single, keep_majority)
    else:
        logger.info("Phase B: Skipped.")

    return current
