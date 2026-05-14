"""
pipeline.py – High-level pipeline orchestrator.

Ties together all stages in the correct order, using PipelineConfig to
drive decisions.  The pipeline can be used programmatically::

    from ncrp.config import PipelineConfig
    from ncrp.pipeline import run_pipeline

    cfg = PipelineConfig(kraken="reads.kraken2", paf="ovlp.paf")
    run_pipeline(cfg)
"""

from __future__ import annotations

import logging
import time
from typing import Dict

from ncrp.algorithms.graph import build_graph
from ncrp.algorithms.propagation import label_propagation
from ncrp.algorithms.refine import refine_labels
from ncrp.algorithms.rescue import rescue_isolated
from ncrp.config import PipelineConfig
from ncrp.core.id_map import IdMap
from ncrp.io.readers import (
    load_labels_and_order,
    parse_edges,
    parse_paf,
)
from ncrp.io.writers import write_tsv
from ncrp.utils.logging import get_logger, stage_timer

logger = get_logger()


def run_pipeline(cfg: PipelineConfig) -> None:
    """
    Execute the full NCRP pipeline according to *cfg*.

    Stages
    ------
    1. Load Kraken2 labels
    2. Parse overlap file (PAF or edges)
    3. Build compact graph
    4. Refine (Phase A correction + Phase B removal)
    5. Label propagation
    6. Merge results
    7. Rescue isolated reads
    8. Write output TSV

    Parameters
    ----------
    cfg: Fully validated PipelineConfig.
    """
    cfg.validate()
    t_start = time.perf_counter()
    logger.info("=" * 50)
    logger.info("NCRP Pipeline Starting")
    logger.info("=" * 50)

    rid_map = IdMap()
    tax_map = IdMap()

    # ------------------------------------------------------------------
    # Stage 1: Load Kraken2 labels
    # ------------------------------------------------------------------
    with stage_timer(logger, "Stage 1: Load Kraken2"):
        init_labels, order = load_labels_and_order(cfg.kraken, rid_map, tax_map)

    # ------------------------------------------------------------------
    # Stage 2: Parse overlap file
    # ------------------------------------------------------------------
    with stage_timer(logger, "Stage 2: Parse Overlaps"):
        if cfg.overlap_type == "paf":
            adj, lengths, _ = parse_paf(cfg.overlap_file, rid_map, cfg.min_overlap)
        else:
            adj, lengths, _ = parse_edges(cfg.overlap_file, rid_map, cfg.min_overlap)

    # ------------------------------------------------------------------
    # Stage 3: Build graph
    # ------------------------------------------------------------------
    with stage_timer(logger, "Stage 3: Build Graph"):
        graph, _ = build_graph(adj, lengths)

    # ------------------------------------------------------------------
    # Stage 4: Refine labels (Phase A + B)
    # ------------------------------------------------------------------
    with stage_timer(logger, "Stage 4: Refine Labels"):
        refined = refine_labels(
            graph,
            init_labels,
            min_support=cfg.local_vote_min_support,
            skip_correction=cfg.skip_correction,
            skip_removal=cfg.skip_removal,
            keep_single=cfg.keep_single,
            keep_majority=cfg.keep_majority,
        )

    # ------------------------------------------------------------------
    # Stage 5: Label propagation
    # ------------------------------------------------------------------
    with stage_timer(logger, "Stage 5: Label Propagation"):
        if not cfg.skip_propagation:
            propagated = label_propagation(graph, refined)
        else:
            logger.info("Label Propagation: Skipped.")
            propagated = refined

    # ------------------------------------------------------------------
    # Stage 6: Merge (prefer propagated; fall back to refined)
    # ------------------------------------------------------------------
    merged: Dict[int, int] = {}
    for r in order:
        if r in propagated:
            merged[r] = propagated[r]
        elif r in refined:
            merged[r] = refined[r]

    # ------------------------------------------------------------------
    # Stage 7: Rescue isolated reads
    # ------------------------------------------------------------------
    with stage_timer(logger, "Stage 7: Rescue"):
        if not cfg.skip_rescue:
            merged, _ = rescue_isolated(
                cfg.overlap_file,
                cfg.overlap_type,
                rid_map,
                merged,
                graph,
                sig_ratio=cfg.rescue_sig_ratio,
                min_score=cfg.rescue_min_score,
            )
        else:
            logger.info("Rescue: Skipped.")

    # ------------------------------------------------------------------
    # Stage 8: Write output
    # ------------------------------------------------------------------
    with stage_timer(logger, "Stage 8: Write Output"):
        write_tsv(cfg.output, order, merged, rid_map, tax_map)

    total = time.perf_counter() - t_start
    logger.info("=" * 50)
    logger.info(f"Pipeline complete. Total time: {total:.2f}s")
    logger.info("=" * 50)
