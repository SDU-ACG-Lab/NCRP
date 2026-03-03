"""
main.py – CLI entry point for NCRP.

Usage
-----
    python main.py --kraken reads.kraken2 --paf overlaps.paf [options]
    python main.py --kraken reads.kraken2 --edges overlaps.edges [options]
"""

import argparse
import sys

from ncrp.config import PipelineConfig
from ncrp.pipeline import run_pipeline
from ncrp.utils.logging import get_logger

logger = get_logger()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="ncrp",
        description="Graph-based taxonomic label correction and propagation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Required ---
    p.add_argument("--kraken", required=True, metavar="FILE",
                   help="Kraken2 classification output file.")

    # --- Overlap source (mutually exclusive) ---
    ovlp = p.add_mutually_exclusive_group(required=True)
    ovlp.add_argument("--paf", metavar="FILE",
                      help="PAF overlap file.")
    ovlp.add_argument("--edges", metavar="FILE",
                      help="3-column edge file (read1 read2 overlap).")

    # --- Output ---
    p.add_argument("--output", default="final.tsv", metavar="FILE",
                   help="Output TSV path.")

    # --- Graph build ---
    p.add_argument("--min-overlap", type=int, default=130,
                   help="Minimum overlap length to include an edge.")

    # --- Refine ---
    p.add_argument("--local-vote-min-support", type=int, default=2,
                   help="Min. neighbor count required for Phase A correction.")

    # --- Ablation: stage skips ---
    ablation = p.add_argument_group("Ablation – stage skips")
    ablation.add_argument("--skip-correction", action="store_true",
                          help="Skip Phase A (correction).")
    ablation.add_argument("--skip-removal", action="store_true",
                          help="Skip Phase B (removal).")
    ablation.add_argument("--skip-propagation", action="store_true",
                          help="Skip label propagation.")
    ablation.add_argument("--skip-rescue", action="store_true",
                          help="Skip rescue step.")

    # --- Ablation: removal strategies ---
    strategies = p.add_argument_group("Ablation – Phase B removal strategies")
    strategies.add_argument("--keep-single", action="store_true",
                            help="Keep label if BFS frontier is single-type "
                                 "(even if it differs from own label).")
    strategies.add_argument("--keep-majority", action="store_true",
                            help="Keep label if self-support > 50%% of frontier votes.")

    # --- Rescue ---
    rescue = p.add_argument_group("Rescue parameters")
    rescue.add_argument("--rescue-sig-ratio", type=float, default=1.2,
                        help="Min. best/second-best score ratio for rescue acceptance.")
    rescue.add_argument("--rescue-min-score", type=int, default=50,
                        help="Min. raw overlap score for rescue acceptance.")

    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    cfg = PipelineConfig(
        kraken=args.kraken,
        paf=args.paf,
        edges=args.edges,
        output=args.output,
        min_overlap=args.min_overlap,
        local_vote_min_support=args.local_vote_min_support,
        skip_correction=args.skip_correction,
        skip_removal=args.skip_removal,
        skip_propagation=args.skip_propagation,
        skip_rescue=args.skip_rescue,
        keep_single=args.keep_single,
        keep_majority=args.keep_majority,
        rescue_sig_ratio=args.rescue_sig_ratio,
        rescue_min_score=args.rescue_min_score,
    )

    try:
        run_pipeline(cfg)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(f"Pipeline error: {exc}")
        sys.exit(1)
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        sys.exit(130)


if __name__ == "__main__":
    main()
