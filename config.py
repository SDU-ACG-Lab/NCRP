"""
config.py – Pipeline configuration dataclass with defaults.
All CLI arguments map 1-to-1 to fields here.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class PipelineConfig:
    # --- I/O ---
    kraken: str = ""
    paf: Optional[str] = None
    edges: Optional[str] = None
    output: str = "final.tsv"

    # --- Graph Build ---
    min_overlap: int = 130

    # --- Refine ---
    local_vote_min_support: int = 2

    # --- Ablation: Stage Skips ---
    skip_correction: bool = False
    skip_removal: bool = False
    skip_propagation: bool = False
    skip_rescue: bool = False

    # --- Ablation: Removal Strategies ---
    keep_single: bool = False    # Phase B: keep if neighbors are single-type
    keep_majority: bool = False  # Phase B: keep if self-support > 50%

    # --- Rescue ---
    rescue_sig_ratio: float = 1.2
    rescue_min_score: int = 50

    def validate(self) -> None:
        """Raise ValueError for invalid combinations."""
        if not self.kraken:
            raise ValueError("--kraken is required.")
        if not self.paf and not self.edges:
            raise ValueError("Either --paf or --edges must be provided.")
        if self.paf and self.edges:
            raise ValueError("--paf and --edges are mutually exclusive.")
        if self.min_overlap < 0:
            raise ValueError("--min-overlap must be non-negative.")
        if self.local_vote_min_support < 1:
            raise ValueError("--local-vote-min-support must be >= 1.")
        if self.rescue_sig_ratio <= 0:
            raise ValueError("--rescue-sig-ratio must be > 0.")

    @property
    def overlap_file(self) -> str:
        return self.paf if self.paf else self.edges  # type: ignore

    @property
    def overlap_type(self) -> str:
        return "paf" if self.paf else "edges"
