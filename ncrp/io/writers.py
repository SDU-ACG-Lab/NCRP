"""
io/writers.py – Output writers for the pipeline.
"""

from __future__ import annotations

import logging
from typing import Dict, List

from ncrp.core.id_map import IdMap

logger = logging.getLogger("ncrp")


def write_tsv(
    output_path: str,
    order: List[int],
    labels: Dict[int, int],
    rid_map: IdMap,
    tax_map: IdMap,
) -> None:
    """
    Write classification results as a 3-column TSV.

    Columns: ``status \\t read_id \\t taxon_id``

    Classified reads use status ``C``; unclassified reads use ``U`` with
    taxon ``0``.

    Parameters
    ----------
    output_path: Destination file path.
    order:       Read IDs in desired output order.
    labels:      Final read_id → taxon_id mapping.
    rid_map:     IdMap to recover original read string IDs.
    tax_map:     IdMap to recover original taxon string IDs.
    """
    logger.info(f"Writing output: {output_path}")
    classified = 0
    unclassified = 0

    with open(output_path, "w") as fh:
        for r in order:
            r_str = rid_map.lookup(r)
            l_int = labels.get(r)
            if l_int is None:
                fh.write(f"U\t{r_str}\t0\n")
                unclassified += 1
            else:
                fh.write(f"C\t{r_str}\t{tax_map.lookup(l_int)}\n")
                classified += 1

    total = classified + unclassified
    pct = 100.0 * classified / total if total else 0.0
    logger.info(
        f"  Classified: {classified}  Unclassified: {unclassified}  "
        f"({pct:.1f}% classified)"
    )
