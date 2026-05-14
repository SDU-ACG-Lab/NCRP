"""
io/readers.py – Parsers for Kraken2, PAF, and edge-list files.

All parsers return plain Python structures; no side effects beyond
populating the provided IdMap instances.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, List, Set, Tuple

from ncrp.core.id_map import IdMap

logger = logging.getLogger("ncrp")

# Type aliases
Labels = Dict[int, int]           # read_id → taxon_id
ReadOrder = List[int]              # read IDs in original file order
AdjDict = Dict[int, Dict[int, int]]  # u → {v: overlap}
Lengths = Dict[int, int]           # read_id → read_length


# ---------------------------------------------------------------------------
# Kraken2 reader
# ---------------------------------------------------------------------------

def load_labels_and_order(
    kraken_file: str,
    rid_map: IdMap,
    tax_map: IdMap,
) -> Tuple[Labels, ReadOrder]:
    """
    Parse a Kraken2 output file.

    Parameters
    ----------
    kraken_file:
        Path to Kraken2 TSV (status, read_id, taxon_id, ...).
    rid_map:
        IdMap for read identifiers (mutated in-place).
    tax_map:
        IdMap for taxon identifiers (mutated in-place).

    Returns
    -------
    labels:
        Mapping read_int_id → taxon_int_id for classified reads.
    order:
        All read_int_ids in file order (used for output ordering).
    """
    labels: Labels = {}
    order: ReadOrder = []
    seen: Set[int] = set()

    logger.info(f"Loading Kraken2 file: {kraken_file}")
    with open(kraken_file, "r", encoding="utf-8", errors="ignore") as fh:
        for lineno, line in enumerate(fh, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                logger.debug(f"  Skipping short line {lineno}: {line!r}")
                continue

            status, rid_s, tax_s = parts[0], parts[1], parts[2]
            rid = rid_map.get_int(rid_s)

            if rid not in seen:
                order.append(rid)
                seen.add(rid)

            if status == "C" and tax_s != "0":
                labels[rid] = tax_map.get_int(tax_s)

    logger.info(f"  Reads total: {len(order)}, Classified: {len(labels)}")
    return labels, order


# ---------------------------------------------------------------------------
# PAF reader
# ---------------------------------------------------------------------------

def parse_paf(
    paf_file: str,
    rid_map: IdMap,
    min_overlap: int = 130,
) -> Tuple[AdjDict, Lengths, dict]:
    """
    Parse a PAF (Pairwise mApping Format) overlap file.

    Only records where column[9] (residue matches) >= *min_overlap* are kept.
    Self-loops are discarded. For duplicate pairs the best overlap is kept.

    Returns
    -------
    adj:    Undirected adjacency {u: {v: best_overlap}}.
    lengths: Read lengths observed in the file.
    stats:  Diagnostic counters.
    """
    adj: AdjDict = defaultdict(dict)
    lengths: Lengths = {}
    stats = {"kept_records": 0, "skipped_low_ov": 0, "parse_errors": 0}

    logger.info(f"Parsing PAF: {paf_file}  (min_overlap={min_overlap})")
    with open(paf_file, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                stats["parse_errors"] += 1
                continue

            try:
                q, t = fields[0], fields[5]
                ov = int(fields[9])
                qlen, tlen = int(fields[1]), int(fields[6])
            except (ValueError, IndexError):
                try:
                    ov = max(0, int(fields[3]) - int(fields[2]))
                    q, t = fields[0], fields[5]
                    qlen = tlen = 0
                except (ValueError, IndexError):
                    stats["parse_errors"] += 1
                    continue

            if ov < min_overlap:
                stats["skipped_low_ov"] += 1
                continue

            u = rid_map.get_int(q)
            v = rid_map.get_int(t)
            if u == v:
                continue

            stats["kept_records"] += 1
            prev = adj[u].get(v)
            if prev is None or ov > prev:
                adj[u][v] = ov
                adj[v][u] = ov
            if qlen:
                lengths[u] = qlen
            if tlen:
                lengths[v] = tlen

    logger.info(
        f"  Kept: {stats['kept_records']}  "
        f"Skipped (low ov): {stats['skipped_low_ov']}  "
        f"Parse errors: {stats['parse_errors']}"
    )
    return adj, lengths, stats


# ---------------------------------------------------------------------------
# Edge-list reader
# ---------------------------------------------------------------------------

def parse_edges(
    edge_file: str,
    rid_map: IdMap,
    min_overlap: int = 130,
) -> Tuple[AdjDict, Lengths, dict]:
    """
    Parse a simple 3-column edge file: ``read1 read2 overlap``.

    Returns the same structure as :func:`parse_paf`.
    """
    adj: AdjDict = defaultdict(dict)
    lengths: Lengths = {}
    stats = {"kept_records": 0, "skipped_low_ov": 0, "parse_errors": 0}

    logger.info(f"Parsing edge file: {edge_file}  (min_overlap={min_overlap})")
    with open(edge_file, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 3:
                stats["parse_errors"] += 1
                continue
            r1, r2, ov_s = parts
            try:
                ov = int(ov_s)
            except ValueError:
                stats["parse_errors"] += 1
                continue

            if ov < min_overlap:
                stats["skipped_low_ov"] += 1
                continue

            u = rid_map.get_int(r1)
            v = rid_map.get_int(r2)
            if u == v:
                continue

            stats["kept_records"] += 1
            prev = adj[u].get(v)
            if prev is None or ov > prev:
                adj[u][v] = ov
                adj[v][u] = ov

    logger.info(
        f"  Kept: {stats['kept_records']}  "
        f"Skipped (low ov): {stats['skipped_low_ov']}  "
        f"Parse errors: {stats['parse_errors']}"
    )
    return adj, lengths, stats
