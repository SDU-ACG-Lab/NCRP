# NCRP

A graph-based taxonomic label correction and propagation pipeline for metagenomic reads.

## Pipeline Stages

1. **Graph Build** – Construct overlap graph from PAF/edges file
2. **Refine Phase A** – Correction via Strict Local Consensus
3. **Refine Phase B** – Removal via Global BFS (with ablation options)
4. **Propagation** – Layered BFS label spreading
5. **Rescue** – Recover unlabeled isolated reads via overlap evidence

## Project Structure

```
ncrp/
├── core/
│   ├── __init__.py
│   └── id_map.py          # ID compression utility
├── io/
│   ├── __init__.py
│   ├── readers.py         # Kraken, PAF, Edges file parsers
│   └── writers.py         # TSV output writer
├── algorithms/
│   ├── __init__.py
│   ├── graph.py           # Graph construction & finalization
│   ├── refine.py          # Phase A (Correction) + Phase B (Removal)
│   ├── propagation.py     # Layered BFS label propagation
│   └── rescue.py          # Rescue step for isolated reads
├── utils/
│   ├── __init__.py
│   └── logging.py         # Structured logger
├── config.py              # Dataclass-based config / defaults
├── pipeline.py            # High-level orchestrator
├── main.py                # CLI entry point
└── README.md
```

## Input Formats

### Kraken2 Classification File (`--kraken`)

Standard Kraken2 output format, tab-separated. At minimum, the first three columns are used:

| Column | Field | Description |
|---|---|---|
| 1 | Status | `C` = classified, `U` = unclassified |
| 2 | Read ID | Unique read identifier (e.g. `S1R0`) |
| 3 | Taxon ID | NCBI taxonomy ID (`0` for unclassified reads) |
| 4+ | — | Ignored by NCRP |

Example:
```
U	S1R0	0	2506	0:581 156889:5 0:1886
C	S1R1	40269	2226	0:187 40269:5 0:1 40269:5 0:1472
C	S1R2	1238	2683	0:985 1238:4 0:1326 1238:12 0:322
```

### PAF Overlap File (`--paf`)

Standard [PAF (Pairwise mApping Format)](https://github.com/lh3/miniasm/blob/master/PAF.md), tab-separated with at least 10 columns. The overlap length is taken from column 10 (number of residue matches). Only records with overlap >= `--min-overlap` (default 130) are kept.

| Column | Field | Used by NCRP |
|---|---|---|
| 1 | Query name | Read ID |
| 6 | Target name | Read ID |
| 10 | Residue matches | Overlap length |

Example:
```
S1R0	2506	911	1791	+	S1R1416196	5135	821	1668	115	880	...
S1R1	2226	902	2213	+	S1R1447390	2793	243	1567	149	1328	...
```

### Edges File (`--edges`)

Simple 3-column space-separated file: `read1 read2 overlap`. Used as an alternative to PAF when you already have a precomputed overlap graph.

Example:
```
S1R0 S1R1416196 880
S1R1 S1R1447390 1328
S1R1 S1R130557 1779
```

## Usage

```bash
# 方式一：直接运行（在仓库根目录执行）
python -m ncrp.main --kraken reads.kraken2 --paf overlaps.paf --output final.tsv

# 方式二：安装后使用命令行
pip install -e .
ncrp --kraken reads.kraken2 --paf overlaps.paf --output final.tsv

# Ablation flags
python -m ncrp.main --kraken reads.kraken2 --paf overlaps.paf \
    --skip-correction \
    --keep-single \
    --keep-majority \
    --skip-rescue
```

## Arguments

| Flag | Default | Description |
|---|---|---|
| `--kraken` | required | Kraken2 classification file |
| `--paf` | — | PAF overlap file (mutually exclusive with --edges) |
| `--edges` | — | Edges file (r1 r2 overlap) |
| `--output` | final.tsv | Output TSV path |
| `--min-overlap` | 130 | Minimum overlap threshold |
| `--local-vote-min-support` | 2 | Min neighbors for Phase A correction |
| `--skip-correction` | False | Skip Phase A |
| `--skip-removal` | False | Skip Phase B |
| `--skip-propagation` | False | Skip propagation |
| `--skip-rescue` | False | Skip rescue |
| `--keep-single` | False | Phase B: keep if neighbors are single-type |
| `--keep-majority` | False | Phase B: keep if self-support > 50% |
| `--rescue-sig-ratio` | 1.2 | Rescue: significance ratio threshold |
| `--rescue-min-score` | 50 | Rescue: minimum overlap score |
