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
