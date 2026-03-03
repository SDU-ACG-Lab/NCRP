# NCRP: a graph-based framework for improving long-read metagenomic classification 

## Installation

### With Conda

We recommend installing NCRP in a dedicated conda environment.

#### From yaml file

Create a file called `environment.yml` in the repository root:

yaml
name: ncrp
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.8
  - matplotlib
Then create and activate the environment:
```
conda env create -f environment.yml
conda activate ncrp
```
### Input file
NCRP expects two main inputs:

#### 1.classification result file

The standard per-read output:
```
C   read_0001   1345   ...
U   read_0002   0       ...
```
Only the three columns are used:

- column 1 : status(C or U)
- column 2 : read ID
- column 3 : taxid(0 for unclassified)
#### 2.overlap graph

Building an overlap graph with minimap2(example). For long reads(ONT/PacBio HiFi), a typical command is : 
```
minimap2 -x ava-ont -t 32 reads.fq reads.fq > overlaps.paf
```
or
```
minimap2 -x ava-pb -t 32 reads.fq reads.fq > overlaps.paf
```

You can provide the graph in two formats:
- minimap2 PAF(--paf) : NCRP parses the file and builds an overlap graph,keeping the best overlap for each read pair.
- Sample edge list(--edges) : a 3-column,tab-delimitered file:
```
readA   readB   350
readB   readC   420
...
```
where the third column is the overlap length. Exactly one of --paf or --edges must be provided.
### Testing your NCRP installation
After installing NCRP and preparing a small test dataset, you can run:
```
cd NCRP
python ncrp_main.py \
  --kraken tests/example.kraken \
  --paf tests/example.paf \
  --output tests/ncrp_labels.tsv
```
If the command finishes successfully and tests/ncrp_labels.tsv is created, your installation works.

To see all available options:
```
python ncrp_main.py --help
```
### Base Usage
```
python ncrp_main.py \
  --kraken path/to/reads.kraken \
  --paf path/to/overlaps.paf \
  --min-overlap 130 \
  --hist-bin 100 \
  --output ncrp_final_labels.tsv
```
### Important command-line options
```
--kraken            Kraken2 per-read output file (required)
--paf               minimap2 PAF overlap file (optional; mutually exclusive with --edges)
--edges             simple edge list "read1 read2 overlap" (optional; mutually exclusive with --paf)
--min-overlap       minimum overlap length (bp) to keep an edge (default: 130)
--hist-bin          bin width (bp) for textual overlap histogram (default: 100)
--output            output file for final labels (default: final_labels.tsv)
```
### Main output files
The main output is ncrp_final_labels.tsv(or the filename you specify with --output), with three columns
```
status   read_id   taxid
```
In addition, NCRP prints to stout:
basic statistics on overlaps(number of raw records, filter edges, unique kept edges, etc.)
