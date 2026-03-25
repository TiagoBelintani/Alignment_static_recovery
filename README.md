# Recover Overall Stats

A Python script to compute per-taxon statistics from multiple sequence alignments.

The program reads a directory containing alignment files and generates summary statistics for each taxon across all alignments.

## Features

- Calculates statistics per taxon across multiple alignments
- Identifies parsimony-informative sites
- Outputs summary tables in CSV format
- Optional generation of summary plots
- Flexible support for different alignment formats

## Requirements

Python ≥ 3.7

Required libraries:

- Biopython

Optional:

- matplotlib (for plots)

Install dependencies:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
Usage

Basic usage:
```bash
python recover_overall_stats.py \
--input alignments_folder \
--output overall_stats.csv
```

Example with plots:
````bash
python recover_overall_stats.py \
--input alignments \
--output stats.csv \
--format fasta \
--plots
```

Arguments

```bash
Argument	Description
--input	Directory containing alignment files
--output	Output CSV file
--format	Alignment format (default: fasta)
--extensions	Accepted file extensions
--missing-chars	Characters treated as gaps/missing
--plots	Generate plots
--plots-dir	Directory for plots
--min-taxa	Minimum taxa per column to test informative sites
```
Output files

Main CSV

Per-taxon statistics including:

number of alignments
mean alignment length
mean gap fraction
informative sites
missing data statistics
Alignment summary CSV

Statistics per alignment:

number of taxa
alignment length
number of informative sites
Plots (optional)

If --plots is enabled, the script generates:

alignments per taxon
mean gap fraction per taxon
mean informative sites per taxon
Example workflow

python recover_overall_stats.py \
--input matrix_95_ocup \
--output stats_geral.csv \
--plots

License

MIT License

Author Tiago Belintani 2025  *Brave the Sun*
