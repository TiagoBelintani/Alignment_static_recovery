#!/usr/bin/env python3
"""
recover_overall_stats.py

Calculates per-taxon statistics from multiple alignments in a folder.

Example usage:
    python recover_overall_stats.py --input input_folder --output overall_stats.csv
    python recover_overall_stats.py --input alignments --output stats.csv --format fasta --plots

Outputs:
- CSV file with per-taxon statistics
- (Optional) PNG plots in a plots directory

Author: general purpose / adaptable
"""

import argparse
import csv
import os
import sys
from pathlib import Path
from collections import defaultdict

from Bio import AlignIO

# matplotlib is optional
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Retrieve overall statistics per taxon from multiple alignments."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Folder containing alignment files."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file. Example: overall_stats.csv"
    )
    parser.add_argument(
        "--format",
        default="fasta",
        help="Alignment format read by Biopython. Default: fasta"
    )
    parser.add_argument(
        "--extensions",
        nargs="+",
        default=[".fa", ".fasta", ".faa", ".fas", ".aln"],
        help="Accepted alignment file extensions. Default: .fa .fasta .faa .fas .aln"
    )
    parser.add_argument(
        "--missing-chars",
        default="-?",
        help="Characters treated as missing/gap. Default: '-?'"
    )
    parser.add_argument(
        "--plots",
        action="store_true",
        help="Generate PNG plots with summary statistics."
    )
    parser.add_argument(
        "--plots-dir",
        default="plots_stats",
        help="Directory where plots will be saved. Default: plots_stats"
    )
    parser.add_argument(
        "--min-taxa",
        type=int,
        default=4,
        help="Minimum number of valid taxa per column to evaluate informative sites. Default: 4"
    )
    return parser.parse_args()


def is_informative_column(column_chars):
    """
    Parsimony-informative site rule:
    at least two states, each present in at least two sequences.
    """
    counts = {}
    for char in column_chars:
        counts[char] = counts.get(char, 0) + 1

    repeated_states = sum(1 for v in counts.values() if v > 1)
    return repeated_states >= 2


def collect_alignment_files(input_dir, extensions):
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Folder not found: {input_dir}")
    if not input_path.is_dir():
        raise NotADirectoryError(f"Path is not a folder: {input_dir}")

    files = [
        f for f in sorted(input_path.iterdir())
        if f.is_file() and f.suffix.lower() in {ext.lower() for ext in extensions}
    ]
    return files


def analyze_alignments(files, aln_format, missing_chars, min_taxa=4):
    taxon_stats = defaultdict(lambda: {
        "count": 0,
        "total_length": 0,
        "total_gap_fraction": 0.0,
        "total_informative_sites_present": 0,
        "total_missing_chars": 0,
        "total_non_missing_chars": 0,
    })

    alignment_summary = []
    total_files = 0
    failed_files = []

    missing_set = set(missing_chars)

    for file_path in files:
        try:
            alignment = AlignIO.read(str(file_path), aln_format)
        except Exception as e:
            failed_files.append((str(file_path), str(e)))
            continue

        total_files += 1
        aln_length = alignment.get_alignment_length()
        n_taxa = len(alignment)

        if aln_length == 0 or n_taxa == 0:
            failed_files.append((str(file_path), "Empty alignment"))
            continue

        informative_columns = set()

        for i in range(aln_length):
            column = [str(rec.seq[i]).upper() for rec in alignment]
            valid_chars = [c for c in column if c not in missing_set]

            if len(valid_chars) < min_taxa:
                continue

            if is_informative_column(valid_chars):
                informative_columns.add(i)

        aln_informative_count = len(informative_columns)

        alignment_summary.append({
            "file": file_path.name,
            "n_taxa": n_taxa,
            "length": aln_length,
            "informative_sites": aln_informative_count,
            "informative_fraction": aln_informative_count / aln_length if aln_length else 0.0
        })

        for rec in alignment:
            seq = str(rec.seq).upper()
            rec_id = rec.id

            missing_count = sum(1 for c in seq if c in missing_set)
            non_missing_count = aln_length - missing_count
            gap_fraction = missing_count / aln_length if aln_length else 0.0
            informative_present = sum(1 for i in informative_columns if seq[i] not in missing_set)

            s = taxon_stats[rec_id]
            s["count"] += 1
            s["total_length"] += aln_length
            s["total_gap_fraction"] += gap_fraction
            s["total_informative_sites_present"] += informative_present
            s["total_missing_chars"] += missing_count
            s["total_non_missing_chars"] += non_missing_count

    return taxon_stats, alignment_summary, total_files, failed_files


def write_csv(output_file, taxon_stats):
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Taxon",
            "Alignments_Present",
            "Mean_Length",
            "Mean_Gaps",
            "Mean_Informative_Sites",
            "Missing_Total",
            "NonMissing_Total",
            "Global_Missing_Fraction"
        ])

        for taxon in sorted(taxon_stats):
            stats = taxon_stats[taxon]
            count = stats["count"]

            if count == 0:
                continue

            total_chars = stats["total_missing_chars"] + stats["total_non_missing_chars"]
            global_missing_fraction = (
                stats["total_missing_chars"] / total_chars if total_chars else 0.0
            )

            writer.writerow([
                taxon,
                count,
                round(stats["total_length"] / count, 2),
                round(stats["total_gap_fraction"] / count, 4),
                round(stats["total_informative_sites_present"] / count, 2),
                stats["total_missing_chars"],
                stats["total_non_missing_chars"],
                round(global_missing_fraction, 4)
            ])


def write_alignment_summary(output_file, alignment_summary):
    summary_path = Path(output_file).with_name(Path(output_file).stem + "_per_alignment.csv")

    with summary_path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "File",
            "Number_Taxa",
            "Length",
            "Informative_Sites",
            "Informative_Site_Fraction"
        ])

        for row in alignment_summary:
            writer.writerow([
                row["file"],
                row["n_taxa"],
                row["length"],
                row["informative_sites"],
                round(row["informative_fraction"], 4)
            ])

    return summary_path


def generate_plots(taxon_stats, plots_dir):
    if not HAS_MATPLOTLIB:
        print("[warning] matplotlib is not installed; plots will not be generated.", file=sys.stderr)
        return

    plots_path = Path(plots_dir)
    plots_path.mkdir(parents=True, exist_ok=True)

    taxa = sorted(taxon_stats.keys())
    counts = [taxon_stats[t]["count"] for t in taxa]
    mean_gaps = [
        taxon_stats[t]["total_gap_fraction"] / taxon_stats[t]["count"]
        if taxon_stats[t]["count"] else 0.0
        for t in taxa
    ]
    mean_informative = [
        taxon_stats[t]["total_informative_sites_present"] / taxon_stats[t]["count"]
        if taxon_stats[t]["count"] else 0.0
        for t in taxa
    ]

    plt.figure(figsize=(12, 6))
    plt.bar(taxa, counts)
    plt.xticks(rotation=90)
    plt.ylabel("Number of alignments")
    plt.xlabel("Taxon")
    plt.title("Alignments per taxon")
    plt.tight_layout()
    plt.savefig(plots_path / "alignments_per_taxon.png", dpi=300)
    plt.close()

    plt.figure(figsize=(12, 6))
    plt.bar(taxa, mean_gaps)
    plt.xticks(rotation=90)
    plt.ylabel("Mean gaps / missing")
    plt.xlabel("Taxon")
    plt.title("Mean gaps per taxon")
    plt.tight_layout()
    plt.savefig(plots_path / "mean_gaps_per_taxon.png", dpi=300)
    plt.close()

    plt.figure(figsize=(12, 6))
    plt.bar(taxa, mean_informative)
    plt.xticks(rotation=90)
    plt.ylabel("Mean informative sites")
    plt.xlabel("Taxon")
    plt.title("Mean informative sites per taxon")
    plt.tight_layout()
    plt.savefig(plots_path / "mean_informative_sites_per_taxon.png", dpi=300)
    plt.close()


def main():
    args = parse_args()

    try:
        files = collect_alignment_files(args.input, args.extensions)
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)

    if not files:
        print("[error] No alignment files found with the provided extensions.", file=sys.stderr)
        sys.exit(1)

    taxon_stats, alignment_summary, total_files, failed_files = analyze_alignments(
        files=files,
        aln_format=args.format,
        missing_chars=args.missing_chars,
        min_taxa=args.min_taxa
    )

    if not taxon_stats:
        print("[error] No valid alignment was processed.", file=sys.stderr)
        if failed_files:
            print("\nFailed files:", file=sys.stderr)
            for fname, err in failed_files:
                print(f"  - {fname}: {err}", file=sys.stderr)
        sys.exit(1)

    write_csv(args.output, taxon_stats)
    aln_summary_path = write_alignment_summary(args.output, alignment_summary)

    if args.plots:
        generate_plots(taxon_stats, args.plots_dir)

    print("\nProcessing completed successfully.")
    print(f"Valid alignments processed: {total_files}")
    print(f"Taxa found: {len(taxon_stats)}")
    print(f"Main CSV: {args.output}")
    print(f"Per-alignment summary: {aln_summary_path}")

    if args.plots:
        if HAS_MATPLOTLIB:
            print(f"Plots saved in: {args.plots_dir}")
        else:
            print("Plots not generated because matplotlib is not installed.")

    if failed_files:
        print(f"\nFailed files: {len(failed_files)}")
        for fname, err in failed_files[:10]:
            print(f"  - {fname}: {err}")
        if len(failed_files) > 10:
            print("  ...")


if __name__ == "__main__":
    main()
