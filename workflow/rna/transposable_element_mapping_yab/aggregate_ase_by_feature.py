#!/usr/bin/env python

import argparse
import sys
from typing import Dict

import pandas as pd
import pyranges as pr


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate GATK ASEReadCounter output to gene/transcript/gene_name "
            "level using a GTF annotation."
        )
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="GTF annotation file with gene_id, transcript_id, and gene_name attributes.",
    )
    parser.add_argument(
        "--ase",
        required=True,
        help="GATK ASEReadCounter output file (tab-delimited).",
    )
    parser.add_argument(
        "--level",
        choices=["gene_id", "transcript_id", "gene_name"],
        default="gene_id",
        help="Feature level for aggregation (default: gene_id).",
    )
    parser.add_argument(
        "--feature-type",
        dest="feature_type",
        choices=["exon", "transcript", "gene"],
        default="exon",
        help="GTF feature type to use for intervals (default: exon).",
    )
    parser.add_argument(
        "--agg",
        choices=["mean", "sum"],
        default="mean",
        help="Aggregation function for counts across SNPs within a feature (default: mean).",
    )
    parser.add_argument(
        "--min-total-count",
        type=int,
        default=0,
        help="Drop SNPs with totalCount < this before aggregation (default: 0).",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        help="Prefix for output files. Will create "
             "<prefix>.snp_by_<level>.tsv and <prefix>.ase_by_<level>.tsv",
    )
    return parser.parse_args()


def parse_gtf_attributes(attr_str: str) -> Dict[str, str]:
    """
    Parse the 9th column of a GTF/GFF-like attributes string into a dict.

    Handles both:
        key "value"; key2 "value2";
    and:
        key=value; key2=value2;
    styles as seen in different annotations.
    """
    attrs: Dict[str, str] = {}
    if attr_str is None:
        return attrs
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        # Prefer space-separated first, fallback to '='
        if " " in field:
            key, val = field.split(" ", 1)
        elif "=" in field:
            key, val = field.split("=", 1)
        else:
            # malformed, skip
            continue
        key = key.strip()
        val = val.strip().strip('"')
        attrs[key] = val
    return attrs


def load_features_from_gtf(gtf_path: str, feature_type: str, level: str) -> pd.DataFrame:
    """
    Load features from a GTF, keep only rows of given feature_type (e.g. exon),
    extract gene_id / transcript_id / gene_name attributes, and build BED-like
    DataFrame suitable for PyRanges.

    Coordinates are converted from 1-based closed (GTF) to 0-based half-open (BED)
    for interval operations:

        start0 = start1 - 1
        end0 = end1
    """
    records = []
    with open(gtf_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, frame, attrs_str = parts
            if ftype != feature_type:
                continue
            try:
                start1 = int(start)
                end1 = int(end)
            except ValueError:
                continue

            attrs = parse_gtf_attributes(attrs_str)
            gene_id = attrs.get("gene_id", None)
            transcript_id = attrs.get("transcript_id", None)
            gene_name = attrs.get("gene_name", None)

            # 0-based, half-open for PyRanges
            start0 = start1 - 1
            end0 = end1

            records.append(
                {
                    "Chromosome": chrom,
                    "Start": start0,
                    "End": end0,
                    "Strand": strand,
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "gene_name": gene_name,
                }
            )

    if not records:
        sys.stderr.write(
            f"[ERROR] No records found in GTF {gtf_path} with feature_type={feature_type}\n"
        )
        sys.exit(1)

    df = pd.DataFrame.from_records(records)

    # Determine feature_id used for grouping / merging
    if level not in df.columns:
        sys.stderr.write(
            f"[ERROR] Requested level '{level}' not found in GTF attributes.\n"
        )
        sys.exit(1)

    df["feature_id"] = df[level]

    # Drop rows where feature_id is missing
    df = df.dropna(subset=["feature_id"])

    return df


def merge_feature_intervals(features_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping/adjacent intervals for each feature_id (e.g. over exons
    of multiple isoforms) to get a union of intervals.

    Returns a DataFrame with columns:
        Chromosome, Start, End, feature_id, gene_id, transcript_id, gene_name
    """
    # Keep mapping from feature_id to a representative gene_id / transcript_id / gene_name
    mapping = (
        features_df.groupby("feature_id")[["gene_id", "transcript_id", "gene_name"]]
        .agg("first")
        .reset_index()
    )

    # Only interval columns + feature_id for merging
    gr = pr.PyRanges(features_df[["Chromosome", "Start", "End", "feature_id"]])

    # Merge intervals by feature_id
    gr_merged = gr.merge(by="feature_id")

    merged_df = gr_merged.as_df()

    # Join back the annotation mapping
    merged_df = merged_df.merge(mapping, on="feature_id", how="left")

    return merged_df


def load_ase_table(ase_path: str, min_total_count: int) -> pd.DataFrame:
    """
    Load GATK ASEReadCounter output as a DataFrame, add 0-based Start/End
    for interval operations, and optionally filter on totalCount.

    Expected columns (GATK ASEReadCounter default):
        contig, position, variantID, refAllele, altAllele,
        refCount, altCount, totalCount, lowMAPQDepth, lowBaseQDepth,
        rawDepth, otherBases, improperPairs
    """
    ase = pd.read_table(ase_path, sep="\t", comment="#")

    # Basic sanity checks
    required_cols = ["contig", "position", "refCount", "altCount", "totalCount"]
    for col in required_cols:
        if col not in ase.columns:
            sys.stderr.write(
                f"[ERROR] Column '{col}' not found in ASEReadCounter table: {ase_path}\n"
            )
            sys.exit(1)

    # Filter on totalCount if requested
    if min_total_count > 0:
        ase = ase.loc[ase["totalCount"] >= min_total_count].copy()

    if ase.empty:
        sys.stderr.write(
            "[WARNING] ASE table is empty after applying min_total_count filter.\n"
        )

    # Add BED-style coordinates for 1bp SNP positions
    ase["Chromosome"] = ase["contig"]
    ase["Start"] = ase["position"] - 1  # 0-based
    ase["End"] = ase["position"]        # half-open

    return ase


def intersect_variants_features(
    ase_df: pd.DataFrame, features_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Intersect per-site ASE counts with feature intervals using PyRanges.

    Returns a DataFrame where each row is an overlapping SNP-feature pair, with:
        Chromosome, Start, End, position, variantID, refCount, altCount, totalCount, ...
        feature_id, gene_id, transcript_id, gene_name
    """
    if ase_df.empty:
        return pd.DataFrame()

    gr_variants = pr.PyRanges(ase_df)
    gr_features = pr.PyRanges(features_df)

    # Intersect SNPs with features (keep only overlapping pairs)
    joined = gr_variants.join(gr_features)  # <- removed how="inner"

    joined_df = joined.as_df()

    return joined_df



def aggregate_by_feature(
    joined_df: pd.DataFrame, level: str, agg_func: str
) -> pd.DataFrame:
    """
    Aggregate SNP-level counts into feature-level counts using mean or sum.

    Output columns include:
        feature_id, gene_id, transcript_id, gene_name, n_snps,
        refCount_<agg>, altCount_<agg>, totalCount_<agg>, ref_fraction
    """
    if joined_df.empty:
        sys.stderr.write(
            "[WARNING] No SNPs overlapped with features. Summary table will be empty.\n"
        )
        return pd.DataFrame(
            columns=[
                "feature_id",
                "gene_id",
                "transcript_id",
                "gene_name",
                "n_snps",
                f"refCount_{agg_func}",
                f"altCount_{agg_func}",
                f"totalCount_{agg_func}",
                "ref_fraction",
            ]
        )

    ref_col_name = f"refCount_{agg_func}"
    alt_col_name = f"altCount_{agg_func}"
    tot_col_name = f"totalCount_{agg_func}"

    # groupby feature_id
    agg_spec = {
        "gene_id": ("gene_id", "first"),
        "transcript_id": ("transcript_id", "first"),
        "gene_name": ("gene_name", "first"),
        "n_snps": ("position", "count"),
        ref_col_name: ("refCount", agg_func),
        alt_col_name: ("altCount", agg_func),
        tot_col_name: ("totalCount", agg_func),
    }

    summary = joined_df.groupby("feature_id").agg(**agg_spec).reset_index()

    # Allelic fraction based on aggregated counts
    denom = summary[ref_col_name] + summary[alt_col_name]
    summary["ref_fraction"] = summary[ref_col_name] / denom.replace(0, pd.NA)

    # This is redundant but nice to have the requested 'level' column explicitly
    summary[level] = summary["feature_id"]

    return summary


def main():
    args = parse_args()

    # 1. Load and process GTF
    sys.stderr.write(
        f"[INFO] Loading GTF from {args.gtf} using feature_type={args.feature_type}, level={args.level}\n"
    )
    features_df = load_features_from_gtf(args.gtf, args.feature_type, args.level)

    sys.stderr.write("[INFO] Merging overlapping intervals per feature_id...\n")
    merged_features_df = merge_feature_intervals(features_df)

    # 2. Load ASE table
    sys.stderr.write(
        f"[INFO] Loading ASEReadCounter table from {args.ase} "
        f"(min_total_count={args.min_total_count})\n"
    )
    ase_df = load_ase_table(args.ase, args.min_total_count)

    # 3. Intersect variants with features
    sys.stderr.write("[INFO] Intersecting SNPs with features...\n")
    joined_df = intersect_variants_features(ase_df, merged_features_df)

    # 4. Write intermediate SNP-feature table
    inter_path = f"{args.output_prefix}.snp_by_{args.level}.tsv"
    sys.stderr.write(f"[INFO] Writing SNP-by-feature table to {inter_path}\n")
    joined_df.to_csv(inter_path, sep="\t", index=False)

    # 5. Aggregate per feature
    sys.stderr.write(
        f"[INFO] Aggregating counts by {args.level} using agg={args.agg}\n"
    )
    summary_df = aggregate_by_feature(joined_df, args.level, args.agg)

    out_path = f"{args.output_prefix}.ase_by_{args.level}.tsv"
    sys.stderr.write(f"[INFO] Writing aggregated ASE table to {out_path}\n")
    summary_df.to_csv(out_path, sep="\t", index=False)

    sys.stderr.write("[INFO] Done.\n")


if __name__ == "__main__":
    main()
