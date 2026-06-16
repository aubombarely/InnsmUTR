#!/usr/bin/env python3
"""
compare_utr_annotation.py  —  InnsmUTR

Benchmarks a UTR-annotated GFF3 (e.g. a PASA-updated gene set) against a
reference / gold-standard GFF3 (e.g. full-length transcript evidence or a
manually curated annotation).

UTR length per transcript is read directly from explicit five_prime_UTR /
three_prime_UTR feature lines (summed per transcript, so spliced/multi-exon
UTRs are handled correctly). This script does not derive UTRs from
mRNA/CDS spans — if a GFF3 only has CDS and exon records, first run AGAT
to add UTR features (see README), then compare the AGAT output.

UTR presence is scored per side as a binary call (present / absent against
the reference), and boundary accuracy is scored among shared-presence
transcripts as the UTR length difference (target − ref), which is already
strand-corrected.

Transcripts are matched between target and reference either by shared ID
(default — use when both files describe the same gene/transcript set,
e.g. before/after running PASA on one annotation) or by CDS-span overlap
(use when target and reference come from independent annotations with
unrelated IDs).

Usage:
    # Match by shared transcript ID
    python compare_utr_annotation.py \\
        -t pasa_updated.gff3 -r reference.gff3 \\
        --outdir comparisons/ --prefix arath --tolerance 20

    # Match by CDS overlap (IDs differ between target and reference)
    python compare_utr_annotation.py \\
        -t pasa_updated.gff3 -r reference.gff3 \\
        --match_by overlap --min_cds_overlap 0.9 \\
        --outdir comparisons/ --prefix arath

Output:
    {outdir}/{prefix}_comparison.tsv   — per-transcript UTR detail
    {outdir}/{prefix}_metrics.tsv      — presence P/R/F1 + boundary accuracy
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

ATTR_RE = re.compile(r"(\w+)=([^;]+)")


def parse_attributes(attr_str: str) -> dict:
    return dict(ATTR_RE.findall(attr_str))


AGAT_HINT = (
    "If this GFF3 only has CDS and exon records, generate UTR features with "
    "AGAT first (see README), then re-run this script on the AGAT output."
)


# ── GFF3 loading ────────────────────────────────────────────────────────────
def load_transcripts(gff3_path: Path, id_attr: str = "ID",
                     mrna_feature: str = "mRNA") -> pd.DataFrame:
    """
    Parse a GFF3 and return one row per coding transcript with its CDS
    span (for matching) and total UTR length per side, read directly from
    five_prime_UTR / three_prime_UTR feature lines (summed per transcript
    so spliced/multi-exon UTRs are counted correctly). Coordinates are
    0-based half-open.

    Columns: id, chrom, strand, cds_start, cds_end, utr5_len, utr3_len
    """
    mrna: dict = {}
    cds: dict = {}
    utr5: dict = {}
    utr3: dict = {}
    n_utr_lines = 0

    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attrs = f[:9]
            attrs_d = parse_attributes(attrs)
            s, e = int(start) - 1, int(end)

            if feature == mrna_feature:
                tid = attrs_d.get(id_attr)
                if tid is None:
                    continue
                mrna[tid] = {"chrom": chrom, "strand": strand}
            elif feature == "CDS":
                parent = attrs_d.get("Parent")
                if parent is None:
                    continue
                if parent not in cds:
                    cds[parent] = [s, e]
                else:
                    cds[parent][0] = min(cds[parent][0], s)
                    cds[parent][1] = max(cds[parent][1], e)
            elif feature == "five_prime_UTR":
                parent = attrs_d.get("Parent")
                if parent is None:
                    continue
                n_utr_lines += 1
                utr5[parent] = utr5.get(parent, 0) + (e - s)
            elif feature == "three_prime_UTR":
                parent = attrs_d.get("Parent")
                if parent is None:
                    continue
                n_utr_lines += 1
                utr3[parent] = utr3.get(parent, 0) + (e - s)

    if n_utr_lines == 0:
        raise ValueError(
            f"{gff3_path}: no five_prime_UTR/three_prime_UTR features found. "
            f"{AGAT_HINT}"
        )

    rows = []
    for tid, m in mrna.items():
        c = cds.get(tid)
        rows.append({
            "id": tid,
            "chrom": m["chrom"],
            "strand": m["strand"],
            "cds_start": c[0] if c else None,
            "cds_end": c[1] if c else None,
            "utr5_len": utr5.get(tid, 0),
            "utr3_len": utr3.get(tid, 0),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        raise ValueError(f"{gff3_path}: no '{mrna_feature}' features found")

    n_no_cds = df["cds_start"].isna().sum()
    if n_no_cds:
        log.warning(f"{gff3_path}: {n_no_cds:,} transcript(s) with no CDS "
                    f"— dropped (non-coding?)")
        df = df.dropna(subset=["cds_start", "cds_end"]).copy()
        df["cds_start"] = df["cds_start"].astype(int)
        df["cds_end"] = df["cds_end"].astype(int)

    log.info(f"{gff3_path}: {len(df):,} coding transcripts loaded")
    return df


# ── Matching ────────────────────────────────────────────────────────────────
def match_by_id(target_df: pd.DataFrame, ref_df: pd.DataFrame) -> pd.DataFrame:
    merged = ref_df.merge(target_df, on="id", how="outer",
                          suffixes=("_ref", "_target"), indicator=True)
    n_ref_only = (merged["_merge"] == "left_only").sum()
    n_target_only = (merged["_merge"] == "right_only").sum()
    if n_ref_only:
        log.warning(f"{n_ref_only:,} reference transcript(s) have no matching "
                    f"target ID — excluded")
    if n_target_only:
        log.warning(f"{n_target_only:,} target transcript(s) have no matching "
                    f"reference ID — excluded")
    matched = merged[merged["_merge"] == "both"].drop(columns="_merge")
    log.info(f"Matched {len(matched):,} transcripts by shared ID")
    return matched


def match_by_overlap(target_df: pd.DataFrame, ref_df: pd.DataFrame,
                     min_cds_overlap: float) -> pd.DataFrame:
    """
    Match each reference transcript to the target transcript with the
    largest CDS-span overlap on the same chrom/strand, above
    min_cds_overlap (as a fraction of the reference CDS span). Uses the
    CDS span (min start, max end), not exon-resolved CDS — adequate for
    matching whole gene models, not for fine boundary comparison.
    """
    import pybedtools

    ref_bed = pybedtools.BedTool.from_dataframe(
        ref_df[["chrom", "cds_start", "cds_end", "id"]]
    )
    tgt_bed = pybedtools.BedTool.from_dataframe(
        target_df[["chrom", "cds_start", "cds_end", "id"]]
    )
    intersect = ref_bed.intersect(tgt_bed, wao=True)

    best: dict = {}
    for feat in intersect:
        ref_id = feat.fields[3]
        tgt_id = feat.fields[7]
        bp = int(feat.fields[-1])
        if tgt_id == "." or bp <= 0:
            continue
        if ref_id not in best or bp > best[ref_id][1]:
            best[ref_id] = (tgt_id, bp)

    ref_indexed = ref_df.set_index("id")
    tgt_indexed = target_df.set_index("id")
    rows = []
    n_unmatched = len(ref_df) - len(best)

    for ref_id, (tgt_id, bp) in best.items():
        ref_row = ref_indexed.loc[ref_id]
        tgt_row = tgt_indexed.loc[tgt_id]
        ref_cds_len = ref_row["cds_end"] - ref_row["cds_start"]
        if ref_row["strand"] != tgt_row["strand"] or ref_cds_len <= 0 or \
           bp / ref_cds_len < min_cds_overlap:
            n_unmatched += 1
            continue
        row = {"id": ref_id}
        for col in ref_df.columns:
            if col == "id":
                continue
            row[f"{col}_ref"] = ref_row[col]
            row[f"{col}_target"] = tgt_row[col]
        rows.append(row)

    log.info(f"Matched {len(rows):,} transcripts by CDS overlap "
             f"(min_cds_overlap={min_cds_overlap}); {n_unmatched:,} reference "
             f"transcripts unmatched")
    pybedtools.cleanup()
    return pd.DataFrame(rows)


# ── Comparison & metrics ───────────────────────────────────────────────────
def classify_presence(ref_has: bool, target_has: bool) -> str:
    if ref_has and target_has:
        return "TP"
    if ref_has and not target_has:
        return "FN"
    if not ref_has and target_has:
        return "FP"
    return "TN"


def build_comparison(matched: pd.DataFrame, tolerance: int) -> pd.DataFrame:
    df = matched.copy()

    for side in ("5", "3"):
        ref_col, tgt_col = f"utr{side}_len_ref", f"utr{side}_len_target"
        df[f"ref_has_{side}utr"] = df[ref_col] > 0
        df[f"target_has_{side}utr"] = df[tgt_col] > 0
        df[f"utr{side}_call"] = [
            classify_presence(r, t)
            for r, t in zip(df[f"ref_has_{side}utr"], df[f"target_has_{side}utr"])
        ]
        # Length difference already encodes the strand-corrected boundary
        # offset, since the inner (CDS) boundary is assumed unchanged.
        df[f"utr{side}_boundary_offset_bp"] = df[tgt_col] - df[ref_col]
        df[f"utr{side}_within_tolerance"] = (
            df[f"utr{side}_call"].eq("TP") &
            (df[f"utr{side}_boundary_offset_bp"].abs() <= tolerance)
        )

    cols = ["id", "chrom_ref", "strand_ref",
           "utr5_len_ref", "utr5_len_target", "utr5_call",
           "utr5_boundary_offset_bp", "utr5_within_tolerance",
           "utr3_len_ref", "utr3_len_target", "utr3_call",
           "utr3_boundary_offset_bp", "utr3_within_tolerance"]
    out = df[cols].rename(columns={"chrom_ref": "chrom", "strand_ref": "strand"})
    return out


def _presence_metrics(df: pd.DataFrame, side: str) -> dict:
    calls = df[f"utr{side}_call"].value_counts()
    tp, fp, fn = calls.get("TP", 0), calls.get("FP", 0), calls.get("FN", 0)
    tn = calls.get("TN", 0)
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
         if (precision + recall) and not pd.isna(precision) and not pd.isna(recall)
         else float("nan"))

    tp_rows = df[df[f"utr{side}_call"] == "TP"]
    offsets = tp_rows[f"utr{side}_boundary_offset_bp"]
    boundary_accuracy = (tp_rows[f"utr{side}_within_tolerance"].mean()
                         if len(tp_rows) else float("nan"))

    return {
        "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "precision": precision, "recall": recall, "f1": f1,
        "boundary_accuracy_rate": boundary_accuracy,
        "mean_signed_offset_bp": offsets.mean() if len(offsets) else float("nan"),
        "mean_abs_offset_bp": offsets.abs().mean() if len(offsets) else float("nan"),
        "median_abs_offset_bp": offsets.abs().median() if len(offsets) else float("nan"),
    }


def write_metrics(comparison: pd.DataFrame, tolerance: int,
                  match_by: str, outpath: Path) -> None:
    m5 = _presence_metrics(comparison, "5")
    m3 = _presence_metrics(comparison, "3")

    with open(outpath, "w") as fh:
        fh.write(f"# InnsmUTR UTR annotation comparison metrics\n")
        fh.write(f"# match_by={match_by}\ttolerance_bp={tolerance}\n")
        fh.write(f"# n_matched_transcripts={len(comparison)}\n")
        fh.write("#\n")
        fh.write("side\tTP\tFP\tFN\tTN\tprecision\trecall\tf1\t"
                "boundary_accuracy_rate\tmean_signed_offset_bp\t"
                "mean_abs_offset_bp\tmedian_abs_offset_bp\n")
        for side_label, m in (("5utr", m5), ("3utr", m3)):
            fh.write(
                f"{side_label}\t{m['TP']}\t{m['FP']}\t{m['FN']}\t{m['TN']}\t"
                f"{m['precision']:.4f}\t{m['recall']:.4f}\t{m['f1']:.4f}\t"
                f"{m['boundary_accuracy_rate']:.4f}\t"
                f"{m['mean_signed_offset_bp']:.2f}\t"
                f"{m['mean_abs_offset_bp']:.2f}\t"
                f"{m['median_abs_offset_bp']:.2f}\n"
            )
    log.info(f"Metrics written to {outpath}")
    log.info(f"5' UTR — precision={m5['precision']:.3f} recall={m5['recall']:.3f} "
             f"f1={m5['f1']:.3f} boundary_accuracy={m5['boundary_accuracy_rate']:.3f}")
    log.info(f"3' UTR — precision={m3['precision']:.3f} recall={m3['recall']:.3f} "
             f"f1={m3['f1']:.3f} boundary_accuracy={m3['boundary_accuracy_rate']:.3f}")


# ── Entry point ────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Benchmark a UTR-annotated GFF3 against a reference GFF3"
    )
    parser.add_argument("-t", "--target", required=True,
                        help="Target GFF3 (e.g. PASA-updated annotation)")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference / gold-standard GFF3")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--prefix", required=True,
                        help="Prefix for output file names")
    parser.add_argument("--match_by", choices=["id", "overlap"], default="id",
                        help="Transcript matching strategy (default: id)")
    parser.add_argument("--id_attr", default="ID",
                        help="GFF3 attribute used as transcript ID (default: ID)")
    parser.add_argument("--mrna_feature", default="mRNA",
                        help="GFF3 feature type for transcripts (default: mRNA)")
    parser.add_argument("--min_cds_overlap", type=float, default=0.9,
                        help="Min fraction of reference CDS span that must be "
                             "covered for an overlap match (default: 0.9)")
    parser.add_argument("--tolerance", type=int, default=20,
                        help="UTR boundary tolerance in bp (default: 20)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    target_df = load_transcripts(Path(args.target), args.id_attr, args.mrna_feature)
    ref_df = load_transcripts(Path(args.reference), args.id_attr, args.mrna_feature)

    if args.match_by == "id":
        matched = match_by_id(target_df, ref_df)
    else:
        matched = match_by_overlap(target_df, ref_df, args.min_cds_overlap)

    if matched.empty:
        log.error("No transcripts could be matched between target and reference")
        sys.exit(1)

    comparison = build_comparison(matched, args.tolerance)

    comparison_path = outdir / f"{args.prefix}_comparison.tsv"
    comparison.to_csv(comparison_path, sep="\t", index=False)
    log.info(f"Per-transcript comparison written to {comparison_path}")

    metrics_path = outdir / f"{args.prefix}_metrics.tsv"
    write_metrics(comparison, args.tolerance, args.match_by, metrics_path)


if __name__ == "__main__":
    main()
