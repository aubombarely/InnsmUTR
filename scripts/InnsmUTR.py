#!/usr/bin/env python3
"""
InnsmUTR.py — Add UTRs to a genome annotation using full-length transcript
evidence.

Aligns transcript sequences (cDNA, IsoSeq, Trinity/StringTie assemblies) to
the genome, identifies regions that extend beyond existing CDS annotations, and
adds five_prime_utr / three_prime_utr features to the GFF3.

Pipeline:
    1. Align transcripts to genome (minimap2 or GMAP)
    2. Parse alignments → per-transcript exon blocks
    3. Match aligned transcripts to gene models by CDS overlap
    4. Extract UTR intervals (transcript exons beyond the CDS boundary)
    5. Write updated GFF3 with UTR features + mRNA/gene coordinate updates
    6. Write per-mRNA stats table

Output directory layout:
    {output}/
    ├── results/
    │   ├── mod01_alignment_{prefix}.stats.tsv   Alignment summary stats
    │   ├── mod02_utr_{prefix}.gff3              Updated annotation with UTRs
    │   ├── mod02_utr_stats_{prefix}.tsv         Per-mRNA UTR stats
    │   └── {prefix}.run_summary.json
    ├── workdir/                                  Intermediate alignment files
    └── logs/
        ├── Run_InnsmUTR.log
        └── {prefix}.emissions.csv (if codecarbon installed)
"""

import argparse
import csv
import getpass
import json
import os
import platform
import re
import resource
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

VERSION = "v0.1.0"

# ── Logging ───────────────────────────────────────────────────────────────────

_LOG_FH = None


def _log(msg: str) -> None:
    ts   = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr)
    if _LOG_FH is not None:
        print(line, file=_LOG_FH, flush=True)


def _banner(title: str) -> None:
    bar = "─" * (len(title) + 4)
    _log(f"┌{bar}┐")
    _log(f"│  {title}  │")
    _log(f"└{bar}┘")


# ── Tool helpers ──────────────────────────────────────────────────────────────

def _checkpoint(path: Path, label: str, force: bool) -> bool:
    if not force and path.exists() and path.stat().st_size > 0:
        _log(f"  [checkpoint] {label} — {path.name} already exists, skipping")
        return True
    return False


def _require_tool(name: str) -> str:
    tool = shutil.which(name)
    if tool is None:
        print(
            f"ERROR: '{name}' not found in PATH.\n"
            f"       Install with:  conda install -c bioconda {name}",
            file=sys.stderr,
        )
        sys.exit(1)
    return tool


def _run(
    cmd: list,
    capture_stdout: bool = False,
    env: dict = None,
    cwd: Path = None,
) -> subprocess.CompletedProcess:
    _log(f"  $ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE if capture_stdout else None,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        cwd=cwd,
    )
    if result.returncode != 0:
        print(
            f"ERROR: command failed (exit {result.returncode}):\n"
            f"{result.stderr[-3000:]}",
            file=sys.stderr,
        )
        sys.exit(1)
    return result


# ── GFF3 parsing ──────────────────────────────────────────────────────────────

def _parse_attrs(attr_str: str) -> dict:
    attrs = {}
    for item in attr_str.strip().split(";"):
        item = item.strip()
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def load_gff3(path: Path) -> tuple:
    """
    Parse GFF3 and return:
      genes          : {gene_id  → {chrom, strand, start, end}}
      mrnas          : {mrna_id  → {gene_id, chrom, strand, start, end,
                                    cds_intervals, exon_intervals}}
      existing_utrs  : set of mrna_ids that already carry UTR features
      records        : list of raw input lines (reproduced verbatim in output)
    """
    genes         = {}
    mrnas         = {}
    cds_by_mrna   = defaultdict(list)
    exon_by_mrna  = defaultdict(list)
    existing_utrs = set()
    records       = []

    with open(path) as fh:
        for line in fh:
            records.append(line)
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attr_str = cols
            start = int(start)
            end   = int(end)
            attrs = _parse_attrs(attr_str)

            if ftype == "gene":
                gid = attrs.get("ID", "")
                if gid:
                    genes[gid] = {
                        "chrom": chrom, "strand": strand,
                        "start": start, "end": end,
                    }

            elif ftype in ("mRNA", "transcript"):
                mid = attrs.get("ID", "")
                gid = attrs.get("Parent", "")
                if mid:
                    mrnas[mid] = {
                        "gene_id": gid,
                        "chrom": chrom, "strand": strand,
                        "start": start, "end": end,
                        "cds_intervals": [],
                        "exon_intervals": [],
                    }

            elif ftype == "CDS":
                parent = attrs.get("Parent", "")
                for pid in parent.split(","):
                    pid = pid.strip()
                    if pid:
                        cds_by_mrna[pid].append((start, end))

            elif ftype == "exon":
                parent = attrs.get("Parent", "")
                for pid in parent.split(","):
                    pid = pid.strip()
                    if pid:
                        exon_by_mrna[pid].append((start, end))

            elif ftype in (
                "five_prime_UTR", "three_prime_UTR",
                "five_prime_utr", "three_prime_utr",
            ):
                parent = attrs.get("Parent", "")
                for pid in parent.split(","):
                    pid = pid.strip()
                    if pid:
                        existing_utrs.add(pid)

    for mid, mrna in mrnas.items():
        mrna["cds_intervals"]  = sorted(cds_by_mrna.get(mid, []))
        mrna["exon_intervals"] = sorted(exon_by_mrna.get(mid, []))

    return genes, mrnas, existing_utrs, records


# ── SAM / CIGAR helpers ───────────────────────────────────────────────────────

def _cigar_exon_blocks(ref_start_0: int, cigar: str) -> list:
    """
    Return exon blocks as (start, end) tuples (1-based, inclusive, GFF3 coords).
    ref_start_0 is the SAM POS field converted to 0-based.
    N operations mark intron boundaries; D stays within the current exon block.
    """
    blocks      = []
    pos         = ref_start_0
    block_start = None

    for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar):
        length = int(m.group(1))
        op     = m.group(2)
        if op in ("M", "=", "X", "D"):
            if block_start is None:
                block_start = pos
            pos += length
        elif op == "N":                          # intron: close current exon block
            if block_start is not None:
                blocks.append((block_start + 1, pos))   # convert to 1-based
                block_start = None
            pos += length
        # I, S, H, P: do not advance reference position

    if block_start is not None:
        blocks.append((block_start + 1, pos))

    return blocks


def _aligned_query_bases(cigar: str) -> int:
    """Read bases consumed by the alignment (M/I/=/X ops), excluding soft clips."""
    return sum(int(m.group(1)) for m in re.finditer(r"(\d+)([MI=X])", cigar))


def _total_query_bases(cigar: str) -> int:
    """Total read length implied by CIGAR (M/I/S/=/X ops)."""
    return sum(int(m.group(1)) for m in re.finditer(r"(\d+)([MIS=X])", cigar))


# ── Module 1: Transcript alignment ───────────────────────────────────────────

def run_minimap2(
    genome: Path, transcripts: Path, workdir: Path, threads: int, force: bool
) -> Path:
    sam = workdir / "transcripts.sam"
    if _checkpoint(sam, "minimap2", force):
        return sam
    _require_tool("minimap2")
    _run(
        [
            "minimap2", "-ax", "splice", "--secondary=no", "-C", "5",
            "-t", str(threads),
            str(genome), str(transcripts),
            "-o", str(sam),
        ]
    )
    _log(f"  minimap2 SAM → {sam.name}")
    return sam


def run_gmap(
    genome: Path, transcripts: Path, workdir: Path, threads: int, force: bool
) -> Path:
    sam     = workdir / "transcripts.sam"
    if _checkpoint(sam, "GMAP", force):
        return sam
    _require_tool("gmap_build")
    _require_tool("gmap")
    db_dir  = workdir / "gmap_db"
    db_name = "genome"

    if force or not db_dir.exists():
        db_dir.mkdir(parents=True, exist_ok=True)
        _run(["gmap_build", "-D", str(db_dir), "-d", db_name, str(genome)])
        _log(f"  GMAP database → {db_dir.name}/")
    else:
        _log("  [checkpoint] GMAP database already exists, skipping build")

    result = _run(
        [
            "gmap", "-D", str(db_dir), "-d", db_name,
            "--format=samse", "--npaths=1",
            "-t", str(threads),
            str(transcripts),
        ],
        capture_stdout=True,
    )
    sam.write_text(result.stdout)
    _log(f"  GMAP SAM → {sam.name}")
    return sam


def parse_sam_alignments(
    sam_path: Path, min_identity: float, min_coverage: float
) -> dict:
    """
    Parse SAM file (minimap2 or GMAP output). Only primary alignments are kept
    (secondary flag 0x100 and supplementary flag 0x800 are discarded).

    Returns alignments_by_chrom: {chrom → [aln_dict, ...]}
    Each aln_dict carries: tx_id, chrom, strand, exon_blocks (1-based
    inclusive), identity, coverage, aln_start, aln_end.
    """
    alignments_by_chrom = defaultdict(list)
    n_total = 0
    n_pass  = 0

    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue

            tx_id = cols[0]
            flag  = int(cols[1])
            chrom = cols[2]
            pos   = int(cols[3]) - 1   # convert to 0-based
            cigar = cols[5]
            seq   = cols[9]

            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue
            if chrom == "*" or cigar == "*":
                continue

            n_total += 1
            strand = "-" if (flag & 0x10) else "+"

            nm = 0
            for tag in cols[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag[5:])
                    break

            exon_blocks   = _cigar_exon_blocks(pos, cigar)
            if not exon_blocks:
                continue

            aln_query_len = _aligned_query_bases(cigar)
            query_len     = len(seq) if seq != "*" else _total_query_bases(cigar)
            if query_len == 0 or aln_query_len == 0:
                continue

            identity = max(0.0, 1.0 - nm / aln_query_len)
            coverage = aln_query_len / query_len

            if identity < min_identity or coverage < min_coverage:
                continue

            n_pass += 1
            alignments_by_chrom[chrom].append(
                {
                    "tx_id":       tx_id,
                    "chrom":       chrom,
                    "strand":      strand,
                    "exon_blocks": exon_blocks,
                    "identity":    round(identity, 4),
                    "coverage":    round(coverage, 4),
                    "aln_start":   exon_blocks[0][0],
                    "aln_end":     exon_blocks[-1][1],
                }
            )

    _log(f"  {n_total} primary alignments parsed, {n_pass} passed filters")
    return dict(alignments_by_chrom)


def write_alignment_stats(
    alignments_by_chrom: dict,
    min_identity: float,
    min_coverage: float,
    path: Path,
) -> None:
    all_aln  = [a for v in alignments_by_chrom.values() for a in v]
    n_aln    = len(all_aln)
    n_chroms = len(alignments_by_chrom)
    mean_id  = round(sum(a["identity"] for a in all_aln) / n_aln, 4) if n_aln else 0
    mean_cov = round(sum(a["coverage"] for a in all_aln) / n_aln, 4) if n_aln else 0

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["alignments_passing_filters",  n_aln])
        w.writerow(["chromosomes_with_alignments",  n_chroms])
        w.writerow(["min_identity_threshold",        min_identity])
        w.writerow(["min_coverage_threshold",        min_coverage])
        w.writerow(["mean_identity",                 mean_id])
        w.writerow(["mean_coverage",                 mean_cov])

    _log(f"  Alignment stats → {path.name}")


# ── Module 2: UTR extraction ──────────────────────────────────────────────────

def _clip_upstream(block: tuple, cds_start: int):
    """Return the part of block strictly before cds_start (1-based), or None."""
    s, e = block
    if s >= cds_start:
        return None
    return (s, min(e, cds_start - 1))


def _clip_downstream(block: tuple, cds_end: int):
    """Return the part of block strictly after cds_end (1-based), or None."""
    s, e = block
    if e <= cds_end:
        return None
    return (max(s, cds_end + 1), e)


def _utr_blocks_for_alignment(
    mrna: dict, aln: dict, min_utr_len: int, max_utr_len: int
) -> dict:
    """
    Given an mRNA and an overlapping transcript alignment, derive candidate
    5' and 3' UTR blocks (filtered by length thresholds).
    """
    cds_start = mrna["cds_intervals"][0][0]
    cds_end   = mrna["cds_intervals"][-1][1]
    strand    = mrna["strand"]

    up_blocks = [
        c for blk in aln["exon_blocks"]
        for c in (_clip_upstream(blk, cds_start),) if c
    ]
    dn_blocks = [
        c for blk in aln["exon_blocks"]
        for c in (_clip_downstream(blk, cds_end),) if c
    ]

    five_blocks  = up_blocks if strand == "+" else dn_blocks
    three_blocks = dn_blocks if strand == "+" else up_blocks

    def _filter(blocks):
        total = sum(e - s + 1 for s, e in blocks)
        return blocks if min_utr_len <= total <= max_utr_len else []

    return {
        "five_prime_utr":  _filter(five_blocks),
        "three_prime_utr": _filter(three_blocks),
    }


def extract_utrs(
    mrnas: dict,
    existing_utrs: set,
    alignments_by_chrom: dict,
    min_utr_len: int,
    max_utr_len: int,
) -> dict:
    """
    For each mRNA (without pre-existing UTR features) find the best-scoring
    overlapping transcript alignment and extract UTR blocks.

    Scoring: identity × coverage + small bonus for alignments that extend
    beyond the CDS (favouring transcripts with UTR evidence over pure
    CDS-matching ones).

    Returns {mrna_id → {'five_prime_utr': [...], 'three_prime_utr': [...],
                         'supporting_tx': str, 'identity': float, 'coverage': float}}
    """
    utrs        = {}
    n_with_5    = 0
    n_with_3    = 0
    n_skipped   = 0

    for mid, mrna in mrnas.items():
        if not mrna["cds_intervals"]:
            continue
        if mid in existing_utrs:
            n_skipped += 1
            continue

        chrom  = mrna["chrom"]
        strand = mrna["strand"]
        cds_s  = mrna["cds_intervals"][0][0]
        cds_e  = mrna["cds_intervals"][-1][1]

        candidates = [
            a for a in alignments_by_chrom.get(chrom, [])
            if a["strand"] == strand
            and a["aln_start"] <= cds_e
            and a["aln_end"]   >= cds_s
        ]
        if not candidates:
            continue

        def _score(a):
            ext = max(0, cds_s - a["aln_start"]) + max(0, a["aln_end"] - cds_e)
            return a["identity"] * a["coverage"] + 1e-4 * ext

        best   = max(candidates, key=_score)
        result = _utr_blocks_for_alignment(mrna, best, min_utr_len, max_utr_len)
        result["supporting_tx"] = best["tx_id"]
        result["identity"]      = best["identity"]
        result["coverage"]      = best["coverage"]
        utrs[mid] = result

        if result["five_prime_utr"]:
            n_with_5 += 1
        if result["three_prime_utr"]:
            n_with_3 += 1

    if n_skipped:
        _log(f"  {n_skipped} mRNAs skipped — already carry UTR features")
    _log(f"  mRNAs with 5' UTR added  : {n_with_5}")
    _log(f"  mRNAs with 3' UTR added  : {n_with_3}")
    return utrs


# ── GFF3 output ───────────────────────────────────────────────────────────────

def _utr_gff3_lines(
    mrna_id: str,
    blocks: list,
    utr_type: str,
    chrom: str,
    strand: str,
    source: str = "InnsmUTR",
) -> list:
    lines = []
    for i, (s, e) in enumerate(sorted(blocks), start=1):
        attrs = f"ID={mrna_id}.{utr_type}.{i};Parent={mrna_id}"
        lines.append(
            "\t".join([chrom, source, utr_type, str(s), str(e),
                       ".", strand, ".", attrs]) + "\n"
        )
    return lines


def write_updated_gff3(
    records: list,
    mrnas: dict,
    genes: dict,
    utrs: dict,
    output_path: Path,
) -> None:
    """
    Reproduce the input GFF3 verbatim, with two modifications:
      1. gene / mRNA coordinates are expanded to encompass new UTR blocks
         when those blocks lie outside the existing feature boundaries.
      2. UTR feature lines are inserted immediately after each mRNA line.
    """
    # --- expanded coordinate ranges ---
    mrna_new_coords: dict = {}
    for mid, utr in utrs.items():
        all_blocks = utr.get("five_prime_utr", []) + utr.get("three_prime_utr", [])
        if not all_blocks:
            continue
        mrna      = mrnas[mid]
        utr_min   = min(s for s, e in all_blocks)
        utr_max   = max(e for s, e in all_blocks)
        new_start = min(mrna["start"], utr_min)
        new_end   = max(mrna["end"],   utr_max)
        if new_start != mrna["start"] or new_end != mrna["end"]:
            mrna_new_coords[mid] = (new_start, new_end)

    gene_new_coords: dict = {}
    for mid, (ns, ne) in mrna_new_coords.items():
        gid = mrnas[mid]["gene_id"]
        if gid and gid in genes:
            prev = gene_new_coords.get(gid, (genes[gid]["start"], genes[gid]["end"]))
            gene_new_coords[gid] = (min(prev[0], ns), max(prev[1], ne))

    # --- build UTR line lists ---
    utr_lines_by_mrna: dict = {}
    for mid, utr in utrs.items():
        mrna  = mrnas[mid]
        lines = []
        for utype in ("five_prime_utr", "three_prime_utr"):
            blocks = utr.get(utype, [])
            if blocks:
                lines.extend(
                    _utr_gff3_lines(
                        mid, blocks, utype, mrna["chrom"], mrna["strand"]
                    )
                )
        if lines:
            utr_lines_by_mrna[mid] = lines

    # --- stream-copy records, patching as we go ---
    with open(output_path, "w") as fh:
        for line in records:
            if line.startswith("#") or not line.strip():
                fh.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                fh.write(line)
                continue

            ftype = cols[2]
            attrs = _parse_attrs(cols[8])
            fid   = attrs.get("ID", "")

            if ftype == "gene" and fid in gene_new_coords:
                cols[3] = str(gene_new_coords[fid][0])
                cols[4] = str(gene_new_coords[fid][1])
                fh.write("\t".join(cols) + "\n")
                continue

            if ftype in ("mRNA", "transcript"):
                if fid in mrna_new_coords:
                    cols[3] = str(mrna_new_coords[fid][0])
                    cols[4] = str(mrna_new_coords[fid][1])
                fh.write("\t".join(cols) + "\n")
                for utr_line in utr_lines_by_mrna.get(fid, []):
                    fh.write(utr_line)
                continue

            fh.write(line)

    _log(f"  Updated GFF3 → {output_path.name}")


# ── Stats ─────────────────────────────────────────────────────────────────────

def write_utr_stats(mrnas: dict, utrs: dict, path: Path) -> None:
    headers = [
        "mrna_id", "gene_id", "chrom", "strand",
        "five_utr_bp", "five_utr_exons",
        "three_utr_bp", "three_utr_exons",
        "supporting_tx", "aln_identity", "aln_coverage",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(headers)
        for mid in sorted(mrnas):
            mrna  = mrnas[mid]
            utr   = utrs.get(mid, {})
            five  = utr.get("five_prime_utr",  [])
            three = utr.get("three_prime_utr", [])
            w.writerow([
                mid,
                mrna["gene_id"],
                mrna["chrom"],
                mrna["strand"],
                sum(e - s + 1 for s, e in five),
                len(five),
                sum(e - s + 1 for s, e in three),
                len(three),
                utr.get("supporting_tx", ""),
                utr.get("identity",      ""),
                utr.get("coverage",      ""),
            ])
    _log(f"  UTR stats → {path.name}")


# ── Run summary ───────────────────────────────────────────────────────────────

def write_run_summary(
    args,
    mrnas: dict,
    utrs: dict,
    path: Path,
    elapsed_s: float,
    peak_mem_mb: float,
    emissions_kg,
) -> None:
    n_with_5 = sum(1 for u in utrs.values() if u.get("five_prime_utr"))
    n_with_3 = sum(1 for u in utrs.values() if u.get("three_prime_utr"))
    summary  = {
        "date":              datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "version":           VERSION,
        "input_gff3":        str(args.gff3),
        "input_genome":      str(args.genome),
        "input_transcripts": str(args.transcripts),
        "n_mrnas_total":     len(mrnas),
        "n_mrnas_with_5utr": n_with_5,
        "n_mrnas_with_3utr": n_with_3,
        "parameters": {
            "aligner":      args.aligner,
            "threads":      args.threads,
            "min_utr_len":  args.min_utr_len,
            "max_utr_len":  args.max_utr_len,
            "min_coverage": args.min_coverage,
            "min_identity": args.min_identity,
        },
        "resource_usage": {
            "wall_clock_s":       round(elapsed_s, 1),
            "peak_mem_mb":        round(peak_mem_mb, 1),
            "emissions_kg_CO2eq": emissions_kg,
        },
    }
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    _log(f"  Run summary → {path.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="InnsmUTR",
        description="Add UTRs to a genome annotation using transcript evidence.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--gff3",        required=True, type=Path,
                    help="Input genome annotation (GFF3)")
    ap.add_argument("--genome",      required=True, type=Path,
                    help="Reference genome FASTA")
    ap.add_argument("--transcripts", required=True, type=Path,
                    help="Transcript FASTA (cDNA / IsoSeq / Trinity / StringTie)")
    ap.add_argument("--output",      required=True,
                    help="Output directory name / run prefix")
    ap.add_argument("--aligner",     default="minimap2",
                    choices=["minimap2", "gmap"],
                    help="Splice-aware aligner  (default: minimap2)")
    ap.add_argument("--threads",     type=int, default=4,
                    help="CPU threads  (default: 4)")
    ap.add_argument("--min_utr_len", type=int, default=10,
                    help="Minimum UTR length in bp to add  (default: 10)")
    ap.add_argument("--max_utr_len", type=int, default=5000,
                    help="Maximum UTR length in bp to add  (default: 5000)")
    ap.add_argument("--min_coverage", type=float, default=0.80,
                    help="Minimum transcript alignment coverage  (default: 0.80)")
    ap.add_argument("--min_identity", type=float, default=0.95,
                    help="Minimum alignment identity  (default: 0.95)")
    ap.add_argument("--disable_co2_tracking", action="store_true",
                    help="Disable carbon footprint tracking even if codecarbon "
                         "is installed")
    ap.add_argument("--force", action="store_true",
                    help="Rerun all steps from scratch even if intermediate "
                         "outputs exist in workdir/")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    # ── Input validation ──────────────────────────────────────────────────────
    for flag, val in [
        ("--gff3",        args.gff3),
        ("--genome",      args.genome),
        ("--transcripts", args.transcripts),
    ]:
        if not val.exists():
            print(f"ERROR: {flag} not found: {val}", file=sys.stderr)
            sys.exit(1)

    # ── Directory layout ──────────────────────────────────────────────────────
    run_dir  = Path(args.output)
    results  = run_dir / "results"
    workdir  = run_dir / "workdir"
    logs_dir = run_dir / "logs"
    for d in (results, workdir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    prefix = run_dir.name

    # ── Run log ───────────────────────────────────────────────────────────────
    global _LOG_FH
    log_path = logs_dir / "Run_InnsmUTR.log"
    _LOG_FH  = open(log_path, "w")
    sep = "=" * 62
    _LOG_FH.write(f"{sep}\n  InnsmUTR {VERSION}  —  Run Log\n{sep}\n")
    _LOG_FH.write(f"Date      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    _LOG_FH.write(f"User      : {getpass.getuser()}\n")
    _LOG_FH.write(f"Server    : {platform.node()}\n")
    _LOG_FH.write(
        f"OS        : {platform.system()} {platform.release()} "
        f"({platform.machine()})\n"
    )
    _LOG_FH.write(f"Directory : {os.getcwd()}\n")
    _LOG_FH.write(f"Command   : {' '.join(sys.argv)}\n")
    _LOG_FH.write(f"{sep}\n\n")
    _LOG_FH.flush()

    # ── Carbon footprint tracker ──────────────────────────────────────────────
    _tracker = None
    if args.disable_co2_tracking:
        _log("  Carbon footprint tracking disabled (--disable_co2_tracking)")
    else:
        try:
            from codecarbon import EmissionsTracker
            _tracker = EmissionsTracker(
                output_dir=str(logs_dir),
                output_file=f"{prefix}.emissions.csv",
                project_name="InnsmUTR",
                log_level="warning",
            )
            _tracker.start()
            _log("  codecarbon tracker started")
        except ImportError:
            _log(
                "  codecarbon not installed — carbon tracking skipped "
                "(conda install -c conda-forge codecarbon)"
            )

    t_start = time.monotonic()

    if args.force:
        _log("--force set: all steps will rerun regardless of existing outputs")
    elif workdir.exists() and any(workdir.iterdir()):
        _log(
            "Existing workdir found — resuming from checkpoints "
            "(use --force to rerun all steps from scratch)"
        )

    # ── Banner ────────────────────────────────────────────────────────────────
    _banner(f"InnsmUTR  {VERSION}")
    _log(f"  GFF3        : {args.gff3}")
    _log(f"  Genome      : {args.genome}")
    _log(f"  Transcripts : {args.transcripts}")
    _log(f"  Output      : {run_dir}/")
    _log(f"  Aligner     : {args.aligner}")
    _log(f"  Threads     : {args.threads}")

    # ── Module 1: Transcript alignment ───────────────────────────────────────
    _banner("Module 1 — Transcript alignment")
    if args.aligner == "minimap2":
        sam_path = run_minimap2(
            args.genome, args.transcripts, workdir, args.threads, args.force
        )
    else:
        sam_path = run_gmap(
            args.genome, args.transcripts, workdir, args.threads, args.force
        )

    _log("  Parsing alignments ...")
    alignments_by_chrom = parse_sam_alignments(
        sam_path, args.min_identity, args.min_coverage
    )
    write_alignment_stats(
        alignments_by_chrom,
        args.min_identity,
        args.min_coverage,
        results / f"mod01_alignment_{prefix}.stats.tsv",
    )

    # ── Module 2: UTR extraction ──────────────────────────────────────────────
    _banner("Module 2 — UTR extraction")
    _log("  Parsing input GFF3 ...")
    genes, mrnas, existing_utrs, records = load_gff3(args.gff3)
    _log(f"  {len(genes)} genes, {len(mrnas)} mRNAs loaded")
    if existing_utrs:
        _log(f"  {len(existing_utrs)} mRNAs already carry UTR features (will be skipped)")

    utrs = extract_utrs(
        mrnas, existing_utrs, alignments_by_chrom,
        args.min_utr_len, args.max_utr_len,
    )

    write_updated_gff3(
        records, mrnas, genes, utrs,
        results / f"mod02_utr_{prefix}.gff3",
    )
    write_utr_stats(
        mrnas, utrs,
        results / f"mod02_utr_stats_{prefix}.tsv",
    )

    # ── Resource usage ────────────────────────────────────────────────────────
    elapsed_s   = time.monotonic() - t_start
    ru          = resource.getrusage(resource.RUSAGE_SELF)
    peak_mem_mb = (
        ru.ru_maxrss / (1024 * 1024)
        if platform.system() == "Darwin"
        else ru.ru_maxrss / 1024
    )

    emissions_kg = None
    if _tracker is not None:
        try:
            emissions_kg = _tracker.stop()
        except Exception:
            pass

    _banner("Resource usage")
    _log(f"  Wall-clock time   : {elapsed_s:.1f} s  ({elapsed_s / 60:.1f} min)")
    _log(f"  Peak memory (RSS) : {peak_mem_mb:.1f} MB")
    if emissions_kg is not None:
        _log(f"  Carbon footprint  : {emissions_kg:.6f} kg CO2eq")
        _log(f"  Emissions log     : {logs_dir}/{prefix}.emissions.csv")

    write_run_summary(
        args, mrnas, utrs,
        results / f"{prefix}.run_summary.json",
        elapsed_s, peak_mem_mb, emissions_kg,
    )

    _banner("Done")
    _log(f"  Results → {results}/")
    _log(f"  Log     → {log_path}")

    if _LOG_FH is not None:
        _LOG_FH.close()


if __name__ == "__main__":
    main()
