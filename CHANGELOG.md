# Changelog

All notable changes to InnsmUTR are documented here.
Dates follow ISO 8601 (YYYY-MM-DD). Changes are grouped by version and type.

---

## [v0.2.0] — 2026-06-24

### Added

- **`InnsmUTR.py`** — new script that *adds* UTRs to a genome annotation
  using full-length transcript evidence (cDNA, IsoSeq, Trinity/StringTie).
  Complements `compare_utr_annotation.py`, which *benchmarks* UTR calls
  against a reference.

  Core workflow:
  1. Align transcripts to the genome with minimap2 (`-ax splice`) or GMAP
     (`--format=samse`); only primary alignments are used.
  2. Parse SAM CIGAR strings to build per-transcript exon blocks.
  3. Match each mRNA's CDS extent against overlapping aligned transcripts;
     pick the best-scoring alignment (identity × coverage, with a bonus for
     extending beyond the CDS).
  4. Clip transcript exon blocks at the CDS boundaries to derive 5' and 3'
     UTR intervals; apply `--min_utr_len` / `--max_utr_len` filters.
  5. Reproduce the input GFF3, inserting `five_prime_utr` / `three_prime_utr`
     features after each mRNA line and updating mRNA/gene coordinates where
     UTRs extend beyond the existing boundaries.
  6. Write a per-mRNA stats TSV with UTR lengths, exon counts, and
     alignment quality metrics.
  7. Write a run summary JSON with parameter and resource usage details.

- **Structured output directory** — `{output}/results/`, `{output}/workdir/`,
  `{output}/logs/`; intermediate SAM files stay in `workdir/`.
- **Module-prefixed output filenames** — `mod01_alignment_*` (alignment
  stats), `mod02_utr_*` (updated GFF3 and stats table).
- **Run log** (`logs/Run_InnsmUTR.log`) — records date, user, server
  hostname, OS, working directory, and full command on every run.
- **Carbon footprint tracking** — automatic when `codecarbon` is installed;
  writes energy and CO₂eq to `logs/{prefix}.emissions.csv`.
- **`--disable_co2_tracking`** — opt out of tracking without uninstalling
  codecarbon.
- **Checkpoint / resume logic** — re-running the same command skips the
  alignment step if `workdir/transcripts.sam` already exists.
- **`--force`** — bypasses all checkpoints and reruns every step from scratch.
- **Resource usage in `run_summary.json`** — wall-clock time (s), peak RSS
  memory (MB), CO₂eq emissions.
- **`--aligner minimap2|gmap`** — choose between minimap2 (default) and GMAP
  as the splice-aware aligner.

---

## [v0.1.0] — 2026-06-20

### Added

- **Initial release** of InnsmUTR.
- **`compare_utr_annotation.py`** — benchmarking tool that compares a
  UTR-annotated GFF3 against a reference GFF3 and reports presence/absence
  metrics (TP/FP/FN/TN, precision, recall, F1) and boundary accuracy
  (fraction of transcripts within a configurable tolerance in bp) for both
  5' and 3' UTRs.
- **Match modes** — `--match_by id` (shared transcript ID, default) and
  `--match_by overlap` (largest CDS-span overlap, requires bedtools /
  pybedtools) for comparing annotations with unrelated IDs.
- **Output tables** — `{prefix}_comparison.tsv` (per-transcript detail) and
  `{prefix}_metrics.tsv` (per-side aggregated metrics).
- Conda environment (`envs/innsmutr.yaml`) with pandas, bedtools, and
  pybedtools.
