# InnsmUTR

<p align="center">
  <img src="assets/innsmutr_logo.svg" width="260" alt="InnsmUTR logo"/>
</p>

Benchmarking tool for UTR annotation. Compares a UTR-annotated GFF3 (e.g. the
output of a PASA UTR-addition run) against a reference / gold-standard GFF3,
and reports presence and boundary-accuracy metrics for 5' and 3' UTRs.

## Why

Adding UTRs to a gene annotation (e.g. with PASA, using full-length
transcript evidence) is easy to run but hard to evaluate — there's no
single "correct" boundary, and call quality varies between under-extension
(missed UTR), over-extension (spurious UTR), and boundary drift. InnsmUTR
gives a repeatable way to score a UTR-annotation run against a reference.

## How it works

For each transcript, the 5' and 3' UTR length is read directly from the
explicit `five_prime_UTR` / `three_prime_UTR` feature lines in the GFF3,
summed per transcript so spliced (multi-exon) UTRs are counted correctly.
**This script does not derive UTRs from `mRNA`/`CDS`/`exon` spans** — if
your GFF3 only has `CDS` and `exon` records (no UTR features), run it
through [AGAT](https://github.com/NBISweden/AGAT) first to generate them,
then compare the AGAT output:

```bash
conda install -c bioconda agat
agat_convert_sp_gxf2gxf.pl -g your_annotation.gff3 -o your_annotation.utrs.gff3
```

(Check AGAT's documentation for the UTR-generation behavior of your
installed AGAT version — recent versions infer and add UTR features as
part of normalizing the gene model hierarchy from `exon`/`CDS` records.)

Transcripts are matched between target and reference either by:
- **shared transcript ID** (default, `--match_by id`) — use when target and
  reference describe the same gene set (e.g. before/after running PASA on
  one annotation)
- **CDS-span overlap** (`--match_by overlap`) — use when target and
  reference come from independent annotations with unrelated IDs; matches
  each reference transcript to the target transcript with the largest
  CDS-span overlap on the same chrom/strand (requires `bedtools`/`pybedtools`)

UTR presence per side is scored as a binary call against the reference
(TP/FP/FN/TN), and boundary accuracy is scored among presence-agreeing
transcripts as the UTR length difference (target − reference), which is
already strand-corrected — a positive offset means the target UTR is
longer (extends further) than the reference.

## Usage

```bash
# Match by shared transcript ID
python scripts/compare_utr_annotation.py \
    -t pasa_updated.gff3 -r reference.gff3 \
    --outdir comparisons/ --prefix arath --tolerance 20

# Match by CDS overlap (IDs differ between target and reference)
python scripts/compare_utr_annotation.py \
    -t pasa_updated.gff3 -r reference.gff3 \
    --match_by overlap --min_cds_overlap 0.9 \
    --outdir comparisons/ --prefix arath
```

### Options

| Flag | Default | Description |
|---|---|---|
| `-t, --target` | required | Target GFF3 (e.g. PASA-updated annotation) |
| `-r, --reference` | required | Reference / gold-standard GFF3 |
| `--outdir` | required | Output directory |
| `--prefix` | required | Prefix for output file names |
| `--match_by` | `id` | `id` or `overlap` |
| `--id_attr` | `ID` | GFF3 attribute used as transcript ID |
| `--mrna_feature` | `mRNA` | GFF3 feature type for transcripts |
| `--min_cds_overlap` | `0.9` | Min fraction of reference CDS span covered, for overlap matching |
| `--tolerance` | `20` | UTR boundary tolerance in bp |

## Output

- `{prefix}_comparison.tsv` — one row per matched transcript: UTR length
  (ref/target), presence call (TP/FP/FN/TN), boundary offset in bp, and
  whether the offset is within tolerance, for each of 5' and 3' UTR.
- `{prefix}_metrics.tsv` — per side (5'/3'): TP/FP/FN/TN counts, precision,
  recall, F1 (presence-based), boundary accuracy rate (fraction of
  presence-agreeing transcripts within tolerance), and mean/median
  boundary offset in bp.

## Requirements

```bash
conda env create -f envs/innsmutr.yaml
conda activate innsmutr
```

`bedtools`/`pybedtools` are only required for `--match_by overlap`; `--match_by id`
(the default) has no dependency beyond `pandas`.
