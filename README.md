# Intragenomic mutation rate variation is shaped by genome architecture

Code and processed data for the analysis of within-genome variation in
spontaneous mutation rates in *Saccharomyces cerevisiae* and *Escherichia coli*
K-12, using selection-screened mutation accumulation (MA) data.

This repository accompanies the paper *"Selection-Screened Mutation Accumulation
Data Shows Intragenomic Mutation Rate Variation is Shaped by Genome Architecture,
not Gene Function"*.

## What this does

Two pipelines:

1. **Selection screening (Ka/Ks)** — screens each MA dataset for residual
   purifying selection using a trinucleotide-context-matched permutation test,
   so that downstream analysis treats retained mutations as an approximately
   neutral sample. See `src/kaks/`.
2. **Gene-level modelling** — builds per-gene feature tables and models mutation
   incidence against functional predictors (protein abundance, mRNA level,
   essentiality) and architectural predictors (GC content, replication timing,
   chromosomal position), controlling for gene length. See `src/modelling/`.

## Repository layout

```
data/
  raw/         # third-party inputs — NOT redistributed here; see data/DATA_SOURCES.md
  processed/   # processed, model-ready data used by this pipeline
src/
  kaks/        # selection-screening (Ka/Ks permutation) pipeline
  modelling/   # feature-table construction and statistical modelling
results/        # model output CSVs
```

## Data availability

Processed, model-ready data used by this pipeline (the cleaned mutation tables
and per-gene feature tables) are included in `data/processed/`, model outputs in
`results/`, and the reference CDS FASTAs in `data/raw/`. Together these let the
two pipelines be run end-to-end.

Note that these processed files are provided as inputs; the scripts that build
them from the original supplements and predictor databases are not included, so
the repository reproduces the analysis from the cleaned data rather than
regenerating that data from scratch. `data/DATA_SOURCES.md` lists every input,
its source, and version/accession.

## Reproducing the analysis

Requires Python 3.11+ and PAML (for `yn00`/`codeml`) available on your `PATH`.

```bash
# 1. set up environment
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. fetch raw inputs (see data/DATA_SOURCES.md for what/where)
#    place them under data/raw/ as described there

# 3. selection screening (per dataset) — see src/kaks/README.md
python src/kaks/kaks_selection.py \
    --csv <mutations.csv> --fasta <reference_cds.fasta> --out result

# 4. gene-level modelling — see src/modelling/README.md
python src/modelling/run_all_models.py \
    --yeast data/processed/sharp_gene_table.csv \
    --ecoli data/processed/M2_gene_table.csv \
    --outdir results
```

PAML (`yn00`/`codeml`) is only required for the Ka/Ks validation options
(`--engine paml` / `--codeml`); the default `fast` NG86 engine needs no PAML.
