# Data sources

This file records each input used in the analysis, its origin, and the
version/accession used.

## What is included in this repository

Included directly:
- `data/processed/sharp_clean.csv`, `data/processed/M2_clean.csv` — cleaned,
  gene-oriented mutation tables (Ka/Ks inputs).
- `data/processed/sharp_gene_table.csv`, `data/processed/M2_gene_table.csv` —
  per-gene feature tables (modelling inputs).
- `data/raw/orf_coding_all_R64-2-1_20150113.fasta` (SGD R64-2-1 ORF CDS) and
  `data/raw/NC_000913_2_cds_nucl.fasta` (NCBI RefSeq NC_000913.2 CDS) —
  reference sequences, committed so the pipeline runs out of the box.

Not included (fetch from the sources below): the original Sharp/Foster mutation
supplements and the predictor databases (PaxDB, PRECISE-1K, PEC, SGD phenotypes,
Berners-Lee timing) used to build the feature tables.

## Mutation accumulation datasets

| Dataset | Organism | Source | Accession / location |
|---|---|---|---|
| Sharp et al. 2018 (wild-type MA) | *S. cerevisiae* | PNAS 115(22):E5046–E5055, doi:10.1073/pnas.1801040115 | Supplementary data |
| Foster et al. 2015 (M2, K-12 MG1655) | *E. coli* | PNAS 112(44):E5990–E5999, doi:10.1073/pnas.1512136112 | Supplementary data |

Datasets evaluated but excluded (see paper §2.3): Zhu et al. 2014
(doi:10.1073/pnas.1323011111); Liu & Zhang 2019 (doi:10.1016/j.cub.2019.03.054);
Sui et al. 2020 (doi:10.1073/pnas.2018633117); Lee et al. 2012
(doi:10.1073/pnas.1210309109); Foster lab ED1a / IAI1 lines; Sane et al. 2025
(doi:10.1371/journal.pbio.3003282).

## Reference sequences / annotations

| Input | Organism | Source | Version |
|---|---|---|---|
| ORF coding set | *S. cerevisiae* | SGD | R64-2-1 |
| Genome / CDS | *E. coli* K-12 MG1655 | NCBI (via EDirect) | NC_000913.2 |

## Predictor databases

| Predictor | Organism | Source | Version / accession |
|---|---|---|---|
| Protein abundance | both | PaxDB (Huang et al. 2025) | v6.0 |
| mRNA abundance | *E. coli* | PRECISE-1K (Lamoureux et al. 2023) | MG1655, glucose, WT baseline |
| mRNA abundance | *S. cerevisiae* | Pelechano et al. 2010 (orig. Nagalakshmi et al. 2008; Miura et al. 2006) | — |
| Essentiality | *E. coli* | PEC (Yamazaki, Niki & Kato 2008) | based on NC_000913.2 |
| Essentiality | *S. cerevisiae* | SGD phenotype annotations (Engel et al. 2024) | systematic deletion collection |
| Replication timing | *S. cerevisiae* | Berners-Lee et al. 2025 model | 1 kb resolution |
| Replication timing | *E. coli* | Derived: normalised distance from oriC | computed from NC_000913.2 |
