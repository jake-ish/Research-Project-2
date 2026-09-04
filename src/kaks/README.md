# Selection-screening (Ka/Ks) pipeline

`kaks_selection.py` — screens an MA dataset for purifying selection using a
trinucleotide-context-matched permutation test. Pure standard library (no
Biopython); optionally calls PAML (`yn00`/`codeml`) for validation.

Two tests are run:
- **Test 1 (Ka/Ks):** observed ω from a concatenated CDS "supergene" vs a null
  built by re-allocating each mutation to a context-matched site (1,000 perms).
- **Test 2 (stops):** observed stop-codon count vs the same context-matched null.

## Inputs
- `--csv`   cleaned, gene-oriented mutation table (in `data/processed/`)
- `--fasta` gene-oriented reference CDS FASTA (in `data/raw/`)

## Run
```bash
# Yeast (Sharp)
python src/kaks/kaks_selection.py \
    --csv data/processed/sharp_clean.csv \
    --fasta data/raw/orf_coding_all_R64-2-1_20150113.fasta \
    --out results/kaks_sharp

# E. coli (Foster M2)
python src/kaks/kaks_selection.py \
    --csv data/processed/M2_clean.csv \
    --fasta data/raw/NC_000913_2_cds_nucl.fasta \
    --out results/kaks_M2

# Validate the fast NG86 engine against real PAML yn00 on N replicates
python src/kaks/kaks_selection.py \
    --csv data/processed/sharp_clean.csv \
    --fasta data/raw/orf_coding_all_R64-2-1_20150113.fasta \
    --engine paml --paml-validate 50
```
Key options: `--nperm` (default 1000), `--seed` (default 1), `--engine fast|paml`,
`--codeml`. PAML is only needed for `--engine paml` / `--codeml`.
