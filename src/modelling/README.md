# Gene-level modelling pipeline

`run_all_models.py` — builds nothing itself; it reads the per-gene feature tables
and fits three complementary models per organism (logistic [primary], negative
binomial with ln(length) offset, Kendall's tau), reports VIFs, and writes a
cross-species comparison.

## Inputs (per-gene feature tables, in `data/processed/`)
- `--yeast` default `sharp_gene_table.csv`
- `--ecoli` default `M2_gene_table.csv`

## Run
```bash
python src/modelling/run_all_models.py \
    --yeast data/processed/sharp_gene_table.csv \
    --ecoli data/processed/M2_gene_table.csv \
    --outdir results
```

## Outputs (to `results/`)
- `yeast_model_results.csv`
- `ecoli_model_results.csv`
- `cross_species_comparison.csv`

Requires: pandas, numpy, scipy, statsmodels.
