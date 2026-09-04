#!/usr/bin/env python3
# =============================================================================
# run_all_models.py
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Answers "Which gene features predict whether a gene accumulates mutations?"
#   for two organisms (yeast Sharp, E. coli Foster M2), then compares them.
#
#   For each organism it fits THREE complementary models so a result rests on
#   agreement between methods rather than on one:
#       1. LOGISTIC REGRESSION  -> P(gene mutated) ~ features        (primary)
#       2. NEGATIVE BINOMIAL    -> mutation_count ~ features, offset=ln(length)
#       3. KENDALL'S TAU        -> rank correlation of density vs each feature
#
#   OUTPUTS (three CSV files):
#       yeast_model_results.csv        - every model coefficient for yeast
#       ecoli_model_results.csv        - every model coefficient for E. coli
#       cross_species_comparison.csv   - shared-core logistic ORs side by side
#   It also prints a readable summary and the collinearity (VIF) diagnostics.
#
# WHY THESE CHOICES (short version)
#   ~92% of genes have zero mutations (max 3), so "hit vs not hit" is the
#   well-powered question -> logistic is primary. ln(CDS length) is always
#   included because a longer gene is a bigger target; without it other
#   predictors would just proxy for length. Protein & RNA are ~0.72 correlated,
#   so we report VIFs and also fit expression measures within the full model.
#
# PREDICTORS (the two organisms differ slightly)
#   Shared core : log_protein, log_rna, gc_content, essential, log_length
#   Yeast adds  : replication_timing   (modelled timing, Berners-Lee et al. 2025)
#   E. coli adds: strand_bin, ori_distance
#   NB: yeast replication_timing and E. coli ori_distance are DIFFERENT measures
#       and are not directly comparable across organisms.
#
# OUTPUT CSV COLUMNS (yeast/ecoli files)
#   organism, model, model_set, predictor, estimate, effect_ratio, p_value, n
#     model        = logistic | negbin | kendall_tau | VIF
#     model_set    = core | full
#     estimate     = coef (logistic/negbin), tau (kendall), or VIF value
#     effect_ratio = odds_ratio (logistic) / rate_ratio (negbin); blank otherwise
#
# USAGE
#   python3 run_all_models.py
#   python3 run_all_models.py --yeast sharp_gene_table.csv --ecoli M2_gene_table.csv --outdir .
#
# REQUIRES: pandas, numpy, scipy, statsmodels
#   install with:  python3 -m pip install pandas numpy scipy statsmodels
# =============================================================================

import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import kendalltau
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor


# -----------------------------------------------------------------------------
# 1. LOAD + PREPARE ONE ORGANISM'S GENE TABLE
# -----------------------------------------------------------------------------
def prepare(path, is_yeast):
    """
    Load a gene table and derive model-ready columns.
      - Yeast: drop non-nuclear genes ('chrMito', '2micron'): different
        mutational process + no nuclear-database predictor coverage.
      - Log10-transform protein, RNA, length (span many orders of magnitude).
      - Essentiality -> numeric 0/1 (NA rows dropped later).
      - E. coli: strand -> 0/1 (1 = minus strand; leading/lagging proxy).
    """
    df = pd.read_csv(path)
    if is_yeast and "chromosome" in df.columns:
        df = df[~df["chromosome"].isin(["chrMito", "2micron"])].copy()
    df["log_protein"] = np.log10(df["protein_abundance"].where(df["protein_abundance"] > 0))
    df["log_rna"]     = np.log10(df["rna_level"].where(df["rna_level"] > 0))
    df["log_length"]  = np.log10(df["CDS_length_bp"])
    df["essential"]   = pd.to_numeric(df["essential"], errors="coerce")
    if not is_yeast and "strand" in df.columns:
        df["strand_bin"] = (df["strand"] == "-").astype(float)
    return df


def complete_cases(df, predictors):
    """Keep only rows with every predictor + response present."""
    need = predictors + ["mutated", "mutation_count", "mutation_density"]
    return df.dropna(subset=[c for c in need if c in df.columns]).copy()


# -----------------------------------------------------------------------------
# 2. COLLINEARITY (VIF) -- printed, and returned as rows
# -----------------------------------------------------------------------------
def vif_rows(df, predictors, organism, model_set):
    X = sm.add_constant(df[predictors].astype(float))
    print(f"\n[Collinearity / VIF] {organism} ({model_set})")
    rows = []
    for i, name in enumerate(X.columns):
        if name == "const":
            continue
        v = variance_inflation_factor(X.values, i)
        print(f"    {name:18s} VIF = {v:5.2f}{'   <-- high' if v > 5 else ''}")
        rows.append(dict(organism=organism, model="VIF", model_set=model_set,
                         predictor=name, estimate=round(v, 3),
                         effect_ratio=np.nan, p_value=np.nan, n=int(len(df))))
    return rows


# -----------------------------------------------------------------------------
# 3. LOGISTIC REGRESSION -> result rows + fitted model
# -----------------------------------------------------------------------------
def logistic_rows(df, predictors, organism, model_set):
    model = smf.logit("mutated ~ " + " + ".join(predictors), data=df).fit(disp=False)
    n = int(model.nobs)
    print(f"\n[Logistic] {organism} ({model_set})  n={n}, mutated={int(df['mutated'].sum())}")
    print(f"    {'predictor':18s}{'coef':>9}{'odds_ratio':>12}{'p_value':>10}")
    rows = []
    for name in model.params.index:
        if name == "Intercept":
            continue
        c, p = model.params[name], model.pvalues[name]
        print(f"    {name:18s}{c:9.3f}{np.exp(c):12.3f}{p:10.3f} {'*' if p < 0.05 else ' '}")
        rows.append(dict(organism=organism, model="logistic", model_set=model_set,
                         predictor=name, estimate=round(c, 4),
                         effect_ratio=round(np.exp(c), 4), p_value=round(p, 4), n=n))
    return rows, model


# -----------------------------------------------------------------------------
# 4. NEGATIVE BINOMIAL -> result rows
# -----------------------------------------------------------------------------
def negbin_rows(df, predictors, organism, model_set):
    try:
        model = smf.glm(
            "mutation_count ~ " + " + ".join(predictors),
            data=df, family=sm.families.NegativeBinomial(),
            offset=np.log(df["CDS_length_bp"].values),
        ).fit()
    except Exception as e:
        print(f"\n[NegBin] {organism} ({model_set}): did not converge ({e})")
        return []
    n = int(model.nobs)
    print(f"\n[NegBin] {organism} ({model_set})  n={n}")
    print(f"    {'predictor':18s}{'coef':>9}{'rate_ratio':>12}{'p_value':>10}")
    rows = []
    for name in model.params.index:
        if name == "Intercept":
            continue
        c, p = model.params[name], model.pvalues[name]
        print(f"    {name:18s}{c:9.3f}{np.exp(c):12.3f}{p:10.3f} {'*' if p < 0.05 else ' '}")
        rows.append(dict(organism=organism, model="negbin", model_set=model_set,
                         predictor=name, estimate=round(c, 4),
                         effect_ratio=round(np.exp(c), 4), p_value=round(p, 4), n=n))
    return rows


# -----------------------------------------------------------------------------
# 5. KENDALL'S TAU -> result rows
# -----------------------------------------------------------------------------
def kendall_rows(df, predictors, organism, model_set):
    print(f"\n[Kendall tau vs mutation_density] {organism} ({model_set})")
    print(f"    {'predictor':18s}{'tau':>9}{'p_value':>10}")
    rows = []
    for p_ in predictors:
        tau, pval = kendalltau(df[p_], df["mutation_density"])
        print(f"    {p_:18s}{tau:9.3f}{pval:10.3f} {'*' if pval < 0.05 else ' '}")
        rows.append(dict(organism=organism, model="kendall_tau", model_set=model_set,
                         predictor=p_, estimate=round(tau, 4),
                         effect_ratio=np.nan, p_value=round(pval, 4), n=int(len(df))))
    return rows


# -----------------------------------------------------------------------------
# 6. RUN EVERYTHING FOR ONE ORGANISM -> all rows + core logistic model
# -----------------------------------------------------------------------------
def analyse(path, is_yeast, organism):
    print("\n" + "=" * 78)
    print(f"ORGANISM: {organism}    (file: {path})")
    print("=" * 78)

    df = prepare(path, is_yeast)
    core = ["log_protein", "log_rna", "gc_content", "essential", "log_length"]
    if is_yeast:
        full = core + (["replication_timing"] if "replication_timing" in df.columns else [])
    else:
        full = core + [c for c in ["strand_bin", "ori_distance"] if c in df.columns]

    dcore = complete_cases(df, core)
    dfull = complete_cases(df, full)
    print(f"\nComplete-case genes: core={len(dcore)} (mutated={int(dcore['mutated'].sum())}) | "
          f"full={len(dfull)} (mutated={int(dfull['mutated'].sum())})")

    rows = []
    rows += vif_rows(dcore, core, organism, "core")

    print("\n----- SHARED-CORE MODELS -----")
    core_rows, core_model = logistic_rows(dcore, core, organism, "core")
    rows += core_rows
    rows += negbin_rows(dcore, core, organism, "core")
    rows += kendall_rows(dcore, core, organism, "core")

    if full != core:
        print("\n----- FULL MODEL (organism-specific predictors incl. timing) -----")
        full_rows, _ = logistic_rows(dfull, full, organism, "full")
        rows += full_rows
        rows += negbin_rows(dfull, full, organism, "full")
        rows += kendall_rows(dfull, full, organism, "full")

    return rows, core_model


# -----------------------------------------------------------------------------
# 7. MAIN: run both, write three CSVs
# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--yeast", default="sharp_gene_table.csv")
    ap.add_argument("--ecoli", default="M2_gene_table.csv")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    col_order = ["organism", "model", "model_set", "predictor",
                 "estimate", "effect_ratio", "p_value", "n"]

    # --- YEAST ---
    yeast_rows, yeast_model = analyse(args.yeast, True, "Yeast (Sharp)")
    yeast_df = pd.DataFrame(yeast_rows)[col_order]
    ypath = os.path.join(args.outdir, "yeast_model_results.csv")
    yeast_df.to_csv(ypath, index=False, lineterminator="\n")

    # --- E. COLI ---
    ecoli_rows, ecoli_model = analyse(args.ecoli, False, "E. coli (Foster M2)")
    ecoli_df = pd.DataFrame(ecoli_rows)[col_order]
    epath = os.path.join(args.outdir, "ecoli_model_results.csv")
    ecoli_df.to_csv(epath, index=False, lineterminator="\n")

    # --- CROSS-SPECIES COMPARISON (shared-core logistic odds ratios) ---
    print("\n" + "=" * 78)
    print("CROSS-SPECIES COMPARISON  (shared-core logistic odds ratios)")
    print("=" * 78)
    core = ["log_protein", "log_rna", "gc_content", "essential", "log_length"]
    comp = []
    print(f"    {'predictor':18s}{'yeast_OR':>10}{'yeast_p':>9}"
          f"{'ecoli_OR':>10}{'ecoli_p':>9}   direction")
    for name in core:
        yb, yp = yeast_model.params.get(name, np.nan), yeast_model.pvalues.get(name, np.nan)
        eb, ep = ecoli_model.params.get(name, np.nan), ecoli_model.pvalues.get(name, np.nan)
        yOR, eOR = np.exp(yb), np.exp(eb)
        direction = "same" if np.sign(yb) == np.sign(eb) else "OPPOSITE"
        both_sig = bool((yp < 0.05) and (ep < 0.05))
        print(f"    {name:18s}{yOR:10.3f}{yp:9.3f}{eOR:10.3f}{ep:9.3f}   {direction}")
        comp.append(dict(
            predictor=name,
            yeast_odds_ratio=round(yOR, 4), yeast_p=round(yp, 4),
            ecoli_odds_ratio=round(eOR, 4), ecoli_p=round(ep, 4),
            direction=direction, significant_in_both=both_sig,
        ))
    comp_df = pd.DataFrame(comp)
    cpath = os.path.join(args.outdir, "cross_species_comparison.csv")
    comp_df.to_csv(cpath, index=False, lineterminator="\n")

    print("\nWrote three CSV files:")
    print(f"   {ypath}")
    print(f"   {epath}")
    print(f"   {cpath}")
    print("\nNote: yeast replication_timing and E. coli ori_distance are different")
    print("measures (in each organism's FULL model rows, not the cross-species")
    print("table) and should not be treated as the same quantity.")


if __name__ == "__main__":
    main()
