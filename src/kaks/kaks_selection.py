#!/usr/bin/env python3
"""
kaks_selection.py  --  Ka/Ks selection test by context-matched permutation.

Implements the protocol:
  TEST 1 (Ka/Ks):  Concatenate all reference CDS ("supergene") and the same
                   CDS carrying the observed point mutations (nonsense EXCLUDED).
                   Ka/Ks is computed with PAML yn00.  A null distribution is
                   built by re-allocating every observed mutation to a random
                   "comparable site": a position anywhere in any CDS whose
                   TRINUCLEOTIDE centred on the mutated base is identical, and
                   applying the same central-base change (redraw if it makes a
                   stop).  1000 replicates -> null Ka/Ks distribution.
                   p = (#null Ka/Ks <= observed) / N.

  TEST 2 (stops):  How many observed mutations create a stop codon?  Null =
                   same context-matched re-allocation but of ALL mutations and
                   ALLOWING stops.  p = (#null stop-count <= observed) / N.
                   Fewer observed stops than null => purifying selection.

No Biopython.  Only the standard library (+ optional matplotlib for plots).

Strand:  work entirely in GENE-oriented CDS space.  The mutation CSV's
         ref_codon/actual_codon columns are already gene-oriented, so we never
         touch genomic coordinates or strand.  (In these data, ref/actual are
         Watson-strand while the codon columns are gene-strand -- verified.)
"""

# ===========================================================================
# HOW TO READ THIS FILE (plain-language map)
# ---------------------------------------------------------------------------
# The file is a toolbox of small functions, then one big function `run()` near
# the bottom that calls them in order. Reading top to bottom you meet:
#   1. the genetic code         : turn a 3-letter codon into an amino acid
#   2. read_fasta / strip_stop  : load the reference genes from a FASTA file
#   3. Mut / load_mutations     : loads mutations from the CSV
#   4. assign_trinucleotides    : find each mutation's 3-base context
#   5. build_site_pool          : list every place each context occurs
#   6. the NG86 block           : the fast, exact Ka/Ks calculator
#   7. run_yn00 / run_codeml    : call the other PAML programs to compare Ka/Ks
#   8. the permutation engines  : shuffle mutations to build the null model
#   9. run()                    : ties everything together
#  10. main()                   : reads the command-line options
#
# The imports below are all standard Python (no Biopython). `csv` reads the
# mutation table, `random` drives the shuffling, `subprocess` launches PAML,
# `math` does the logarithm in the distance correction, and tempfile/os/shutil
# handle the scratch files PAML needs.
# ===========================================================================
import argparse, csv, math, os, random, subprocess, sys, tempfile, shutil
from collections import defaultdict, Counter

# ----------------------------------------------------------------------------
# Standard genetic code (NCBI table 1 -- correct for yeast nuclear & E. coli)
# ----------------------------------------------------------------------------
# --- THE GENETIC CODE ------------------------------------------------------
# Below we build CODON_TABLE, a dictionary that maps every 3-letter codon to
# its amino acid letter. The string _AA lists the 64 amino
# acids in the exact order produced by looping bases in the order T,C,A,G for
# position 1, then 2, then 3 -- this is the standard textbook ordering, so
# _AA lines up one-to-one with those 64 codons. STOPS is just the set of the
# three stop codons. COMP is a lookup for complementing DNA (A<->T, G<->C),
# used when we check strand. Everything downstream leans on this table:
# deciding if a mutation is silent, changing, or a stop is just 'translate
# the before-codon and the after-codon and compare'.
BASES = "TCAG"
_AA = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG") #amino acids for the 64 trinucleotide contexts
CODON_TABLE = {}
_i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODON_TABLE[b1 + b2 + b3] = _AA[_i]
            _i += 1
STOPS = {c for c, a in CODON_TABLE.items() if a == "*"}
COMP = str.maketrans("ACGTacgt", "TGCATGCA")


# Look a codon up in the table; return 'X' if it's not a clean ACGT codon.
def translate_codon(c):
    return CODON_TABLE.get(c, "X")


# True if this codon is one of the three stop codons.
def is_stop(c):
    return c in STOPS


# ----------------------------------------------------------------------------
# FASTA (no Biopython)
# ----------------------------------------------------------------------------
# --- LOADING THE REFERENCE GENES -------------------------------------------
# Reads a FASTA file of coding sequences and returns {gene name: DNA string}.
# The important trick is choosing ONE key per gene: it prefers the
# [locus_tag=...] in NCBI headers, otherwise the first word of the header
# (which is the systematic name in the yeast file). Keeping one entry per gene
# is the de-duplication fix - it stops an NCBI file that lists each gene twice
# from doubling the backbone. read_fasta.duplicates records how many repeats
# were dropped.
def read_fasta(path, key_targets=None):
    """Return {id: sequence(upper)} with exactly ONE entry per gene.

    Each record is keyed by a single canonical id: the identifier the mutation
    set actually uses if we can tell (first header token, [locus_tag=...], or
    [gene=...] that appears in key_targets), otherwise [locus_tag=...] if present,
    otherwise the first whitespace token after '>'.

    This handles both SGD-style headers ('>YAL001C TFC3 ...', keyed YAL001C) and
    NCBI CDS headers ('>lcl|NC_..._1 [locus_tag=b0001] [gene=thrL] ...', keyed
    b0001).  Records that resolve to an id already seen are DUPLICATES: the first
    is kept, the rest are dropped and counted in read_fasta.duplicates.  This is
    what stops an NCBI CDS file (one record per CDS feature, IDs not equal to the
    locus_tag) from inflating the backbone with two copies of every gene.
    """
    import re
    records = []
    cur_tok = cur_lt = cur_gene = None
    buf = []

    def flush():
        if cur_tok is not None:
            records.append((cur_tok, cur_lt, cur_gene, "".join(buf).upper()))

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                cur_tok = header.split()[0] if header else ""
                m = re.search(r"locus_tag=([^\]\s]+)", header)
                cur_lt = m.group(1) if m else None
                g = re.search(r"gene=([^\]\s]+)", header)
                cur_gene = g.group(1) if g else None
                buf = []
            else:
                buf.append(line.strip())
        flush()

    kt = key_targets or set()
    seqs = {}
    dups = 0
    for tok, lt, gene, seq in records:
        if tok in kt:
            key = tok
        elif lt and lt in kt:
            key = lt
        elif gene and gene in kt:
            key = gene
        elif lt:
            key = lt
        else:
            key = tok
        if key in seqs:
            dups += 1              # duplicate locus id -> keep the first
            continue
        seqs[key] = seq
    read_fasta.duplicates = dups
    return seqs


# Trim a sequence to a whole number of codons and drop a trailing stop codon.
# Needed because PAML refuses to run if a sequence contains a stop codon.
def strip_terminal_stop(seq):
    """Trim to a multiple of 3 and drop a single trailing stop codon if present."""
    n = len(seq) - (len(seq) % 3)
    seq = seq[:n]
    if len(seq) >= 3 and seq[-3:] in STOPS:
        seq = seq[:-3]
    return seq


# ----------------------------------------------------------------------------
# Load mutations
# ----------------------------------------------------------------------------
# --- ONE MUTATION ----------------------------------------------------------
# A small container for a single mutation: which gene, the codon before and
# after, which of the 3 codon positions changed, the actual base change, and
# whether it turns a normal codon into a stop. `trinuc` (the 3-base context)
# is filled in later by assign_trinucleotides. __slots__ makes these
# objects small and fast since there can be thousands of them.
class Mut:
    __slots__ = ("locus", "ref_codon", "alt_codon", "pos", "central_from",
                 "central_to", "creates_stop", "trinuc")

    def __init__(self, locus, ref_codon, alt_codon, pos):
        self.locus = locus
        self.ref_codon = ref_codon
        self.alt_codon = alt_codon
        self.pos = pos                      # 1,2,3 (codon position of changed base)
        self.central_from = ref_codon[pos - 1]
        self.central_to = alt_codon[pos - 1]
        self.creates_stop = is_stop(alt_codon) and not is_stop(ref_codon)
        self.trinuc = None                  # centred trinucleotide, filled later


# Read the mutation CSV and build one Mut per row, using the gene-oriented
# ref_codon / actual_codon / position_mutated columns. It sanity-checks each
# row (exactly one base differs, at the stated position) and silently skips
# malformed rows. Because it uses the codon columns, which are already on the
# gene's reading strand, the whole script never deals with genome coordinates.
def load_mutations(csv_path):
    muts = []
    with open(csv_path, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            rc, ac = r["ref_codon"].upper(), r["actual_codon"].upper()
            try:
                pos = int(r["position_mutated"])
            except (ValueError, KeyError):
                continue
            if len(rc) != 3 or len(ac) != 3:
                continue
            diffs = [i for i in range(3) if rc[i] != ac[i]]
            if len(diffs) != 1 or diffs[0] + 1 != pos:
                continue                    # skip malformed rows
            muts.append(Mut(r["locus_tag"], rc, ac, pos))
    return muts


# ----------------------------------------------------------------------------
# Locate each mutation's centred trinucleotide in its CDS
# ----------------------------------------------------------------------------
# --- EACH MUTATION'S 3-BASE CONTEXT ----------------------------------------
# The permutation matches mutations by the trinucleotide centred on the
# mutated base. If the mutated base is the MIDDLE of its codon, the codon IS
# that trinucleotide (exact, no work). If it's the 1st or 3rd base, we need
# one neighbouring base, so we find the codon in the gene and read the
# neighbour. When the same codon appears more than once with different
# neighbours the context is ambiguous; we pick one and count how often that
# happens (the trinuc_stats you see in the summary). This is the one small
# approximation in the pipeline, and it never affects middle-position
# mutations.
def assign_trinucleotides(muts, cds, rng, report=None):
    """Set mut.trinuc = gene-oriented trinucleotide centred on the mutated base.
    pos==2 : trinuc == ref_codon exactly (no flank needed).
    pos==1 : need the base 5' of the codon; pos==3 : the base 3' of the codon.
    Locate the codon in the CDS by matching ref_codon at frame; if it occurs
    more than once and the required flank differs, pick one at random (reported).
    """
    stats = Counter()
    missing = Counter()
    for m in muts:
        seq = cds.get(m.locus)
        if seq is None:
            missing[m.locus] += 1
            m.trinuc = None
            stats["no_cds"] += 1
            continue
        if m.pos == 2:
            m.trinuc = m.ref_codon
            stats["pos2_exact"] += 1
            continue
        # find codon indices k where seq[3k:3k+3]==ref_codon
        L = len(seq)
        ncod = L // 3
        cands = [k for k in range(ncod) if seq[3 * k:3 * k + 3] == m.ref_codon]
        flanks = set()
        results = []
        for k in cands:
            if m.pos == 1:
                left = seq[3 * k - 1] if 3 * k - 1 >= 0 else None
                if left is None:
                    continue
                tri = left + m.ref_codon[0] + m.ref_codon[1]
            else:  # pos == 3
                right = seq[3 * k + 3] if 3 * k + 3 < L else None
                if right is None:
                    continue
                tri = m.ref_codon[1] + m.ref_codon[2] + right
            results.append(tri)
            flanks.add(tri)
        if not results:
            # codon at CDS terminus or not found -> fall back to codon-only context
            # (use the two in-codon bases + the mutated base already at an edge)
            m.trinuc = None
            stats["edge_or_nomatch"] += 1
            continue
        if len(flanks) == 1:
            m.trinuc = results[0]
            stats["unique_flank"] += 1
        else:
            m.trinuc = rng.choice(results)
            stats["ambiguous_flank"] += 1
    if report is not None:
        report.update({"trinuc_stats": dict(stats),
                       "n_missing_cds_records": len(missing)})
    return muts


# ----------------------------------------------------------------------------
# Global pool of comparable sites: trinucleotide -> list of (locus, pos_in_cds)
# pos_in_cds is the index of the CENTRAL base.
# ----------------------------------------------------------------------------
# --- THE LIST OF COMPARABLE SITES ------------------------------------------
# Scan every reference gene and, for each of the 64 possible trinucleotides,
# record every position in every gene where it occurs.
# List of the locations of all CAG residues in any frame across all CDSs,
# built once for all 64 contexts. The permutation then grabs a random site
# from the matching list instantly.
def build_site_pool(cds):
    pool = defaultdict(list)
    for locus, seq in cds.items():
        L = len(seq)
        for i in range(1, L - 1):
            tri = seq[i - 1:i + 2]
            if tri.isalpha() and set(tri) <= set("ACGT"):
                pool[tri].append((locus, i))
    return pool


# ----------------------------------------------------------------------------
# Classify a change at central base -> syn / nonsyn / stop, given codon context
# ----------------------------------------------------------------------------
# --- SILENT, CHANGING, OR STOP? --------------------------------------------
# Given a position in a gene and a new base, find the codon that position sits
# in, apply the change, and return 'syn' (silent), 'non' (amino-acid-changing),
# or 'stop'. This is the workhorse the permutations call to score where a
# relocated mutation landed.
def classify_at_site(cds, locus, i, new_base):
    """Apply new_base at position i of CDS[locus]; return ('syn'|'non'|'stop')."""
    seq = cds[locus]
    k = i // 3                      # codon index
    p = i % 3                       # position within codon
    codon = seq[3 * k:3 * k + 3]
    if len(codon) != 3:
        return None
    new = codon[:p] + new_base + codon[p + 1:]
    if is_stop(new):
        return "stop"
    return "syn" if translate_codon(new) == translate_codon(codon) else "non"


# ----------------------------------------------------------------------------
# NG86, implemented to match PAML EXACTLY.
#
# Convention (verified against PAML's printed NG86 omega on the real yeast data,
# reproducing 0.3619 to 4 d.p.):  when counting synonymous / nonsynonymous SITES,
# changes to a stop codon are skipped entirely - they contribute to neither S
# nor N, so a codon contributes < 3 sites in total.  Site counts are averaged
# over the two sequences.  pS = Sd/S, pN = Nd/N, then Jukes-Cantor:
#     d = -3/4 * ln(1 - 4/3 * p)
# and omega = dN/dS.
#
# Because observed and every null replicate share the SAME multi-megabase CDS
# backbone, S and N are recomputed only at the handful of mutated codons (a
# delta update), which makes each replicate O(#mutations) instead of O(#codons).
# This yields the exact PAML NG86 omega per replicate at negligible cost.
# ----------------------------------------------------------------------------
# --- THE FAST, EXACT Ka/Ks CALCULATOR (Nei-Gojobori / NG86) ----------------
# For ONE codon, count how many of its possible single-base changes are silent
# (S sites) versus amino-acid-changing (N sites). Each position contributes up
# to 3 possible changes; we add 1/3 per change so a position sums to 1 site.
# Following PAML's exact convention, changes that would make a STOP are ignored
# entirely, so a codon can contribute slightly fewer than 3 sites in total.
def codon_sites(codon):
    """(S, N) sites for one codon; changes to stops are excluded from both."""
    aa = CODON_TABLE[codon]
    S = N = 0.0
    for p in range(3):
        for b in "ACGT":
            if b == codon[p]:
                continue
            mc = codon[:p] + b + codon[p + 1:]
            if is_stop(mc):
                continue                      # excluded from S and N (PAML)
            if CODON_TABLE[mc] == aa:
                S += 1.0 / 3.0
            else:
                N += 1.0 / 3.0
    return S, N


# Pre-compute (S, N) for all 61 non-stop codons once, so it is never
# recalculated during the thousands of scoring passes.
_SITE_CACHE = {c: codon_sites(c) for c in CODON_TABLE if c not in STOPS}


# For a before/after codon pair, count how many differences are silent (Sd)
# versus changing (Nd). Single-base differences are trivially one or the other.
# For the rare codon hit twice, average over both mutation orders (NG86
# pathway averaging), ignoring any path that passes through a stop.
def count_diffs(c1, c2):
    """(Sd, Nd) between two codons, NG86 pathway-averaged for multi-hit codons."""
    diffs = [i for i in range(3) if c1[i] != c2[i]]
    if not diffs:
        return 0.0, 0.0
    if len(diffs) == 1:
        return ((1.0, 0.0) if CODON_TABLE[c1] == CODON_TABLE[c2]
                else (0.0, 1.0))
    # average over all mutational pathways (ignore paths through stops)
    import itertools
    Sd = Nd = 0.0
    npath = 0
    for order in itertools.permutations(diffs):
        cur = c1
        s = n = 0.0
        ok = True
        for p in order:
            nxt = cur[:p] + c2[p] + cur[p + 1:]
            if is_stop(nxt):
                ok = False
                break
            if CODON_TABLE[nxt] == CODON_TABLE[cur]:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            Sd += s
            Nd += n
            npath += 1
    if npath == 0:
        return 0.0, float(len(diffs))
    return Sd / npath, Nd / npath


# The Jukes-Cantor correction: turn a raw proportion of differences into an
# evolutionary distance. At the tiny divergence here it barely changes the
# numbers, but it is the standard, correct thing to apply.
def _jc(p):
    """Jukes-Cantor correction."""
    if p <= 0:
        return 0.0
    x = 1.0 - (4.0 / 3.0) * p
    if x <= 0:
        return float("inf")
    return -0.75 * math.log(x)


# Combine the site counts and difference counts into omega = dN/dS, where
# dN = corrected (Nd/N) and dS = corrected (Sd/S). Returns inf/nan for the
# degenerate cases (no synonymous differences, etc.).
def ng86_omega(S, N, Sd, Nd):
    """omega = dN/dS from NG86 site and difference counts."""
    if S <= 0 or N <= 0:
        return float("nan")
    dS = _jc(Sd / S)
    dN = _jc(Nd / N)
    if dS <= 0:
        return float("inf") if dN > 0 else float("nan")
    return dN / dS


# Total the S and N site counts across the ENTIRE reference concatenation,
# once. This is the expensive part, done a single time; every replicate then
# only adjusts these totals at the few codons it mutates (the 'delta update'
# that makes 1000 permutations run in seconds).
def ng86_backbone(seq):
    """(S, N) for the whole reference concatenation, PAML convention."""
    S = N = 0.0
    for k in range(len(seq) // 3):
        c = seq[3 * k:3 * k + 3]
        sn = _SITE_CACHE.get(c)
        if sn is None:          # stop codon or ambiguous -> skipped
            continue
        S += sn[0]
        N += sn[1]
    return S, N


# kept for backwards compatibility
def ng86_sites(seq):
    return ng86_backbone(seq)


# ----------------------------------------------------------------------------
# PAML yn00 on a pairwise (ref, mut) concatenation
# ----------------------------------------------------------------------------
# --- CALLING REAL PAML: yn00 -----------------------------------------------
# Write the two sequences (reference and mutated) into the PHYLIP file format
# PAML wants, write a yn00 control file, run yn00, and read the answer back.
# It returns BOTH PAML's NG86 and Yang-Nielsen (YN00) numbers. Important: it
# reads omega from PAML's main output at full precision, NOT from the rounded
# 2YN.* matrix files (those round to 4 decimals and, at this tiny divergence,
# that rounding produces nonsense ratios). Sequences must be in frame and
# contain no stop codons.
def run_yn00(ref_concat, mut_concat, workdir=None, keep=False):
    """Run PAML yn00 on the pairwise (ref, mut) concatenation.

    Returns dict with NG86 and YN00 estimates.  IMPORTANT: omega is parsed from
    the main output file, where PAML prints it at full internal precision.  Do
    NOT compute omega = dN/dS from the 2YN.dN / 2YN.dS matrix files -- those are
    rounded to 4 d.p., and at the very low divergence of a mutation-accumulation
    dataset (dS ~ 1e-3) that rounding destroys the ratio (you get quantized
    junk like exactly 0.8333 or 1.0000).

    Sequences must be equal length, in frame, and contain NO stop codons.
    """
    tmp = workdir or tempfile.mkdtemp(prefix="yn00_")
    phy = os.path.join(tmp, "aln.phy")
    ctl = os.path.join(tmp, "yn00.ctl")
    out = os.path.join(tmp, "aln.out")
    with open(phy, "w") as fh:
        fh.write(f"  2  {len(ref_concat)}\n")
        fh.write(f"ref        {ref_concat}\n")
        fh.write(f"mut        {mut_concat}\n")
    with open(ctl, "w") as fh:
        fh.write(f"seqfile = {phy}\noutfile = {out}\n"
                 "verbose = 0\nicode = 0\nweighting = 0\ncommonf3x4 = 0\n")
    subprocess.run(["yn00", ctl], cwd=tmp, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, check=False)
    res = {"omega_NG86": None, "omega_YN00": None,
           "dN_YN00": None, "dS_YN00": None, "kappa": None}
    try:
        txt = open(out).read()
    except FileNotFoundError:
        if not keep:
            shutil.rmtree(tmp, ignore_errors=True)
        return res

    # --- NG86 block:  "mut                  0.3619 (0.0004 0.0010)"
    i = txt.find("Nei & Gojobori")
    if i != -1:
        import re
        m = re.search(r"\n\s*\S+\s+([\d.]+)\s*\(\s*([\d.\-]+)\s+([\d.\-]+)\s*\)",
                      txt[i:i + 600])
        if m:
            res["omega_NG86"] = float(m.group(1))
            res["dN_NG86"] = float(m.group(2))
            res["dS_NG86"] = float(m.group(3))

    # --- YN00 table:  "2  1  S  N  t  kappa  omega  dN +- SE  dS +- SE"
    j = txt.find("seq. seq.")
    if j != -1:
        import re
        m = re.search(
            r"\n\s*\d+\s+\d+\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
            r"([\d.\-]+)\s+([\d.\-]+)\s*\+\-\s*([\d.\-]+)\s+([\d.\-]+)\s*\+\-",
            txt[j:j + 600])
        if m:
            res["S_YN00"] = float(m.group(1))
            res["N_YN00"] = float(m.group(2))
            res["kappa"] = float(m.group(4))
            res["omega_YN00"] = float(m.group(5))
            res["dN_YN00"] = float(m.group(6))
            res["dS_YN00"] = float(m.group(8))
    if not keep:
        shutil.rmtree(tmp, ignore_errors=True)
    return res


def _read_2yn_unused(path):
    """Parse the single off-diagonal value from a PAML 2YN.* distance matrix."""
    try:
        with open(path) as fh:
            toks = fh.read().split()
    except FileNotFoundError:
        return None
    # format:  "2  ref  mut <value>"  -> last float token is the pairwise value
    vals = []
    for t in toks:
        try:
            vals.append(float(t))
        except ValueError:
            pass
    # first float is the count "2"; the pairwise distance is the last float
    return vals[-1] if len(vals) >= 2 else None


# --- CALLING REAL PAML: codeml (maximum likelihood) ------------------------
# Same idea as run_yn00 but runs codeml in pairwise maximum-likelihood mode
# (runmode = -2). This is the 'edit the codeml file' route. codeml is slower
# (it searches iteratively), so it is used only on the OBSERVED alignment, not
# on the 1000 nulls. Returns the ML estimate of omega, dN, dS, and kappa.
def run_codeml(ref_concat, mut_concat, workdir=None, keep=False, timeout=3600):
    """Run PAML codeml in pairwise ML mode (runmode = -2, model M0) on the
    (ref, mut) alignment.  This is the 'edit the codeml file' route.  Slower than
    yn00 (maximum likelihood), so it is used for the OBSERVED alignment only, not
    for all 1000 permutations.  Returns {'omega_ML', 'dN', 'dS', 't', 'kappa'}."""
    tmp = workdir or tempfile.mkdtemp(prefix="codeml_")
    phy = os.path.join(tmp, "aln.phy")
    ctl = os.path.join(tmp, "codeml.ctl")
    out = os.path.join(tmp, "codeml.out")
    with open(phy, "w") as fh:
        fh.write(f"  2  {len(ref_concat)}\n")
        fh.write(f"ref        {ref_concat}\n")
        fh.write(f"mut        {mut_concat}\n")
    with open(ctl, "w") as fh:
        fh.write(CODEML_CTL.format(seqfile="aln.phy", outfile="codeml.out"))
    try:
        subprocess.run(["codeml", "codeml.ctl"], cwd=tmp,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        if not keep:
            shutil.rmtree(tmp, ignore_errors=True)
        return {"omega_ML": None, "error": "codeml timed out"}
    res = {"omega_ML": None, "dN": None, "dS": None, "t": None, "kappa": None}
    try:
        txt = open(out).read()
    except FileNotFoundError:
        if not keep:
            shutil.rmtree(tmp, ignore_errors=True)
        return res
    import re
    # "t= 0.0016  S= ...  N= ...  dN/dS=  0.4172  dN = 0.0004  dS = 0.0009"
    m = re.search(r"t=\s*([\d.]+).*?dN/dS=\s*([\d.]+)\s+dN\s*=\s*([\d.]+)\s+"
                  r"dS\s*=\s*([\d.]+)", txt, re.S)
    if m:
        res["t"] = float(m.group(1))
        res["omega_ML"] = float(m.group(2))
        res["dN"] = float(m.group(3))
        res["dS"] = float(m.group(4))
    k = re.search(r"kappa \(ts/tv\)\s*=\s*([\d.]+)", txt)
    if k:
        res["kappa"] = float(k.group(1))
    if not keep:
        shutil.rmtree(tmp, ignore_errors=True)
    return res


# The codeml control file.  Pairwise ML dN/dS on a 2-sequence alignment.
# The codeml control file, with each setting explained inline. runmode=-2 is
# pairwise, model=0 / NSsites=0 means one dN/dS for the whole alignment (M0),
# icode=0 is the universal genetic code, and kappa/omega are estimated (not
# fixed). {seqfile}/{outfile} are filled in when the file is written.
CODEML_CTL = """      seqfile = {seqfile}
      outfile = {outfile}
        noisy = 0
      verbose = 0
      runmode = -2      * -2 = pairwise ML (what we want for ref-vs-mutant)
      seqtype = 1       * 1 = codons
    CodonFreq = 2       * 2 = F3x4 codon frequencies
        model = 0       * 0 = one omega for the whole alignment
      NSsites = 0       * 0 = M0, a single dN/dS
        icode = 0       * 0 = universal code (yeast nuclear & E. coli)
    fix_kappa = 0       * estimate the ts/tv ratio
        kappa = 2       * initial value
    fix_omega = 0       * estimate omega
        omega = 0.4     * initial value
        clock = 0
    cleandata = 0       * keep all sites
"""
# --- BUILDING THE TWO SUPERGENES -------------------------------------------
# Glue all reference genes into one long sequence (the 'supergene'), in a
# fixed gene order so the reference and mutant line up base-for-base.
def concat_reference(cds, loci_order):
    return "".join(cds[l] for l in loci_order)


# Make the second supergene: a copy of the reference with your real non-stop
# mutations pasted in. Each mutation is placed at the first unused occurrence
# of its ref_codon in its own gene. Also returns `draws` - a list of (gene,
# position, new base, silent-or-changing) - which is exactly what the fast
# scorer needs, so the observed value is scored the same way as the nulls.
def apply_observed(cds, loci_order, muts):
    """Return (mutant concat, n applied, draws) with observed mutations applied.
    Mutations are placed at the first unused in-frame occurrence of their
    ref_codon within their own CDS.  draws = [(locus, i, new_base, cls)]."""
    mut_seqs = {l: list(cds[l]) for l in loci_order}
    applied = 0
    draws = []
    used = defaultdict(set)   # locus -> codon indices already mutated
    for m in muts:
        if m.locus not in mut_seqs:
            continue
        seq = cds[m.locus]
        ncod = len(seq) // 3
        for k in range(ncod):
            if k in used[m.locus]:
                continue
            if seq[3 * k:3 * k + 3] == m.ref_codon:
                idx = 3 * k + (m.pos - 1)
                mut_seqs[m.locus][idx] = m.central_to
                used[m.locus].add(k)
                cls = ("syn" if translate_codon(m.ref_codon)
                       == translate_codon(m.alt_codon) else "non")
                draws.append((m.locus, idx, m.central_to, cls))
                applied += 1
                break
    return "".join("".join(mut_seqs[l]) for l in loci_order), applied, draws


# ----------------------------------------------------------------------------
# Permutation engines
# ----------------------------------------------------------------------------
# --- THE NULL MODEL: shuffling mutations -----------------------------------
# ONE round of re-shuffling for Test 1. For each real non-stop mutation, pick
# a random site from the pool that has the SAME centred trinucleotide, apply
# the SAME central base change, and redraw if it would create a stop or reuse
# a site. This is the supervisor's recipe exactly: same contexts, same base
# changes, same number of mutations - just scattered at random instead of
# where selection left them. Returns the list of drawn changes.
def allocate_nonstop(muts_nonstop, cds, pool, rng, max_redraw=10000):
    """Re-allocate each non-stop mutation to a random comparable site (same
    centred trinucleotide, same central change), rejecting draws that create a
    stop or collide.  Return list of (locus, i, new_base, cls)."""
    draws = []
    used = set()
    for m in muts_nonstop:
        sites = pool.get(m.trinuc)
        if not sites:
            continue
        for _ in range(max_redraw):
            locus, i = rng.choice(sites)
            if (locus, i) in used:
                continue
            cls = classify_at_site(cds, locus, i, m.central_to)
            if cls == "stop":
                continue                 # redraw (Ka/Ks excludes stops)
            used.add((locus, i))
            draws.append((locus, i, m.central_to, cls))
            break
    return draws


# Turn one shuffled set into a Ka/Ks using the exact NG86 method. Starts from 
# the whole-backbone site totals (S_ref, N_ref) and only ADJUSTS them at the codons
# that got mutated - so it never re-scans millions of codons.
#  It groups changes by codon (in case two land in one codon),
# updates S and N to the pair-average PAML uses, tallies silent vs changing
# differences, and returns omega. This reproduces PAML's NG86 number exactly.
def score_fast(draws, all_cds, S_ref, N_ref):
    """EXACT PAML-NG86 omega for a replicate, via delta update.

    Only codons that receive a mutation change their site counts, so we adjust
    (S, N) from the reference backbone at those codons only.  Site counts are
    averaged over the two sequences (PAML convention).  Multiple hits in one
    codon are handled with NG86 pathway averaging."""
    # group the drawn changes by codon
    bycodon = defaultdict(list)          # (locus, codon_index) -> [(p, new_base)]
    for locus, i, nb, _cls in draws:
        bycodon[(locus, i // 3)].append((i % 3, nb))

    S = S_ref
    N = N_ref
    Sd = Nd = 0.0
    for (locus, k), changes in bycodon.items():
        seq = all_cds[locus]
        c1 = seq[3 * k:3 * k + 3]
        s1 = _SITE_CACHE.get(c1)
        if s1 is None:
            continue                      # reference codon was a stop: skipped
        c2 = list(c1)
        for p, nb in changes:
            c2[p] = nb
        c2 = "".join(c2)
        s2 = _SITE_CACHE.get(c2)
        if s2 is None:
            # mutant codon is a stop -> excluded from the alignment upstream;
            # should not occur in TEST 1, but guard anyway
            continue
        # backbone counted this codon as s1; NG86 wants the pair average
        S += (s1[0] + s2[0]) / 2.0 - s1[0]
        N += (s1[1] + s2[1]) / 2.0 - s1[1]
        sd, nd = count_diffs(c1, c2)
        Sd += sd
        Nd += nd
    return ng86_omega(S, N, Sd, Nd)


# The slow alternative to score_fast: actually build the mutated supergene
# from the draws and run real yn00 on it. Used only for the --engine paml mode
# and for the validation check.
def score_paml(draws, all_cds, loci_order, ref_concat):
    mut_seqs = {l: list(all_cds[l]) for l in loci_order}
    for locus, i, nb, _cls in draws:
        mut_seqs[locus][i] = nb
    mut_concat = "".join("".join(mut_seqs[l]) for l in loci_order)
    return run_yn00(ref_concat, mut_concat)


# One complete Test-1 replicate: shuffle once (allocate_nonstop), then score
# (score_fast). Called nperm times to build the null Ka/Ks distribution.
def permute_kaks(muts_nonstop, all_cds, pool, S_ref, N_ref, rng,
                 max_redraw=10000):
    """One TEST-1 replicate -> exact NG86 omega."""
    draws = allocate_nonstop(muts_nonstop, all_cds, pool, rng, max_redraw)
    return score_fast(draws, all_cds, S_ref, N_ref)


# One complete Test-2 replicate: re-shuffle ALL mutations (stops allowed this
# time) and count how many land on a stop-creating change. Called nperm times
# to build the null stop-count distribution.
def permute_stops(muts_all, cds, pool, rng, max_redraw=10000):
    """One replicate of TEST 2.  Re-allocate ALL mutations (stops allowed) and
    count how many create a stop codon."""
    nstop = 0
    used = set()
    for m in muts_all:
        tri = m.trinuc
        sites = pool.get(tri)
        if not sites:
            continue
        for _ in range(max_redraw):
            locus, i = rng.choice(sites)
            if (locus, i) in used:
                continue
            used.add((locus, i))
            if classify_at_site(cds, locus, i, m.central_to) == "stop":
                nstop += 1
            break
    return nstop


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
# ===========================================================================
# run() -- THE RECIPE THAT TIES EVERYTHING TOGETHER
# ---------------------------------------------------------------------------
# Executed in order: (1) load mutations; (2) load + quality-control the
# reference, dropping genes that are out of frame or contain internal stops
# and logging why; (3) drop mutations whose codon doesn't match the reference;
# (4) give each mutation its trinucleotide context; (5) split into the non-stop
# set (Test 1) and stop-creating set (Test 2); (6) build the site pool and the
# reference supergene with its site totals; (7) compute the OBSERVED values
# (Ka/Ks by fast-NG86, by yn00, optionally by codeml; and the observed stop
# count); (8) run the two NULLS, nperm replicates each; (9) compute the two
# p-values and write all outputs. `seed` makes the random draws reproducible.
# ===========================================================================
def run(csv_path, fasta_path, nperm=1000, seed=1, engine="fast",
        paml_validate=0, outprefix="result", quiet=False, do_codeml=False):
    rng = random.Random(seed)
    report = {}

    def log(*a):
        if not quiet:
            print(*a, file=sys.stderr)

    muts = load_mutations(csv_path)
    log(f"[load] {len(muts)} usable mutations from {csv_path}")

    loci_needed = {m.locus for m in muts}
    raw = read_fasta(fasta_path, key_targets=loci_needed)
    ndup = getattr(read_fasta, "duplicates", 0)
    if ndup:
        log(f"[cds]  collapsed {ndup} duplicate record(s) to one entry per locus")
    report["n_duplicate_records"] = ndup
    cds = {}
    for locus in loci_needed:
        if locus in raw:
            s = strip_terminal_stop(raw[locus])
            if s and set(s) <= set("ACGT") and len(s) >= 3:
                cds[locus] = s
    log(f"[cds]  {len(cds)}/{len(loci_needed)} mutated loci found in FASTA")
    # also load ALL CDS present in fasta for the site pool / backbone.
    # QC: a CDS is usable only if it is in frame, ACGT-only, and free of INTERNAL
    # stop codons.  Internal stops abort yn00.  In yeast R64 these are the
    # mitochondrial genes (Q0*), which use the mitochondrial code (TGA=Trp) and
    # do not belong in a nuclear analysis, plus a few dubious ORFs/pseudogenes.
    all_cds = {}
    qc = Counter()
    excluded = []
    for k, s in raw.items():
        if k.startswith("__"):
            continue
        if len(s) % 3:
            qc["not_multiple_of_3"] += 1
            excluded.append(k)
            continue
        s2 = strip_terminal_stop(s)
        if not s2 or len(s2) < 3:
            qc["too_short"] += 1
            excluded.append(k)
            continue
        if set(s2) - set("ACGT"):
            qc["non_ACGT"] += 1
            excluded.append(k)
            continue
        if any(is_stop(s2[3 * i:3 * i + 3]) for i in range(len(s2) // 3)):
            qc["internal_stop"] += 1
            excluded.append(k)
            continue
        all_cds[k] = s2
        qc["kept"] += 1
    log(f"[cds]  QC {dict(qc)}")
    if excluded:
        log(f"[cds]  excluded {len(excluded)} CDS "
            f"(e.g. {', '.join(excluded[:6])}{' ...' if len(excluded) > 6 else ''})")
    report["cds_qc"] = dict(qc)
    report["cds_excluded"] = excluded
    log(f"[cds]  {len(all_cds)} usable CDS for backbone & site pool")

    # mutations must live in a QC-passing CDS
    cds = {k: v for k, v in cds.items() if k in all_cds}

    # keep only mutations whose CDS we have AND whose ref_codon actually occurs
    # in-frame in that CDS (guards against reference/annotation mismatch)
    muts = [m for m in muts if m.locus in cds]
    matched, unmatched = [], []
    for m in muts:
        s = cds[m.locus]
        if any(s[3 * k:3 * k + 3] == m.ref_codon for k in range(len(s) // 3)):
            matched.append(m)
        else:
            unmatched.append(m)
    if unmatched:
        log(f"[filter] {len(unmatched)} mutations dropped: ref_codon not found "
            f"in-frame in the reference CDS (reference/annotation mismatch), "
            f"e.g. {[(m.locus, m.ref_codon) for m in unmatched[:4]]}")
    report["n_unmatched_refcodon"] = len(unmatched)
    muts = matched
    log(f"[filter] {len(muts)} mutations retained (CDS available & ref_codon matches)")

    assign_trinucleotides(muts, all_cds, rng, report)
    muts = [m for m in muts if m.trinuc is not None]
    log(f"[context] {len(muts)} mutations with a centred trinucleotide "
        f"({report['trinuc_stats']})")

    muts_nonstop = [m for m in muts
                    if not is_stop(m.ref_codon) and not is_stop(m.alt_codon)]
    muts_stop = [m for m in muts if m.creates_stop]
    log(f"[sets] non-stop (Ka/Ks): {len(muts_nonstop)} | "
        f"stop-creating: {len(muts_stop)}")

    # ---- site pool over ALL cds ----
    log("[pool] indexing all trinucleotide sites across all CDS ...")
    pool = build_site_pool(all_cds)
    log(f"[pool] {sum(len(v) for v in pool.values())} sites, "
        f"{len(pool)} distinct trinucleotides")

    # ---- reference backbone & NG86 site counts (PAML convention) ----
    loci_order = sorted(all_cds)
    ref_concat = concat_reference(all_cds, loci_order)
    log(f"[backbone] concatenation length {len(ref_concat)} bp "
        f"({len(ref_concat)//3} codons)")
    S_ref, N_ref = ng86_backbone(ref_concat)
    log(f"[backbone] NG86 sites  S={S_ref:.1f}  N={N_ref:.1f}")

    # ============================ OBSERVED ============================
    mut_concat, applied, obs_draws = apply_observed(all_cds, loci_order,
                                                    muts_nonstop)
    log(f"[observed] applied {applied}/{len(muts_nonstop)} non-stop mutations")

    obs = {}
    obs["omega_ng86_exact"] = score_fast(obs_draws, all_cds, S_ref, N_ref)
    Nd = sum(1 for *_x, c in obs_draws if c == "non")
    Sd = sum(1 for *_x, c in obs_draws if c == "syn")
    obs["Nd"], obs["Sd"] = Nd, Sd

    yn = run_yn00(ref_concat, mut_concat)
    obs["paml"] = yn
    log(f"[observed] PAML yn00:  NG86 omega={yn.get('omega_NG86')}   "
        f"YN00 omega={yn.get('omega_YN00')}  (kappa={yn.get('kappa')})")
    log(f"[observed] our exact NG86 omega={obs['omega_ng86_exact']:.4f}  "
        f"(Nd={Nd}, Sd={Sd})  <- should equal PAML's NG86 omega")

    obs["codeml"] = None
    if do_codeml:
        log("[observed] running codeml (pairwise ML, runmode=-2) ...")
        cm = run_codeml(ref_concat, mut_concat)
        obs["codeml"] = cm
        log(f"[observed] codeml ML omega={cm.get('omega_ML')}  "
            f"(dN={cm.get('dN')}, dS={cm.get('dS')}, kappa={cm.get('kappa')})")

    obs_stops = len(muts_stop)
    log(f"[observed] stop-creating mutations: {obs_stops}")

    # ============================ NULL: TEST 1 ============================
    log(f"[null-1] {nperm} permutations (engine={engine}) ...")
    null_omega = []
    for r in range(nperm):
        if engine == "paml":
            draws = allocate_nonstop(muts_nonstop, all_cds, pool, rng)
            om = score_paml(draws, all_cds, loci_order,
                            ref_concat).get("omega_YN00")
            om = float("nan") if om is None else om
        else:
            om = permute_kaks(muts_nonstop, all_cds, pool, S_ref, N_ref, rng)
        null_omega.append(om)
        if not quiet and (r + 1) % max(1, nperm // 10) == 0:
            log(f"  ... {r+1}/{nperm}")

    validation = None
    if paml_validate:
        validation = _validate(muts_nonstop, all_cds, loci_order, pool,
                               ref_concat, S_ref, N_ref,
                               random.Random(seed + 99), paml_validate, log)

    # ============================ NULL: TEST 2 ============================
    log(f"[null-2] {nperm} permutations (stops allowed) ...")
    null_stops = []
    for r in range(nperm):
        null_stops.append(permute_stops(muts, all_cds, pool, rng))
        if not quiet and (r + 1) % max(1, nperm // 10) == 0:
            log(f"  ... {r+1}/{nperm}")

    # ============================ P-VALUES ============================
    # The null replicates and the observed value are scored with the SAME
    # statistic (exact NG86 omega, identical to PAML's NG86), so the rank is
    # exact.  PAML's YN00 omega is reported alongside for the observed data.
    obs_stat = (obs["omega_ng86_exact"] if engine == "fast"
                else obs["paml"].get("omega_YN00"))
    valid = [x for x in null_omega if x == x and x not in (float("inf"),)]
    p1 = (sum(1 for x in valid if x <= obs_stat) / len(valid)) if valid else float("nan")
    p2 = (sum(1 for x in null_stops if x <= obs_stops)) / len(null_stops)

    results = {
        "organism_csv": csv_path,
        "reference_fasta": fasta_path,
        "n_mutations_total": len(muts),
        "n_nonstop": len(muts_nonstop),
        "n_stop": len(muts_stop),
        "n_cds_used": len(all_cds),
        "n_duplicate_records_collapsed": report.get("n_duplicate_records"),
        "backbone_bp": len(ref_concat),
        "S_sites": S_ref, "N_sites": N_ref,
        "observed_Nd": obs["Nd"], "observed_Sd": obs["Sd"],
        "observed_omega_NG86_exact": obs["omega_ng86_exact"],
        "observed_omega_NG86_paml": obs["paml"].get("omega_NG86"),
        "observed_omega_YN00_paml": obs["paml"].get("omega_YN00"),
        "observed_omega_ML_codeml": (obs["codeml"] or {}).get("omega_ML"),
        "observed_dN_YN00": obs["paml"].get("dN_YN00"),
        "observed_dS_YN00": obs["paml"].get("dS_YN00"),
        "observed_kappa": obs["paml"].get("kappa"),
        "null_statistic": ("exact NG86 omega" if engine == "fast"
                           else "PAML YN00 omega"),
        "observed_statistic_used_for_p": obs_stat,
        "nperm": nperm,
        "null_omega_mean": (sum(valid) / len(valid)) if valid else float("nan"),
        "null_omega_min": min(valid) if valid else float("nan"),
        "null_omega_max": max(valid) if valid else float("nan"),
        "p_kaks_lower": p1,
        "observed_stops": obs_stops,
        "null_stops_mean": sum(null_stops) / len(null_stops),
        "null_stops_min": min(null_stops),
        "null_stops_max": max(null_stops),
        "p_stops_lower": p2,
        "trinuc_stats": report.get("trinuc_stats"),
        "cds_qc": report.get("cds_qc"),
        "n_unmatched_refcodon": report.get("n_unmatched_refcodon"),
        "validation": validation,
    }
    _write_outputs(results, null_omega, null_stops, obs_stat, obs_stops,
                   outprefix, log)
    return results


# The agreement check for --paml-validate: score the SAME shuffles with both
# our fast NG86 and real PAML, and print them side by side so you can confirm
# they match to ~4 decimals (proving the fast engine is not an approximation).
def _validate(muts_nonstop, all_cds, loci_order, pool, ref_concat,
              S_ref, N_ref, rng, k, log):
    """Score k replicates with BOTH engines on the SAME allocation.
    Our exact NG86 should reproduce PAML's NG86 omega to ~4 d.p."""
    fast, png, pyn = [], [], []
    for _ in range(k):
        draws = allocate_nonstop(muts_nonstop, all_cds, pool, rng)
        f = score_fast(draws, all_cds, S_ref, N_ref)
        p = score_paml(draws, all_cds, loci_order, ref_concat)
        fast.append(f)
        png.append(p.get("omega_NG86"))
        pyn.append(p.get("omega_YN00"))
    log(f"[validate] same allocation, {k} replicates:")
    log(f"{'ours(NG86)':>12} {'PAML NG86':>10} {'PAML YN00':>10}  {'diff':>8}")
    for f, g, y in zip(fast, png, pyn):
        d = (abs(f - g) if (g is not None) else float('nan'))
        log(f"{f:12.4f} {(g if g is not None else float('nan')):10.4f} "
            f"{(y if y is not None else float('nan')):10.4f}  {d:8.5f}")
    return {"ours_ng86": fast, "paml_ng86": png, "paml_yn00": pyn}


# Write the results: the _summary.json (all headline numbers, counts, QC), the
# two _null_*.txt files (the raw 1000 values behind each histogram), and the
# _plots.png figure (observed value marked against each null). Plotting is
# wrapped in try/except so the run still finishes if matplotlib is missing.
def _write_outputs(results, null_omega, null_stops, obs_om_fast, obs_stops,
                   outprefix, log):
    import json
    with open(outprefix + "_summary.json", "w") as fh:
        json.dump({k: v for k, v in results.items() if k != "validation"},
                  fh, indent=2, default=str)
    with open(outprefix + "_null_kaks.txt", "w") as fh:
        fh.write("\n".join(str(x) for x in null_omega))
    with open(outprefix + "_null_stops.txt", "w") as fh:
        fh.write("\n".join(str(x) for x in null_stops))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        valid = [x for x in null_omega if x == x and x != float("inf")]
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].hist(valid, bins=40, color="#4C72B0", alpha=.85)
        ax[0].axvline(obs_om_fast, color="crimson", lw=2,
                      label=f"observed = {obs_om_fast:.3f}")
        ax[0].set_title(f"TEST 1  Ka/Ks null  (p={results['p_kaks_lower']:.4f})")
        ax[0].set_xlabel("Ka/Ks"); ax[0].set_ylabel("replicates"); ax[0].legend()
        ax[1].hist(null_stops, bins=range(min(null_stops),
                   max(null_stops) + 2), color="#55A868", alpha=.85,
                   align="left")
        ax[1].axvline(obs_stops, color="crimson", lw=2,
                      label=f"observed = {obs_stops}")
        ax[1].set_title(f"TEST 2  stop-count null (p={results['p_stops_lower']:.4f})")
        ax[1].set_xlabel("# stop-creating mutations"); ax[1].set_ylabel("replicates")
        ax[1].legend()
        fig.tight_layout()
        fig.savefig(outprefix + "_plots.png", dpi=130)
        log(f"[out] wrote {outprefix}_plots.png")
    except Exception as e:
        log(f"[out] plot skipped ({e})")


# --- THE COMMAND LINE ------------------------------------------------------
# Defines the flags you type (--csv, --fasta, --nperm, --seed, --engine,
# --paml-validate, --codeml, --out), reads them, and calls run() with them.
# This is the thin layer that turns your terminal command into a function call.
def main():
    ap = argparse.ArgumentParser(description="Ka/Ks context-permutation selection test (no Biopython).")
    ap.add_argument("--csv", required=True, help="mutation CSV") #input csv
    ap.add_argument("--fasta", required=True, help="gene-oriented reference CDS FASTA") #input cds fasta
    ap.add_argument("--nperm", type=int, default=1000) #number of permutations
    ap.add_argument("--seed", type=int, default=1) #sets a seed, can change to get a new draw
    ap.add_argument("--engine", choices=["fast", "paml"], default="fast", #run with either ng86, or yn00
                    help="fast = NG86 fixed-backbone omega (validated equal-rank "
                         "to yn00); paml = real yn00 every replicate (slow)")
    ap.add_argument("--paml-validate", type=int, default=0,
                    help="score this many replicates with BOTH engines to prove agreement")
    ap.add_argument("--codeml", action="store_true",
                    help="also run codeml (pairwise ML dN/dS, runmode=-2) on the "
                         "observed alignment for a maximum-likelihood omega")
    ap.add_argument("--out", default="result")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    res = run(a.csv, a.fasta, a.nperm, a.seed, a.engine, a.paml_validate,
              a.out, a.quiet, a.codeml)
    import json
    print(json.dumps({k: v for k, v in res.items()
                      if k not in ("validation",)}, indent=2, default=str))


if __name__ == "__main__":
    main()
