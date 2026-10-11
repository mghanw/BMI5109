"""Pseudobulk target-vs-NTC expression contrasts on sparse log1p CP10k matrices."""

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from .dispersion import bh


def indicator(labels, levels):
    """Sparse groups x cells indicator for the given label levels (other cells ignored)."""
    labels = list(labels)
    position = {level: i for i, level in enumerate(levels)}
    rows = np.array([position.get(v, -1) for v in labels])
    keep = rows >= 0
    return sp.csr_matrix((np.ones(keep.sum()), (rows[keep], np.flatnonzero(keep))),
                         shape=(len(levels), len(labels)))


def group_moments(X, labels, levels):
    """Per group: cell count, mean, mean of squares and detection rate of every column."""
    X = sp.csr_matrix(X)
    G = indicator(labels, levels)
    n = np.asarray(G.sum(axis=1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        inv = sp.diags(np.where(n > 0, 1 / n, 0))
        squares = X.copy()
        squares.data = squares.data.astype(np.float64) ** 2
        detected = X.copy()
        detected.data = (detected.data > 0).astype(np.float64)
        mean = np.asarray((inv @ G @ X).todense())
        mean_sq = np.asarray((inv @ G @ squares).todense())
        detection = np.asarray((inv @ G @ detected).todense())
    return n, mean, mean_sq, detection


def detected_genes(X, meta, min_detection, conditions=("DMSO", "Vem")):
    """Columns detected (> 0) in at least min_detection of NTC cells, per condition."""
    labels = np.where(meta.target.eq("NTC"), meta.condition, None)
    _, _, _, detection = group_moments(X, labels, list(conditions))
    return pd.DataFrame(detection.T >= min_detection, columns=[f"detected_ntc_{c}" for c in conditions])


def welch(mean_a, var_a, n_a, mean_b, var_b, n_b):
    va, vb = var_a / n_a, var_b / n_b
    se = np.sqrt(va + vb)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (mean_a - mean_b) / se
        df = (va + vb) ** 2 / (va ** 2 / (n_a - 1) + vb ** 2 / (n_b - 1))
        p = 2 * stats.t.sf(np.abs(t), df)
    p = np.where(se > 0, p, np.nan)
    return t, p


def expression_contrasts(X, genes, meta, targets, target_gene_ids, conditions=("DMSO", "Vem")):
    """Long table of target - NTC differences in mean log1p CP10k per (condition, target, gene).

    X: cells x genes (already restricted to the analysed genes), genes: DataFrame with gene_id, symbol.
    target_gene_ids: target -> perturbed gene_id; that gene is excluded from the target's rows.
    BH within each (condition, target).
    """
    levels = [(c, t) for c in conditions for t in ["NTC", *targets]]
    labels = list(zip(meta.condition, meta.target))
    n, mean, mean_sq, _ = group_moments(X, labels, levels)
    var = (mean_sq - mean ** 2) * (n / np.maximum(n - 1, 1))[:, None]
    base = lambda g: g.split(".")[0]
    gene_bases = genes.gene_id.map(base).to_numpy()
    frames = []
    for condition in conditions:
        r = levels.index((condition, "NTC"))
        for target in targets:
            i = levels.index((condition, target))
            keep = gene_bases != base(target_gene_ids[target])
            t, p = welch(mean[i], var[i], n[i], mean[r], var[r], n[r])
            frame = pd.DataFrame({"condition": condition, "target": target, "gene_id": genes.gene_id.to_numpy(),
                                  "symbol": genes.symbol.to_numpy(), "mean_target": mean[i], "mean_ntc": mean[r],
                                  "diff": mean[i] - mean[r], "welch_t": t, "p_value": p,
                                  "n_target": int(n[i]), "n_ntc": int(n[r])})[keep]
            frame["q_value"] = bh(frame.p_value)
            frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def summarize_expression(contrasts, cutoff, q=0.05, top=10):
    """Counts of genes with |diff| >= cutoff (all / q < q) and the top genes per (condition, target)."""
    big = contrasts["diff"].abs().ge(cutoff)
    sig = contrasts.q_value.lt(q)
    summary = contrasts.assign(big=big, big_sig=big & sig, sig=sig, up=big & sig & contrasts["diff"].gt(0)).groupby(
        ["condition", "target"]).agg(n_genes_tested=("gene_id", "size"), n_q_lt_005=("sig", "sum"),
                                     n_absdiff_ge_cutoff=("big", "sum"), n_absdiff_ge_cutoff_q_lt_005=("big_sig", "sum"),
                                     n_up_ge_cutoff_q_lt_005=("up", "sum"), n_target=("n_target", "first"),
                                     n_ntc=("n_ntc", "first")).reset_index()
    summary["n_down_ge_cutoff_q_lt_005"] = summary.n_absdiff_ge_cutoff_q_lt_005 - summary.n_up_ge_cutoff_q_lt_005
    summary["absdiff_cutoff"] = cutoff
    ranked = contrasts[sig].assign(absdiff=contrasts["diff"].abs()).sort_values(
        ["condition", "target", "absdiff", "symbol"], ascending=[True, True, False, True])
    tops = ranked.groupby(["condition", "target"]).head(top).copy()
    tops["rank"] = tops.groupby(["condition", "target"]).cumcount() + 1
    tops = tops[["condition", "target", "rank", "symbol", "gene_id", "diff", "mean_target", "mean_ntc", "q_value"]]
    return summary, tops.reset_index(drop=True)
