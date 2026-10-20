"""Network gene panel for cell-specific CSCN graphs, selected from NTC cells only.

The panel is the forced genes (eligible targets and the primary Fufa program genes) followed by fill
genes ranked by a Seurat-style mean.var.plot standardization, computed on log1p CP10k with
equal-frequency mean bins. Expression of target-perturbed cells never enters any statistic: the only
entry point that accepts all cells (build_panel) drops non-NTC rows before any arithmetic, and every
other function raises when given non-NTC cells.
"""

import re

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .gene_sets import scoring_genes, unversioned

NTC = "NTC"
EXCLUDE_PATTERN = r"^MT-|^(RP[LS]\d+([AXY]|AL|L)?\d*|RPLP[0-2]|RPSA)(P\d+)?$"


def ntc_only(meta):
    """Raise unless every cell is a non-targeting control."""
    if not meta.target.eq(NTC).all():
        raise ValueError("Network panel and graph inputs must contain NTC cells only")


def log_cp10k(counts, total):
    """log1p counts per 10k; same arithmetic as scripts/aim1.py::log_cp10k (float32, sparse CSR)."""
    X = sp.csr_matrix(counts, dtype=np.float32, copy=True)
    X.data = np.log1p(X.data / np.repeat(np.asarray(total), np.diff(X.indptr)).astype(np.float32)
                      * np.float32(1e4))
    return X


def ntc_gene_stats(X, meta, conditions):
    """Per condition: detection (> 0), mean and variance (ddof=1) of every column over NTC cells."""
    ntc_only(meta)
    X = sp.csr_matrix(X)
    out = {}
    for condition in conditions:
        rows = np.flatnonzero(meta.condition.eq(condition).to_numpy())
        if len(rows) < 2:
            raise ValueError(f"Fewer than two NTC cells in {condition}")
        sub = X[rows]
        n = len(rows)
        total = np.asarray(sub.sum(axis=0, dtype=np.float64)).ravel()
        squares = sub.copy()
        squares.data = squares.data.astype(np.float64) ** 2
        total_sq = np.asarray(squares.sum(axis=0)).ravel()
        mean = total / n
        out[f"detection_{condition}"] = np.asarray((sub > 0).sum(axis=0)).ravel() / n
        out[f"mean_{condition}"] = mean
        out[f"var_{condition}"] = np.maximum(total_sq - n * mean ** 2, 0) / (n - 1)
    return pd.DataFrame(out)


def excluded_symbol(symbols, pattern=EXCLUDE_PATTERN):
    """Mitochondrial and ribosomal-protein symbols (RPL/RPS, RPLP0-2, RPSA and their pseudogenes)."""
    rule = re.compile(pattern)
    return np.array([bool(rule.search(str(s))) for s in symbols])


def dispersion_z(mean, var, reference, n_bins=20):
    """Standardized dispersion log(var / mean) within equal-frequency mean bins.

    Bin edges are mean quantiles of the reference genes, and each bin's centre and SD come from
    the reference genes in it. Genes with zero mean or variance get NaN; a bin with SD 0 gives z 0.
    """
    mean, var, reference = np.asarray(mean, float), np.asarray(var, float), np.asarray(reference, bool)
    defined = (mean > 0) & (var > 0)
    if (reference & ~defined).any():
        raise ValueError("Reference genes need positive mean and variance")
    d = np.full(len(mean), np.nan)
    d[defined] = np.log(var[defined] / mean[defined])
    edges = np.quantile(mean[reference], np.linspace(0, 1, n_bins + 1))
    bins = np.clip(np.searchsorted(edges[1:-1], mean, side="right"), 0, n_bins - 1)
    z = np.full(len(mean), np.nan)
    for b in range(n_bins):
        ref = reference & (bins == b)
        if not ref.any():
            continue
        centre = d[ref].mean()
        sd = d[ref].std(ddof=1) if ref.sum() > 1 else 0.0
        members = defined & (bins == b)
        z[members] = (d[members] - centre) / sd if sd > 0 else 0.0
    return z


def rank_fill_genes(stats, universe, conditions):
    """Rank universe genes by the mean of their within-condition dispersion ranks (1 = most dispersed).

    stats: one row per gene with gene_id and dispersion_z_<condition>. Ties: mean z, then gene_id.
    """
    sub = stats[np.asarray(universe, bool)].copy()
    for condition in conditions:
        order = np.lexsort((sub.gene_id.to_numpy(), -sub[f"dispersion_z_{condition}"].to_numpy()))
        rank = np.empty(len(sub), int)
        rank[order] = np.arange(1, len(sub) + 1)
        sub[f"dispersion_rank_{condition}"] = rank
    sub["mean_rank"] = sub[[f"dispersion_rank_{c}" for c in conditions]].mean(axis=1)
    mean_z = sub[[f"dispersion_z_{c}" for c in conditions]].mean(axis=1)
    order = np.lexsort((sub.gene_id.to_numpy(), -mean_z.to_numpy(), sub.mean_rank.to_numpy()))
    sub = sub.iloc[order].reset_index(drop=True)
    sub["fill_rank"] = np.arange(1, len(sub) + 1)
    return sub


def panel_of_size(forced, ranked, size):
    """Forced genes followed by the top size - len(forced) ranked fill genes."""
    forced, ranked = list(forced), list(ranked)
    if size < len(forced):
        raise ValueError(f"Panel size {size} is smaller than the {len(forced)} forced genes")
    if size - len(forced) > len(ranked):
        raise ValueError(f"Panel size {size} exceeds forced plus ranked fill genes")
    if set(forced) & set(ranked):
        raise ValueError("Forced and fill genes overlap")
    return forced + ranked[:size - len(forced)]


def build_panel(counts, meta, var, forced, params):
    """Panel statistics from all cells' raw counts, using NTC rows only.

    counts: cells x genes raw counts (sparse); meta: condition, target, total_counts (whole
    transcriptome) per cell; var: indexed by gene_id with symbol, gene_type, is_par_y; forced:
    DataFrame with gene_id (as in var.index), symbol, role, target, in_fufa in panel order.
    Returns forced (with statistics and the low_detection flag), ranked fill genes and genes (the
    statistics of every scoring gene).
    """
    is_ntc = meta.target.eq(NTC).to_numpy()
    counts = sp.csr_matrix(counts)[is_ntc]
    meta = meta[is_ntc]
    ntc_only(meta)
    conditions = list(params["conditions"])
    keep = scoring_genes(var)
    genes = var[keep].reset_index().rename(columns={var.index.name or "index": "gene_id"})
    missing = set(forced.gene_id) - set(genes.gene_id)
    if missing or forced.gene_id.duplicated().any():
        raise ValueError(f"Forced genes must be unique scoring genes; missing: {sorted(missing)}")
    X = log_cp10k(counts[:, np.flatnonzero(keep)], meta.total_counts.to_numpy(float))
    stats = pd.concat([genes[["gene_id", "symbol"]], ntc_gene_stats(X, meta, conditions)], axis=1)
    is_forced = stats.gene_id.isin(forced.gene_id).to_numpy()
    detected = np.all([stats[f"detection_{c}"].to_numpy() >= params["fill_min_detection"]
                       for c in conditions], axis=0)
    excluded = excluded_symbol(stats.symbol, params["exclude_pattern"])
    universe = detected & ~excluded & ~is_forced
    for condition in conditions:
        stats[f"dispersion_z_{condition}"] = dispersion_z(
            stats[f"mean_{condition}"], stats[f"var_{condition}"], universe, params["n_bins"])
    ranked = rank_fill_genes(stats, universe, conditions)
    forced_stats = forced.merge(stats.drop(columns="symbol"), on="gene_id", how="left")
    forced_stats["low_detection"] = np.all(
        [forced_stats[f"detection_{c}"].to_numpy() < params["forced_low_detection"] for c in conditions], axis=0)
    return dict(forced=forced_stats, ranked=ranked, genes=stats,
                n_universe=int(universe.sum()), n_detected=int(detected.sum()),
                n_excluded_detected=int((detected & excluded & ~is_forced).sum()),
                n_ntc={c: int(meta.condition.eq(c).sum()) for c in conditions})


def forced_genes(targets, target_gene_ids, fufa_genes, var):
    """Forced panel genes in panel order: targets (by name), then Fufa genes not already present.

    targets: eligible target names; target_gene_ids: target -> gene_id; fufa_genes: program rows
    (symbol, gene_id, duplicate_in_program) from config/programs.json; var: indexed by versioned
    gene_id with symbol. Gene IDs are matched unversioned against the scoring genes.
    """
    keep = scoring_genes(var)
    lookup = pd.Series(var.index[keep], index=[unversioned(g) for g in var.index[keep]])
    symbol = var.symbol
    rows = []
    for target in sorted(targets):
        gene_id = lookup.get(unversioned(target_gene_ids[target]))
        if gene_id is None:
            raise ValueError(f"Target {target} has no scoring gene")
        rows.append(dict(gene_id=gene_id, symbol=symbol[gene_id], role="target", target=target))
    fufa_ids = []
    for g in fufa_genes:
        if g.get("duplicate_in_program") or not isinstance(g.get("gene_id"), str):
            continue
        gene_id = lookup.get(unversioned(g["gene_id"]))
        if gene_id is None:
            raise ValueError(f"Fufa gene {g['symbol']} has no scoring gene")
        fufa_ids.append(gene_id)
    present = {r["gene_id"] for r in rows}
    rows += [dict(gene_id=g, symbol=symbol[g], role="Fufa", target="") for g in fufa_ids if g not in present]
    frame = pd.DataFrame(rows)
    frame["in_fufa"] = frame.gene_id.isin(fufa_ids)
    if frame.gene_id.duplicated().any():
        raise ValueError("Two targets map to the same gene")
    return frame


def load_ntc_expression(path, gene_ids, conditions):
    """NTC cells' log1p CP10k for the given genes (dense float32) and their metadata.

    Non-NTC rows are dropped before normalization; total_counts is the whole-transcriptome depth.
    """
    import anndata as ad

    adata = ad.read_h5ad(path)
    obs = adata.obs
    keep = (obs.target.astype(str).eq(NTC) & obs.condition.astype(str).isin(conditions)).to_numpy()
    meta = obs[keep].copy()
    for col in ("condition", "target", "lig1_well"):
        meta[col] = meta[col].astype(str)
    position = pd.Series(np.arange(adata.n_vars), index=adata.var.index)
    columns = position.reindex(list(gene_ids))
    if columns.isna().any():
        raise ValueError(f"Genes missing from the expression matrix: {list(columns[columns.isna()].index)}")
    counts = sp.csr_matrix(adata.X[np.flatnonzero(keep)])
    ntc_only(meta)
    X = log_cp10k(counts[:, columns.to_numpy(int)], meta.total_counts.to_numpy(float))
    return np.asarray(X.todense(), dtype=np.float32), meta
