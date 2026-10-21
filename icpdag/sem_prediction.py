"""Linear SEM coefficients and knockdown simulation on fold-specific DAGs.

Equation per node j: x_j = a_j + sum_{l in pa(j)} b_jl x_l + e_j on log1p CP10k. Local coefficients
are ridge-shrunk towards the pooled regression with the same parent set:
    b = argmin sum_{i in N} (x_ij - a - b'x_i,pa)^2 + lambda |N| ||b - b_pooled||^2.
A knockdown of g by delta gives Delta x = (I - B)^-1 e_g delta, computed by forward substitution.
"""

import networkx as nx
import numpy as np

from .sem_orientation import local_ida_choices, orient, undirected_neighbours

MODELS = ["cell_specific_cscn_sem", "shared_cscn_sem", "cluster_level_sem", "local_association_network_sem",
          "unweighted_propagation_ablation", "zero_downstream_baseline"]


def _centre(X, y):
    return X - X.mean(axis=0), y - y.mean()


def pooled_ols(X, y):
    """OLS slopes with an intercept (minimum-norm least squares when collinear)."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    if X.shape[1] == 0:
        return np.zeros(0)
    Xc, yc = _centre(X, y)
    return np.linalg.lstsq(Xc, yc, rcond=None)[0]


def ridge_to_pooled(X, y, b_pooled, lam):
    """Local slopes shrunk towards b_pooled; the intercept is unpenalised (data are centred).

    lam = 0 gives the local OLS fit, lam = inf gives b_pooled.
    """
    X, y, b_pooled = np.asarray(X, float), np.asarray(y, float), np.asarray(b_pooled, float)
    if X.shape[1] == 0:
        return np.zeros(0)
    if np.isinf(lam):
        return b_pooled.copy()
    if lam == 0:
        return pooled_ols(X, y)
    Xc, yc = _centre(X, y)
    n = len(y)
    A = Xc.T @ Xc + lam * n * np.eye(X.shape[1])
    return np.linalg.solve(A, Xc.T @ yc + lam * n * b_pooled)


def intercept(X, y, b):
    X, y = np.asarray(X, float), np.asarray(y, float)
    return float(y.mean() - (X.mean(axis=0) @ b if X.shape[1] else 0.0))


class PooledCoefficients:
    """OLS on a fixed set of cells (all NTC cells of a condition, a calibration fold or a cluster)."""

    def __init__(self, values, index):
        self.values, self.index, self.cache = values, index, {}

    def __call__(self, node, parents):
        key = (node, parents)
        if key not in self.cache:
            cols = [self.index[p] for p in parents]
            self.cache[key] = pooled_ols(self.values[:, cols], self.values[:, self.index[node]])
        return self.cache[key]


class LocalCoefficients:
    """Ridge-to-pooled coefficients on one cell's neighbourhood rows."""

    def __init__(self, values, index, rows, pooled, lam):
        self.values, self.index, self.rows, self.pooled, self.lam = values, index, rows, pooled, lam
        self.cache = {}

    def __call__(self, node, parents):
        key = (node, parents)
        if key not in self.cache:
            cols = [self.index[p] for p in parents]
            local = self.values[self.rows]
            self.cache[key] = ridge_to_pooled(local[:, cols], local[:, self.index[node]],
                                              self.pooled(node, parents), self.lam)
        return self.cache[key]


def unweighted_coefficients(node, parents):
    """Unweighted propagation ablation: b_jl = 1 / |pa(j)|, no fitting, no signs."""
    return np.full(len(parents), 1.0 / len(parents)) if parents else np.zeros(0)


def sorted_parents(dag, node, index):
    return tuple(sorted(dag.predecessors(node), key=index.__getitem__))


def knockdown_response(dag, source, coefficients, index):
    """Per-unit response (I - B)^-1 e_source by forward substitution over source's descendants."""
    delta = np.zeros(len(index))
    delta[index[source]] = 1.0
    descendants = nx.descendants(dag, source)
    if not descendants:
        return delta
    for node in nx.topological_sort(dag.subgraph(descendants | {source})):
        if node == source:
            continue
        parents = sorted_parents(dag, node, index)
        b = coefficients(node, parents)
        delta[index[node]] = float(b @ delta[[index[p] for p in parents]])
    return delta


def coefficient_matrix(dag, coefficients, index):
    """Dense B (B[j, l] = b_jl) for checks against forward substitution."""
    B = np.zeros((len(index), len(index)))
    for node in dag.nodes:
        parents = sorted_parents(dag, node, index)
        if parents:
            B[index[node], [index[p] for p in parents]] = coefficients(node, parents)
    return B


def cell_response(dags, source, coefficients, index):
    """Equal-weight mean over the local IDA choices; None when no valid choice exists."""
    if not dags:
        return None
    return np.mean([knockdown_response(d, source, coefficients, index) for d in dags], axis=0)


def zero_response(source, index):
    """Zero-downstream baseline: only the perturbed gene changes (per unit delta)."""
    delta = np.zeros(len(index))
    delta[index[source]] = 1.0
    return delta


def ida_dags(oriented, target, model):
    """DAG choices for a held-out target from an orient() result, with local IDA statistics."""
    if undirected_neighbours(oriented["pdag"], target):
        out = local_ida_choices(oriented["pdag"], target, model)
        return out["dags"], {k: out[k] for k in ("n_undirected", "n_candidates", "n_failed", "n_fallback")}
    dags = [oriented["dag"]] if oriented["valid"] else []
    return dags, dict(n_undirected=0, n_candidates=1, n_failed=int(not dags), n_fallback=int(oriented["fallback"]))


def graph_responses(pdag, seeds, targets, model, coefficients, index):
    """Per held-out target: (mean per-unit response over valid choices or None, IDA stats, DAG choices)."""
    oriented = orient(pdag, seeds, model)
    out = {}
    for target in targets:
        dags, stats = ida_dags(oriented, target, model)
        out[target] = (cell_response(dags, target, coefficients, index), stats, dags)
    return out


def group_responses(pdags, labels, values, seeds, targets, model, index):
    """Pooled-coefficient models (shared: one group; cluster: one per label).

    pdags: group label -> starting PDAG; coefficients are OLS on that group's cells with the group
    graph's own parent sets. Returns group label -> graph_responses output.
    """
    labels = np.asarray(labels)
    return {group: graph_responses(pdag, seeds, targets, model, PooledCoefficients(values[labels == group], index),
                                   index)
            for group, pdag in pdags.items()}


def node_prediction_errors(dag, x_cell, local_values, pooled, lams, index):
    """Squared errors of predicting each node with >= 1 parent from its parents, per lambda.

    local_values: the neighbourhood rows (training calibration fold); pooled: b_pooled provider on
    the training calibration fold. Returns (sum of squared errors per lambda, number of nodes).
    """
    sse, n_nodes = np.zeros(len(lams)), 0
    for node in dag.nodes:
        parents = sorted_parents(dag, node, index)
        if not parents:
            continue
        cols = [index[p] for p in parents]
        X, y = local_values[:, cols], local_values[:, index[node]]
        b_pooled = pooled(node, parents)
        for i, lam in enumerate(lams):
            b = ridge_to_pooled(X, y, b_pooled, lam)
            prediction = intercept(X, y, b) + x_cell[cols] @ b
            sse[i] += (x_cell[index[node]] - prediction) ** 2
        n_nodes += 1
    return sse, n_nodes


def delta_values(target_gene, target_folds, ntc_gene, ntc_folds, folds=(0, 1)):
    """delta(c) = mean of the perturbed gene in target cells of calib_fold c minus NTC cells of calib_fold c.

    Takes only the perturbed gene's values, so no other expression of target cells can enter.
    """
    target_gene, ntc_gene = np.asarray(target_gene, float), np.asarray(ntc_gene, float)
    target_folds, ntc_folds = np.asarray(target_folds), np.asarray(ntc_folds)
    if target_gene.ndim != 1 or ntc_gene.ndim != 1:
        raise ValueError("delta uses the perturbed gene's values only")
    out = {}
    for c in folds:
        t, n = target_gene[target_folds == c], ntc_gene[ntc_folds == c]
        if not len(t) or not len(n):
            raise ValueError(f"No cells in calib_fold {c}")
        out[c] = float(t.mean() - n.mean())
    return out


def target_delta(values, meta, target, node, condition):
    """delta per calib_fold for one held-out target: reads only the perturbed gene's column."""
    target_rows = (meta.condition.eq(condition) & meta.target.eq(target)).to_numpy()
    ntc_rows = (meta.condition.eq(condition) & meta.target.eq("NTC")).to_numpy()
    column = values[node].to_numpy()
    return delta_values(column[target_rows], meta.calib_fold.to_numpy()[target_rows],
                        column[ntc_rows], meta.calib_fold.to_numpy()[ntc_rows])


def fufa_weights(program_genes, controls, scale, nodes):
    """Panel weights w with predicted score shift = w . Delta x under the frozen linear scoring.

    score = (mean(program genes) - mean(background controls) - centre) / scale; genes outside the
    panel have a predicted change of 0 for every model.
    """
    position = {n: i for i, n in enumerate(nodes)}
    w = np.zeros(len(nodes))
    for g in program_genes:
        if g in position:
            w[position[g]] += 1.0 / len(program_genes)
    for g in controls:
        if g in position:
            w[position[g]] -= 1.0 / len(controls)
    return w / scale, dict(program_in_panel=sum(g in position for g in program_genes),
                           n_program=len(program_genes), controls_in_panel=sum(g in position for g in controls),
                           n_controls=len(controls))
