"""Comparator graphs for the SEM models: density-matched population PC, clusters, local association.

The shared and cluster graphs reuse the pinned CSCN PC-stable skeleton search and unshielded-collider
orientation; a subclass replaces only the CI test with a population Spearman partial-correlation
test (Fisher z). The threshold is chosen on structure only: the sparsest graph in which at least as
many targets have an incident edge as have an edge in >= 10% of the same cells' CSCN graphs
(target-coverage matching). Limitation: the test family differs from the cell-specific model (local
CSCN count test vs population partial correlation).
"""

from contextlib import redirect_stdout
from io import StringIO
import time

import numpy as np
from scipy.stats import norm, rankdata
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

from ._vendor.cscn.core import CSCN


def spearman_matrix(values):
    """Rank correlation matrix; columns without variance get NaN rows and columns."""
    ranks = np.apply_along_axis(rankdata, 0, np.asarray(values, dtype=float))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(ranks, rowvar=False)
    constant = np.asarray(values).std(axis=0) == 0
    corr[constant, :] = np.nan
    corr[:, constant] = np.nan
    return corr, np.flatnonzero(constant)


class PopulationPC(CSCN):
    """Vendored PC-stable search with a population partial-correlation CI test.

    Independence iff |atanh(r_xy.Z)| * sqrt(n - |Z| - 3) < z_crit, with r_xy.Z from the inverse of
    the rank-correlation submatrix of {x, y} and Z (pseudo-inverse if singular). Undefined
    correlations count as independent. z statistics are cached across thresholds.
    """

    def __init__(self, nodes, corr, n_obs, z_crit, max_conditioning, cache=None, max_seconds=None):
        super().__init__(variable_names=list(nodes), significance_level=0.5, max_cond_vars=max_conditioning,
                         query_engine="bitmap", pc_var="stable")
        self.corr, self.n_obs, self.z_crit = corr, int(n_obs), float(z_crit)
        self.z_cache = {} if cache is None else cache
        self.deadline = None if max_seconds is None else time.monotonic() + max_seconds
        self.run_core(np.zeros((2, len(nodes)), dtype=np.float32), use_nmf=False)

    def z_statistic(self, i, j, given):
        key = (i, j, given)
        if key not in self.z_cache:
            index = [i, j, *given]
            sub = self.corr[np.ix_(index, index)]
            if not np.isfinite(sub).all():
                self.z_cache[key] = 0.0
            else:
                precision = np.linalg.pinv(sub)
                denominator = np.sqrt(precision[0, 0] * precision[1, 1])
                r = -precision[0, 1] / denominator if denominator > 0 else 0.0
                r = float(np.clip(r, -0.999999999, 0.999999999))
                self.z_cache[key] = abs(np.arctanh(r)) * np.sqrt(max(self.n_obs - len(given) - 3, 1))
        return self.z_cache[key]

    def conditional_independence_test(self, X, Y, Z, data=None, independencies=None, key_cell_idx=0,
                                      significance_level=None, sigmoid_score=None):
        i, j = sorted((self._coerce_gene_id(X), self._coerce_gene_id(Y)))
        given = tuple(sorted(self._coerce_gene_id(z) for z in Z))
        if self.deadline is not None and time.monotonic() > self.deadline:
            raise RuntimeError("Population PC fit exceeded the runtime guard")
        return bool(self.z_statistic(i, j, given) < self.z_crit)

    def fit(self):
        diagnostics = StringIO()
        with redirect_stdout(diagnostics):
            skeleton, separators = self._estimate_pc_skeleton(0)
        if "ICT ERROR:" in diagnostics.getvalue():
            raise RuntimeError(diagnostics.getvalue())
        return skeleton, self._orient_unshielded_colliders(skeleton, separators)


def coverage_matched_graph(values, nodes, target_nodes, min_targets, max_conditioning, max_seconds=600,
                           iterations=40):
    """Sparsest population PC graph in which >= min_targets of target_nodes have an incident edge.

    Searches z_crit (larger = sparser): starts at the level-0 bound (the largest threshold at which
    min_targets targets still have a marginal edge), lowers it by factors of 0.9 until the PC skeleton
    covers enough targets, then bisects; returns the largest evaluated z_crit with enough coverage.
    Runtime guard: raises if any single PC fit takes longer than max_seconds. Structure only.
    """
    values = np.asarray(values, dtype=float)
    corr, constant = spearman_matrix(values)
    n = len(values)
    index = {name: i for i, name in enumerate(nodes)}
    targets = [t for t in target_nodes if t in index]
    cache, results = {}, {}
    with np.errstate(invalid="ignore"):
        z = np.abs(np.arctanh(np.clip(np.nan_to_num(corr), -0.999999999, 0.999999999))) * np.sqrt(n - 3)
    np.fill_diagonal(z, 0.0)
    upper = z[np.triu_indices(len(nodes), 1)]

    def evaluate(threshold):
        if threshold not in results:
            skeleton, pdag = PopulationPC(nodes, corr, n, threshold, max_conditioning, cache, max_seconds).fit()
            covered = sum(skeleton.degree(t) > 0 for t in targets)
            results[threshold] = (covered, skeleton, pdag)
        return results[threshold][0]

    if min_targets <= 0:
        threshold = float(upper.max()) + 1.0
        evaluate(threshold)
    else:
        if min_targets > len(targets):
            raise ValueError("More targets requested than available")
        best_marginal = np.sort(z[[index[t] for t in targets]].max(axis=1))[::-1]
        high = low = float(best_marginal[min_targets - 1])
        while evaluate(low) < min_targets:
            high = low
            low *= 0.9
            if low < 1e-6:
                raise RuntimeError(f"Coverage of {min_targets} targets is not reachable")
        for _ in range(iterations):
            if high - low < 1e-6:
                break
            mid = 0.5 * (low + high)
            if evaluate(mid) >= min_targets:
                low = mid
            else:
                high = mid
        threshold = max(z_ for z_, r in results.items() if r[0] >= min_targets)
    covered, skeleton, pdag = results[threshold]
    return dict(skeleton=skeleton, pdag=pdag, z_crit=threshold, edges=skeleton.number_of_edges(),
                targets_with_edge=covered, min_targets=int(min_targets), n_obs=n,
                log10_alpha=float((np.log(2) + norm.logsf(threshold)) / np.log(10)),
                zero_variance=[nodes[i] for i in constant], iterations=len(results),
                level0_edges=int((upper >= threshold).sum()))


def principal_components(values, n_components, seed):
    """Centred (not scaled) PCA scores, full SVD."""
    return PCA(n_components=n_components, svd_solver="full", random_state=seed).fit_transform(
        np.asarray(values, dtype=float))


def nearest_neighbours(scores, size, candidates=None):
    """Indices of the `size` nearest cells (Euclidean in scores) for every cell, itself included.

    candidates: optional positions the neighbours are drawn from (default all cells). Stable order.
    """
    scores = np.asarray(scores, dtype=float)
    pool = np.arange(len(scores)) if candidates is None else np.asarray(candidates)
    if size > len(pool):
        raise ValueError("Neighbourhood larger than the candidate pool")
    out = np.empty((len(scores), size), dtype=np.int32)
    for start in range(0, len(scores), 256):
        block = scores[start:start + 256]
        d = ((block[:, None, :] - scores[pool][None, :, :]) ** 2).sum(axis=2)
        order = np.argsort(d, axis=1, kind="stable")[:, :size]
        out[start:start + len(block)] = pool[order]
    return out


def clusters(scores, k, seed):
    return KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(np.asarray(scores, dtype=float))


def association_pairs(values, n_edges):
    """Top n_edges gene pairs by absolute Pearson correlation; ties by node order; constant genes excluded."""
    values = np.asarray(values, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(values, rowvar=False)
    i, j = np.triu_indices(values.shape[1], 1)
    strength = np.abs(corr[i, j])
    keep = np.isfinite(strength)
    i, j, strength = i[keep], j[keep], strength[keep]
    order = np.lexsort((j, i, -strength))[:n_edges]
    return list(zip(i[order].tolist(), j[order].tolist()))
