"""Distribution contrasts of fixed outcome scores: dispersion, Wasserstein and KS.

All statistics compare a target group with the NTC group of the same condition. Bootstrap draws
resample cells within each (condition, group) stratum with indices seeded by the stratum label,
so the same NTC draws are shared by every target, readout and contrast type.
"""

import zlib

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ks_2samp

CHUNK = 500


# ---------------------------------------------------------------- point statistics
def log_var_ratio(x, y):
    return float(np.log(np.var(x, ddof=1) / np.var(y, ddof=1)))


def iqr(values, axis=-1):
    q1, q3 = np.quantile(values, [0.25, 0.75], axis=axis)
    return q3 - q1


def iqr_diff(x, y):
    return float(iqr(x) - iqr(y))


def _grid(n, m):
    """Segments of [0, 1] on which both empirical quantile functions are constant."""
    u = np.union1d(np.arange(n + 1) / n, np.arange(m + 1) / m)
    mid = (u[:-1] + u[1:]) / 2
    return np.diff(u), np.floor(mid * n).astype(int), np.floor(mid * m).astype(int)


def wasserstein_sorted(xs, ys, grid=None):
    """Exact W1 between rows of sorted samples xs (B x n) and ys (B x m)."""
    xs, ys = np.atleast_2d(xs), np.atleast_2d(ys)
    width, ix, iy = grid or _grid(xs.shape[1], ys.shape[1])
    out = np.empty(len(xs))
    for s in range(0, len(xs), CHUNK):
        out[s:s + CHUNK] = np.abs(xs[s:s + CHUNK][:, ix] - ys[s:s + CHUNK][:, iy]) @ width
    return out


def ks_sorted(xs, ys):
    """Exact two-sample KS statistic between rows of xs (B x n) and ys (B x m), ties handled."""
    xs, ys = np.atleast_2d(xs), np.atleast_2d(ys)
    n, m = xs.shape[1], ys.shape[1]
    out = np.empty(len(xs))
    for s in range(0, len(xs), CHUNK):
        values = np.concatenate([xs[s:s + CHUNK], ys[s:s + CHUNK]], axis=1)
        order = np.argsort(values, axis=1, kind="stable")
        merged = np.take_along_axis(values, order, axis=1)
        cx = np.cumsum(order < n, axis=1)
        cy = np.arange(1, n + m + 1) - cx
        gap = np.abs(cx * m - cy * n) / (n * m)
        last = np.ones_like(merged, dtype=bool)
        last[:, :-1] = merged[:, 1:] != merged[:, :-1]
        out[s:s + CHUNK] = np.where(last, gap, 0).max(axis=1)
    return out


def wasserstein1(x, y):
    return float(wasserstein_sorted(np.sort(x), np.sort(y))[0])


def centered_wasserstein1(x, y):
    return wasserstein1(np.asarray(x) - np.mean(x), np.asarray(y) - np.mean(y))


def ks_stat(x, y):
    return float(ks_sorted(np.sort(x), np.sort(y))[0])


# ---------------------------------------------------------------- resampling
def stratum_rng(seed, *labels):
    """Generator seeded by the config seed and a stable hash of the stratum labels."""
    key = zlib.crc32("|".join(map(str, labels)).encode())
    return np.random.default_rng([seed, key])


def bootstrap_indices(n, repeats, seed, *labels):
    return stratum_rng(seed, "cell_bootstrap", *labels).integers(0, n, (repeats, n), dtype=np.int32)


class GroupDraws:
    """Bootstrap draws of one (condition, group) sample: sorted values and summary statistics."""

    def __init__(self, values, repeats, seed, condition, group):
        self.values = np.asarray(values, float)
        idx = bootstrap_indices(len(self.values), repeats, seed, condition, group)
        draws = self.values[idx]
        draws.sort(axis=1)
        self.sorted = draws
        self.mean = draws.mean(axis=1)
        self.var = draws.var(axis=1, ddof=1)
        self.iqr = iqr(draws, axis=1)

    @property
    def n(self):
        return len(self.values)


def bootstrap_p(samples):
    """Two-sided percentile bootstrap p-value for H0: statistic = 0."""
    samples = np.asarray(samples)
    samples = samples[np.isfinite(samples)]
    tail = min((samples <= 0).sum(), (samples >= 0).sum())
    return float(min(1.0, 2 * (1 + tail) / (len(samples) + 1)))


def permutation_p(x, y, n_perm, rng, center=False):
    """Label-permutation p-value for W1 (center=True: W1 of group-centered values, shape only)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if center:
        x, y = x - x.mean(), y - y.mean()
    observed = wasserstein1(x, y)
    pooled = np.concatenate([x, y])
    n, grid = len(x), _grid(len(x), len(y))
    exceed = 0
    for s in range(0, n_perm, CHUNK):
        block = rng.permuted(np.tile(pooled, (min(CHUNK, n_perm - s), 1)), axis=1)
        a, b = np.sort(block[:, :n], axis=1), np.sort(block[:, n:], axis=1)
        if center:
            a, b = a - a.mean(axis=1, keepdims=True), b - b.mean(axis=1, keepdims=True)
        exceed += int((wasserstein_sorted(a, b, grid) >= observed - 1e-12).sum())
    return float((1 + exceed) / (n_perm + 1))


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        q[ok] = false_discovery_control(p[ok], method="bh")
    return q


# ---------------------------------------------------------------- contrasts
def _row(statistic, kind, target, condition, readout, estimate, samples, p, method, n_t, n_r):
    samples = np.asarray(samples, float)
    finite = samples[np.isfinite(samples)]
    low, high = np.quantile(finite, [0.025, 0.975]) if len(finite) else (np.nan, np.nan)
    return dict(statistic=statistic, kind=kind, target=target, condition=condition, readout=readout,
                estimate=float(estimate), ci_low=float(low), ci_high=float(high), p_value=p, p_method=method,
                n_target=n_t, n_reference=n_r)


def dispersion_contrasts(scores, meta, targets, repeats, n_perm, seed, conditions=("DMSO", "Vem"),
                         progress=lambda _: None):
    """Target-vs-NTC distribution contrasts, drug effects and gene x drug interactions.

    scores: cells x readouts (frozen units); meta: condition, target (same index).
    Returns one row per (statistic, kind, target, condition, readout); BH within families of
    (statistic, kind, condition, readout) across targets.
    """
    meta = meta.reindex(scores.index)
    rows = []
    for readout in scores:
        values = scores[readout].to_numpy(float)
        draws = {}
        for condition in conditions:
            for group in ["NTC", *targets]:
                mask = (meta.condition.eq(condition) & meta.target.eq(group)).to_numpy()
                if mask.sum() < 2:
                    raise ValueError(f"Insufficient cells for {condition}/{group}")
                draws[condition, group] = GroupDraws(values[mask], repeats, seed, condition, group)
            ref = draws[condition, "NTC"]
            ref_centered = ref.sorted - ref.mean[:, None]
            for target in targets:
                d = draws[condition, target]
                x, y, n_t, n_r = d.values, ref.values, d.n, ref.n
                grid = _grid(n_t, n_r)
                rng = stratum_rng(seed, "permutation", readout, condition, target)
                lvr = np.log(d.var / ref.var)
                rows.append(_row("log_var_ratio", "gene", target, condition, readout, log_var_ratio(x, y), lvr,
                                 bootstrap_p(lvr), "bootstrap", n_t, n_r))
                diff = d.iqr - ref.iqr
                rows.append(_row("iqr_diff", "gene", target, condition, readout, iqr_diff(x, y), diff,
                                 bootstrap_p(diff), "bootstrap", n_t, n_r))
                w = wasserstein_sorted(d.sorted, ref.sorted, grid)
                rows.append(_row("wasserstein", "gene", target, condition, readout, wasserstein1(x, y), w,
                                 permutation_p(x, y, n_perm, rng), "permutation", n_t, n_r))
                wc = wasserstein_sorted(d.sorted - d.mean[:, None], ref_centered, grid)
                rows.append(_row("centered_wasserstein", "gene", target, condition, readout,
                                 centered_wasserstein1(x, y), wc,
                                 permutation_p(x, y, n_perm, rng, center=True), "permutation", n_t, n_r))
                rows.append(_row("ks", "gene", target, condition, readout, ks_stat(x, y),
                                 ks_sorted(d.sorted, ref.sorted), float(ks_2samp(x, y).pvalue), "ks_2samp",
                                 n_t, n_r))
            progress(f"dispersion {readout} / {condition}")
        if len(conditions) != 2:
            continue
        first, second = conditions
        for group in ["NTC", *targets]:
            a, b = draws[first, group], draws[second, group]
            label = f"{second}_minus_{first}"
            lvr = np.log(b.var / a.var)
            rows.append(_row("log_var_ratio", "drug", group, label, readout, log_var_ratio(b.values, a.values),
                             lvr, bootstrap_p(lvr), "bootstrap", b.n, a.n))
            diff = b.iqr - a.iqr
            rows.append(_row("iqr_diff", "drug", group, label, readout, iqr_diff(b.values, a.values), diff,
                             bootstrap_p(diff), "bootstrap", b.n, a.n))
            rows.append(_row("wasserstein", "drug", group, label, readout, wasserstein1(b.values, a.values),
                             wasserstein_sorted(b.sorted, a.sorted), np.nan, "descriptive", b.n, a.n))
            rows.append(_row("ks", "drug", group, label, readout, ks_stat(b.values, a.values),
                             ks_sorted(b.sorted, a.sorted), float(ks_2samp(b.values, a.values).pvalue),
                             "ks_2samp", b.n, a.n))
            if group == "NTC":
                continue
            r1, r2 = draws[first, "NTC"], draws[second, "NTC"]
            est = (log_var_ratio(b.values, a.values) - log_var_ratio(r2.values, r1.values))
            inter = np.log(b.var / a.var) - np.log(r2.var / r1.var)
            rows.append(_row("log_var_ratio", "interaction", group, label, readout, est, inter,
                             bootstrap_p(inter), "bootstrap", b.n + a.n, r1.n + r2.n))
            est = iqr_diff(b.values, a.values) - iqr_diff(r2.values, r1.values)
            inter = (b.iqr - a.iqr) - (r2.iqr - r1.iqr)
            rows.append(_row("iqr_diff", "interaction", group, label, readout, est, inter,
                             bootstrap_p(inter), "bootstrap", b.n + a.n, r1.n + r2.n))
        progress(f"dispersion {readout} drug / interaction")
    frame = pd.DataFrame(rows)
    frame["q_value"] = np.nan
    for _, idx in frame.groupby(["statistic", "kind", "condition", "readout"]).groups.items():
        frame.loc[idx, "q_value"] = bh(frame.loc[idx, "p_value"])
    return frame


# ---------------------------------------------------------------- depth sensitivity
def residualize_depth(score, log_umi, fit_mask):
    """Remove the linear log-UMI trend fitted in fit_mask cells; units are unchanged."""
    score, log_umi, fit_mask = np.asarray(score, float), np.asarray(log_umi, float), np.asarray(fit_mask, bool)
    xc = log_umi[fit_mask] - log_umi[fit_mask].mean()
    slope = float(xc @ (score[fit_mask] - score[fit_mask].mean()) / (xc @ xc))
    fitted = slope * xc
    resid = score[fit_mask] - score[fit_mask].mean() - fitted
    r2 = float(1 - resid.var() / score[fit_mask].var())
    adjusted = score - slope * (log_umi - log_umi[fit_mask].mean())
    return adjusted, slope, r2


# ---------------------------------------------------------------- well cluster bootstrap
def cluster_bootstrap(values, meta, targets, repeats, seed, conditions=("DMSO", "Vem")):
    """Resample whole wells (with replacement) within each condition; draws shared by all groups.

    values: 1-D scores aligned with meta (columns condition, target, well).
    Returns per-condition draws {condition: dict(mean=B x G, var=B x G, groups=[...], wells=n)}.
    """
    meta = meta.assign(value=np.asarray(values, float))
    groups = ["NTC", *targets]
    out = {}
    for condition in conditions:
        sub = meta[meta.condition.eq(condition) & meta.target.isin(groups)]
        wells = np.sort(meta.loc[meta.condition.eq(condition), "well"].unique())
        stats = sub.assign(sq=sub.value ** 2).groupby(["well", "target"]).agg(
            n=("value", "size"), s1=("value", "sum"), s2=("sq", "sum"))
        shape = (len(wells), len(groups))
        n, s1, s2 = (stats[c].unstack("target").reindex(index=wells, columns=groups).fillna(0).to_numpy()
                     .reshape(shape) for c in ("n", "s1", "s2"))
        rng = stratum_rng(seed, "cluster_bootstrap", condition)
        counts = rng.multinomial(len(wells), np.full(len(wells), 1 / len(wells)), size=repeats).astype(float)
        nn, t1, t2 = counts @ n, counts @ s1, counts @ s2
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(nn >= 1, t1 / nn, np.nan)
            var = np.where(nn >= 2, (t2 - t1 ** 2 / nn) / (nn - 1), np.nan)
        out[condition] = dict(mean=mean, var=var, n=nn, groups=groups, wells=len(wells))
    return out


def cluster_contrasts(draws, targets, readout, conditions=("DMSO", "Vem")):
    """KD - NTC mean, log variance ratio (per condition) and mean gene x drug interaction."""
    rows = []

    def add(statistic, kind, target, condition, samples):
        samples = np.asarray(samples, float)
        ok = np.isfinite(samples)
        low, high = np.quantile(samples[ok], [0.025, 0.975]) if ok.any() else (np.nan, np.nan)
        rows.append(dict(statistic=statistic, kind=kind, target=target, condition=condition, readout=readout,
                         cluster_ci_low=float(low), cluster_ci_high=float(high),
                         cluster_draws_dropped=int((~ok).sum())))

    for condition in conditions:
        d = draws[condition]
        for target in targets:
            j, r = d["groups"].index(target), 0
            add("mean", "gene", target, condition, d["mean"][:, j] - d["mean"][:, r])
            add("log_var_ratio", "gene", target, condition, np.log(d["var"][:, j] / d["var"][:, r]))
    first, second = conditions
    a, b = draws[first], draws[second]
    for target in targets:
        j = a["groups"].index(target)
        inter = (b["mean"][:, j] - a["mean"][:, j]) - (b["mean"][:, 0] - a["mean"][:, 0])
        add("mean", "interaction", target, f"{second}_minus_{first}", inter)
    return pd.DataFrame(rows)
