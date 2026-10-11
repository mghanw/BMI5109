"""Split-half reliability and power simulation for the Aim 1 primary dispersion call.

The power simulation uses NTC cells only: a pseudo-target of n NTC cells is drawn without replacement,
the remaining NTC cells are its control, and mean-preserving heterogeneity is injected into the
pseudo-target. The test is the Aim 1 primary call (icpdag/dispersion.py): percentile-bootstrap p of
the log variance ratio, together with a mean-difference bootstrap CI that includes 0. Within a
replicate every scenario shares the split, the bootstrap draws and the mixture cell assignment
(common random numbers).
"""

import numpy as np
import pandas as pd

from icpdag import dispersion as ds


# ---------------------------------------------------------------- reliability
def split_halves(genes, rng):
    """Random partition of genes into floor(k/2) and ceil(k/2) genes (each half in input order)."""
    genes = list(genes)
    if len(genes) < 2:
        raise ValueError("At least two genes are needed for a split")
    order = rng.permutation(len(genes))
    half = len(genes) // 2
    return [genes[i] for i in sorted(order[:half])], [genes[i] for i in sorted(order[half:])]


def spearman_brown(r):
    """Reliability of the full score from the correlation of two half-scores."""
    r = np.asarray(r, float)
    return 2 * r / (1 + r)


# ---------------------------------------------------------------- heterogeneity injection
def inject_mixture(x, pi, d, order):
    """Shift round(pi*n) cells by +d and another round(pi*n) cells by -d; the mean is unchanged.

    order: a permutation of range(n); the first k cells go up, the next k go down, so assignments
    are nested across pi for a fixed order.
    """
    x = np.asarray(x, float)
    k = int(round(pi * len(x)))
    if 2 * k > len(x):
        raise ValueError("pi must be at most 0.5")
    out = x.copy()
    out[order[:k]] += d
    out[order[k:2 * k]] -= d
    return out


def inject_scale(x, s):
    """Scale deviations from the sample mean by s; the sample log variance ratio changes by 2 log s."""
    x = np.asarray(x, float)
    m = x.mean()
    return m + s * (x - m)


def mixture_true_lvr(pi, d_sd):
    """Expected log variance ratio of the symmetric mixture, d in units of the NTC SD."""
    return float(np.log1p(2 * pi * d_sd ** 2))


# ---------------------------------------------------------------- bootstrap
def bootstrap_counts(idx, n):
    """B x n multiplicity matrix of bootstrap index draws (rows of idx)."""
    idx = np.asarray(idx)
    repeats = len(idx)
    flat = (np.arange(repeats, dtype=np.int64)[:, None] * n + idx).ravel()
    return np.bincount(flat, minlength=repeats * n).reshape(repeats, n).astype(float)


def bootstrap_mean_var(values, counts):
    """Bootstrap means and ddof=1 variances from a multiplicity matrix.

    values: n or n x S (several samples sharing the same draws). Equal to GroupDraws.mean / .var for
    the same indices, up to floating-point error.
    """
    values = np.asarray(values, float)
    n = counts.shape[1]
    s1 = counts @ values
    s2 = counts @ values ** 2
    mean = s1 / n
    var = (s2 - n * mean ** 2) / (n - 1)
    return mean, var


def check_resolution(repeats, alphas):
    """A percentile-bootstrap p is at least 2/(B+1); every alpha must be reachable."""
    floor = 2 / (repeats + 1)
    bad = [a for a in alphas if floor >= a]
    if bad:
        raise ValueError(f"B={repeats} cannot reach alpha {bad}: minimum bootstrap p is {floor:.2e}")


# ---------------------------------------------------------------- simulation
def make_scenarios(pis, d_sds, scale_lvrs):
    rows = [dict(model="no_injection", pi=np.nan, d_sd=np.nan, s=1.0, true_lvr=0.0)]
    rows += [dict(model="mixture", pi=pi, d_sd=d, s=np.nan, true_lvr=mixture_true_lvr(pi, d))
             for pi in pis for d in d_sds]
    rows += [dict(model="scale", pi=np.nan, d_sd=np.nan, s=float(np.exp(v / 2)), true_lvr=float(v))
             for v in scale_lvrs]
    return rows


def simulate(ntc, n, scenarios, reps, repeats, seed, condition, alphas):
    """Per-replicate primary-test results for each scenario.

    ntc: NTC scores of one condition. Returns a long DataFrame (scenario x replicate) with the
    estimated log variance ratio, its bootstrap p and whether the mean-difference CI includes 0.
    """
    check_resolution(repeats, alphas)
    ntc = np.asarray(ntc, float)
    total = len(ntc)
    if not 2 <= n <= total - 2:
        raise ValueError("Pseudo-target size must leave at least two control cells")
    sd = np.std(ntc, ddof=1)
    rows = []
    for rep in range(reps):
        rng = ds.stratum_rng(seed, "power", condition, n, rep)
        perm = rng.permutation(total)
        x, y = ntc[perm[:n]], ntc[perm[n:]]
        order = rng.permutation(n)
        samples = []
        for sc in scenarios:
            if sc["model"] == "mixture":
                samples.append(inject_mixture(x, sc["pi"], sc["d_sd"] * sd, order))
            else:
                samples.append(inject_scale(x, sc["s"]))
        xs = np.column_stack(samples)
        tc = bootstrap_counts(ds.bootstrap_indices(n, repeats, seed, "power", condition, n, rep, "target"), n)
        cc = bootstrap_counts(ds.bootstrap_indices(total - n, repeats, seed, "power", condition, n, rep,
                                                   "control"), total - n)
        tm, tv = bootstrap_mean_var(xs, tc)
        cm, cv = bootstrap_mean_var(y, cc)
        lvr_boot = np.log(tv / cv[:, None])
        diff = tm - cm[:, None]
        low, high = np.quantile(diff, [0.025, 0.975], axis=0)
        var_y = np.var(y, ddof=1)
        for j, sc in enumerate(scenarios):
            rows.append(dict(scenario=j, replicate=rep, estimate=float(np.log(np.var(xs[:, j], ddof=1) / var_y)),
                             p_value=ds.bootstrap_p(lvr_boot[:, j]),
                             mean_ci_includes_0=bool(low[j] <= 0 <= high[j])))
    return pd.DataFrame(rows)


def summarize(results, scenarios, alphas):
    """Power per scenario: lvr test alone and the full primary call, with Monte Carlo SE."""
    out = []
    for j, sc in enumerate(scenarios):
        r = results[results.scenario.eq(j)]
        reps = len(r)
        row = dict(sc, reps=reps, mean_estimate=r.estimate.mean(), sd_estimate=r.estimate.std(ddof=1),
                   rate_mean_ci_includes_0=r.mean_ci_includes_0.mean())
        for label, alpha in alphas.items():
            for kind, hit in (("lvr", r.p_value.lt(alpha)), ("call", r.p_value.lt(alpha) & r.mean_ci_includes_0)):
                power = hit.mean()
                row[f"power_{kind}_{label}"] = power
                row[f"mcse_{kind}_{label}"] = np.sqrt(power * (1 - power) / reps)
        out.append(row)
    return pd.DataFrame(out)


def min_detectable(true_lvr, power, target=0.8):
    """Smallest true log variance ratio with power >= target, by linear interpolation; NaN if never."""
    true_lvr, power = np.asarray(true_lvr, float), np.asarray(power, float)
    order = np.argsort(true_lvr)
    x, p = true_lvr[order], power[order]
    hits = np.flatnonzero(p >= target)
    if not len(hits):
        return np.nan
    i = hits[0]
    if i == 0:
        return float(x[0])
    return float(x[i - 1] + (target - p[i - 1]) * (x[i] - x[i - 1]) / (p[i] - p[i - 1]))


def signal_scale_lvr(observed_lvr, reliability):
    """Log variance ratio of the reproducible signal implied by an observed one.

    Assumes measurement-noise variance is equal in target and NTC cells:
    exp(observed) = 1 + reliability * (R - 1). Undefined (NaN) when reliability <= 0.
    """
    observed_lvr, reliability = np.broadcast_arrays(np.asarray(observed_lvr, float), np.asarray(reliability, float))
    out = np.full(observed_lvr.shape, np.nan)
    ok = reliability > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        out[ok] = np.log1p(np.expm1(observed_lvr[ok]) / reliability[ok])
    return out[()] if out.ndim == 0 else out
