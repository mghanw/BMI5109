"""Marginal factorial effects with shared stratified bootstrap draws."""

import json

import numpy as np
import pandas as pd


def bootstrap_means(values, repeats, rng):
    out = np.empty((repeats, values.shape[1]))
    for start in range(0, repeats, 50):
        n = min(50, repeats - start)
        indices = rng.integers(0, len(values), (n, len(values)))
        out[start:start+n] = values[indices].mean(axis=1)
    return out


CONDITIONS = ("DMSO", "Vem")


def define_contrasts(targets, comparisons):
    """Encode unique group coefficients so shared controls cancel exactly."""
    rows = []

    def add(kind, target, comparator, condition, terms):
        rows.append(dict(contrast_id=f"{kind}__{target}__{comparator}__{condition}",
                         kind=kind, target=target, comparator=comparator,
                         condition=condition, terms=terms))

    for target in targets:
        for condition in CONDITIONS:
            add("gene", target, "NTC", condition,
                {(condition, target): 1, (condition, "NTC"): -1})
        add("interaction", target, "NTC", "Vem_minus_DMSO",
            {("Vem", target): 1, ("DMSO", target): -1,
             ("Vem", "NTC"): -1, ("DMSO", "NTC"): 1})
    for target in ["NTC", *targets]:
        add("drug", target, "DMSO", "Vem_minus_DMSO",
            {("Vem", target): 1, ("DMSO", target): -1})
    for left, right in comparisons:
        if left == right or left not in targets or right not in targets:
            raise ValueError("Pair comparisons require two distinct included targets")
        for condition in CONDITIONS:
            add("pair", left, right, condition,
                {(condition, left): 1, (condition, right): -1})
        add("pair_interaction", left, right, "Vem_minus_DMSO",
            {("Vem", left): 1, ("DMSO", left): -1,
             ("Vem", right): -1, ("DMSO", right): 1})
    return rows


def estimate_contrasts(scores, metadata, contrasts, repeats, seed, progress=lambda _: None):
    """Estimate contrasts for fixed outcome scores standardized to DMSO NTC SD."""
    if scores.index.has_duplicates or metadata.index.has_duplicates:
        raise ValueError("Cell identifiers must be unique")
    if set(scores.index) != set(metadata.index):
        raise ValueError("Scores and metadata must contain identical cells")
    metadata = metadata.reindex(scores.index)
    if not np.isfinite(scores.to_numpy(dtype=float)).all():
        raise ValueError("Outcomes must be finite")
    if any(c.startswith("match_") for c in scores):
        raise ValueError("Post-treatment matching variables are not ATE outcomes")
    if repeats < 100:
        raise ValueError("At least 100 bootstrap draws are required")
    rng = np.random.default_rng(seed)
    groups, draws, group_rows = {}, {}, []
    for key, ids in metadata.groupby(["condition", "target"], sort=True).groups.items():
        values = scores.loc[ids].to_numpy(dtype=float)
        if len(values) < 2:
            raise ValueError(f"Insufficient cells for group {key}")
        groups[key] = dict(mean=values.mean(axis=0), variance=values.var(axis=0, ddof=1), n=len(values))
        draws[key] = bootstrap_means(values, repeats, rng)
        for j, outcome in enumerate(scores):
            group_rows.append(dict(condition=key[0], target=key[1], outcome=outcome,
                                   n_cells=len(values), mean=groups[key]["mean"][j],
                                   sd=np.sqrt(groups[key]["variance"][j])))
        progress(f"ATE group means: {key[0]} / {key[1]}, n={len(values)}")
    rows = []
    for contrast in contrasts:
        terms = contrast["terms"]
        if set(terms) - set(groups):
            raise ValueError(f"Missing comparison groups: {set(terms) - set(groups)}")
        estimate = sum(weight * groups[key]["mean"] for key, weight in terms.items())
        samples = sum(weight * draws[key] for key, weight in terms.items())
        low, high = np.quantile(samples, [0.025, 0.975], axis=0)
        variance = sum(weight ** 2 * groups[key]["variance"] / groups[key]["n"]
                       for key, weight in terms.items())
        base = {k: v for k, v in contrast.items() if k != "terms"}
        counts = json.dumps({f"{c}/{t}": groups[c, t]["n"] for c, t in terms}, sort_keys=True)
        for j, outcome in enumerate(scores):
            rows.append(dict(**base, outcome=outcome, estimate=estimate[j], ci_low=low[j], ci_high=high[j],
                             iid_cell_se=np.sqrt(variance[j]), group_cell_counts=counts,
                             units="frozen_DMSO_NTC_SD", interval_scope="observed_pool_frozen_scores_cell_bootstrap",
                             causal_status="ATE_only_under_additional_identification_assumptions"))
    return pd.DataFrame(rows), pd.DataFrame(group_rows), groups


def reference_sensitivity(scores, metadata, contrasts, groups):
    """Change NTC composition, without treating constructs as screen replicates."""
    metadata = metadata.reindex(scores.index)
    per_guide, rows = {}, []
    for (condition, target, guide), ids in metadata.groupby(["condition", "target", "guide"]).groups.items():
        values = scores.loc[ids].to_numpy()
        per_guide[condition, target, guide] = dict(mean=values.mean(axis=0), n=len(ids))
        for j, outcome in enumerate(scores):
            rows.append(dict(condition=condition, target=target, guide=guide, outcome=outcome,
                             n_cells=len(ids), mean=values[:, j].mean()))
    guide_sets = [{g for c, t, g in per_guide if c == condition and t == "NTC"} for condition in CONDITIONS]
    if not guide_sets[0] or guide_sets[0] != guide_sets[1]:
        raise ValueError("Reference sensitivity needs the same observed NTC constructs in both conditions")
    guides = sorted(guide_sets[0])
    if len(guides) < 2:
        raise ValueError("At least two NTC constructs are needed for reference sensitivity")
    variants = [("equal_construct_weight", "all", guides, False)]
    variants += [("leave_one_out", guide, [g for g in guides if g != guide], True) for guide in guides]
    variants += [("single_construct", guide, [guide], True) for guide in guides]
    effects = []
    for scheme, changed, included, cell_weighted in variants:
        means = {key: val["mean"] for key, val in groups.items()}
        for condition in CONDITIONS:
            sub = [per_guide[condition, "NTC", guide] for guide in included]
            means[condition, "NTC"] = np.average([a["mean"] for a in sub], axis=0,
                                                 weights=[a["n"] for a in sub] if cell_weighted else None)
        for contrast in contrasts:
            if not any(target == "NTC" for _, target in contrast["terms"]):
                continue
            value = sum(w * means[key] for key, w in contrast["terms"].items())
            for j, outcome in enumerate(scores):
                effects.append(dict(contrast_id=contrast["contrast_id"], outcome=outcome,
                                    reference_scheme=scheme, changed_construct=changed,
                                    estimate=value[j], n_reference_constructs=len(included)))
    frame = pd.DataFrame(effects)
    summary = frame.groupby(["contrast_id", "outcome", "reference_scheme"]).estimate.agg(
        minimum="min", maximum="max").reset_index()
    return frame, summary, pd.DataFrame(rows)
