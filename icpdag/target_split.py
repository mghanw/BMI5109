"""Held-out target split frozen before any response statistic is computed.

The split uses only target identities, eligibility and knockdown strength. Knockdown strength is
the remaining expression of the perturbed gene itself, so it carries no response information.
"""

import re

import numpy as np
import pandas as pd

RESPONSE_PATTERN = re.compile(r"score|estimate|program|fufa|kanai|tsoi|response|contrast", re.I)


def eligible_targets(qc):
    """Targets eligible in every condition of the Phase 1 QC table."""
    rows = qc[qc.target.ne("NTC")].assign(eligible=lambda d: d.eligible.eq(True))
    flags = rows.pivot(index="target", columns="condition", values="eligible")
    flags = flags.fillna(False).astype(bool)
    return sorted(flags.index[flags.all(axis=1)])


def make_split(qc, n_folds, seed):
    """Seeded fold assignment, stratified by knockdown strength tertile.

    qc: Phase 1 target QC with columns condition, target, eligible, knockdown_remaining_ratio.
    Returns one row per eligible target: target, kd_strength, kd_tertile (0 = strongest), fold.
    """
    bad = [c for c in qc.columns if RESPONSE_PATTERN.search(c)]
    if bad:
        raise ValueError(f"Response-like columns are not allowed in the split input: {bad}")
    targets = eligible_targets(qc)
    if len(targets) < n_folds:
        raise ValueError("Fewer eligible targets than folds")
    ratio = qc[qc.target.isin(targets)].pivot(index="target", columns="condition",
                                               values="knockdown_remaining_ratio")
    if ratio.isna().any().any():
        raise ValueError("Missing knockdown ratio for an eligible target")
    strength = ratio.mean(axis=1).rename("kd_strength").reset_index()
    strength = strength.sort_values(["kd_strength", "target"], kind="stable").reset_index(drop=True)
    tertile = np.empty(len(strength), int)
    for k, block in enumerate(np.array_split(np.arange(len(strength)), 3)):
        tertile[block] = k
    strength["kd_tertile"] = tertile
    rng = np.random.default_rng(seed)
    fold, counter = {}, 0
    for k in range(3):
        members = strength.loc[strength.kd_tertile.eq(k), "target"].to_numpy()
        for target in rng.permutation(members):
            fold[target] = counter % n_folds
            counter += 1
    strength["fold"] = strength.target.map(fold)
    return strength.sort_values("target").reset_index(drop=True)[["target", "kd_strength", "kd_tertile", "fold"]]


def load_split(path):
    """Target -> (fold, kd_tertile) table from config/target_split.json."""
    import json

    with open(path, encoding="utf-8") as f:
        spec = json.load(f)
    return pd.DataFrame([dict(target=t, **v) for t, v in spec["targets"].items()])


def add_fold(frame, split, column="target"):
    """Attach fold and kd_tertile; NTC and non-eligible rows get NA."""
    lookup = split.set_index("target")[["fold", "kd_tertile"]]
    out = frame.drop(columns=[c for c in ("fold", "kd_tertile") if c in frame], errors="ignore")
    out = out.join(lookup, on=column)
    out["fold"] = out.fold.astype("Int64")
    out["kd_tertile"] = out.kd_tertile.astype("Int64")
    return out
