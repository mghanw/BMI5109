"""Expression-matched program scores on a frozen DMSO NTC scale."""

import numpy as np
import pandas as pd


def program_scores(log_expr, symbols, annotation, meta, programs, cfg):
    lookup = {g: i for i, g in enumerate(symbols)}
    ntc = (meta.target.eq("NTC") & meta.condition.eq("DMSO")).to_numpy()
    mean = np.asarray(log_expr[ntc].mean(axis=0)).ravel()
    detected = np.asarray((log_expr[ntc] > 0).mean(axis=0)).ravel()
    coding = set(annotation.loc[annotation.gene_type.eq("protein_coding"), "symbol"])
    forbidden = set(cfg["network_genes"]) | set(cfg["gene_readouts"])
    for info in programs.values():
        forbidden.update(info["genes"])
    candidates = np.array([i for i, g in enumerate(symbols)
                           if g in coding and g not in forbidden and detected[i] >= 0.01])
    order = candidates[np.argsort(mean[candidates], kind="stable")]
    bins = np.array_split(order, 24)
    centers = np.array([mean[b].mean() for b in bins])
    scores, raw_scores, audit = {}, {}, []
    rng = np.random.default_rng(cfg["seed"])
    for name, info in programs.items():
        genes = [g for g in info["genes"] if g in lookup]
        if len(genes) < max(1, 0.8 * len(info["genes"])):
            raise ValueError(f"Less than 80% gene mapping coverage for {name}")
        indices = [lookup[g] for g in genes]
        controls = set()
        for i in indices:
            pool = bins[np.argmin(np.abs(centers - mean[i]))]
            controls.update(rng.choice(pool, min(50, len(pool)), replace=False).tolist())
        controls = sorted(controls)
        raw = np.asarray(log_expr[:, indices].mean(axis=1)).ravel()
        if info.get("role") == "matching":
            # Matching genes must be disjoint from readout genes, including
            # background-control genes used by readout scores.
            controls = []
            background = np.zeros(len(meta))
        else:
            background = np.asarray(log_expr[:, controls].mean(axis=1)).ravel()
        value = raw - background
        center, scale = value[ntc].mean(), value[ntc].std(ddof=1)
        raw_center, raw_scale = raw[ntc].mean(), raw[ntc].std(ddof=1)
        if min(scale, raw_scale) <= 1e-10:
            raise ValueError(f"Unidentifiable score in DMSO NTC: {name}")
        scores[name] = (value - center) / scale
        raw_scores[name] = (raw - raw_center) / raw_scale
        audit.append(dict(program=name, **info, mapped_genes=genes,
                          missing_genes=sorted(set(info["genes"]) - set(genes)),
                          n_requested=len(info["genes"]), n_mapped=len(genes),
                          detectable_DMSO_NTC=int(sum(detected[i] > 0 for i in indices)),
                          controls=[symbols[i] for i in controls], center=float(center), scale=float(scale),
                          raw_center=float(raw_center), raw_scale=float(raw_scale)))
    return pd.DataFrame(scores, index=meta.index), pd.DataFrame(raw_scores, index=meta.index), audit
