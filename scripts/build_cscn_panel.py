"""Build (or verify) the network gene panel from NTC cells.

Usage:
    python scripts/build_cscn_panel.py [--config config/cscn_panel.json] [--freeze-size P]

Inputs: NTC cell expression of rna_singlets.h5ad (both conditions), the identities of the eligible
targets (config/target_split.json) and the mapped Fufa program genes (config/programs.json). No
expression of target-perturbed cells, no response and no Aim 1 table is read.

Writes the forced genes, the full ranked fill-gene list and the input hashes into the config, and
results/summary/cscn_panel_genes.csv (the frozen panel, or the largest candidate panel before
freezing). --freeze-size P fixes the panel size. A frozen config is only verified: the panel is
recomputed and any difference raises; the file is never overwritten.
"""
import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from icpdag.cscn_panel import build_panel, forced_genes  # noqa: E402
from icpdag.target_split import load_split  # noqa: E402

COMPARED = ["inputs", "n_ntc_cells", "n_fill_universe", "forced_genes", "fill_genes_ranked"]


def sha256(path, text=False):
    data = Path(path).read_bytes()
    if text:  # LF-normalised so the digest does not depend on the checkout's line endings
        data = data.replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def load_inputs(params):
    import anndata as ad

    adata = ad.read_h5ad(ROOT / params["expression"])
    obs = adata.obs
    meta = pd.DataFrame({"condition": obs.condition.astype(str).to_numpy(),
                         "target": obs.target.astype(str).to_numpy(),
                         "total_counts": obs.total_counts.to_numpy(float)})
    targets = sorted(load_split(ROOT / params["target_split"]).target)
    # Identity metadata only: which gene each eligible target's guides perturb
    ids = obs.loc[obs.target.astype(str).isin(targets), ["target", "target_gene_id"]].astype(str).drop_duplicates()
    if ids.target.duplicated().any() or set(ids.target) != set(targets):
        raise ValueError("Each eligible target needs exactly one target_gene_id")
    with open(ROOT / params["programs"], encoding="utf-8") as f:
        fufa = json.load(f)["programs"][params["fufa_program"]]["genes"]
    forced = forced_genes(targets, dict(zip(ids.target, ids.target_gene_id)), fufa, adata.var)
    return adata.X, meta, adata.var, forced


def build_spec(params):
    counts, meta, var, forced = load_inputs(params)
    result = build_panel(counts, meta, var, forced, params)
    f = result["forced"]
    spec = {
        "description": "Network gene panel for cell-specific CSCN graphs. Selected from NTC cells only "
                       "(both conditions); target identities and the Fufa gene list are the only other "
                       "inputs. A panel of size P is the forced genes followed by the top P - "
                       "n_forced fill genes.",
        "parameters": params,
        "method": {
            "expression": "log1p CP10k with whole-transcriptome total_counts (as scripts/aim1.py)",
            "fill_universe": "scoring genes (one per symbol, no PAR_Y) detected (count > 0) in >= "
                             "fill_min_detection of NTC cells in each condition, excluding exclude_pattern "
                             "and forced genes",
            "dispersion": "Seurat-style mean.var.plot standardization, computed on log1p CP10k with "
                          "equal-frequency mean bins: per condition, d = log(variance / mean) over NTC "
                          "cells; z-score of d within n_bins mean-quantile bins of the fill universe",
            "ranking": "mean of the two within-condition z ranks (1 = most dispersed); ties by mean z, "
                       "then gene_id",
            "low_detection_flag": "forced gene detected in < forced_low_detection of NTC cells in both "
                                  "conditions; forced genes are kept regardless",
        },
        "inputs": {
            "expression_sha256": sha256(ROOT / params["expression"]),
            "programs_sha256": sha256(ROOT / params["programs"], text=True),
            "target_split_sha256": sha256(ROOT / params["target_split"], text=True),
        },
        "n_ntc_cells": result["n_ntc"],
        "n_fill_universe": result["n_universe"],
        "n_forced": len(f),
        "forced_genes": [dict(gene_id=r.gene_id, symbol=r.symbol, role=r.role, target=r.target,
                              in_fufa=bool(r.in_fufa), low_detection=bool(r.low_detection))
                         for r in f.itertuples()],
        "fill_genes_ranked": [[r.gene_id, r.symbol] for r in result["ranked"].itertuples()],
    }
    return spec, result


def dump(spec, path):
    text = json.dumps(spec, indent=2)
    text = re.sub(r'\[\s+("[^"\n]*"),\s+("[^"\n]*")\s+\]', r"[\1, \2]", text)  # one line per fill gene
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text + "\n")


def genes_table(result, spec, size, conditions):
    forced = result["forced"].assign(fill_rank=pd.NA)
    n_fill = size - len(forced)
    fill = result["ranked"].head(n_fill).assign(role="variable", target="", in_fufa=False, low_detection=False)
    table = pd.concat([forced, fill], ignore_index=True)
    sizes = sorted(spec["parameters"]["candidate_sizes"])
    table["smallest_candidate_size"] = [
        sizes[0] if pd.isna(r) else next((s for s in sizes if s - len(forced) >= r), pd.NA)
        for r in table.fill_rank]
    cols = ["gene_id", "symbol", "role", "target", "in_fufa", "fill_rank", "smallest_candidate_size"]
    for c in conditions:
        cols += [f"detection_{c}", f"mean_{c}", f"dispersion_z_{c}"]
    table = table[cols + ["low_detection"]]
    table["fill_rank"] = table.fill_rank.astype("Int64")
    return table


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "cscn_panel.json"))
    parser.add_argument("--freeze-size", type=int)
    parser.add_argument("--summary", default=str(ROOT / "results" / "summary" / "cscn_panel_genes.csv"))
    args = parser.parse_args(argv)
    with open(args.config, encoding="utf-8") as f:
        existing = json.load(f)
    params = existing["parameters"]
    spec, result = build_spec(params)
    if existing.get("frozen"):
        diff = [k for k in COMPARED if existing.get(k) != json.loads(json.dumps(spec[k]))]
        if diff:
            raise ValueError(f"Frozen panel differs from the recomputed panel in: {diff}")
        if args.freeze_size and args.freeze_size != existing["panel_size"]:
            raise ValueError(f"Panel already frozen at size {existing['panel_size']}")
        print(f"verified frozen panel (size {existing['panel_size']}): identical")
        size = existing["panel_size"]
    else:
        size = args.freeze_size
        if size is not None:
            if size not in params["candidate_sizes"]:
                raise ValueError(f"Panel size must be one of {params['candidate_sizes']}")
            if "forced_genes" in existing:
                diff = [k for k in COMPARED if existing.get(k) != json.loads(json.dumps(spec[k]))]
                if diff:
                    raise ValueError(f"Recomputed panel differs from the provisional config in: {diff}")
        spec.update(panel_size=size, frozen=size is not None)
        dump(spec, args.config)
        print(f"wrote {args.config} (frozen: {size is not None})")
    table = genes_table(result, spec, size or max(params["candidate_sizes"]), params["conditions"])
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.summary, index=False, lineterminator="\n")
    f = result["forced"]
    print(f"NTC cells: {result['n_ntc']}; forced genes: {len(f)} "
          f"({f.role.value_counts().to_dict()}, in Fufa: {int(f.in_fufa.sum())})")
    print(f"detected >= {params['fill_min_detection']:.0%} in both conditions: {result['n_detected']}; "
          f"excluded MT/ribosomal: {result['n_excluded_detected']}; fill universe: {result['n_universe']}")
    low = f[f.low_detection]
    print(f"low-detection forced genes: {len(low)} by role {low.role.value_counts().to_dict()}: "
          f"{', '.join(low.symbol)}")
    print(f"wrote {args.summary} ({len(table)} genes)")


if __name__ == "__main__":
    main()
