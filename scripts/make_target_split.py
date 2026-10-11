"""Freeze the held-out target split for Phase 2 and later phases.

Usage:
    python scripts/make_target_split.py [--qc results/summary/phase1_target_qc.csv]
                                        [--out config/target_split.json] [--folds 5] [--seed 20261011]

Inputs are Phase 1 eligibility and knockdown remaining ratios only (the perturbed gene's own
expression). No response statistic is read.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from icpdag.target_split import make_split  # noqa: E402

QC_COLUMNS = ["condition", "target", "eligible", "knockdown_remaining_ratio"]


def build_spec(qc_path, n_folds, seed):
    qc = pd.read_csv(qc_path, usecols=QC_COLUMNS)
    split = make_split(qc, n_folds, seed)
    digest = hashlib.sha256(Path(qc_path).read_bytes()).hexdigest()
    return {
        "description": "Held-out target split over eligible target identities, frozen before any "
                       "response-based step. No response information was used.",
        "source_table": Path(qc_path).as_posix(),
        "source_sha256": digest,
        "assignment_rule": "paper_gt5_max0.5_sec0.3",
        "eligibility": "> 50 singlets and adequate knockdown (remaining ratio <= 0.5, BH q < 0.05) in both conditions",
        "n_folds": n_folds,
        "seed": seed,
        "stratification": "kd_strength = mean of DMSO and Vem knockdown remaining ratios; tertiles by rank "
                          "(0 = strongest knockdown); within each tertile targets are permuted with "
                          "numpy default_rng(seed) and dealt round-robin to folds with a running counter",
        "n_targets": len(split),
        "fold_sizes": {str(k): int(v) for k, v in split.fold.value_counts().sort_index().items()},
        "targets": {r.target: {"fold": int(r.fold), "kd_tertile": int(r.kd_tertile),
                               "kd_strength": round(float(r.kd_strength), 6)} for r in split.itertuples()},
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--qc", default="results/summary/phase1_target_qc.csv")
    parser.add_argument("--out", default="config/target_split.json")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20261011)
    args = parser.parse_args(argv)
    spec = build_spec(args.qc, args.folds, args.seed)
    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(spec, f, indent=2)
        f.write("\n")
    print(f"wrote {args.out}: {spec['n_targets']} targets, fold sizes {spec['fold_sizes']}")


if __name__ == "__main__":
    main()
