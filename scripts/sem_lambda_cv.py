"""NTC-only cross-validation of the ridge-to-pooled shrinkage lambda (one lambda per condition).

Usage:
    python scripts/sem_lambda_cv.py [--config config/sem_prediction.json] [--conditions DMSO,Vem]
                                    [--workers W] [--cells N]

For c in {0, 1}: every NTC cell of calib_fold 1-c uses its fold-specific completed DAG (stored PDAG ->
training-fold seeds -> Meek -> completion, for each of the 5 outer folds), a neighbourhood of the
round(0.35 n_c) nearest calib_fold-c NTC cells (same 20 PCs) and b_pooled from calib_fold c, and
predicts every node with >= 1 parent from its parents. The criterion is the mean squared error over
nodes and cells, averaged over both directions and the 5 DAG sets; lambda = argmin over the grid.
lambda = infinity (pure b_pooled) is reported as a reference row and cannot be selected. Seed
evidence uses training-fold targets only; no response enters the criterion. Writes
outputs/sem_prediction/lambda_cv/<condition>.json and results/summary/sem_prediction_lambda_cv.csv
(with --cells: outputs/sem_prediction_lambda_cv_smoke/ only).
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from aim1 import assert_no_cell_ids  # noqa: E402
from icpdag.cscn_networks import reference_model  # noqa: E402
from icpdag.sem_graphs import nearest_neighbours  # noqa: E402
from icpdag.sem_inputs import condition_inputs, fold_targets, load_panel_cells, target_nodes  # noqa: E402
from icpdag.sem_orientation import fold_seeds, orient  # noqa: E402
from icpdag.sem_prediction import PooledCoefficients, node_prediction_errors  # noqa: E402
from icpdag.target_split import load_split  # noqa: E402

_W = {}


def _init(state):
    _W.update(state)
    _W["model"] = reference_model(state["values"], state["nodes"], state["cscn"])
    _W["pooled"] = {c: PooledCoefficients(state["values"][state["folds"] == c], state["index"]) for c in (0, 1)}


def _run(cells):
    values, lams = _W["values"], _W["lams"]
    sse = np.zeros((len(_W["seeds"]), len(lams)))
    nodes = np.zeros(len(_W["seeds"]), dtype=np.int64)
    skipped = np.zeros(len(_W["seeds"]), dtype=np.int64)
    for k in cells:
        train = 1 - int(_W["folds"][k])  # the cell is evaluated; fit on the other calibration fold
        local = values[_W["neighbours"][train][k]]
        start = _W["graphs"].pdag(k)
        for f, seeds in enumerate(_W["seeds"]):
            oriented = orient(start, seeds, _W["model"])
            if not oriented["valid"]:
                skipped[f] += 1
                continue
            e, n = node_prediction_errors(oriented["dag"], values[k], local, _W["pooled"][train], lams, _W["index"])
            sse[f] += e
            nodes[f] += n
    return sse, nodes, skipped


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "sem_prediction.json"))
    parser.add_argument("--conditions")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--cells", type=int, help="Smoke test: seeded subset of evaluation cells")
    args = parser.parse_args(argv)
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    networks = json.loads((ROOT / cfg["inputs"]["networks"]).read_text(encoding="utf-8"))
    cscn = dict(window=networks["window"], ci_alpha=networks["ci_alpha"], max_conditioning=networks["max_conditioning"])
    grid = list(cfg["coefficients"]["lambda_grid"])
    lams = grid + [np.inf]
    split = load_split(ROOT / cfg["inputs"]["target_split"])
    folds = sorted(split.fold.unique())
    conditions = args.conditions.split(",") if args.conditions else cfg["conditions"]
    out_dir = ROOT / ("outputs/sem_prediction_lambda_cv_smoke" if args.cells else cfg["outputs_dir"] + "/lambda_cv")
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for condition in conditions:
        t0 = time.time()
        inp = condition_inputs(ROOT, cfg, condition)
        tnodes = target_nodes(inp["panel"])
        panel_values, panel_meta = load_panel_cells(ROOT / cfg["inputs"]["expression"], inp["ids"], [condition],
                                                    split.target)
        panel_values = pd.DataFrame(panel_values, index=panel_meta.index, columns=inp["nodes"])
        seeds = [fold_seeds(panel_values, panel_meta, fold_targets(split, f)[1], tnodes, condition)[0] for f in folds]
        cal = inp["meta"].calib_fold.to_numpy().astype(int)
        neighbours = {}
        for c in (0, 1):
            pool = np.flatnonzero(cal == c)
            neighbours[c] = nearest_neighbours(inp["pcs"], int(round(0.35 * len(pool))), candidates=pool)
        cells = np.arange(len(cal))
        if args.cells:
            cells = np.sort(np.random.default_rng(cfg["seed"]).choice(cells, args.cells, replace=False))
        state = dict(values=inp["values"], nodes=inp["nodes"], index=inp["index"], graphs=inp["graphs"],
                     neighbours=neighbours, folds=cal, seeds=seeds, lams=lams, cscn=cscn)
        sse = np.zeros((len(folds), len(lams)))
        n_nodes = np.zeros(len(folds), dtype=np.int64)
        skipped = np.zeros(len(folds), dtype=np.int64)
        chunks = [cells[i:i + 16].tolist() for i in range(0, len(cells), 16)]
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn"),
                                 initializer=_init, initargs=(state,)) as pool:
            for e, n, s in pool.map(_run, chunks):
                sse += e
                n_nodes += n
                skipped += s
        mse_by_fold = sse / n_nodes[:, None]
        mse = mse_by_fold.mean(axis=0)
        best = int(np.argmin(mse[:len(grid)]))
        chosen = grid[best]
        boundary = chosen in (grid[0], grid[-1])
        for i, lam in enumerate(lams):
            rows.append(dict(condition=condition, **{"lambda": "inf" if np.isinf(lam) else lam},
                             selectable=not np.isinf(lam), mse=mse[i],
                             **{f"mse_fold{f}": mse_by_fold[j, i] for j, f in enumerate(folds)},
                             selected=(not np.isinf(lam)) and lam == chosen, at_grid_boundary=boundary,
                             n_cells=len(cells), n_node_predictions=int(n_nodes.sum()),
                             cells_without_dag=int(skipped.sum())))
        (out_dir / f"{condition}.json").write_text(json.dumps(dict(
            condition=condition, lambda_selected=chosen, at_grid_boundary=boundary, grid=grid,
            mse={("inf" if np.isinf(l) else str(l)): float(m) for l, m in zip(lams, mse)},
            n_cells=len(cells), smoke_test=bool(args.cells)), indent=2), encoding="utf-8")
        print(f"{condition}: lambda {chosen}{' (grid boundary)' if boundary else ''}; MSE "
              f"{dict(zip(['inf' if np.isinf(l) else l for l in lams], np.round(mse, 5)))} ({time.time() - t0:.0f}s)")
    table = pd.DataFrame(rows)
    assert_no_cell_ids(table)
    path = (out_dir if args.cells else ROOT / "results" / "summary") / "sem_prediction_lambda_cv.csv"
    table.to_csv(path, index=False, lineterminator="\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
