"""Timing smoke test for cell-specific CSCN graphs over panel size x maximum conditioning order.

Usage:
    python scripts/cscn_panel_timing.py [--config config/cscn_networks.json] [--workers W]
                                        [--sizes 150,200] [--orders 2,3] [--focal 20] [--run-name NAME]
    Window x CI alpha sensitivity at one panel size and order:
    python scripts/cscn_panel_timing.py --sizes 250 --orders 2 --windows 0.35,0.5,0.65 --alphas 0.05,0.1

Reference set = all NTC cells of the condition; focal cells = the same seeded NTC cells per
condition in every configuration. Configurations run cheapest first (panel size, then order); a
configuration is skipped when a cheaper one in both dimensions exceeded the threshold (size x order
grid only). Within a configuration, after each round of completed cells, pending cells are cancelled
once the running median exceeds the threshold, and workers are stopped when no cell finishes within
the hard timeout. Writes outputs/<run>/ (config.json, log, per-cell timings, grid, per-target edge
shares); the default full size x order grid is also copied to results/summary/cscn_panel_timing.csv.
Hyperparameters are compared on graph diagnostics and runtime only, never prediction performance.
"""
import argparse
from concurrent.futures import FIRST_COMPLETED, wait
import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from icpdag.cscn_networks import peak_rss_mb, run_chunk, worker_pool  # noqa: E402
from icpdag.cscn_panel import load_ntc_expression  # noqa: E402
from infer_cscn_networks import focal_subset, panel_nodes  # noqa: E402

LOG = logging.getLogger("cscn_timing")
KEYS = ["panel_size", "max_conditioning", "window", "ci_alpha"]


def stop_workers(pool):
    for process in list(getattr(pool, "_processes", {}).values()):
        process.terminate()
    pool.shutdown(wait=False, cancel_futures=True)


def time_configuration(references, focal, nodes, n_targets, cfg, workers, threshold, hard_timeout):
    tasks = [dict(condition=c, locals=[int(i)], strict=True) for pair in zip(*[
        [(c, i) for i in focal[c]] for c in references]) for c, i in pair]
    results, status, start = [], "ok", time.perf_counter()
    pool = worker_pool(references, nodes, cfg, workers)
    try:
        futures = {pool.submit(run_chunk, t): t for t in tasks}
        pending, last = set(futures), time.perf_counter()
        while pending:
            done, pending = wait(pending, timeout=5, return_when=FIRST_COMPLETED)
            for future in done:
                if future.cancelled():
                    continue
                result = future.result()
                result["condition"] = futures[future]["condition"]
                results.append(result)
                last = time.perf_counter()
                if len(results) % workers == 0:
                    median = np.median([r["rows"][0]["seconds"] for r in results])
                    if median > threshold and status == "ok":
                        status = "exceeded"
                        for f in pending:
                            f.cancel()
            pending = {f for f in pending if not f.cancelled()}
            if pending and time.perf_counter() - last > hard_timeout:
                status = "timeout"
                stop_workers(pool)
                break
    finally:
        if status != "timeout":
            pool.shutdown(wait=True)
    seconds = np.array([r["rows"][0]["seconds"] for r in results])
    if status == "ok" and len(seconds) and np.median(seconds) > threshold:
        status = "exceeded"
    rows = []
    for r in results:
        row = dict(r["rows"][0], condition=r["condition"], pid=r["pid"], setup_seconds=r["setup_seconds"],
                   peak_rss_mb=r["peak_rss_mb"])
        with_edge = sorted((set(r["edges"]["source"].tolist()) | set(r["edges"]["destination"].tolist()))
                           & set(range(n_targets)))
        row["target_fraction"] = len(with_edge) / n_targets
        row["targets_with_edge"] = ";".join(map(str, with_edge))
        rows.append(row)
    return pd.DataFrame(rows), status, time.perf_counter() - start


def summarize(cells, status, key, workers, n_total, wall):
    base = dict(key, status=status, n_cells_timed=len(cells))
    if not len(cells):
        return dict(base, config_wall_seconds=wall)
    s = cells.seconds
    setup = cells.drop_duplicates(["pid", "condition"]).setup_seconds
    worker_peak = cells.groupby("pid").peak_rss_mb.max()
    mean = s.mean()
    return dict(base, seconds_mean=mean, seconds_median=s.median(), seconds_p90=s.quantile(0.9), seconds_max=s.max(),
                pc_seconds_median=cells.seconds_pc.median(), setup_seconds_median=setup.median(),
                setup_seconds_max=setup.max(), peak_rss_worker_mb_max=worker_peak.max(),
                peak_rss_total_mb_est=workers * worker_peak.max() + peak_rss_mb(),
                edges_per_cell_median=cells.n_edges.median(), undirected_share_median=(
                    cells.n_undirected / cells.n_edges.where(cells.n_edges > 0)).median(),
                target_fraction_with_edge_mean=cells.target_fraction.mean(),
                ci_tests_per_cell_median=cells.n_ci_tests.median(), meek_valid_rate=cells.meek_valid.mean(),
                projected_wall_hours=(n_total * mean / workers + setup.max()) / 3600,
                config_wall_seconds=wall)


def target_edge_share(cells, targets):
    """Per configuration, target and condition: share of focal cells where the target has >= 1 edge."""
    rows = []
    for key, group in cells.groupby(KEYS + ["condition"], sort=True):
        present = [set(map(int, s.split(";"))) if s else set()
                   for s in group.targets_with_edge.fillna("").astype(str)]
        for i, t in enumerate(targets):
            rows.append(dict(zip(KEYS + ["condition"], key), target=t["target"], symbol=t["symbol"],
                             n_cells=len(group), share_cells_with_edge=np.mean([i in p for p in present])))
    return pd.DataFrame(rows)


def floats(text):
    return [float(v) for v in text.split(",")] if text else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "cscn_networks.json"))
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--sizes")
    parser.add_argument("--orders")
    parser.add_argument("--windows", help="Sensitivity: comma-separated CSCN windows (sigmoid_score)")
    parser.add_argument("--alphas", help="Sensitivity: comma-separated CI significance levels")
    parser.add_argument("--focal", type=int)
    parser.add_argument("--run-name")
    args = parser.parse_args(argv)
    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)
    with open(ROOT / cfg["paths"]["panel"], encoding="utf-8") as f:
        panel = json.load(f)
    t = cfg["timing"]
    sizes = [int(v) for v in args.sizes.split(",")] if args.sizes else t["panel_sizes"]
    orders = [int(v) for v in args.orders.split(",")] if args.orders else t["max_conditioning"]
    windows = floats(args.windows) or [cfg["window"]]
    alphas = floats(args.alphas) or [cfg["ci_alpha"]]
    sensitivity = args.windows is not None or args.alphas is not None
    if sensitivity and (len(sizes) != 1 or len(orders) != 1):
        raise ValueError("A window x alpha sensitivity grid needs one panel size and one order")
    n_focal = args.focal or t["focal_per_condition"]
    full = (not sensitivity and sizes == t["panel_sizes"] and orders == t["max_conditioning"]
            and n_focal == t["focal_per_condition"])
    default_run = "cscn_panel_sensitivity" if sensitivity else ("cscn_panel_timing" if full else "smoke_cscn_panel_timing")
    run = args.run_name or default_run
    out = ROOT / cfg["paths"]["outputs_dir"] / run
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(out / "log.txt", mode="w", encoding="utf-8")])
    with open(out / "config.json", "w", encoding="utf-8") as f:
        json.dump(dict(cfg, sizes=sizes, orders=orders, windows=windows, alphas=alphas, focal_per_condition=n_focal,
                       workers=args.workers, cpu_count=os.cpu_count(), panel_frozen=panel.get("frozen"),
                       argv=sys.argv[1:]), f, indent=2)

    gene_ids, nodes = panel_nodes(panel, max(sizes))
    targets = [g for g in panel["forced_genes"] if g["role"] == "target"]
    n_targets = len(targets)
    if any(g["role"] != "target" for g in panel["forced_genes"][:n_targets]):
        raise ValueError("Target nodes must come first in the panel")
    values, meta = load_ntc_expression(ROOT / panel["parameters"]["expression"], gene_ids, cfg["conditions"])
    n_total = len(meta)
    rows_by = {c: np.flatnonzero(meta.condition.eq(c).to_numpy()) for c in cfg["conditions"]}
    focal = {c: focal_subset(len(r), n_focal, cfg["seed"], k) for k, (c, r) in enumerate(rows_by.items())}
    LOG.info("NTC cells %s; focal %d per condition; workers %d of %d CPUs",
             {c: len(r) for c, r in rows_by.items()}, n_focal, args.workers, os.cpu_count())

    grid, exceeded, all_cells = [], [], []
    configurations = [dict(panel_size=P, max_conditioning=m, window=w, ci_alpha=a)
                      for P in sorted(sizes) for m in sorted(orders) for w in windows for a in alphas]
    for key in configurations:
        P, m = key["panel_size"], key["max_conditioning"]
        if not sensitivity and any(p <= P and o <= m for p, o in exceeded):
            grid.append(dict(key, status="skipped", n_cells_timed=0))
            LOG.info("%s skipped (a cheaper configuration exceeded %ds)", key, t["threshold_seconds"])
            continue
        references = {c: values[r][:, :P] for c, r in rows_by.items()}
        net_cfg = dict(window=key["window"], ci_alpha=key["ci_alpha"], max_conditioning=m)
        cells, status, wall = time_configuration(references, focal, nodes[:P], n_targets, net_cfg, args.workers,
                                                 t["threshold_seconds"], t["hard_timeout_seconds"])
        if status != "ok":
            exceeded.append((P, m))
        summary = summarize(cells, status, key, args.workers, n_total, wall)
        grid.append(summary)
        all_cells.append(cells.assign(**key))
        LOG.info("%s %s: n=%d mean %.1fs median %.1fs p90 %.1fs max %.1fs; edges %.0f; targets %.2f; "
                 "projected %.2f h", key, status, len(cells), summary.get("seconds_mean", np.nan),
                 summary.get("seconds_median", np.nan), summary.get("seconds_p90", np.nan),
                 summary.get("seconds_max", np.nan), summary.get("edges_per_cell_median", np.nan),
                 summary.get("target_fraction_with_edge_mean", np.nan), summary.get("projected_wall_hours", np.nan))
        frame = pd.DataFrame(grid).assign(cpu_count=os.cpu_count(), workers=args.workers, n_ntc_cells=n_total)
        frame.to_csv(out / "timing_grid.csv", index=False)
        pd.concat(all_cells).to_csv(out / "cell_timings.csv.gz", index=False)
    frame = pd.DataFrame(grid).assign(cpu_count=os.cpu_count(), workers=args.workers, n_ntc_cells=n_total)
    frame.to_csv(out / "timing_grid.csv", index=False)
    if all_cells:
        shares = target_edge_share(pd.concat(all_cells, ignore_index=True), targets)
        shares.to_csv(out / "target_edge_share.csv", index=False)
    if full:
        summary_path = ROOT / cfg["paths"]["summary_dir"] / "cscn_panel_timing.csv"
        frame.to_csv(summary_path, index=False, lineterminator="\n")
        LOG.info("wrote %s", summary_path)


if __name__ == "__main__":
    main()
