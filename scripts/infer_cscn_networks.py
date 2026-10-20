"""Cell-specific CSCN graphs for all NTC cells of one condition (resumable, shardable), and the merge.

Usage:
    python scripts/infer_cscn_networks.py --condition DMSO [--shard 0/4] [--workers W]
    python scripts/infer_cscn_networks.py --condition DMSO --merge
    Smoke test (seeded subset; writes only to outputs/<run>/, panel need not be frozen):
    python scripts/infer_cscn_networks.py --condition DMSO --cells 200 --panel-size 250 \
        --max-conditioning 2 [--shard 0/2] [--run-name cscn_panel_networks_smoke]

The reference set is all NTC cells of the condition. A full run needs a frozen panel and
max_conditioning set in config/cscn_networks.json. Each shard writes
outputs/<run>/<condition>/shard_<i>_of_<n>/ (config.json, chunk files, done); rerunning skips
finished chunks. --merge checks the shards and writes <condition>/cells.csv.gz and edges.npz.
Paths are pathlib-only and workers use the spawn start method, so the same command runs as a PBS
array job on Linux.
"""
import argparse
import hashlib
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from icpdag.cscn_networks import merge_shards, run_shard  # noqa: E402
from icpdag.cscn_panel import load_ntc_expression, panel_of_size  # noqa: E402

LOG = logging.getLogger("cscn_networks")


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def panel_nodes(panel, size):
    """Panel gene IDs and node names (symbols) of the size-P panel."""
    forced = [g["gene_id"] for g in panel["forced_genes"]]
    ranked = [g for g, _ in panel["fill_genes_ranked"]]
    symbol = {g["gene_id"]: g["symbol"] for g in panel["forced_genes"]}
    symbol.update({g: s for g, s in panel["fill_genes_ranked"]})
    ids = panel_of_size(forced, ranked, size)
    return ids, [symbol[g] for g in ids]


def focal_subset(n, size, seed, condition_index):
    rng = np.random.default_rng([seed, condition_index])
    return np.sort(rng.choice(n, size=size, replace=False))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "cscn_networks.json"))
    parser.add_argument("--condition", required=True)
    parser.add_argument("--shard", default="0/1", help="i/n: run block i of n")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--cells", type=int, help="Smoke test: seeded subset of focal cells")
    parser.add_argument("--panel-size", type=int, help="Smoke test only (full runs use the frozen size)")
    parser.add_argument("--max-conditioning", type=int, help="Smoke test only (full runs use the config)")
    parser.add_argument("--run-name")
    args = parser.parse_args(argv)
    with open(args.config, encoding="utf-8") as f:
        cfg = json.load(f)
    with open(ROOT / cfg["paths"]["panel"], encoding="utf-8") as f:
        panel = json.load(f)
    if args.condition not in cfg["conditions"]:
        raise ValueError(f"Unknown condition {args.condition}")
    smoke = args.cells is not None
    if smoke:
        size, order = args.panel_size, args.max_conditioning
        if size is None or order is None:
            raise ValueError("A smoke test needs --panel-size and --max-conditioning")
    else:
        if args.panel_size is not None or args.max_conditioning is not None:
            raise ValueError("--panel-size / --max-conditioning are for smoke tests; full runs use the configs")
        if not panel.get("frozen") or cfg["max_conditioning"] is None:
            raise ValueError("A full run needs a frozen panel and max_conditioning in the network config")
        size, order = panel["panel_size"], cfg["max_conditioning"]
    run = args.run_name or ("cscn_panel_networks_smoke" if smoke else "cscn_panel_networks")
    out = ROOT / cfg["paths"]["outputs_dir"] / run / args.condition
    out.mkdir(parents=True, exist_ok=True)
    shard, n_shards = (int(v) for v in args.shard.split("/"))
    if not 0 <= shard < n_shards:
        raise ValueError("--shard must be i/n with 0 <= i < n")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(out / ("merge_log.txt" if args.merge else f"log_shard_{shard:03d}_of_{n_shards:03d}.txt"),
                            mode="a", encoding="utf-8")])

    gene_ids, nodes = panel_nodes(panel, size)
    values, meta = load_ntc_expression(ROOT / panel["parameters"]["expression"], gene_ids, [args.condition])
    LOG.info("%s: %d NTC reference cells x %d panel genes", args.condition, *values.shape)
    if args.merge:
        cells, edges = merge_shards(out, meta)
        LOG.info("merged %d cells, %d edges, %d CI errors", len(cells), len(edges["source"]),
                 int(cells.error.astype(str).str.len().gt(0).sum()))
        return
    network_cfg = dict(window=cfg["window"], ci_alpha=cfg["ci_alpha"], max_conditioning=order,
                       query_engine=cfg["query_engine"], pc_variant=cfg["pc_variant"])
    focal = None
    if smoke:
        focal = focal_subset(len(values), args.cells, cfg["seed"], cfg["conditions"].index(args.condition))
    provenance = dict(panel_size=size, panel_sha256=file_sha256(ROOT / cfg["paths"]["panel"]),
                      expression_sha256=panel["inputs"]["expression_sha256"], smoke_test=smoke)
    run_config = dict(network_cfg, **provenance, conditions=cfg["conditions"], run_name=run)
    run_path = out.parent / "config.json"
    if run_path.exists() and json.loads(run_path.read_text(encoding="utf-8")) != run_config:
        raise ValueError(f"{run_path} exists with a different configuration")
    run_path.write_text(json.dumps(run_config, indent=2), encoding="utf-8")
    directory = run_shard(values, meta, nodes, network_cfg, out, shard, n_shards, focal=focal,
                          workers=args.workers, chunk_size=cfg["chunk_size"], provenance=provenance, log=LOG.info)
    with open(directory / "run_info.jsonl", "a", encoding="utf-8") as f:
        f.write(json.dumps(dict(git_commit=git_commit(), argv=sys.argv[1:], workers=args.workers,
                                cpu_count=os.cpu_count())) + "\n")
    LOG.info("finished %s", directory)


if __name__ == "__main__":
    main()
