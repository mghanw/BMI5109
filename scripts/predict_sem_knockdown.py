"""Knockdown predictions of every SEM model for the held-out targets of one (condition, outer fold).

Usage:
    python scripts/predict_sem_knockdown.py --condition DMSO --fold 0 [--shard 0/1] [--workers W]
    python scripts/predict_sem_knockdown.py --condition DMSO --fold 0 --merge
    Smoke test (subset of held-out targets and reference cells; outputs/<run>/ only):
    python scripts/predict_sem_knockdown.py --condition DMSO --fold 0 --targets A,B,C --cells 200 \
        [--run-name sem_prediction_smoke]

Definitions: config/sem_prediction.json (must be frozen). Needs scripts/build_sem_graphs.py and
scripts/sem_lambda_cv.py first. Response data read: training-fold targets (orientation seeds) and
the held-out target's own perturbed gene per calibration fold (delta); nothing else.

Run phase (resumable by chunk, shardable over reference cells): per reference cell, the cell-level
models (cell_specific_cscn_sem, unweighted_propagation_ablation, local_association_network_sem).
Merge phase: checks the shards, adds the pooled-coefficient models (shared_cscn_sem,
cluster_level_sem) and zero_downstream_baseline, and writes per-unit and delta-scaled pseudobulk
predictions, per-cell Fufa score shifts and statistics to <out>/<condition>/fold<f>/predictions.npz.
"""
import argparse
from collections import Counter
import hashlib
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
from icpdag.cscn_networks import reference_model  # noqa: E402
from icpdag.sem_graphs import association_pairs  # noqa: E402
from icpdag.sem_inputs import condition_inputs, fold_targets, load_panel_cells, target_nodes  # noqa: E402
from icpdag.sem_orientation import fold_seeds, undirected_pdag  # noqa: E402
from icpdag.sem_prediction import (MODELS, LocalCoefficients, PooledCoefficients, cell_response,  # noqa: E402
                                   fufa_weights, graph_responses, group_responses, target_delta,
                                   unweighted_coefficients, zero_response)
from icpdag.target_split import load_split  # noqa: E402

CELL_MODELS = ["cell_specific_cscn_sem", "unweighted_propagation_ablation", "local_association_network_sem"]
STAT_KEYS = ("n_undirected", "n_candidates", "n_failed", "n_fallback")
_W = {}


# ---------------------------------------------------------------- frozen Fufa weights
def frozen_fufa_weights(cfg, nodes, cache):
    """Panel weights of the frozen Aim 1 Fufa score (same genes, controls, centre and scale)."""
    if cache.exists():
        spec = json.loads(cache.read_text(encoding="utf-8"))
    else:
        import aim1
        aim1_cfg = json.loads((ROOT / cfg["inputs"]["program_scoring"]).read_text(encoding="utf-8"))
        counts, var_all, obs = aim1.load_data(aim1_cfg)
        keep = aim1.gs.scoring_genes(var_all)
        var = var_all[keep]
        X = aim1.log_cp10k(counts[:, keep], obs.total_counts.to_numpy(float))
        meta = obs[["condition", "target"]].copy()
        library_symbols, _, _ = aim1.library_genes(aim1_cfg, var, sorted(load_split(ROOT / cfg["inputs"]["target_split"]).target))
        _, mapping = aim1.load_programs(ROOT / aim1_cfg["paths"]["programs"])
        programs = aim1.gs.programs_for_scoring(mapping)
        _, _, audit = aim1.score_cells(X, var, meta, programs, library_symbols, aim1_cfg["seed"])
        fufa = next(a for a in audit if a["program"] == aim1_cfg["primary_readout"])
        frozen = pd.read_csv(ROOT / "results" / "summary" / "aim1_program_audit.csv").set_index("program")
        row = frozen.loc[aim1_cfg["primary_readout"]]
        if not (np.isclose(row.center, fufa["center"]) and np.isclose(row.scale, fufa["scale"])
                and int(row.n_background_controls) == len(fufa["controls"])):
            raise RuntimeError("Recomputed Fufa scoring does not match the frozen Aim 1 audit")
        spec = dict(program=aim1_cfg["primary_readout"], program_genes=fufa["mapped_genes"],
                    controls=fufa["controls"], center=fufa["center"], scale=fufa["scale"])
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(spec, indent=2), encoding="utf-8")
    w, counts = fufa_weights(spec["program_genes"], spec["controls"], spec["scale"], nodes)
    return w, counts


# ---------------------------------------------------------------- worker
def _init(state):
    _W.update(state)
    _W["model"] = reference_model(state["values"], state["nodes"], state["cscn"])
    _W["pooled"] = PooledCoefficients(state["values"], state["index"])


def _run(cells):
    values, index, held = _W["values"], _W["index"], _W["held"]
    T, P = len(held), len(index)
    sums = {m: np.zeros((T, P)) for m in CELL_MODELS}
    counts = {m: np.zeros(T, dtype=np.int64) for m in CELL_MODELS}
    fufa = {m: np.full((len(cells), T), np.nan) for m in CELL_MODELS}
    stats, seconds = Counter(), Counter()
    for row, k in enumerate(cells):
        local = LocalCoefficients(values, index, _W["neighbours"][k], _W["pooled"], _W["lam"])
        t0 = time.perf_counter()
        cell = graph_responses(_W["graphs"].pdag(k), _W["seeds"], held, _W["model"], local, index)
        t1 = time.perf_counter()
        pairs = association_pairs(values[_W["neighbours"][k]], _W["graphs"].n_edges(k))
        start = undirected_pdag(_W["nodes"], [(_W["nodes"][i], _W["nodes"][j]) for i, j in pairs])
        assoc = graph_responses(start, _W["seeds"], held, _W["model"], local, index)  # same rows, same cache
        t2 = time.perf_counter()
        unweighted = {g: cell_response(cell[g][2], g, unweighted_coefficients, index) for g in held}
        t3 = time.perf_counter()
        for t, g in enumerate(held):
            for m, (response, s) in {"cell_specific_cscn_sem": cell[g][:2],
                                     "unweighted_propagation_ablation": (unweighted[g], cell[g][1]),
                                     "local_association_network_sem": assoc[g][:2]}.items():
                stats[(m, g, "units")] += 1
                for key in STAT_KEYS:
                    stats[(m, g, key)] += s[key]
                if response is None:
                    stats[(m, g, "excluded")] += 1
                    continue
                sums[m][t] += response
                counts[m][t] += 1
                fufa[m][row, t] = float(_W["w"] @ response)
        seconds["cell_specific_cscn_sem"] += t1 - t0
        seconds["local_association_network_sem"] += t2 - t1
        seconds["unweighted_propagation_ablation"] += t3 - t2
    return dict(cells=list(cells), sums=sums, counts=counts, fufa=fufa, stats=dict(stats), seconds=dict(seconds))


# ---------------------------------------------------------------- helpers
def digest(array):
    return hashlib.sha256(np.asarray(array, dtype=np.int64).tobytes()).hexdigest()


def save_chunk(directory, k, result):
    tmp = directory / f"chunk_{k:05d}.tmp.npz"
    arrays = {f"{kind}_{m}": result[kind][m] for kind in ("sums", "counts", "fufa") for m in CELL_MODELS}
    with open(tmp, "wb") as f:
        np.savez_compressed(f, cells=np.asarray(result["cells"]), **arrays)
    meta = dict(stats=[[list(map(str, key)), int(v)] for key, v in result["stats"].items()], seconds=result["seconds"])
    (directory / f"chunk_{k:05d}.json").write_text(json.dumps(meta), encoding="utf-8")
    os.replace(tmp, directory / f"chunk_{k:05d}.npz")


def load_inputs(cfg, args):
    inp = condition_inputs(ROOT, cfg, args.condition)
    split = load_split(ROOT / cfg["inputs"]["target_split"])
    held_targets, train = fold_targets(split, args.fold)
    if args.targets:
        chosen = args.targets.split(",")
        if set(chosen) - set(held_targets):
            raise ValueError(f"Not held out in fold {args.fold}: {sorted(set(chosen) - set(held_targets))}")
        held_targets = sorted(chosen)
    tnodes = target_nodes(inp["panel"])
    panel_values, panel_meta = load_panel_cells(ROOT / cfg["inputs"]["expression"], inp["ids"], [args.condition],
                                                split.target)
    panel_values = pd.DataFrame(panel_values, index=panel_meta.index, columns=inp["nodes"])
    seeds, _ = fold_seeds(panel_values, panel_meta, train, tnodes, args.condition)
    deltas = np.array([[target_delta(panel_values, panel_meta, t, tnodes[t], args.condition)[c] for c in (0, 1)]
                       for t in held_targets])
    del panel_values, panel_meta
    lam_spec = json.loads((ROOT / cfg["outputs_dir"] / "lambda_cv" / f"{args.condition}.json").read_text(encoding="utf-8"))
    if lam_spec.get("smoke_test"):
        raise ValueError("lambda comes from a smoke cross-validation")
    w, w_counts = frozen_fufa_weights(cfg, inp["nodes"], ROOT / cfg["outputs_dir"] / "fufa_weights.json")
    refs = np.arange(len(inp["values"]))
    if args.cells:
        refs = np.sort(np.random.default_rng(cfg["seed"]).choice(refs, args.cells, replace=False))
    return dict(inp=inp, held_targets=held_targets, held=[tnodes[t] for t in held_targets], seeds=seeds,
                deltas=deltas, lam=float(lam_spec["lambda_selected"]), w=w, w_counts=w_counts, refs=refs)


# ---------------------------------------------------------------- merge
def merge(cfg, args, data, directory, networks):
    inp, held, refs = data["inp"], data["held"], data["refs"]
    shards = sorted(directory.glob("shard_*_of_*"))
    configs = [json.loads((s / "config.json").read_text(encoding="utf-8")) for s in shards]
    common = [{k: v for k, v in c.items() if k != "shard"} for c in configs]
    if not shards or any(c != common[0] for c in common) or \
            sorted(c["shard"] for c in configs) != list(range(common[0]["n_shards"])):
        raise ValueError("Missing shards or shards with different configurations")
    if any(not (s / "done").exists() for s in shards):
        raise ValueError("Unfinished shards")
    T, P, n_ref = len(held), len(inp["nodes"]), len(refs)
    position = {int(k): i for i, k in enumerate(refs)}
    per_unit = np.zeros((len(MODELS), T, P))
    included = np.zeros((len(MODELS), T), dtype=np.int64)
    fufa = np.full((len(MODELS), n_ref, T), np.nan)
    stats, seconds, seen = Counter(), Counter(), []
    for s in shards:
        for k in range(int((s / "done").read_text())):
            with np.load(s / f"chunk_{k:05d}.npz") as z:
                cells = z["cells"]
                seen.extend(cells.tolist())
                rows = [position[int(c)] for c in cells]
                for m in CELL_MODELS:
                    i = MODELS.index(m)
                    per_unit[i] += z[f"sums_{m}"]
                    included[i] += z[f"counts_{m}"]
                    fufa[i, rows] = z[f"fufa_{m}"]
            meta = json.loads((s / f"chunk_{k:05d}.json").read_text(encoding="utf-8"))
            stats.update({tuple(key): v for key, v in meta["stats"]})
            seconds.update(meta["seconds"])
    if sorted(seen) != sorted(position) or len(seen) != len(set(seen)):
        raise ValueError("Chunks do not cover the reference cells exactly once")
    model = reference_model(inp["values"], inp["nodes"], dict(window=networks["window"], ci_alpha=networks["ci_alpha"],
                                                              max_conditioning=networks["max_conditioning"]))
    t0 = time.perf_counter()
    groups = {"shared_cscn_sem": group_responses({0: inp["structure"]["shared"]}, np.zeros(len(inp["values"])),
                                                 inp["values"], data["seeds"], held, model, inp["index"])}
    seconds["shared_cscn_sem"] += time.perf_counter() - t0
    t0 = time.perf_counter()
    labels = inp["labels"]
    groups["cluster_level_sem"] = group_responses(
        {c: inp["structure"][f"cluster_{c + 1}"] for c in np.unique(labels)}, labels, inp["values"], data["seeds"],
        held, model, inp["index"])
    seconds["cluster_level_sem"] += time.perf_counter() - t0
    ref_labels = {"shared_cscn_sem": np.zeros(n_ref, dtype=int), "cluster_level_sem": labels[refs]}
    for m, by_group in groups.items():
        i = MODELS.index(m)
        for t, g in enumerate(held):
            for group, responses in by_group.items():
                response, s, _ = responses[g]
                members = ref_labels[m] == group
                stats[(m, g, "units")] += int(members.sum())
                for key in STAT_KEYS:
                    stats[(m, g, key)] += int(members.sum()) * s[key]
                if response is None:
                    stats[(m, g, "excluded")] += int(members.sum())
                    continue
                per_unit[i, t] += members.sum() * response
                included[i, t] += members.sum()
                fufa[i, members, t] = float(data["w"] @ response)
    i = MODELS.index("zero_downstream_baseline")
    for t, g in enumerate(held):
        response = zero_response(g, inp["index"])
        per_unit[i, t] = n_ref * response
        included[i, t] = n_ref
        fufa[i, :, t] = float(data["w"] @ response)
        stats[("zero_downstream_baseline", g, "units")] += n_ref
    with np.errstate(invalid="ignore", divide="ignore"):
        per_unit /= included[:, :, None]
    scaled = per_unit[:, :, None, :] * data["deltas"][None, :, :, None]
    out = directory / "predictions.npz"
    np.savez_compressed(out, models=np.asarray(MODELS), targets=np.asarray(data["held_targets"]),
                        nodes=np.asarray(inp["nodes"]), per_unit=per_unit, scaled=scaled, delta=data["deltas"],
                        included=included, fufa_per_unit=fufa,
                        fufa_scaled=fufa[:, :, :, None] * data["deltas"][None, None, :, :],
                        reference_cells=refs, lam=data["lam"])
    rows = []
    for m in MODELS:
        for t, (g, target) in enumerate(zip(held, data["held_targets"])):
            row = dict(model=m, condition=args.condition, fold=args.fold, target=target, node=g,
                       units=stats[(m, g, "units")], excluded=stats[(m, g, "excluded")])
            row.update({key: stats[(m, g, key)] for key in STAT_KEYS})
            rows.append(row)
    pd.DataFrame(rows).to_csv(directory / "prediction_stats.csv", index=False)
    (directory / "merge.json").write_text(json.dumps(dict(seconds=dict(seconds), n_reference=n_ref,
                                                          w_counts=data["w_counts"]), indent=2), encoding="utf-8")
    return out, seconds


# ---------------------------------------------------------------- main
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(ROOT / "config" / "sem_prediction.json"))
    parser.add_argument("--condition", required=True)
    parser.add_argument("--fold", type=int, required=True)
    parser.add_argument("--shard", default="0/1")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--targets", help="Smoke test: comma-separated held-out targets")
    parser.add_argument("--cells", type=int, help="Smoke test: seeded subset of reference cells")
    parser.add_argument("--run-name")
    parser.add_argument("--chunk-size", type=int, default=8)
    args = parser.parse_args(argv)
    cfg = json.loads(Path(args.config).read_text(encoding="utf-8"))
    if not cfg.get("frozen"):
        raise ValueError("config/sem_prediction.json must be frozen (committed) before any prediction run")
    networks = json.loads((ROOT / cfg["inputs"]["networks"]).read_text(encoding="utf-8"))
    smoke = bool(args.targets or args.cells)
    run = args.run_name or ("sem_prediction_smoke" if smoke else "sem_prediction/predictions")
    directory = ROOT / "outputs" / run / args.condition / f"fold{args.fold}"
    directory.mkdir(parents=True, exist_ok=True)
    t_load = time.time()
    data = load_inputs(cfg, args)
    print(f"{args.condition} fold {args.fold}: {len(data['held'])} held-out targets, {len(data['refs'])} reference "
          f"cells, {len(data['seeds'])} seed proposals, lambda {data['lam']} ({time.time() - t_load:.0f}s)", flush=True)
    if args.merge:
        out, seconds = merge(cfg, args, data, directory, networks)
        print(f"wrote {out}; model seconds {json.dumps({k: round(v, 1) for k, v in seconds.items()})}")
        return
    shard, n_shards = (int(v) for v in args.shard.split("/"))
    mine = np.array_split(data["refs"], n_shards)[shard]
    sdir = directory / f"shard_{shard:03d}_of_{n_shards:03d}"
    sdir.mkdir(parents=True, exist_ok=True)
    config = dict(condition=args.condition, fold=args.fold, shard=shard, n_shards=n_shards,
                  targets=data["held_targets"], lam=data["lam"], reference_sha256=digest(data["refs"]),
                  chunk_size=args.chunk_size, config_sha256=hashlib.sha256(Path(args.config).read_bytes()
                                                                           .replace(b"\r\n", b"\n")).hexdigest())
    path = sdir / "config.json"
    if path.exists() and json.loads(path.read_text(encoding="utf-8")) != config:
        raise ValueError(f"{path} exists with a different configuration")
    path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    chunks = [mine[i:i + args.chunk_size].tolist() for i in range(0, len(mine), args.chunk_size)]
    pending = [k for k in range(len(chunks)) if not (sdir / f"chunk_{k:05d}.npz").exists()]
    inp = data["inp"]
    state = dict(values=inp["values"], nodes=inp["nodes"], index=inp["index"], graphs=inp["graphs"],
                 neighbours=inp["neighbours"], seeds=data["seeds"], held=data["held"], lam=data["lam"], w=data["w"],
                 cscn=dict(window=networks["window"], ci_alpha=networks["ci_alpha"],
                           max_conditioning=networks["max_conditioning"]))
    t0 = time.time()
    print(f"shard {shard}/{n_shards}: {len(mine)} cells, {len(pending)}/{len(chunks)} chunks to run", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=mp.get_context("spawn"),
                             initializer=_init, initargs=(state,)) as pool:
        for k, result in zip(pending, pool.map(_run, [chunks[k] for k in pending])):
            save_chunk(sdir, k, result)
    (sdir / "done").write_text(str(len(chunks)), encoding="utf-8")
    print(f"shard done ({time.time() - t0:.0f}s wall)")


if __name__ == "__main__":
    main()
