"""Cell-specific CSCN skeletons and collider-oriented PDAGs on the network panel (NTC cells only).

Uses the pinned official CSCN (icpdag/_vendor) as icpdag/gene_ckm.py does: bitmap queries, stable
PC, no NMF, then unshielded-collider orientation. Intervention seeding is fold-specific and happens
later, so the stored graph is the PDAG after collider orientation. Unseeded Meek propagation plus
the official DAG completion is computed as a diagnostic only.

The reference set is all NTC cells of one condition; every cell can be a focal cell. Inference of a
focal cell is deterministic and independent of the other focal cells, so results do not depend on
chunking, sharding or worker count.
"""

from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
import hashlib
from io import StringIO
import json
import multiprocessing as mp
import os
from pathlib import Path
import time

import networkx as nx
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from ._vendor.cscn.core import CSCN
from .cscn_panel import ntc_only
from .intervention_orientation import valid_extension

EDGE_KEYS = ("cell_index", "source", "destination", "directed")
CELL_COLUMNS = ["cell_index", "n_edges", "n_directed", "n_undirected", "n_ci_tests", "seconds",
                "seconds_pc", "error", "meek_valid", "meek_error"]


def peak_rss_mb():
    """Peak resident memory of this process in MB (peak working set on Windows)."""
    try:
        import psutil

        info = psutil.Process().memory_info()
        if hasattr(info, "peak_wset"):
            return info.peak_wset / 2 ** 20
    except ImportError:
        pass
    import resource

    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def reference_model(values, nodes, cfg):
    """Official CSCN on one condition's reference cells (cells x panel nodes)."""
    model = CSCN(variable_names=list(nodes), sigmoid_score=cfg["window"], significance_level=cfg["ci_alpha"],
                 max_cond_vars=cfg["max_conditioning"], query_engine="bitmap", pc_var="stable",
                 enable_query_stats=True)
    model.run_core(np.ascontiguousarray(values, dtype=np.float32), use_nmf=False)
    return model


def encode_edges(pdag, nodes):
    """(source, destination, directed) per adjacency; undirected edges are stored with source < destination."""
    index = {n: i for i, n in enumerate(nodes)}
    rows = []
    for a, b in pdag.edges:
        i, j = index[a], index[b]
        if pdag.has_edge(b, a):
            if i < j:
                rows.append((i, j, False))
        else:
            rows.append((i, j, True))
    rows.sort(key=lambda r: (min(r[0], r[1]), max(r[0], r[1])))
    src = np.array([r[0] for r in rows], dtype=np.int16)
    dst = np.array([r[1] for r in rows], dtype=np.int16)
    return src, dst, np.array([r[2] for r in rows], dtype=bool)


def decode_edges(source, destination, directed, nodes):
    """PDAG (networkx DiGraph; undirected edges in both directions) from encoded edges."""
    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    for i, j, d in zip(source, destination, directed):
        graph.add_edge(nodes[i], nodes[j])
        if not d:
            graph.add_edge(nodes[j], nodes[i])
    return graph


def infer_cell(model, local, nodes, strict=True):
    """Skeleton and collider PDAG of one focal cell, plus the unseeded Meek / DAG-completion diagnostic.

    An 'ICT ERROR' from the official CI test raises when strict; otherwise it is recorded and the cell
    gets no edges.
    """
    start = time.perf_counter()
    model.clear_task_cache()
    model.reset_query_stats()
    diagnostics = StringIO()
    with redirect_stdout(diagnostics):
        skeleton, separators = model._estimate_pc_skeleton(int(local))
    stats = model.get_query_stats()
    n_tests = int(stats["ci_queries"] - stats["ci_cache_hits"])
    empty = (np.zeros(0, np.int16), np.zeros(0, np.int16), np.zeros(0, bool))
    row = dict(n_ci_tests=n_tests, error="", meek_valid=False, meek_error="")
    if "ICT ERROR:" in diagnostics.getvalue():
        if strict:
            raise RuntimeError(diagnostics.getvalue())
        row.update(error=diagnostics.getvalue().strip()[:500], meek_error="skipped: CI test error",
                   n_edges=0, n_directed=0, n_undirected=0, seconds_pc=time.perf_counter() - start,
                   seconds=time.perf_counter() - start)
        return row, empty
    pdag = model._orient_unshielded_colliders(skeleton, separators)
    row["seconds_pc"] = time.perf_counter() - start
    edges = encode_edges(pdag, nodes)
    adjacency = {frozenset((nodes[i], nodes[j])) for i, j in zip(edges[0], edges[1])}
    if adjacency != {frozenset(e) for e in skeleton.edges}:
        raise RuntimeError("Collider orientation changed the skeleton adjacencies")
    try:
        meek = model._apply_meek_rules(pdag.copy(), apply_r4=False)
        dag = model._build_dag_representation(meek)
        row["meek_valid"] = bool(valid_extension(meek, dag))
        if not row["meek_valid"]:
            row["meek_error"] = "DAG completion inconsistent with PDAG"
    except (ValueError, RuntimeError) as exc:
        row["meek_error"] = str(exc)[:500]
    row.update(n_edges=len(edges[0]), n_directed=int(edges[2].sum()), n_undirected=int((~edges[2]).sum()),
               seconds=time.perf_counter() - start)
    return row, edges


# ---------------------------------------------------------------- worker processes
_WORKER = {}


def _init_worker(references, nodes, cfg):
    """references: condition -> reference matrix. One model per condition, built once per worker."""
    _WORKER["limits"] = threadpool_limits(limits=1)
    _WORKER.update(nodes=list(nodes), models={}, setup={})
    for condition, values in references.items():
        start = time.perf_counter()
        _WORKER["models"][condition] = reference_model(values, nodes, cfg)
        _WORKER["setup"][condition] = time.perf_counter() - start


def run_chunk(task):
    """task: dict(condition, locals, strict). Returns per-cell rows and encoded edges."""
    model, nodes = _WORKER["models"][task["condition"]], _WORKER["nodes"]
    rows, edges = [], []
    for local in task["locals"]:
        row, e = infer_cell(model, local, nodes, task.get("strict", False))
        rows.append(dict(cell_index=int(local), **row))
        edges.append(np.full(len(e[0]), int(local), dtype=np.int32))
        edges.extend(e)
    cat = [np.concatenate(edges[k::4]) if edges else np.zeros(0) for k in range(4)]
    return dict(rows=rows, edges=dict(zip(EDGE_KEYS, cat)), pid=os.getpid(),
                setup_seconds=_WORKER["setup"][task["condition"]], peak_rss_mb=peak_rss_mb())


def worker_pool(references, nodes, cfg, workers):
    """Spawn-context process pool whose workers hold the reference models (identical on any OS)."""
    return ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"),
                               initializer=_init_worker, initargs=(references, list(nodes), cfg))


def merge_edges(parts):
    if not parts:
        return {k: np.zeros(0, dtype=t) for k, t in zip(EDGE_KEYS, (np.int32, np.int16, np.int16, bool))}
    out = {k: np.concatenate([p[k] for p in parts]) for k in EDGE_KEYS}
    order = np.lexsort((np.maximum(out["source"], out["destination"]),
                        np.minimum(out["source"], out["destination"]), out["cell_index"]))
    return {k: v[order].astype(t) for (k, v), t in zip(out.items(), (np.int32, np.int16, np.int16, bool))}


def infer_networks(values, meta, nodes, cfg, focal=None, workers=1, chunk_size=32, strict=False):
    """Graphs for focal cells (positions in values; default all) with all cells of values as reference.

    values: reference cells x nodes (one condition, NTC only); meta: the same cells' metadata.
    Returns (cells, edges): one row per focal cell and the concatenated encoded edges.
    """
    ntc_only(meta)
    if meta.condition.nunique() != 1 or len(meta) != len(values):
        raise ValueError("One condition per reference set; metadata must match the expression rows")
    condition = str(meta.condition.iloc[0])
    focal = np.arange(len(values)) if focal is None else np.sort(np.asarray(focal, dtype=int))
    tasks = [dict(condition=condition, locals=focal[i:i + chunk_size].tolist(), strict=strict)
             for i in range(0, len(focal), chunk_size)]
    references = {condition: values}
    if workers == 1:
        _init_worker(references, nodes, cfg)
        results = [run_chunk(t) for t in tasks]
    else:
        with worker_pool(references, nodes, cfg, workers) as pool:
            results = list(pool.map(run_chunk, tasks))
    cells = pd.DataFrame([r for res in results for r in res["rows"]], columns=CELL_COLUMNS)
    return cells.sort_values("cell_index").reset_index(drop=True), merge_edges([r["edges"] for r in results])


# ---------------------------------------------------------------- shards
def focal_digest(focal):
    return hashlib.sha256(np.asarray(focal, dtype=np.int64).tobytes()).hexdigest()


def _atomic_write(path, writer):
    tmp = path.with_name(path.name + ".tmp")
    writer(tmp)
    os.replace(tmp, path)


def run_shard(values, meta, nodes, cfg, out_dir, shard, n_shards, *, focal=None, workers=1,
              chunk_size=32, provenance=None, log=print):
    """Infer one shard of the focal cells into out_dir/shard_<i>_of_<n>/; resumable by chunk.

    meta: reference cells' metadata (condition, calib_fold, lig1_well; index = cell id). The shard
    is block `shard` of np.array_split over the sorted focal positions. Finished chunks are skipped.
    """
    ntc_only(meta)
    condition = str(meta.condition.iloc[0])
    focal = np.arange(len(values)) if focal is None else np.sort(np.asarray(focal, dtype=int))
    mine = np.array_split(focal, n_shards)[shard]
    directory = Path(out_dir) / f"shard_{shard:03d}_of_{n_shards:03d}"
    directory.mkdir(parents=True, exist_ok=True)
    config = dict(cfg, condition=condition, nodes=list(nodes), shard=shard, n_shards=n_shards,
                  n_reference=len(values), n_focal=len(focal), focal_sha256=focal_digest(focal),
                  chunk_size=chunk_size, **(provenance or {}))
    path = directory / "config.json"
    if path.exists():
        old = json.loads(path.read_text(encoding="utf-8"))
        if old != json.loads(json.dumps(config)):
            raise ValueError(f"{path} exists with a different configuration")
    else:
        _atomic_write(path, lambda p: p.write_text(json.dumps(config, indent=2), encoding="utf-8"))
    chunks = [mine[i:i + chunk_size] for i in range(0, len(mine), chunk_size)]
    pending = [(k, c) for k, c in enumerate(chunks) if not (directory / f"chunk_{k:05d}.npz").exists()]
    log(f"{condition} shard {shard}/{n_shards}: {len(mine)} cells, {len(pending)}/{len(chunks)} chunks to run")
    tasks = [dict(condition=condition, locals=c.tolist(), strict=False) for _, c in pending]

    def save(k, result):
        cells = pd.DataFrame(result["rows"], columns=CELL_COLUMNS)

        def write_edges(p):
            with open(p, "wb") as f:
                np.savez_compressed(f, **result["edges"])

        _atomic_write(directory / f"chunk_{k:05d}.csv.gz",
                      lambda p: cells.to_csv(p, index=False, compression="gzip"))
        _atomic_write(directory / f"chunk_{k:05d}.npz", write_edges)

    if pending:
        references = {condition: values}
        if workers == 1:
            _init_worker(references, nodes, cfg)
            for (k, _), task in zip(pending, tasks):
                save(k, run_chunk(task))
        else:
            with worker_pool(references, nodes, cfg, workers) as pool:
                for (k, _), result in zip(pending, pool.map(run_chunk, tasks)):
                    save(k, result)
    (directory / "done").write_text(str(len(chunks)), encoding="utf-8")
    return directory


def read_shard(directory):
    directory = Path(directory)
    n = int((directory / "done").read_text(encoding="utf-8"))
    cells = [pd.read_csv(directory / f"chunk_{k:05d}.csv.gz", keep_default_na=False,
                         dtype={"error": str, "meek_error": str}) for k in range(n)]
    edges = []
    for k in range(n):
        with np.load(directory / f"chunk_{k:05d}.npz") as z:
            edges.append({key: z[key] for key in EDGE_KEYS})
    frame = pd.concat(cells, ignore_index=True) if cells else pd.DataFrame(columns=CELL_COLUMNS)
    return frame, merge_edges(edges)


def merge_shards(out_dir, meta):
    """Merge all shards of one condition; checks configs, completeness and focal coverage.

    meta: the reference cells' metadata in reference order (index = cell id). Writes cells.csv.gz
    (with cell_id, condition, calib_fold, lig1_well) and edges.npz plus nodes to out_dir.
    """
    out_dir = Path(out_dir)
    shards = sorted(out_dir.glob("shard_*_of_*"))
    if not shards:
        raise ValueError(f"No shards in {out_dir}")
    configs = [json.loads((s / "config.json").read_text(encoding="utf-8")) for s in shards]
    common = [{k: v for k, v in c.items() if k != "shard"} for c in configs]
    if any(c != common[0] for c in common):
        raise ValueError("Shards were run with different configurations")
    n_shards = common[0]["n_shards"]
    if sorted(c["shard"] for c in configs) != list(range(n_shards)):
        raise ValueError(f"Expected shards 0..{n_shards - 1}, found {[c['shard'] for c in configs]}")
    missing = [s.name for s in shards if not (s / "done").exists()]
    if missing:
        raise ValueError(f"Unfinished shards: {missing}")
    parts = [read_shard(s) for s in shards]
    cells = pd.concat([p[0] for p in parts], ignore_index=True).sort_values("cell_index").reset_index(drop=True)
    edges = merge_edges([p[1] for p in parts])
    if cells.cell_index.duplicated().any() or focal_digest(cells.cell_index) != common[0]["focal_sha256"]:
        raise ValueError("Merged cells do not cover the focal set exactly once")
    if len(meta) != common[0]["n_reference"]:
        raise ValueError("Metadata does not match the reference set")
    info = meta.iloc[cells.cell_index.to_numpy()]
    cells.insert(1, "cell_id", info.index.to_numpy())
    for col in ("condition", "calib_fold", "lig1_well"):
        cells.insert(cells.columns.get_loc("n_edges"), col, info[col].to_numpy())
    counts = np.bincount(np.searchsorted(cells.cell_index.to_numpy(), edges["cell_index"]), minlength=len(cells))
    if not np.array_equal(counts, cells.n_edges.to_numpy()):
        raise RuntimeError("Edge lists do not match the per-cell edge counts")
    cells.to_csv(out_dir / "cells.csv.gz", index=False, compression="gzip")
    np.savez_compressed(out_dir / "edges.npz", nodes=np.asarray(common[0]["nodes"], dtype=str), **edges)
    (out_dir / "config.json").write_text(json.dumps(common[0], indent=2), encoding="utf-8")
    return cells, edges
