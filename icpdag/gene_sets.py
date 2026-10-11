"""External program gene sets: parsing, GENCODE v38 symbol mapping and target overlap.

Source files and their provenance are listed in data/ref/SOURCES.md (scripts/fetch_refs.py).
"""

import re

import numpy as np
import pandas as pd

FUFA = "Fufa2019_AZD_up_SOX10_bound_pigmentation"
KANAI_UP = "Kanai2024_YAP_TAZ_up"
KANAI_DOWN = "Kanai2024_YAP_TAZ_down"
TSOI = {"Tsoi_undifferentiated": "Undifferentiated", "Tsoi_neural_crest_like": "Neural crest-like",
        "Tsoi_transitory": "Transitory", "Tsoi_melanocytic": "Melanocytic"}
# Table S3 labels genes up in one subtype or shared by two adjacent subtypes.
TSOI_LABELS = {
    "Undifferentiated": ["Undifferentiated"],
    "Undifferentiated-Neural crest-like": ["Undifferentiated", "Neural crest-like"],
    "Neural crest-like": ["Neural crest-like"],
    "Neural crest-like-Transitory": ["Neural crest-like", "Transitory"],
    "Transitory": ["Transitory"],
    "Transitory-Melanocytic": ["Transitory", "Melanocytic"],
    "Melanocytic": ["Melanocytic"],
}
# Order of the earlier analysis (fixes the background-control random draws in program_scores).
ORDER = [*TSOI, KANAI_UP, KANAI_DOWN, FUFA]
ROLES = {FUFA: "primary", KANAI_UP: "secondary", KANAI_DOWN: "exploratory",
         **{name: "secondary" for name in TSOI}}


def _clean(values):
    return [str(v).strip() for v in values if isinstance(v, str) and str(v).strip()]


def parse_fufa(table):
    """Fufa 2019 Table S2: pigmentation genes up under AZD6244 with SOX10 binding under AZD6244.

    table: the sheet read with header row 'Gene Symbol', 'Log2 Fold Change', ..., 'SOX10 binding ...'.
    """
    frame = table.iloc[:, :5].copy()
    frame.columns = ["symbol", "log2fc", "h3k27ac_loss", "h3k27ac_gain", "sox10"]
    frame = frame.dropna(subset=["symbol"])
    bound = frame.sox10.astype(str).str.strip().str.upper().str.startswith("X")
    keep = frame.log2fc.astype(float).gt(0) & bound
    return _clean(frame.loc[keep, "symbol"])


def parse_kanai(table):
    """Kanai 2024 Table S3 tab 6: 132-gene signature; columns 13/14 hold the up/down totals.

    table: the sheet read without a header (two header rows, then genes).
    """
    body = table.iloc[2:]
    up_total = pd.to_numeric(body[13], errors="coerce").fillna(0)
    down_total = pd.to_numeric(body[14], errors="coerce").fillna(0)
    if (up_total.gt(0) & down_total.gt(0)).any():
        raise ValueError("A Kanai signature gene is counted as both up and down")
    return {KANAI_UP: _clean(body.loc[up_total.gt(0), 0]), KANAI_DOWN: _clean(body.loc[down_total.gt(0), 0])}


def parse_tsoi(table):
    """Tsoi 2018 Table S3: a subtype program is the genes labelled with that subtype only.

    Genes shared by two adjacent subtypes (e.g. 'Undifferentiated-Neural crest-like') are not used.
    This own-label definition reproduces the earlier report's Tsoi contrasts.
    """
    unknown = set(table.Signature.dropna()) - set(TSOI_LABELS)
    if unknown:
        raise ValueError(f"Unknown Tsoi signature labels: {sorted(unknown)}")
    return {name: _clean(table.loc[table.Signature.eq(subtype), "Gene"]) for name, subtype in TSOI.items()}


def read_sources(paths):
    """Read the three supplementary tables; paths: dict with keys fufa, kanai, tsoi."""
    fufa = parse_fufa(pd.read_excel(paths["fufa"], sheet_name="Table S2", header=1))
    kanai = parse_kanai(pd.read_excel(paths["kanai"], sheet_name="Tab 6-YAPTAZ Signature", header=None))
    tsoi = parse_tsoi(pd.read_excel(paths["tsoi"], sheet_name="Table S3", header=1))
    return {name: {**tsoi, **kanai, FUFA: fufa}[name] for name in ORDER}


def scoring_genes(var):
    """Mask of matrix genes used for scoring: no PAR_Y rows, one gene per symbol.

    A duplicated symbol keeps its protein-coding gene (then the lowest gene_id).
    """
    frame = var.reset_index().rename(columns={var.index.name or "index": "gene_id"})
    frame["row"] = np.arange(len(frame))
    frame = frame[~frame.is_par_y.astype(bool)]
    frame = frame.assign(rank=frame.gene_type.ne("protein_coding").astype(int))
    keep = frame.sort_values(["symbol", "rank", "gene_id"]).drop_duplicates("symbol").row
    mask = np.zeros(len(var), bool)
    mask[keep.to_numpy()] = True
    return mask


def _split(value):
    return [] if not isinstance(value, str) else [s.strip() for s in value.split("|") if s.strip()]


def hgnc_index(hgnc):
    """Lookups from a symbol to HGNC rows, by approved/previous symbol and by alias."""
    by_name, by_alias = {}, {}
    for i, row in enumerate(hgnc.itertuples()):
        for name in [row.symbol, *_split(row.prev_symbol)]:
            by_name.setdefault(name, set()).add(i)
        for name in _split(row.alias_symbol):
            by_alias.setdefault(name, set()).add(i)
    return by_name, by_alias


def map_symbols(source, genes, hgnc):
    """Map source symbols to scoring genes (GENCODE v38).

    genes: DataFrame with gene_id (unversioned or versioned), symbol, for scoring genes only.
    hgnc: HGNC complete set with symbol, prev_symbol, alias_symbol, ensembl_gene_id.
    Exact symbol match first; otherwise an HGNC approved/previous-symbol match, then an alias
    match, accepted only when it identifies exactly one scoring gene (by Ensembl ID or symbol).
    """
    genes = genes.assign(base=genes.gene_id.str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True))
    by_symbol = genes.set_index("symbol")
    by_base = genes.set_index("base")
    by_name, by_alias = hgnc_index(hgnc)
    rows = []
    for symbol in source:
        if symbol in by_symbol.index:
            rows.append(dict(source_symbol=symbol, symbol=symbol, gene_id=by_symbol.loc[symbol, "gene_id"],
                             match_type="exact"))
            continue
        hit = None
        for kind, index in (("hgnc_previous", by_name), ("hgnc_alias", by_alias)):
            candidates = set()
            for i in index.get(symbol, ()):
                record = hgnc.iloc[i]
                if isinstance(record.ensembl_gene_id, str) and record.ensembl_gene_id in by_base.index:
                    candidates.add(by_base.loc[record.ensembl_gene_id, "symbol"])
                elif record.symbol in by_symbol.index:
                    candidates.add(record.symbol)
            if len(candidates) == 1:
                hit = (kind, candidates.pop())
                break
            if len(candidates) > 1:
                hit = ("ambiguous", None)
                break
        if hit and hit[1] is not None:
            rows.append(dict(source_symbol=symbol, symbol=hit[1], gene_id=by_symbol.loc[hit[1], "gene_id"],
                             match_type=hit[0]))
        else:
            rows.append(dict(source_symbol=symbol, symbol=None, gene_id=None,
                             match_type="ambiguous" if hit else "missing"))
    return pd.DataFrame(rows, columns=["source_symbol", "symbol", "gene_id", "match_type"])


def mapping_table(sources, genes, hgnc):
    """One row per (program, source gene); later duplicates of a mapped symbol are flagged."""
    frames = []
    for name, source in sources.items():
        mapped = map_symbols(source, genes, hgnc).assign(program=name, readout_role=ROLES[name])
        mapped["duplicate_in_program"] = mapped.symbol.notna() & mapped.symbol.duplicated()
        frames.append(mapped)
    return pd.concat(frames, ignore_index=True)[["program", "readout_role", "source_symbol", "symbol",
                                                 "gene_id", "match_type", "duplicate_in_program"]]


def programs_for_scoring(mapping):
    """program_scores input: mapped v38 symbol, or the source symbol when unmapped (counted missing)."""
    programs = {}
    for name in ORDER:
        rows = mapping[mapping.program.eq(name) & ~mapping.duplicate_in_program]
        genes = [s if isinstance(s, str) else src for s, src in zip(rows.symbol, rows.source_symbol)]
        programs[name] = {"genes": genes, "readout_role": ROLES[name]}
    return programs


def summarize_mapping(mapping):
    out = mapping.groupby(["program", "readout_role"]).agg(
        n_source=("source_symbol", "size"),
        n_exact=("match_type", lambda s: int(s.eq("exact").sum())),
        n_hgnc_previous=("match_type", lambda s: int(s.eq("hgnc_previous").sum())),
        n_hgnc_alias=("match_type", lambda s: int(s.eq("hgnc_alias").sum())),
        n_missing=("match_type", lambda s: int(s.isin(["missing", "ambiguous"]).sum())),
        n_duplicate=("duplicate_in_program", "sum"),
        missing_genes=("source_symbol", lambda s: ";".join(sorted(
            s[mapping.loc[s.index, "match_type"].isin(["missing", "ambiguous"])])))).reset_index()
    out["n_mapped"] = out.n_exact + out.n_hgnc_previous + out.n_hgnc_alias - out.n_duplicate
    out["coverage"] = (out.n_exact + out.n_hgnc_previous + out.n_hgnc_alias) / out.n_source
    return out.set_index("program").loc[ORDER].reset_index()


def target_overlap(programs, target_symbols):
    """(target, program) pairs where the perturbed gene is a program member."""
    rows = [dict(target=t, program=name, target_symbol=s)
            for t, s in target_symbols.items() for name, info in programs.items() if s in info["genes"]]
    return pd.DataFrame(rows, columns=["target", "program", "target_symbol"])


def unversioned(gene_id):
    return re.sub(r"\.\d+(_PAR_Y)?$", "", str(gene_id))
