"""Build the external program gene sets for Aim 1 scoring.

Usage:
    python -I scripts/build_programs.py [--out config/programs.json]

Reads the supplementary tables recorded in data/ref/SOURCES.md, maps symbols to GENCODE v38
(exact, then unique HGNC previous/alias symbols), and writes config/programs.json (public gene
lists with provenance) and results/summary/aim1_program_mapping.csv. Stops if the Fufa program
does not have 45 source genes or any program maps below 80% of its source genes.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from icpdag import gene_sets as gs  # noqa: E402

REF = ROOT / "data" / "ref"
PATHS = {"fufa": REF / "fufa2019" / "13072_2019_297_MOESM2_ESM.xlsx",
         "kanai": REF / "kanai2024" / "cancers-16-00852-s001_Table_S3.xlsx",
         "tsoi": REF / "tsoi2018" / "1-s2.0-S1535610818301223-mmc4.xlsx"}
HGNC = REF / "hgnc" / "hgnc_complete_set.txt"
PROVENANCE = {
    gs.FUFA: dict(paper="Fufa et al. 2019, Epigenetics Chromatin 12:50 (PMC6688322)", source="fufa",
                  sheet="Table S2", rule="Log2 Fold Change > 0 AND 'SOX10 binding under AZD6244' marked X or X*"),
    gs.KANAI_UP: dict(paper="Kanai et al. 2024, Cancers 16:852 (PMC10930532)", source="kanai",
                      sheet="Tab 6-YAPTAZ Signature", rule="Total Up > 0"),
    gs.KANAI_DOWN: dict(paper="Kanai et al. 2024, Cancers 16:852 (PMC10930532)", source="kanai",
                        sheet="Tab 6-YAPTAZ Signature", rule="Total Down > 0"),
    **{name: dict(paper="Tsoi et al. 2018, Cancer Cell 33:890 (PMC5953834)", source="tsoi", sheet="Table S3",
                  rule=f"Signature label is exactly '{sub}' (labels shared by two adjacent subtypes are not "
                       f"used; this reproduces the earlier report)")
       for name, sub in gs.TSOI.items()},
}
EXPECTED = {gs.FUFA: 45, gs.KANAI_UP: 80, gs.KANAI_DOWN: 52}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=str(ROOT / "config" / "programs.json"))
    parser.add_argument("--summary", default=str(ROOT / "results" / "summary" / "aim1_program_mapping.csv"))
    parser.add_argument("--rna", default=str(ROOT / "data" / "processed" / "rna_singlets.h5ad"))
    args = parser.parse_args(argv)
    import anndata as ad

    var = ad.read_h5ad(args.rna, backed="r").var
    genes = var[gs.scoring_genes(var)].reset_index().rename(columns={var.index.name or "index": "gene_id"})
    hgnc = pd.read_csv(HGNC, sep="\t", dtype=str,
                       usecols=["symbol", "prev_symbol", "alias_symbol", "ensembl_gene_id"])
    sources = gs.read_sources(PATHS)
    for name, n in EXPECTED.items():
        if len(sources[name]) != n:
            raise SystemExit(f"STOP: {name} has {len(sources[name])} source genes, expected {n}")
    mapping = gs.mapping_table(sources, genes, hgnc)
    summary = gs.summarize_mapping(mapping)
    print(summary.drop(columns="missing_genes").to_string(index=False))
    for row in summary.itertuples():
        if row.missing_genes:
            print(f"  {row.program} missing: {row.missing_genes}")
    if summary.coverage.lt(0.8).any():
        raise SystemExit("STOP: a program maps below 80% of its source genes")
    renamed = mapping[mapping.match_type.str.startswith("hgnc")]
    print("renamed:", "; ".join(f"{r.source_symbol}->{r.symbol}" for r in renamed.drop_duplicates("source_symbol").itertuples()))
    spec = {
        "description": "External expression programs (public supplementary tables), mapped to GENCODE v38 "
                       "symbols of data/processed/rna_singlets.h5ad. Scoring order fixes background draws.",
        "order": gs.ORDER,
        "hgnc_sha256": sha256(HGNC),
        "programs": {},
    }
    for name in gs.ORDER:
        info = PROVENANCE[name]
        rows = mapping[mapping.program.eq(name)]
        spec["programs"][name] = {
            "readout_role": gs.ROLES[name], "paper": info["paper"], "sheet": info["sheet"], "rule": info["rule"],
            "source_file": PATHS[info["source"]].relative_to(ROOT).as_posix(),
            "source_sha256": sha256(PATHS[info["source"]]),
            "genes": rows[["source_symbol", "symbol", "gene_id", "match_type", "duplicate_in_program"]]
            .astype(object).where(rows[["source_symbol", "symbol", "gene_id", "match_type",
                                        "duplicate_in_program"]].notna(), None).to_dict("records"),
        }
    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(spec, f, indent=1)
        f.write("\n")
    summary.to_csv(args.summary, index=False)
    print(f"wrote {args.out} and {args.summary}")


if __name__ == "__main__":
    main()
