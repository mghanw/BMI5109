"""Fetch public reference files into data/ref/ and record provenance in data/ref/SOURCES.md.

Usage:
    python scripts/fetch_refs.py [--download] [--gtf-check data/processed/rna_mtx/rows.txt]

Without --download only existing files are checksummed. --gtf-check reports, for every GTF in
data/ref/gencode/, how many matrix gene IDs match (versioned and unversioned).
"""
import argparse
import datetime as dt
import gzip
import hashlib
import re
import urllib.request
from pathlib import Path

REF = Path("data/ref")
REPO_COMMIT = "bafdfff28a6fc36854a5defe32440abc98eaabc0"
REPO_FILES = [
    "README.md",
    "Dual_omics_processing_scripts/barcode_files/full_1st_lig_barcodes_96.txt",
    "Dual_omics_processing_scripts/barcode_files/full_2nd_lig_barcodes_384.txt",
    "Dual_omics_processing_scripts/barcode_files/dualKD_perturbfate_sgRNA_info_df.txt",
    "bulk_screen_processing/perturbfate_dual_KD_library_with_NTC.csv",
    "Dual_omics_cell_calling/dual_omics_combine_multi_modal_cells_functions.R",
    "Dual_omics_processing_scripts/RNA_modality/perturbFate_dual_omics_RNA.main.sh",
    "Dual_omics_processing_scripts/RNA_modality/gene_count_processing_sciRNAseq_exon_intron.R",
    "Dual_omics_processing_scripts/example_files/example_RNA_sample_id.txt",
    "Dual_omics_processing_scripts/sgRNA_modality/perturbfate_2lig_sgRNA_sc_counting.10xCS.58bp_ver.py",
]
SOURCES = [
    ("GEO series SOFT (sample annotations)",
     "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291147/soft/GSE291147_family.soft.gz",
     REF / "GSE291147_family.soft.gz"),
    ("Paper full text, open access (Europe PMC PMC13233327; doi 10.1038/s41586-026-10367-0)",
     "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13233327/fullTextXML",
     REF / "PMC13233327_fulltext.xml"),
    ("GENCODE v38 primary-assembly GTF (named in the authors' RNA pipeline)",
     "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz",
     REF / "gencode" / "gencode.v38.primary_assembly.annotation.gtf.gz"),
] + [(f"PerturbFate code repository @ {REPO_COMMIT[:7]}",
      f"https://raw.githubusercontent.com/JunyueCaoLab/PerturbFate/{REPO_COMMIT}/{f}",
      REF / f"perturbfate_repo_{REPO_COMMIT}" / f.replace("/", "_")) for f in REPO_FILES]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def gtf_gene_ids(path):
    with gzip.open(path, "rt") as f:
        return {re.search(r'gene_id "([^"]+)"', line).group(1)
                for line in f if not line.startswith("#") and line.split("\t", 3)[2] == "gene"}


def gtf_check(rows_path):
    ids = set(Path(rows_path).read_text().split())
    unversion = lambda s: {re.sub(r"\.\d+", "", x) for x in s}
    for gtf in sorted((REF / "gencode").glob("*.gtf.gz")):
        G = gtf_gene_ids(gtf)
        print(f"{gtf.name}: GTF genes {len(G)}, matrix {len(ids)}, versioned match {len(ids & G)}, "
              f"unversioned match {len(unversion(ids) & unversion(G))}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--gtf-check")
    args = parser.parse_args()
    lines = ["# Reference sources", "", "Public, read-only files used by scripts/prepare_data.py. "
             "Data only: never executed.", "",
             "| Description | URL | File | Accessed | Bytes | SHA256 |", "|---|---|---|---|---|---|"]
    for desc, url, path in SOURCES:
        if not path.exists():
            if not args.download:
                print(f"missing: {path}")
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            urllib.request.urlretrieve(url, path)
        accessed = dt.datetime.fromtimestamp(path.stat().st_mtime).date().isoformat()
        lines.append(f"| {desc} | {url} | {path.relative_to(REF).as_posix()} | {accessed} | "
                     f"{path.stat().st_size} | {sha256(path)} |")
    (REF / "SOURCES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {REF / 'SOURCES.md'}")
    if args.gtf_check:
        gtf_check(args.gtf_check)


if __name__ == "__main__":
    main()
