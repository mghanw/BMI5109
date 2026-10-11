"""Convert R-exported MatrixMarket to a cells x features AnnData, oriented by metadata barcodes.

Usage:
    python scripts/mtx_to_h5ad.py <mtx_dir> <metadata.csv.gz> <out.h5ad> [barcode_column]

If barcode_column is omitted, every metadata column is tried and the one matching the most
matrix row/column names is used. If nothing matches exactly, names are compared with the
"_RNA_" token removed (sgRNA libraries are named without it) and matched cells are renamed
to the metadata cell names.
"""
import sys

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse


def main(mtx_dir, meta_path, out_path, barcode_column=None):
    m = scipy.sparse.csr_matrix(scipy.io.mmread(f"{mtx_dir}/matrix.mtx"))
    rows = pd.read_csv(f"{mtx_dir}/rows.txt", header=None, dtype=str, keep_default_na=False)[0]
    cols = pd.read_csv(f"{mtx_dir}/cols.txt", header=None, dtype=str, keep_default_na=False)[0]
    meta = pd.read_csv(meta_path, dtype=str)
    print("Metadata columns:", list(meta.columns))
    print(f"Matrix shape (rows x cols): {m.shape}")

    candidates = [barcode_column] if barcode_column else list(meta.columns)
    best = None
    for normalize in (False, True):
        key = (lambda x: x.replace("_RNA_", "_")) if normalize else (lambda x: x)
        row_set, col_set = set(rows.map(key)), set(cols.map(key))
        for c in candidates:
            names = set(meta[c].map(key))
            hit = (c, normalize, len(names & col_set), len(names & row_set))
            if best is None or max(hit[2:]) > max(best[2:]):
                best = hit
        if max(best[2:]) > 0:
            break
    column, normalized, in_cols, in_rows = best
    print(f"Barcode column: {column!r} (matches: {in_cols} in columns, {in_rows} in rows; "
          f"'_RNA_' token ignored: {normalized})")
    if max(in_cols, in_rows) == 0:
        sys.exit("No metadata column matches the matrix names; inspect barcode formats manually.")

    if in_cols >= in_rows:
        m, cells, features = m.T.tocsr(), cols, rows  # R stored features x cells
    else:
        cells, features = rows, cols
    if cells.duplicated().any() or features.duplicated().any():
        print("Warning: duplicated cell or feature names; making feature names unique.")
    if normalized:
        # Rename matrix cells to the metadata (RNA) cell names so all files share one ID.
        lookup = {name.replace("_RNA_", "_"): name for name in meta[column]}
        cells = cells.map(lambda x: lookup.get(x.replace("_RNA_", "_"), x))
    adata = ad.AnnData(X=m, obs=pd.DataFrame(index=cells.values), var=pd.DataFrame(index=features.values))
    adata.var_names_make_unique()
    matched = adata.obs_names.isin(meta[column])
    print(f"Cells matched to metadata: {int(matched.sum())} / {adata.n_obs}")
    adata.uns["barcode_column"] = column
    adata.write_h5ad(out_path, compression="gzip")
    print(adata)


if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    main(*sys.argv[1:])
