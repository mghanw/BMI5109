# Download the three GSE291147 files used by the project into data\raw.
# Run from the repository root:  powershell -ExecutionPolicy Bypass -File scripts\download_data.ps1
$ErrorActionPreference = "Stop"
$base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291147/suppl"
$files = @(
    "GSE291147_Dual_omics_RNA_gene_count_matrix.RDS",
    "GSE291147_Dual_omics_RNA_cell_metadata.csv.gz",
    "GSE291147_Dual_omics_sgRNA_count_matrix.RDS"
)
New-Item -ItemType Directory -Force data\raw, data\processed, docs | Out-Null
foreach ($f in $files) {
    $dest = "data\raw\$f"
    if (Test-Path $dest) { Write-Host "Exists, skipping: $f"; continue }
    Write-Host "Downloading $f ..."
    curl.exe -L --fail -o $dest "$base/$f"
    if ($LASTEXITCODE -ne 0) { throw "Download failed: $f" }
}
Get-ChildItem data\raw | Format-Table Name, Length
