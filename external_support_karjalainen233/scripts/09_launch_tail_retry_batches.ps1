$python = "D:\conda\python.exe"
$script = "D:\codex\metabolic_repo\external_support_karjalainen233\scripts\03_standardize_external_gwas_to_txt.py"
$retryDir = "D:\metabolic\233\retry_tail"
$gwasDir = "D:\metabolic\233\GWAS"
$logDir = "D:\metabolic\233\logs"

New-Item -ItemType Directory -Force -Path $retryDir | Out-Null
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$batches = Get-ChildItem $retryDir -Filter "external233_tail_batch_*.tsv" | Sort-Object Name
if (-not $batches) {
    throw "No tail retry batch files found in $retryDir"
}

foreach ($batch in $batches) {
    $rows = Import-Csv -Delimiter "`t" $batch.FullName
    foreach ($row in $rows) {
        $tmp = Join-Path $gwasDir ($row.study_accession + ".tmp.txt")
        if (Test-Path $tmp) {
            Remove-Item $tmp -Force -ErrorAction SilentlyContinue
        }
    }
}

foreach ($batch in $batches) {
    $base = [System.IO.Path]::GetFileNameWithoutExtension($batch.Name)
    $log = Join-Path $logDir ("tail_" + $base + ".log")
    $err = Join-Path $logDir ("tail_" + $base + ".err.log")
    Start-Process -WindowStyle Hidden -FilePath $python -ArgumentList @(
        $script,
        "--accessions-file", $batch.FullName
    ) -RedirectStandardOutput $log -RedirectStandardError $err | Out-Null
    Write-Host ("Started tail retry: {0}" -f $batch.Name)
    Start-Sleep -Seconds 2
}
