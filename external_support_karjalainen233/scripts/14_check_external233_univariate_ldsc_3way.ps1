$ErrorActionPreference = "Stop"

$statusDir = "D:\metabolic\233\ldsc_univariate\status"
$resultsDir = "D:\metabolic\233\ldsc_univariate\results"
$sumstatsDir = "D:\metabolic\233\sumstats"

$sumstats = (Get-ChildItem $sumstatsDir -Filter "*.sumstats.gz" -File -ErrorAction SilentlyContinue | Measure-Object).Count
$logs = (Get-ChildItem $resultsDir -Filter "*_h2.log" -File -ErrorAction SilentlyContinue | Measure-Object).Count
$status = Get-ChildItem $statusDir -Filter "*.tsv" -File -ErrorAction SilentlyContinue | Sort-Object LastWriteTime -Descending
$py = Get-Process python -ErrorAction SilentlyContinue | Select-Object Id,StartTime,CPU | Sort-Object StartTime -Descending

Write-Output "sumstats_gz=$sumstats"
Write-Output "h2_logs=$logs"
Write-Output ""
Write-Output "Recent status files:"
$status | Select-Object -First 6 Name,LastWriteTime,Length | Format-Table -AutoSize
Write-Output ""
Write-Output "Active python processes:"
$py | Format-Table -AutoSize
