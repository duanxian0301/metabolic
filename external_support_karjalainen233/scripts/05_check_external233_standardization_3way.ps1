$gwasDir = "D:\metabolic\233\GWAS"
$qcDir = "D:\metabolic\233\QC"
$logDir = "D:\metabolic\233\logs"

Write-Host "GWAS txt count:"
(Get-ChildItem -Path $gwasDir -Filter *.txt -ErrorAction SilentlyContinue | Measure-Object).Count

Write-Host "`nQC summary count:"
(Get-ChildItem -Path $qcDir -Filter *.qc_summary.tsv -ErrorAction SilentlyContinue | Measure-Object).Count

Write-Host "`nRecent logs:"
Get-ChildItem -Path $logDir -Filter *.log -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 6 Name, Length, LastWriteTime |
    Format-Table -AutoSize

Write-Host "`nRunning python processes:"
Get-Process python -ErrorAction SilentlyContinue |
    Select-Object Id, CPU, StartTime, Path |
    Sort-Object StartTime |
    Format-Table -AutoSize
