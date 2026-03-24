param(
    [int]$Groups = 3,
    [string]$AccessionsFile = ""
)

$ErrorActionPreference = "Stop"

$pythonExe = "D:\conda\python.exe"
$script = "D:\codex\metabolic_repo\external_support_karjalainen233\scripts\11_run_external233_univariate_ldsc.py"
$manifest = "D:\metabolic\233\ldsc_univariate\external233_ldsc_manifest.tsv"
$logDir = "D:\metabolic\233\ldsc_univariate\launcher_logs"

New-Item -ItemType Directory -Force -Path $logDir | Out-Null

$rows = if ([string]::IsNullOrWhiteSpace($AccessionsFile)) {
    (Import-Csv -Delimiter "`t" $manifest).Count
} else {
    (Import-Csv -Delimiter "`t" $AccessionsFile).Count
}
if ($rows -lt 1) {
    throw "Input list is empty."
}

$groupSize = [math]::Ceiling($rows / [double]$Groups)
$ranges = @()
for ($i = 1; $i -le $Groups; $i++) {
    $start = (($i - 1) * $groupSize) + 1
    $end = [math]::Min($i * $groupSize, $rows)
    $ranges += @{ Name = ("grp{0}" -f $i); Start = $start; End = $end }
}

foreach ($r in $ranges) {
    if ($r.Start -gt $rows) { continue }
    $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $log = Join-Path $logDir ("{0}_{1}.log" -f $r.Name, $stamp)
    if ([string]::IsNullOrWhiteSpace($AccessionsFile)) {
        $args = @($script, "--start", [string]$r.Start, "--end", [string]$r.End)
    } else {
        $args = @($script, "--accessions-file", $AccessionsFile, "--start", [string]$r.Start, "--end", [string]$r.End)
    }
    Start-Process -FilePath $pythonExe -ArgumentList $args -RedirectStandardOutput $log -RedirectStandardError ($log + ".err") -WindowStyle Hidden | Out-Null
    Write-Output ("Launched {0}: rows {1}-{2} -> {3}" -f $r.Name, $r.Start, $r.End, $log)
}
