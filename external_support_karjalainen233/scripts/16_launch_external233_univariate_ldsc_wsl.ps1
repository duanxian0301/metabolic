param(
    [int]$Groups = 8,
    [string]$AccessionsFile = ""
)

$ErrorActionPreference = "Stop"

$script = "/mnt/d/codex/metabolic_repo/external_support_karjalainen233/scripts/15_run_external233_univariate_ldsc_wsl.py"
$manifest = "D:\metabolic\233\ldsc_univariate\external233_ldsc_manifest.tsv"
$logDir = "D:\metabolic\233\ldsc_univariate\launcher_logs_wsl"

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
for ($i = 1; $i -le $Groups; $i++) {
    $start = (($i - 1) * $groupSize) + 1
    $end = [math]::Min($i * $groupSize, $rows)
    if ($start -gt $rows) { continue }
    $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $log = Join-Path $logDir ("wsl_grp{0}_{1}.log" -f $i, $stamp)
    if ([string]::IsNullOrWhiteSpace($AccessionsFile)) {
        $bashCmd = "/usr/bin/python3 $script --start $start --end $end"
    } else {
        $accessionsWsl = $AccessionsFile.Replace('\','/').Replace('D:','/mnt/d')
        $bashCmd = "/usr/bin/python3 $script --accessions-file '$accessionsWsl' --start $start --end $end"
    }
    Start-Process -FilePath "wsl.exe" -ArgumentList "-e","bash","-lc",$bashCmd -RedirectStandardOutput $log -RedirectStandardError ($log + ".err") -WindowStyle Hidden | Out-Null
    Write-Output ("Launched wsl_grp{0}: rows {1}-{2} -> {3}" -f $i, $start, $end, $log)
}
