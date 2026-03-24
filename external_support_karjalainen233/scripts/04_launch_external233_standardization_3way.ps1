$manifest = "D:\metabolic\233\manifests\external233_processing_manifest.tsv"
$script = "D:\codex\metabolic_repo\external_support_karjalainen233\scripts\03_standardize_external_gwas_to_txt.py"
$python = "D:\conda\python.exe"
$logDir = "D:\metabolic\233\logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

if (-not (Test-Path $manifest)) {
    throw "Manifest not found: $manifest"
}

$n = (Import-Csv -Delimiter "`t" $manifest).Count
if ($n -le 0) {
    throw "No rows in manifest: $manifest"
}

$chunk = [math]::Ceiling($n / 3.0)
$groups = @(
    @{ Name = "grp1"; Start = 1; End = [math]::Min($chunk, $n) },
    @{ Name = "grp2"; Start = $chunk + 1; End = [math]::Min($chunk * 2, $n) },
    @{ Name = "grp3"; Start = $chunk * 2 + 1; End = $n }
)

foreach ($g in $groups) {
    if ($g.Start -gt $g.End) { continue }
    $log = Join-Path $logDir ("external233_{0}_{1}_{2}.log" -f $g.Name, $g.Start, $g.End)
    $err = Join-Path $logDir ("external233_{0}_{1}_{2}.err.log" -f $g.Name, $g.Start, $g.End)
    Start-Process -WindowStyle Hidden -FilePath $python -ArgumentList @(
        $script,
        "--start", $g.Start,
        "--end", $g.End
    ) -RedirectStandardOutput $log -RedirectStandardError $err | Out-Null
    Start-Sleep -Milliseconds 300
    Write-Host ("Started {0}: {1}-{2}" -f $g.Name, $g.Start, $g.End)
    Write-Host ("Log: {0}" -f $log)
}
