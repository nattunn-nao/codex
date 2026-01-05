param(
    [string]$LogPath,
    [string]$OutPath
)

if (-not $LogPath -or -not $OutPath) {
    Write-Error "Usage: .\gms_progress.ps1 LOGFILE OUTFILE"
    exit 1
}

$culture = [System.Globalization.CultureInfo]::InvariantCulture

# LOG生成待ち
while (-not (Test-Path $LogPath)) {
    Start-Sleep -Seconds 1
}

# ヘッダ
if (-not (Test-Path $OutPath)) {
    "# step  energy(H)        grad_max        rms" |
        Out-File -FilePath $OutPath -Encoding ascii
}

$seenSteps = @{}
$jobStart = $null
$jobEnd   = $null

while ($true) {

    if (-not (Test-Path $LogPath)) {
        Start-Sleep -Seconds 1
        continue
    }

    $lines = Get-Content -Path $LogPath

    foreach ($line in $lines) {

        # 開始時刻
        if (-not $jobStart -and $line -match 'EXECUTION OF GAMESS BEGUN\s+(\d{2}:\d{2}:\d{2})\s+(\d{2}-[A-Z]{3}-\d{4})') {
            $timeStr = $matches[1]
            $dateStr = $matches[2]
            $dtStr   = "$dateStr $timeStr"
            $jobStart = [datetime]::ParseExact($dtStr, 'dd-MMM-yyyy HH:mm:ss', $culture)
            continue
        }

        # 終了時刻
        if (-not $jobEnd -and $line -match 'EXECUTION OF GAMESS TERMINATED.*\s(\d{2}:\d{2}:\d{2})\s+(\d{2}-[A-Z]{3}-\d{4})') {
            $timeStr = $matches[1]
            $dateStr = $matches[2]
            $dtStr   = "$dateStr $timeStr"
            $jobEnd  = [datetime]::ParseExact($dtStr, 'dd-MMM-yyyy HH:mm:ss', $culture)
            continue
        }

        # NSERCH: ... E= ... GRAD. MAX= ... R.M.S.= ...
        if ($line -match 'NSERCH:\s*([0-9]+)\s+E=\s*([-\d\.]+).*GRAD\. MAX=\s*([0-9\.]+)\s+R\.M\.S\.\s*=\s*([0-9\.]+)') {

            $nserch = [int]$matches[1]
            $E      = $matches[2]
            $gmax   = $matches[3]
            $rms    = $matches[4]

            if ($seenSteps.ContainsKey($nserch)) { continue }

            $stepIndex = $nserch + 1
            $lineOut = ("{0,3}  {1,16}  {2,12}  {3,12}" -f `
                        $stepIndex, $E, $gmax, $rms)
            $lineOut | Out-File -FilePath $OutPath -Append -Encoding ascii
            $seenSteps[$nserch] = $true
            continue
        }
    }

    if ($jobEnd) { break }

    Start-Sleep -Seconds 1
}

# 実行時間の追記
if ($jobStart -and $jobEnd) {
    $wall = New-TimeSpan -Start $jobStart -End $jobEnd

    "# =========================================" | Out-File -FilePath $OutPath -Append -Encoding ascii
    ("# start_time  {0}" -f $jobStart.ToString("yyyy-MM-dd HH:mm:ss")) | Out-File -FilePath $OutPath -Append -Encoding ascii
    ("# end_time    {0}" -f $jobEnd.ToString("yyyy-MM-dd HH:mm:ss"))   | Out-File -FilePath $OutPath -Append -Encoding ascii
    ("# wall_time   {0}" -f $wall.ToString())                          | Out-File -FilePath $OutPath -Append -Encoding ascii
}
