param(
    [int[]]$Sizes = @(1000, 2000, 3000, 3798, 4500, 5000, 6000),
    [switch]$ForceRerun,
    [double]$MaxWorkingSetFraction = 0.80,
    [double]$MinAvailableGB = 1.5
)

$ErrorActionPreference = "Stop"

# Paths are relative to this repository, not to the author's workstation.
# Track only the launched process and its descendants. Unrelated jobs are excluded.
function Get-BenchmarkProcesses([int]$RootProcessId) {
    if ($script:TrackedRoot -ne $RootProcessId) {
        $script:TrackedRoot = $RootProcessId
        $script:TrackedIds = @($RootProcessId)
    }
    $pending = @($script:TrackedIds)
    while ($pending.Count -gt 0) {
        $newIds = @()
        foreach ($parentId in $pending) {
            $children = @(Get-CimInstance Win32_Process -Filter "ParentProcessId = $parentId" -ErrorAction Stop)
            foreach ($child in $children) {
                if ($script:TrackedIds -notcontains [int]$child.ProcessId) {
                    $script:TrackedIds += [int]$child.ProcessId
                    $newIds += [int]$child.ProcessId
                }
            }
        }
        $pending = $newIds
    }
    return @(Get-Process -Id $script:TrackedIds -ErrorAction SilentlyContinue)
}
# Verify process-tree access before starting a long-running experiment.
$null = Get-CimInstance Win32_Process -Filter "ProcessId = $PID" -ErrorAction Stop

$repoRoot = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
function Resolve-Program([string]$Configured, [string]$Fallback) {
    if ($Configured) {
        if (-not (Test-Path -LiteralPath $Configured -PathType Leaf)) { throw "Executable not found: $Configured" }
        return (Resolve-Path -LiteralPath $Configured).Path
    }
    return (Get-Command $Fallback -ErrorAction Stop).Source
}
$outputRoot = Join-Path $repoRoot "results/HBC_ELNMF_scaling"
$workRoot = Join-Path $outputRoot "work"
$logRoot = Join-Path $outputRoot "logs"
New-Item -ItemType Directory -Force -Path $outputRoot, $workRoot, $logRoot | Out-Null
$benchmarkR = Join-Path $PSScriptRoot "HBC_ELNMF_scaling_benchmark.R"
$rscript = Resolve-Program $env:ELNMF_RSCRIPT "Rscript"
if ($MaxWorkingSetFraction -le 0 -or $MaxWorkingSetFraction -gt 1 -or $MinAvailableGB -lt 0) {
    throw "Invalid memory safety settings"
}
# To reproduce the-spot full-run safety policy: -Sizes -MaxWorkingSetFraction 0.97 -MinAvailableGB 1.5

Add-Type @"
using System;
using System.Runtime.InteropServices;
public class ScalingMemoryStatusNative {
    [StructLayout(LayoutKind.Sequential, CharSet = CharSet.Auto)]
    public class MEMORYSTATUSEX {
        public uint dwLength;
        public uint dwMemoryLoad;
        public ulong ullTotalPhys;
        public ulong ullAvailPhys;
        public ulong ullTotalPageFile;
        public ulong ullAvailPageFile;
        public ulong ullTotalVirtual;
        public ulong ullAvailVirtual;
        public ulong ullAvailExtendedVirtual;
        public MEMORYSTATUSEX() { dwLength = (uint)Marshal.SizeOf(typeof(MEMORYSTATUSEX)); }
    }
    [DllImport("kernel32.dll", CharSet = CharSet.Auto, SetLastError = true)]
    public static extern bool GlobalMemoryStatusEx([In, Out] MEMORYSTATUSEX lpBuffer);
}
"@

$memoryStatus = New-Object ScalingMemoryStatusNative+MEMORYSTATUSEX
if (-not [ScalingMemoryStatusNative]::GlobalMemoryStatusEx($memoryStatus)) {
    throw "Unable to read physical memory information."
}
$totalPhysicalBytes = [double]$memoryStatus.ullTotalPhys
$stopBytes = $MaxWorkingSetFraction * $totalPhysicalBytes

$csvPath = Join-Path $outputRoot "HBC_ELNMF_scaling_benchmark.csv"
$results = @()
if (Test-Path -LiteralPath $csvPath) { $results = @(Import-Csv -LiteralPath $csvPath) }

foreach ($n in $Sizes) {
    $existing = @($results | Where-Object { [int]$_.N_spots -eq $n -and $_.Status -eq "Completed" })
    if (-not $ForceRerun -and $existing.Count -gt 0) {
        Write-Host "SKIP N=$n (already completed)"
        continue
    }

    $runName = "N{0:D5}" -f $n
    $runDir = Join-Path $workRoot $runName
    New-Item -ItemType Directory -Force -Path $runDir | Out-Null
    $startMarker = Join-Path $runDir "TIMED_REGION_STARTED.marker"
    $endMarker = Join-Path $runDir "TIMED_REGION_COMPLETED.marker"
    $errorMarker = Join-Path $runDir "ERROR.txt"
    $runtimePath = Join-Path $runDir "runtime_seconds.txt"
    Remove-Item -LiteralPath $startMarker, $endMarker, $errorMarker, $runtimePath -Force -ErrorAction SilentlyContinue

    $stdoutPath = Join-Path $logRoot ($runName + "_stdout.log")
    $stderrPath = Join-Path $logRoot ($runName + "_stderr.log")
    Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue

    Write-Host ("START N={0} at {1}" -f $n, (Get-Date -Format "yyyy-MM-dd HH:mm:ss"))
    $processLaunch = Get-Date
    $process = Start-Process -FilePath $rscript `
        -ArgumentList @("--vanilla", ('"' + $benchmarkR + '"'), $n, ('"' + $outputRoot + '"')) `
        -WorkingDirectory $repoRoot `
        -RedirectStandardOutput $stdoutPath `
        -RedirectStandardError $stderrPath `
        -WindowStyle Hidden `
        -PassThru

    $peakWorkingSet = 0L
    $currentWorkingSet = 0L
    $timedStart = $null
    $terminatedForMemory = $false
    $nextReport = 30.0
    $wallClock = [System.Diagnostics.Stopwatch]::StartNew()

    while (-not $process.HasExited) {
        if ($null -eq $timedStart -and (Test-Path -LiteralPath $startMarker)) {
            $timedStart = (Get-Item -LiteralPath $startMarker).LastWriteTime
            Write-Host ("TIMED REGION N={0} started at {1}" -f $n, $timedStart)
        }

        if ($null -ne $timedStart) {
            try {
                $activeR = @(Get-BenchmarkProcesses $process.Id)
                $sum = ($activeR | Measure-Object -Property WorkingSet64 -Sum).Sum
                $currentWorkingSet = if ($null -eq $sum) { 0L } else { [long]$sum }
                if ($currentWorkingSet -gt $peakWorkingSet) { $peakWorkingSet = $currentWorkingSet }

                [ScalingMemoryStatusNative]::GlobalMemoryStatusEx($memoryStatus) | Out-Null
                if ($currentWorkingSet -ge $stopBytes -or [double]$memoryStatus.ullAvailPhys -le $MinAvailableGB * 1GB) {
                    $terminatedForMemory = $true
                    Write-Host ("MEMORY LIMIT N={0}: {1:N2} GB ({2:N2}%). Stopping this run." -f `
                        $n, ($currentWorkingSet / 1GB), (100.0 * $currentWorkingSet / $totalPhysicalBytes))
                    foreach ($rp in $activeR) {
                        Stop-Process -Id $rp.Id -Force -ErrorAction SilentlyContinue
                    }
                    break
                }
            } catch {
            }
        }

        if ($wallClock.Elapsed.TotalSeconds -ge $nextReport) {
            $phase = if ($null -eq $timedStart) { "input preparation" } else { "timed method" }
            Write-Host ("RUNNING N={0}: phase={1}, wall={2:N1} min, current={3:N2} GB, peak={4:N2} GB" -f `
                $n, $phase, $wallClock.Elapsed.TotalMinutes, ($currentWorkingSet / 1GB), ($peakWorkingSet / 1GB))
            $nextReport += 30.0
        }
        Start-Sleep -Milliseconds 500
        $process.Refresh()
    }

    try { $process.WaitForExit() } catch {}
    $wallClock.Stop()

    $runtimeSeconds = $null
    if ((Test-Path -LiteralPath $runtimePath) -and -not $terminatedForMemory) {
        $runtimeSeconds = [double]((Get-Content -LiteralPath $runtimePath -Raw).Trim())
    } elseif ($null -ne $timedStart) {
        $runtimeSeconds = ((Get-Date) - $timedStart).TotalSeconds
    }

    $retained = $null
    $summaryPath = Join-Path $runDir "run_summary.csv"
    if (Test-Path -LiteralPath $summaryPath) {
        $summaryRow = Import-Csv -LiteralPath $summaryPath | Select-Object -First 1
        $retained = $summaryRow.Retained_clusterings
    }

    $exitCode = if ($terminatedForMemory) { -9 } else { $process.ExitCode }
    $status = if ($terminatedForMemory) {
        "Terminated_for_memory_safety"
    } elseif ($exitCode -eq 0 -and (Test-Path -LiteralPath $endMarker)) {
        "Completed"
    } else {
        "Failed"
    }

    $results = @($results | Where-Object { [int]$_.N_spots -ne $n })
    $results += [pscustomobject]@{
        Dataset = "HBC"
        N_spots = $n
        Run = 1
        Seed = 1234
        Runtime_sec = if ($null -eq $runtimeSeconds) { $null } else { [math]::Round($runtimeSeconds, 3) }
        Runtime_min = if ($null -eq $runtimeSeconds) { $null } else { [math]::Round($runtimeSeconds / 60.0, 3) }
        Peak_RAM_GB = [math]::Round($peakWorkingSet / 1GB, 3)
        Total_RAM_GB = [math]::Round($totalPhysicalBytes / 1GB, 3)
        Peak_RAM_percent = [math]::Round(100.0 * $peakWorkingSet / $totalPhysicalBytes, 2)
        Max_working_set_fraction = $MaxWorkingSetFraction
        Minimum_available_RAM_GB_setting = $MinAvailableGB
        Status = $status
        Retained_clusterings = $retained
        Input_mode = if ($n -lt 3798) { "stratified_subsample" } elseif ($n -eq 3798) { "original" } else { "original_plus_resampled_spots" }
        Measurement = "Timed preprocessing through final labels; process WorkingSet64 sampled every 0.5 s; one run"
    }
    $results = @($results | Sort-Object { [int]$_.N_spots })
    $results | Export-Csv -LiteralPath $csvPath -NoTypeInformation -Encoding UTF8

    Write-Host ("FINISH N={0}: status={1}, runtime={2}, peak={3:N2} GB ({4:N2}%)" -f `
        $n, $status, $runtimeSeconds, ($peakWorkingSet / 1GB), (100.0 * $peakWorkingSet / $totalPhysicalBytes))

    if ($status -ne "Completed") {
        Write-Host "Stopping the scaling sequence; larger sizes will not be run."
        break
    }
    [System.GC]::Collect()
    Start-Sleep -Seconds 5
}

Write-Host "RESULTS $csvPath"

