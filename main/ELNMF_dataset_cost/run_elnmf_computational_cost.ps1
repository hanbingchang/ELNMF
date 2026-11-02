param(
    [string[]]$OnlyDatasets = @()
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
$outputRoot = Join-Path $repoRoot "results/ELNMF_computational_cost"
$scriptRoot = $PSScriptRoot
$logRoot = Join-Path $outputRoot "logs"
New-Item -ItemType Directory -Force -Path $outputRoot, $logRoot | Out-Null

Add-Type @"
using System;
using System.Runtime.InteropServices;
public class MemoryStatusNative {
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

$memoryStatus = New-Object MemoryStatusNative+MEMORYSTATUSEX
if (-not [MemoryStatusNative]::GlobalMemoryStatusEx($memoryStatus)) {
    throw "Unable to read physical memory information."
}
$totalPhysicalBytes = [double]$memoryStatus.ullTotalPhys

# Set ELNMF_RSCRIPT to select the historical R version (4.5.1).
$rscript = Resolve-Program $env:ELNMF_RSCRIPT "Rscript"
$datasets = @(
    [pscustomobject]@{ Dataset = "MVC"; Script = "MVC_benchmark.R" },
    [pscustomobject]@{ Dataset = "HBC"; Script = "HBC_benchmark.R" },
    [pscustomobject]@{ Dataset = "DLPFC 151672"; Script = "DLPFC_151672_benchmark.R" }
)

if ($OnlyDatasets.Count -gt 0) {
    $datasets = @($datasets | Where-Object { $OnlyDatasets -contains $_.Dataset })
    if ($datasets.Count -eq 0) {
        throw "No matching datasets were selected."
    }
}

$csvPath = Join-Path $outputRoot "ELNMF_runtime_peak_memory.csv"
$results = @()
if (Test-Path -LiteralPath $csvPath) {
    $results = @(Import-Csv -LiteralPath $csvPath)
}

foreach ($item in $datasets) {
    $safeName = ($item.Dataset -replace '[^A-Za-z0-9]+', '_').Trim('_')
    $benchmarkScript = Join-Path $scriptRoot $item.Script
    if (-not (Test-Path -LiteralPath $benchmarkScript)) { throw "Missing script: $benchmarkScript" }

    $stdoutPath = Join-Path $logRoot ($safeName + "_stdout.log")
    $stderrPath = Join-Path $logRoot ($safeName + "_stderr.log")
    Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue

    Write-Host ("START {0} at {1}" -f $item.Dataset, (Get-Date -Format "yyyy-MM-dd HH:mm:ss"))
    $benchmarkStart = Get-Date
    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
    $process = Start-Process -FilePath $rscript `
        -ArgumentList @("--vanilla", ('"' + $benchmarkScript + '"')) `
        -WorkingDirectory $repoRoot `
        -RedirectStandardOutput $stdoutPath `
        -RedirectStandardError $stderrPath `
        -WindowStyle Hidden `
        -PassThru

    $peakWorkingSet = 0L
    $currentWorkingSet = 0L
    $nextReport = 30.0
    while (-not $process.HasExited) {
        try {
            $activeRProcesses = @(Get-BenchmarkProcesses $process.Id)
            $currentWorkingSet = [long](($activeRProcesses | Measure-Object -Property WorkingSet64 -Sum).Sum)
            if ($currentWorkingSet -gt $peakWorkingSet) {
                $peakWorkingSet = $currentWorkingSet
            }
        } catch {
        }
        if ($stopwatch.Elapsed.TotalSeconds -ge $nextReport) {
            Write-Host ("RUNNING {0}: elapsed={1:N1} min, current={2:N2} GB, peak={3:N2} GB" -f `
                $item.Dataset,
                $stopwatch.Elapsed.TotalMinutes,
                ($currentWorkingSet / 1GB),
                ($peakWorkingSet / 1GB))
            $nextReport += 30.0
        }
        Start-Sleep -Milliseconds 500
    }
    $process.WaitForExit()
    $stopwatch.Stop()

    $logText = ""
    if (Test-Path -LiteralPath $stdoutPath) {
        $logText = Get-Content -LiteralPath $stdoutPath -Raw
    }
    $spots = $null
    $genes = $null
    $dimMatch = [regex]::Match($logText, '(?:高变异基因表达矩阵维度|过滤后X维度).*?(\d+)\s*[×x]\s*(\d+)')
    if ($dimMatch.Success) {
        $genes = [int]$dimMatch.Groups[1].Value
        $spots = [int]$dimMatch.Groups[2].Value
    }
    $ari = $null
    $nmi = $null
    $metricMatches = [regex]::Matches($logText, 'ARI=([0-9.]+),\s*NMI=([0-9.]+)')
    if ($metricMatches.Count -gt 0) {
        $lastMetric = $metricMatches[$metricMatches.Count - 1]
        $ari = [double]$lastMetric.Groups[1].Value
        $nmi = [double]$lastMetric.Groups[2].Value
    }

    $status = if ($process.ExitCode -eq 0) { "Completed" } else { "Failed" }
    $results = @($results | Where-Object { $_.Dataset -ne $item.Dataset })
    $results += [pscustomobject]@{
        Dataset = $item.Dataset
        Run = 1
        Spots = $spots
        Genes = $genes
        Runtime_seconds = [math]::Round($stopwatch.Elapsed.TotalSeconds, 3)
        Runtime_minutes = [math]::Round($stopwatch.Elapsed.TotalMinutes, 3)
        Peak_RAM_GB = [math]::Round($peakWorkingSet / 1GB, 3)
        Total_physical_RAM_GB = [math]::Round($totalPhysicalBytes / 1GB, 3)
        Peak_RAM_percent = [math]::Round(100.0 * $peakWorkingSet / $totalPhysicalBytes, 2)
        ARI = $ari
        NMI = $nmi
        Exit_code = $process.ExitCode
        Status = $status
        Measurement = "Fresh R process; WorkingSet64 sampled every 0.5 s; complete pipeline including preprocessing and saved intermediate V matrices"
    }
    $results | Export-Csv -LiteralPath $csvPath -NoTypeInformation -Encoding UTF8
    Write-Host ("FINISH {0}: status={1}, runtime={2:N2} min, peak={3:N2} GB ({4:N2}%)" -f `
        $item.Dataset,
        $status,
        $stopwatch.Elapsed.TotalMinutes,
        ($peakWorkingSet / 1GB),
        (100.0 * $peakWorkingSet / $totalPhysicalBytes))

    if ($process.ExitCode -ne 0) {
        Write-Host ("See error log: {0}" -f $stderrPath)
    }
    [System.GC]::Collect()
    Start-Sleep -Seconds 3
}

$results | Export-Csv -LiteralPath $csvPath -NoTypeInformation -Encoding UTF8
Write-Host ("RESULTS {0}" -f $csvPath)
