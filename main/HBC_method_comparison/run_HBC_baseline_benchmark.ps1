param(
    [string[]]$OnlyMethods = @()
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
$resultRoot = Join-Path $repoRoot "results"
$artifactRoot = Join-Path $resultRoot "HBC_baseline_computational_benchmark"
$scriptRoot = $PSScriptRoot
$logRoot = Join-Path $artifactRoot "logs"
$labelRoot = Join-Path $artifactRoot "labels"
$csvPath = Join-Path $artifactRoot "HBC_computational_benchmark.csv"
New-Item -ItemType Directory -Force -Path $artifactRoot, $logRoot, $labelRoot | Out-Null
# Configure interpreter paths through environment variables; resolve only selected methods.
# Historical versions: CGNMF R 4.5.2; Seurat R 4.5.1. Python methods use separate environments.
$methods = @(
    [pscustomobject]@{
        Method = "CGNMF"
        Executable = $env:CGNMF_RSCRIPT
        Fallback = "Rscript"
        Arguments = @("--vanilla", (Join-Path $scriptRoot "HBC_CGNMF_benchmark.R"))
        WorkDir = $repoRoot
        Device = "CPU"
        ProcessNames = @("Rscript", "Rterm")
    },
    [pscustomobject]@{
        Method = "stHGC"
        Executable = $env:STHGC_PYTHON
        Fallback = "python"
        Arguments = @("-u", (Join-Path $scriptRoot "HBC_stHGC_benchmark.py"))
        WorkDir = $repoRoot
        Device = "CPU"
        ProcessNames = @("python")
    },
    [pscustomobject]@{
        Method = "MNMST"
        Executable = $env:MNMST_PYTHON
        Fallback = "python"
        Arguments = @("-u", (Join-Path $scriptRoot "HBC_MNMST_benchmark.py"))
        WorkDir = $repoRoot
        Device = "CPU"
        ProcessNames = @("python")
    },
    [pscustomobject]@{
        Method = "Seurat"
        Executable = $env:SEURAT_RSCRIPT
        Fallback = "Rscript"
        Arguments = @("--vanilla", (Join-Path $scriptRoot "HBC_Seurat_benchmark.R"))
        WorkDir = $repoRoot
        Device = "CPU"
        ProcessNames = @("Rscript", "Rterm")
    }
)

if ($OnlyMethods.Count -gt 0) {
    $methods = @($methods | Where-Object { $OnlyMethods -contains $_.Method })
    if ($methods.Count -eq 0) { throw "No matching methods selected." }
}

$rows = @()
if (Test-Path -LiteralPath $csvPath) {
    $rows = @(Import-Csv -LiteralPath $csvPath)
}

$oldOutput = $env:HBC_BENCH_OUTPUT
$oldPythonUnbuffered = $env:PYTHONUNBUFFERED
$env:HBC_BENCH_OUTPUT = $labelRoot
$env:PYTHONUNBUFFERED = "1"

try {
    foreach ($item in $methods) {
        $item.Executable = Resolve-Program $item.Executable $item.Fallback

        $stdoutPath = Join-Path $logRoot ($item.Method + "_stdout.log")
        $stderrPath = Join-Path $logRoot ($item.Method + "_stderr.log")
        Remove-Item -LiteralPath $stdoutPath, $stderrPath -Force -ErrorAction SilentlyContinue

        $quotedArgs = @($item.Arguments | ForEach-Object {
            if ($_ -match '\s') { '"' + $_ + '"' } else { $_ }
        })

        Write-Host ("START {0}: preparing HBC input" -f $item.Method)
        $launchTime = Get-Date
        $process = Start-Process -FilePath $item.Executable `
            -ArgumentList $quotedArgs `
            -WorkingDirectory $item.WorkDir `
            -RedirectStandardOutput $stdoutPath `
            -RedirectStandardError $stderrPath `
            -WindowStyle Hidden `
            -PassThru

        $markerSeen = $false
        $endSeen = $false
        $peakWorkingSet = 0L
        $runtimeSec = $null
        $clusterCount = $null
        $lastReport = Get-Date

        while ($true) {
            $activeProcesses = @(Get-BenchmarkProcesses $process.Id)
            if ($activeProcesses.Count -eq 0) { break }

            $logText = if (Test-Path -LiteralPath $stdoutPath) {
                [string](Get-Content -LiteralPath $stdoutPath -Raw -ErrorAction SilentlyContinue)
            } else { "" }
            if ($null -eq $logText) { $logText = "" }

            if (-not $markerSeen -and $logText -match ("BENCHMARK_START\|" + [regex]::Escape($item.Method))) {
                $markerSeen = $true
                Write-Host ("MEASURE {0}: method computation started" -f $item.Method)
            }

            if ($markerSeen) {
                try {
                    $workingSet = [long](($activeProcesses | Measure-Object -Property WorkingSet64 -Sum).Sum)
                    if ($workingSet -gt $peakWorkingSet) { $peakWorkingSet = $workingSet }
                } catch { }
            }

            $endPattern = "BENCHMARK_END\|" + [regex]::Escape($item.Method) + "\|runtime_sec=([0-9.]+)\|clusters=([0-9]+)"
            $endMatch = [regex]::Match($logText, $endPattern)
            if ($endMatch.Success) {
                $runtimeSec = [double]$endMatch.Groups[1].Value
                $clusterCount = [int]$endMatch.Groups[2].Value
                $endSeen = $true
            }

            if (((Get-Date) - $lastReport).TotalSeconds -ge 30) {
                if ($markerSeen) {
                    Write-Host ("RUNNING {0}: peak host RAM={1:N2} GB" -f $item.Method, ($peakWorkingSet / 1GB))
                } else {
                    Write-Host ("PREPARING {0}: waiting for timed-region marker" -f $item.Method)
                }
                $lastReport = Get-Date
            }
            Start-Sleep -Milliseconds 500
        }
        $process.WaitForExit()

        $finalLog = if (Test-Path -LiteralPath $stdoutPath) {
            [string](Get-Content -LiteralPath $stdoutPath -Raw -ErrorAction SilentlyContinue)
        } else { "" }
        if ($null -eq $finalLog) { $finalLog = "" }
        if (-not $endSeen) {
            $endPattern = "BENCHMARK_END\|" + [regex]::Escape($item.Method) + "\|runtime_sec=([0-9.]+)\|clusters=([0-9]+)"
            $endMatch = [regex]::Match($finalLog, $endPattern)
            if ($endMatch.Success) {
                $runtimeSec = [double]$endMatch.Groups[1].Value
                $clusterCount = [int]$endMatch.Groups[2].Value
                $endSeen = $true
            }
        }

        if ($process.ExitCode -ne 0 -or -not $markerSeen -or -not $endSeen) {
            Write-Host ("FAILED {0}: exit={1}; see {2}" -f $item.Method, $process.ExitCode, $stderrPath)
            continue
        }

        $rows = @($rows | Where-Object { $_.Method -ne $item.Method })
        $rows += [pscustomobject]@{
            Method = $item.Method
            Run = 1
            Runtime_sec = [math]::Round($runtimeSec, 3)
            Peak_RAM_GB = [math]::Round($peakWorkingSet / 1GB, 3)
            Peak_GPU_GB = "NA"
            Device = $item.Device
        }
        $rows | Export-Csv -LiteralPath $csvPath -NoTypeInformation -Encoding UTF8
        Write-Host ("FINISH {0}: runtime={1:N3} sec, peak RAM={2:N3} GB, clusters={3}" -f `
            $item.Method, $runtimeSec, ($peakWorkingSet / 1GB), $clusterCount)
        [System.GC]::Collect()
        Start-Sleep -Seconds 3
    }
} finally {
    $env:HBC_BENCH_OUTPUT = $oldOutput
    $env:PYTHONUNBUFFERED = $oldPythonUnbuffered
}

$rows | Export-Csv -LiteralPath $csvPath -NoTypeInformation -Encoding UTF8
Write-Host ("RESULTS {0}" -f $csvPath)
