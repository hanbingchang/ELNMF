param(
    [int]$N = 6000,
    [double]$MinAvailableGB = 0.75,
    [double]$MaxWorkingSetFraction = 0.97
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
$runDir = Join-Path $outputRoot ("work/N{0:D5}" -f $N)
$vDir = Join-Path $runDir "V_matrices"
$scriptR = Join-Path $PSScriptRoot "HBC_ELNMF_scaling_benchmark.R"
if (-not (Test-Path -LiteralPath $vDir)) { throw "Run the scaling benchmark first to generate the saved V matrices: $vDir" }
New-Item -ItemType Directory -Force -Path (Join-Path $outputRoot "logs") | Out-Null

$vFiles = @(Get-ChildItem -LiteralPath $vDir -Filter "V_lambda*_k*.rds" -File)
if ($vFiles.Count -ne 50) { throw "Expected 50 V matrices, found $($vFiles.Count)" }

Add-Type @"
using System;
using System.Runtime.InteropServices;
public class HBCConsensusMemoryNative {
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

function Get-MemoryStatus {
    $status = New-Object HBC6000ConsensusMemoryNative+MEMORYSTATUSEX
    if (-not [HBC6000ConsensusMemoryNative]::GlobalMemoryStatusEx($status)) {
        throw "Unable to read physical memory information."
    }
    return $status
}

$initialMemory = Get-MemoryStatus
$totalBytes = [double]$initialMemory.ullTotalPhys
$maxWorkingSetBytes = $MaxWorkingSetFraction * $totalBytes
$minAvailableBytes = $MinAvailableGB * 1GB
$initialAvailableGB = [math]::Round([double]$initialMemory.ullAvailPhys / 1GB, 3)

$startMarker = Join-Path $runDir "CONSENSUS_RESUME_STARTED.marker"
$endMarker = Join-Path $runDir "CONSENSUS_RESUME_COMPLETED.marker"
$runtimePath = Join-Path $runDir "consensus_resume_runtime_seconds.txt"
$stdout = Join-Path $outputRoot ("logs\N{0:D5}_consensus_memory_resume_stdout.log" -f $N)
$stderr = Join-Path $outputRoot ("logs\N{0:D5}_consensus_memory_resume_stderr.log" -f $N)
Remove-Item -LiteralPath $startMarker,$endMarker,$runtimePath,$stdout,$stderr -Force -ErrorAction SilentlyContinue

$rscript = Resolve-Program $env:ELNMF_RSCRIPT "Rscript"
$launch = Get-Date
Write-Host ("START CONSENSUS N={0}; initial available={1:N2} GB; safety floor={2:N2} GB" -f $N,$initialAvailableGB,$MinAvailableGB)
$p = Start-Process -FilePath $rscript `
    -ArgumentList @("--vanilla",('"'+$scriptR+'"'),$N,('"'+$outputRoot+'"'), 'consensus') `
    -WorkingDirectory $repoRoot `
    -RedirectStandardOutput $stdout `
    -RedirectStandardError $stderr `
    -WindowStyle Hidden `
    -PassThru

$peakWorkingSet = 0L
$peakPrivate = 0L
$minAvailable = [long]$initialMemory.ullAvailPhys
$nextReport = 30.0
$sw = [Diagnostics.Stopwatch]::StartNew()
$terminated = $false
$terminationReason = ""

while (-not $p.HasExited) {
    $active = @(Get-BenchmarkProcesses $p.Id)
    $wsSum = ($active | Measure-Object WorkingSet64 -Sum).Sum
    $privateSum = ($active | Measure-Object PrivateMemorySize64 -Sum).Sum
    $currentWS = if ($null -eq $wsSum) { 0L } else { [long]$wsSum }
    $currentPrivate = if ($null -eq $privateSum) { 0L } else { [long]$privateSum }
    if ($currentWS -gt $peakWorkingSet) { $peakWorkingSet = $currentWS }
    if ($currentPrivate -gt $peakPrivate) { $peakPrivate = $currentPrivate }

    $memoryNow = Get-MemoryStatus
    $availableNow = [long]$memoryNow.ullAvailPhys
    if ($availableNow -lt $minAvailable) { $minAvailable = $availableNow }

    if ($currentWS -ge $maxWorkingSetBytes -or $availableNow -le $minAvailableBytes) {
        $terminated = $true
        $terminationReason = if ($currentWS -ge $maxWorkingSetBytes) {
            "process_working_set_reached_safety_fraction"
        } else {
            "system_available_ram_reached_safety_floor"
        }
        Write-Host ("MEMORY SAFETY STOP: working set={0:N2} GB, private={1:N2} GB, available={2:N2} GB" -f `
            ($currentWS/1GB),($currentPrivate/1GB),($availableNow/1GB))
        foreach ($rp in $active) { Stop-Process -Id $rp.Id -Force -ErrorAction SilentlyContinue }
        break
    }

    if ($sw.Elapsed.TotalSeconds -ge $nextReport) {
        Write-Host ("CONSENSUS N={0}: elapsed={1:N1} min, working set={2:N2} GB, private={3:N2} GB, peak WS={4:N2} GB, available={5:N2} GB" -f `
            $N,$sw.Elapsed.TotalMinutes,($currentWS/1GB),($currentPrivate/1GB),($peakWorkingSet/1GB),($availableNow/1GB))
        $nextReport += 30
    }
    Start-Sleep -Milliseconds 500
    $p.Refresh()
}

try { $p.WaitForExit() } catch {}
$sw.Stop()
$status = if ($terminated) {
    "Terminated_for_memory_safety"
} elseif ($p.ExitCode -eq 0 -and (Test-Path -LiteralPath $endMarker)) {
    "Completed"
} else {
    "Failed"
}
$consensusSeconds = if (Test-Path -LiteralPath $runtimePath) {
    [double]((Get-Content -LiteralPath $runtimePath -Raw).Trim())
} else {
    $sw.Elapsed.TotalSeconds
}

$result = [pscustomobject]@{
    Dataset = "HBC"
    N_spots = $N
    Run = 1
    Seed = 1234
    Peak_WorkingSet_GB = [math]::Round($peakWorkingSet/1GB,3)
    Peak_WorkingSet_percent = [math]::Round(100*$peakWorkingSet/$totalBytes,2)
    Peak_PrivateMemory_GB = [math]::Round($peakPrivate/1GB,3)
    Total_RAM_GB = [math]::Round($totalBytes/1GB,3)
    Initial_available_RAM_GB = $initialAvailableGB
    Minimum_available_RAM_GB = [math]::Round($minAvailable/1GB,3)
    Consensus_runtime_sec = [math]::Round($consensusSeconds,3)
    Status = $status
    Termination_reason = $terminationReason
    Measurement = "Consensus-only checkpoint rerun using the saved 50 V matrices; process memory sampled every 0.5 s"
}

$resultPath = Join-Path $outputRoot ("HBC_ELNMF_{0}_consensus_memory_resume.csv" -f $N)
$result | Export-Csv -LiteralPath $resultPath -NoTypeInformation -Encoding UTF8
Write-Host ("FINISH: status={0}, peak WS={1:N3} GB ({2:N2}%), peak private={3:N3} GB, minimum available={4:N3} GB" -f `
    $status,($peakWorkingSet/1GB),(100*$peakWorkingSet/$totalBytes),($peakPrivate/1GB),($minAvailable/1GB))
Write-Host "RESULT $resultPath"
