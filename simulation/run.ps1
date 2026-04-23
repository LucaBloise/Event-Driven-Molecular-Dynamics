#!/usr/bin/env pwsh
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$SimulationArgs = @()
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$srcDir = Join-Path $scriptDir "src/main/java"
$binDir = Join-Path $scriptDir "bin"
$mainClass = "ar.edu.itba.sds.tp3.simulation.SimulationMain"

$javacCommand = Get-Command javac -ErrorAction SilentlyContinue
if (-not $javacCommand) {
    $javacCommand = Get-Command javac.exe -ErrorAction SilentlyContinue
}

if (-not $javacCommand) {
    $jdksDir = Join-Path $HOME ".jdks"
    if (-not (Test-Path $jdksDir) -and $env:USERPROFILE) {
        $jdksDir = Join-Path $env:USERPROFILE ".jdks"
    }

    if (Test-Path $jdksDir) {
        $fallbackJavac = Get-ChildItem -Path $jdksDir -Recurse -File -Filter "javac.exe" |
            Select-Object -First 1
        if ($fallbackJavac) {
            $javacCommand = [PSCustomObject]@{ Source = $fallbackJavac.FullName }
        }
    }
}

if (-not $javacCommand) {
    Write-Error "Error: javac was not found in PATH. Install a JDK (17+) and retry."
    exit 1
}

$javacPath = $javacCommand.Source
$javaPath = $null

$javaCommand = Get-Command java -ErrorAction SilentlyContinue
if (-not $javaCommand) {
    $javaCommand = Get-Command java.exe -ErrorAction SilentlyContinue
}

if ($javaCommand) {
    $javaPath = $javaCommand.Source
} elseif ($javacPath.ToLower().EndsWith("javac.exe")) {
    $javaPath = Join-Path (Split-Path -Parent $javacPath) "java.exe"
} else {
    $javaPath = "java"
}

if (-not (Test-Path $binDir)) {
    New-Item -ItemType Directory -Path $binDir | Out-Null
}

$javaFiles = Get-ChildItem -Path $srcDir -Recurse -File -Filter "*.java" |
    Sort-Object FullName |
    Select-Object -ExpandProperty FullName

if (-not $javaFiles -or $javaFiles.Count -eq 0) {
    Write-Error "Error: no Java sources were found under $srcDir"
    exit 1
}

& $javacPath -d $binDir @javaFiles
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}

& $javaPath -cp $binDir $mainClass @SimulationArgs
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}