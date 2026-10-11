param([string]$Python = 'python')
$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
& $Python -X utf8 (Join-Path $here 'audit_cache_safety.py')
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host 'Review safety_audit_results.json: BLOCKED means the candidate must remain OFF.'
