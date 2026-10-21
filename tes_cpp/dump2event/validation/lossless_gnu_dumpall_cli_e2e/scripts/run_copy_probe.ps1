$ErrorActionPreference = 'Stop'
$compiler = (Get-Command gfortran -ErrorAction Stop).Source
$root = Split-Path -Parent $PSScriptRoot
$src = Join-Path $PSScriptRoot 'copy_formatted_records_probe.f90'
$work = Join-Path $env:TEMP 'gnu_lossless_copy_probe'
New-Item -ItemType Directory -Force -Path $work | Out-Null
$exe = Join-Path $work 'copy_records.exe'
& $compiler -std=f2008 -Wall -Wextra -O2 $src -o $exe
if ($LASTEXITCODE -ne 0) { throw "gfortran failed: $LASTEXITCODE" }
$cases = @(
  @{ Name='short'; Length=17; FinalLf=$true },
  @{ Name='exact80'; Length=80; FinalLf=$true },
  @{ Name='long81'; Length=81; FinalLf=$true },
  @{ Name='hundreds'; Length=600; FinalLf=$true },
  @{ Name='multi_k'; Length=9000; FinalLf=$true },
  @{ Name='empty'; Length=0; FinalLf=$true },
  @{ Name='no_final_lf'; Length=81; FinalLf=$false }
)
$rows = foreach ($case in $cases) {
  $input = Join-Path $work ($case.Name + '.in')
  $output = Join-Path $work ($case.Name + '.out')
  $payload = 'X' * $case.Length
  $suffix = if ($case.FinalLf) { "`n" } else { '' }
  [IO.File]::WriteAllText($input, $payload + $suffix, [Text.Encoding]::ASCII)
  & $exe $input $output
  if ($LASTEXITCODE -ne 0) { throw "copy failed for $($case.Name): $LASTEXITCODE" }
  $actual = [Text.Encoding]::ASCII.GetString([IO.File]::ReadAllBytes($output)).TrimEnd("`r", "`n")
  $inputBytes = [IO.File]::ReadAllBytes($input)
  $outputBytes = [IO.File]::ReadAllBytes($output)
  [pscustomobject]@{
    case = $case.Name
    input_bytes = $inputBytes.Length
    output_bytes = $outputBytes.Length
    payload_equal = ($actual -ceq $payload)
    input_has_final_lf = ($inputBytes.Length -gt 0 -and $inputBytes[$inputBytes.Length-1] -eq 10)
    output_has_final_lf = ($outputBytes.Length -gt 0 -and $outputBytes[$outputBytes.Length-1] -eq 10)
  }
}
$rows | Format-Table -AutoSize
