param(
  [Parameter(Mandatory=$true)][string]$ElmerSource
)
$ErrorActionPreference = 'Stop'
$expectedCommit = '5a8de867068be0568f09af40fb90fee300dfbede'
$expected = @{
  'fem/src/modules/HeatSolve.F90' = 'fddcaed084563ede6b534dd5bb658d8dbad424fb95f98bb8bd890ba50498e819'
  'fem/src/Lists.F90' = '4f641d128d80672dc43888460bbeb2c7519719cfc3225ccb8ff688d95a2d4dcc'
  'matc/src/parser.c' = '8ddf30c20de8ace2eaddd5ec2f1634dea690a6f755e9d0383394cba9d9bc926a'
  'matc/src/eval.c' = '594fef16510a9291c64c7551d18941c3a9a70504d4ee990703bfb64f224d1c9d'
}
$actualCommit = (git -C $ElmerSource rev-parse HEAD).Trim()
if ($actualCommit -ne $expectedCommit) { throw "Unexpected Elmer revision: $actualCommit" }
foreach ($relative in $expected.Keys) {
  $path = Join-Path $ElmerSource $relative
  $actual = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
  if ($actual -ne $expected[$relative]) { throw "Source hash mismatch: $relative ($actual)" }
  Write-Output "PASS $relative $actual"
}
Write-Output "PASS official-source provenance $actualCommit"
Write-Output 'NOTE the Phase24 candidate patch is not applicable to pristine release-26.1; no patch or build is attempted.'
