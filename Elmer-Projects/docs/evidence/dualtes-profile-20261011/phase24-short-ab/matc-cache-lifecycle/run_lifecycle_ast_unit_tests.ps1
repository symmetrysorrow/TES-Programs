param(
  [string]$MatcSource = 'C:\Users\Nao\AppData\Local\Temp\phase24-03-build-source\matc\src',
  [string]$MatcDll = 'C:\Users\Nao\AppData\Local\Temp\phase24-03-release-build\matc\src\libmatc.dll',
  [string]$Gcc = 'G:\Elmer-Migration\msys64\ucrt64\bin\gcc.exe'
)
$ErrorActionPreference = 'Stop'
$test = Join-Path $PSScriptRoot 'matc_ast_purity_test.c'
$exe = Join-Path $PSScriptRoot 'matc_ast_purity_test.exe'
$log = Join-Path $PSScriptRoot 'matc_ast_unit.log'
$env:PATH = "$(Split-Path $Gcc);$(Split-Path $MatcDll);$env:PATH"
& $Gcc $test "-I$MatcSource" "-L$(Split-Path $MatcDll)" -lmatc -o $exe
if ($LASTEXITCODE) { throw 'MATC AST unit-test compile failed.' }
& $exe 2>&1 | Tee-Object -FilePath $log
if ($LASTEXITCODE) { throw 'MATC AST unit tests failed.' }
Write-Output 'MATC_AST_UNIT_TEST_PASS'

