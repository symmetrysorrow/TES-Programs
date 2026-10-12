$ErrorActionPreference = 'Continue'
$gcc = 'G:\Elmer-Migration\msys64\ucrt64\bin\gcc.exe'
$root = $PSScriptRoot
$probe = Join-Path $root 'asan_capability_probe.c'
$exe = Join-Path $root 'asan_capability_probe.exe'
$log = Join-Path $root 'asan_capability_probe.log'
& $gcc -fsanitize=address $probe -o $exe 2>&1 | Tee-Object -FilePath $log
$rc = $LASTEXITCODE
"gcc_exit_code=$rc" | Add-Content $log
if ($rc -eq 0) { & $exe 2>&1 | Tee-Object -FilePath $log -Append; "probe_exit_code=$LASTEXITCODE" | Add-Content $log }
else { 'ASAN_UNAVAILABLE' | Add-Content $log }
