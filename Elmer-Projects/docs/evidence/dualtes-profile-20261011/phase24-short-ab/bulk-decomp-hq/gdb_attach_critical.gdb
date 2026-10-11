set pagination off
set breakpoint pending on
set print pretty on
file C:/Users/Nao/AppData/Local/Temp/phase24-debug-runtime/bin/ElmerSolver_mpi.exe
attach 5064
sharedlibrary
break C:/Users/Nao/AppData/Local/Temp/phase24-debug-src/fem/src/CRSMatrix.F90:1174
commands
silent
printf "\\n=== PHASE24 CRITICAL ENTRY ===\\n"
info sharedlibrary
info address '.gomp_critical_user_phase24_position_cache_setup'
p/x &gmtr
p gmtr%numberofrows
p/x gmtr%rows
p/x gmtr%cols
p/x gmtr%values
p n
p ndofs
p elementindex
p/x &indices
p/x &lmtr
p/x indices(1)
p/x indices(2)
p/x indices(3)
p/x indices(4)
bt 5
continue
end
continue
