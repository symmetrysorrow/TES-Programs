set pagination off
set print pretty on
set breakpoint pending on
break CRSMatrix.F90:1174
commands
silent
printf "\n=== PHASE24 CRITICAL ENTRY ===\n"
info sharedlibrary
info address '.gomp_critical_user_phase24_position_cache_setup'
p/x &gmtr
p gmtr%numberofrows
p gmtr%rows
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
x/8bx '.gomp_critical_user_phase24_position_cache_setup'
bt 5
continue
end
run
