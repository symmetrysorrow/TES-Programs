# Phase22 native HeatSolve timing contract

The project-side parser accepts the following marker from the rebuilt native
`HeatSolve` module:

```text
HeatSolveWallBreakdown: step=<n> iter=<n> total_wall_s=<s> element_traversal_wall_s=<s> local_fem_wall_s=<s> local_stiffness_wall_s=<s> local_mass_wall_s=<s> nonlinear_material_wall_s=<s> global_insertion_wall_s=<s> rhs_assembly_wall_s=<s> matrix_conversion_wall_s=<s> boundary_assembly_wall_s=<s>
```

The current Phase22 native integration source has wall timers around
`DiffuseConvective(Gen)Compose`, `DefaultUpdateMass`/
`Default1stOrderTime`, `DefaultUpdateEquations`, the bulk traversal residual,
and the boundary assembly loop. Local stiffness formation is not independently
separable inside `Compose`; RHS and sparse-format conversion remain in the
documented residual until lower-level hooks are available.

This source change is in the sibling native integration worktree at
`D:/github/TES-Programs/tools/elmer-phase20-native-integration` and remains
unbuilt because the installed module files are not available in the current
WSL include path. Until the rebuilt solver is installed, the production report
keeps `assembly_wall` uninstrumented and places that time in
`unclassified_wall`; it never substitutes the legacy CPU assembly counter.
