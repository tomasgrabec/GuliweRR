"""
    Green

Aggregator for GULIWERR's Green's-function surface-response solver: `GreenSolverTypes`
(the result type), `GreenSextic` (k3 sextic root-finding and the polarisation pre-pass),
and `GreenSolver` (boundary-condition assembly, the velocity loop, and the layered stack).

Nested exactly two levels below `GuliwerrJ` (`GuliwerrJ.Solvers.Green`), the same depth as
`GuliwerrJ.Solvers.Christoffel`/`GuliwerrJ.Solvers.Ritz` — anything here that needs
`CoreTypes` climbs to it via `using ...CoreTypes: ...` (three dots); anything that needs
`Christoffel` (a SIBLING under `Solvers`, not an ancestor) reaches it via
`using ..Christoffel: ...` (two dots) — see `GreenSextic.jl`'s own header for the worked
reasoning. This creates a real load-order constraint: `Solvers.jl` must `include` this file
AFTER `christoffel/Christoffel.jl`.

None of the files below declare their own `module` — one module per top-level folder, not
per file (see `CoreTypes.jl`'s header).

**Name-collision note for `GreenSolver.jl`:** `Ritz` already exports `calc`. If
`Green` needs a `calc` too, use `import ..Ritz: calc` and add methods to the existing
generic function rather than defining a fresh binding — `using .Ritz` + `using .Green` both
exporting a distinct `calc` binding is an unresolvable ambiguity when the bare name is used
(a binding clash, not a dispatch one — different argument types do not help).
"""
module Green

include("GreenSolverTypes.jl")
include("GreenSextic.jl")

export GreenSolver, nvel, has_minus, has_half_space_full_output, has_layered_full_output
export GreenSexticPrepass, green_sextic_prepass, green_layer_setup, green_layers, GreenLayer,
       calc_detw_by_christoffel, solve_k3_allvel, canonicalize_k3,
       check_identical_records_allvel,
       calc_u_allvel, substrate_wave_indices, kvec_cols, outer_cols,
       OMEGA_TOLERANCE, ZERO_TOLERANCE, LIMITING_ABSORPTION_ETA,
       NEAR_ISO_ANGLE_STEP, NEAR_ISO_WARN_THRESH

include("GreenSolver.jl")
export calc_traction_cols, calc_gfull_from_partial_waves, build_bc_matrix,
       green_component, green_gfull, DISP_FIELD_NDEPTH, DISP_FIELD_MAX_DEPTH, GREEN_CHUNK

end # module Green
