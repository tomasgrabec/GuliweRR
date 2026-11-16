# src/inverseCalcs/Inv.jl -- the `Inv` launchers
#
# EXPERIMENTAL (as in MATLAB): validated on synthetic data only; the interface may change.
#
# Port of `guliwerr/inverseCalcs/@Inv/Inv.m` -- a MATLAB STATIC NAMESPACE
# (`classdef Inv`, `methods (Static)` only, never instantiated), the same
# pattern as `Christoffel`/`CalcCommon`/`ParallelPolicy`/`Calc`: a namespace becomes free
# functions, no type at all. `Inv` has NO same-named type anywhere in this port, so this one
# COULD have been a nested `module Inv` with no collision risk -- but `Inv.jl` is `include`d
# into the EXISTING `module InverseCalcs`, not wrapped in its own module (one module per
# top-level FOLDER, not per file, is this port's convention: a per-file module named after
# the file's main type would collide with that type name; `inverseCalcs/`'s one module is
# `InverseCalcs`, declared in `InverseCalcs.jl`). So this file, like `Calc.jl`, is the
# THINNEST POSSIBLE LAUNCHER: four one-line aliases, no logic, no defaults of their
# own (a default duplicated here could silently drift from the real constructor's own
# default -- `Inv.m`'s own docblock says exactly this: "Thin varargin passthrough only -- no
# defaults live here, and none should be added").
#
# NAMING: MATLAB's `Inv.measured(...)`/`Inv.association(...)`/`Inv.sensitivity(...)`/
# `Inv.elasticity(...)` become plain top-level functions `measured`/`association`/
# `sensitivity`/`elasticity` in `InverseCalcs` (exactly how `Calc.angSweepHS` became the bare
# function `angSweepHS`) -- no collision with any STRUCT FIELD of the same name
# (`ElasticityInversion.measured`/`.association`/`.sensitivity`, `MeasuredVelocities`'s own
# internal state, etc.): Julia field access (`.field`) and top-level function names are
# different namespaces entirely, so a type having a field called `measured` and a module
# having a function called `measured` do not conflict. `sensitivity(...)` forwards to the
# shared `calc` generic (already imported one level up in `InverseCalcs.jl`); Julia's
# multiple dispatch picks `SensitivityAnalysis.calc`'s method from the argument types exactly
# as it already does for every other caller of that generic -- no special-casing needed here.
# EXPORTED (matching `ForwardCalcs.jl`'s own `Calc.jl`-launcher precedent: `singleHS`,
# `angSweepHS`, etc. are exported plain names, not qualified-only). Note that this export
# surface (`measured`, `association`, `sensitivity`, `elasticity`) has generic-sounding names.

"""
    measured(args...; kwargs...) -> MeasuredVelocities

Thin `Inv` namespace launcher for [`MeasuredVelocities`](@ref). Port of
`guliwerr/inverseCalcs/@Inv/Inv.m`'s `Inv.measured`.
"""
measured(args...; kwargs...) = MeasuredVelocities(args...; kwargs...)

"""
    association(args...; kwargs...) -> ModeAssociation

Thin `Inv` namespace launcher for [`fromSweep`](@ref). Port of `Inv.m`'s `Inv.association`.
"""
association(args...; kwargs...) = fromSweep(args...; kwargs...)

"""
    sensitivity(args...; kwargs...) -> SensitivityAnalysis

Thin `Inv` namespace launcher for `SensitivityAnalysis`'s own [`calc`](@ref) method. Port of
`Inv.m`'s `Inv.sensitivity`.
"""
sensitivity(args...; kwargs...) = calc(args...; kwargs...)

"""
    elasticity(args...; kwargs...) -> ElasticityInversion

Thin `Inv` namespace launcher for [`ElasticityInversion`](@ref). Port of `Inv.m`'s
`Inv.elasticity`.
"""
elasticity(args...; kwargs...) = ElasticityInversion(args...; kwargs...)

export measured, association, sensitivity, elasticity
