# =========================================================================================
# RitzSolverSAW.jl -- port of guliwerr/solvers/ritz/RitzSolverSAW.m
#
# The half-space (and layered half-space) leg of the Ritz family. The solver only SOLVES and
# reports three per-mode diagnostics; existence/classification is RitzScore's / GreenScore's
# job.
#
# Shape: immutable struct, `RitzSolver` abstract supertype, shared fields duplicated (see
# RitzSolver.jl's header). `withDiagnostics` is a `Bool` (MATLAB validates 0/1 at run time;
# Julia's signature does it).
# =========================================================================================

using LinearAlgebra

using ...CoreTypes: Structure, CalcParams, Material, AddLayer, nlay

export RitzSolverSAW, RitzDiagnostics, modeDiagnostics, calcWithGamma

"""
    RitzDiagnostics

Per-mode diagnostics of a solved [`RitzSolverSAW`](@ref) (MATLAB: the `diagnostics` struct,
same field names). Empty vectors until `calc(s, true)` / [`modeDiagnostics`](@ref) runs.

- `freqDisp` -- `exp(-10 ((v - vg)/v)^2)`: 1 for a non-dispersive mode (a true half-space
  surface wave), -> 0 for a dispersive one
- `oopSurfDispl` -- `|u3|` at the top surface (x3 = surface normal)
- `depthInsensitivity` -- the solver's `depthInsensitivity` minus its minimum: large for
  modes whose velocity barely changes with the artificial boundary depth (substrate-only for
  a layered structure), i.e. modes confined near the surface
"""
struct RitzDiagnostics
    freqDisp::Vector{Float64}
    oopSurfDispl::Vector{Float64}
    depthInsensitivity::Vector{Float64}
end
RitzDiagnostics() = RitzDiagnostics(Float64[], Float64[], Float64[])

Base.:(==)(a::RitzDiagnostics, b::RitzDiagnostics) =
    a.freqDisp == b.freqDisp && a.oopSurfDispl == b.oopSurfDispl &&
    a.depthInsensitivity == b.depthInsensitivity

"""
    RitzSolverSAW

Ritz-Rayleigh solver for the half space / layered half space. Port of
`guliwerr/solvers/ritz/RitzSolverSAW.m`. The first 13 fields are [`RitzSolver`](@ref)'s;
own field: `diagnostics::RitzDiagnostics`.

```julia
s = RitzSolverSAW(structure = str, angleDeg = 30.0)
s = calc(s)            # solve + diagnostics (MATLAB default withDiagnostics = true)
s = calc(s, false)     # solve only
```
"""
struct RitzSolverSAW <: RitzSolver
    structure::Structure
    angle::Float64
    nCompModes::Int
    freqs::Vector{Float64}
    vels::Vector{Float64}
    velsGroup::Vector{Float64}
    eigVecs_con::Matrix{ComplexF64}
    eigVals::Vector{Float64}
    nullConstrMat::Matrix{Float64}
    basisVecL_con::Int
    dVel_dDepth::Vector{Float64}
    depthInsensitivity::Vector{Float64}
    diagnostics::RitzDiagnostics
end

"""
    RitzSolverSAW(; structure = nothing, angleDeg = nothing) -> RitzSolverSAW

Unsolved solver. Both arguments optional (`Structure()`, `0.0`); a non-finite angle is an
`ArgumentError` `RitzSolver:constructor:badAngle` (a wrong-typed structure is a
`MethodError`/`TypeError`, Julia's analogue of `RitzSolver:constructor:badStructure`).
"""
function RitzSolverSAW(; structure::Union{Structure,Nothing} = nothing,
                          angleDeg::Union{Real,Nothing} = nothing)
    return RitzSolverSAW(
        structure === nothing ? Structure() : structure, _check_angle(angleDeg),
        0, Float64[], Float64[], Float64[],
        zeros(ComplexF64, 0, 0), Float64[],
        zeros(Float64, 0, 0), 0,
        Float64[], Float64[], RitzDiagnostics())
end

function _rebuild(s::RitzSolverSAW;
                   structure          = s.structure,
                   angle              = s.angle,
                   nCompModes         = s.nCompModes,
                   freqs              = s.freqs,
                   vels               = s.vels,
                   velsGroup          = s.velsGroup,
                   eigVecs_con        = s.eigVecs_con,
                   eigVals            = s.eigVals,
                   nullConstrMat      = s.nullConstrMat,
                   basisVecL_con      = s.basisVecL_con,
                   dVel_dDepth        = s.dVel_dDepth,
                   depthInsensitivity = s.depthInsensitivity,
                   diagnostics        = s.diagnostics)
    return RitzSolverSAW(structure, angle, nCompModes, freqs, vels, velsGroup,
                          eigVecs_con, eigVals, nullConstrMat, basisVecL_con,
                          dVel_dDepth, depthInsensitivity, diagnostics)
end

"""
    calc(s::RitzSolverSAW, withDiagnostics::Bool = true, precomp = NamedTuple()) -> RitzSolverSAW

Calculate the modes of a half-space (or layered-on-half-space) domain. Port of
`RitzSolverSAW.calc`.

`precomp` may carry (each read lazily, [`precomp_field`](@ref)):
- `nullConstrMat` -- constraint null space ([`buildConstraintMatrix`](@ref) with `SAWBoundary()`)
- `dKdC_cell`, `d2KdCdD_cell` -- from [`derivsByDepth`](@ref), together or not at all
  (a half-supplied pair is an `ArgumentError` from `GammaMatrix`'s `calc`)
- `ritzStruct` -- [`ritzStructurePrecompute`](@ref) of this structure and `nullConstrMat`
  (the angle-independent `M_con` and real-gauge back-map); ignored unless built from the
  same `nullConstrMat` object

Steps: assemble K, dK/dD, M; real gauge if it applies; solve on the Hermitian / real
symmetric path; group velocity (layer-weighted); `depthInsensitivity` from the substrate
block alone for `nLay > 1`; diagnostics when `withDiagnostics`.
"""
function calc(s::RitzSolverSAW, withDiagnostics::Bool = true,
               precomp::NamedTuple = NamedTuple())
    return first(calcWithGamma(s, withDiagnostics, precomp))
end

"""
    calcWithGamma(s::RitzSolverSAW, withDiagnostics = true, precomp = NamedTuple())
        -> (s, gm::GammaMatrix)

[`calc`](@ref), also returning the assembled `GammaMatrix` (`K`, `dK/dD`, `M`) the solve used,
for a caller that evaluates Rayleigh quotients or derivatives on the same matrices (the solver
does not store it, see the [`RitzSolver`](@ref) field table).
"""
function calcWithGamma(s::RitzSolverSAW, withDiagnostics::Bool = true,
                       precomp::NamedTuple = NamedTuple())
    nullConstrMat = precomp_field(precomp, :nullConstrMat,
        () -> buildConstraintMatrix(s.structure.calcParams, SAWBoundary()))
    dKdC_cell    = precomp_field(precomp, :dKdC_cell)
    d2KdCdD_cell = precomp_field(precomp, :d2KdCdD_cell)

    # K, dK/dD, M: a local, not stored (see the RitzSolver field table)
    gm = calc(s.structure, SAWBoundary(), s.angle;
               dKdC_cell = dKdC_cell, d2KdCdD_cell = d2KdCdD_cell)
    out = _rebuild(s; nullConstrMat = nullConstrMat,
                       diagnostics = RitzDiagnostics())

    sp = precomp_field(precomp, :ritzStruct)
    K_con, M_con, gauge = _project_for_solve(out, gm, sp)
    out = solveEigenvalueProblem(out, K_con, M_con, gauge, sp)
    out = calcGroupVelocity(out, gm, nlay(out.structure) > 1)

    if withDiagnostics && out.nCompModes > 0
        out = modeDiagnostics(out)
    end
    return out, gm
end

"""
    modeDiagnostics(s::RitzSolverSAW) -> RitzSolverSAW

Fill `s.diagnostics` (see [`RitzDiagnostics`](@ref)) for every computed mode. Port of
`RitzSolverSAW.modeDiagnostics`. Diagnostics, not a classification. `oopSurfDispl` is the
absolute value of the top-surface u3 row of [`calcDisplacement`](@ref), evaluated alone.
"""
function modeDiagnostics(s::RitzSolverSAW)
    nM = s.nCompModes
    v  = s.vels[1:nM]
    freqDisp = abs.(exp.(-10 .* ((v .- s.velsGroup[1:nM]) ./ max.(v, eps(Float64))) .^ 2))
    oopSurfDispl = abs.(_top_surface_u3(s, nM))
    di = s.depthInsensitivity[1:nM]
    depthInsens = di .- minimum(di)
    return _rebuild(s; diagnostics = RitzDiagnostics(freqDisp, oopSurfDispl, depthInsens))
end

