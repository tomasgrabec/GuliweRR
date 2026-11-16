# =========================================================================================
# CalcIsofreqSingleLayered.jl -- port of guliwerr/forwardCalcs/CalcIsofreqSingleLayered.m:
# single-angle, single-frequency container for the isoFreq
# Green's-function scan -- fixes one angle and one real reference frequency `omega0`, then
# sweeps wavevector magnitude `kPar` via `GreenSolver.calc`'s `scanMode = :isoFreq`. No Ritz
# and no phase split -- Ritz fixes `|kPar| = 1` by construction.
#
# NO per-file module (see Calc.jl's header for why).
#
# MUTABLE: a `mutable struct`, like the other Calc types, although the MATLAB class is a plain
# value class: every step mutates the object and returns it, so both `calc(x)` and
# `x = calc(x)` work; `copy(x)` is a shallow copy.
#
# DOUBLE-ROTATION INVARIANT -- identical contract to the other single-angle Calc
# types: `structure` is rotated ONCE at construction, and `calc` below passes `0.0`, never
# `cs.angle`, to `GreenSolver`. A rotated anisotropic case (a non-zero, non-symmetric angle) is
# cross-checked against `CalcSingleHS.calcGreen` in the tests.
# =========================================================================================

"""
    CalcIsofreqSingleLayered

Single-angle, single-frequency container for the isoFreq Green's-function scan.

Port of `guliwerr/forwardCalcs/CalcIsofreqSingleLayered.m`. **Mutable struct** (see this file's header). Fixes one
propagation angle and one reference frequency `omega0`, then
sweeps wavevector magnitude `kParRange` via [`calc`](@ref).

# Fields
- `structure::Structure` -- the ROTATED structure for this angle; any `nlay >= 1`
  (half-space or layered -- `GreenSolver` is layered-generic).
- `angle::Float64` -- propagation angle [deg], reference only.
- `omega0::Float64` -- fixed reference frequency; `NaN` (unset) only on the empty/default
  object -- REQUIRED at construction otherwise (no default; a `CalcIsofreqSingleLayered`
  with no frequency reference is meaningless).
- `kParRange::Vector{Float64}` -- wavevector-magnitude axis to sweep; REQUIRED at
  construction otherwise. Must be positive; need not be sorted, but SHOULD be ascending if
  `GreenPeaks.find` will be used on the result (it needs increasing velocity/x-axis).
- `absG::Vector{Float64}` -- [`calc`](@ref)'s result: `abs(G33)` at each `kParRange[i]`,
  indexed by `kParRange`, NOT by the internally-derived velocity range.
"""
mutable struct CalcIsofreqSingleLayered
    structure::Structure
    angle::Float64
    omega0::Float64
    kParRange::Vector{Float64}

    absG::Vector{Float64}
end

"""
    CalcIsofreqSingleLayered(; structure = nothing, angle = nothing, omega0 = nothing,
                                kParRange = nothing) -> CalcIsofreqSingleLayered

Construct `CalcIsofreqSingleLayered` for one propagation angle and one fixed reference
frequency.

Port of `CalcIsofreqSingleLayered.m`. `structure === nothing` (MATLAB's `nargin == 0`)
returns a well-formed empty object (`omega0 = NaN`, `kParRange`/`absG` empty) with no
validation attempted -- matching MATLAB's early `return`.

Otherwise `omega0` and `kParRange` are BOTH required, with no default (MATLAB's own
docstring: "a `CalcIsofreqSingleLayered` with no frequency reference is meaningless"):
- `omega0 === nothing` throws `ArgumentError`
  (`"CalcIsofreqSingleLayered:constructor:missingOmega0"`); otherwise must be real, finite,
  positive.
- `kParRange === nothing` or empty throws `ArgumentError`
  (`"CalcIsofreqSingleLayered:constructor:missingKParRange"`); otherwise must be a real,
  finite, positive vector.

`structure` is then rotated ONCE (see this file's header).

Keyword-only, not positional -- see `CalcSingleHS`'s constructor for why.
"""
function CalcIsofreqSingleLayered(; structure::Union{Structure,Nothing} = nothing,
                                     angle::Union{Real,Nothing} = nothing,
                                     omega0::Union{Real,Nothing} = nothing,
                                     kParRange = nothing)
    if structure === nothing
        return CalcIsofreqSingleLayered(Structure(), 0.0, NaN, Float64[], Float64[])
    end

    ang = angle === nothing ? 0.0 : Float64(angle)

    omega0 === nothing &&
        throw(ArgumentError("CalcIsofreqSingleLayered:constructor:missingOmega0 -- omega0 " *
            "is required (no default) -- a CalcIsofreqSingleLayered with no frequency " *
            "reference is meaningless."))
    (isfinite(omega0) && omega0 > 0) ||
        throw(ArgumentError("CalcIsofreqSingleLayered: omega0 must be a real, finite, " *
            "positive scalar; got $omega0"))

    (kParRange === nothing || isempty(kParRange)) &&
        throw(ArgumentError("CalcIsofreqSingleLayered:constructor:missingKParRange -- " *
            "kParRange is required (no default) -- the wavevector-magnitude axis to sweep."))
    kpr = Float64.(vec(collect(kParRange)))
    (all(isfinite, kpr) && all(>(0), kpr)) ||
        throw(ArgumentError("CalcIsofreqSingleLayered: kParRange must be a real, finite, " *
            "positive vector"))

    structRot = rotateStructure(structure, ang)
    return CalcIsofreqSingleLayered(structRot, ang, Float64(omega0), kpr, Float64[])
end

"""
    _rebuild(cs::CalcIsofreqSingleLayered; <field> = cs.<field>, ...) -> CalcIsofreqSingleLayered

In-place field update (the Calc types are mutable, MATLAB handle semantics). Internal; not exported.
"""
function _rebuild(x::CalcIsofreqSingleLayered; kw...)
    for (k, v) in kw
        setproperty!(x, k, v)
    end
    return x
end

"""
    copy(x::CalcIsofreqSingleLayered)

Shallow copy (MATLAB `copy` of a handle object, one level).
"""
Base.copy(x::CalcIsofreqSingleLayered) = CalcIsofreqSingleLayered((getfield(x, k) for k in fieldnames(CalcIsofreqSingleLayered))...)

# -----------------------------------------------------------------------------------------
# calc -- CalcIsofreqSingleLayered.m
# -----------------------------------------------------------------------------------------

"""
    calc(cs::CalcIsofreqSingleLayered) -> CalcIsofreqSingleLayered

Run the isoFreq Green's-function scan over `cs.kParRange`.

Port of `CalcIsofreqSingleLayered.calc`. Converts
`kParRange` into the `velRange` argument `calc(GreenSolver, ...)`'s `:isoFreq` scan mode
expects: `|kPar| = omega0/v`, so `v = omega0/|kPar|`. `GreenSolver`'s result is in the same
element order as the `velRange` it was given, so `absG` comes back already indexed by
`kParRange`'s own order (regardless of whether `kParRange` is sorted).

`angle = 0.0`, not `cs.angle` -- `cs.structure` is already rotated (see this file's header,
"DOUBLE-ROTATION INVARIANT"). `forceDir = displDir = 3` (G33) with `computeMinus = false`,
since `G_eff` is analytically zero for G33 regardless of `k`-sign (same reasoning as
[`greenSpectrum`](@ref)'s own `computeMinus = false`).
"""
function calc(cs::CalcIsofreqSingleLayered)
    velRangeForGreen = cs.omega0 ./ cs.kParRange
    greenResult = calc(GreenSolver, cs.structure, 0.0, velRangeForGreen;
                        forceDir = 3, displDir = 3, scanMode = :isoFreq,
                        omega0 = cs.omega0, computeMinus = false)
    absG = vec(abs.(greenResult.GFull[3, 3, :]))
    return _rebuild(cs; absG = absG)
end
