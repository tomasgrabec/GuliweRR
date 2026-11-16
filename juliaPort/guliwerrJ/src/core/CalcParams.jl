# src/core/CalcParams.jl
#
# Calculation parameters for the Ritz-Rayleigh solver: Legendre integration matrices and the
# sampled Legendre basis. Port of `guliwerr/core/@CalcParams/CalcParams.m` (a MATLAB *value*
# class — every `Set*` method returns a new object rather than mutating `obj` in place).
#
# Translation pattern: a value class whose setters
# already return `obj = obj.SetX(...)` maps onto an **immutable Julia `struct`** plus plain
# setter *functions* that return a new value — `CalcParams` never has a mutating method, by
# construction, so there is no way to hand back an object whose `legOrder` and `legDiag`/
# `legNondiag`/`polyLeg` have drifted out of sync (the guarantee MATLAB's
# `RecomputeLegBasis` enforced only by being "the sole owner of that derivation", i.e. by
# convention/discipline, not by the type system).
#
# Plain top-level definitions, not `module CalcParams` — see MatlabCompat.jl's header
# for why: attaching a module docstring above
# `module CalcParams` while it also contains `struct CalcParams` throws `MethodError: no
# method matching doc!(...)` at load time (Julia's doc system resolves the shadowed binding
# to the type, not the module, once the struct is defined), and separately, `using .CalcParams`
# one level up cannot override the auto-bound module name with the struct's — bare
# `CalcParams(...)` fails with `objects of type Module are not callable`. Both problems
# disappear once there is no per-file submodule to collide with in the first place; this file
# is `include`d directly into `CoreTypes`, and `legendre_matrix_diag_nondiag`/
# `basis_decomposition_2d` from `Legendre.jl` are already in scope as sibling definitions in
# that same module (Legendre.jl is `include`d first — see `CoreTypes.jl`).

export CalcParams, set_rel_depth, set_real_thickness, set_leg_order

# -----------------------------------------------------------------------------------------
# The type
# -----------------------------------------------------------------------------------------

"""
    CalcParams(; relDepth = 20.0, legOrder = 50, realThickness = NaN)

Legendre integration matrices and sampled basis for the Ritz-Rayleigh solver.
Port of `guliwerr/core/@CalcParams/CalcParams.m`.

# Fields
- `relDepth::Float64` — domain depth relative to the wavelength, which is fixed at 1
  throughout the solvers (`|k_par| = 1`).
- `legOrder::Int` — max order of the Legendre polynomials.
- `realThickness::Float64` — real thickness, same units `Structure.ApplyWavelength` is
  called with; `NaN` means unset (MATLAB's own default/sentinel, unchanged here).
- `legDiag::Matrix{Float64}` — `D_mn = ∫ P̃_m' P̃_n'`, `legOrder × legOrder`.
- `legNondiag::Matrix{Float64}` — `N_mn = ∫ P̃_m' P̃_n`, `legOrder × legOrder`.
- `basisVecL::Int` — length of the full basis vector, `3 * legOrder`.
- `polyLeg::Matrix{Float64}` — basis sampled on the depth mesh, `legOrder × nMesh`
  (`nMesh = length(-1:0.01:1) = 201`, since the constructor always uses mesh spacing
  `0.01` — see [`Legendre.basis_decomposition_2d`](@ref)'s docstring on why that must not be
  unified with that function's own `0.02` default).

**Immutable: every derived field (`legDiag`, `legNondiag`, `basisVecL`, `polyLeg`) is
rebuilt in this constructor**, so a `CalcParams` value can never hold a mesh stale for its
own `legOrder` — the invariant `CalcParams.RecomputeLegBasis` documented itself as "the
sole owner of". Because the *only* way to build one is through this keyword method (the
seven-field positional constructor Julia generates automatically is not part of the public
API — nothing stops calling it directly, but nothing here does, and doing so bypasses the
derivation), there is no analogue of MATLAB's `SetAccess = protected`/private-method
discipline to maintain: the type itself enforces it.

See [`legendre_matrix_diag_nondiag`](@ref) and [`basis_decomposition_2d`](@ref) (both in
`Legendre.jl`, same module) for the derivation.

Throws `ArgumentError` if `legOrder < 1`, or if `relDepth` / `realThickness` is not a
positive finite scalar or `NaN` (MATLAB validates the same:
`CalcParams.checkLegOrder`, `CalcParams:SetRelDepth|SetRealThickness:badInput`).
"""
struct CalcParams
    relDepth::Float64
    legOrder::Int
    realThickness::Float64
    legDiag::Matrix{Float64}
    legNondiag::Matrix{Float64}
    basisVecL::Int
    polyLeg::Matrix{Float64}
end

"""
    _check_positive_or_nan(x, methodName, argName)

MATLAB `CalcParams.checkPositiveOrNaN`: a positive finite scalar, or `NaN` (unset;
a film's `relDepth` placeholder until `ApplyWavelength` derives it).
"""
function _check_positive_or_nan(x::Real, methodName::String, argName::String)
    (isnan(x) || (isfinite(x) && x > 0)) ||
        throw(ArgumentError("CalcParams:$methodName:badInput -- $argName must be a positive " *
                             "finite scalar, or NaN (unset); got $x"))
    return nothing
end

function CalcParams(; relDepth::Real = 20.0, legOrder::Integer = 50,
                       realThickness::Real = NaN)
    _check_positive_or_nan(relDepth, "SetRelDepth", "relDepth")
    _check_positive_or_nan(realThickness, "SetRealThickness", "realThickness")
    legOrder >= 1 || throw(ArgumentError("CalcParams: legOrder must be a positive integer, got $legOrder"))
    d, n = legendre_matrix_diag_nondiag(legOrder)
    return CalcParams(Float64(relDepth), Int(legOrder), Float64(realThickness),
                       d, n, 3 * legOrder,
                       basis_decomposition_2d(legOrder, 0.01))
end

# -----------------------------------------------------------------------------------------
# Setters — return a new value, exactly as MATLAB's `obj = obj.SetX(...)` did
# -----------------------------------------------------------------------------------------

"""
    set_rel_depth(cp::CalcParams, v::Real) -> CalcParams

Return a copy of `cp` with a new relative depth. Port of `CalcParams.SetRelDepth`. Nothing
derived (`legDiag`, `legNondiag`, `basisVecL`, `polyLeg`) depends on `relDepth`, so those
fields are carried over unchanged — no recomputation, just as `SetRelDepth` never called
`RecomputeLegBasis`.
"""
function set_rel_depth(cp::CalcParams, v::Real)
    _check_positive_or_nan(v, "SetRelDepth", "relDepth")
    return CalcParams(Float64(v), cp.legOrder, cp.realThickness,
                      cp.legDiag, cp.legNondiag, cp.basisVecL, cp.polyLeg)
end

"""
    set_real_thickness(cp::CalcParams, v::Real) -> CalcParams

Return a copy of `cp` with a new real thickness. Port of `CalcParams.SetRealThickness`; like
[`set_rel_depth`](@ref), nothing derived depends on this field.
"""
function set_real_thickness(cp::CalcParams, v::Real)
    _check_positive_or_nan(v, "SetRealThickness", "realThickness")
    return CalcParams(cp.relDepth, cp.legOrder, Float64(v),
                      cp.legDiag, cp.legNondiag, cp.basisVecL, cp.polyLeg)
end

"""
    set_leg_order(cp::CalcParams, legO::Integer) -> CalcParams

Return a copy of `cp` at a new Legendre order. Port of `CalcParams.SetLegOrder`. Goes
through the keyword constructor (not the raw 7-field one), so **every** derived field is
rebuilt from scratch — the `RecomputeLegBasis` guarantee — and an invalid `legO` throws the
same `ArgumentError` the constructor does. `cp` itself is left untouched: this returns a new
value, it does not and cannot mutate `cp` (`CalcParams` is immutable).
"""
set_leg_order(cp::CalcParams, legO::Integer) =
    CalcParams(; relDepth = cp.relDepth, legOrder = legO,
                  realThickness = cp.realThickness)
