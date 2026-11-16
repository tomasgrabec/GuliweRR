# =========================================================================================
# Structure.jl
#
# Port of guliwerr/core/@Structure/Structure.m: the layer stack (Material + CalcParams per
# layer) that every solver in the codebase is ultimately built on.
#
# NO per-file module (see CoreTypes.jl's header):
# `struct Structure` lives directly in `CoreTypes`, exactly like `Material`/`ElastConsts`.
#
# MUTABLE vs IMMUTABLE: `classdef Structure` has no `< handle` -- a MATLAB *value* class,
# and its own header comment says so explicitly: "Value class: AddLayer, RemoveLayer and
# ModifyLayer return the modified object." All four mutating-looking methods
# (`AddLayer`, `RemoveLayer`, `ModifyLayer`, `ApplyWavelength`) follow the
# `obj = obj.SetX(...)` rebuild-and-return idiom and never write to `obj`'s properties without
# also returning `obj`. Nothing in the MATLAB code assigns to a Structure instance's
# properties from outside Structure.m (the other same-named assignments, e.g.
# `CalcAngSweepHS.mater`, belong to GUI-form structs and other classes' own properties).
# Conclusion: **immutable struct**, same choice as ElastConsts/Material/CalcParams.
#
# DERIVED FIELDS BECOME FUNCTIONS --
# THE ONE STRUCTURAL DEVIATION FROM A LITERAL FIELD-FOR-FIELD PORT:
# MATLAB stores `nLay` and `strBasisVecL` as properties that `AddLayer`/`RemoveLayer`/
# `ModifyLayer` each update by hand in lockstep with the mutation -- three separate places
# that must all agree, and a fourth import point
# (`ApplyWavelength`, which calls `ModifyLayer` and so inherits the bookkeeping for free) that
# would silently go stale if a future edit ever bypassed `ModifyLayer`. In Julia both become
# plain functions computed on demand from `mater`/`calcParams` directly, so there is nothing
# to keep in sync and nothing that can go stale:
#   nlay(s)           = length(s.mater)                                    -- was `obj.nLay`
#   str_basis_vec_l(s) = sum(cp.basisVecL for cp in s.calcParams)          -- was `obj.strBasisVecL`
# Named exactly `nlay`/`str_basis_vec_l` (lower_snake_case, the naming for a NEW port-only
# entity with no direct MATLAB *property* counterpart -- there is no
# `obj.nlay` to grep for and match casing against, unlike e.g. `AddLayer`, which stays
# PascalCase because it names an existing MATLAB *method*).
#
# NAMING NOTE (methods, as opposed to the two derived quantities above) -- same deliberate
# convention as Material.jl/CalcParams.jl: `AddLayer`, `RemoveLayer`, `ModifyLayer`,
# `ApplyWavelength`, `isMirrorSymmetric`, `isSagittalMirrorSymmetric` keep their MATLAB method names verbatim (receiver
# `obj` becomes an explicit first argument) for greppability against the original.
#
# CONSTRUCTOR SHAPE -- keyword-only, not positional-with-defaults (a positional outer
# constructor with the same arity as the struct's field count can break package
# precompilation outright, as for `Material`):
# `Structure` has 4 fields (`name`, `nMaxModes`, `mater`, `calcParams`); a positional outer
# constructor mirroring MATLAB's `Structure(materIn,calcParamsIn,name,nMaxModes)` one-for-one
# would have the SAME arity as the struct's own auto-generated positional inner constructor
# and hit exactly that precompilation bug. The keyword constructor below uses `materIn`/
# `calcParamsIn` -- not `mater`/`calcParams` -- as its keyword names, deliberately: those
# keywords accept a single first-layer Material/CalcParams (mirroring MATLAB's own constructor
# argument names, Structure.m), not the Vector each field actually stores, so reusing the
# field names would misleadingly suggest they take vectors.
# =========================================================================================

"""
    Structure

The GuliweRR computational domain -- a stack of layers, each described by a `Material` and a
`CalcParams`.

Port of `guliwerr/core/@Structure/Structure.m`. **Immutable struct** -- see the file-level
comment in `Structure.jl` for why.

Fields (kept camelCase, to match the MATLAB properties one-for-one):

- `name`       :: structure name
- `nMaxModes`  :: maximum number of computed modes
- `mater`      :: `Vector{Material}`, layer 1 = bottom … layer end = top
- `calcParams` :: `Vector{CalcParams}`, same layer order as `mater`

`nLay` and `strBasisVecL` are NOT fields here -- they are the functions [`nlay`](@ref) and
[`str_basis_vec_l`](@ref) (see this file's header comment for why: computed on demand, they
cannot go stale the way MATLAB's hand-updated properties could).

Construct via the keyword-optional form [`Structure(; materIn, calcParamsIn, name,
nMaxModes)`](@ref) (every argument optional, mirroring `Structure.m`'s `nargin`-guarded
constructor); do not call the raw 4-field positional form directly outside the package -- it
performs none of the constructor's "adding `materIn` goes through `AddLayer`" logic.
"""
struct Structure
    name::String
    nMaxModes::Int
    mater::Vector{Material}
    calcParams::Vector{CalcParams}
end

"""
    nlay(s::Structure) -> Int

Number of layers in the structure. Port of the MATLAB property `Structure.nLay`, turned into
a function -- see `Structure.jl`'s header comment.
"""
nlay(s::Structure) = length(s.mater)

"""
    str_basis_vec_l(s::Structure) -> Int

Sum of `basisVecL` across all layers. Port of the MATLAB property `Structure.strBasisVecL`,
turned into a function -- see `Structure.jl`'s header comment. `0` for an empty structure
(matches `Structure.m`'s own `strBasisVecL = 0` default).
"""
str_basis_vec_l(s::Structure) = sum(cp.basisVecL for cp in s.calcParams; init = 0)

"""
    Structure(; materIn=nothing, calcParamsIn=nothing, name=nothing, nMaxModes=nothing) -> Structure

Construct a `Structure`, optionally adding a first layer.

INPUT (all optional, all keyword -- see this file's CONSTRUCTOR SHAPE note):
`materIn` -- `Material` for the first layer; omitted/`nothing` adds no layer, exactly as
MATLAB's `isempty(materIn)` guard does (Structure.m: "An empty materIn is what
suppresses AddLayer, not the argument count"); `calcParamsIn` -- `CalcParams` for that first
layer (a default `CalcParams()` is used if `materIn` is given but `calcParamsIn` is not);
`name` -- structure name; `nMaxModes` -- maximum number of modes to compute (MATLAB default
500, applied by the struct's own default construction below, not by this keyword's default).

Port of `guliwerr/core/@Structure/Structure.m`.
"""
function Structure(; materIn::Union{Material,Nothing} = nothing,
                      calcParamsIn::Union{CalcParams,Nothing} = nothing,
                      name::Union{AbstractString,Nothing} = nothing,
                      nMaxModes::Union{Integer,Nothing} = nothing)
    obj = Structure("str 1", 500, Material[], CalcParams[])
    if materIn !== nothing
        cp0 = calcParamsIn === nothing ? CalcParams() : calcParamsIn
        obj = AddLayer(obj, materIn, cp0)
    end
    if name !== nothing
        obj = Structure(String(name), obj.nMaxModes, obj.mater, obj.calcParams)
    end
    if nMaxModes !== nothing
        obj = Structure(obj.name, Int(nMaxModes), obj.mater, obj.calcParams)
    end
    return obj
end

# -----------------------------------------------------------------------------------------
# AddLayer / RemoveLayer / ModifyLayer (Structure.m) -- each returns a NEW Structure,
# since the struct is immutable; none mutates `obj`, and none mutates `obj`'s existing
# `mater`/`calcParams` VECTORS in place either (new vectors are built every time) -- both
# matter for value semantics, since a `Vector` field only stores a pointer.
# -----------------------------------------------------------------------------------------

"""
    AddLayer(obj::Structure, mater0=nothing, calcParams0=nothing) -> Structure

Append a layer to the top of the structure, returning a new [`Structure`](@ref).

INPUT (both optional): `mater0` -- `Material` for the new layer (default `Material()` if
omitted); `calcParams0` -- `CalcParams` for the new layer (default `CalcParams()` if
omitted).

Port of `guliwerr/core/@Structure/Structure.m`.
"""
function AddLayer(obj::Structure, mater0::Union{Material,Nothing} = nothing,
                   calcParams0::Union{CalcParams,Nothing} = nothing)
    m0  = mater0 === nothing ? Material() : mater0
    cp0 = calcParams0 === nothing ? CalcParams() : calcParams0
    return Structure(obj.name, obj.nMaxModes, [obj.mater; m0], [obj.calcParams; cp0])
end

"""
    RemoveLayer(obj::Structure, ind) -> Structure

Remove a layer from the structure by index (1 = bottom, `nlay(obj)` = top), returning a new
[`Structure`](@ref).

Throws `ArgumentError` if `ind` is not a real numeric scalar, or if it is out of range
(mirroring `Structure.m`'s two distinct error identifiers, `Structure:RemoveLayer:badInput`
and `Structure:RemoveLayer:outOfRange`, as two distinct message prefixes -- Julia has no
error-identifier system).

Port of `guliwerr/core/@Structure/Structure.m`.
"""
function RemoveLayer(obj::Structure, ind)
    (ind isa Real) ||
        throw(ArgumentError("Structure:RemoveLayer:badInput -- Layer index must be a " *
            "numeric scalar, got $(typeof(ind))"))
    n = nlay(obj)
    (ind < 1 || ind > n) &&
        throw(ArgumentError("Structure:RemoveLayer:outOfRange -- Layer index $ind is out " *
            "of range (structure has $n layers)"))

    idx = Int(ind)
    newMater      = copy(obj.mater);      deleteat!(newMater, idx)
    newCalcParams = copy(obj.calcParams); deleteat!(newCalcParams, idx)
    return Structure(obj.name, obj.nMaxModes, newMater, newCalcParams)
end

"""
    ModifyLayer(obj::Structure, indLay, mater0=nothing, calcParams0=nothing) -> Structure

Replace the material and/or calcParams of an existing layer, returning a new
[`Structure`](@ref).

INPUT: `indLay` -- layer index to modify; `nothing` -> top-most layer (mirroring MATLAB's
`isempty(indLay)` guard); `mater0` -- new `Material`; `nothing` -> keep existing (mirroring
MATLAB's `isempty(mater0)` guard); `calcParams0` -- new `CalcParams`; `nothing` -> keep
existing.

Port of `guliwerr/core/@Structure/Structure.m`.
"""
function ModifyLayer(obj::Structure, indLay::Union{Integer,Nothing},
                      mater0::Union{Material,Nothing} = nothing,
                      calcParams0::Union{CalcParams,Nothing} = nothing)
    n = nlay(obj)
    idx = indLay === nothing ? n : indLay
    (idx < 1 || idx > n) &&
        throw(ArgumentError("Structure:ModifyLayer:outOfRange -- Layer index $idx is out " *
            "of range (structure has $n layers)"))

    m0  = mater0 === nothing ? obj.mater[idx] : mater0
    cp0 = calcParams0 === nothing ? obj.calcParams[idx] : calcParams0

    newMater      = copy(obj.mater);      newMater[idx]      = m0
    newCalcParams = copy(obj.calcParams); newCalcParams[idx] = cp0
    return Structure(obj.name, obj.nMaxModes, newMater, newCalcParams)
end

"""
    ApplyWavelength(obj::Structure, wavelengthReal::Real) -> Structure

Derive `relDepth = realThickness / wavelengthReal` for every FILM layer (layers `2:nlay`),
returning a new [`Structure`](@ref). The substrate (layer 1) is skipped -- see
`Structure.m`'s comment (`legOrder` stays fixed across a wavelength sweep, so a
derived substrate `relDepth` would grow arbitrarily large at a long-wavelength point and
undermesh the substrate's Legendre basis; set it deliberately with [`set_rel_depth`](@ref) if
needed).

Throws `ArgumentError` if `wavelengthReal` is not finite and positive, or if any film layer
has no `realThickness` set (`isnan`).

Port of `guliwerr/core/@Structure/Structure.m`.
"""
function ApplyWavelength(obj::Structure, wavelengthReal::Real)
    (isfinite(wavelengthReal) && wavelengthReal > 0) ||
        throw(ArgumentError("Structure:ApplyWavelength:expectedPositive -- wavelengthReal " *
            "must be a finite, positive real scalar, got $wavelengthReal"))

    result = obj
    for ilay in 2:nlay(obj)   # layer 1 = substrate, skipped
        cp0 = result.calcParams[ilay]
        isnan(cp0.realThickness) &&
            throw(ArgumentError("Structure:ApplyWavelength:missingRealThickness -- Layer " *
                "$ilay has no realThickness set -- cannot derive relDepth."))
        cp0 = set_rel_depth(cp0, cp0.realThickness / wavelengthReal)
        result = ModifyLayer(result, ilay, nothing, cp0)
    end
    return result
end

"""
    isMirrorSymmetric(obj::Structure, tol::Real=1e-10) -> Bool

Is the plate invariant under the mid-plane mirror `z -> D-z`? True iff, for every `i` and
`j = nlay(obj)+1-i` (the middle layer of an odd stack paired with itself), `relDepth` and
density match and layer `j`'s standard-frame tensor equals the MIRROR IMAGE of layer `i`'s:
`C_j = M C_i M`, `M = diag(1,1,-1)`. `tol` is relative.

This is the exact condition for every mode to have a definite mid-plane parity. Pairing a
layer with its mirror image, not an identical copy, makes one test cover both requirements:
the stack must be symmetric AND each crystal must look the same from either face -- a single
cubic (111) layer fails, (001)/(110) and isotropic layers pass; a stack of two mirror-twin
layers passes. Angle-independent (the `x3` mirror commutes with in-plane rotations).
Conservative: a plate symmetric only after merging identical adjacent layers returns `false`.
Exact parity does not make the Lamb/TH split exact -- see [`isSagittalMirrorSymmetric`](@ref).

Port of `Structure.isMirrorSymmetric`.
"""
function isMirrorSymmetric(obj::Structure, tol::Real = 1e-10)
    signZ = _mirrorSign(3)
    n = nlay(obj)
    for i in 1:cld(n, 2)
        j = n + 1 - i
        Ci = GetStandardElastTensor(obj.mater[i], 0.0)
        Cj = GetStandardElastTensor(obj.mater[j], 0.0)
        dI, dJ = obj.calcParams[i].relDepth, obj.calcParams[j].relDepth
        rI, rJ = obj.mater[i].rho, obj.mater[j].rho
        if abs(dI - dJ) > tol * max(abs(dI), abs(dJ)) ||
                abs(rI - rJ) > tol * max(abs(rI), abs(rJ)) ||
                maximum(abs.(Cj .- signZ .* Ci)) > tol * maximum(abs.(Ci))
            return false
        end
    end
    return true
end

"""
    isSagittalMirrorSymmetric(obj::Structure, angleDeg::Real, tol::Real=1e-10) -> Bool

Is the plate invariant under the sagittal mirror `x2 -> -x2` at this propagation angle? True
iff every layer's standard-frame tensor at `angleDeg` satisfies `C = M C M`,
`M = diag(1,-1,1)` (`tol` relative). Where it holds, SH (`u2`) decouples from the sagittal
(`u1`, `u3`) motion, so the Lamb/TH split of the mode labels is exact; with
[`isMirrorSymmetric`](@ref) as well, all four labels are. Angle-dependent (Cu(001): `true`
along <100> and <110>, `false` at 30 deg).

Port of `Structure.isSagittalMirrorSymmetric`.
"""
function isSagittalMirrorSymmetric(obj::Structure, angleDeg::Real, tol::Real = 1e-10)
    signY = _mirrorSign(2)
    for i in 1:nlay(obj)
        C = GetStandardElastTensor(obj.mater[i], angleDeg)
        maximum(abs.(C .- signY .* C)) > tol * maximum(abs.(C)) && return false
    end
    return true
end

# Sign pattern of a mirror normal to axis `ax` on a 4th-order tensor:
# (M C M)_ijkl = s_ijkl C_ijkl, s_ijkl = (-1)^(how many of i,j,k,l equal ax).
function _mirrorSign(ax::Integer)
    v = ones(3); v[ax] = -1.0
    return [v[i] * v[j] * v[k] * v[l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
end

export Structure, nlay, str_basis_vec_l, AddLayer, RemoveLayer, ModifyLayer,
       ApplyWavelength, isMirrorSymmetric, isSagittalMirrorSymmetric
