# =========================================================================================
# Material.jl
#
# Port of guliwerr/core/@Material/Material.m (the type, constructor, four Set* methods, and
# four instance-level rotation wrappers: GetRotMat, GetRotElastTensor, GetRotElastVoigt,
# GetRotDirMat). The other 12 files under core/@Material/ (Euler/DCM/Miller/dirMat
# free-standing conversions, plus the two private LatticeParamSpec helpers) live in
# Orientation.jl, `include`d directly afterwards into the same `CoreTypes` module -- see
# that file's own header comment for why the split, and CoreTypes.jl for the include order
# (Material.jl before Orientation.jl; the two are mutually reachable by plain function calls
# either way, since Julia resolves free-function calls at call time, not include time --
# only STRUCT FIELD TYPES must be defined before use, and neither file's struct depends on
# a type from the other).
#
# NO per-file module (see CoreTypes.jl's header):
# `struct Material` lives directly in `CoreTypes`, exactly like `ElastConsts`.
#
# NAMING NOTE -- same convention as ElastConsts.jl: every function name below
# (`SetElastConsts`, `SetDirMat`, `SetDensity`, `SetName`, `GetRotMat`, `GetRotElastTensor`,
# `GetRotElastVoigt`, `GetRotDirMat`) is kept IDENTICAL to its MATLAB instance-method name,
# with the receiver (`obj`) becoming an explicit first argument (here there is only one
# type, so this is a plain free function, not multiple dispatch; the naming is kept for
# greppability against the MATLAB source).
#
# MUTABLE vs IMMUTABLE: `classdef Material` has no `< handle` -- a MATLAB *value* class, and
# its own header comment says so explicitly ("Value class. Set*/Get* methods take the
# object; other methods are static."). Its four Set* methods all follow the
# `obj = obj.SetX(...)` rebuild-and-return idiom and never mutate `obj` in place. Nothing in
# the MATLAB code assigns to a Material instance's properties outside its own
# constructor/setters (the one superficially similar assignment, `obj.panels{k}.dirMat` in
# `guliwerr/gui/StructureControlGUI.m`, writes to a plain GUI-form struct that holds display
# state, not to a `Material`).
# Conclusion: **immutable struct**, same choice as ElastConsts and CalcParams.
#
# TRANSLATION DECISIONS (language-shape differences, not behaviour changes):
# 1. MATLAB's `elConsts` constructor argument accepts an ElastConsts object, a numeric
#    vector/matrix, OR a cell array `{consts}` / `{consts, symmNo}` (Material.m,
#    SetElastConsts). A MATLAB cell array literal is not idiomatic Julia; the natural
#    equivalent for "an ordered, heterogeneous, fixed-length bundle" is a `Tuple`, so
#    `{consts, symmNo}` becomes `(consts, symmNo)` and `{consts}` becomes `(consts,)`. This
#    is a syntax mapping only -- the three-way dispatch (object / vector-or-matrix / bundle)
#    is preserved exactly, including the "cell with 1 vs >1 elements" branch
#    (SetElastConsts's `if length(elConsts) > 1`).
# 2. MATLAB's `name` property is declared `string` (Material.m); Julia's `String` is the
#    direct equivalent and is what every `SetName`/constructor `name` argument is coerced to.
# 3. The public constructor takes `elConsts`/`rho`/`dirMat`/`name` as KEYWORD arguments
#    (`Material(; elConsts = ..., rho = ..., ...)`), not MATLAB-mirroring positional-optional
#    ones (MATLAB's `nargin > 0 && ~isempty(x)` guards become keyword arguments with
#    defaults). This is not just a style choice here: a POSITIONAL outer
#    constructor with all 4 arguments defaulted collides with the struct's own
#    auto-generated positional inner constructor (also arity 4, one per field) in a
#    Julia/Pkg-specific way: package precompilation fails outright with "Method
#    overwriting is not permitted during Module precompilation" in this package's layout,
#    even though the identical pattern precompiles cleanly as a plain top-level script. The keyword form
#    has a different underlying arity (0 positional arguments) from the struct's 4-argument
#    positional constructor, so the two can never collide.
# =========================================================================================

"""
    Material

Elastic constants, density and crystallographic orientation for one layer of a GULIWERR
structure.

Port of `guliwerr/core/@Material/Material.m`. **Immutable struct** -- see the file-level
comment in `Material.jl` for why.

Fields (kept camelCase, to match the MATLAB properties one-for-one):

- `name`     :: material name
- `elConsts` :: [`ElastConsts`](@ref) -- elastic constants
- `rho`      :: density [g*cm^-3]
- `dirMat`   :: 3x3 orientation matrix, row vectors (right-hand: row 3 = cross(row 1, row 2)):
    - row 1 -- surface normal (plane)
    - row 2 -- propagation direction
    - row 3 -- in-plane completion

  Two frames derive from `dirMat`, and they differ (Material.m): [`CalcRotMat`](@ref)
  gives x1 = propagation, x2 = plane, x3 = completion; [`GetStandardDirMat`](@ref) gives
  x1 = propagation, x2 = completion, x3 = plane -- the literature convention the solvers use.

Construct via the keyword-optional form [`Material(; elConsts, rho, dirMat,
name)`](@ref) (every argument optional, mirroring `Material.m`'s `nargin`-guarded
constructor -- see this file's TRANSLATION DECISIONS note 3 for why keyword rather than
positional); do not call the raw 4-field positional form (`name` first) directly outside
the package -- it performs none of the constructor's coercion (`SetElastConsts`'s
object/vector/tuple dispatch, `SetDirMat`'s right-hand completion).
"""
struct Material
    name::String
    elConsts::ElastConsts
    rho::Float64
    dirMat::Matrix{Float64}
end

"""
    Material(; elConsts=nothing, rho=nothing, dirMat=nothing, name=nothing) -> Material

Construct a `Material`, starting from MATLAB's own default (`name = "mater"`,
`elConsts = ElastConsts([100, 30], 1)`, `rho = 1`,
`dirMat = [1 0 0; 0 0 1; 0 -1 0]`, Material.m) and applying each argument that is not
`nothing`, in the same order as `Material.m`'s constructor body (elConsts, rho, dirMat,
name -- Material.m).

INPUT (all optional, all keyword -- see this file's TRANSLATION DECISIONS note 3):
`elConsts` -- an [`ElastConsts`](@ref) object, a numeric vector or matrix, or a `Tuple`
`(consts,)` / `(consts, symmNo)` (see [`SetElastConsts`](@ref) and TRANSLATION DECISIONS
note 1 for the MATLAB-cell-to-Tuple mapping); `rho` -- density [g*cm^-3]; `dirMat` -- 2x3,
3x2 or 3x3 matrix of direction row vectors (see [`SetDirMat`](@ref)/[`CompleteDirMat`](@ref));
`name` -- material name.

Port of `guliwerr/core/@Material/Material.m`.
"""
function Material(; elConsts = nothing, rho = nothing, dirMat = nothing, name = nothing)
    obj = Material("mater", ElastConsts([100.0, 30.0], 1), 1.0,
                    [1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 -1.0 0.0])
    if elConsts !== nothing
        obj = SetElastConsts(obj, elConsts)
    end
    if rho !== nothing
        obj = SetDensity(obj, rho)
    end
    if dirMat !== nothing
        obj = SetDirMat(obj, dirMat)
    end
    if name !== nothing
        obj = SetName(obj, name)
    end
    return obj
end

# -----------------------------------------------------------------------------------------
# Set methods (Material.m) -- each returns a NEW Material, per the immutability
# finding above; none mutates `obj`.
# -----------------------------------------------------------------------------------------

"""
    SetElastConsts(obj::Material, elConsts) -> Material

Set the elastic constants of a material, returning a new [`Material`](@ref).

`elConsts` may be an [`ElastConsts`](@ref) object, a numeric vector/matrix (passed straight
to the [`ElastConsts`](@ref) constructor, which infers the symmetry class), or a `Tuple`
`(consts,)` / `(consts, symmNo)` -- the Julia equivalent of MATLAB's cell-array form
`{consts}` / `{consts, symmNo}` (this file's TRANSLATION DECISIONS header note).

Port of `guliwerr/core/@Material/Material.m`.
"""
function SetElastConsts(obj::Material, elConsts)
    ec = if elConsts isa ElastConsts
        elConsts
    elseif elConsts isa Tuple
        length(elConsts) > 1 ? ElastConsts(elConsts[1], elConsts[2]) : ElastConsts(elConsts[1])
    elseif elConsts isa AbstractArray{<:Real}
        ElastConsts(elConsts)
    else
        throw(ArgumentError("Material:SetElastConsts:badInput -- elConsts must be an " *
            "ElastConsts object, numeric vector, or Tuple (consts[, symmNo])"))
    end
    return Material(obj.name, ec, obj.rho, obj.dirMat)
end

"""
    SetDirMat(obj::Material, matr) -> Material

Set the orientation matrix; the result is always right-hand orthonormal
(via [`CompleteDirMat`](@ref)). `matr` is a 2x3, 3x2 or 3x3 matrix of direction row vectors.

Port of `guliwerr/core/@Material/Material.m`.
"""
function SetDirMat(obj::Material, matr)
    return Material(obj.name, obj.elConsts, obj.rho, CompleteDirMat(matr))
end

"""
    SetDensity(obj::Material, rho) -> Material

Set the density of the material [g*cm^-3] (positive finite scalar, else `ArgumentError`,
as MATLAB's `validateattributes`).

Port of `guliwerr/core/@Material/Material.m`.
"""
function SetDensity(obj::Material, rho::Real)
    (isfinite(rho) && rho > 0) ||
        throw(ArgumentError("Material.SetDensity: rho must be a positive finite scalar, got $rho"))
    return Material(obj.name, obj.elConsts, Float64(rho), obj.dirMat)
end

"""
    SetName(obj::Material, name) -> Material

Set the material name.

Port of `guliwerr/core/@Material/Material.m`.
"""
function SetName(obj::Material, name::AbstractString)
    return Material(String(name), obj.elConsts, obj.rho, obj.dirMat)
end

# -----------------------------------------------------------------------------------------
# Rotation methods (Material.m)
# -----------------------------------------------------------------------------------------

"""
    GetRotMat(obj::Material, angleDeg, relativeToggle=false) -> Matrix{Float64}

Rotation matrix (3x2) for a given in-plane rotation angle.

`relativeToggle`: `false`/`0` (default) -- rotate within `obj.dirMat`; `true`/`1` -- rotate
in the absolute `[100]`/`[010]` frame. Anything else throws
(`Material:GetRotMat:badInput`).

Port of `Material.GetRotMat` (logical toggle; 0/1 also accepted).
"""
function GetRotMat(obj::Material, angleDeg::Real, relativeToggle::Union{Nothing,Bool,Integer} = nothing)
    rt = relativeToggle === nothing ? false : relativeToggle
    (rt isa Bool || rt == 0 || rt == 1) ||
        throw(ArgumentError("Material:GetRotMat:badInput -- relativeToggle must be false " *
            "(relative to dirMat) or true (absolute frame)"))
    dMat = Bool(rt) ? CompleteDirMat([1.0 0.0 0.0; 0.0 1.0 0.0]) : obj.dirMat
    return CalcRotMat(dMat, angleDeg)
end

"""
    GetRotElastTensor(obj::Material, angleDeg, relativeToggle=false) -> Array{Float64,4}

Elastic tensor (3x3x3x3) rotated by the given in-plane angle. The angle is required
(`0` for the dirMat frame; `obj.elConsts.tensor` is the unrotated tensor): without it this
throws `Material:GetRotElastTensor:noAngle` (as in MATLAB).

Port of `Material.GetRotElastTensor`.
"""
function GetRotElastTensor(obj::Material, angleDeg::Union{Nothing,Real} = nothing,
                            relativeToggle::Union{Nothing,Bool,Integer} = nothing)
    angleDeg === nothing &&
        throw(ArgumentError("Material:GetRotElastTensor:noAngle -- An angle is required (0 for " *
            "the dirMat frame; obj.elConsts.tensor is the unrotated tensor)."))
    rotMat = GetRotMat(obj, angleDeg, relativeToggle)
    return RotateTensor(obj.elConsts.tensor, rotMat)
end

"""
    GetRotElastVoigt(obj::Material, angleDeg, relativeToggle=false) -> Matrix{Float64}

Voigt matrix (6x6) of the elastic tensor rotated by the given in-plane angle. The angle is
required (`Material:GetRotElastVoigt:noAngle`, as [`GetRotElastTensor`](@ref)).

Port of `Material.GetRotElastVoigt`.
"""
function GetRotElastVoigt(obj::Material, angleDeg::Union{Nothing,Real} = nothing,
                           relativeToggle::Union{Nothing,Bool,Integer} = nothing)
    angleDeg === nothing &&
        throw(ArgumentError("Material:GetRotElastVoigt:noAngle -- An angle is required (0 for " *
            "the dirMat frame; obj.elConsts.voigt is the unrotated matrix)."))
    return Tensor2Voigt(GetRotElastTensor(obj, angleDeg, relativeToggle))
end

"""
    GetStandardElastTensor(obj::Material, angleDeg=0.0) -> Array{Float64,4}

Elastic tensor in the standard calculation frame (x1 propagation, x3 surface normal) at
in-plane angle `angleDeg`: `GetStandardDirMat(obj.dirMat, angleDeg)` then
`GetRotElastTensor(..., 0)`. The one place that chain lives; the
limiting velocities, the Green solver, the spot maps and the Ritz score all start from it.

Port of `Material.GetStandardElastTensor`.
"""
function GetStandardElastTensor(obj::Material, angleDeg::Real = 0.0)
    matStd = SetDirMat(obj, GetStandardDirMat(obj.dirMat, angleDeg))
    return GetRotElastTensor(matStd, 0.0)
end

"""
    GetRotDirMat(obj::Material, angleDeg) -> Matrix{Float64}

Direction matrix rotated in-plane by the given angle (row 1, the surface normal, is the
pivot and stays unchanged).

Port of `guliwerr/core/@Material/Material.m`.
"""
function GetRotDirMat(obj::Material, angleDeg::Real)
    return RotDirMatFromAngle(obj.dirMat, angleDeg)
end

export Material, SetElastConsts, SetDirMat, SetDensity, SetName,
       GetRotMat, GetRotElastTensor, GetRotElastVoigt, GetRotDirMat, GetStandardElastTensor
