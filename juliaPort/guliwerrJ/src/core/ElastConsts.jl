# =========================================================================================
# ElastConsts.jl
#
# Port of guliwerr/core/@ElastConsts/ElastConsts.m (the type + constructor) plus its
# Voigt/vector/tensor conversion methods: Voigt2Vector.m, Vector2Voigt.m, Any2Tensor.m,
# Tensor2Voigt.m, Tensor2Vector.m, AsColumnVec.m.
#
# Symmetry-class bookkeeping (SymmClassLibrary, FindType, FindSymmClass_input,
# IndepVoigtPositions, GetIndOf2DReduct_15of21, the six
# {Cubic,Hexagonal,Tetragonal}{,De}conversion functions) lives in ElastConstsSymmetry.jl;
# rotation and derived-modulus functions (RotateTensor, ElastRotReductMat, RightHandFrame,
# Tensor2Youngs, Voigt2Bulk, HillAverages, VoigtAverages, ReussAverages) live in
# ElastConstsDerived.jl. `CalcAnisotropyRatio.m` is EXCLUDED from all three files -- it
# lives in `solvers/christoffel/AnisotropyRatio.jl`, because it calls
# `Christoffel.phaseVels` and `Christoffel` in turn depends on `ElastConsts`
# (SymmClassLibrary, Any2Tensor, Tensor2Vector, the constructor) -- keeping it here would
# make a dependency cycle between core and solvers. `ElastConsts.GetAnisotropyRatio` (the
# instance-method wrapper that dispatches to CalcAnisotropyRatio for symmClassNo >= 3) is
# left out for the same reason.
#
# All three files are `include`d directly into `module CoreTypes` (see src/core/CoreTypes.jl)
# -- there is no `module ElastConsts` wrapper (that would collide with the `ElastConsts`
# *type* name defined below, and a MATLAB value class is not a namespace, so it gets no
# module of its own; only genuine static-namespace classdefs like `ParallelPolicy` become
# Julia modules). Every name below is therefore a plain type/function living directly in
# the `CoreTypes` module -- there is no "ElastConsts." prefix as in MATLAB's
# `ElastConsts.Any2Tensor(...)` static-method call syntax.
#
# NAMING NOTE: every function name below (`Voigt2Vector`, `Any2Tensor`, `Tensor2Voigt`,
# `Tensor2Vector`, `AsColumnVec`) is kept IDENTICAL to its MATLAB original, rather than
# renamed to lower_snake_case, so the port stays greppable against the MATLAB source.
#
# MUTABLE vs IMMUTABLE: `classdef ElastConsts` has no `< handle` -- it is a MATLAB
# *value* class. Its methods that look like mutators (`SetNewSymmClass`) follow the
# `obj = obj.Method(...)`
# rebuild-and-return idiom; they never mutate `obj` in place and hand back the same
# object. Nothing in the MATLAB code assigns directly to the properties of
# an ElastConsts instance outside its own constructor (the only assignments to identically
# named properties, `obj.nIndepConsts` and `obj.symmClassNo` in
# guliwerr/inverseCalcs/SensitivityAnalysis.m, belong to that class's own properties).
# Conclusion: **immutable struct**, the default mapping for a value class with no in-place
# mutation.
#
# ALIASING NOTE -- a genuine MATLAB-to-Julia semantics gap, not a MATLAB bug: MATLAB
# arrays are value types (assignment copies); Julia arrays are reference types
# (assignment aliases the same backing memory). `Any2Tensor`'s tensor-input branch is a
# straight passthrough in MATLAB (`outTen = inp;`), which is a real, independent copy
# there. The Julia port copies explicitly (`Array{Float64}(inp)`) at that one spot so that
# a caller's later in-place mutation of their own tensor array can never retroactively
# change a tensor already stored inside a (supposedly immutable) `ElastConsts` value.
# =========================================================================================

"""
    ElastConsts

Elastic constants in the three formats GULIWERR uses: independent-constant vector, Voigt
6x6 matrix, and full 3x3x3x3 tensor, kept in sync with a symmetry class.

Port of `guliwerr/core/@ElastConsts/ElastConsts.m`. **Immutable struct** -- see the
file-level comment in `ElastConsts.jl` for why.

Fields mirror the MATLAB properties one-for-one (kept camelCase, to ease
side-by-side comparison with the MATLAB source):

- `vector`       :: independent-constant vector, length `nIndepConsts`
- `voigt`        :: 6x6 Voigt matrix
- `tensor`       :: 3x3x3x3 elastic tensor
- `symmClassNo`  :: symmetry class number, 1-10 (see [`SymmClassLibrary`](@ref))
- `nIndepConsts` :: number of independent constants (`== length(vector)`)
- `shearAnisotropyRatio` :: `(v_shear,max / v_shear,min)^2` over ALL directions, `>= 1`,
  rotation-invariant; stored only where a closed form exists (isotropic `1`, cubic
  `max(A_Z, 1/A_Z)` with the Zener ratio `A_Z = 2 c44 / (c11 - c12)`), `NaN` otherwise
  (no closed form; `calc_anisotropy_ratio` scans). A cache, read by `GreenSolver.calc`'s
  near-isotropy check.

Construct via [`ElastConsts(consts, symmNo)`](@ref) or the zero-argument default; do not
call this positional-6-field form directly outside the package (it performs none of the
constructor's validation or derivation).
"""
struct ElastConsts
    vector::Vector{Float64}
    voigt::Matrix{Float64}
    tensor::Array{Float64,4}
    symmClassNo::Int
    nIndepConsts::Int
    shearAnisotropyRatio::Float64
end

"""
    ZenerRatio(cubicVec) -> Float64

Zener anisotropy ratio `A = 2 c44 / (c11 - c12)` of a cubic `[c11, c12, c44]` vector (NOT
the shear anisotropy ratio, which is `max(A, 1/A)`). Port of `ElastConsts.ZenerRatio`.
"""
ZenerRatio(cubicVec::AbstractVector{<:Real}) = 2 * cubicVec[3] / (cubicVec[1] - cubicVec[2])

"""Closed-form `shearAnisotropyRatio` (see [`ElastConsts`](@ref)); `NaN` without one."""
function _shear_anisotropy_ratio(vector::AbstractVector{<:Real}, symmNo::Integer)
    symmNo == 1 && return 1.0
    if symmNo == 2
        zener = ZenerRatio(vector)
        (isfinite(zener) && zener > 0) && return max(zener, 1 / zener)
    end
    return NaN
end

"""
    ElastConsts(consts, symmNo=nothing) -> ElastConsts
    ElastConsts() -> ElastConsts

Construct from a vector of independent constants, a 6x6 Voigt matrix, or a 3x3x3x3
tensor. `symmNo` (see [`SymmClassLibrary`](@ref)) is inferred via
[`FindSymmClass_input`](@ref) when omitted -- see that function's own docstring for when
inference is unreliable (prefer passing `symmNo` explicitly for Voigt/tensor input, and
always for a 6- or 7-constant vector, which is ambiguous between tetragonal and trigonal).

The zero-argument form gives MATLAB's default: a blank triclinic (`symmNo = 10`) object of
NaNs, matching `ElastConsts()` in `ElastConsts.m`.

Port of `guliwerr/core/@ElastConsts/ElastConsts.m` (the constructor body only).
`SetNewSymmClass` and every `Get*` instance-method wrapper
(`GetReducedElastConsts`, `GetYoungsModulus`, `GetVoigtAverages`,
`GetReussAverages`, `GetHillAverages`, `GetAnisotropyRatio`, `GetRotElastTensor`,
`GetRotElastVoigt`, `GetCubicConversion`, `GetTetragonalConversion`,
`GetHexagonalConversion`) are **not ported** -- in Julia these wrappers add nothing over calling the underlying free
function directly on the struct's own fields, e.g.
`ElastRotReductMat(rotMat, ec.symmClassNo) * ec.vector`
in place of `ec.GetReducedElastConsts(rotMat)`. `GetAnisotropyRatio` additionally cannot
be ported here at all, cubic/isotropic branch included, without either duplicating or
forward-referencing `CalcAnisotropyRatio` -- see this file's header comment; the
closed-form part lives in the `shearAnisotropyRatio` field. (MATLAB no longer has
`SetElastVector`, `GetBulkModulus` and `GetShearAnisotropyRatio`;
their replacements are `ElastConsts(vector, ec.symmClassNo)`, `Voigt2Bulk(ec.voigt)` and
`ec.shearAnisotropyRatio`.)
"""
function ElastConsts(consts, symmNo::Union{Nothing,Integer} = nothing)
    sN = symmNo === nothing ? FindSymmClass_input(consts) : Int(symmNo)

    _, expectedNumberOfConsts = SymmClassLibrary(sN)

    # For vector input, validate the count before Vector2Voigt silently under-reads it
    # (ElastConsts.m: "validate count before Vector2Voigt silently truncates").
    # NOTE: Julia's `ndims` differs from MATLAB's for a 1-D `Vector` (ndims == 1, not 2);
    # `FindType` (ElastConstsSymmetry.jl) already accounts for that language difference,
    # so reusing it here keeps this guard equivalent to MATLAB's `isvector(consts)`.
    if FindType(consts) == "vector" && length(consts) > 1 &&
       length(consts) != expectedNumberOfConsts
        throw(ArgumentError("ElastConsts:constructor:wrongCount -- number of input " *
                             "constants does not match the symmetry class."))
    end

    tensor = Any2Tensor(consts, sN)
    voigt  = Tensor2Voigt(tensor)
    vector = Tensor2Vector(tensor, sN)

    nIndepConsts = length(vector)
    if nIndepConsts != expectedNumberOfConsts
        throw(ArgumentError("ElastConsts:constructor:wrongCount -- number of input " *
                             "constants does not match the symmetry class."))
    end

    vecF = Float64.(collect(vector))
    return ElastConsts(vecF, Float64.(collect(voigt)), Float64.(collect(tensor)), sN,
                        Int(nIndepConsts), _shear_anisotropy_ratio(vecF, sN))
end

ElastConsts() = ElastConsts(fill(NaN, 6, 6), 10)

# -----------------------------------------------------------------------------------------
# Voigt2Vector.m
# -----------------------------------------------------------------------------------------

"""
    Voigt2Vector(cij, symmNo=0) -> Vector{Float64}

Transform a Voigt matrix to the vector of independent constants, in
[`IndepVoigtPositions`](@ref) order. `symmNo = 0` (the default, matching MATLAB's
`nargin == 1 || symmNo == 0`) infers the class via [`FindSymmClass_input`](@ref).

Port of `guliwerr/core/@ElastConsts/Voigt2Vector.m`.
"""
function Voigt2Vector(cij::AbstractMatrix{<:Real}, symmNo::Integer = 0)
    sN = symmNo == 0 ? FindSymmClass_input(cij) : symmNo
    if !(1 <= sN <= 10)
        throw(ArgumentError("ElastConsts:Voigt2Vector:unknownSymm -- unknown symmetry " *
                             "class number: $sN (must be 1-10)"))
    end
    posList = IndepVoigtPositions(sN)
    return [Float64(cij[r, c]) for (r, c) in posList]
end

# -----------------------------------------------------------------------------------------
# Vector2Voigt.m
# -----------------------------------------------------------------------------------------

"""
    Vector2Voigt(inpVec, symmNo=0) -> Matrix{Float64}

Transform a vector of independent elastic constants into a 6x6 Voigt matrix. The order in
which `inpVec` is read per symmetry class mirrors [`IndepVoigtPositions`](@ref) exactly
(the vector-to-Voigt-to-vector round trip depends on the two staying in step).

Port of `guliwerr/core/@ElastConsts/Vector2Voigt.m`.
"""
function Vector2Voigt(inpVec, symmNo::Integer = 0)
    sN = symmNo == 0 ? FindSymmClass_input(inpVec) : symmNo
    v = vec(inpVec)   # accept a Vector, or a 1xN / Nx1 Matrix, exactly like MATLAB's inpVec(k)

    if sN == 1                                   # isotropic
        c11, c44 = v[1], v[2]
        c12 = c11 - 2c44
        return Float64[c11 c12 c12 0 0 0
                        c12 c11 c12 0 0 0
                        c12 c12 c11 0 0 0
                        0   0   0  c44 0  0
                        0   0   0   0 c44 0
                        0   0   0   0  0 c44]
    elseif sN == 2                                # cubic
        c11, c12, c44 = v[1], v[2], v[3]
        return Float64[c11 c12 c12 0 0 0
                        c12 c11 c12 0 0 0
                        c12 c12 c11 0 0 0
                        0   0   0  c44 0  0
                        0   0   0   0 c44 0
                        0   0   0   0  0 c44]
    elseif sN == 3                                # tetragonal, type I
        c11, c12, c13, c33, c44, c66 = v[1], v[2], v[3], v[4], v[5], v[6]
        return Float64[c11 c12 c13 0 0 0
                        c12 c11 c13 0 0 0
                        c13 c13 c33 0 0 0
                        0   0   0  c44 0  0
                        0   0   0   0 c44 0
                        0   0   0   0  0 c66]
    elseif sN == 4                                # tetragonal, type II (c16 != 0)
        c11, c12, c13, c16, c33, c44, c66 = v[1], v[2], v[3], v[4], v[5], v[6], v[7]
        return Float64[c11  c12  c13 0   0   c16
                        c12  c11  c13 0   0  -c16
                        c13  c13  c33 0   0    0
                        0    0    0  c44  0    0
                        0    0    0   0  c44   0
                        c16 -c16   0   0   0  c66]
    elseif sN == 5                                # trigonal, type I
        c11, c12, c13, c14, c33, c44 = v[1], v[2], v[3], v[4], v[5], v[6]
        c66 = (c11 - c12) / 2
        return Float64[c11  c12  c13  c14  0   0
                        c12  c11  c13 -c14  0   0
                        c13  c13  c33  0    0   0
                        c14 -c14   0  c44   0   0
                        0    0     0   0  c44  c14
                        0    0     0   0  c14  c66]
    elseif sN == 6                                # trigonal, type II (c15 != 0)
        c11, c12, c13, c14, c15, c33, c44 = v[1], v[2], v[3], v[4], v[5], v[6], v[7]
        c66 = (c11 - c12) / 2
        return Float64[c11  c12  c13  c14  c15  0
                        c12  c11  c13 -c14 -c15  0
                        c13  c13  c33  0    0    0
                        c14 -c14   0  c44   0  -c15
                        c15 -c15   0   0   c44  c14
                        0    0     0 -c15  c14  c66]
    elseif sN == 7                                # hexagonal
        c11, c12, c13, c33, c44 = v[1], v[2], v[3], v[4], v[5]
        c66 = (c11 - c12) / 2
        return Float64[c11 c12 c13 0 0 0
                        c12 c11 c13 0 0 0
                        c13 c13 c33 0 0 0
                        0   0   0  c44 0  0
                        0   0   0   0 c44 0
                        0   0   0   0  0 c66]
    elseif sN == 8                                # orthorhombic
        c11, c12, c13, c22, c23, c33, c44, c55, c66 =
            v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9]
        return Float64[c11 c12 c13 0 0 0
                        c12 c22 c23 0 0 0
                        c13 c23 c33 0 0 0
                        0   0   0  c44 0  0
                        0   0   0   0 c55 0
                        0   0   0   0  0 c66]
    elseif sN == 9                                # monoclinic
        c11, c12, c13, c15, c22, c23, c25, c33, c35, c44, c46, c55, c66 =
            v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13]
        return Float64[c11 c12 c13  0  c15  0
                        c12 c22 c23  0  c25  0
                        c13 c23 c33  0  c35  0
                        0   0   0  c44  0  c46
                        c15 c25 c35  0  c55  0
                        0   0   0  c46  0  c66]
    elseif sN == 10                               # triclinic
        c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66 =
            v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],
            v[15],v[16],v[17],v[18],v[19],v[20],v[21]
        return Float64[c11 c12 c13 c14 c15 c16
                        c12 c22 c23 c24 c25 c26
                        c13 c23 c33 c34 c35 c36
                        c14 c24 c34 c44 c45 c46
                        c15 c25 c35 c45 c55 c56
                        c16 c26 c36 c46 c56 c66]
    else
        throw(ArgumentError("ElastConsts:Vector2Voigt:unknownSymm -- unknown symmetry " *
                             "class number: $sN (must be 1-10)"))
    end
end

# -----------------------------------------------------------------------------------------
# Any2Tensor.m
# -----------------------------------------------------------------------------------------

"""
    Any2Tensor(inp, symmNo=0) -> Array{Float64,4}

Transform any form of elastic constants (vector / Voigt 6x6 / tensor 3x3x3x3) to the full
3x3x3x3 tensor.

Port of `guliwerr/core/@ElastConsts/Any2Tensor.m`.

**Not `@inferred`-stable** (Julia cannot infer a single concrete return type): the
"vector" branch recurses into this same function with a Voigt `Matrix{Float64}`
argument, and whether THAT call takes the "voigt" branch (6x6 shape) or the "vector"
branch (a degenerate 1xN/Nx1 shape) depends on `size(...)`, a runtime value the compiler
cannot see through the recursion -- so Julia infers `Union{Array{Float64,4},
Matrix{Float64}}` for that path even though every real call in this package always
bottoms out at a concrete `Array{Float64,4}`. This mirrors MATLAB's own runtime
`FindType`-based dispatch faithfully; it is not a hot-loop entry point (unlike
`RotateTensor`/`Tensor2Youngs`, which are `@inferred`-stable), so the Union return is left as-is
rather than restructured purely for inference's sake.
"""
function Any2Tensor(inp, symmNo::Integer = 0)
    t = FindType(inp)
    if t == "tensor"
        # Explicit copy, not a passthrough alias -- see this file's ALIASING NOTE header.
        return Array{Float64}(inp)
    elseif t == "voigt"
        # Voigt index pair -> tensor index pair (Julia matrix-literal syntax already
        # matches MATLAB's `[1 6 5; 6 2 4; 5 4 3]` element-for-element -- no row/col-major
        # translation needed since this table happens to be symmetric).
        p = [1 6 5; 6 2 4; 5 4 3]
        outTen = zeros(Float64, 3, 3, 3, 3)
        for l in 1:3, k in 1:3, j in 1:3, i in 1:3
            outTen[i, j, k, l] = inp[p[i, j], p[k, l]]
        end
        return outTen
    elseif t == "vector"
        cij = Vector2Voigt(inp, symmNo)
        return Any2Tensor(cij, symmNo)
    else
        throw(ArgumentError("ElastConsts:Any2Tensor:badInput -- input must be a vector, " *
                             "Voigt (6x6), or elastic tensor (3x3x3x3)"))
    end
end

# -----------------------------------------------------------------------------------------
# Tensor2Voigt.m
# -----------------------------------------------------------------------------------------

"""
    Tensor2Voigt(inpTen) -> Matrix{Float64}

Transform the elastic tensor (3x3x3x3) to a Voigt matrix (6x6).

Port of `guliwerr/core/@ElastConsts/Tensor2Voigt.m`.
"""
function Tensor2Voigt(inpTen::AbstractArray{<:Real,4})
    T = inpTen
    return Float64[
        T[1,1,1,1] T[1,1,2,2] T[1,1,3,3] T[1,1,2,3] T[1,1,1,3] T[1,1,1,2]
        T[1,1,2,2] T[2,2,2,2] T[2,2,3,3] T[2,2,2,3] T[2,2,1,3] T[2,2,1,2]
        T[1,1,3,3] T[2,2,3,3] T[3,3,3,3] T[3,3,2,3] T[3,3,1,3] T[3,3,1,2]
        T[1,1,2,3] T[2,2,2,3] T[3,3,2,3] T[2,3,2,3] T[2,3,1,3] T[2,3,1,2]
        T[1,1,1,3] T[2,2,1,3] T[3,3,1,3] T[2,3,1,3] T[1,3,1,3] T[1,3,1,2]
        T[1,1,1,2] T[2,2,1,2] T[3,3,1,2] T[2,3,1,2] T[1,3,1,2] T[1,2,1,2]
    ]
end

# -----------------------------------------------------------------------------------------
# Tensor2Vector.m
# -----------------------------------------------------------------------------------------

"""
    Tensor2Vector(inpTen, symmNo=0) -> Vector{Float64}

Transform the elastic tensor to the vector of independent constants. `symmNo = 0` (the
default) infers the class via [`FindSymmClass_input`](@ref).

Port of `guliwerr/core/@ElastConsts/Tensor2Vector.m`.
"""
function Tensor2Vector(inpTen::AbstractArray{<:Real,4}, symmNo::Integer = 0)
    sN = symmNo == 0 ? FindSymmClass_input(inpTen) : symmNo
    c = Tensor2Voigt(inpTen)
    return Voigt2Vector(c, sN)
end

# -----------------------------------------------------------------------------------------
# AsColumnVec.m -- MATLAB `methods (Static, Access = private)`; kept unexported here.
# -----------------------------------------------------------------------------------------

"""
    AsColumnVec(vecIn, nExpected, fname) -> Vector{Float64}

Normalise an optional elastic-constant vector input to an `nExpected`-element column.
Empty input (or, in this port, `nothing`) gives an all-zeros vector -- how the
`*Conversion` functions ([`CubicConversion`](@ref) and friends) are called when only
their transformation matrix (their second return value) is wanted.

Port of `guliwerr/core/@ElastConsts/AsColumnVec.m`. Private in MATLAB
(`methods (Static, Access = private)`); Julia has no access control, so
this is simply left unexported rather than renamed.
"""
function AsColumnVec(vecIn, nExpected::Integer, fname::AbstractString)
    if vecIn === nothing || isempty(vecIn)
        return zeros(Float64, nExpected)
    end
    ok = (vecIn isa AbstractVector) ||
         (vecIn isa AbstractMatrix && any(==(1), size(vecIn)))
    if !ok || length(vecIn) != nExpected
        throw(ArgumentError("ElastConsts:AsColumnVec:wrongSize -- ElastConsts.$fname: " *
                             "wrong input size -- expected a $nExpected-element vector."))
    end
    return Float64.(vec(vecIn))
end

export ElastConsts, Voigt2Vector, Vector2Voigt, Any2Tensor, Tensor2Voigt, Tensor2Vector, ZenerRatio
