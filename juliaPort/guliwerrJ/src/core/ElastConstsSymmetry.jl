# =========================================================================================
# ElastConstsSymmetry.jl
#
# Port of guliwerr/core/@ElastConsts/SymmClassLibrary.m, FindType.m, FindSymmClass_input.m,
# IndepVoigtPositions.m, GetIndOf2DReduct_15of21.m, and the six
# {Cubic,Hexagonal,Tetragonal}{,De}conversion.m files. `include`d directly into
# `module CoreTypes` alongside ElastConsts.jl/ElastConstsDerived.jl -- see the header
# comment at the top of ElastConsts.jl for the module-layout and naming-convention
# rationale shared by all three files (no `module ElastConstsSymmetry` wrapper; MATLAB
# names kept verbatim; immutable-struct decision recorded there, not repeated here).
#
# `using LinearAlgebra` here (needed for `inv` and the matrix left-division `\` inside the
# three *Deconversion functions) -- Project.toml already lists LinearAlgebra as a `[deps]`
# entry.
# =========================================================================================

using LinearAlgebra

# -----------------------------------------------------------------------------------------
# FindType.m -- MATLAB `methods (Static, Access = private)`; kept unexported here.
# -----------------------------------------------------------------------------------------

"""
    FindType(input) -> String

Identify the format of elastic-constant input from its dimensions: `"tensor"` (3x3x3x3),
`"voigt"` (6x6), `"vector"` (any vector shape, including a bare scalar), or `"unknown"`.

Port of `guliwerr/core/@ElastConsts/FindType.m`. Private in MATLAB; left unexported here
(Julia has no access control, so this is a naming convention, not enforcement).

**Language-difference note**: MATLAB's `ndims` is always >= 2 (a MATLAB "vector" is a
1xN or Nx1 array with `ndims == 2`; a MATLAB scalar is a 1x1 array, `ndims == 2` too).
Julia's `ndims` is 1 for a plain `Vector` and 0 for a bare `Number` (not an array at
all). This is handled here by dispatching separately on `Number` (-> `"vector"`, matching
MATLAB's scalar-is-a-1x1-array behaviour) and by treating Julia's `ndims == 1` the same as
MATLAB's "one of the two dims is 1" vector case, rather than by trying to force Julia
arrays to imitate MATLAB's always->=2-dimensional convention.
"""
FindType(input::Number) = "vector"

function FindType(input::AbstractArray)
    nd = ndims(input)
    if nd == 4
        return "tensor"
    elseif nd == 1
        return "vector"
    elseif nd == 2
        sz = size(input)
        if any(==(1), sz)
            return "vector"
        elseif sz == (6, 6)
            return "voigt"
        else
            return "unknown"
        end
    else
        return "unknown"
    end
end

FindType(input) = "unknown"   # anything else (matches MATLAB's implicit "otherwise")

# -----------------------------------------------------------------------------------------
# FindSymmClass_input.m
# -----------------------------------------------------------------------------------------

"""
    FindSymmClass_input(inp) -> Int

Try to identify the symmetry class number from the input's structure.

Port of `guliwerr/core/@ElastConsts/FindSymmClass_input.m`.

**GOTCHA (carried verbatim from the MATLAB docstring)**: reliable only from a vector of
independent constants. Voigt (and tensor, converted to Voigt first) input counts
`unique(cij)` entries, which miscounts whenever two independent constants happen to be
numerically equal. Prefer supplying `symmNo` explicitly.

**Language-difference note**: MATLAB's `unique` treats every `NaN` as distinct from every
other `NaN` (so `unique([NaN,NaN])` keeps both); Julia's `unique` uses `isequal`, under
which `NaN` equals itself, so `unique([NaN,NaN])` collapses to one entry. This only
matters if this function is ever called directly on NaN-filled input (it never is from
the `ElastConsts` constructor, which always resolves `symmNo` before touching the
constructor's default-NaN path); not worth reproducing MATLAB's distinct-NaN quirk here.
"""
function FindSymmClass_input(inp)
    isVoigt = true
    local NUnique::Int
    local cij

    t = FindType(inp)
    if t == "tensor"
        cij = Tensor2Voigt(inp)
    elseif t == "voigt"
        cij = inp
    elseif t == "vector"
        isVoigt = false
        NUnique = length(inp)
    else
        throw(ArgumentError("ElastConsts:FindSymmClass_input:badInput -- input must be " *
                             "a vector, Voigt (6x6), or elastic tensor (3x3x3x3)"))
    end

    if isVoigt
        CUnique = sort(unique(vec(cij)))
        if !isempty(CUnique) && CUnique[1] == 0
            CUnique = CUnique[2:end]   # the structural zeros of every class but triclinic
        end
        NUnique = length(CUnique)
    end

    local symmNo::Int
    if NUnique == 2
        symmNo = 1
        if isVoigt && cij[1,2] != (cij[1,1] - 2*cij[4,4])
            throw(ArgumentError("ElastConsts:FindSymmClass_input:notSymmetric -- symmetry check " *
                                 "failed for assumed class $symmNo (isotropic)"))
        end
    elseif NUnique == 3
        symmNo = 2
    elseif NUnique == 5
        symmNo = 7
    elseif NUnique == 6
        if isVoigt
            if cij[1,4] == 0
                symmNo = 3          # tetragonal
            else
                symmNo = 5          # trigonal
                if cij[6,6] != (cij[1,1] - cij[1,2]) / 2
                    throw(ArgumentError("ElastConsts:FindSymmClass_input:notSymmetric -- symmetry " *
                                         "check failed for assumed class $symmNo (trigonal I)"))
                end
            end
        else
            throw(ArgumentError("ElastConsts:FindSymmClass_input:ambiguousClass -- symmetry class not " *
                "given and cannot be determined (6 constants: tetragonal(3) or trigonal(5))"))
        end
    elseif NUnique == 7
        if isVoigt && cij[1,4] == 0
            symmNo = 4
        elseif isVoigt
            symmNo = 6
            if cij[6,6] != (cij[1,1] - cij[1,2]) / 2
                throw(ArgumentError("ElastConsts:FindSymmClass_input:notSymmetric -- symmetry check " *
                                     "failed for assumed class $symmNo (trigonal II)"))
            end
        else
            throw(ArgumentError("ElastConsts:FindSymmClass_input:ambiguousClass -- symmetry class not " *
                "given and cannot be determined (7 constants: tetragonal(4) or trigonal(6))"))
        end
    elseif NUnique == 9
        symmNo = 8
    elseif NUnique == 13
        symmNo = 9
    elseif NUnique == 21
        symmNo = 10
    else
        throw(ArgumentError("ElastConsts:FindSymmClass_input:ambiguous -- cannot " *
                             "determine symmetry class: $NUnique unique constants found"))
    end

    return symmNo
end

# -----------------------------------------------------------------------------------------
# SymmClassLibrary.m
# -----------------------------------------------------------------------------------------

const _SYMM_CLASS_NAMES = String[
    "isotropic: [c11,c44]",
    "cubic: [c11,c12,c44]",
    "tetragonal, type I: [c11,c12,c13,c33,c44,c66]",
    "tetragonal, type II  (c_16 ~= 0): [c11,c12c,c13,c16,c33,c44,c66]",
    "trigonal, type I: [c11,c12,c13,c14,c33,c44]",
    "trigonal, type II  (c_15 ~= 0): [c11,c12,c13,c14,c15,c33,c44]",
    "hexagonal: [c11,c12,c13,c33,c44]",
    "orthogonal: [c11,c12,c13,c22,c23,c33,c44,c55,c66]",
    "monoclinic: [c11,c12,c13,c15,c22,c23,c25,c33,c35,c44,c46,c55,c66]",
    "triclinic: [c11,c12,...,c16,c22,c23,...,c26,c33,...]",
]
const _SYMM_CLASS_COUNTS = Int[2, 3, 6, 7, 6, 7, 5, 9, 13, 21]

"""
    SymmClassLibrary() -> (names::Vector{String}, counts::Vector{Int})
    SymmClassLibrary(symmNo) -> (name::String, count)

Library of symmetry classes and their number of independent constants. With no argument,
prints (via `println`, in place of MATLAB's `disp`) and returns all ten. With `symmNo`,
returns that class's name and count; anything outside 1-10 throws
`ElastConsts:SymmClassLibrary:unknownSymm` (as in MATLAB).

Port of `guliwerr/core/@ElastConsts/SymmClassLibrary.m`.
"""
function SymmClassLibrary()
    for ii in 1:10
        println(lpad(ii, 2, '0'), "         ", lpad(_SYMM_CLASS_COUNTS[ii], 2, '0'),
                "               ", _SYMM_CLASS_NAMES[ii])
    end
    return _SYMM_CLASS_NAMES, _SYMM_CLASS_COUNTS
end

function SymmClassLibrary(symmNo::Real)
    (isinteger(symmNo) && 1 <= symmNo <= 10) ||
        throw(ArgumentError("ElastConsts:SymmClassLibrary:unknownSymm -- Unknown symmetry " *
                             "class number (must be 1-10)."))
    return _SYMM_CLASS_NAMES[Int(symmNo)], _SYMM_CLASS_COUNTS[Int(symmNo)]
end

# -----------------------------------------------------------------------------------------
# IndepVoigtPositions.m
# -----------------------------------------------------------------------------------------

"""
    IndepVoigtPositions(symmNo) -> Vector{Tuple{Int,Int}}

Grid positions (`(row, col)`, upper-triangle 6x6 Voigt positions) of each independent
Voigt constant for symmetry class `symmNo`, one per [`Vector2Voigt`](@ref) argument entry,
in the SAME order `Vector2Voigt` reads `inpVec[k]` for that class. Must stay in step with
`Vector2Voigt`'s branches -- the vector-to-Voigt round trip fails if the two drift apart.

Port of `guliwerr/core/@ElastConsts/IndepVoigtPositions.m` (MATLAB cell array of `[row,col]`
pairs here becomes a `Vector` of `(row,col)` tuples). MATLAB DERIVES the list from
`Vector2Voigt` (the first upper-triangle position, row-major, fed by one constant alone;
cached per class) instead of hard-coding it; the result is identical to this table for all
ten classes, so the table is kept -- it is the cheaper of the two equal answers.
"""
function IndepVoigtPositions(symmNo::Integer)
    if symmNo == 1
        return [(1,1), (4,4)]
    elseif symmNo == 2
        return [(1,1), (1,2), (4,4)]
    elseif symmNo == 3
        return [(1,1), (1,2), (1,3), (3,3), (4,4), (6,6)]
    elseif symmNo == 4
        return [(1,1), (1,2), (1,3), (1,6), (3,3), (4,4), (6,6)]
    elseif symmNo == 5
        return [(1,1), (1,2), (1,3), (1,4), (3,3), (4,4)]
    elseif symmNo == 6
        return [(1,1), (1,2), (1,3), (1,4), (1,5), (3,3), (4,4)]
    elseif symmNo == 7
        return [(1,1), (1,2), (1,3), (3,3), (4,4)]
    elseif symmNo == 8
        return [(1,1), (1,2), (1,3), (2,2), (2,3), (3,3), (4,4), (5,5), (6,6)]
    elseif symmNo == 9
        return [(1,1), (1,2), (1,3), (1,5), (2,2), (2,3), (2,5), (3,3), (3,5), (4,4),
                (4,6), (5,5), (6,6)]
    elseif symmNo == 10
        return [(1,1), (1,2), (1,3), (1,4), (1,5), (1,6),
                (2,2), (2,3), (2,4), (2,5), (2,6),
                (3,3), (3,4), (3,5), (3,6),
                (4,4), (4,5), (4,6),
                (5,5), (5,6),
                (6,6)]
    else
        throw(ArgumentError("ElastConsts:IndepVoigtPositions:unknownSymm -- unknown " *
                             "symmetry class number: $symmNo (must be 1-10)"))
    end
end

# -----------------------------------------------------------------------------------------
# GetIndOf2DReduct_15of21.m
# -----------------------------------------------------------------------------------------

"""
    GetIndOf2DReduct_15of21() -> (indSet::Vector{Int}, nReduct::Int)

Indices (into a 21-element triclinic vector, [`Vector2Voigt`](@ref) symmNo=10 order) of
the 15 (of 21) elastic constants surviving the 2D reduction. See
[`ElastRotReductMat`](@ref) (`ElastConstsDerived.jl`) for which Voigt pairs are dropped
and why.

Port of `guliwerr/core/@ElastConsts/GetIndOf2DReduct_15of21.m`.
"""
function GetIndOf2DReduct_15of21()
    indSet = [1, 3, 4, 5, 6, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    return indSet, length(indSet)
end

# -----------------------------------------------------------------------------------------
# Cubic/Hexagonal/Tetragonal (De)conversion pairs
# -----------------------------------------------------------------------------------------

"""
    CubicConversion(vecIn=Float64[]) -> (vecOut, transfMat)

Convert cubic elastic constants `[c11,c12,c44]` to `[B, c', c44]`, where `B` is the bulk
modulus and `c' = (c11-c12)/2`. Omit `vecIn` (or pass an empty vector) to get only the
3x3 transformation matrix, with `vecOut` an all-zeros placeholder ([`AsColumnVec`](@ref)).

Port of `guliwerr/core/@ElastConsts/CubicConversion.m`.
"""
function CubicConversion(vecIn = Float64[])
    v = AsColumnVec(vecIn, 3, "CubicConversion")
    transfMat = [1/3  2/3  0.0     # B
                 0.5 -0.5  0.0     # c'
                 0.0  0.0  1.0]    # c44
    return transfMat * v, transfMat
end

"""
    CubicDeconversion(vecIn=Float64[]) -> (vecOut, deTransfMat)

Inverse of [`CubicConversion`](@ref): recover `[c11,c12,c44]` from `[B,c',c44]`.

Port of `guliwerr/core/@ElastConsts/CubicDeconversion.m`.
"""
function CubicDeconversion(vecIn = Float64[])
    v = AsColumnVec(vecIn, 3, "CubicDeconversion")
    _, transfMat = CubicConversion()
    deTransfMat = inv(transfMat)
    return transfMat \ v, deTransfMat
end

"""
    HexagonalConversion(vecIn=Float64[]) -> (vecOut, transfMat)

Convert hexagonal elastic constants `[c11,c12,c13,c33,c44]` to a physical-mode basis
`[C11,C33,G1,G2,G3]`: uniaxial elongations `C11` (basal plane) and `C33` (c-axis), plus
shear moduli `G1 = c44` (axial shear), `G2 = (c11-c12)/2` (basal shear),
`G3 = (c11+c33-2*c13)/4` (breathing shear, volume-preserving c/a distortion). `c12` and
`c13` carry no mode of their own -- their Poisson coupling is absorbed into `G2`/`G3`.
`G1 = G2 = G3` in the isotropic limit.

Port of `guliwerr/core/@ElastConsts/HexagonalConversion.m`.
"""
function HexagonalConversion(vecIn = Float64[])
    v = AsColumnVec(vecIn, 5, "HexagonalConversion")
    transfMat = [1.0    0.0    0.0    0.0   0.0     # C11
                 0.0    0.0    0.0    1.0   0.0     # C33
                 0.0    0.0    0.0    0.0   1.0     # G1 == axial shear (c44)
                 0.5   -0.5    0.0    0.0   0.0     # G2 == basal shear
                 0.25   0.0   -0.5   0.25   0.0]    # G3 == breathing shear
    return transfMat * v, transfMat
end

"""
    HexagonalDeconversion(vecIn=Float64[]) -> (vecOut, deTransfMat)

Inverse of [`HexagonalConversion`](@ref): recover `[c11,c12,c13,c33,c44]` from
`[C11,C33,G1,G2,G3]`.

Port of `guliwerr/core/@ElastConsts/HexagonalDeconversion.m`.
"""
function HexagonalDeconversion(vecIn = Float64[])
    v = AsColumnVec(vecIn, 5, "HexagonalDeconversion")
    _, transfMat = HexagonalConversion()
    deTransfMat = inv(transfMat)
    return transfMat \ v, deTransfMat
end

"""
    TetragonalConversion(vecIn=Float64[]) -> (vecOut, transfMat)

Calculate equivalent tetragonal constants `[B, c33, c', c'', c44, c66]` from
`[c11,c12,c13,c33,c44,c66]` (tetragonal type I only), where `c''` is the volume-conserving
"breathing" mode.

Port of `guliwerr/core/@ElastConsts/TetragonalConversion.m`.
"""
function TetragonalConversion(vecIn = Float64[])
    v = AsColumnVec(vecIn, 6, "TetragonalConversion")
    transfMat = [2/9   2/9   4/9   1/9   0.0  0.0     # B
                 0.0   0.0   0.0   1.0   0.0  0.0     # c33
                 0.5  -0.5   0.0   0.0   0.0  0.0     # c'
                 1/6   1/6  -4/6   2/6   0.0  0.0     # c'' == breathing mode
                 0.0   0.0   0.0   0.0   1.0  0.0     # c44
                 0.0   0.0   0.0   0.0   0.0  1.0]    # c66
    return transfMat * v, transfMat
end

"""
    TetragonalDeconversion(vecIn=Float64[]) -> (vecOut, deTransfMat)

Inverse of [`TetragonalConversion`](@ref): recover `[c11,c12,c13,c33,c44,c66]` from
`[B,c33,c',c'',c44,c66]`.

Port of `guliwerr/core/@ElastConsts/TetragonalDeconversion.m`.
"""
function TetragonalDeconversion(vecIn = Float64[])
    v = AsColumnVec(vecIn, 6, "TetragonalDeconversion")
    _, transfMat = TetragonalConversion()
    deTransfMat = inv(transfMat)
    return transfMat \ v, deTransfMat
end

export SymmClassLibrary, FindSymmClass_input, IndepVoigtPositions, GetIndOf2DReduct_15of21,
       CubicConversion, CubicDeconversion, HexagonalConversion, HexagonalDeconversion,
       TetragonalConversion, TetragonalDeconversion
