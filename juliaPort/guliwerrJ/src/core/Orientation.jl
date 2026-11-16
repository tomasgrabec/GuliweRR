# =========================================================================================
# Orientation.jl
#
# Port of the 12 free-standing files under guliwerr/core/@Material/ that are NOT Material.m
# itself: Euler2DCM.m, DCM2Euler.m, Euler2DirMat.m, DirMat2Euler.m, Miller2DirMat.m,
# DirMat2Miller.m, CompleteDirMat.m, GetStandardDirMat.m, RotDirMatFromAngle.m,
# CalcRotMat.m, LatticeParamSpec.m, LatticeParamSpecReport.m. `include`d directly into
# `module CoreTypes` (see CoreTypes.jl and Material.jl's header for why there is no
# `module Orientation` wrapper), after Material.jl --
# though nothing here actually requires that order: every cross-reference between the two
# files is a plain function call, resolved at call time, and neither file's `struct`
# depends on a type defined in the other (only struct field types must be defined before
# use).
#
# All of MATLAB's `cosd`/`sind`/`acosd`/`atan2d` map directly onto Julia Base's own
# degree-argument trig functions (`cosd`, `sind`, `acosd`, and the two-argument `atand(y,x)`
# for `atan2d`) -- no hand-rolled degree conversion needed anywhere in this file.
#
# NAMING NOTE -- same convention as ElastConsts.jl/Material.jl: every function
# name below is kept IDENTICAL to its MATLAB original (`Euler2DCM`, `DCM2Euler`, ...).
#
# PRIVACY NOTE: `LatticeParamSpec`/`LatticeParamSpecReport` are
# `methods (Static, Access = private)` in MATLAB (Material.m). Julia has no access
# control; following the AsColumnVec precedent in ElastConsts.jl, both
# are simply left UNEXPORTED rather than renamed -- reachable from outside CoreTypes via the
# fully-qualified `GuliwerrJ.CoreTypes.LatticeParamSpec` path, as `AsColumnVec`/`FindType`
# are.
#
# TRANSLATION DECISIONS (language-shape differences, not behaviour changes):
# 1. `Miller2DirMat`/`DirMat2Miller`'s optional `latticeParams` MATLAB `struct` (fields
#    b_over_a, c_over_a, alpha, beta, gamma) becomes a Julia `NamedTuple`, e.g.
#    `(c_over_a = 1.633,)` for `struct('c_over_a',1.633)` -- the natural Julia analogue of
#    MATLAB's ad-hoc `struct()`. `LatticeParamSpec`'s `spec.defaults`/the merged parameter
#    set `p` are NamedTuples for the same reason; `spec.required` stays a `Vector{String}` (not
#    `Vector{Symbol}`) so `spec.required == ["b_over_a","c_over_a"]` reads the same as
#    MATLAB's `{'b_over_a','c_over_a'}`, and so `join(spec.required, ", ")` needs
#    no `String.()` conversion.
# 2. `Miller2DirMat`'s 1x3 ROW-vector output (`dcm`) and `DirMat2Miller`'s row-vector input
#    (`dirCos`) are plain `Vector{Float64}` here, not a 1x3 `Matrix` -- MATLAB's row/column
#    distinction for a single direction-cosine triple carries no meaning once ported (every
#    consumer here is `dot`/`cross`/`norm`, which don't care), and a bare `Vector` is what
#    every other "one direction" quantity here is (`RightHandFrame`'s rows,
#    once sliced out, `RotDirMatFromAngle`'s row arithmetic).
# 3. Both requirements-listing entry points (`Miller2DirMat(symmNo)` /
#    `Miller2DirMat(nothing, symmNo)`, and the `DirMat2Miller` equivalents) are reachable in
#    Julia via the SAME multi-purpose function as the real computation, exactly mirroring
#    MATLAB's own `nargin`/`isempty`-gated dual behaviour (one function, two unrelated return
#    shapes depending on input). This is a genuine, intentional type instability -- see the
#    docstrings below and the precedent set by `Any2Tensor` in ElastConsts.jl:
#    ported faithfully, and deliberately NOT type-stable (`@inferred`).
# =========================================================================================

using LinearAlgebra

_is_scalarlike(x) = x isa Number || (x isa AbstractArray && length(x) == 1)
_as_scalar(x) = x isa Number ? x : first(x)

# -----------------------------------------------------------------------------------------
# Euler2DCM.m
# -----------------------------------------------------------------------------------------

"""
    Euler2DCM(phi1, Phi, phi2) -> Matrix{Float64}

Build a direction-cosine matrix from Bunge (ZXZ) Euler angles [degrees].

Row/column 3 is Bunge's own distinguished ("tilt") axis -- NOT GULIWERR's surface normal.
This is the raw, generic angle<->matrix conversion with no [`Material`](@ref) `dirMat`
semantics at all. Use [`Euler2DirMat`](@ref) instead for a ready-to-use `Material.dirMat`
(row 1 = surface normal).

Exact inverse of [`DCM2Euler`](@ref) (see there for the sign/frame convention).

Port of `guliwerr/core/@Material/Euler2DCM.m`.
"""
function Euler2DCM(phi1::Real, Phi::Real, phi2::Real)
    c1 = cosd(phi1); s1 = sind(phi1)
    C  = cosd(Phi);  S  = sind(Phi)
    c2 = cosd(phi2); s2 = sind(phi2)

    g = Float64[c1*c2 - s1*s2*C   s1*c2 + c1*s2*C   s2*S
               -c1*s2 - s1*c2*C  -s1*s2 + c1*c2*C   c2*S
                s1*S             -c1*S              C   ]

    return Matrix(g')
end

# -----------------------------------------------------------------------------------------
# DCM2Euler.m
# -----------------------------------------------------------------------------------------

"""
    DCM2Euler(dcm) -> (phi1, Phi, phi2)

Extract Bunge (ZXZ) Euler angles [degrees] from a direction-cosine matrix. `dcm` is a 3x3
matrix (rows normalised internally, handedness preserved AS GIVEN -- rows are read
directly, never reconstructed via a right-hand cross product), or a 2x3/3x2 vector pair
`[n1; n2]`, right-hand completed via [`CompleteDirMat`](@ref).

`phi1 in [0,360)`, `Phi in [0,180]`.

CONVENTION: `g = Rz(phi2)*Rx(Phi)*Rz(phi1)` transforms SPECIMEN-frame components into
CRYSTAL-frame components; `dcm = g'`, so `dcm`'s rows are the specimen axes in crystal
components. Row/column 3 is Bunge's distinguished ("tilt") axis.

This is the raw, generic matrix<->angle conversion with no [`Material`](@ref) `dirMat`
semantics at all. Use [`DirMat2Euler`](@ref) instead if you have a genuine `Material.dirMat`
(row 1 = surface normal).

Do NOT call this on [`GetStandardDirMat`](@ref)'s output: that `[3,2,1]` reorder is an odd
permutation, giving an improper (`det = -1`) matrix, which has no Bunge decomposition at
all. A number still comes back, but it is not a physical rotation and [`Euler2DCM`](@ref)
will not round-trip it.

Same convention as MTEX (`orientation.byEuler(...).matrix == dcm`) and as the Nottingham
SRAS orientation output.

Port of `guliwerr/core/@Material/DCM2Euler.m`.
"""
function DCM2Euler(dcm::AbstractMatrix{<:Real})
    d = Matrix{Float64}(dcm)
    if size(d, 1) == 2 || size(d, 2) == 2
        d = CompleteDirMat(d)
    else
        for jj in 1:3
            d[jj, :] ./= norm(d[jj, :])
        end
    end

    tol = 1e-9
    Phi = acosd(d[3, 3])
    phi1, phi2 = if abs(sind(Phi)) < tol
        # gimbal lock (Phi = 0 or 180 deg): only phi1+/-phi2 is determined, so by
        # convention the whole rotation goes into phi1
        atand(d[2, 1], d[1, 1]), 0.0
    else
        atand(d[1, 3], -d[2, 3]), atand(d[3, 1], d[3, 2])
    end

    return mod(phi1, 360), Phi, mod(phi2, 360)
end

# -----------------------------------------------------------------------------------------
# Euler2DirMat.m
# -----------------------------------------------------------------------------------------

"""
    Euler2DirMat(phi1, Phi, phi2) -> Matrix{Float64}

Build a `Material.dirMat` orientation from Bunge (ZXZ) Euler angles [degrees]: 3x3
Convention-A direction matrix (row 1 = surface normal, row 2 = in-plane propagation
reference, row 3 = completion) -- exact inverse of [`DirMat2Euler`](@ref), ready to pass
straight into [`Material`](@ref)'s constructor.

Builds Bunge's own direction-cosine matrix (row/col 3 = his distinguished tilt axis, via
[`Euler2DCM`](@ref)) and reorders it `[3,1,2]` onto GULIWERR's row-1-is-surface-normal
convention. Need the raw, unreordered form instead? Use [`Euler2DCM`](@ref).

Port of `guliwerr/core/@Material/Euler2DirMat.m`.
"""
function Euler2DirMat(phi1::Real, Phi::Real, phi2::Real)
    dcm = Euler2DCM(phi1, Phi, phi2)
    return dcm[[3, 1, 2], :]
end

# -----------------------------------------------------------------------------------------
# DirMat2Euler.m
# -----------------------------------------------------------------------------------------

"""
    DirMat2Euler(dirMat) -> (phi1, Phi, phi2)

Extract Bunge (ZXZ) Euler angles [degrees] from a `Material.dirMat` orientation.
`dirMat` is a 3x3 Convention-A direction matrix (row 1 = surface normal), or a 2x3/3x2
vector pair, right-hand completed via [`CompleteDirMat`](@ref). Exact inverse of
[`Euler2DirMat`](@ref).

Right-hand completes the input (a no-op for a `dirMat` that is already a proper
right-handed frame, which every genuine `Material.dirMat` is) then reorders it `[2,3,1]`
onto Bunge's own row/col-3-is-tilt-axis frame before extracting angles ([`DCM2Euler`](@ref)).

Port of `guliwerr/core/@Material/DirMat2Euler.m`.
"""
function DirMat2Euler(dirMat::AbstractMatrix{<:Real})
    dm = CompleteDirMat(dirMat)
    return DCM2Euler(dm[[2, 3, 1], :])
end

# -----------------------------------------------------------------------------------------
# CompleteDirMat.m
# -----------------------------------------------------------------------------------------

"""
    CompleteDirMat(inp) -> Matrix{Float64}

Build a right-hand orthonormal coordinate system from two direction vectors, via
[`RightHandFrame`](@ref).

INPUT: `inp` -- 2x3 or 3x2 matrix containing two row vectors `[n1; n2]` (a 3x3 is accepted
too -- only its first two rows are read, so its third row is REBUILT right-handed, not
preserved).
OUTPUT: 3x3 matrix of row vectors: row 1 -- `n1` (preserved); row 2 -- `n2` (adjusted to be
perpendicular to `n1`); row 3 -- `n3 = cross(n1,n2)` (right-hand completion).

Port of `guliwerr/core/@Material/CompleteDirMat.m`.
"""
function CompleteDirMat(inp::AbstractMatrix{<:Real})
    return RightHandFrame(inp, "dirMat recalculated based on its first vector (input vectors not perpendicular).")
end

# -----------------------------------------------------------------------------------------
# GetStandardDirMat.m
# -----------------------------------------------------------------------------------------

"""
    GetStandardDirMat(dirMat, angleDeg=0.0) -> Matrix{Float64}

Reorder (and optionally rotate) a direction matrix into the standard x1=propagation,
x2=in-plane, x3=surface-normal convention of the literature (Every, Maznev), from the
internal row order `[plane; direction; completion]`.

OUTPUT: 3x3 direction matrix, rows `[completion; direction; plane]`.

The angle is applied here, BEFORE the reorder, because [`RotDirMatFromAngle`](@ref) mixes
rows 2 & 3 -- rotating an already-reordered matrix would swing the propagation direction
towards the surface normal. So pass angle 0 downstream, e.g. `CalcRotMat(stdDirMat, 0)`.

Port of `guliwerr/core/@Material/GetStandardDirMat.m`.
"""
function GetStandardDirMat(dirMat::AbstractMatrix{<:Real}, angleDeg::Real = 0.0)
    rotDirMat = RotDirMatFromAngle(dirMat, angleDeg)
    return rotDirMat[[3, 2, 1], :]
end

# -----------------------------------------------------------------------------------------
# RotDirMatFromAngle.m
# -----------------------------------------------------------------------------------------

"""
    RotDirMatFromAngle(dirMat, angleDeg) -> Matrix{Float64}

Rotate a direction matrix in-plane by the given angle [degrees]. `dirMat` is a 3x3 direction
matrix (row vectors; complete right-hand system).

OUTPUT: 3x3 direction matrix, rows 2 and 3 rotated by `angleDeg` from row 2 towards row 3
(counter-clockwise seen looking down the surface normal, row 1, onto the surface); with
`Euler2DirMat` this is from specimen X towards Y. Row 1 is the pivot and stays unchanged. Sole definition of the in-plane rotation's sign
convention.

Port of `guliwerr/core/@Material/RotDirMatFromAngle.m`.
"""
function RotDirMatFromAngle(dirMat::AbstractMatrix{<:Real}, angleDeg::Real)
    c = cosd(angleDeg); s = sind(angleDeg)
    row1 = Float64.(dirMat[1, :])
    row2 = dirMat[2, :] .* c .+ dirMat[3, :] .* s
    row3 = dirMat[2, :] .* (-s) .+ dirMat[3, :] .* c
    return permutedims(hcat(row1, row2, row3))   # stack as ROWS
end

# -----------------------------------------------------------------------------------------
# CalcRotMat.m
# -----------------------------------------------------------------------------------------

"""
    CalcRotMat(dirMat, angleDeg=0.0) -> Matrix{Float64}

Build a rotation matrix (3x2) from a direction matrix and an angle [degrees], for
[`RotateTensor`](@ref): column 1 -- propagation direction, rotated within the surface plane
(the plane spanned by `dirMat` rows 2 and 3); column 2 -- surface normal (`dirMat` row 1),
the pivot, unchanged.

Port of `guliwerr/core/@Material/CalcRotMat.m`.
"""
function CalcRotMat(dirMat::AbstractMatrix{<:Real}, angleDeg::Real = 0.0)
    rotDirMat = RotDirMatFromAngle(dirMat, angleDeg)
    return hcat(rotDirMat[2, :], rotDirMat[1, :])
end

# -----------------------------------------------------------------------------------------
# LatticeParamSpec.m / LatticeParamSpecReport.m -- private in MATLAB, unexported here
# -----------------------------------------------------------------------------------------

"""
    LatticeParamSpecResult

Return type of [`LatticeParamSpec`](@ref): `defaults` (a `NamedTuple` with all of
`b_over_a, c_over_a, alpha, beta, gamma` at their symmetry-implied values), `required`
(the field names the caller must supply, as `String`s -- see Orientation.jl's TRANSLATION
DECISIONS header note for why `String` and not `Symbol`), and `className` (the symmetry
class name, from [`SymmClassLibrary`](@ref)).
"""
struct LatticeParamSpecResult
    defaults::NamedTuple{(:b_over_a, :c_over_a, :alpha, :beta, :gamma), NTuple{5,Float64}}
    required::Vector{String}
    className::String
end

"""
    LatticeParamSpec(symmNo) -> LatticeParamSpecResult

Default lattice-shape parameters and which are required, by symmetry class.

CONVENTIONS (carried from MATLAB): trigonal (5/6) in the hexagonal axis setting, not
rhombohedral; monoclinic (9) unique-axis-b (`alpha = gamma = 90`, `beta` free).

Port of `guliwerr/core/@Material/LatticeParamSpec.m`. Private in MATLAB
(`methods (Static, Access = private)`); left unexported here, per this file's header note.
"""
function LatticeParamSpec(symmNo::Integer)
    defaults = (b_over_a = 1.0, c_over_a = 1.0, alpha = 90.0, beta = 90.0, gamma = 90.0)
    required = if symmNo in (1, 2)             # isotropic, cubic
        String[]
    elseif symmNo in (3, 4)                    # tetragonal, type I/II
        ["c_over_a"]
    elseif symmNo in (5, 6, 7)                 # trigonal I/II (hex-axis setting), hexagonal
        defaults = merge(defaults, (gamma = 120.0,))
        ["c_over_a"]
    elseif symmNo == 8                         # orthorhombic
        ["b_over_a", "c_over_a"]
    elseif symmNo == 9                         # monoclinic, unique axis b
        ["b_over_a", "c_over_a", "beta"]
    elseif symmNo == 10                        # triclinic
        ["b_over_a", "c_over_a", "alpha", "beta", "gamma"]
    else
        throw(ArgumentError("Material:LatticeParamSpec:badSymmNo -- Unknown symmetry " *
            "class number $symmNo (see SymmClassLibrary)."))
    end
    className = SymmClassLibrary(symmNo)[1]
    return LatticeParamSpecResult(defaults, required, className)
end

"""
    LatticeParamSpecReport(fname, symmNo) -> LatticeParamSpecResult

Print (via `println`, in place of MATLAB's `disp`) and return which `latticeParams` a
symmetry class requires. `fname` names the calling function, for the printed prefix.

Backs the requirements-listing mode of [`Miller2DirMat`](@ref) and [`DirMat2Miller`](@ref).

Port of `guliwerr/core/@Material/LatticeParamSpecReport.m`. Private in MATLAB; left
unexported here, per this file's header note.
"""
function LatticeParamSpecReport(fname::AbstractString, symmNo::Integer)
    spec = LatticeParamSpec(symmNo)
    msg = if isempty(spec.required)
        "$(spec.className): no latticeParams needed (angles/ratios are fixed by symmetry)."
    else
        "$(spec.className): latticeParams must include: $(join(spec.required, ", "))"
    end
    println("!(N) Material.$fname: $msg")
    return spec
end

# -----------------------------------------------------------------------------------------
# Miller2DirMat.m
# -----------------------------------------------------------------------------------------

"""
    Miller2DirMat(hkl, symmNo=nothing, latticeParams=nothing, indexType=nothing) -> Vector{Float64}
    Miller2DirMat(symmNo) -> LatticeParamSpecResult
    Miller2DirMat(nothing, symmNo) -> LatticeParamSpecResult

Convert a Miller (or Miller-Bravais) index to a direction-cosine vector, OR (requirements-
listing mode, when `hkl` is scalar-like and `symmNo` is omitted, or when `hkl` is `nothing`
or empty) list the `latticeParams` fields `symmNo` requires.

**Genuinely dual-typed return** (`Vector{Float64}` in the normal case, a
[`LatticeParamSpecResult`](@ref) in listing mode), mirroring MATLAB's own `nargin`/`isempty`
dynamic dispatch exactly -- see this file's header TRANSLATION DECISIONS note 3. Not
type-stable (`@inferred`), same precedent as `Any2Tensor` in
ElastConsts.jl.

INPUT: `hkl` -- Miller index `[h,k,l]`, or Miller-Bravais `[h,k,i,l]` (`i = -(h+k)`;
trigonal/hexagonal PLANE NORMALS only); `symmNo` -- symmetry class number (see
[`SymmClassLibrary`](@ref)); `latticeParams` -- `NamedTuple` of lattice-shape parameters
(ratios and angles only -- absolute scale never affects a direction cosine; call with
`hkl = nothing` to list the fields `symmNo` requires; unrequired fields are ignored):
`b_over_a, c_over_a` (axis-length ratios), `alpha, beta, gamma` (interaxial angles
[degrees]: `alpha = angle(b,c)`, `beta = angle(a,c)`, `gamma = angle(a,b)`); `indexType` --
`"direction"` (real-space `[uvw]`, DEFAULT) or `"planeNormal"` (reciprocal-space `(hkl)`).
NOT parallel except for cubic symmetry -- a surface normal is a planeNormal, a propagation
direction is a direction.

Port of `guliwerr/core/@Material/Miller2DirMat.m`.
"""
function Miller2DirMat(hkl, symmNo::Union{Nothing,Real} = nothing,
                        latticeParams::Union{Nothing,NamedTuple} = nothing,
                        indexType::Union{Nothing,AbstractString} = nothing)
    # requirements-listing mode: Miller2DirMat(symmNo) or Miller2DirMat(nothing, symmNo)
    if symmNo === nothing
        _is_scalarlike(hkl) || throw(ArgumentError("Material:Miller2DirMat:needSymmNo -- " *
            "Give a symmetry class number: Miller2DirMat(symmNo) or " *
            "Miller2DirMat(nothing, symmNo)."))
        return LatticeParamSpecReport("Miller2DirMat", Int(_as_scalar(hkl)))
    end
    if hkl === nothing || isempty(hkl)
        return LatticeParamSpecReport("Miller2DirMat", Int(symmNo))
    end

    lp    = latticeParams === nothing ? NamedTuple() : latticeParams
    itype = indexType === nothing ? "direction" : indexType

    spec = LatticeParamSpec(Int(symmNo))
    for f in spec.required
        sym = Symbol(f)
        if !(haskey(lp, sym) && getfield(lp, sym) !== nothing)
            throw(ArgumentError("Material:Miller2DirMat:missingParam -- $(spec.className) " *
                "requires latticeParams.$f (missing). Full requirement list: " *
                "$(join(spec.required, ", "))"))
        end
    end
    p = spec.defaults
    for f in spec.required
        sym = Symbol(f)
        p = merge(p, NamedTuple{(sym,)}((Float64(getfield(lp, sym)),)))
    end

    hklv = Float64.(vec(hkl))

    # Miller-Bravais 4-index [h k i l] -> 3-index [h k l]; plane normals only, since the
    # real-space zone-axis form [UVTW] is a different transform and is not implemented
    if length(hklv) == 4
        if lowercase(itype) == "direction"
            throw(ArgumentError("Material:Miller2DirMat:zoneAxisNotSupported -- 4-index " *
                "Miller-Bravais [h k i l] is a PLANE-NORMAL notation; the equivalent " *
                "4-index zone-axis notation for a DIRECTION ([UVTW]) uses a different, " *
                "non-trivial transform not implemented here. Supply the 3-index direction " *
                "[u v w] instead."))
        end
        if abs(hklv[3] + hklv[1] + hklv[2]) > 1e-9
            throw(ArgumentError("Material:Miller2DirMat:badMillerBravais -- Miller-Bravais " *
                "index must satisfy i = -(h+k); got h=$(hklv[1]), k=$(hklv[2]), " *
                "i=$(hklv[3])."))
        end
        hklv = hklv[[1, 2, 4]]
    end

    # direct lattice vectors; a fixed to 1, as absolute scale cannot affect a direction cosine
    a = 1.0; b = p.b_over_a; c = p.c_over_a
    al = deg2rad(p.alpha); be = deg2rad(p.beta); ga = deg2rad(p.gamma)

    aVec = [a, 0.0, 0.0]
    bVec = [b * cos(ga), b * sin(ga), 0.0]
    cx = c * cos(be)
    cy = c * (cos(al) - cos(be) * cos(ga)) / sin(ga)
    cz2 = c^2 - cx^2 - cy^2
    if cz2 < -1e-9
        throw(ArgumentError("Material:Miller2DirMat:badLatticeParams -- latticeParams do " *
            "not describe a valid unit cell (negative cell volume^2)."))
    end
    cVec = [cx, cy, sqrt(max(cz2, 0.0))]

    it = lowercase(itype)
    vCart = if it == "direction"
        hklv[1] * aVec + hklv[2] * bVec + hklv[3] * cVec
    elseif it == "planenormal"
        V = dot(aVec, cross(bVec, cVec))
        aStar = cross(bVec, cVec) / V
        bStar = cross(cVec, aVec) / V
        cStar = cross(aVec, bVec) / V
        hklv[1] * aStar + hklv[2] * bStar + hklv[3] * cStar
    else
        throw(ArgumentError("Material:Miller2DirMat:badIndexType -- indexType must be " *
            "'direction' or 'planeNormal'."))
    end

    nrm = norm(vCart)
    if nrm < 1e-12
        throw(ArgumentError("Material:Miller2DirMat:zeroVector -- hkl maps to the zero vector."))
    end
    return vCart ./ nrm
end

# -----------------------------------------------------------------------------------------
# DirMat2Miller.m
# -----------------------------------------------------------------------------------------

"""
    DirMat2Miller(dirCos, symmNo=nothing, latticeParams=nothing, indexType=nothing, maxIndex=nothing)
        -> (hkl, residualDeg)
    DirMat2Miller(symmNo) -> (LatticeParamSpecResult, Float64[])
    DirMat2Miller(nothing, symmNo) -> (LatticeParamSpecResult, Float64[])

Find the smallest-integer Miller index whose direction best matches a direction cosine, by
brute force over reduced integer triples scored against [`Miller2DirMat`](@ref) (exact for
all 7 crystal systems, but `O(maxIndex^3)` -- not for a hot loop).

Indexing has no exact general solution -- an arbitrary direction need not lie along any
low-index crystallographic direction. **ALWAYS check `residualDeg`.**

Requirements-listing mode (see [`Miller2DirMat`](@ref)'s docstring for the same pattern)
always returns a 2-tuple here too, `(spec, Float64[])`, so the normal and listing-mode calls
share one return shape (unlike MATLAB, where `nargout` silently drops the second output) --
callers wanting just the spec write `spec, _ = DirMat2Miller(symmNo)`.

INPUT: `dirCos` -- direction-cosine vector (need not be pre-normalised); `symmNo`,
`latticeParams`, `indexType` -- see [`Miller2DirMat`](@ref) (call with `dirCos = nothing` to
see `latticeParams` requirements); `maxIndex` -- search `|h|,|k|,|l|` up to this integer
[default 12].
OUTPUT: `hkl` -- best-fit small-integer Miller index `[h,k,l]`, or Miller-Bravais
`[h,k,i,l]` (`i=-(h+k)`) for trigonal/hexagonal with `indexType = "planeNormal"`;
`residualDeg` -- angular misfit [degrees] between `hkl`'s direction and `dirCos`
(`0` = exact match).

Port of `guliwerr/core/@Material/DirMat2Miller.m`.
"""
function DirMat2Miller(dirCos, symmNo::Union{Nothing,Real} = nothing,
                        latticeParams::Union{Nothing,NamedTuple} = nothing,
                        indexType::Union{Nothing,AbstractString} = nothing,
                        maxIndex::Union{Nothing,Integer} = nothing)
    if symmNo === nothing
        _is_scalarlike(dirCos) || throw(ArgumentError("Material:DirMat2Miller:needSymmNo -- " *
            "Give a symmetry class number: DirMat2Miller(symmNo) or " *
            "DirMat2Miller(nothing, symmNo)."))
        return LatticeParamSpecReport("DirMat2Miller", Int(_as_scalar(dirCos))), Float64[]
    end
    if dirCos === nothing || isempty(dirCos)
        return LatticeParamSpecReport("DirMat2Miller", Int(symmNo)), Float64[]
    end

    lp    = latticeParams === nothing ? NamedTuple() : latticeParams
    itype = indexType === nothing ? "direction" : indexType
    maxIdx = maxIndex === nothing ? 12 : Int(maxIndex)

    dv = Float64.(vec(dirCos))
    dv = dv ./ norm(dv)

    bestDot = -Inf
    bestH = [1, 0, 0]
    for h in -maxIdx:maxIdx, k in -maxIdx:maxIdx, l in -maxIdx:maxIdx
        (h == 0 && k == 0 && l == 0) && continue
        gcd(gcd(abs(h), abs(k)), abs(l)) > 1 && continue   # skip non-reduced triples
        cand = Miller2DirMat([h, k, l], Int(symmNo), lp, itype)
        d = dot(cand, dv)
        if d > bestDot
            bestDot = d
            bestH = [h, k, l]
        end
    end

    residualDeg = acosd(min(bestDot, 1.0))
    hkl = bestH
    if Int(symmNo) in (5, 6, 7) && lowercase(itype) == "planenormal"
        hkl = [hkl[1], hkl[2], -(hkl[1] + hkl[2]), hkl[3]]
    end
    return hkl, residualDeg
end

export Euler2DCM, DCM2Euler, Euler2DirMat, DirMat2Euler, Miller2DirMat, DirMat2Miller,
       CompleteDirMat, GetStandardDirMat, RotDirMatFromAngle, CalcRotMat
