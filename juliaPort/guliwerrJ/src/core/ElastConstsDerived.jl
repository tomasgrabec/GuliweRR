# =========================================================================================
# ElastConstsDerived.jl
#
# Port of guliwerr/core/@ElastConsts/RotateTensor.m, ElastRotReductMat.m,
# RightHandFrame.m, Tensor2Youngs.m, Voigt2Bulk.m, HillAverages.m, VoigtAverages.m,
# ReussAverages.m. `include`d directly into `module CoreTypes` alongside ElastConsts.jl
# and ElastConstsSymmetry.jl -- see the header comment at the top of ElastConsts.jl for
# the shared module-layout / naming-convention / immutability rationale.
#
# `using LinearAlgebra` here for `norm`, `dot`, `cross` (RightHandFrame, RotateTensor's
# orthogonality check) and `inv` (Tensor2Youngs, Voigt2Bulk, ReussAverages). Already a
# `[deps]` entry in Project.toml.
#
# `RotateTensor` and `Tensor2Youngs` are the two "hot" entry points from this file:
# `Vis.plotVelSurface`, `Vis.plotYoungsSurface` and `Vis.plotGroupVelSurface` each call the
# corresponding ElastConsts/Christoffel function inside a 61x121-point double loop (7381
# evaluations), so both are written to be type-stable (`@inferred`).
# =========================================================================================

using LinearAlgebra

# -----------------------------------------------------------------------------------------
# RightHandFrame.m
# -----------------------------------------------------------------------------------------

"""
    RightHandFrame(twoVecs, noteMsg=default) -> Matrix{Float64}

Build a right-hand orthonormal 3x3 frame (ROW vectors) from two direction vectors.

INPUT: `twoVecs` -- 2x3 or 3x2 matrix of two direction vectors; a 3x3 is also accepted and
only its first two ROWS are read (its third row is recomputed, so handedness is imposed,
not preserved).
`noteMsg` -- text of the `ElastConsts:RightHandFrame:notPerpendicular` warning issued when
the two input vectors are not perpendicular, `|dot(n1, n2)| > 1e-3` (MATLAB issues a
warning, for an obtuse pair too).

OUTPUT: 3x3 matrix of row vectors -- row 1 = `n1`, the normalised FIRST input vector
(preserved exactly); row 2 = `n2`, the second vector made exactly perpendicular to `n1`;
row 3 = `n3 = cross(n1,n2)`, the right-hand completion.

**ROLE-AGNOSTIC by design** (carried from the MATLAB docstring): callers disagree on what
the two vectors mean (`Material.CompleteDirMat` passes `[surface normal; propagation]`,
`RotateTensor` passes `[propagation; surface normal]`), so no role-specific behaviour is
added here.

Port of `guliwerr/core/@ElastConsts/RightHandFrame.m`.
"""
function RightHandFrame(twoVecs::AbstractMatrix{<:Real},
                         noteMsg::AbstractString =
                             "Coordinate system recalculated based on the first vector (input vectors not perpendicular).")
    tv = twoVecs
    if size(tv, 1) > size(tv, 2) && size(tv, 1) == 3
        tv = permutedims(tv)
    end
    if size(tv, 1) < 2 || size(tv, 2) != 3
        throw(ArgumentError("ElastConsts:RightHandFrame:badInput -- input must hold two " *
            "direction vectors as rows (2x3, 3x2 or 3x3); got $(size(tv,1))x$(size(tv,2))."))
    end

    n1 = Float64.(tv[1, :]) ./ norm(tv[1, :])
    n2 = Float64.(tv[2, :]) ./ norm(tv[2, :])

    if abs(dot(n1, n2)) > 1e-3   # inputs not perpendicular: n2 is about to be adjusted
        @warn "ElastConsts:RightHandFrame:notPerpendicular -- $noteMsg"
    end

    n3 = cross(n1, n2)
    n2 = cross(n3, n1)      # make n2 exactly perpendicular to n1

    return permutedims(hcat(n1, n2, n3))   # stack as ROWS
end

# -----------------------------------------------------------------------------------------
# RotateTensor.m
# -----------------------------------------------------------------------------------------

"""
    RotateTensor(inpTen, rotMat) -> Array{Float64,4}

Rotate an elastic tensor by a rotation matrix. `rotMat` may be 2x3/3x2 (two direction
vectors, completed to a right-hand frame via [`RightHandFrame`](@ref)) or 3x3 (used as
given, after an orthogonality check and per-row renormalisation).

`C'_{i1i2i3i4} = sum_{j1j2j3j4} R_{i1j1} R_{i2j2} R_{i3j3} R_{i4j4} C_{j1j2j3j4}`,
evaluated as four sequential mode-n matrix products -- this reshape/permute trick relies
on MATLAB and Julia sharing the SAME column-major array layout (unlike e.g. a NumPy port),
so it translates directly with no index-order correction.

Port of `guliwerr/core/@ElastConsts/RotateTensor.m`.
"""
function RotateTensor(inpTen::AbstractArray{<:Real,4}, rotMat::AbstractMatrix{<:Real})
    n = length(rotMat)
    local R::Matrix{Float64}
    if n == 6                                  # 2x3 or 3x2: complete to a right-hand frame
        R = RightHandFrame(rotMat,
            "Rotation matrix recalculated based on its first vector (input vectors not perpendicular).")
    elseif n == 9
        if size(rotMat) == (3, 3)
            Rm = Matrix{Float64}(rotMat)
            if any(abs.(Rm' * Rm - I) .> 1e-6)
                throw(ArgumentError("ElastConsts:RotateTensor:notOrthogonal -- Rotation matrix not orthogonal."))
            end
            for jj in 1:3
                Rm[jj, :] ./= norm(Rm[jj, :])
            end
            R = Rm
        else
            throw(ArgumentError("ElastConsts:RotateTensor:badFormat -- Unknown rotation matrix " *
                                 "format (set a 2x3 or 3x3 matrix)."))
        end
    else
        throw(ArgumentError("ElastConsts:RotateTensor:badFormat -- Unknown rotation matrix " *
                             "format (set a 2x3 or 3x3 matrix)."))
    end

    rotTen = Array{Float64}(inpTen)
    for d in 1:4
        perm = Int[d; setdiff(1:4, d)]     # bring dimension d to front, keep the rest ascending
        iperm = invperm(perm)
        rotTen = permutedims(rotTen, perm)
        rotTen = reshape(R * reshape(rotTen, 3, :), (3, 3, 3, 3))
        rotTen = permutedims(rotTen, iperm)
    end
    return rotTen
end

# -----------------------------------------------------------------------------------------
# ElastRotReductMat.m
# -----------------------------------------------------------------------------------------

"""
    ElastRotReductMat(rotMat, symmNo) -> Matrix{Float64}

Map from independent constants to the 15 reduced constants of the 2D RR (Rayleigh-Ritz)
system, for a given rotation. Output is `15 x noIndepConsts`.

The reduction drops every Voigt pair touching index 2 (the completion axis under the
standard x1=propagation, x2=completion, x3=surface-normal convention): `e_22` vanishes for
a field with no completion-axis dependence, so `C_I2`/`C_2J` terms add nothing to the
strain energy.

Port of `guliwerr/core/@ElastConsts/ElastRotReductMat.m`.
"""
function ElastRotReductMat(rotMat::AbstractMatrix{<:Real}, symmNo::Integer)
    indReducedEl, nReductElast = GetIndOf2DReduct_15of21()
    _, noIndepConstsRaw = SymmClassLibrary(symmNo)
    noIndepConsts = Int(noIndepConstsRaw)

    indep2RedElastMat = fill(NaN, nReductElast, noIndepConsts)

    for const_no in 1:noIndepConsts
        dummyConsts = zeros(Float64, noIndepConsts)
        dummyConsts[const_no] = 1.0
        cTen0   = Any2Tensor(dummyConsts, symmNo)
        cTenRot = RotateTensor(cTen0, rotMat)
        cVec    = Tensor2Vector(cTenRot, 10)

        indep2RedElastMat[:, const_no] = cVec[indReducedEl]
    end

    return indep2RedElastMat
end

# -----------------------------------------------------------------------------------------
# Tensor2Youngs.m
# -----------------------------------------------------------------------------------------

"""
    Tensor2Youngs(elTen, rotMat=nothing) -> Float64

Young's modulus along the x1 axis of the (optionally rotated) elastic tensor.

Port of `guliwerr/core/@ElastConsts/Tensor2Youngs.m`.
"""
function Tensor2Youngs(elTen::AbstractArray{<:Real,4}, rotMat = nothing)
    ten = rotMat === nothing ? elTen : RotateTensor(elTen, rotMat)
    elVoigt = Tensor2Voigt(ten)
    compVoigt = inv(elVoigt)   # compliance matrix; its (1,1) element is 1/E
    return 1 / compVoigt[1, 1]
end

# -----------------------------------------------------------------------------------------
# Voigt2Bulk.m
# -----------------------------------------------------------------------------------------

"""
    Voigt2Bulk(elVoigt) -> Float64

Bulk modulus from the Voigt matrix.

Port of `guliwerr/core/@ElastConsts/Voigt2Bulk.m`.
"""
function Voigt2Bulk(elVoigt::AbstractMatrix{<:Real})
    compVoigt = inv(elVoigt)   # compliance matrix; its 3x3 block sums to 1/B
    return 1 / sum(@view compVoigt[1:3, 1:3])
end

# -----------------------------------------------------------------------------------------
# VoigtAverages.m
# -----------------------------------------------------------------------------------------

"""
    VoigtAverages(elVoigt) -> (c11_voigt, c44_voigt)

Voigt averages -- moduli of the equivalent statistically isotropic material under the
Voigt (1928) uniform-STRAIN assumption, the upper VRH bound. `c11_voigt` is M (P-wave /
longitudinal modulus), `c44_voigt` is G (S-wave / transverse modulus). Kube & Turner 2015,
Eq. (11), with its typographical errors corrected.

Port of `guliwerr/core/@ElastConsts/VoigtAverages.m`.
"""
function VoigtAverages(elVoigt::AbstractMatrix{<:Real})
    c = elVoigt
    c11_voigt = (3c[1,1] + 2c[1,2] + 2c[1,3] + 3c[2,2] + 2c[2,3] + 3c[3,3] +
                 4c[4,4] + 4c[5,5] + 4c[6,6]) / 15
    c44_voigt = (c[1,1] - c[1,2] - c[1,3] + c[2,2] - c[2,3] + c[3,3] +
                 3c[4,4] + 3c[5,5] + 3c[6,6]) / 15
    return c11_voigt, c44_voigt
end

# -----------------------------------------------------------------------------------------
# ReussAverages.m
# -----------------------------------------------------------------------------------------

"""
    ReussAverages(elVoigt) -> (c11_reuss, c44_reuss)

Reuss averages -- moduli of the equivalent statistically isotropic material under the
Reuss (1929) uniform-STRESS assumption, the lower VRH bound. Standard general-symmetry
Reuss formulas for `K` and `G` from the compliance matrix `S = inv(elVoigt)`, combined via
`M = K + 4/3*G`. Reuss's `K` here is the same `1/sum(compliance 3x3 block)` formula as
[`Voigt2Bulk`](@ref).

Port of `guliwerr/core/@ElastConsts/ReussAverages.m`.
"""
function ReussAverages(elVoigt::AbstractMatrix{<:Real})
    s = inv(elVoigt)
    K_reuss = 1 / (s[1,1] + s[2,2] + s[3,3] + 2*(s[1,2] + s[1,3] + s[2,3]))
    G_reuss = 15 / (4*(s[1,1] + s[2,2] + s[3,3]) - 4*(s[1,2] + s[1,3] + s[2,3]) +
                     3*(s[4,4] + s[5,5] + s[6,6]))
    c11_reuss = K_reuss + 4/3*G_reuss
    c44_reuss = G_reuss
    return c11_reuss, c44_reuss
end

# -----------------------------------------------------------------------------------------
# HillAverages.m
# -----------------------------------------------------------------------------------------

"""
    HillAverages(elVoigt) -> (c11_hill, c44_hill)

Hill averages -- arithmetic mean of the Voigt and Reuss bounds (Hill 1952), the standard
practical estimate of a polycrystal's isotropic elastic moduli.

Port of `guliwerr/core/@ElastConsts/HillAverages.m`.
"""
function HillAverages(elVoigt::AbstractMatrix{<:Real})
    c11_voigt, c44_voigt = VoigtAverages(elVoigt)
    c11_reuss, c44_reuss = ReussAverages(elVoigt)
    return (c11_voigt + c11_reuss) / 2, (c44_voigt + c44_reuss) / 2
end

export RotateTensor, ElastRotReductMat, RightHandFrame, Tensor2Youngs, Voigt2Bulk,
       HillAverages, VoigtAverages, ReussAverages
