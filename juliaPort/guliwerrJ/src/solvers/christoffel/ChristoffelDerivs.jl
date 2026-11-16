# =========================================================================================
# ChristoffelDerivs.jl -- `include`d directly into `module Christoffel` (see
# Christoffel.jl's header for the "one module per folder" rationale and the
# `using ...CoreTypes` dependency this whole folder shares).
#
# Port of guliwerr/solvers/christoffel/@Christoffel/phaseVelAngleDeriv.m and
# groupVelVectorsForDir.m, plus Christoffel.m's private `GroupVelForDir`.
#
# Hot path: `phaseVelAngleDeriv` is the inner kernel of `LimitingVels.findRoots`' coarse
# scan (~314 calls per tensor) and of every fminbnd iteration. `_pvad_perm` below takes the
# permuted 9x9 tensor matrix (`_kronperm`) hoisted by the caller, so the 81-entry
# `permutedims` is done once per tensor instead of twice per angle.
# =========================================================================================

"""
    _kronperm(elTensor) -> Matrix{Float64}

The 9x9 matrix `GammaFromKron` contracts with `kron(n, n)`:
`reshape(permutedims(elTensor, (1, 3, 2, 4)), 9, 9)`. Hoist it out of any per-direction
loop over one tensor.
"""
_kronperm(elTensor::AbstractArray{<:Real,4}) =
    reshape(permutedims(Float64.(elTensor), (1, 3, 2, 4)), 9, 9)

"""
    _pvad_perm(M9, angleRad, rho) -> (vels, dVels, eigVecs)

[`phaseVelAngleDeriv`](@ref) on a hoisted [`_kronperm`](@ref) matrix; also returns the
eigenvectors (columns, descending), which `LimitingVels` needs for the root polarisation.
Same operations, in the same order, as MATLAB's `phaseVelAngleDeriv` -> `matrix` ->
`GammaFromKron` -> `SortedEig` chain (direction normalised for the eigenproblem only;
`dGamma` from the unnormalised `n`, `dn`, exactly as MATLAB).
"""
function _pvad_perm(M9::AbstractMatrix{Float64}, angleRad::Real, rho::Real)
    c, s = cos(angleRad), sin(angleRad)
    n  = (c, 0.0, s)
    dn = (-s, 0.0, c)
    nrm = sqrt(c * c + s * s)
    nh = (c / nrm, 0.0, s / nrm)            # matrix() normalises the direction

    nn  = Vector{Float64}(undef, 9)
    dnn = Vector{Float64}(undef, 9)
    @inbounds for b in 1:3, a in 1:3        # kron index (b-1)*3 + a  <->  x[b] * y[a]
        nn[(b - 1) * 3 + a]  = nh[b] * nh[a]
        dnn[(b - 1) * 3 + a] = dn[b] * n[a] + n[b] * dn[a]
    end
    G  = reshape(M9 * nn, 3, 3) ./ rho
    dG = reshape(M9 * dnn, 3, 3) ./ rho

    eigVals, eigVecs = SortedEig(G)
    vels = sqrt.(eigVals)
    dVels = Vector{Float64}(undef, 3)
    @inbounds for iwt in 1:3
        g = @view eigVecs[:, iwt]
        dVels[iwt] = dot(g, dG, g) / (2 * vels[iwt])
    end
    return vels, dVels, eigVecs
end

# -----------------------------------------------------------------------------------------
# GroupVelForDir -- MATLAB `methods (Static, Access = private)` (Christoffel.m)
# -----------------------------------------------------------------------------------------

"""
    GroupVelForDir(normC, normCperm, n) -> (vG, vP, eigVecs)

Phase/group velocities and polarisations for ONE direction.

INPUT: `normC` -- density-NORMALISED 3x3x3x3 elastic tensor; `normCperm` --
`reshape(permutedims(normC,(1,4,2,3)), 9, 9)`, hoisted by the caller (no direction
dependence); `n` -- propagation direction, 3-vector.
OUTPUT: `vG` -- 3x3, columns are group velocity vectors `[L,T1,T2]`; `vP` -- 3-vector phase
velocities, descending; `eigVecs` -- 3x3, columns are polarisation eigenvectors.

`vG_i = C(i,j,k,l)*g(j)*g(k)*s(l)`, with `s = n/vP` the slowness vector (Royer-Dieulesaint,
*Elastic Waves in Solids I*, Eq. 4.24).

Port of `guliwerr/solvers/christoffel/@Christoffel/Christoffel.m` (`GroupVelForDir`).
"""
function GroupVelForDir(normC::AbstractArray{<:Real,4}, normCperm::AbstractMatrix{<:Real},
                         n::AbstractVector{<:Real})
    eigVals, eigVecs, _ = matrix(normC, n, 1.0)
    vP = sqrt.(eigVals)

    vG = zeros(Float64, 3, 3)
    for iwt in 1:3
        g = eigVecs[:, iwt]
        s = n ./ vP[iwt]                        # slowness vector
        P = g * g'                              # 3x3 outer product g⊗g
        T = reshape(normCperm * vec(P), 3, 3)   # T(i,l) = sum_jk normC(i,j,k,l)*P(j,k)
        vG[:, iwt] = T * s
    end
    return vG, vP, eigVecs
end

# -----------------------------------------------------------------------------------------
# groupVelVectorsForDir.m
# -----------------------------------------------------------------------------------------

"""
    groupVelVectorsForDir(elTensor, dirVec, rho=1.0) -> (vG, vP, eigVecsOut)

Bulk group-velocity vectors for an arbitrary 3D direction (e.g. full-sphere ray surfaces).

INPUT: `elTensor` -- 3x3x3x3 tensor; `dirVec` -- propagation direction (need not be
normalised, must be nonzero); `rho` -- density (default `1`, positive finite).
OUTPUT: `vG` -- 3x3, columns are group velocity vectors `[qL, qT1, qT2]`; `vP` -- 3-vector
phase velocities, descending; `eigVecsOut` -- 3x3, columns are polarisation eigenvectors.

Royer-Dieulesaint Eq. 4.24 contraction, `vG_i = C(i,j,k,l) g(j) g(k) s(l) / rho` (`g` the
polarisation, `s = n / vP` the slowness). See `docs/reference/architecture/christoffel.md`.

Errors (MATLAB ids): `Christoffel:groupVelVectorsForDir:badRho`,
`Christoffel:groupVelVectorsForDir:badDirVec`.

Port of `guliwerr/solvers/christoffel/@Christoffel/groupVelVectorsForDir.m`.
"""
function groupVelVectorsForDir(elTensor::AbstractArray{<:Real,4}, dirVec::AbstractVector{<:Real},
                                rho::Real = 1.0)
    (isfinite(rho) && rho > 0) ||
        throw(ArgumentError("Christoffel:groupVelVectorsForDir:badRho -- rho must be a positive finite scalar."))
    (length(dirVec) == 3 && any(!iszero, dirVec)) ||
        throw(ArgumentError("Christoffel:groupVelVectorsForDir:badDirVec -- dirVec must be a nonzero 3-vector."))
    normC = elTensor ./ rho
    n = Float64.(vec(dirVec)) ./ norm(dirVec)
    normC_perm = reshape(permutedims(normC, (1, 4, 2, 3)), 9, 9)
    return GroupVelForDir(normC, normC_perm, n)
end

# -----------------------------------------------------------------------------------------
# phaseVelAngleDeriv.m
# -----------------------------------------------------------------------------------------

"""
    phaseVelAngleDeriv(elTensor, angleRad, rho=1.0) -> (vels, dVels)

Phase velocities and their first derivative w.r.t. sagittal-plane angle (direction
`[cos; 0; sin]`).

`angleRad` is IN RADIANS, unlike [`phaseVelScan`](@ref)'s degrees; angle-derivative
formulas are radian-native, and mixing the two costs a silent factor of `180/pi`.

OUTPUT: `vels` -- `[vL; vT1; vT2]`, descending; `dVels` -- `d(vels)/d(angleRad)`, same order
(Hellmann-Feynman: `dv = g' dGamma g / (2 v)`). Derivation: `docs/reference/architecture/christoffel.md`.

Port of `guliwerr/solvers/christoffel/@Christoffel/phaseVelAngleDeriv.m`.
"""
function phaseVelAngleDeriv(elTensor::AbstractArray{<:Real,4}, angleRad::Real, rho::Real = 1.0)
    size(elTensor) == (3, 3, 3, 3) ||
        throw(ArgumentError("Christoffel:phaseVelAngleDeriv:badInput -- elTensor must be a 3x3x3x3 tensor"))
    vels, dVels, _ = _pvad_perm(_kronperm(elTensor), angleRad, rho)
    return vels, dVels
end

export phaseVelAngleDeriv, groupVelVectorsForDir
