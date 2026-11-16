# =========================================================================================
# Christoffel.jl
#
# Port of guliwerr/solvers/christoffel/@Christoffel/ (matrix.m, GammaFromKron -- a private
# static method of Christoffel.m itself --, phaseVels.m, phaseVelScan.m). The other two
# @Christoffel/ files (phaseVelAngleDeriv, groupVelVectorsForDir) live in
# ChristoffelDerivs.jl; the relocated `ElastConsts.CalcAnisotropyRatio` lives in
# AnisotropyRatio.jl; `LimitingVels.m` lives in LimitingVels.jl. All three are `include`d at
# the bottom of this file, directly into THIS module -- one module per top-level folder
# (`christoffel/`), not one per file (see CoreTypes.jl's header for the "no per-file
# submodules" rationale). Their functions therefore need no `Christoffel.` prefix to call
# each other, exactly like MatlabCompat.jl/Legendre.jl/etc. inside `module CoreTypes`.
#
# `Christoffel` genuinely IS a namespace with no type of the same name -- `classdef
# Christoffel` in MATLAB is `methods (Static)` only (Christoffel.m), never instantiated.
# So `module Christoffel` here does NOT collide with anything (there is no `struct
# Christoffel` to fight over the name), unlike `ElastConsts.jl`/`Material.jl`, where the
# file's main export shares the file's name and therefore CANNOT be wrapped in a same-named
# module.
#
# MODULE NESTING / CROSS-MODULE DEPENDENCY -- load-bearing assumption: this file needs
# `ElastConsts`, `Material` and several of their free functions from `CoreTypes`, which is a
# SIBLING of `Solvers` (this module's parent), not an ancestor reachable by plain lexical
# nesting. The nesting is exactly
#     module GuliwerrJ
#       include("core/CoreTypes.jl");  using .CoreTypes      # CoreTypes is a DIRECT child
#       include("solvers/Solvers.jl"); using .Solvers        # Solvers is a DIRECT child
#         # Solvers.jl: module Solvers; include("christoffel/Christoffel.jl"); using .Christoffel; end
#     end
# i.e. `Christoffel` sits exactly TWO levels below `GuliwerrJ` (GuliwerrJ -> Solvers ->
# Christoffel). Julia's relative-import dot-count climbs the LEXICAL module-nesting tree,
# not the `using` graph, so from inside `module Christoffel`, `..` reaches `Solvers` (its
# direct parent) and `...` reaches `GuliwerrJ` (its grandparent) -- and `CoreTypes` is bound
# directly in `GuliwerrJ`'s own namespace (via that first `include`), so `using ...CoreTypes`
# below resolves correctly REGARDLESS of what `Solvers.jl` itself does or does not `using`.
# This works precisely BECAUSE Julia submodules do not automatically inherit a parent
# module's own `using`-imported names -- each module must re-import what it needs itself --
# so relying on `Solvers.jl` to have already brought `CoreTypes` into scope would be fragile;
# climbing to the common ancestor by dot-count is the robust idiom instead. **If Solvers.jl
# wraps Christoffel at a different nesting depth than exactly two levels below GuliwerrJ,
# this `using ...CoreTypes` line (and its twin in `ritz/GammaMatrix.jl`) needs its dot-count
# adjusted to match.**
# =========================================================================================

module Christoffel

using LinearAlgebra

using ...CoreTypes: ElastConsts, Material, SymmClassLibrary, Any2Tensor, Tensor2Vector,
                     GetRotDirMat, GetStandardDirMat, SetDirMat, GetRotElastTensor,
                     fminbnd_brent

export matrix, GammaFromKron, phaseVels, phaseVelScan

# -----------------------------------------------------------------------------------------
# Small shared helper (used by this file and by ChristoffelDerivs.jl / LimitingVels.jl,
# `include`d below into this same module): MATLAB's `nargin<2||isempty(angleDeg)` /
# `isscalar(angleDeg)` idiom for "one angle or a vector of angles, default 0" recurs across
# phaseVelScan.
# -----------------------------------------------------------------------------------------

"""Coerce a MATLAB-style `angleDeg` argument (scalar OR vector, degrees) to `Vector{Float64}`."""
_as_angle_vec(x::Real) = [Float64(x)]
_as_angle_vec(x) = Float64.(collect(x))

# -----------------------------------------------------------------------------------------
# GammaFromKron.m -- MATLAB `methods (Static, Access = private)` inside Christoffel.m
# (Christoffel.m). Julia has no access control; kept unexported only informally
# (not listed in the module `export` line above the way `matrix`/`phaseVels` etc. are) --
# it is reachable as `Christoffel.GammaFromKron` regardless of export, same as MATLAB's private
# static method is reachable only via `Christoffel.GammaFromKron` from inside the classdef.
# -----------------------------------------------------------------------------------------

"""
    GammaFromKron(elTensor, nnKron, rho=1.0) -> Matrix{Float64}

Christoffel-style contraction `Gamma(i,k) = C(i,j,k,l)*nn(j,l) / rho`.

`nnKron` is `kron(n,n)` (9-vector) for the ordinary Christoffel matrix, or its total
differential `kron(dn,n) + kron(n,dn)` for the angle derivative
([`phaseVelAngleDeriv`](@ref)). The `[1,3,2,4]` permutation is what pairs `(j,l)` with
`kron`'s index order -- see `Christoffel.m`'s own comment.

Port of `guliwerr/solvers/christoffel/@Christoffel/Christoffel.m` (`GammaFromKron`,
a private static method of the classdef itself, not a separate file).
"""
function GammaFromKron(elTensor::AbstractArray{<:Real,4}, nnKron::AbstractVector{<:Real},
                        rho::Real = 1.0)
    M9 = reshape(permutedims(elTensor, (1, 3, 2, 4)), 9, 9)
    return reshape(M9 * nnKron, 3, 3) ./ rho
end

# -----------------------------------------------------------------------------------------
# SortedEig -- MATLAB `methods (Static, Access = private)` (Christoffel.m)
# -----------------------------------------------------------------------------------------

"""
    SortedEig(M) -> (eigVals, eigVecs)

Eigendecomposition of a real 3x3 Christoffel matrix, sorted DESCENDING by eigenvalue --
the `[L, T1, T2]` wave-type convention used throughout this module.

`M` is symmetrised first, `(M + M')/2`, exactly as MATLAB's `SortedEig` does: a Christoffel
matrix is symmetric physically, but one built from a ROTATED elastic tensor is not, at
round-off level, and MATLAB's `eig` then takes its general solver (an exactly degenerate isotropic shear pair came back as a complex pair, fused silica at
150 deg). Julia's `Symmetric` wrap alone would read only the upper triangle, which differs
from the symmetrised matrix at round-off; symmetrising first keeps the input bit-identical
to MATLAB's, and the wrap still gives a concrete (type-stable) `Eigen{Float64}` via the
symmetric LAPACK solver both languages use.

Eigenvector SIGN is not cross-language reproducible (see [`matrix`](@ref)); compare
eigenvectors as rays (up to sign) or, for a degenerate pair, as subspaces.

Port of `guliwerr/solvers/christoffel/@Christoffel/Christoffel.m` (`SortedEig`).
"""
function SortedEig(M::AbstractMatrix{<:Real})
    F = eigen(Symmetric((M .+ M') ./ 2))
    ord = sortperm(F.values; rev = true)
    return F.values[ord], F.vectors[:, ord]
end

# -----------------------------------------------------------------------------------------
# matrix.m
# -----------------------------------------------------------------------------------------

"""
    matrix(elTensor, dirVec, rho=1.0, normalizeDir=true) -> (eigVals, eigVecs, GammaMat)

Christoffel matrix and its eigendecomposition.

INPUT: `elTensor` -- 3x3x3x3 elastic tensor; `dirVec` -- propagation direction (3-vector);
`rho` -- density (default `1`, i.e. pass `1` for density-normalised input); `normalizeDir`
-- normalise `dirVec` before use (default `true`).
OUTPUT: `eigVals` -- eigenvalues sorted descending (`= rho*v_P^2` for L, T1, T2); `eigVecs`
-- corresponding eigenvectors as columns (polarisations); `GammaMat` -- 3x3 Christoffel
matrix.

**Degenerate case:** for an isotropic material the two shear eigenvalues coincide and
`eigVecs`' corresponding two columns span an arbitrary orthonormal basis of that 2-D
eigenspace -- driver-dependent between MATLAB and Julia even though both call the same
LAPACK family. Compare those two columns as a subspace, never element-wise or as individual
rays.

**Non-degenerate case -- sign is ALSO not cross-language reproducible.** One might expect
that away from a degeneracy both languages' `syev`-family solvers agree on the eigenvector
sign, but they do not: on the Cu(001) test fixture (19 angles, all three eigenvalues well
separated, no near-degeneracy), 15 of 19 angles have at least one eigenvector column exactly
sign-flipped between MATLAB and Julia, even though the eigenvalues themselves agree to
`rtol = 1e-12` and each flipped column is still the correct ray (agreement up to sign at
`rtol = 1e-10` throughout). So the eigensolvers carry no portable sign convention in general.
**Always compare `matrix`'s `eigVecs` (and anything built from them without an intervening
`abs(...)`) as rays (up to sign), never element-wise, even away from any eigenvalue
degeneracy.**

Port of `guliwerr/solvers/christoffel/@Christoffel/matrix.m`.
"""
function matrix(elTensor::AbstractArray{<:Real,4}, dirVec::AbstractVector{<:Real},
                 rho::Real = 1.0, normalizeDir::Bool = true)
    size(elTensor) == (3, 3, 3, 3) ||
        throw(ArgumentError("Christoffel:matrix:badInput -- elTensor must be a 3x3x3x3 tensor"))
    length(dirVec) == 3 ||
        throw(ArgumentError("Christoffel:matrix:badInput -- dirVec must have 3 elements"))

    n = Float64.(vec(dirVec))
    normalizeDir && (n = n ./ norm(n))

    # Not the Ritz GammaMatrix class (guliwerr/solvers/ritz/@GammaMatrix) -- same name,
    # different physics object (Christoffel.m's own comment).
    GammaMat = GammaFromKron(elTensor, kron(n, n), rho)

    eigVals, eigVecs = SortedEig(GammaMat)
    return eigVals, eigVecs, GammaMat
end

# -----------------------------------------------------------------------------------------
# phaseVels.m
# -----------------------------------------------------------------------------------------

"""
    phaseVels(elTensor, dirVec, rho=1.0) -> Vector{Float64}

Phase velocities `[vL; vT1; vT2]` for a given tensor and direction, sorted descending.

Port of `guliwerr/solvers/christoffel/@Christoffel/phaseVels.m`.

**Divergence from MATLAB, deliberately not reproduced:** MATLAB's `sqrt` of a negative
real silently returns a complex value; Julia's `sqrt` throws `DomainError`. Every eigenvalue
of a Christoffel matrix built from a physically stable (positive semi-definite) elastic
tensor is `>= 0`, so this never fires for a physically valid material -- recorded here
rather than papered over with a `max(..., 0)` clamp that would silently hide a genuinely
unstable elastic tensor instead of erroring on it.
"""
function phaseVels(elTensor::AbstractArray{<:Real,4}, dirVec::AbstractVector{<:Real},
                    rho::Real = 1.0)
    eigVals, _, _ = matrix(elTensor, dirVec, rho)
    return sqrt.(eigVals)
end

# -----------------------------------------------------------------------------------------
# phaseVelScan.m
# -----------------------------------------------------------------------------------------

"""
    phaseVelScan(mater::Material, angleDeg=0.0) -> Matrix{Float64}

Phase velocities scanned over rotation angles for a [`Material`](@ref). `vels(wavetype,
angle)` -- `[vL; vT1; vT2]` at each angle.

Port of `guliwerr/solvers/christoffel/@Christoffel/phaseVelScan.m`. MATLAB's
`~isa(mater,'Material')` guard (`Christoffel:phaseVelScan:badInput`) becomes ordinary Julia
dispatch on the `mater::Material` annotation -- a non-`Material` first argument is a
`MethodError` at the call site rather than a runtime `error(...)`.
"""
function phaseVelScan(mater::Material, angleDeg = 0.0)
    angs = _as_angle_vec(angleDeg)
    vels = zeros(Float64, 3, length(angs))
    for (ia, a) in pairs(angs)
        rotDirMat = GetRotDirMat(mater, a)
        dirVec = rotDirMat[2, :]    # propagation direction for this angle
        vels[:, ia] = phaseVels(mater.elConsts.tensor, dirVec, mater.rho)
    end
    return vels
end

# -----------------------------------------------------------------------------------------
# The other three files in this folder -- plain top-level definitions `include`d directly
# into THIS module, per the "one module per top-level folder" convention (see this file's
# header). Order matters here exactly as it does in CoreTypes.jl: ChristoffelDerivs.jl uses
# `matrix`/`SortedEig` defined above; AnisotropyRatio.jl uses `phaseVels` defined above;
# LimitingVels.jl uses `_kronperm`/`_pvad_perm` (ChristoffelDerivs.jl) and `fminbnd_brent`
# (`using` above).
# -----------------------------------------------------------------------------------------

include("ChristoffelDerivs.jl")
include("AnisotropyRatio.jl")
include("LimitingVels.jl")

end # module Christoffel
