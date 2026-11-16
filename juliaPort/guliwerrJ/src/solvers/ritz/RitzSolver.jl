# =========================================================================================
# RitzSolver.jl -- port of guliwerr/solvers/ritz/RitzSolver.m
#
# Base of the one real GULIWERR inheritance chain (`RitzSolverSAW < RitzSolver`,
# `RitzSolverLamb < RitzSolver`). Shared eigen-solve, real gauge, group velocity, constraint
# matrix, displacement, angular slope.
#
# SHAPE: `RitzSolver` is an abstract type; each concrete subtype is an IMMUTABLE struct declaring the shared fields
# itself, plus a `_rebuild(s; field = value, ...)` copy-with-replacements method that the
# shared functions below use where MATLAB writes `obj.field = ...; return obj`. No per-file
# `module` (a `module RitzSolver` would collide with the abstract type).
#
# THE EIGEN PATH (as in MATLAB):
#
#  1. `solveConstrainedEig` symmetrises K_con (and M_con) and solves on the Hermitian /
#     Cholesky path -- `eigen(Hermitian(K))` / `eigen(Hermitian(K), Hermitian(M))` here,
#     `eig` on an exactly Hermitian input in MATLAB. Eigenvalues come back real and sorted.
#     There is no need to avoid the `Hermitian` wrap to preserve a general-eigensolver
#     (zgeev) phase convention: MATLAB does not rely on one either.
#  2. `realGauge`: when every layer is invariant under a 2-fold rotation about the surface
#     normal, D^H K D with D = diag(1, 1, i) per Legendre coefficient is REAL symmetric, and
#     the solve runs on the real symmetric path (2-3x faster). The eigenvectors are
#     mapped back, alpha = N^T (d .* (N beta)), so `eigVecs_con` keeps its meaning.
#  3. `phaseFix` (in `calcDisplacement`) makes each mode's largest in-plane Legendre
#     coefficient real and positive. The full-space vector N * alpha is independent of which
#     orthonormal basis N of the null space was returned (N N^T is the projector onto that
#     space), so after phaseFix `calcDisplacement` is reproducible ACROSS implementations
#     (MATLAB `null` vs Julia `nullspace`) up to round-off, for every non-degenerate mode --
#     signed values agree element-wise with MATLAB. `eigVecs_con` itself stays
#     basis-dependent (MATLAB keeps the solver's phase there too); only quantities built from
#     N * eigVecs_con, or phase-invariant products in one basis, are comparable between the
#     languages.
#
# EFFICIENCY (the Ritz solve is the hot path of every sweep and inversion):
#  - K is block-diagonal by layer: the gauge test, the real-gauge matrix and the projection
#    N^T K N are done block by block, with real BLAS (a complex K is split into its real and
#    imaginary parts, so no complex-by-real generic matmul is ever hit).
#  - Group velocity, dVel_dDepth and depthInsensitivity are Rayleigh quotients
#    alpha^H (N^T dK N) alpha = (N alpha)^H dK (N alpha). MATLAB projects dK/dD
#    (strBVL^2 * nCon) and then forms the full nCon x nCon product V' * dK_con * V only to
#    take its diagonal; here the quotients are taken per layer block on the full-space
#    eigenvectors (bVL^2 * nModes per block). Same numbers to round-off, a fraction of the
#    work, and the projected dK/dD is never formed.
#  - A single-layer structure has M = I exactly (normFac = 1), so M_con = N^T N = I to
#    round-off and MATLAB's identity test always selects the standard solve; the projection
#    of M is skipped and the standard path taken directly.
#  - `modeDiagnostics` evaluates only the top-surface u3 row of the displacement (all it
#    reads), not the whole depth profile.
# =========================================================================================

using LinearAlgebra

using ...CoreTypes: Structure, CalcParams, Material, nlay, str_basis_vec_l,
                     basis_decomposition_2d

export RitzSolver, solveConstrainedEig, realGauge, gaugeVector, isGaugeReal, phaseFix,
       REAL_GAUGE_TOL, precomp_field,
       buildConstraintMatrix, solveEigenvalueProblem,
       calcGroupVelocity, calcDisplacement, angleSlope, pruneEigVecs, ritzStructurePrecompute,
       ritzGramBlocks, gaugeBackMap

# -----------------------------------------------------------------------------------------
# The abstract type
# -----------------------------------------------------------------------------------------

"""
    RitzSolver

Abstract supertype of GULIWERR's Ritz-Rayleigh forward solvers:
[`RitzSolverSAW`](@ref) (half space / layered half space) and [`RitzSolverLamb`](@ref)
(free plate). Port of `guliwerr/solvers/ritz/RitzSolver.m`.

Every subtype declares these fields (MATLAB's `RitzSolver` properties):

| field | type | meaning |
|---|---|---|
| `structure` | `Structure` | computational domain |
| `angle` | `Float64` | propagation angle [deg] w.r.t. `dirMat` row 2 |
| `nCompModes` | `Int` | number of stored modes |
| `freqs` | `Vector{Float64}` | `sqrt(eigVals)/(2pi)` |
| `vels` | `Vector{Float64}` | phase velocities (`= freqs`, `relLength = 1`) |
| `velsGroup` | `Vector{Float64}` | group velocities `vp + sum_i d_i dvp/dd_i` |
| `eigVecs_con` | `Matrix{ComplexF64}` | eigenvectors in the constrained basis, M-orthonormal |
| `eigVals` | `Vector{Float64}` | eigenvalues `omega^2` (normalised), ascending |
| `nullConstrMat` | `Matrix{Float64}` | `N`, orthonormal basis of the constraint null space |
| `basisVecL_con` | `Int` | rows of `eigVecs_con` (columns of `N`) |
| `dVel_dDepth` | `Vector{Float64}` | `dvp/dD` (every layer `+1` at once for `nLay > 1`) |
| `depthInsensitivity` | `Vector{Float64}` | `log10(1/abs(dvp/dD))` (substrate-only for layered SAW) |

plus a `_rebuild(s; field = s.field, ...)` keyword method covering every field.

The assembled `GammaMatrix` (`K`, `dK/dD`, `M`) is **not** a field: `calc` keeps it as a local
and hands it to [`_project_for_solve`](@ref) and [`calcGroupVelocity`](@ref). Nothing reads it
after the solve, and stored per angle it was most of a sweep's memory (~3 MB per angle at
strBVL 240, more than the full eigenvector set). Rebuild it with `calc(s.structure,
SAWBoundary(), s.angle)` (`LambBoundary()` for a plate), ~1 ms, as MATLAB does.
"""
abstract type RitzSolver end

"""
    _rebuild(s::T; <field> = s.<field>, ...) -> T

Copy of `s` with the named fields replaced (immutable-struct analogue of MATLAB's
`obj.field = value` inside a value-class method). One method per concrete subtype.
"""
function _rebuild end

"MATLAB `RitzSolver.REAL_GAUGE_TOL`: `||Im(D'KD)||_F / ||D'KD||_F` must be at round-off level."
const REAL_GAUGE_TOL = 1e-13

# -----------------------------------------------------------------------------------------
# precompField -- MATLAB `methods (Static, Access = protected)`
# -----------------------------------------------------------------------------------------

"""
    precomp_field(precomp::NamedTuple, field::Symbol, default = nothing)

Read one field of a caller-supplied precompute bundle, or compute a fallback lazily.
Port of `RitzSolver.precompField`. `default` is a zero-argument closure, called only when
`field` is absent, `nothing` or an empty array -- the fallbacks are expensive (a null-space
solve, a reference-plate solve) and must not run when the caller supplied the value.
"""
function precomp_field(precomp::NamedTuple, field::Symbol, default = nothing)
    if haskey(precomp, field)
        val = getfield(precomp, field)
        if val !== nothing && !(val isa AbstractArray && isempty(val))
            return val
        end
    end
    return default === nothing ? nothing : default()
end

# -----------------------------------------------------------------------------------------
# Layer bookkeeping
# -----------------------------------------------------------------------------------------

"Row range of each layer's block in the full (unconstrained) basis."
function _layer_ranges(calcParams::AbstractVector{CalcParams})
    out = Vector{UnitRange{Int}}(undef, length(calcParams))
    r0 = 0
    for (il, cp) in pairs(calcParams)
        out[il] = (r0 + 1):(r0 + cp.basisVecL)
        r0 += cp.basisVecL
    end
    return out
end

"True when `N` is the identity returned for an unconstrained (single-layer Lamb) problem."
_is_unconstrained(N::AbstractMatrix) = size(N, 1) == size(N, 2)

# -----------------------------------------------------------------------------------------
# buildConstraintMatrix -- MATLAB `methods (Static)`
# -----------------------------------------------------------------------------------------

"""
    buildConstraintMatrix(calcParams, ::SAWBoundary) -> Matrix{Float64}
    buildConstraintMatrix(calcParams, ::LambBoundary) -> Matrix{Float64}

Orthonormal basis `N` of the null space of the boundary-condition constraint matrix,
`strBasisVecL x (strBasisVecL - nConstraints)`. Port of `RitzSolver.buildConstraintMatrix`
(`bcType = 'saw' | 'lamb'` becomes dispatch on `SAWBoundary()`/`LambBoundary()`).

- `SAWBoundary`: rigid bottom (`u_i = 0` at the bottom of layer 1) plus inter-layer
  displacement continuity.
- `LambBoundary`: inter-layer continuity only; a single layer has no constraints and `N`
  is exactly `I`.

The basis is NOT the one MATLAB's `null` returns (the null singular values are exactly
zero, so LAPACK's basis of that subspace is arbitrary), but it spans the same subspace. No
output depends on the choice except `eigVecs_con` itself: eigenvalues, Rayleigh quotients,
`N * eigVecs_con` (hence `calcDisplacement`) and Hermitian overlaps within one basis are all
basis-independent.
"""
buildConstraintMatrix(calcParams::CalcParams, bc::Union{SAWBoundary,LambBoundary}) =
    buildConstraintMatrix([calcParams], bc)
buildConstraintMatrix(calcParams::AbstractVector{CalcParams}, ::SAWBoundary) =
    _build_constraint_matrix(calcParams, true)
buildConstraintMatrix(calcParams::AbstractVector{CalcParams}, ::LambBoundary) =
    _build_constraint_matrix(calcParams, false)

function _build_constraint_matrix(calcParams::AbstractVector{CalcParams}, rigidBottom::Bool)
    noLays = length(calcParams)
    legO   = [cp.legOrder  for cp in calcParams]
    bVL    = [cp.basisVecL for cp in calcParams]
    strBVL = sum(bVL)

    constrMat = zeros(Float64, strBVL, strBVL)   # zero rows are ignored by the null space
    constrNo  = 1

    if rigidBottom      # u_i = 0 at the bottom (Legendre argument -1) of layer 1
        for iiD in 1:3
            for iiL in 1:bVL[1]
                iL = mod(iiL - 1, legO[1])
                iD = cld(iiL, legO[1])
                constrMat[constrNo, iiL] = (iD == iiD) * ((-1)^iL) * sqrt((2 * iL + 1) / 2)
            end
            constrNo += 1
        end
    end

    if noLays > 1       # inter-layer displacement continuity
        for ilay in 1:(noLays - 1)
            legIntLow  = ilay == 1 ? (1:bVL[ilay]) :
                         ((sum(bVL[1:(ilay - 1)]) + 1):sum(bVL[1:ilay]))
            legIntHigh = (sum(bVL[1:ilay]) + 1):sum(bVL[1:(ilay + 1)])
            for iiD in 1:3
                for iiL in legIntLow          # top of the lower layer: P_m(+1)
                    k  = iiL - first(legIntLow) + 1
                    iL = mod(k - 1, legO[ilay])
                    iD = cld(k, legO[ilay])
                    constrMat[constrNo, iiL] = (iD == iiD) * sqrt((2 * iL + 1) / 2)
                end
                for iiL in legIntHigh         # bottom of the upper layer: P_m(-1)
                    k  = iiL - last(legIntLow)
                    iL = mod(k - 1, legO[ilay + 1])
                    iD = cld(k, legO[ilay + 1])
                    constrMat[constrNo, iiL] = -(iD == iiD) * ((-1)^iL) * sqrt((2 * iL + 1) / 2)
                end
                constrNo += 1
            end
        end
    end

    constrNo == 1 && return Matrix{Float64}(I, strBVL, strBVL)   # single-layer Lamb
    return nullspace(constrMat)
end

# -----------------------------------------------------------------------------------------
# Real gauge -- MATLAB `methods (Static)`: realGauge, gaugeVector, isGaugeReal
# -----------------------------------------------------------------------------------------

"""
    gaugeVector(calcParams) -> Vector{ComplexF64}

Diagonal `d` of the gauge `D = diag(1, 1, i)` per Legendre coefficient: `1` on every layer's
`u1`, `u2` coefficients, `im` on its `u3` coefficients. Port of `RitzSolver.gaugeVector`.
"""
gaugeVector(cp::CalcParams) = gaugeVector([cp])
function gaugeVector(calcParams::AbstractVector{CalcParams})
    d = ComplexF64[]
    for cp in calcParams
        L = cp.legOrder
        append!(d, ones(ComplexF64, 2L))
        append!(d, fill(ComplexF64(im), L))
    end
    return d
end

"""
    isGaugeReal(Kg, scale = norm(Kg)) -> Bool

Is a gauged matrix real to round-off: `norm(imag(Kg)) <= REAL_GAUGE_TOL * scale`
(Frobenius norms). Port of `RitzSolver.isGaugeReal`. Works on any array (the inversion
passes per-parameter derivative stacks).
"""
function isGaugeReal(Kg::AbstractArray, scale::Real = norm(Kg))
    s2 = 0.0
    @inbounds for x in Kg
        s2 += abs2(imag(x))
    end
    return sqrt(s2) <= REAL_GAUGE_TOL * scale
end

"""
    realGauge(K, calcParams) -> (Kout, d)

Real-symmetric gauge `u3 -> i u3` of the full stiffness matrix, if valid. Port of
`RitzSolver.realGauge`.

If every layer's elastic tensor is invariant under a 2-fold rotation about the surface
normal (isotropic, cubic (001)/(110), hexagonal, tetragonal/orthorhombic with a principal
axis normal, ...), the imaginary part of `K` sits only in the in-plane/u3 blocks and
`D^H K D` is real. Returns `(real(D^H K D)::Matrix{Float64}, d::Vector{ComplexF64})` when
the test [`isGaugeReal`](@ref) passes, else `(K, nothing)` (MATLAB: `d = []`). The test is
at round-off level on purpose: for cubic (111) the residual is ~3e-5 of `norm(K)` and
forcing it real shifts eigenvalues by 9-19 %.

Elementwise, `conj(d_i) K_ij d_j` multiplies by exactly `1` or `+-im`, so the gauged entries
are bit-identical to MATLAB's `conj(d) .* K .* d.'`.
"""
function realGauge(K::AbstractMatrix, calcParams)
    cps = calcParams isa CalcParams ? [calcParams] : calcParams
    Kr = _real_gauge_blocks(K, _layer_ranges(cps), cps)
    Kr === nothing && return K, nothing
    return Kr, gaugeVector(cps)
end

"""
    _real_gauge_blocks(K, ranges, calcParams) -> Union{Matrix{Float64}, Nothing}

`real(D^H K D)` if the gauge test passes, else `nothing`. `K` must be block-diagonal by
layer (as every `GammaMatrix` is); only the diagonal blocks are read, and the off-block
zeros contribute nothing to either norm. One pass, no complex temporary.
"""
function _real_gauge_blocks(K::AbstractMatrix{<:Complex}, ranges, calcParams)
    n  = size(K, 1)
    Kr = zeros(Float64, n, n)
    im2 = 0.0
    all2 = 0.0
    @inbounds for (il, r) in pairs(ranges)
        L  = calcParams[il].legOrder
        r0 = first(r) - 1
        nb = length(r)
        for jj in 1:nb
            jU3 = jj > 2L
            for ii in 1:nb
                iU3 = ii > 2L
                z = K[r0 + ii, r0 + jj]
                # conj(d_i) * z * d_j with d in {1, im}
                g = iU3 == jU3 ? z : (iU3 ? complex(imag(z), -real(z)) :   # -im * z
                                            complex(-imag(z), real(z)))    #  im * z
                Kr[r0 + ii, r0 + jj] = real(g)
                im2  += abs2(imag(g))
                all2 += abs2(g)
            end
        end
    end
    return sqrt(im2) <= REAL_GAUGE_TOL * sqrt(all2) ? Kr : nothing
end
_real_gauge_blocks(K::AbstractMatrix{<:Real}, ranges, calcParams) = Matrix{Float64}(K)

# -----------------------------------------------------------------------------------------
# Projection into the constrained basis
# -----------------------------------------------------------------------------------------

"""
    _sandwich(N, K, ranges) -> N^T K N

`N^T K N` for a block-diagonal `K`, summed block by block with real BLAS. A complex `K` is
split into real and imaginary parts (`N` is real), so `N^T K N` is two real sandwiches --
the plain transpose of MATLAB's `N.'`, never a conjugation of `K`.
"""
function _sandwich(N::AbstractMatrix{Float64}, K::AbstractMatrix{Float64}, ranges)
    m = size(N, 2)
    out = zeros(Float64, m, m)
    for r in ranges
        whole = r == axes(N, 1)          # one layer: the block is the whole matrix, no copies
        Nb = whole ? N : N[r, :]
        KN = (whole ? K : K[r, r]) * Nb
        mul!(out, transpose(Nb), KN, 1.0, 1.0)
    end
    return out
end
function _sandwich(N::AbstractMatrix{Float64}, K::AbstractMatrix{ComplexF64}, ranges)
    m  = size(N, 2)
    re = zeros(Float64, m, m)
    ii = zeros(Float64, m, m)
    for r in ranges
        Nb = r == axes(N, 1) ? N : N[r, :]
        Kb = view(K, r, r)
        mul!(re, transpose(Nb), real.(Kb) * Nb, 1.0, 1.0)
        mul!(ii, transpose(Nb), imag.(Kb) * Nb, 1.0, 1.0)
    end
    return complex.(re, ii)
end

"""
    _project_for_solve(s::RitzSolver, gm::GammaMatrix, sp = nothing) -> (K_con, M_con, gauge)

The hot-path projection used by the subclasses' `calc`: `K_con` in the real gauge when it
applies (`Matrix{Float64}`, `gauge = d`), else complex (`gauge = nothing`); `M_con` is
`nothing` for a single layer (M = I exactly, see the file header), else `N^T M N` --
taken from `sp` ([`ritzStructurePrecompute`](@ref)) when it matches.
"""
function _project_for_solve(s::RitzSolver, gm::GammaMatrix, sp = nothing)
    cps    = s.structure.calcParams
    ranges = _layer_ranges(cps)
    N      = s.nullConstrMat
    K      = gm.gammaMatAll
    Kr     = _real_gauge_blocks(K, ranges, cps)
    unc    = _is_unconstrained(N)
    if Kr === nothing
        K_con = unc ? Matrix{ComplexF64}(K) : _sandwich(N, K, ranges)
        gauge = nothing
    else
        K_con = unc ? Kr : _sandwich(N, Kr, ranges)
        gauge = gaugeVector(cps)
    end
    M_con = if length(cps) == 1
        nothing
    elseif _structure_cache_valid(sp, N, gm.normMatAll, ranges)
        sp.M_con
    else
        _project_mass(N, gm.normMatAll, ranges)
    end
    return K_con, M_con, gauge
end

# -----------------------------------------------------------------------------------------
# Per-structure (angle-independent) part of the projection -- cached by sweeps
# -----------------------------------------------------------------------------------------

"Per-layer `normFac` read off a block-diagonal `M` (one scaled identity per layer block)."
_block_norm_facs(M::AbstractMatrix, ranges) = Float64[M[first(r), first(r)] for r in ranges]

"""
    _project_mass(N, M, ranges) -> Matrix{Float64}

`M_con = N^T M N` for the Ritz mass matrix, whose layer blocks are scaled identities
(`normFac * I`): `sum_b normFac_b * N_b^T N_b`, one real BLAS product per block instead of
the dense sandwich with the (diagonal) block. Same matrix to round-off.
"""
function _project_mass(N::AbstractMatrix{Float64}, M::AbstractMatrix, ranges)
    _is_unconstrained(N) && return Matrix{Float64}(M)
    m   = size(N, 2)
    out = zeros(Float64, m, m)
    for r in ranges
        Nb = N[r, :]
        mul!(out, transpose(Nb), Nb, M[first(r), first(r)], 1.0)
    end
    return out
end

"""
    ritzStructurePrecompute(structure, nullConstrMat) -> NamedTuple

The angle-independent pieces of a Ritz solve that a sweep can share between its angles
(they depend only on `calcParams`, the layer depths and densities, and `N`), returned as the
precompute fields `ritzStruct = (; N, normFacs, M_con, gaugeBack)`:

- `M_con = N^T M N` (layered structures only; `nothing` for one layer, where `M = I`), with
  the per-layer `normFacs` it was built from;
- `gaugeBack` -- the real-gauge back-map `N^T (d .* (N B))` in compact form
  ([`_gauge_back_map`](@ref)), or `nothing` when it does not apply.

[`calc`](@ref) uses them only when `N` is the very matrix they were built from (`===`) and,
for `M_con`, when the solve's `M` has the same `normFacs`; otherwise it recomputes. Results
equal the uncached path to round-off (the products are associated differently).

`gaugeBack` and `gram` (keywords) take the parts that depend on `N` alone from a wavelength
sweep ([`ritzGramBlocks`](@ref), `ritzWavelengthSweep`): the back-map is reused as given and
`M_con` is the weighted sum `sum_b normFac_b * gram[b]` (round-off-equal to the products).
"""
function ritzStructurePrecompute(structure::Structure, N::AbstractMatrix{Float64};
                                 gaugeBack = :compute, gram = nothing)
    cps    = structure.calcParams
    ranges = _layer_ranges(cps)
    gb = gaugeBack === :compute ? _gauge_back_map(N, gaugeVector(cps)) : gaugeBack
    if gram !== nothing                  # wavelength-sweep path, shared with ritzSAWPrecompute!
        m = size(N, 2); layered = length(cps) > 1
        return _structure_from_gram!(zeros(length(cps)),
                                     layered ? Matrix{Float64}(undef, m, m) : nothing,
                                     layered ? Matrix{Float64}(undef, m, m) : nothing,
                                     structure, N, gb, gram)
    end
    normFacs = Float64[(cps[i].relDepth * structure.mater[i].rho) /
                       (cps[1].relDepth * structure.mater[1].rho) for i in eachindex(cps)]
    M_con = nothing
    if length(cps) > 1
        if _is_unconstrained(N)
            M_con = zeros(Float64, size(N, 2), size(N, 2))
            for (il, r) in pairs(ranges), k in r
                M_con[k, k] = normFacs[il]
            end
        else
            M_con = zeros(Float64, size(N, 2), size(N, 2))
            for (il, r) in pairs(ranges)
                Nb = N[r, :]
                mul!(M_con, transpose(Nb), Nb, normFacs[il], 1.0)
            end
        end
    end
    return (; N, normFacs, M_con, gaugeBack = gb, M_conU = nothing)
end

"""
    _structure_from_gram!(normFacs, M_con, M_conU, structure, N, gaugeBack, gram) -> NamedTuple

The wavelength-sweep path of [`ritzStructurePrecompute`](@ref), written into the given arrays
(fresh ones there, a [`RitzWorkspace`](@ref)'s in `ritzSAWPrecompute!`): the per-layer
`normFacs`, `M_con = sum_b normFac_b * gram[b]` and `M_conU`, the upper Cholesky factor of the
symmetrised `M_con` (angle-independent; the lowest-mode solve of `ritzSensitivities` uses it,
as `solveConstrainedEig` would build it). `M_con`/`M_conU` are `nothing` for one layer.
"""
function _structure_from_gram!(normFacs::AbstractVector{Float64}, M_con, M_conU,
                               structure::Structure, N::AbstractMatrix{Float64}, gaugeBack, gram)
    cps = structure.calcParams
    for i in eachindex(cps)
        normFacs[i] = (cps[i].relDepth * structure.mater[i].rho) /
                      (cps[1].relDepth * structure.mater[1].rho)
    end
    length(cps) > 1 ||
        return (; N, normFacs, M_con = nothing, gaugeBack, M_conU = nothing)
    M_con .= normFacs[1] .* gram[1]
    for il in 2:length(cps)
        M_con .+= normFacs[il] .* gram[il]
    end
    M_conU .= (M_con .+ transpose(M_con)) ./ 2
    cholesky!(Hermitian(M_conU))                 # upper triangle <- U (potrf, as cholesky())
    n = size(M_conU, 1)
    @inbounds for j in 1:n, i in (j + 1):n       # as Matrix(F.U): zeros below the diagonal
        M_conU[i, j] = 0.0
    end
    return (; N, normFacs, M_con, gaugeBack, M_conU)
end

"""
    ritzGramBlocks(calcParams, N) -> Vector{Matrix{Float64}}

Per-layer Gram blocks `N_b^T N_b` of the constraint null-space basis (`N_b` = the layer's
rows of `N`), so that `M_con = N^T M N = sum_b normFac_b * N_b^T N_b` for any set of layer
depths and densities. They depend on the Legendre orders alone (through `N`), so a sweep over
wavelengths, where only `normFac_b` changes, builds them once.
"""
function ritzGramBlocks(calcParams::AbstractVector{CalcParams}, N::AbstractMatrix{Float64})
    return [(Nb = N[r, :]; transpose(Nb) * Nb) for r in _layer_ranges(calcParams)]
end

"""
    gaugeBackMap(calcParams, N) -> Union{Nothing, NamedTuple}

The precomputed real-gauge back-map of [`ritzStructurePrecompute`](@ref) (its `gaugeBack`
field) for the constraint basis `N`; depends on `N` and the Legendre orders only.
"""
gaugeBackMap(calcParams::AbstractVector{CalcParams}, N::AbstractMatrix{Float64}) =
    _gauge_back_map(N, gaugeVector(calcParams))

"""
    _gauge_back_map(N, d) -> Union{Nothing, NamedTuple}

Precomputed real-gauge back-map. With `Pi = N^T diag(isU3) N` and `N^T N = I`,
`N^T (d .* (N B)) = (I - Pi) B + im Pi B`. `Pi` is non-zero only on the constrained
coordinates that carry u3 (the constraints act per displacement component, so `N`'s null
space splits by component up to round-off); the map is stored as that index set `idx` and the
dense block `Pi[idx, idx]`, so each angle costs one `|idx|^2 x nModes` product instead of
three `strBVL x nCon x nModes` ones. Returns `nothing` (use the plain back-map) when `N` is the
identity, `N^T N` is not `I` to 1e-12, or `Pi` has entries above 1e-13 outside the block.
"""
function _gauge_back_map(N::AbstractMatrix{Float64}, d::AbstractVector)
    _is_unconstrained(N) && return nothing
    isU3 = imag.(d) .!= 0
    Pi = transpose(N) * (N .* isU3)
    n  = size(Pi, 1)
    norm(transpose(N) * N - I) <= 1e-12 || return nothing
    tol = 1e-13
    idx = [j for j in 1:n if any(i -> abs(Pi[i, j]) > tol, 1:n)]
    inIdx = falses(n); inIdx[idx] .= true
    @inbounds for j in 1:n, i in 1:n
        (inIdx[i] && inIdx[j]) || abs(Pi[i, j]) <= tol || return nothing
    end
    return (; idx, Pib = Pi[idx, idx])
end

"`N^T (d .* (N B))` from a [`_gauge_back_map`](@ref) bundle: `(I - Pi) B + im Pi B`."
function _gauge_back_cached(gb, B::AbstractMatrix{Float64})
    idx = gb.idx
    PB  = gb.Pib * B[idx, :]
    out = Matrix{ComplexF64}(B)
    @inbounds for j in axes(B, 2), (k, i) in pairs(idx)
        out[i, j] = complex(B[i, j] - PB[k, j], PB[k, j])
    end
    return out
end

function _structure_cache_valid(sp, N, M, ranges)
    sp === nothing && return false
    sp.N === N || return false
    sp.M_con === nothing && return false
    return sp.normFacs == _block_norm_facs(M, ranges)
end

# -----------------------------------------------------------------------------------------
# solveConstrainedEig -- MATLAB `methods (Static)`, THE eigen-solve
# -----------------------------------------------------------------------------------------

"""
    solveConstrainedEig(K_con, M_con = nothing; nev = nothing) -> (eigVals, V)

Sorted eigenpairs of `K_con a = omega^2 M_con a` on the Hermitian path. Port of
`RitzSolver.solveConstrainedEig`; the one place the Ritz eigen-solve happens (the solvers'
`solveEigenvalueProblem` and `ElasticityInversion`'s fit loop both call it).

- `K_con` -- `n x n`, Hermitian up to round-off; real (`Matrix{Float64}`, the real gauge)
  or complex. Symmetrised, `(K + K')/2`, then solved as `Hermitian`.
- `M_con` -- `n x n` real symmetric positive definite up to round-off, or `nothing`/empty
  for the identity. A matrix within `1e-8` (max-abs) of `I` is treated as the identity
  (MATLAB's test: `N^T I N` carries ~1e-15 round-off). Otherwise symmetrised and solved by
  the Cholesky path (`eigen(Hermitian(K), Hermitian(M))`); if `M_con` is not numerically
  positive definite, the general QZ solver is used instead (MATLAB falls back the same way)
  and the eigenvectors are renormalised.
- `nev` -- optional: compute only the lowest `nev` eigenpairs (standard path only; the
  generalised path always computes all). Eigenvalues are identical to round-off to the
  first `nev` of the full solve -- used by the solvers when `nMaxModes < n`.
- `lowest` -- with `nev`: compute only the lowest `nev` eigenpairs on the generalised path
  too, by bisection and inverse iteration on the Cholesky-reduced matrix (`dsyevr`/`zheevr`
  range `'I'`) with `abstol = 2*floatmin`. Pays for small `nev` (4-layer, n = 588: 10 ms for
  8 modes against 16 ms for all) and loses against the full solve beyond a few tens of modes.
  The tight `abstol` matters: LAPACK's default (`eps*norm(T)`) left the lowest eigenvalues of
  a stack with a 3 nm film 3e-8 off, against 1e-11 for the full solve and this path.
- `cholM` -- with `lowest`: the upper Cholesky factor of the symmetrised `M_con`
  (`ritzStruct.M_conU` of a wavelength-sweep precompute), so a sweep factorises it once.

Returns `eigVals::Vector{Float64}` ascending and `V` (`Matrix{Float64}` for a real `K_con`,
`Matrix{ComplexF64}` otherwise), M-orthonormal (`V' M V = I`). Each column's phase (sign,
for real `V`) is the solver's choice; see [`phaseFix`](@ref).
"""
function solveConstrainedEig(K_con::AbstractMatrix{T}, M_con = nothing;
                              nev::Union{Integer,Nothing} = nothing,
                              lowest::Bool = false,
                              cholM = nothing) where {T<:Union{Float64,ComplexF64}}
    n  = size(K_con, 1)
    Kh = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        Kh[i, j] = (K_con[i, j] + conj(K_con[j, i])) / 2
    end
    useStandard = M_con === nothing || isempty(M_con) || _near_identity(M_con)
    sub = lowest && nev !== nothing && nev < n
    if useStandard && sub
        λ, V = LinearAlgebra.LAPACK.syevr!('V', 'I', 'U', Kh, 0.0, 0.0, 1, Int(nev),
                                           2 * floatmin(Float64))
    elseif useStandard
        if nev !== nothing && nev < n
            F = eigen!(Hermitian(Kh), 1:Int(nev))
        else
            F = eigen!(Hermitian(Kh))
        end
        λ = F.values::Vector{Float64}
        V = F.vectors::Matrix{T}
    else
        Ms = Matrix{Float64}(undef, n, n)
        @inbounds for j in 1:n, i in 1:n
            Ms[i, j] = real(M_con[i, j] + M_con[j, i]) / 2
        end
        local λ::Vector{Float64}, V::Matrix{T}
        try
            λ, V = sub ? _generalised_hermitian_eig_lowest(Kh, Ms, Int(nev), cholM) :
                         _generalised_hermitian_eig(Kh, Ms)
        catch e
            e isa PosDefException || rethrow()
            F  = eigen(Matrix{ComplexF64}(Kh), Matrix{ComplexF64}(Ms), sortby = nothing)
            λ  = real.(F.values)
            V  = T <: Real ? real.(F.vectors) : Matrix{ComplexF64}(F.vectors)
        end
        # renormalise so a^H M a = 1 (the Cholesky path already does; guards the QZ fallback)
        # (real BLAS on the real/imaginary parts: `Ms * V` real x complex is a generic matmul)
        MV = T <: Real ? Ms * V : complex.(Ms * real.(V), Ms * imag.(V))
        @inbounds for j in axes(V, 2)
            s = 0.0
            for i in axes(V, 1)
                s += real(conj(V[i, j]) * MV[i, j])
            end
            sc = 1 / sqrt(s)
            for i in axes(V, 1)
                V[i, j] *= sc
            end
        end
    end
    if !issorted(λ)
        ord = sortperm(λ)
        λ = λ[ord]
        V = V[:, ord]
    end
    return λ, V
end

"""
    _generalised_hermitian_eig(Kh, Ms) -> (λ, V)

`K v = λ M v` for Hermitian `K` (`Kh`, full storage) and real symmetric positive definite `M`
(`PosDefException` otherwise). Real `K`: LAPACK `dsygvd` (`eigen(Hermitian, Hermitian)`).
Complex `K` (non-gauged stacks): the same Cholesky reduction done explicitly -- `M = U^T U`,
`zhegst` for `C = U^-H K U^-1`, then the MRRR solver `zheevr` (`eigen!(Hermitian)`) and
`V = U^-1 W` with real triangular solves. OpenBLAS's `zhegvd` spends most of its time in the
divide-and-conquer stage; on a 4-layer (111) stack (n = 588, 1 thread) this path is 30 %
faster (467-569 -> 327-388 ms), eigenvalues equal to 5e-12 relative (round-off of the
reduction, cond ~1e6), eigenvectors to 4e-15 in the M inner product. For a real `K` MRRR was
slower than `dsygvd` (50-77 vs 107-115 ms), so that path is unchanged.
"""
function _generalised_hermitian_eig(Kh::Matrix{Float64}, Ms::Matrix{Float64})
    F = eigen(Hermitian(Kh), Hermitian(Ms))
    return F.values, F.vectors
end
function _generalised_hermitian_eig(Kh::Matrix{ComplexF64}, Ms::Matrix{Float64})
    U  = cholesky!(Hermitian(copy(Ms))).U            # real; PosDefException if not SPD
    C  = _zhegst_upper!(Kh, Matrix{ComplexF64}(U))   # Kh is overwritten (caller's scratch)
    F  = eigen!(Hermitian(C, :U))
    W  = F.vectors
    Ut = UpperTriangular(Matrix(U))
    Vr = ldiv!(Ut, real.(W))
    Vi = ldiv!(Ut, imag.(W))
    return F.values, complex.(Vr, Vi)
end

"""
    _generalised_hermitian_eig_lowest(Kh, Ms, k) -> (λ, V)

Lowest `k` eigenpairs of `K v = λ M v`: Cholesky `M = U^T U`, `?hegst`/`?sygst`, then
`?syevr`/`?heevr` over index range `1:k` with `abstol = 2*floatmin` (see
[`solveConstrainedEig`](@ref)'s `lowest`), `V = U^-1 W`.
"""
function _generalised_hermitian_eig_lowest(Kh::Matrix{T}, Ms::Matrix{Float64}, k::Int,
                                           cholU = nothing) where {T}
    Um = cholU === nothing ? Matrix(cholesky!(Hermitian(copy(Ms))).U) : cholU  # PosDefException if not SPD
    C  = T <: Real ? _dsygst_upper!(Kh, Um) : _zhegst_upper!(Kh, Matrix{ComplexF64}(Um))
    λ, W = LinearAlgebra.LAPACK.syevr!('V', 'I', 'U', C, 0.0, 0.0, 1, k, 2 * floatmin(Float64))
    Ut = UpperTriangular(Um)
    T <: Real && return λ, ldiv!(Ut, W)
    return λ, complex.(ldiv!(Ut, real.(W)), ldiv!(Ut, imag.(W)))
end

"LAPACK `dsygst` (itype 1, upper): `A <- U^-T A U^-1` for the Cholesky factor `B = U`."
function _dsygst_upper!(A::Matrix{Float64}, B::Matrix{Float64})
    n = LinearAlgebra.checksquare(A)
    info = Ref{LinearAlgebra.BlasInt}()
    ccall((LinearAlgebra.BLAS.@blasfunc(dsygst_), LinearAlgebra.libblastrampoline), Cvoid,
          (Ref{LinearAlgebra.BlasInt}, Ref{UInt8}, Ref{LinearAlgebra.BlasInt}, Ptr{Float64},
           Ref{LinearAlgebra.BlasInt}, Ptr{Float64}, Ref{LinearAlgebra.BlasInt},
           Ptr{LinearAlgebra.BlasInt}, Clong),
          1, 'U', n, A, max(1, stride(A, 2)), B, max(1, stride(B, 2)), info, 1)
    LinearAlgebra.LAPACK.chklapackerror(info[])
    return A
end

"LAPACK `zhegst` (itype 1, upper): `A <- U^-H A U^-1` for the Cholesky factor `B = U`."
function _zhegst_upper!(A::Matrix{ComplexF64}, B::Matrix{ComplexF64})
    n = LinearAlgebra.checksquare(A)
    info = Ref{LinearAlgebra.BlasInt}()
    ccall((LinearAlgebra.BLAS.@blasfunc(zhegst_), LinearAlgebra.libblastrampoline), Cvoid,
          (Ref{LinearAlgebra.BlasInt}, Ref{UInt8}, Ref{LinearAlgebra.BlasInt}, Ptr{ComplexF64},
           Ref{LinearAlgebra.BlasInt}, Ptr{ComplexF64}, Ref{LinearAlgebra.BlasInt},
           Ptr{LinearAlgebra.BlasInt}, Clong),
          1, 'U', n, A, max(1, stride(A, 2)), B, max(1, stride(B, 2)), info, 1)
    LinearAlgebra.LAPACK.chklapackerror(info[])
    return A
end

function _near_identity(M::AbstractMatrix)
    @inbounds for j in axes(M, 2), i in axes(M, 1)
        abs(M[i, j] - (i == j)) < 1e-8 || return false
    end
    return true
end


# -----------------------------------------------------------------------------------------
# solveEigenvalueProblem -- MATLAB `methods (Access = protected)`
# -----------------------------------------------------------------------------------------

"""
    _gauge_back(N, d, B) -> Matrix{ComplexF64}

Map gauged eigenvectors `B` (constrained basis) back to the original basis,
`N^T (d .* (N B))` -- MATLAB's `N.' * (gauge .* (N * V))`. `range(N)` is invariant under
`D` (each constraint acts on one displacement component), so this recovers the coordinates
exactly. For real `B` the work is three real products (`d` is `1` or `im`).
"""
function _gauge_back(N::AbstractMatrix{Float64}, d::AbstractVector, B::AbstractMatrix{Float64})
    if _is_unconstrained(N)
        return d .* B
    end
    U  = N * B
    isU3 = imag.(d) .!= 0
    Ur = U .* .!isU3
    Ui = U .* isU3
    return complex.(transpose(N) * Ur, transpose(N) * Ui)
end
_gauge_back(N::AbstractMatrix{Float64}, d::AbstractVector, B::AbstractMatrix) =
    _is_unconstrained(N) ? d .* B : Matrix{ComplexF64}(transpose(N) * (d .* (N * B)))

"""
    solveEigenvalueProblem(s::RitzSolver, K_con, M_con, gauge = nothing, sp = nothing) -> typeof(s)

Solve `K_con a = omega^2 M_con a` and return `s` with `eigVals`, `eigVecs_con`, `freqs`,
`vels`, `nCompModes`, `basisVecL_con` filled. Port of `RitzSolver.solveEigenvalueProblem`.
`gauge` is `d` from [`realGauge`](@ref) when `K_con` is the real gauged matrix; the
eigenvectors are then mapped back to the original basis (complex, u1/u2 real, u3
imaginary). Only the lowest `nMaxModes` eigenpairs are computed when that is fewer than
`n` (standard path). `sp` -- optional [`ritzStructurePrecompute`](@ref) bundle: its
`gaugeBack` replaces the per-angle back-map when it was built from `s.nullConstrMat`.
"""
function solveEigenvalueProblem(s::RitzSolver, K_con::AbstractMatrix,
                                 M_con::Union{AbstractMatrix,Nothing} = nothing,
                                 gauge::Union{AbstractVector,Nothing} = nothing,
                                 sp = nothing)
    nmax = s.structure.nMaxModes
    λ, V = solveConstrainedEig(K_con, M_con; nev = nmax)
    nC   = min(nmax, length(λ))
    eigVals = nC == length(λ) ? λ : λ[1:nC]
    Vk = nC == size(V, 2) ? V : V[:, 1:nC]        # V is the solver's own fresh array
    eigVecs_con = if gauge === nothing
        Vk isa Matrix{ComplexF64} ? Vk : Matrix{ComplexF64}(Vk)
    elseif sp !== nothing && sp.N === s.nullConstrMat && sp.gaugeBack !== nothing &&
           Vk isa Matrix{Float64}
        _gauge_back_cached(sp.gaugeBack, Vk)
    else
        _gauge_back(s.nullConstrMat, gauge, Vk)
    end
    # MATLAB's abs(sqrt(eigVals)) for a round-off-negative eigenvalue = sqrt(abs(.))
    freqs = sqrt.(abs.(eigVals)) ./ (2π)
    return _rebuild(s;
        nCompModes    = nC,
        eigVals       = eigVals,
        eigVecs_con   = eigVecs_con,
        freqs         = freqs,
        vels          = copy(freqs),          # v = freq (relLength = 1, |k| = 1)
        basisVecL_con = size(eigVecs_con, 1))
end

# -----------------------------------------------------------------------------------------
# Full-space eigenvectors and per-layer Rayleigh quotients
# -----------------------------------------------------------------------------------------

"""
    _full_eigvecs(N, A) -> Matrix{ComplexF64}

`N * A` for real `N` and complex `A`, with real BLAS (real and imaginary parts separately).
"""
function _full_eigvecs(N::AbstractMatrix{Float64}, A::AbstractMatrix{ComplexF64})
    _is_unconstrained(N) && return Matrix{ComplexF64}(A)
    return complex.(N * real.(A), N * imag.(A))
end

"""
    _block_quotients(U, Kfull, ranges) -> Matrix{Float64}  (nLay x nModes)

`q[il, m] = real(U_b[:, m]^H K_b U_b[:, m])` per layer block `b` of a block-diagonal `K`.
"""
function _block_quotients(U::AbstractMatrix{ComplexF64}, K::AbstractMatrix, ranges)
    nM = size(U, 2)
    q  = zeros(Float64, length(ranges), nM)
    for (il, r) in pairs(ranges)
        whole = r == axes(U, 1)          # one layer: no block copies
        Ub = whole ? U : U[r, :]
        KU = (whole ? K : K[r, r]) * Ub
        @inbounds for m in 1:nM
            acc = 0.0
            for i in axes(Ub, 1)
                acc += real(conj(Ub[i, m]) * KU[i, m])
            end
            q[il, m] = acc
        end
    end
    return q
end

# -----------------------------------------------------------------------------------------
# calcGroupVelocity -- MATLAB `methods (Access = protected)`
# -----------------------------------------------------------------------------------------

"""
    calcGroupVelocity(s::RitzSolver, gm::GammaMatrix, substrateOnlyInsens::Bool = false) -> typeof(s)

Group velocity from the Rayleigh-quotient depth derivatives; fills `dVel_dDepth`,
`velsGroup`, `depthInsensitivity`. Port of `RitzSolver.calcGroupVelocity`:

- `dVel_dDepth = Re(a^H dK/dD a) / v / (8 pi^2)` with `dK/dD` the block-diagonal sum over
  layers (every layer `+1` at once for `nLay > 1`).
- `velsGroup = vp + sum_i d_i dvp/dd_i`: each layer's block weighted by its OWN `relDepth`
  (`vp + relDepth * dVel_dDepth` for one layer). For `nLay > 1`, `vp + D_total *
  dVel_dDepth` would over-weight thin films by `D_total/d_film`; this agrees with a
  finite-difference `d omega/dk`.
- `depthInsensitivity = log10(1/max(|dvp/dD|, eps))`, from the substrate block alone when
  `substrateOnlyInsens` (layered SAW: only the substrate's `relDepth` is an artificial
  truncation depth), else from all blocks. `gm` is the solve's [`GammaMatrix`](@ref) (its
  `dGammaDMatAll` is read).

MATLAB's arguments `(dKdD_con, dKdD_con_depthOnly)` -- projected matrices -- are replaced by
the flag: the quotients are evaluated per block on `N * eigVecs_con`, which is the same
number without projecting `dK/dD` (see the file header).
"""
function calcGroupVelocity(s::RitzSolver, gm::GammaMatrix, substrateOnlyInsens::Bool = false)
    nM     = s.nCompModes
    cps    = s.structure.calcParams
    ranges = _layer_ranges(cps)
    U      = _full_eigvecs(s.nullConstrMat, view(s.eigVecs_con, :, 1:nM))
    q      = _block_quotients(U, gm.dGammaDMatAll, ranges)
    v      = s.vels[1:nM]
    c      = 8 * π^2

    qSum        = vec(sum(q, dims = 1))
    dVel_dDepth = qSum ./ v ./ c
    if length(cps) == 1
        velsGroup = v .+ cps[1].relDepth .* dVel_dDepth
    else
        relDepths = [cp.relDepth for cp in cps]
        velsGroup = v .+ vec(sum(relDepths .* q, dims = 1)) ./ v ./ c
    end
    qIns = substrateOnlyInsens ? q[1, :] : qSum
    depthInsensitivity = log10.(1 ./ max.(abs.(qIns ./ v ./ c), eps(Float64)))

    return _rebuild(s;
        dVel_dDepth        = dVel_dDepth,
        velsGroup          = velsGroup,
        depthInsensitivity = depthInsensitivity)
end

# -----------------------------------------------------------------------------------------
# phaseFix -- MATLAB `methods (Static, Access = protected)`
# -----------------------------------------------------------------------------------------

"""
    phaseFix(eigVecs, calcParams) -> Vector{ComplexF64}

Per-mode unit phase factor making each mode's largest in-plane (u1/u2) Legendre
coefficient real and positive (the largest coefficient overall for a mode with no in-plane
part, `max|in-plane| <= 1e-8 max|all|`). Port of `RitzSolver.phaseFix`. `eigVecs` are FULL
basis eigenvectors (`N * eigVecs_con`), `strBVL x nModes`; `eigVecs .* transpose(ph)` has
the fixed phase. Taken on the coefficients, so it does not depend on the evaluation grid.
Ties in magnitude pick the first row, as MATLAB's `max` does.
"""
phaseFix(eigVecs::AbstractMatrix, calcParams) =
    _phase_fix((i, m) -> eigVecs[i, m], size(eigVecs, 1), size(eigVecs, 2), calcParams)

# phaseFix on element `get(i, m)` of an `nR x nM` eigenvector set (calcDisplacement passes
# `complex(Er[i, m], Ei[i, m])`, the real and imaginary parts held separately)
function _phase_fix(get, nR::Integer, nM::Integer, calcParams)
    cps = calcParams isa CalcParams ? [calcParams] : calcParams
    isInPlane = falses(nR)
    r0 = 0
    for cp in cps
        isInPlane[(r0 + 1):(r0 + 2 * cp.legOrder)] .= true
        r0 += cp.basisVecL
    end
    ph = Vector{ComplexF64}(undef, nM)
    @inbounds for m in 1:nM
        mxIn = -1.0; iIn = 0; mxAll = -1.0; iAll = 0
        for i in 1:nR
            a = abs(get(i, m))
            if a > mxAll; mxAll = a; iAll = i; end
            if isInPlane[i] && a > mxIn; mxIn = a; iIn = i; end
        end
        refRow = mxIn <= 1e-8 * mxAll ? iAll : iIn
        ref = get(refRow, m)
        ph[m] = conj(ref) / max(abs(ref), floatmin(Float64))
    end
    return ph
end

# -----------------------------------------------------------------------------------------
# calcDisplacement -- MATLAB public method
# -----------------------------------------------------------------------------------------

"""
    calcDisplacement(s::RitzSolver, lastModeNo = nothing, dscrtSize = 0.01)
        -> Array{ComplexF64,3}

Depth-profile displacement of the first `lastModeNo` modes (`nothing`: all),
`[depthPoints*nLay x lastModeNo x 3]`, COMPLEX with a fixed phase per mode
([`phaseFix`](@ref)). Port of `RitzSolver.calcDisplacement`.

When every layer has a 2-fold axis along the surface normal (the real gauge applies), the
in-plane components come out real and SIGNED and u3 as `im` x real; otherwise u1, u2, u3
carry genuine relative phases. Callers needing a magnitude apply `abs`. `dscrtSize = 0.01`
reuses each layer's precomputed `polyLeg`.
"""
function calcDisplacement(s::RitzSolver, lastModeNo::Union{Integer,Nothing} = nothing,
                           dscrtSize::Real = 0.01)
    lm = lastModeNo === nothing ? s.nCompModes : Int(lastModeNo)
    displSize = round(Int, 2 / dscrtSize, RoundNearestTiesAway) + 1   # MATLAB round
    cps = s.structure.calcParams
    nL  = length(cps)

    # full-basis eigenvectors N * A as real and imaginary parts (what _full_eigvecs computes,
    # without the complex copy), phase-fixed in place with Base's complex product formula
    N  = s.nullConstrMat
    A  = view(s.eigVecs_con, :, 1:lm)
    Er = _is_unconstrained(N) ? real.(A) : N * real.(A)
    Ei = _is_unconstrained(N) ? imag.(A) : N * imag.(A)
    ph = _phase_fix((i, m) -> complex(Er[i, m], Ei[i, m]), size(Er, 1), lm, cps)
    @inbounds for m in 1:lm
        pr = real(ph[m]); pim = imag(ph[m])
        for i in axes(Er, 1)
            er = Er[i, m]; ei = Ei[i, m]
            Er[i, m] = er * pr - ei * pim
            Ei[i, m] = er * pim + ei * pr
        end
    end

    displDir = zeros(ComplexF64, displSize * nL, lm, 3)
    row0 = 0
    bVL_prev = 0
    for ilay in 1:nL
        cp   = cps[ilay]
        pLeg = dscrtSize == 0.01 ? cp.polyLeg : basis_decomposition_2d(cp.legOrder, dscrtSize)
        nP   = size(pLeg, 2)
        legO = cp.legOrder
        Pr   = Matrix{Float64}(undef, nP, lm)
        Pim  = Matrix{Float64}(undef, nP, lm)
        for iDir in 1:3
            rows = (bVL_prev + (iDir - 1) * legO + 1):(bVL_prev + iDir * legO)
            mul!(Pr, transpose(pLeg), view(Er, rows, :))
            mul!(Pim, transpose(pLeg), view(Ei, rows, :))
            @inbounds for m in 1:lm, i in 1:nP
                displDir[row0 + i, m, iDir] = complex(Pr[i, m], Pim[i, m])
            end
        end
        row0 += nP
        bVL_prev += cp.basisVecL
    end
    return displDir
end

"""
    _top_surface_u3(s::RitzSolver, nM) -> Vector{ComplexF64}

The last row of `calcDisplacement(s, nM)[:, :, 3]` (top surface, x3 = normal) without
evaluating the rest of the profile and without the phase fix (callers take `abs`).
"""
function _top_surface_u3(s::RitzSolver, nM::Integer)
    cps  = s.structure.calcParams
    cp   = cps[end]
    bVL_prev = sum(c.basisVecL for c in cps) - cp.basisVecL
    legO = cp.legOrder
    rows = (bVL_prev + 2 * legO + 1):(bVL_prev + 3 * legO)
    p    = view(cp.polyLeg, :, size(cp.polyLeg, 2))
    N    = s.nullConstrMat
    A    = view(s.eigVecs_con, :, 1:nM)
    # p^T (N_rows A) = (N_rows^T p)^T A: one real mat-vec, then one complex one, no copies
    w    = _is_unconstrained(N) ? nothing : transpose(view(N, rows, :)) * p
    return w === nothing ? vec(transpose(p) * view(A, rows, :)) : transpose(A) * w
end

# -----------------------------------------------------------------------------------------
# angleSlope -- MATLAB public method
# -----------------------------------------------------------------------------------------

"""
    angleSlope(s::RitzSolver, modeNos, dKdC_cell = nothing) -> Vector{Float64}

Angular slope `dv/dtheta` [km/s per deg] of the selected modes at this angle
(Hellmann-Feynman, no extra solve). Port of `RitzSolver.angleSlope`.

`K` is linear in each layer's 15 rotated reduced constants `c`, so
`dK/dtheta = sum_c (dc/dtheta) dK/dc` exactly; only `dc/dtheta` is a central difference
(`+-0.01` deg, on the constants alone). Then
`dv/dtheta = Re(a^H dK_con/dtheta a) / (8 pi^2 v)`; exact for a non-degenerate mode.
`dKdC_cell` -- [`derivsByDepth`](@ref)'s first output, as the solve used (built if omitted).
Needs the full eigenvector set (error before `pruneEigVecs`).
"""
function angleSlope(s::RitzSolver, modeNos, dKdC_cell = nothing)
    size(s.eigVecs_con, 2) < s.nCompModes &&
        error("RitzSolver:angleSlope:noEigVecs -- the full Ritz eigenvector set is needed " *
              "(compute before retentionLevel pruning).")
    modes = collect(Int, modeNos isa Integer ? [modeNos] : modeNos)
    isempty(modes) && return Float64[]
    dKdC = dKdC_cell === nothing ? first(derivsByDepth(s.structure)) : dKdC_cell
    h  = 0.01
    cP = RedElastCellFor(s.structure, s.angle + h)
    cM = RedElastCellFor(s.structure, s.angle - h)
    cps    = s.structure.calcParams
    ranges = _layer_ranges(cps)
    U  = _full_eigvecs(s.nullConstrMat, view(s.eigVecs_con, :, modes))
    acc = zeros(Float64, length(modes))
    for (il, r) in pairs(ranges)
        bVL = cps[il].basisVecL
        dc  = (cP[il] .- cM[il]) ./ (2h)
        dKl = reshape(reshape(dKdC[il], bVL * bVL, 15) * dc, bVL, bVL)
        dKl .*= LayerNormFac(s.structure, il)
        Ub  = r == axes(U, 1) ? U : U[r, :]
        KU  = dKl * Ub
        for (k, _) in pairs(modes)
            acc[k] += real(dot(view(Ub, :, k), view(KU, :, k)))
        end
    end
    return acc ./ s.vels[modes] ./ (8 * π^2)
end

# -----------------------------------------------------------------------------------------
# pruneEigVecs -- MATLAB public method
# -----------------------------------------------------------------------------------------

"""
    pruneEigVecs(s::RitzSolver, keepCols = nothing) -> typeof(s)

Copy of `s` keeping only the named columns of `eigVecs_con` (`nothing`/empty: clear it).
`vels`, `eigVals` and the per-mode diagnostics are untouched. Port of
`RitzSolver.pruneEigVecs`; returns a new object (immutable struct, value semantics as in
MATLAB). `ArgumentError` for a non-integer or non-positive index.
"""
function pruneEigVecs(s::RitzSolver, keepCols = nothing)
    if keepCols === nothing || (keepCols isa AbstractArray && isempty(keepCols))
        return _rebuild(s; eigVecs_con = zeros(ComplexF64, 0, 0))
    end
    cols = keepCols isa Integer ? [Int(keepCols)] : collect(keepCols)
    all(c -> c isa Integer || isinteger(c), cols) ||
        throw(ArgumentError("RitzSolver.pruneEigVecs: keepCols must contain integers"))
    all(c -> c >= 1, cols) ||
        throw(ArgumentError("RitzSolver.pruneEigVecs: keepCols must be positive"))
    return _rebuild(s; eigVecs_con = s.eigVecs_con[:, Int.(cols)])
end

"Validate the shared constructor arguments (MATLAB `RitzSolver:constructor:badAngle`)."
function _check_angle(angleDeg)
    angleDeg === nothing && return 0.0
    (angleDeg isa Real && isfinite(angleDeg)) ||
        throw(ArgumentError("RitzSolver:constructor:badAngle -- angleDeg must be a real finite scalar."))
    return Float64(angleDeg)
end
