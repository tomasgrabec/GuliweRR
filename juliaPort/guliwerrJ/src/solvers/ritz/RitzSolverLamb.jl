# =========================================================================================
# RitzSolverLamb.jl -- port of guliwerr/solvers/ritz/RitzSolverLamb.m
#
# Port of guliwerr/solvers/ritz/RitzSolverLamb.m: the free-plate (Lamb-wave) leg of the
# `RitzSolverSAW < RitzSolver` / `RitzSolverLamb < RitzSolver` inheritance chain.
#
# NO per-file module, immutable struct, abstract-supertype-plus-duplicated-fields -- all three
# decisions are documented once in `RitzSolver.jl`'s header.
#
# TRANSLATION DECISIONS
# 1. **`withModeTypes` defaults to `false`.**
# 2. **`modeRecognition` classifies by mid-plane parity fractions** ([`parityFractions`](@ref)).
# 3. **`locateZGV` is the same Newton refinement as MATLAB** ([`ZgvNewton`](@ref)).
# 4. **`fastCalcGroupVel` always returns the `(velGroup, solver)` pair.** MATLAB switches on
#    `nargout` to skip building the second output; Julia has no `nargout`, and here there is
#    nothing to skip anyway -- the solver object IS the thing the group velocity is read off,
#    so returning it costs one pointer.
# =========================================================================================

using LinearAlgebra

using ...CoreTypes: Structure, CalcParams, Material, AddLayer, ModifyLayer, nlay,
                     set_rel_depth

export RitzSolverLamb, modeRecognition, parityFractions, fastCalcGroupVel, locateZGV,
       calcTransGroupVelocity, kPolynomial

# -----------------------------------------------------------------------------------------
# The type
# -----------------------------------------------------------------------------------------

"""
    RitzSolverLamb

Ritz-Rayleigh solver for the free plate (Lamb waves).

Port of `guliwerr/solvers/ritz/RitzSolverLamb.m`. **Immutable struct**, subtype of the
abstract [`RitzSolver`](@ref); the first 13 fields are `RitzSolver.m`'s own properties,
declared here because Julia abstract types carry no data (see `RitzSolver.jl`'s header).

Own fields beyond the shared ones:
- `modesType::Vector{Int}` -- `1 = Asymm`, `2 = Symm`, `3 = TH-Asymm`, `4 = TH-Symm`: the
  largest of the four mid-plane parity fractions ([`parityFractions`](@ref))
- `modeTypeConfidence::Vector{Float64}` -- top parity fraction minus the runner-up, in `[0,1]`.
  1 = a pure type; low = a mixed mode: parity mixing (a stack or crystal cut without a
  mid-plane mirror) or, within one parity class, SH-Lamb polarisation coupling.
- `realThickness::Vector{Float64}` -- physical plate thickness (per layer)
- `realWavVec::Float64` -- physical wavevector `k` [rad/m, or any consistent unit]
- `realWavelength::Float64` -- `2pi/realWavVec` (`NaN` when `realWavVec == 0`)
- `freqsPhysical::Vector{Float64}` -- real frequencies `vels * realWavVec / (2pi)`
- `velsGroupTrans::Vector{Float64}` -- transverse group velocity `c_gtheta = (1/k) domega/dtheta
  = dvp/dtheta` [same unit as `vels`] per mode: the in-plane group-velocity component
  perpendicular to k, positive towards increasing angle (`velsGroup` is the one along k);
  zero on a direction with a sagittal mirror ([`calcTransGroupVelocity`](@ref))
"""
struct RitzSolverLamb <: RitzSolver
    # --- RitzSolver.m's own properties (duplicated, not inherited: Julia abstract types carry no data)
    structure::Structure
    angle::Float64
    nCompModes::Int
    freqs::Vector{Float64}
    vels::Vector{Float64}
    velsGroup::Vector{Float64}
    eigVecs_con::Matrix{ComplexF64}
    eigVals::Vector{Float64}
    nullConstrMat::Matrix{Float64}
    basisVecL_con::Int
    dVel_dDepth::Vector{Float64}
    depthInsensitivity::Vector{Float64}
    # --- RitzSolverLamb.m's own
    modesType::Vector{Int}
    modeTypeConfidence::Vector{Float64}
    realThickness::Vector{Float64}
    realWavVec::Float64
    realWavelength::Float64
    freqsPhysical::Vector{Float64}
    velsGroupTrans::Vector{Float64}
end

"""
    RitzSolverLamb(; structure = nothing, angleDeg = nothing,
                     realThickness = nothing, realWavVec = nothing) -> RitzSolverLamb

Construct an unsolved `RitzSolverLamb`.

Port of `RitzSolverLamb.m` (and, through it, `RitzSolver.m`). All four arguments
are optional, mirroring MATLAB's `nargin`/`~isempty` guards.

**Supplying `realWavVec` (non-zero) rescales the structure.** `relDepth` is redefined so that
`k*d = 2*pi*relDepth`, i.e. `relDepth = d/lambda`: every layer's `relDepth` is multiplied by
`(sum(realThickness) * realWavVec / 2pi) / sum(relDepth)`, and the rescale is skipped entirely
when that factor is within `1e-14` of 1 (`RitzSolverLamb.m`). This is what makes a
wavevector sweep work by constructing one solver per `k`.

Keyword-only, not positional -- see [`RitzSolverSAW`](@ref)'s constructor for why.
"""
function RitzSolverLamb(; structure::Union{Structure,Nothing} = nothing,
                           angleDeg::Union{Real,Nothing} = nothing,
                           realThickness = nothing,
                           realWavVec::Union{Real,Nothing} = nothing)
    rt = realThickness === nothing ? Float64[] :
         (realThickness isa Real ? [Float64(realThickness)] :
                                    Float64.(collect(realThickness)))

    s = RitzSolverLamb(
        structure === nothing ? Structure() : structure,
        _check_angle(angleDeg),
        0, Float64[], Float64[], Float64[],
        zeros(ComplexF64, 0, 0), Float64[],
        zeros(Float64, 0, 0), 0,
        Float64[], Float64[],
        Int[], Float64[], rt, 0.0, NaN, Float64[], Float64[])

    if realWavVec !== nothing && realWavVec != 0
        isempty(s.realThickness) &&
            throw(ArgumentError("RitzSolverLamb:constructor:needThickness -- a physical " *
                                "wavevector needs the physical thickness (realThickness)."))
        k = Float64(realWavVec)
        s = _rebuild(s; realWavVec = k, realWavelength = 2π / k)

        # scale relDepth so that k*d = 2*pi*relDepth  (relDepth = d/lambda)
        if nlay(s.structure) > 0
            totalThick    = sum(s.realThickness)
            totalRelDepth = sum(cp.relDepth for cp in s.structure.calcParams)
            scaleFactor   = (totalThick * k / (2π)) / totalRelDepth
            if abs(scaleFactor - 1) > 1e-14
                newStr = s.structure
                for ilay in 1:nlay(newStr)
                    newRD = newStr.calcParams[ilay].relDepth * scaleFactor
                    newCP = set_rel_depth(newStr.calcParams[ilay], newRD)
                    newStr = ModifyLayer(newStr, ilay, nothing, newCP)
                end
                s = _rebuild(s; structure = newStr)
            end
        end
    end
    return s
end

"""
    _rebuild(s::RitzSolverLamb; <field> = s.<field>, ...) -> RitzSolverLamb

Copy-with-replacements (see `RitzSolver.jl`'s header). Internal; not exported.
"""
function _rebuild(s::RitzSolverLamb;
                   structure          = s.structure,
                   angle              = s.angle,
                   nCompModes         = s.nCompModes,
                   freqs              = s.freqs,
                   vels               = s.vels,
                   velsGroup          = s.velsGroup,
                   eigVecs_con        = s.eigVecs_con,
                   eigVals            = s.eigVals,
                   nullConstrMat      = s.nullConstrMat,
                   basisVecL_con      = s.basisVecL_con,
                   dVel_dDepth        = s.dVel_dDepth,
                   depthInsensitivity = s.depthInsensitivity,
                   modesType          = s.modesType,
                   modeTypeConfidence = s.modeTypeConfidence,
                   realThickness      = s.realThickness,
                   realWavVec         = s.realWavVec,
                   realWavelength     = s.realWavelength,
                   freqsPhysical      = s.freqsPhysical,
                   velsGroupTrans     = s.velsGroupTrans)
    return RitzSolverLamb(structure, angle, nCompModes, freqs, vels, velsGroup,
                           eigVecs_con, eigVals, nullConstrMat, basisVecL_con,
                           dVel_dDepth, depthInsensitivity,
                           modesType, modeTypeConfidence, realThickness, realWavVec,
                           realWavelength, freqsPhysical, velsGroupTrans)
end

# -----------------------------------------------------------------------------------------
# calc -- RitzSolverLamb.m
# -----------------------------------------------------------------------------------------

"""
    calc(s::RitzSolverLamb, withModeTypes::Bool = false, precomp = NamedTuple()) -> RitzSolverLamb

Calculate the Lamb modes of a free plate. Port of `RitzSolverLamb.calc`.

`precomp` may carry `nullConstrMat`, `rawByDepthCell` (from [`RawCellFor`](@ref); legOrder
only, shareable across a whole sweep) and `gammaByDepthCell` (from [`elastByDepth`](@ref)),
each read lazily. `rawByDepthCell` is built here rather than inside the Lamb `calc` of
`GammaMatrix` because [`calcTransGroupVelocity`](@ref) needs it too.

A single layer has no constraints (`N = I`): the projection is skipped, the real gauge is
tested on the full `K`, and the standard (not generalised) solve is used, as in MATLAB.
"""
function calc(s::RitzSolverLamb, withModeTypes::Bool = false,
               precomp::NamedTuple = NamedTuple())
    nullConstrMat = precomp_field(precomp, :nullConstrMat,
        () -> buildConstraintMatrix(s.structure.calcParams, LambBoundary()))
    rawByDepthCell   = precomp_field(precomp, :rawByDepthCell, () -> RawCellFor(s.structure))
    gammaByDepthCell = precomp_field(precomp, :gammaByDepthCell,
        () -> elastByDepth(s.structure, s.angle, rawByDepthCell))

    # K, dK/dD, M: a local, not stored (see the RitzSolver field table)
    gm = calc(s.structure, LambBoundary(), s.angle; gammaByDepthCell = gammaByDepthCell)
    out = _rebuild(s; nullConstrMat = nullConstrMat)

    K_con, M_con, gauge = _project_for_solve(out, gm)
    out = solveEigenvalueProblem(out, K_con, M_con, gauge)
    out = calcGroupVelocity(out, gm, false)
    out = calcTransGroupVelocity(out, rawByDepthCell)

    if out.realWavVec != 0
        out = _rebuild(out; freqsPhysical = out.vels .* out.realWavVec ./ (2π))
    end

    if withModeTypes && out.nCompModes > 0
        out = modeRecognition(out)
    end
    return out
end

"""
    calcTransGroupVelocity(s::RitzSolverLamb, rawByDepthCell = nothing) -> RitzSolverLamb

Transverse group velocity `c_gtheta` of every computed mode (Hellmann-Feynman in angle), filling
`velsGroupTrans`. Port of `RitzSolverLamb.calcTransGroupVelocity`.

`c_gtheta = (1/k) domega/dtheta = dvp/dtheta` (theta in rad, k fixed): the in-plane
group-velocity component perpendicular to k, positive towards increasing angle. The same
derivative as [`angleSlope`](@ref) (`dv/dtheta = a^H (dK/dtheta) a / (8 pi^2 v)`, `dK/dtheta`
from the `+-0.01` deg central difference of the rotated constants alone), built the Lamb path's
way: K is linear in the 15 reduced constants, so [`elastByDepth`](@ref) contracted with
`dc/dtheta` gives the three depth pages of `dK/dtheta` in one contraction; then the depth
factors, density and layer normalisation of the Lamb `calc`. No extra solve. Agrees with a
`+-0.01` deg finite difference of the solved velocities and, to four decimals on a Cu(001)
plate, with GEWtool's transverse energy velocity.
`rawByDepthCell` -- [`RawCellFor`](@ref)`(s.structure)`, built if omitted.
"""
function calcTransGroupVelocity(s::RitzSolverLamb, rawByDepthCell = nothing)
    raw = rawByDepthCell === nothing ? RawCellFor(s.structure) : rawByDepthCell
    nM = s.nCompModes
    nM == 0 && return _rebuild(s; velsGroupTrans = Float64[])
    h  = 0.01   # deg, as angleSlope
    cP = RedElastCellFor(s.structure, s.angle + h)
    cM = RedElastCellFor(s.structure, s.angle - h)
    dc = [(cP[il] .- cM[il]) ./ (2h) for il in 1:nlay(s.structure)]
    dPages = elastByDepth(s.structure, s.angle, raw, dc)

    cps    = s.structure.calcParams
    ranges = _layer_ranges(cps)
    U   = _full_eigvecs(s.nullConstrMat, s.eigVecs_con[:, 1:nM])
    acc = zeros(Float64, nM)
    for (il, r) in pairs(ranges)
        pg = dPages[il]
        dK, _ = ContractDepth(pg[:, :, 1], pg[:, :, 2], pg[:, :, 3], cps[il].relDepth)
        dK = dK ./ s.structure.mater[il].rho .* LayerNormFac(s.structure, il)
        Ub = U[r, :]
        KU = dK * Ub
        for k in 1:nM
            acc[k] += real(dot(view(Ub, :, k), view(KU, :, k)))
        end
    end
    slopePerDeg = acc ./ s.vels[1:nM] ./ (8 * π^2)
    return _rebuild(s; velsGroupTrans = slopePerDeg .* (180 / π))
end

"""
    modeRecognition(s::RitzSolverLamb) -> RitzSolverLamb

Classify every computed mode as Asymm (1) / Symm (2) / TH-Asymm (3) / TH-Symm (4) by mid-plane
parity, filling `modesType` (argmax of [`parityFractions`](@ref)) and `modeTypeConfidence`
(largest fraction minus the second largest).

Port of `RitzSolverLamb.modeRecognition`. Physically defined, so no reference plate: a mode is
labelled by where its kinetic-energy-weighted displacement sits once each component is split
into its even and odd parts about the plate's mid-plane (`docs/reference/architecture/plate.md`,
"Mode-type classification").
"""
function modeRecognition(s::RitzSolverLamb)
    nM = s.nCompModes
    V  = s.eigVecs_con[:, 1:nM]
    U  = nlay(s.structure) > 1 ? s.nullConstrMat * V : V    # one layer: N = I
    fr = parityFractions(U, s.structure)
    modesType          = zeros(Int, nM)
    modeTypeConfidence = zeros(Float64, nM)
    for im in 1:nM
        f = @view fr[:, im]
        i1 = argmax(f)                              # first maximum, as MATLAB's sort/max
        f2 = maximum(f[j] for j in 1:4 if j != i1)
        modesType[im]          = i1
        modeTypeConfidence[im] = max(f[i1] - f2, 0.0)
    end
    return _rebuild(s; modesType = modesType, modeTypeConfidence = modeTypeConfidence)
end

# -----------------------------------------------------------------------------------------
# parityFractions -- MATLAB `methods (Static)`
# -----------------------------------------------------------------------------------------

"""
    parityFractions(U, structure) -> Matrix{Float64}   # [4 x nModes]

Energy fractions of the four mid-plane parity/polarisation classes of every mode, rows in
`modesType` order, each column summing to 1:
`1 Asymm = E(u1 odd) + E(u3 even)`, `2 Symm = E(u1 even) + E(u3 odd)`, `3 TH-Asymm = E(u2 odd)`,
`4 TH-Symm = E(u2 even)`, with `E(.) = int rho_s |.|^2 dz / int rho_s |u|^2 dz`, even/odd about
the mid-plane and `rho_s(z) = (rho(z) + rho(D - z))/2`.

Port of `RitzSolverLamb.parityFractions`. `U` holds the eigenvectors in the full Legendre
basis (`nullConstrMat * eigVecs_con`).

The symmetrised density keeps the even and odd parts orthogonal, so the fractions sum to 1
exactly; it is the kinetic-energy weight itself for a mirror-symmetric stack (and for one
layer, where `rho` cancels). One layer: the normalised Legendre polynomials are orthonormal on
the layer and `P_n(-x) = (-1)^n P_n(x)` about its centre, so the even/odd energies are sums of
`|coefficient|^2` over even/odd `n` -- no depth grid. Several layers: the displacement at depth
`z` and at its mirror `D - z` on a Gauss-Legendre grid over the lower half of the plate whose
sub-intervals end at every interface and mirrored interface; `max(legOrder)` nodes per
sub-interval integrate the (polynomial) products exactly.
"""
function parityFractions(U::AbstractMatrix, structure::Structure)
    cps  = structure.calcParams
    nLay = length(cps)
    nM   = size(U, 2)
    E    = zeros(6, nM)    # rows: u1e u1o u2e u2o u3e u3o
    if nLay == 1
        L = cps[1].legOrder
        @inbounds for m in 1:nM, ic in 1:3, n in 1:L
            a2 = abs2(U[(ic - 1) * L + n, m])
            E[2ic - (isodd(n) ? 1 : 0), m] += a2       # n = 1, 3, ... is degree 0, 2, ...: even
        end
    else
        zLow, w = _mirrorQuadrature(cps)
        zMir = 1 .- zLow
        tB   = vcat(0.0, cumsum([cp.relDepth for cp in cps])) ./ sum(cp.relDepth for cp in cps)
        rhoL = [m.rho for m in structure.mater]
        layLow = _layerOf(zLow, tB)
        layMir = _layerOf(zMir, tB)
        w = w .* (rhoL[layLow] .+ rhoL[layMir]) ./ 2
        row0 = vcat(0, cumsum([cp.basisVecL for cp in cps]))
        nQ = length(zLow)
        uLow = zeros(eltype(U), nQ, nM)
        uMir = zeros(eltype(U), nQ, nM)
        for ic in 1:3
            fill!(uLow, 0); fill!(uMir, 0)
            for ilay in 1:nLay
                L    = cps[ilay].legOrder
                rows = row0[ilay] .+ (ic - 1) * L .+ (1:L)
                q = findall(==(ilay), layLow)
                isempty(q) || (uLow[q, :] = _normLegendre(L, _localCoord(zLow[q], tB, ilay)) * U[rows, :])
                q = findall(==(ilay), layMir)
                isempty(q) || (uMir[q, :] = _normLegendre(L, _localCoord(zMir[q], tB, ilay)) * U[rows, :])
            end
            for m in 1:nM
                se = 0.0; so = 0.0
                @inbounds for k in 1:nQ
                    se += w[k] * abs2(uLow[k, m] + uMir[k, m])
                    so += w[k] * abs2(uLow[k, m] - uMir[k, m])
                end
                E[2ic - 1, m] = se
                E[2ic, m]     = so
            end
        end
    end
    fr = vcat(E[2:2, :] .+ E[5:5, :], E[1:1, :] .+ E[6:6, :], E[4:4, :], E[3:3, :])
    return fr ./ max.(sum(fr, dims = 1), floatmin(Float64))
end

# Gauss-Legendre nodes/weights on the lower half [0, 1/2] of the normalised thickness,
# sub-intervals ending at every interface and mirrored interface (MATLAB mirrorQuadrature).
function _mirrorQuadrature(cps)
    t = vcat(0.0, cumsum([cp.relDepth for cp in cps])) ./ sum(cp.relDepth for cp in cps)
    b = sort(vcat(min.(t, 1 .- t), 0.5))
    b = b[vcat(true, diff(b) .> 1e-12)]            # merge coincident breakpoints
    n = maximum(cp.legOrder for cp in cps)
    beta = [k / sqrt(4k^2 - 1) for k in 1:n-1]     # Golub-Welsch
    F  = eigen(SymTridiagonal(zeros(n), beta))
    x  = F.values
    wx = 2 .* F.vectors[1, :] .^ 2
    h  = diff(b)
    z  = vec(b[1:end-1]' .+ (x .+ 1) ./ 2 .* h')   # [n x nSub], column-major as MATLAB
    w  = vec(wx ./ 2 .* h')
    return z, w
end

# Layer index of each normalised depth z (tB = layer boundaries, bottom = 0)
_layerOf(z, tB) = [1 + count(j -> zi >= tB[j], 2:length(tB)-1) for zi in z]

_localCoord(z, tB, ilay) = 2 .* (z .- tB[ilay]) ./ (tB[ilay + 1] - tB[ilay]) .- 1

# Normalised Legendre polynomials sqrt((2n+1)/2) P_n(x), n = 0..L-1, [length(x) x L]
function _normLegendre(L::Integer, x::AbstractVector)
    P = ones(length(x), L)
    L > 1 && (P[:, 2] .= x)
    for n in 1:L-2
        @views P[:, n + 2] .= ((2n + 1) .* x .* P[:, n + 1] .- n .* P[:, n]) ./ (n + 1)
    end
    return P .* sqrt.(((0:L-1) .+ 0.5))'
end

# -----------------------------------------------------------------------------------------
# fastCalcGroupVel / locateZGV -- MATLAB `methods (Static)`, RitzSolverLamb.m
# -----------------------------------------------------------------------------------------

"""
    fastCalcGroupVel(structure, angleDeg, realThickness, modeNoIn, wavVecIn,
                      precomp = NamedTuple()) -> (velGroup, solver)

Group velocity of one Lamb branch at one wavevector: build a [`RitzSolverLamb`](@ref) at
`wavVecIn`, solve it with mode recognition off, and read `velsGroup[modeNoIn]`. A thin
wrapper for checks.

Port of `RitzSolverLamb.fastCalcGroupVel`. MATLAB returns the
solver object only when asked for a second output; Julia has no `nargout`, and the object is
built regardless, so it is always returned (this file's TRANSLATION DECISIONS note 4).
"""
function fastCalcGroupVel(structure::Structure, angleDeg::Real, realThickness,
                           modeNoIn::Integer, wavVecIn::Real,
                           precomp::NamedTuple = NamedTuple())
    fc = RitzSolverLamb(structure = structure, angleDeg = angleDeg,
                         realThickness = realThickness, realWavVec = wavVecIn)
    fc = calc(fc, false, precomp)
    return fc.velsGroup[modeNoIn], fc
end

"""
    locateZGV(structure, angleDeg, realThickness, modeNoIn, wavVec0,
               precomp = NamedTuple()) -> (wavVecZGV, freqZGV, solver)

Locate the zero-group-velocity (ZGV) point nearest a seed on one Lamb mode, by Newton
refinement ([`ZgvNewton`](@ref)`.solve`) of the plate operator written as a polynomial in k
([`kPolynomial`](@ref)), seeded with raw mode `modeNoIn`'s frequency and eigenvector at
`wavVec0`. Port of `RitzSolverLamb.locateZGV`.

No mode number is held fixed, so a branch crossing between the seed and the ZGV does not
matter; `freqZGV` is the Newton frequency, and the mode it belongs to in `solver` (solved at
`wavVecZGV`) need not be `modeNoIn`. Returns `NaN, NaN, solver-at-wavVec0` when Newton does not
converge or lands outside `[wavVec0/10, wavVec0*10]` (the search interval).

`precomp` -- `nullConstrMat`, `rawByDepthCell` (built if absent).
"""
function locateZGV(structure::Structure, angleDeg::Real, realThickness,
                    modeNoIn::Integer, wavVec0::Real,
                    precomp::NamedTuple = NamedTuple())
    N   = precomp_field(precomp, :nullConstrMat,
        () -> buildConstraintMatrix(structure.calcParams, LambBoundary()))
    raw = precomp_field(precomp, :rawByDepthCell, () -> RawCellFor(structure))
    pc  = (; nullConstrMat = N, rawByDepthCell = raw)
    rs0 = calc(RitzSolverLamb(structure = structure, angleDeg = angleDeg,
                               realThickness = realThickness, realWavVec = wavVec0), false, pc)
    poly = kPolynomial(structure, angleDeg, realThickness, pc)
    kZ, omegaZ, _, ok, _ = ZgvNewton.solve(poly.L2, poly.L1, poly.L0, poly.M,
                                           Float64(wavVec0), 2π * rs0.freqsPhysical[modeNoIn])
    if !ok || kZ < wavVec0 / 10 || kZ > wavVec0 * 10
        return NaN, NaN, rs0
    end
    objOut = calc(RitzSolverLamb(structure = structure, angleDeg = angleDeg,
                                  realThickness = realThickness, realWavVec = kZ), false, pc)
    return kZ, omegaZ / (2π), objOut
end

"""
    kPolynomial(structure, angleDeg, realThickness, precomp = NamedTuple())
        -> (; L2, L1, L0, M, N)

The plate's Ritz eigenproblem as a quadratic polynomial in the wavevector:
`(k^2 L2 + k L1 + L0) a = omega^2 M a` at every physical wavevector `k` with `omega = 2 pi f`
(`f` = `freqsPhysical`), in the constrained basis `N` (`nullConstrMat`). All Hermitian
(`L2`, `L0`, `M` real in every case checked, `L1` imaginary). Port of
`RitzSolverLamb.kPolynomial`.

At fixed physical layer thicknesses `h_i` (only `sum(realThickness)` is used, the layers keep
the relDepth ratios of `structure.calcParams`, as in the constructor), [`ContractDepth`](@ref)'s
`K = K0 + K1/(relDepth/2) + K2/(relDepth/2)^2` rescaled from `(2 pi v)^2` to `omega^2` splits
per layer into `L2 = K0/(2pi)^2`, `L1 = K1/(pi h_i)`, `L0 = 4 K2/h_i^2`, each times
`normFac_i/rho_i`, and `M` = the `normFac`-weighted identity.
"""
function kPolynomial(structure::Structure, angleDeg::Real, realThickness,
                      precomp::NamedTuple = NamedTuple())
    N   = precomp_field(precomp, :nullConstrMat,
        () -> buildConstraintMatrix(structure.calcParams, LambBoundary()))
    raw = precomp_field(precomp, :rawByDepthCell, () -> RawCellFor(structure))
    pages = elastByDepth(structure, angleDeg, raw)
    cps  = structure.calcParams
    n    = nlay(structure)
    relD = [cp.relDepth for cp in cps]
    h    = sum(realThickness) .* relD ./ sum(relD)
    rho  = [structure.mater[il].rho for il in 1:n]
    nTot = sum(cp.basisVecL for cp in cps)
    b2 = zeros(ComplexF64, nTot, nTot); b1 = zeros(ComplexF64, nTot, nTot)
    b0 = zeros(ComplexF64, nTot, nTot); bm = zeros(Float64, nTot, nTot)
    for (il, r) in pairs(_layer_ranges(cps))
        nf = (h[il] * rho[il]) / (h[1] * rho[1])
        w  = nf / rho[il]
        pg = pages[il]
        b2[r, r] .= pg[:, :, 1] .* (w / (2π)^2)
        b1[r, r] .= pg[:, :, 2] .* (w / (π * h[il]))
        b0[r, r] .= pg[:, :, 3] .* (w * 4 / h[il]^2)
        for j in r; bm[j, j] = nf; end
    end
    sym(A) = (A + A') / 2
    return (; L2 = sym(transpose(N) * b2 * N), L1 = sym(transpose(N) * b1 * N),
              L0 = sym(transpose(N) * b0 * N), M = real(sym(transpose(N) * bm * N)), N = N)
end
