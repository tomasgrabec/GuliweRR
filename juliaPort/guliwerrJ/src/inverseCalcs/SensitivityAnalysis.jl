# src/inverseCalcs/SensitivityAnalysis.jl
#
# EXPERIMENTAL (as in MATLAB): validated on synthetic data only; the interface may change.
#
# Port of `guliwerr/inverseCalcs/SensitivityAnalysis.m`: the sensitivity /
# Jacobian matrix G, its SVD-based eigendecomposition (transformMatrix T, sensitivity
# eigenvalues psi = 2 s.^2), and the cutoff deciding how many transformed constants are
# optimised.
#
# Row k of G is one data point (an associated Ritz mode, a limiting bulk velocity, a
# pulse-echo phase velocity, or the bulk-modulus pseudo-datum); column i one independent
# constant of the inverted layer (ElastConsts.vector order); G_ki = dv_k/dc_i at the guess.
#
# - Ritz rows: Rayleigh quotients Re[alpha' (N' dK/dc_i N) alpha] / (8 pi^2 v)
#   (`ritzVelocityDerivative`). The Julia port evaluates them on the inverted layer's
#   full-space eigenvector block U = N[rows, :] * alpha (U' dG U), never projecting dK/dc:
#   same numbers to round-off, O(bVL^2) instead of O(nCon^2) per constant.
# - Limiting / phase / bulk rows: ONE central-difference loop (`channelDerivative`) over the
#   SAME forward evaluators the inversion's residual uses (`channelValue`), with
#   `limitingStep` for the limiting channel and `phaseBulkStep` for the other two. A limiting
#   root lost inside a difference is an error here (`limitingRootLost`).
#
# MATLAB is a value class built once by the static `calc`; Julia: an immutable struct and a
# method of the shared `calc` generic. `association` is duck-typed on `.points` (see
# `ModeAssociation.jl`'s `pointsAtAngle`).
#
# As in MATLAB: no `sensitivityHessian` or `invertedLayerNo` fields, validated options, one Ritz-class set (`activeRitzClasses`), one
# dataset check (`checkDatasets`).

using LinearAlgebra: svd, Diagonal, mul!

using ..CoreTypes: Structure, Material, nlay, ElastConsts, SetElastConsts, SetDirMat,
                    GetRotDirMat, GetStandardDirMat, GetRotElastTensor, Any2Tensor, Voigt2Bulk

import ..Solvers: calc
using ..Solvers: dummyForLayer, matrix, refineRoot, LimitingRoot, buildConstraintMatrix,
                 derivsByDepth, SAWBoundary

using ..ForwardCalcs: limitingAndBulk, resolveParallel

export SensitivityAnalysis, DEFAULT_SENSITIVITY_CUTOFF, waveTypeIndex, phaseVelocityAt,
       bulkVelocityAt, limitingStep, phaseBulkStep, limitingRootsAt, limitingVelocityNear,
       limitingVelocityTracked, rotatedTensorFor, channelDerivative, channelValue,
       ritzVelocityDerivative, checkDatasets

"""Default of the `sensitivityCutoff` option (error amplification <= 160x; empirical,
Grabec et al. 2024)."""
const DEFAULT_SENSITIVITY_CUTOFF = 1 / 160

"""
    SensitivityAnalysis

The sensitivity matrix `G` (after the point weights), `psi`, `T`, the normalised roots
`sqrt(psi/psi[1])` and the automatic cutoff `nOptimisedConstants`. Fields as MATLAB's
properties: `nIndepConsts, symmClassNo, nDatasets, initialVoigtVec, initialTransformedVec
(= T' c0), sensitivityMatrix, sensitivityEigenvalues, transformMatrix,
normalisedSensitivityRoots, sensitivityCutoff, nOptimisedConstants, pointWeights,
includePseudoSAW, nRitzRows, nLimitingRows, nPhaseRows, hasBulkModulusRow`.
"""
struct SensitivityAnalysis
    nIndepConsts::Int
    symmClassNo::Int
    nDatasets::Int

    initialVoigtVec::Vector{Float64}
    initialTransformedVec::Vector{Float64}

    sensitivityMatrix::Matrix{Float64}
    sensitivityEigenvalues::Vector{Float64}
    transformMatrix::Matrix{Float64}
    normalisedSensitivityRoots::Vector{Float64}

    sensitivityCutoff::Float64
    nOptimisedConstants::Int

    pointWeights::Vector{Float64}
    includePseudoSAW::Bool

    nRitzRows::Int
    nLimitingRows::Int
    nPhaseRows::Int
    hasBulkModulusRow::Bool
end

"""Empty `SensitivityAnalysis` (MATLAB's no-argument constructor)."""
SensitivityAnalysis() = SensitivityAnalysis(0, 0, 1, Float64[], Float64[], zeros(0, 0),
    Float64[], zeros(0, 0), Float64[], DEFAULT_SENSITIVITY_CUTOFF, 0, Float64[], true,
    0, 0, 0, false)

# -----------------------------------------------------------------------------------------
# Shared forward-model evaluators of the non-Ritz channels (also the inversion's residual)
# -----------------------------------------------------------------------------------------

"""
    waveTypeIndex(waveType) -> Int

Christoffel branch index `iwt` of a limiting/phase-velocity wave type (descending
eigenvalue order): longitudinal 1, fastShear 2, slowShear 3.
"""
function waveTypeIndex(waveType::AbstractString)
    waveType == "longitudinal" && return 1
    waveType == "fastShear" && return 2
    waveType == "slowShear" && return 3
    throw(ArgumentError("SensitivityAnalysis:waveTypeIndex:badWaveType -- \"$waveType\" has " *
        "no Christoffel branch index (expected longitudinal, fastShear or slowShear)."))
end

"""
    phaseVelocityAt(cVec, symmNo, dirVec, rho, waveType) -> Float64

Bulk phase velocity along `dirVec`, given in the inverted layer's CRYSTAL frame.
"""
function phaseVelocityAt(cVec::AbstractVector{<:Real}, symmNo::Integer,
                          dirVec::AbstractVector{<:Real}, rho::Real, waveType::AbstractString)
    idx = waveTypeIndex(waveType)
    eigVals, _, _ = matrix(Any2Tensor(cVec, symmNo), dirVec, rho)
    return sqrt(eigVals[idx])
end

"""
    bulkVelocityAt(cVec, symmNo, rho) -> Float64

Bulk-modulus pseudo-velocity `sqrt(B/rho)`.
"""
function bulkVelocityAt(cVec::AbstractVector{<:Real}, symmNo::Integer, rho::Real)
    B = Voigt2Bulk(ElastConsts(cVec, symmNo).voigt)
    return sqrt(B / rho)
end

"""
    limitingStep(cAbs) -> Float64

Central-difference step of the LIMITING channel [GPa], `max(1e-4|c|, 1e-3)`: with the
bracket held fixed and `refineRoot` at TolX 1e-12 the limiting velocity is smooth over a
wide range of h, but `1e-6|c|` sits at its noisy end.
"""
limitingStep(cAbs::Real) = max(1e-4 * abs(cAbs), 1e-3)

"""
    phaseBulkStep(cAbs) -> Float64

Central-difference step of the pulse-echo and bulk-modulus channels [GPa],
`max(1e-6|c|, 1e-6)` (both are analytic functions of the constants).
"""
phaseBulkStep(cAbs::Real) = max(1e-6 * abs(cAbs), 1e-6)

"""
    ritzVelocityDerivative(quadForms, vels)

Ritz-channel `dv/dc` from the Rayleigh-quotient numerators `Re[alpha' dK alpha]`
(M-orthonormal `alpha`): `quadForms ./ (8 pi^2 v)`.
"""
ritzVelocityDerivative(quadForms, vels) = quadForms ./ (8 * pi^2 .* vels)

"""
    rotatedTensorFor(matLay, symmNo, cVec, angleDeg) -> Array{Float64,4}

The rotated elastic tensor `limitingAndBulk` builds internally (propagation rotation of
the dirMat, then the standard-frame tensor at angle 0 -- the double-rotation rule); must
stay bit-identical to that chain, or a refined and a rescanned root would not agree.
"""
function rotatedTensorFor(matLay::Material, symmNo::Integer, cVec::AbstractVector{<:Real},
                           angleDeg::Real)
    matPerturbed = SetElastConsts(matLay, ElastConsts(cVec, symmNo))
    matRot = SetDirMat(matPerturbed, GetRotDirMat(matPerturbed, angleDeg))
    matStd = SetDirMat(matRot, GetStandardDirMat(matRot.dirMat))
    return GetRotElastTensor(matStd, 0)
end

"""
    limitingRootsAt(matLay, symmNo, cVec, angleDeg) -> Vector{LimitingRoot}

The full LABELLED Christoffel limiting-root set at one angle (every accepted sagittal root
of all three wave types), so that a caller can re-identify the ONE root it tracks. Never
reduce `shearLimVels` by max()/min() to pick a branch.
"""
function limitingRootsAt(matLay::Material, symmNo::Integer, cVec::AbstractVector{<:Real},
                          angleDeg::Real)
    matPerturbed = SetElastConsts(matLay, ElastConsts(cVec, symmNo))
    matRot = SetDirMat(matPerturbed, GetRotDirMat(matPerturbed, angleDeg))
    lim = limitingAndBulk(matRot)
    nShear = length(lim.shearLimVels)
    roots = Vector{LimitingRoot}(undef, nShear + 1)
    roots[1] = LimitingRoot(1, lim.longLimVel, lim.longLimTheta, lim.longLimPol)
    for q in 1:nShear
        roots[q + 1] = LimitingRoot(lim.shearLimIwt[q], lim.shearLimVels[q],
                                     lim.shearLimTheta[q], lim.shearLimPols[:, q])
    end
    return roots
end

"""
    limitingVelocityNear(matLay, symmNo, cVec, angleDeg, thetaAnchor, iwt, halfWidth = nothing)
        -> (v, theta, pol, ok)

ONE anchored limiting root re-found locally by `refineRoot` (bracket centred on the
accepted iterate's theta*, independent of `cVec`); `ok = false` = the root no longer exists.
"""
function limitingVelocityNear(matLay::Material, symmNo::Integer, cVec::AbstractVector{<:Real},
                               angleDeg::Real, thetaAnchor::Real, iwt::Integer,
                               halfWidth::Union{Nothing,Real} = nothing)
    rotElTensor = rotatedTensorFor(matLay, symmNo, cVec, angleDeg)
    v, theta, pol, _, ok = halfWidth === nothing ?
        refineRoot(rotElTensor, matLay.rho, thetaAnchor, iwt) :
        refineRoot(rotElTensor, matLay.rho, thetaAnchor, iwt, halfWidth)
    return v, theta, pol, ok
end

"""
    limitingVelocityTracked(matLay, symmNo, cVec, angleDeg, pt, halfWidth = nothing)
        -> (v, theta, ok)

Anchored local refinement with the mandatory full-rescan fallback (re-identification by
[`limitingRootScore`](@ref)); `ok = false` only if both failed. `pt` is one limiting point
(`limRootTheta`, `limBranchIwt`, `limAnchorPol`, `calculatedVelocity`).
"""
function limitingVelocityTracked(matLay::Material, symmNo::Integer, cVec::AbstractVector{<:Real},
                                  angleDeg::Real, pt, halfWidth::Union{Nothing,Real} = nothing)
    v, theta, _, ok = limitingVelocityNear(matLay, symmNo, cVec, angleDeg, pt.limRootTheta,
                                            Int(pt.limBranchIwt), halfWidth)
    ok && return v, theta, true
    roots = limitingRootsAt(matLay, symmNo, cVec, angleDeg)
    isempty(roots) && return NaN, Float64(pt.limRootTheta), false
    bestScore, q = findmax(limitingRootScore(pt, roots))
    bestScore <= 0 && return NaN, Float64(pt.limRootTheta), false
    return roots[q].vel, roots[q].theta, true
end

"""
    channelValue(kind::Symbol, cVec, ctx) -> (v, ok)

One non-Ritz data point's calculated velocity at the Voigt vector `cVec`. `kind`:
`:limiting` (ctx: `symmNo, matLay, pt`), `:phase` (ctx: `symmNo, rho, pv`), `:bulk` (ctx:
`symmNo, rho`). `ok = false` only for a limiting point whose root was lost (`v = NaN`).
"""
function channelValue(kind::Symbol, cVec::AbstractVector{<:Real}, ctx)
    if kind === :limiting
        v, _, ok = limitingVelocityTracked(ctx.matLay, ctx.symmNo, cVec, ctx.pt.angleDeg, ctx.pt)
        return v, ok
    elseif kind === :phase
        return phaseVelocityAt(cVec, ctx.symmNo, ctx.pv.propagationDir, ctx.rho,
                               ctx.pv.waveType), true
    elseif kind === :bulk
        return bulkVelocityAt(cVec, ctx.symmNo, ctx.rho), true
    end
    throw(ArgumentError("SensitivityAnalysis:channelValue:badKind -- \"$kind\" is not a " *
        "data channel (expected limiting, phase or bulk)."))
end

"""
    channelDerivative(kind::Symbol, x, nCols, toVoigt, ctx) -> (row::Vector{Float64}, ok)

Central-difference `+dv/dx_i`, `i = 1..nCols`, of one non-Ritz data point -- the one
finite-difference loop of `calc` (`x` = Voigt constants, `toVoigt = identity`) and of the
inversion's Jacobian (`x` = transformed constants, `toVoigt = c -> T*c`; its Jacobian row
is `-row`). Step `limitingStep` (`:limiting`) or `phaseBulkStep` (`:phase`, `:bulk`).
`ok = false` if a perturbed limiting evaluation lost its root (that entry is NaN).
"""
function channelDerivative(kind::Symbol, x::AbstractVector{<:Real}, nCols::Integer, toVoigt, ctx)
    kind in (:limiting, :phase, :bulk) || throw(ArgumentError(
        "SensitivityAnalysis:channelDerivative:badKind -- \"$kind\" is not a data channel " *
        "(expected limiting, phase or bulk)."))
    row = zeros(nCols)
    ok = true
    xP = Vector{Float64}(x)
    for i in 1:nCols
        h = kind === :limiting ? limitingStep(x[i]) : phaseBulkStep(x[i])
        xP[i] = x[i] + h
        vP, okP = channelValue(kind, toVoigt(xP), ctx)
        xP[i] = x[i] - h
        vM, okM = channelValue(kind, toVoigt(xP), ctx)
        xP[i] = x[i]
        row[i] = (vP - vM) / (2h)
        ok = ok && okP && okM
    end
    return row, ok
end

"""
    checkDatasets(structure, association, measured, invertedLayerNo, idPrefix)
        -> (nDatasets, symmNoAll, LAll)

Dataset-array checks shared by `calc` and the `ElasticityInversion` constructor: equal
lengths (`<idPrefix>:datasetCountMismatch`), `invertedLayerNo` in range for every dataset
(`<idPrefix>:badLayer`); returns every dataset's inverted-layer symmClassNo / nIndepConsts.
Arguments are vectors (one entry per dataset).
"""
function checkDatasets(structure::AbstractVector, association::AbstractVector,
                        measured::AbstractVector, invertedLayerNo::Integer, idPrefix::AbstractString)
    nDatasets = length(structure)
    (length(association) == nDatasets && length(measured) == nDatasets) || throw(ArgumentError(
        "$idPrefix:datasetCountMismatch -- structure ($nDatasets), association " *
        "($(length(association))) and measured ($(length(measured))) must all supply the " *
        "same number of datasets."))
    (invertedLayerNo >= 1 && all(nlay(s) >= invertedLayerNo for s in structure)) ||
        throw(ArgumentError("$idPrefix:badLayer -- invertedLayerNo $invertedLayerNo is out of " *
            "range for at least one dataset (nLay = [$(join(nlay.(structure), " "))])."))
    symmNoAll = [s.mater[invertedLayerNo].elConsts.symmClassNo for s in structure]
    LAll = [s.mater[invertedLayerNo].elConsts.nIndepConsts for s in structure]
    return nDatasets, symmNoAll, LAll
end

# -----------------------------------------------------------------------------------------
# Ritz rows
# -----------------------------------------------------------------------------------------

"""
    layerRows(structure, invertedLayerNo) -> UnitRange{Int}

Full-space rows (Legendre coefficients) of the inverted layer.
"""
function layerRows(structure::Structure, invertedLayerNo::Integer)
    below = sum((structure.calcParams[i].basisVecL for i in 1:(invertedLayerNo - 1)); init = 0)
    return (below + 1):(below + structure.calcParams[invertedLayerNo].basisVecL)
end

"""
    ritzRowsAtAngle(structure, invertedLayerNo, dKdC_cell, N, angleDeg, association,
                    ritzIdxGlobal, includeSet) -> (rows, localIdx)

Sensitivity rows of every Ritz-bound point at one angle. `includeSet` must be the set the
caller built `ritzIdxGlobal` with. Quadratic forms on the layer block of the full-space
eigenvectors (see the file header).
"""
function ritzRowsAtAngle(structure::Structure, invertedLayerNo::Integer, dKdC_cell,
                          N::AbstractMatrix{Float64}, angleDeg::Real, association,
                          ritzIdxGlobal::AbstractVector{<:Integer}, includeSet)
    globalIdxAtAngle = pointsAtAngle(association, angleDeg, includeSet)
    localIdx = Int.(indexin(globalIdxAtAngle, ritzIdxGlobal))

    dG_str, _ = dummyForLayer(structure, invertedLayerNo, dKdC_cell, angleDeg)
    L = size(dG_str, 3)
    rows = layerRows(structure, invertedLayerNo)
    pts = association.points
    nPts = length(globalIdxAtAngle)
    A = Matrix{ComplexF64}(undef, size(N, 2), nPts)
    vels = Vector{Float64}(undef, nPts)
    for k in 1:nPts
        A[:, k] = pts[globalIdxAtAngle[k]].anchorEigVec_con
        vels[k] = pts[globalIdxAtAngle[k]].calculatedVelocity
    end
    U = N[rows, :] * A                      # full-space eigenvectors, inverted-layer block
    out = zeros(nPts, L)
    dGU = similar(U)
    for i in 1:L
        mul!(dGU, view(dG_str, rows, rows, i), U)
        for k in 1:nPts
            q = 0.0
            @inbounds for r in axes(U, 1)
                q += real(conj(U[r, k]) * dGU[r, k])
            end
            out[k, i] = ritzVelocityDerivative(q, vels[k])
        end
    end
    return out, localIdx
end

"""
    datasetSensitivityRows(structure, association, invertedLayerNo, measured, cVec0, symmNo,
                           L, includeSet, useParallel) -> (G_d, counts)

ONE dataset's block: Ritz rows (threaded over the unique angles, each angle writing its own
slot), then limiting, phase-velocity and bulk-modulus rows.
"""
function datasetSensitivityRows(structure::Structure, association, invertedLayerNo::Integer,
                                 measured::MeasuredVelocities, cVec0::Vector{Float64},
                                 symmNo::Integer, L::Integer, includeSet, useParallel::Symbol)
    matLay = structure.mater[invertedLayerNo]
    rho = matLay.rho
    pts = association.points

    ritzIdxGlobal = findall(p -> p.bindingClass in includeSet, pts)
    nRitz = length(ritzIdxGlobal)
    G_ritz = zeros(nRitz, L)
    if nRitz > 0
        N = buildConstraintMatrix(structure.calcParams, SAWBoundary())
        dKdC_cell, _ = derivsByDepth(structure)
        uniqueAngles = sort(unique(pts[i].angleDeg for i in ritzIdxGlobal))
        nA = length(uniqueAngles)
        rowsVec = Vector{Matrix{Float64}}(undef, nA)
        idxVec = Vector{Vector{Int}}(undef, nA)
        job(ia) = ritzRowsAtAngle(structure, invertedLayerNo, dKdC_cell, N, uniqueAngles[ia],
                                  association, ritzIdxGlobal, includeSet)
        if resolveParallel(useParallel, nA)
            with_blas_threads(1) do
                Threads.@threads for ia in 1:nA
                    rowsVec[ia], idxVec[ia] = with_task_blas_single(() -> job(ia))
                end
            end
        else
            for ia in 1:nA
                rowsVec[ia], idxVec[ia] = job(ia)
            end
        end
        for ia in 1:nA
            G_ritz[idxVec[ia], :] .= rowsVec[ia]
        end
    end

    # Limiting rows: anchored local refinement from each point's tracked theta*/branch at
    # the wider limiting step (a full rescan per perturbation, or the 1e-6|c| step, gives
    # WRONG derivatives, not merely slow ones).
    limIdxGlobal = findall(p -> p.bindingClass == "limiting", pts)
    nLim = length(limIdxGlobal)
    G_lim = zeros(nLim, L)
    for k in 1:nLim
        pt = pts[limIdxGlobal[k]]
        row, ok = channelDerivative(:limiting, cVec0, L, identity,
                                    (symmNo = symmNo, matLay = matLay, pt = pt))
        ok || throw(ArgumentError("SensitivityAnalysis:calc:limitingRootLost -- Limiting " *
            "point $(limIdxGlobal[k]) ($(pt.waveType) at $(pt.angleDeg) deg): the tracked " *
            "Christoffel root vanished within the finite-difference step, so its sensitivity " *
            "row is undefined. Remove the point or start from another initial guess."))
        G_lim[k, :] .= row
    end

    nPhase = length(measured.phaseVelocities)
    G_phase = zeros(nPhase, L)
    for k in 1:nPhase
        G_phase[k, :] .= first(channelDerivative(:phase, cVec0, L, identity,
            (symmNo = symmNo, rho = rho, pv = measured.phaseVelocities[k])))
    end

    hasBulk = !isnan(measured.bulkModulus)
    G_bulk = hasBulk ?
        permutedims(first(channelDerivative(:bulk, cVec0, L, identity, (symmNo = symmNo, rho = rho)))) :
        zeros(0, L)

    return vcat(G_ritz, G_lim, G_phase, G_bulk),
           (nRitz = nRitz, nLim = nLim, nPhase = nPhase, hasBulk = hasBulk)
end

# -----------------------------------------------------------------------------------------
# calc
# -----------------------------------------------------------------------------------------

"""
    calc(structure, association, invertedLayerNo, measured;
         pointWeights = nothing, includePseudoSAW = true, useParallel = :auto,
         sensitivityCutoff = DEFAULT_SENSITIVITY_CUTOFF) -> SensitivityAnalysis

Build the sensitivity matrix and its eigendecomposition. `structure`/`association`/
`measured`: scalars (one dataset) or vectors (a JOINT inversion of several cuts of the SAME
crystal: same inverted-layer index, symmetry class and constant count; the per-dataset
blocks are stacked before ONE SVD; row order dataset by dataset, each Ritz / limiting /
phase / bulk). The initial guess is dataset 1's inverted-layer constants.
`pointWeights` (one per row, that order; default unweighted) multiply G's rows;
`includePseudoSAW` (default true) decides whether pseudo-Rayleigh points contribute rows.
`psi`/`T` come from `svd(G)`, never `eig(G'G)` (squaring the condition number loses the
small eigenvalues the cutoff depends on); each column of `T` is sign-pinned (largest entry
positive).
"""
function calc(structure::Union{Structure,AbstractVector{<:Structure}}, association,
              invertedLayerNo::Integer,
              measured::Union{MeasuredVelocities,AbstractVector{<:MeasuredVelocities}};
              pointWeights::Union{Nothing,AbstractVector{<:Real}} = nothing,
              includePseudoSAW::Bool = true,
              useParallel::Symbol = :auto,
              sensitivityCutoff::Real = DEFAULT_SENSITIVITY_CUTOFF)
    if pointWeights !== nothing
        all(w -> isfinite(w) && w >= 0, pointWeights) || throw(ArgumentError(
            "SensitivityAnalysis:calc:badPointWeights -- PointWeights must be finite and " *
            "nonnegative."))
    end
    (isfinite(sensitivityCutoff) && 0 < sensitivityCutoff <= 1) || throw(ArgumentError(
        "SensitivityAnalysis:calc:badSensitivityCutoff -- SensitivityCutoff must be in (0, 1]."))

    structureVec = structure isa AbstractVector ? structure : [structure]
    associationVec = association isa AbstractVector ? association : [association]
    measuredVec = measured isa AbstractVector ? measured : [measured]
    nDatasets, symmNoAll, LAll = checkDatasets(structureVec, associationVec, measuredVec,
                                               invertedLayerNo, "SensitivityAnalysis:calc")
    if any(symmNoAll .!= symmNoAll[1]) || any(LAll .!= LAll[1])
        throw(ArgumentError("SensitivityAnalysis:calc:symmetryMismatch -- All jointly-" *
            "inverted datasets must share the same symmetry class and number of independent " *
            "constants for the inverted layer (it is the same crystal): symmClassNo = " *
            "[$(join(symmNoAll, " "))], nIndepConsts = [$(join(LAll, " "))]."))
    end
    symmNo = symmNoAll[1]
    L = LAll[1]
    cVec0 = Float64.(vec(structureVec[1].mater[invertedLayerNo].elConsts.vector))
    includeSet = activeRitzClasses(includePseudoSAW)

    Gcell = Vector{Matrix{Float64}}(undef, nDatasets)
    counts = Vector{NamedTuple}(undef, nDatasets)
    for d in 1:nDatasets
        Gcell[d], counts[d] = datasetSensitivityRows(structureVec[d], associationVec[d],
            invertedLayerNo, measuredVec[d], cVec0, symmNo, L, includeSet, useParallel)
    end
    G = reduce(vcat, Gcell)
    nData = size(G, 1)

    if pointWeights === nothing
        w = ones(nData)
    else
        w = Float64.(vec(pointWeights))
        length(w) == nData || throw(ArgumentError("SensitivityAnalysis:calc:badWeights -- " *
            "PointWeights has $(length(w)) entries; expected $nData (one per data row)."))
    end
    G = Diagonal(w) * G

    F = svd(G)
    T = copy(F.V)
    psi = 2 .* F.S .^ 2
    for pCol in 1:size(T, 2)               # sign pinning, reproducible run to run
        iMax = argmax(abs.(view(T, :, pCol)))
        T[iMax, pCol] < 0 && (T[:, pCol] .= .-T[:, pCol])
    end
    normRoots = sqrt.(psi ./ psi[1])
    nOpt = count(>=(sensitivityCutoff), normRoots)

    return SensitivityAnalysis(Int(L), Int(symmNo), Int(nDatasets), cVec0, T' * cVec0,
                                G, psi, T, normRoots, Float64(sensitivityCutoff), Int(nOpt),
                                w, includePseudoSAW,
                                sum(c.nRitz for c in counts), sum(c.nLim for c in counts),
                                sum(c.nPhase for c in counts), any(c.hasBulk for c in counts))
end
