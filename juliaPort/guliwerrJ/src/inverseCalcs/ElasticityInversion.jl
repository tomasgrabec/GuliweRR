# src/inverseCalcs/ElasticityInversion.jl
#
# EXPERIMENTAL (as in MATLAB): validated on synthetic data only; the interface may change.
#
# Port of `guliwerr/inverseCalcs/ElasticityInversion.m`: the optimisation
# orchestrator. Owns the Structure(s), the inverted-layer index, the MeasuredVelocities, the
# LIVE ModeAssociation(s), the SensitivityAnalysis, the optimiser configuration and the
# results; `calc!` fits the leading P transformed elastic constants (the rest frozen at the
# initial guess), growing P stage by stage (`stepwise = true`) or in one stage.
#
# The objective (`_cost_and_jacobian!`) runs three passes per dataset:
#   PASS 1  Ritz / pseudo-SAW channel, per unique associated angle: K_con by LINEAR
#           SUPERPOSITION of precomputed blocks (K is linear in the constants; no forward
#           sweep), `solveConstrainedEig`, re-identification by M-weighted MAC (`rebind`),
#           optionally STEERED by a live Green peak for pseudo-SAW points; residual and
#           analytic Jacobian from the Rayleigh quotients at the CURRENT eigenvectors.
#   PASS 2  limiting channel: labelled root set, re-identification (`rebindLimiting`),
#           Jacobian by central difference with anchored local refinement.
#   PASS 3  pulse-echo phase velocities and the bulk-modulus pseudo-datum.
# The rebound association is promoted to the live state ONLY when the pooled cost improved,
# so a rejected trial step cannot corrupt the tracking. A point whose re-identification
# collapsed keeps its ROW with residual and Jacobian zeroed (never deleted).
#
# Real gauge: per angle, when the fixed part and every per-parameter
# derivative are real in the gauge u3 -> i u3, the fit loop solves the REAL symmetric
# problem; the Jacobian quotients use the solve basis `Vs`, `rebind` gets the original basis
# `V = gaugeBack * Vs` (unitary map, quotients agree).
#
# Julia-specific (deliberate divergences):
# - PASS 1 is evaluated angle-parallel (Julia threads, BLAS 1 thread each): every angle
#   solves its own eigenproblem and re-binds its OWN points (`rebind` at one angle reads and
#   writes only the points at that angle), writing into its own slot; the slots are merged
#   and the rows assembled in angle order afterwards -- bitwise what the serial loop gives.
# - The optimiser is a projected Levenberg-Marquardt on (r, J) (`_lm_box`), not MATLAB's
#   `lsqnonlin` trust-region-reflective: same Gauss-Newton family and the same tolerances,
#   native box bounds by projection. Chosen over Fminbox(LBFGS) (first order:
#   several times more objective evaluations).
# - The projections of K are block-wise (only the inverted layer's block of dK/dT is
#   nonzero): round-off-level differences from MATLAB's full N.' * K * N.
#
# Mutable struct: the objective is called back by the optimiser an unknown number of times
# and accumulates state across calls (MATLAB handle class).

using LinearAlgebra: Hermitian, Symmetric, cholesky, issuccess, dot, mul!

using ..Solvers: solveConstrainedEig, gaugeVector, isGaugeReal, buildConstraintMatrix,
                 derivsByDepth, SAWBoundary, GammaMatrix

using ..ForwardCalcs: rotateStructure, limitingAndBulk, GreenScore, greenScoreSpectrum,
                      surfaceVRef, defaultThreshold, RAYLEIGH_LABEL_TOL

export ElasticityInversion, calc!, recalcForwardAtResult, pseudoSAWGreenPeaks

# -----------------------------------------------------------------------------------------
# Per-dataset precompute
# -----------------------------------------------------------------------------------------

"""One unique angle's fixed stiffness and per-transformed-constant derivatives, projected
(`T = Float64` in the real gauge, `ComplexF64` otherwise)."""
struct _AngleBlock{T<:Union{Float64,ComplexF64}}
    angleDeg::Float64
    K_con_fixed::Matrix{T}
    dKdT::Array{T,3}          # [nCon x nCon x L], dK/d(transformed p)
end

"""One dataset's iterate-independent precompute (MATLAB `iterState(d)`)."""
struct _IterState
    uniqueAngles::Vector{Float64}
    M_con::Matrix{Float64}
    gaugeBack::Matrix{ComplexF64}          # N' * (d .* N)
    blocks::Vector{Union{_AngleBlock{Float64},_AngleBlock{ComplexF64}}}
end

gaugeOn(b::_AngleBlock{Float64}) = true
gaugeOn(b::_AngleBlock{ComplexF64}) = false

# -----------------------------------------------------------------------------------------
# The type
# -----------------------------------------------------------------------------------------

"""
    ElasticityInversion

Fields as MATLAB's properties: `structure, invertedLayerNo, nDatasets, measured, association,
sensitivity, optimisedConstantRange, stepwise, anchorPolicy (:fixed | :chained),
rebindVelocityWindow, includePseudoSAW, pointWeights, toleranceLevel, lowerBound,
upperBound, macWarnThreshold, limMatchScoreThresh, limitingRescanCadence,
pseudoSAWRefinement (:none | :green), pseudoSAWGreenCadence, costFunctionValue, result`, and
the protected bookkeeping `iterState, evalCount, promotedEvalCount, nGreenRescans,
warnedKeys, initialStructure, initialAssociation`.
"""
mutable struct ElasticityInversion
    structure::Vector{Structure}
    invertedLayerNo::Int
    nDatasets::Int
    measured::Vector{MeasuredVelocities}
    association::Vector{ModeAssociation}
    sensitivity::SensitivityAnalysis

    optimisedConstantRange::Union{Nothing,Tuple{Int,Int}}
    stepwise::Bool
    anchorPolicy::Symbol
    rebindVelocityWindow::Float64
    includePseudoSAW::Bool
    pointWeights::Vector{Float64}
    toleranceLevel::Float64
    lowerBound::Vector{Float64}
    upperBound::Vector{Float64}

    macWarnThreshold::Float64
    limMatchScoreThresh::Float64
    limitingRescanCadence::Int
    pseudoSAWRefinement::Symbol
    pseudoSAWGreenCadence::Int

    costFunctionValue::Float64
    result::NamedTuple

    iterState::Vector{_IterState}
    evalCount::Int
    promotedEvalCount::Int
    nGreenRescans::Int
    warnedKeys::Set{Tuple{String,Int,Int,Int}}
    initialStructure::Vector{Structure}
    initialAssociation::Vector{ModeAssociation}
end

_normalize_range(r::Nothing) = nothing
_normalize_range(r::Integer) = (Int(r), Int(r))
function _normalize_range(r)
    length(r) == 2 || throw(ArgumentError(
        "ElasticityInversion:calc:badRange -- optimisedConstantRange must be a scalar N or a " *
        "[first last] pair; got a value of length $(length(r))."))
    return (Int(r[1]), Int(r[2]))
end

"""
    ElasticityInversion(structure, measured, association, sensitivity, invertedLayerNo;
        anchorPolicy = :chained, rebindVelocityWindow = DEFAULT_REBIND_WINDOW,
        includePseudoSAW = <the sensitivity's>, pointWeights = <the sensitivity's>,
        toleranceLevel = 2, lowerBound = [], upperBound = [], optimisedConstantRange = nothing,
        stepwise = true, macWarnThreshold = 0.7, limMatchScoreThresh = 0.5,
        limitingRescanCadence = 1, pseudoSAWRefinement = :green, pseudoSAWGreenCadence = 5)

Construct the orchestrator. `structure`/`measured`/`association`: scalars or vectors (one per
dataset of a JOINT inversion); `sensitivity`: one, built from the same dataset arrays.
`includePseudoSAW` and `pointWeights` decide which rows G has and how they are weighted, so
both default to the sensitivity's own and a contradicting `includePseudoSAW` errors
(`includePseudoSAWMismatch`); explicit weights must have one entry per data row of the
sensitivity, in its row order. Bounds act on the transformed constants (one per independent
constant). `optimisedConstantRange`: `nothing` (= `(1, nOptimisedConstants)`), `N` (= `(N, N)`)
or `(first, last)`.
"""
function ElasticityInversion(structure, measured, association,
                               sensitivity::SensitivityAnalysis, invertedLayerNo::Integer;
                               anchorPolicy::Symbol = :chained,
                               rebindVelocityWindow::Real = DEFAULT_REBIND_WINDOW,
                               includePseudoSAW::Union{Nothing,Bool} = nothing,
                               pointWeights::Union{Nothing,AbstractVector{<:Real}} = nothing,
                               toleranceLevel::Real = 2,
                               lowerBound::AbstractVector{<:Real} = Float64[],
                               upperBound::AbstractVector{<:Real} = Float64[],
                               optimisedConstantRange = nothing,
                               stepwise::Bool = true,
                               macWarnThreshold::Real = DEFAULT_TRACK_OPTS.macWarnThreshold,
                               limMatchScoreThresh::Real = DEFAULT_TRACK_OPTS.limMatchScoreThresh,
                               limitingRescanCadence::Integer = 1,
                               pseudoSAWRefinement::Symbol = :green,
                               pseudoSAWGreenCadence::Integer = 5)
    structureVec = structure isa AbstractVector ? Vector{Structure}(collect(structure)) : Structure[structure]
    measuredVec = measured isa AbstractVector ? Vector{MeasuredVelocities}(collect(measured)) :
                                                  MeasuredVelocities[measured]
    associationVec = association isa AbstractVector ?
        Vector{ModeAssociation}(collect(association)) : ModeAssociation[association]

    nDatasets, symmNoAll, LAll = checkDatasets(structureVec, associationVec, measuredVec,
                                               invertedLayerNo, "ElasticityInversion:constructor")
    (all(==(sensitivity.symmClassNo), symmNoAll) && all(==(sensitivity.nIndepConsts), LAll)) ||
        throw(ArgumentError("ElasticityInversion:constructor:symmetryMismatch -- All jointly-" *
            "inverted datasets' inverted layer must share the symmetry class and independent-" *
            "constant count the supplied sensitivity was built with (it is the same crystal): " *
            "symmClassNo = [$(join(symmNoAll, " "))] (sensitivity: $(sensitivity.symmClassNo)), " *
            "nIndepConsts = [$(join(LAll, " "))] (sensitivity: $(sensitivity.nIndepConsts))."))

    pseudoSAWRefinement in (:none, :green) || throw(ArgumentError(
        "ElasticityInversion:constructor:badPseudoSAWRefinement -- PseudoSAWRefinement must " *
        "be :none or :green; got :$pseudoSAWRefinement."))
    anchorPolicy in (:fixed, :chained) || throw(ArgumentError(
        "ElasticityInversion:constructor:badAnchorPolicy -- AnchorPolicy must be :fixed or " *
        ":chained; got :$anchorPolicy."))
    (isfinite(rebindVelocityWindow) && rebindVelocityWindow > 0) || throw(ArgumentError(
        "ElasticityInversion:constructor:badRebindVelocityWindow -- RebindVelocityWindow must " *
        "be finite and positive."))
    isfinite(toleranceLevel) || throw(ArgumentError(
        "ElasticityInversion:constructor:badToleranceLevel -- ToleranceLevel must be finite."))
    L = sensitivity.nIndepConsts
    for (name, b) in (("LowerBound", lowerBound), ("UpperBound", upperBound))
        isempty(b) && continue
        any(isnan, b) && throw(ArgumentError(
            "ElasticityInversion:constructor:badBounds -- $name must not contain NaN."))
        length(b) == L || throw(ArgumentError("ElasticityInversion:constructor:badBounds -- " *
            "$name has $(length(b)) entries; expected $L (one per transformed constant)."))
    end

    if includePseudoSAW === nothing
        includePseudoSAW = sensitivity.includePseudoSAW
    elseif includePseudoSAW != sensitivity.includePseudoSAW
        throw(ArgumentError("ElasticityInversion:constructor:includePseudoSAWMismatch -- " *
            "IncludePseudoSAW = $includePseudoSAW contradicts the sensitivity, which was built " *
            "with IncludePseudoSAW = $(sensitivity.includePseudoSAW): G's rows and the " *
            "residual's rows would describe different datasets."))
    end
    nDataRows = size(sensitivity.sensitivityMatrix, 1)
    if pointWeights === nothing
        w = any(!=(1), sensitivity.pointWeights) ? copy(sensitivity.pointWeights) : Float64[]
    else
        w = Float64.(collect(pointWeights))
        if !isempty(w)
            all(x -> isfinite(x) && x >= 0, w) || throw(ArgumentError(
                "ElasticityInversion:constructor:badPointWeights -- PointWeights must be " *
                "finite and nonnegative."))
            length(w) == nDataRows || throw(ArgumentError(
                "ElasticityInversion:constructor:badWeights -- PointWeights has $(length(w)) " *
                "entries; expected $nDataRows (one per data row of the sensitivity)."))
        end
    end
    (0 <= macWarnThreshold <= 1) || throw(ArgumentError(
        "ElasticityInversion:constructor:badMacWarnThreshold -- MacWarnThreshold must be in [0,1]."))
    (0 <= limMatchScoreThresh <= 1) || throw(ArgumentError(
        "ElasticityInversion:constructor:badLimMatchScoreThresh -- LimMatchScoreThresh must " *
        "be in [0,1]."))
    limitingRescanCadence >= 1 || throw(ArgumentError(
        "ElasticityInversion:constructor:badLimitingRescanCadence -- LimitingRescanCadence " *
        "must be an integer >= 1."))
    pseudoSAWGreenCadence >= 1 || throw(ArgumentError(
        "ElasticityInversion:constructor:badPseudoSAWGreenCadence -- PseudoSAWGreenCadence " *
        "must be an integer >= 1."))

    return ElasticityInversion(structureVec, Int(invertedLayerNo), nDatasets, measuredVec,
        associationVec, sensitivity,
        _normalize_range(optimisedConstantRange), stepwise, anchorPolicy,
        Float64(rebindVelocityWindow), includePseudoSAW, w, Float64(toleranceLevel),
        Float64.(collect(lowerBound)), Float64.(collect(upperBound)),
        Float64(macWarnThreshold), Float64(limMatchScoreThresh), Int(limitingRescanCadence),
        pseudoSAWRefinement, Int(pseudoSAWGreenCadence),
        NaN, NamedTuple(),
        _IterState[], 0, 0, 0, Set{Tuple{String,Int,Int,Int}}(),
        copy(structureVec), copy(associationVec))
end

# -----------------------------------------------------------------------------------------
# prepareIteration
# -----------------------------------------------------------------------------------------

# N_b' * B * N_b for a real row block N_b (B real or complex), real BLAS only
_proj(Nb::Matrix{Float64}, B::AbstractMatrix{Float64}) = transpose(Nb) * (B * Nb)
_proj(Nb::Matrix{Float64}, B::AbstractMatrix{ComplexF64}) =
    complex.(transpose(Nb) * (real.(B) * Nb), transpose(Nb) * (imag.(B) * Nb))

"""
    _prepare_angle(structure, invLay, T, N, ranges, rowsLay, d, dKdC_cell, d2KdCdD_cell, angleDeg)

One angle of `prepareIteration`: fixed stiffness (inverted layer's block removed), the
derivatives w.r.t. the transformed constants, the per-angle real-gauge decision (fixed part
AND every derivative real in the gauge), projections into the constrained basis.
"""
function _prepare_angle(structure::Structure, invLay::Int, T::Matrix{Float64},
                        N::Matrix{Float64}, ranges, rowsLay::UnitRange{Int},
                        d::Vector{ComplexF64}, dKdC_cell, d2KdCdD_cell, angleDeg::Float64)
    L = size(T, 1)
    gmA = calc(structure, SAWBoundary(), angleDeg; dKdC_cell = dKdC_cell,
               d2KdCdD_cell = d2KdCdD_cell)
    K_all = gmA.gammaMatAll
    K_full = copy(K_all)
    K_full[rowsLay, rowsLay] .= 0

    dG_str, _ = dummyForLayer(structure, invLay, dKdC_cell, angleDeg)
    nb = length(rowsLay)
    dGT = reshape(reshape(dG_str[rowsLay, rowsLay, :], :, L) * T, nb, nb, L)   # layer block only

    dl = d[rowsLay]
    Kg = conj.(d) .* K_full .* transpose(d)
    dGTg = conj.(dl) .* dGT .* transpose(dl)
    on = isGaugeReal(Kg, norm(K_all))   # |d| = 1: norm of the gauged full K = norm(K)
    for p in 1:L
        on = on && isGaugeReal(view(dGTg, :, :, p))
    end
    Nb = N[rowsLay, :]
    nCon = size(N, 2)
    if on
        Kr = real.(Kg)
        Kc = zeros(nCon, nCon)
        for r in ranges
            Kc .+= _proj(N[r, :], view(Kr, r, r))
        end
        dK = Array{Float64}(undef, nCon, nCon, L)
        for p in 1:L
            dK[:, :, p] = _proj(Nb, real.(view(dGTg, :, :, p)))
        end
        return _AngleBlock{Float64}(angleDeg, Kc, dK)
    else
        Kc = zeros(ComplexF64, nCon, nCon)
        for r in ranges
            Kc .+= _proj(N[r, :], view(K_full, r, r))
        end
        dK = Array{ComplexF64}(undef, nCon, nCon, L)
        for p in 1:L
            dK[:, :, p] = _proj(Nb, view(dGT, :, :, p))
        end
        return _AngleBlock{ComplexF64}(angleDeg, Kc, dK)
    end
end

"""
    _prepare_iteration!(inv)

Per dataset: `M_con` (elasticity- and angle-independent), the gauge back-map, and per unique
associated angle the `_AngleBlock` (threaded over angles). Once per `calc!`.
"""
function _prepare_iteration!(inv::ElasticityInversion)
    T = inv.sensitivity.transformMatrix
    invLay = inv.invertedLayerNo
    states = Vector{_IterState}(undef, inv.nDatasets)
    for d in 1:inv.nDatasets
        s = inv.structure[d]
        N = buildConstraintMatrix(s.calcParams, SAWBoundary())
        dKdC_cell, d2KdCdD_cell = derivsByDepth(s)
        uniqueAngles = uniqueRitzAngles(inv.association[d])
        gm0 = calc(s, SAWBoundary(), 0.0; dKdC_cell = dKdC_cell, d2KdCdD_cell = d2KdCdD_cell)
        M_con = Matrix{Float64}(transpose(N) * (real.(gm0.normMatAll) * N))
        dg = gaugeVector(s.calcParams)
        gaugeBack = Matrix{ComplexF64}(transpose(N) * (dg .* N))
        ranges = _layer_ranges_of(s)
        rowsLay = layerRows(s, invLay)
        nA = length(uniqueAngles)
        blocks = Vector{Union{_AngleBlock{Float64},_AngleBlock{ComplexF64}}}(undef, nA)
        job(ia) = _prepare_angle(s, invLay, T, N, ranges, rowsLay, dg, dKdC_cell,
                                 d2KdCdD_cell, uniqueAngles[ia])
        _foreach_angle(nA) do ia
            blocks[ia] = job(ia)
        end
        states[d] = _IterState(uniqueAngles, M_con, gaugeBack, blocks)
    end
    inv.iterState = states
    return nothing
end

# full-space row range of every layer
function _layer_ranges_of(s::Structure)
    out = UnitRange{Int}[]
    r0 = 0
    for cp in s.calcParams
        push!(out, (r0 + 1):(r0 + cp.basisVecL))
        r0 += cp.basisVecL
    end
    return out
end

"""Run `f(ia)` for `ia = 1:n`, threaded over the angles when Julia has threads and `n > 1`;
every `ia` must write only its own slot. BLAS runs single-threaded in both cases: the
per-angle matrices (nCon of tens to a few hundred) are far too small for threaded BLAS to
pay, and on a loaded machine its threads make the serial path 10-100x slower (measured)."""
function _foreach_angle(f, n::Integer)
    with_blas_threads(1) do
        if Threads.nthreads() > 1 && n > 1
            Threads.@threads :dynamic for ia in 1:n
                with_task_blas_single(() -> f(ia))
            end
        else
            for ia in 1:n
                f(ia)
            end
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------------------
# Tracking options, limiting roots, Green steering
# -----------------------------------------------------------------------------------------

"""This object's re-identification thresholds, as `rebind`/`rebindLimiting` consume them."""
_current_track_opts(inv::ElasticityInversion) =
    merge(DEFAULT_TRACK_OPTS, (macWarnThreshold = inv.macWarnThreshold,
                               limMatchScoreThresh = inv.limMatchScoreThresh,
                               useGreenSteering = inv.pseudoSAWRefinement === :green))

"""
    _track_limiting_roots(inv, assoc, matLay, symmNo, voigtVec, updateAnchor, trackOpts)

Re-identify every tracked limiting root at the current iterate: the full labelled root set
every `limitingRescanCadence`-th evaluation (default: always), anchored local refinement in
between (a failed refinement forces the full scan of that angle).
"""
function _track_limiting_roots(inv::ElasticityInversion, assoc::ModeAssociation, matLay::Material,
                                 symmNo::Integer, voigtVec::AbstractVector{<:Real},
                                 updateAnchor::Bool, trackOpts::NamedTuple)
    limAngles = uniqueLimitingAngles(assoc)
    isempty(limAngles) && return assoc
    fullScanDue = inv.limitingRescanCadence <= 1 ||
                  mod(inv.evalCount - 1, inv.limitingRescanCadence) == 0
    # the full root scans of all angles are independent: computed angle-parallel up front
    rootsAt = Vector{Vector{LimitingRoot}}(undef, length(limAngles))
    if fullScanDue
        _foreach_angle(length(limAngles)) do ia
            rootsAt[ia] = limitingRootsAt(matLay, symmNo, voigtVec, limAngles[ia])
        end
    end
    for (ia, angleDeg) in enumerate(limAngles)
        idxAtAngle = pointsAtAngle(assoc, angleDeg, ("limiting",))
        doFullScan = fullScanDue
        if !doFullScan
            vRef = zeros(length(idxAtAngle))
            thRef = zeros(length(idxAtAngle))
            for (kk, gi) in enumerate(idxAtAngle)
                pt = assoc.points[gi]
                v, th, _, ok = limitingVelocityNear(matLay, symmNo, voigtVec, angleDeg,
                                                     pt.limRootTheta, Int(pt.limBranchIwt))
                if !ok
                    doFullScan = true
                    break
                end
                vRef[kk] = v
                thRef[kk] = th
            end
            if !doFullScan
                for (kk, gi) in enumerate(idxAtAngle)
                    assoc = setLimitingRefined(assoc, gi, vRef[kk], thRef[kk])
                end
            end
        end
        if doFullScan
            roots = isassigned(rootsAt, ia) ? rootsAt[ia] :
                    limitingRootsAt(matLay, symmNo, voigtVec, angleDeg)
            assoc = rebindLimiting(assoc, angleDeg, roots, updateAnchor, trackOpts)
        end
    end
    return assoc
end

"""
    _structure_at_elasticity(inv, d, voigtVec) -> Structure

Dataset `d`'s Structure with the inverted layer at a candidate elasticity (UNROTATED).
"""
function _structure_at_elasticity(inv::ElasticityInversion, d::Integer,
                                    voigtVec::AbstractVector{<:Real})
    elC = ElastConsts(voigtVec, inv.sensitivity.symmClassNo)
    matLayNew = SetElastConsts(inv.structure[d].mater[inv.invertedLayerNo], elC)
    return ModifyLayer(inv.structure[d], inv.invertedLayerNo, matLayNew, nothing)
end

"""
    pseudoSAWGreenPeaks(structureD, angleDeg) -> Vector{GreenScore.Peak}

The pseudo-SAW peaks the forward 'green' classifier reports at one angle: G13_eff at
`SawSpectrum.ETA` on `GreenScore.bandAxis` at the default step, `GreenScore.direction`
(pole-fitted candidates, probe slope, score), `GreenScore.accepted` at the default 'high'
threshold, and the pseudoRayleigh side `vel >= vS (1 + RAYLEIGH_LABEL_TOL)`. `structureD`
is UNROTATED (rotated once here; the solver gets angle 0).
"""
function pseudoSAWGreenPeaks(structureD::Structure, angleDeg::Real)
    structRot = rotateStructure(structureD, angleDeg)
    lim = limitingAndBulk(structRot.mater[end])
    v = GreenScore.bandAxis(lim.vS, lim.vTf, GreenScore.DEFAULT_STEP)
    g = greenScoreSpectrum(structRot, v, "ElasticityInversion.pseudoSAWGreenPeaks")
    pk = GreenScore.direction(structureD, angleDeg, g.G13Eff, g.velAxis, lim.vS, lim.vTf,
                              surfaceVRef(structureD))
    isempty(pk) && return pk
    keep = GreenScore.accepted(pk, defaultThreshold(:green, false)) .&
           ([p.vel for p in pk] .>= lim.vS * (1 + RAYLEIGH_LABEL_TOL))
    return pk[keep]
end

"""
    _refresh_green_peaks_at_angle(assoc, structureD, angleDeg) -> (assoc, nSolves)

One Green solve at ONE angle; store on every pseudo-SAW point bound there the peak nearest
its OWN previous peak (continuity in the peak's sequence, seeded by `anchorGreenPeakVel`),
NaN if there is none. Never by nearness to `calculatedVelocity` (that would mirror the
Ritz tracking it is meant to check).
"""
function _refresh_green_peaks_at_angle(assoc::ModeAssociation, structureD::Structure,
                                        angleDeg::Real)
    idxP = pointsAtAngle(assoc, angleDeg, ("pseudoSAW",))
    isempty(idxP) && return assoc, 0
    peaks = pseudoSAWGreenPeaks(structureD, angleDeg)
    pv = [p.vel for p in peaks]
    for gi in idxP
        if isempty(pv)
            assoc = setGreenPeak(assoc, gi, NaN)
            continue
        end
        p = assoc.points[gi]
        vRef = if isfinite(p.greenPeakVel) && p.greenPeakVel > 0
            p.greenPeakVel
        elseif isfinite(p.anchorGreenPeakVel) && p.anchorGreenPeakVel > 0
            p.anchorGreenPeakVel
        else
            p.calculatedVelocity
        end
        assoc = setGreenPeak(assoc, gi, pv[argmin(abs.(pv .- vRef))])
    end
    return assoc, 1
end

"""
    _maybe_refresh_green_peaks!(inv, voigtVec)

Re-measure every pseudo-SAW point's Green peak at the iterate just PROMOTED: on the first
promotion, then every `pseudoSAWGreenCadence`-th. A rejected trial never burns a solve here.
"""
function _maybe_refresh_green_peaks!(inv::ElasticityInversion, voigtVec::AbstractVector{<:Real})
    inv.pseudoSAWRefinement === :green || return nothing
    due = inv.promotedEvalCount == 1 || mod(inv.promotedEvalCount, inv.pseudoSAWGreenCadence) == 0
    due || return nothing
    for d in 1:inv.nDatasets
        assoc = inv.association[d]
        idxP = findall(p -> p.bindingClass == "pseudoSAW", assoc.points)
        isempty(idxP) && continue
        structureD = _structure_at_elasticity(inv, d, voigtVec)
        anglesP = sort(unique(assoc.points[i].angleDeg for i in idxP))
        res = Vector{Tuple{ModeAssociation,Int}}(undef, length(anglesP))
        # angle-parallel: each angle's refresh touches only its own points
        _foreach_angle(length(anglesP)) do ia
            res[ia] = _refresh_green_peaks_at_angle(assoc, structureD, anglesP[ia])
        end
        points = copy(assoc.points)
        for (ia, a) in enumerate(anglesP)
            for gi in pointsAtAngle(assoc, a, ("pseudoSAW",))
                points[gi] = res[ia][1].points[gi]
            end
            inv.nGreenRescans += res[ia][2]
        end
        inv.association[d] = ModeAssociation(points)
    end
    return nothing
end

_is_green_corroborated(p::AssocPoint, opts::NamedTuple) =
    opts.useGreenSteering && p.bindingClass == "pseudoSAW" &&
    isfinite(p.greenPeakVel) && p.greenPeakVel > 0 &&
    isfinite(p.greenAgreement) && p.greenAgreement <= opts.greenFilterRel

_pinned_suffix(p::AssocPoint) = p.pinned ? " (pinned -- kept in the fit anyway)" : ""

"""
    _emit_tracking_warnings!(inv, P)

Warn ONCE per (dataset, point, stage) about a collapsed re-identification
(`ElasticityInversion:tracking:limRootLost | modeCollision | lowMAC`, the MATLAB id leading
the message). A Green-corroborated pseudo-SAW match is not low-confidence for its MAC.
"""
function _emit_tracking_warnings!(inv::ElasticityInversion, P::Integer)
    opts = _current_track_opts(inv)
    for d in 1:inv.nDatasets
        for (k, p) in enumerate(inv.association[d].points)
            if p.bindingClass == "limiting"
                p.limRootLost || continue
                id = "ElasticityInversion:tracking:limRootLost"
                msg = "Dataset $d, point $k ($(p.waveType) at $(p.angleDeg) deg, stage P = $P): " *
                    "the tracked Christoffel root was lost at this iterate (no root cleared " *
                    "limMatchScoreThresh $(inv.limMatchScoreThresh), last score " *
                    "$(p.limMatchScore), or a finite-difference step lost it); the point's " *
                    "residual and Jacobian row are zeroed$(_pinned_suffix(p))."
            elseif p.modeCollision
                id = "ElasticityInversion:tracking:modeCollision"
                msg = "Dataset $d, point $k ($(p.waveType) at $(p.angleDeg) deg, stage P = $P): " *
                    "the one-to-one assignment gave this point no Ritz mode at all. Its " *
                    "residual and Jacobian row are zeroed$(_pinned_suffix(p))."
            elseif p.macValue < inv.macWarnThreshold && !_is_green_corroborated(p, opts)
                id = "ElasticityInversion:tracking:lowMAC"
                msg = "Dataset $d, point $k ($(p.waveType) at $(p.angleDeg) deg, stage P = $P): " *
                    "MAC $(p.macValue) is below macWarnThreshold ($(inv.macWarnThreshold)). " *
                    "Its residual and Jacobian row are zeroed and its anchor was not " *
                    "updated$(_pinned_suffix(p))."
            else
                continue
            end
            key = (id, d, k, P)
            key in inv.warnedKeys && continue
            push!(inv.warnedKeys, key)
            @warn "$id -- $msg"
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------------------
# The objective
# -----------------------------------------------------------------------------------------

"""Result of one angle of PASS 1."""
struct _AngleResult
    points::Vector{Pair{Int,AssocPoint}}   # the angle's Ritz points after re-binding
    rowIdx::Vector{Int}                    # their positions in the sensitivity's row order
    r::Vector{Float64}
    J::Matrix{Float64}                     # [nRows x P]
    nGreen::Int
end

# K = K_fixed + sum_p chat_p dKdT_p
function _assemble_K(b::_AngleBlock{T}, chatFull::Vector{Float64}) where {T}
    K = copy(b.K_con_fixed)
    n2 = length(K)
    @inbounds for p in eachindex(chatFull)
        c = chatFull[p]
        off = (p - 1) * n2
        @simd for i in 1:n2
            K[i] += c * b.dKdT[off + i]
        end
    end
    return K
end

# Re[alpha' dKdT_p alpha] for p = 1:P
function _quad_forms!(out::AbstractVector{Float64}, b::_AngleBlock{T}, alpha::AbstractVector{T},
                      tmp::Vector{T}) where {T}
    for p in eachindex(out)
        mul!(tmp, view(b.dKdT, :, :, p), alpha)
        out[p] = real(dot(alpha, tmp))
    end
    return out
end

"""
    _pass1_angle(inv, d, b, chatFull, voigtVec, P, scratchAssoc, trackOpts, updateAnchor,
                 ritzRank, rowOffset) -> _AngleResult

One angle of PASS 1: assemble, solve, map back from the gauge, re-bind; if Green steering
just failed at this angle (a pseudo-SAW point with no candidate inside the Green filter:
`modeCollision`), re-measure the peak at the current elasticity and redo the match from
the PRE-rebind state -- the one sanctioned Green solve on the trial path. Then the rows.
"""
function _pass1_angle(inv::ElasticityInversion, d::Int, b::_AngleBlock{T},
                      M_con::Matrix{Float64}, gaugeBack::Matrix{ComplexF64},
                      chatFull::Vector{Float64}, voigtVec::Vector{Float64}, P::Int,
                      scratchAssoc::ModeAssociation, trackOpts::NamedTuple, updateAnchor::Bool,
                      ritzRank::Vector{Int}, rowOffset::Int) where {T}
    angleDeg = b.angleDeg
    K = _assemble_K(b, chatFull)
    eigVals, Vs = solveConstrainedEig(K, M_con)
    V = T === Float64 ? gaugeBack * Vs : Vs
    vels = sqrt.(abs.(eigVals)) ./ (2π)

    assocBefore = scratchAssoc
    assoc = rebind(scratchAssoc, angleDeg, V, M_con, vels, inv.rebindVelocityWindow,
                   updateAnchor, trackOpts)
    nGreen = 0
    if trackOpts.useGreenSteering
        idxPS = pointsAtAngle(assoc, angleDeg, ("pseudoSAW",))
        if any(assoc.points[i].modeCollision for i in idxPS)
            sG = _structure_at_elasticity(inv, d, voigtVec)
            assocBefore, nGreen = _refresh_green_peaks_at_angle(assocBefore, sG, angleDeg)
            assoc = rebind(assocBefore, angleDeg, V, M_con, vels, inv.rebindVelocityWindow,
                           updateAnchor, trackOpts)
        end
    end

    idxPts = pointsAtAngle(assoc, angleDeg, RITZ_CLASSES)
    pts = Pair{Int,AssocPoint}[gi => assoc.points[gi] for gi in idxPts]
    rowIdx = Int[]
    r = Float64[]
    Jrows = Vector{Float64}[]
    q = zeros(P)
    tmp = Vector{T}(undef, size(Vs, 1))
    for gi in idxPts
        pt = assoc.points[gi]
        (pt.bindingClass == "pseudoSAW" && !inv.includePseudoSAW) && continue
        push!(rowIdx, rowOffset + ritzRank[gi])
        if !pt.activeThisEvaluation
            push!(r, 0.0)
            push!(Jrows, zeros(P))
            continue
        end
        push!(r, pt.measuredVelocity - pt.calculatedVelocity)
        _quad_forms!(q, b, view(Vs, :, Int(pt.ritzModeNo)), tmp)
        push!(Jrows, .-ritzVelocityDerivative(q, pt.calculatedVelocity))
    end
    J = isempty(Jrows) ? zeros(0, P) : permutedims(reduce(hcat, Jrows))
    return _AngleResult(pts, rowIdx, r, J, nGreen)
end

"""
    _cost_and_jacobian!(inv, x, P) -> (r, J)

The objective: residuals (measured - calculated) and analytic Jacobian, POOLED over every
dataset, rows weighted by `pointWeights` in the sensitivity's row order. Promotes every
dataset's rebound association TOGETHER, updates `costFunctionValue`, emits the tracking
warnings and runs the cadenced Green rescan, ONLY when the pooled cost improved on the best
of the current stage.
"""
function _cost_and_jacobian!(inv::ElasticityInversion, x::AbstractVector{<:Real}, P::Integer)
    sens = inv.sensitivity
    T = sens.transformMatrix
    chatFull = copy(sens.initialTransformedVec)
    chatFull[1:P] .= x
    voigtVec = T * chatFull

    inv.evalCount += 1
    trackOpts = _current_track_opts(inv)
    updateAnchor = inv.anchorPolicy === :chained
    activeCls = activeRitzClasses(inv.includePseudoSAW)

    rAll = Float64[]
    Jparts = Matrix{Float64}[]
    wRow = Int[]
    rowOffset = 0
    scratchAssocArr = copy(inv.association)
    nGreenTrial = 0

    for d in 1:inv.nDatasets
        scratchAssoc = scratchAssocArr[d]
        isRitzRow = [p.bindingClass in activeCls for p in scratchAssoc.points]
        ritzRank = cumsum(isRitzRow)
        nRitzRows = count(isRitzRow)

        # ================= PASS 1 (angle-parallel) =================
        st = inv.iterState[d]
        nA = length(st.blocks)
        res = Vector{_AngleResult}(undef, nA)
        base = scratchAssoc
        _foreach_angle(nA) do ia
            res[ia] = _pass1_angle(inv, d, st.blocks[ia], st.M_con, st.gaugeBack, chatFull,
                                   voigtVec, Int(P), base, trackOpts, updateAnchor, ritzRank,
                                   rowOffset)
        end
        points = copy(scratchAssoc.points)
        for ia in 1:nA
            for (gi, pt) in res[ia].points
                points[gi] = pt
            end
            append!(rAll, res[ia].r)
            push!(Jparts, res[ia].J)
            append!(wRow, res[ia].rowIdx)
            nGreenTrial += res[ia].nGreen
        end
        scratchAssoc = ModeAssociation(points)

        # ================= PASS 2: limiting =================
        matLay = inv.structure[d].mater[inv.invertedLayerNo]
        symmNo = sens.symmClassNo
        rho = matLay.rho
        scratchAssoc = _track_limiting_roots(inv, scratchAssoc, matLay, symmNo, voigtVec,
                                               updateAnchor, trackOpts)
        limIdx = findall(p -> p.bindingClass == "limiting", scratchAssoc.points)
        toVoigt = c -> T * c
        # the points' finite-difference rows are independent: computed point-parallel
        fd = Vector{Tuple{Vector{Float64},Bool}}(undef, length(limIdx))
        _foreach_angle(length(limIdx)) do kk
            pt = scratchAssoc.points[limIdx[kk]]
            fd[kk] = pt.activeThisEvaluation ?
                channelDerivative(:limiting, chatFull, P, toVoigt,
                                  (symmNo = symmNo, matLay = matLay, pt = pt)) :
                (zeros(P), true)
        end
        for (kk, gi) in enumerate(limIdx)
            pt = scratchAssoc.points[gi]
            push!(wRow, rowOffset + nRitzRows + kk)
            if !pt.activeThisEvaluation
                push!(rAll, 0.0)
                push!(Jparts, zeros(1, P))
                continue
            end
            row, okFD = fd[kk]
            if okFD
                push!(rAll, pt.measuredVelocity - pt.calculatedVelocity)
                push!(Jparts, -permutedims(row))
            else
                # a perturbed evaluation lost the tracked root: the row is undefined,
                # the point is excluded like a root rebindLimiting found gone
                scratchAssoc = setLimitingLost(scratchAssoc, gi)
                push!(rAll, 0.0)
                push!(Jparts, zeros(1, P))
            end
        end

        # ================= PASS 3: pulse-echo and bulk modulus =================
        meas = inv.measured[d]
        nPhase = length(meas.phaseVelocities)
        for (kk, pv) in enumerate(meas.phaseVelocities)
            push!(wRow, rowOffset + nRitzRows + length(limIdx) + kk)
            vCalc = phaseVelocityAt(voigtVec, symmNo, pv.propagationDir, rho, pv.waveType)
            push!(rAll, pv.velocity - vCalc)
            row, _ = channelDerivative(:phase, chatFull, P, toVoigt,
                                       (symmNo = symmNo, rho = rho, pv = pv))
            push!(Jparts, -permutedims(row))
        end
        hasBulk = !isnan(meas.bulkModulus)
        if hasBulk
            push!(wRow, rowOffset + nRitzRows + length(limIdx) + nPhase + 1)
            push!(rAll, sqrt(meas.bulkModulus / rho) - bulkVelocityAt(voigtVec, symmNo, rho))
            row, _ = channelDerivative(:bulk, chatFull, P, toVoigt, (symmNo = symmNo, rho = rho))
            push!(Jparts, -permutedims(row))
        end

        scratchAssocArr[d] = scratchAssoc
        rowOffset += nRitzRows + length(limIdx) + nPhase + Int(hasBulk)
    end
    inv.nGreenRescans += nGreenTrial

    r = rAll
    J = isempty(Jparts) ? zeros(0, P) : reduce(vcat, Jparts)
    if !isempty(inv.pointWeights)
        length(inv.pointWeights) == length(r) || throw(ArgumentError(
            "ElasticityInversion:costAndJacobian:badWeights -- PointWeights has " *
            "$(length(inv.pointWeights)) entries but the residual has $(length(r)) rows: the " *
            "association differs from the one the sensitivity was built from."))
        w = inv.pointWeights[wRow]
        r = w .* r
        J = w .* J
    end

    cost = sum(abs2, r)
    if isnan(inv.costFunctionValue) || cost < inv.costFunctionValue
        inv.costFunctionValue = cost
        inv.association = scratchAssocArr      # ALL datasets promoted together
        inv.promotedEvalCount += 1
        _emit_tracking_warnings!(inv, P)
        _maybe_refresh_green_peaks!(inv, voigtVec)
    end
    return r, J
end

# -----------------------------------------------------------------------------------------
# Optimiser: projected Levenberg-Marquardt
# -----------------------------------------------------------------------------------------

"""
    _lm_box(fJ, x0, lo, hi; stepTol, funTol, optTol, maxIter = 400)
        -> (x, resnorm, exitflag, info)

Levenberg-Marquardt with Marquardt scaling on `fJ(x) -> (r, J)`, steps projected onto the
box `lo <= x <= hi`. A trial is accepted when it lowers `sum(r.^2)` (then lambda / 3), else
lambda * 4 and the step is recomputed from the same `J` (no new evaluation). Stops (MATLAB
lsqnonlin's tolerance roles) on: first-order optimality `norm(projected J'r, Inf) <=
optTol * max(1, initial)` (exitflag 1), relative cost change of an accepted step `<= funTol`
(3), step `<= stepTol (stepTol + norm(x))` (2), `maxIter` (0). `info`: iterations,
evaluations, firstOrderOpt.
"""
function _lm_box(fJ, x0::Vector{Float64}, lo::Vector{Float64}, hi::Vector{Float64};
                 stepTol::Float64, funTol::Float64, optTol::Float64, maxIter::Int = 400)
    x = clamp.(x0, lo, hi)
    r, J = fJ(x)
    nEval = 1
    cost = sum(abs2, r)
    g = J' * r
    pg(x, g) = [((x[i] <= lo[i] && g[i] > 0) || (x[i] >= hi[i] && g[i] < 0)) ? 0.0 : g[i]
                for i in eachindex(x)]
    opt0 = max(1.0, norm(pg(x, g), Inf))
    lambda = 1e-3
    exitflag = 0
    it = 0
    while it < maxIter
        it += 1
        if norm(pg(x, g), Inf) <= optTol * opt0
            exitflag = 1
            break
        end
        A = J' * J
        dA = [max(A[i, i], 1e-12 * max(1.0, maximum(abs, A))) for i in axes(A, 1)]
        accepted = false
        stepSmall = false
        local xNew, rNew, JNew, costNew
        for _ in 1:30
            M = copy(A)
            for i in axes(M, 1)
                M[i, i] += lambda * dA[i]
            end
            delta = -(Symmetric(M) \ g)
            xNew = clamp.(x .+ delta, lo, hi)
            step = norm(xNew .- x)
            if step <= stepTol * (stepTol + norm(x))
                stepSmall = true
                break
            end
            rNew, JNew = fJ(xNew)
            nEval += 1
            costNew = sum(abs2, rNew)
            if costNew < cost
                accepted = true
                break
            end
            lambda *= 4
        end
        if stepSmall
            exitflag = 2
            break
        end
        accepted || break                     # no descent found at any damping
        relDrop = (cost - costNew) / max(cost, eps())
        x, r, J, cost = xNew, rNew, JNew, costNew
        g = J' * r
        lambda = max(lambda / 3, 1e-12)
        if relDrop <= funTol
            exitflag = 3
            break
        end
    end
    return x, cost, exitflag,
           (iterations = it, funcCount = nEval, firstOrderOpt = norm(pg(x, g), Inf))
end

# -----------------------------------------------------------------------------------------
# calc!
# -----------------------------------------------------------------------------------------

"""
    calc!(inv::ElasticityInversion) -> inv

Run the schedule (step-wise growth of P, or one stage at `range[2]`), populate `inv.result`
and write the recovered elasticity into every dataset's inverted layer. Repeatable: every
call restarts from the constructor's structure and association.
"""
function calc!(inv::ElasticityInversion)
    L = inv.sensitivity.nIndepConsts
    first_, last_ = inv.optimisedConstantRange === nothing ?
        (1, inv.sensitivity.nOptimisedConstants) : inv.optimisedConstantRange
    (1 <= first_ <= last_ <= L) || throw(ArgumentError(
        "ElasticityInversion:calc:badRange -- optimisedConstantRange must be a scalar N or a " *
        "[first last] pair with 1 <= first <= last <= $L (L = $L independent constants); " *
        "got ($first_, $last_)."))

    inv.structure = copy(inv.initialStructure)
    inv.association = copy(inv.initialAssociation)
    _prepare_iteration!(inv)
    inv.costFunctionValue = NaN
    inv.evalCount = 0
    inv.promotedEvalCount = 0
    inv.nGreenRescans = 0

    tol = inv.toleranceLevel
    stepTol, funTol, optTol = 10.0^(-(3 + tol)), 10.0^(-(5 + tol)), 10.0^(-(5 + tol))
    xPrev = Float64[]
    convergenceLog = NamedTuple[]
    Plist = inv.stepwise ? collect(first_:last_) : [last_]
    for P in Plist
        inv.costFunctionValue = NaN     # every stage: no gating against the previous optimum
        empty!(inv.warnedKeys)          # warn once per (point, stage)
        x0 = P == Plist[1] ? inv.sensitivity.initialTransformedVec[1:P] :
             vcat(xPrev, inv.sensitivity.initialTransformedVec[P])
        lo = isempty(inv.lowerBound) ? fill(-Inf, P) : inv.lowerBound[1:P]
        hi = isempty(inv.upperBound) ? fill(Inf, P) : inv.upperBound[1:P]
        xSol, resnorm, exitflag, output = _lm_box(x -> _cost_and_jacobian!(inv, x, P), x0, lo, hi;
                                                  stepTol, funTol, optTol)
        xPrev = xSol
        push!(convergenceLog, (P = P, exitflag = exitflag, output = output, resnorm = resnorm))
    end
    _finalize_result!(inv, xPrev, last_, convergenceLog)
    return inv
end

"""
    recalcForwardAtResult(inv; angles = nothing, vMax = nothing, calcMode = :full,
                          retentionLevel = :full)

Fresh forward sweep at the inverted elasticity (dataset 1), for display/checking; never
called inside `calc!`. `angles` default: dataset 1's associated Ritz angles. `vMax`
(layered only) optionally extends the displayed Green axis.
"""
function recalcForwardAtResult(inv::ElasticityInversion; angles = nothing, vMax = nothing,
                                 calcMode::Symbol = :full, retentionLevel::Symbol = :full)
    angs = angles === nothing ? uniqueRitzAngles(inv.association[1]) : angles
    structure1 = inv.structure[1]
    sweep = if nlay(structure1) == 1
        angSweepHS(; structure = structure1, angles = angs, calcMode, retentionLevel)
    elseif vMax === nothing
        angSweepLayered(; structure = structure1, angles = angs, calcMode, retentionLevel)
    else
        angSweepLayered(; structure = structure1, angles = angs, vMax, calcMode, retentionLevel)
    end
    return calc(sweep)
end

# -----------------------------------------------------------------------------------------
# Result
# -----------------------------------------------------------------------------------------

"""
    _finalize_result!(inv, transformedVecFinal, P, convergenceLog)

Reverse-transform, errors in both spaces (Grabec et al. 2024, Eqs. 13-14:
`delta_p = sqrt(2 Fmin / (N psi_p))` for ALL L combinations, `sigma_i = sqrt(sum_p T_ip^2
delta_p^2)`), the SAME elasticity written into every dataset, `inv.result`.
"""
function _finalize_result!(inv::ElasticityInversion, transformedVecFinal::AbstractVector{<:Real},
                             P::Integer, convergenceLog::Vector{<:NamedTuple})
    sens = inv.sensitivity
    L = sens.nIndepConsts
    T = sens.transformMatrix
    chatFull = copy(sens.initialTransformedVec)
    chatFull[1:P] .= transformedVecFinal
    voigtVec = T * chatFull

    elConsts = ElastConsts(voigtVec, sens.symmClassNo)
    for d in 1:inv.nDatasets
        matFinal = SetElastConsts(inv.structure[d].mater[inv.invertedLayerNo], elConsts)
        inv.structure[d] = ModifyLayer(inv.structure[d], inv.invertedLayerNo, matFinal, nothing)
    end
    Fmin = convergenceLog[end].resnorm

    activeCls = activeRitzClasses(inv.includePseudoSAW)
    nData = 0
    for d in 1:inv.nDatasets
        pts = inv.association[d].points
        nData += count(p -> p.bindingClass in activeCls, pts) +
                 count(p -> p.bindingClass == "limiting", pts) +
                 length(inv.measured[d].phaseVelocities) + Int(!isnan(inv.measured[d].bulkModulus))
    end
    deltaHat = sqrt.(2 .* Fmin ./ (nData .* sens.sensitivityEigenvalues))
    frozen = falses(L)
    frozen[(P + 1):end] .= true
    voigtErr = sqrt.((T .^ 2) * (deltaHat .^ 2))

    inv.result = (transformedElasticVec = chatFull,
                  transformedElasticVecErrors = deltaHat,
                  sensitivityEigenvalues = sens.sensitivityEigenvalues,
                  transformMatrix = T,
                  frozenCombinations = frozen,
                  voigtElasticVec = voigtVec,
                  voigtElasticVecErrors = voigtErr,
                  elConsts = elConsts,
                  costFunctionValue = Fmin,
                  nDataPoints = nData,
                  convergence = convergenceLog,
                  tracking = _build_tracking_report(inv))
    return nothing
end

"""Per-class re-identification health, for `inv.result.tracking`."""
function _build_tracking_report(inv::ElasticityInversion)
    classes = ("ritzMode", "pseudoSAW", "limiting")
    minConf = Dict(c => NaN for c in classes)
    nPerCls = Dict(c => 0 for c in classes)
    nLowCls = Dict(c => 0 for c in classes)
    nCollision = nRootLost = nPinned = nExcluded = 0
    greenAgree = Float64[]
    lowCounts = Vector{Vector{Int}}(undef, inv.nDatasets)
    opts = _current_track_opts(inv)
    for d in 1:inv.nDatasets
        assoc = inv.association[d]
        pts = assoc.points
        if isempty(pts)
            lowCounts[d] = Int[]
            continue
        end
        conf = confidence(assoc)
        lowCounts[d] = [p.lowConfidenceCount for p in pts]
        corroborated = [_is_green_corroborated(p, opts) for p in pts]
        for c in classes
            m = findall(p -> p.bindingClass == c, pts)
            isempty(m) && continue
            nPerCls[c] += length(m)
            cVals = filter(!isnan, conf[m])
            isempty(cVals) || (minConf[c] = isnan(minConf[c]) ? minimum(cVals) :
                                                                min(minConf[c], minimum(cVals)))
            thr = c == "limiting" ? inv.limMatchScoreThresh : inv.macWarnThreshold
            for i in m
                isnan(conf[i]) && continue
                gated = c == "limiting" ? conf[i] < thr : (conf[i] < thr && !corroborated[i])
                gated && (nLowCls[c] += 1)
            end
        end
        nCollision += count(p -> p.modeCollision, pts)
        nRootLost += count(p -> p.limRootLost, pts)
        nPinned += count(p -> p.pinned, pts)
        nExcluded += count(p -> !p.activeThisEvaluation, pts)
        for p in pts
            (p.bindingClass == "pseudoSAW" && isfinite(p.greenAgreement)) &&
                push!(greenAgree, p.greenAgreement)
        end
    end
    return (minConfidence = NamedTuple(Symbol(k) => v for (k, v) in minConf),
            nPointsPerClass = NamedTuple(Symbol(k) => v for (k, v) in nPerCls),
            nBelowThresholdPerClass = NamedTuple(Symbol(k) => v for (k, v) in nLowCls),
            nModeCollisions = nCollision,
            nLimitingRootsLost = nRootLost,
            nPinnedPoints = nPinned,
            nExcludedAtFinalEval = nExcluded,
            lowConfidenceCounts = lowCounts,
            pseudoSAWRefinement = inv.pseudoSAWRefinement,
            nGreenRescans = inv.nGreenRescans,
            greenAgreementMax = isempty(greenAgree) ? NaN : maximum(greenAgree),
            greenAgreementMedian = isempty(greenAgree) ? NaN : median(greenAgree),
            nObjectiveEvaluations = inv.evalCount,
            nPromotedEvaluations = inv.promotedEvalCount)
end
