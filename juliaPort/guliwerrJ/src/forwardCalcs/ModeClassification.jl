# =========================================================================================
# ModeClassification.jl -- port of guliwerr/forwardCalcs/ModeClassification.m (the Green-primary
# classification, its three classification methods and the CLASSIFICATION_VERSION stamp).
#
# EXPERIMENTAL (as in MATLAB): checked against the measured maps published in Grabec et al.
# 2024 and the UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
#
# EXISTENCE is decided by one score per candidate, computed when the angle was calculated
# and stored on the Calc object (threshold-free):
#   'green'   half-space  GreenScore          greenResult.greenScore        (recognizeFull)
#             layered     GreenScoreLayered   greenResult.greenScoreLayered (recognizeLayered)
#   'refined' half-space  GreenSpotScore      greenResult.refined           (recognizeRefined)
#   'ritz'    both        RitzScore           ritzScore.cands               (recognizeRitzOnlyScore)
# IDENTITY (which Ritz mode represents a Green peak) comes from Ritz: the mode nearest the
# single-direction pole (matchRitzToGreenPeak; a 'refined' spot peak via its nearest pole,
# GreenSpotScore.identify, matchSpotPeak only as fallback).
# The label (Rayleigh / pseudoRayleigh, or layered 'SAW' + aboveSubstrateST) comes last.
#
# Julia shape: `ModeClassification` is one immutable struct holding the settings (threshold,
# velMatchTol, matchScoreThresh, associateLimiting, carryEntries) AND the outputs (mode,
# method, significantModes, ritzScorePeaks, greenScorePeaks + greenAccepted), as the MATLAB
# value class does; `calc(mc, ritz, structure, christoffelResult, greenResult, ritzScore)`
# returns a new one. The MATLAB input structs become NamedTuples with the same field names.
#
# Deliberate representation choices (no numerical effect):
# - `type`, `method`, `grade` are Symbols; MATLAB's '' (limiting entries) is `:none`.
# - `aboveSubstrateST` is `Union{Bool,Missing}` (MATLAB NaN = `missing`).
# - The Green candidates keep their own record type (`GreenScore.Peak` /
#   `GreenScoreLayered.Peak`); MATLAB's per-record `accepted` field is the parallel vector
#   `greenAccepted`. `ritzScorePeaks` is `Vector{RitzScore.Peak}` with `accepted` filled.
# =========================================================================================

"""
    RAYLEIGH_LABEL_TOL

Relative tolerance of the half-space Rayleigh/pseudoRayleigh label: a Green peak below
`minShear * (1 + RAYLEIGH_LABEL_TOL)` is 'Rayleigh'. Equal to the Green damping of the
existence solve (`SawSpectrum.ETA`), which broadens every pole to about `eta * v`.
"""
const RAYLEIGH_LABEL_TOL = 3e-3

"""
    CLASSIFICATION_VERSION

Classification stamp, bumped whenever a classification rule or a stored candidate field
changes. A Calc result records the value it was classified with (`classificationVersion`);
another value makes it stale (it cannot be reclassified; recompute it).
"""
const CLASSIFICATION_VERSION = 3

"""
    GREEN_TIE_REL

Relative tie window of [`matchRitzToGreenPeak`](@ref): Ritz modes within
`GREEN_TIE_REL * vGreen` of the nearest mode's distance count as equally near, and the
identity weight decides between them (the 'green' peak's own placement error against a fine
G13_eff peak: p95 0.4-5 m/s).
"""
const GREEN_TIE_REL = 5e-4

"""
    LIMIT_ORTH_TOL

Orthogonality tolerance of [`degenerateLimitGroups`](@ref) (|cos| between two limiting
polarisations): `findRoots` makes a degenerate group's polarisations orthonormal to
round-off; a mirror twin's are not.
"""
const LIMIT_ORTH_TOL = 1e-6

"""
    SignificantMode

One `significantModes` entry (MATLAB `emptyEntries` field set): `type` (`:Rayleigh`,
`:pseudoRayleigh`, `:limiting` half-space; `:SAW`, `:limiting` layered), `vel` (the velocity
of the method that decided existence: `greenPeakVel` on 'green' (the pole) / 'refined' (the
spot peak), the supporting mode's eigenvelocity on 'ritz', the associated mode's
eigenvelocity on limiting entries), `modeNo`, `modeVel` (the eigenvelocity of Ritz mode
`modeNo`, every entry), `matchScore`, `greenPeakVel` (NaN on
'ritz' / limiting entries), `oopDisplScore`, `dispersionScore`, `limVelAssoc` (limiting entries only),
`aboveSubstrateST` (layered SAW entries only, else `missing`), `relPeakHeight`, `method`
(`:green` / `:refined` / `:ritz`, `:none` on limiting entries), `score` (the existence score
on that method's own scale), `scoreDb = 20 score`, `grade` (`:strong` passes the method's
'low' standard, `:weak` only the active one, `:none` on limiting entries).
"""
struct SignificantMode
    type::Symbol
    vel::Float64
    modeNo::Int
    modeVel::Float64
    matchScore::Float64
    greenPeakVel::Float64
    oopDisplScore::Float64
    dispersionScore::Float64
    limVelAssoc::Float64
    aboveSubstrateST::Union{Bool,Missing}
    relPeakHeight::Float64
    method::Symbol
    score::Float64
    scoreDb::Float64
    grade::Symbol
end

"""
    ModeClassification(; threshold = :high, velMatchTol = 0.02, matchScoreThresh = 0.05,
                         associateLimiting = true, carryEntries = SignificantMode[])

Settings of a classification (unclassified: empty outputs). Run it with [`calc`](@ref).
`threshold` is `:high` | `:low` (Strings accepted) | a number on the method's own scale.
"""
struct ModeClassification
    mode::Symbol
    method::Symbol
    significantModes::Vector{SignificantMode}
    velMatchTol::Float64
    matchScoreThresh::Float64
    threshold::Union{Symbol,Float64}
    associateLimiting::Bool
    carryEntries::Vector{SignificantMode}
    ritzScorePeaks::Vector{RitzScore.Peak}
    greenScorePeaks::Union{Vector{GreenScore.Peak},Vector{GreenScoreLayered.Peak}}
    greenAccepted::Vector{Bool}
end

_thresholdArg(t::Real) = Float64(t)
_thresholdArg(t::Union{Symbol,AbstractString}) = Symbol(lowercase(String(t)))

function ModeClassification(; threshold = :high, velMatchTol::Real = 0.02,
                            matchScoreThresh::Real = 0.05, associateLimiting::Bool = true,
                            carryEntries::AbstractVector{SignificantMode} = SignificantMode[])
    return ModeClassification(:full, :green, SignificantMode[], Float64(velMatchTol),
                              Float64(matchScoreThresh), _thresholdArg(threshold),
                              associateLimiting, collect(carryEntries), RitzScore.Peak[],
                              GreenScore.Peak[], Bool[])
end

_withOutputs(mc::ModeClassification; mode, method, sm, ritzPeaks = RitzScore.Peak[],
             greenPeaks = GreenScore.Peak[], greenAcc = Bool[]) =
    ModeClassification(mode, method, sm, mc.velMatchTol, mc.matchScoreThresh, mc.threshold,
                       mc.associateLimiting, mc.carryEntries, ritzPeaks, greenPeaks, greenAcc)

# ---- entry builders (MATLAB makeSAWEntry / makeLimitingEntry) -----------------------------

# vel: greenVel on the Green paths (the 'green' pole, the 'refined' spot peak), the mode's
# eigenvelocity on the Ritz-only path (greenVel NaN); modeVel always the mode's eigenvelocity.
function makeSAWEntry(type::Symbol, ritz, modeNo::Integer, matchScore::Real, greenVel::Real,
                      relPeakHeight::Real, aboveSubstrateST, score::Real, method::Symbol,
                      strong::Bool)
    vMode = ritz.vels[modeNo]
    return SignificantMode(type, isfinite(greenVel) ? Float64(greenVel) : vMode, modeNo, vMode,
                           matchScore, greenVel,
                           ritz.diagnostics.oopSurfDispl[modeNo],
                           ritz.diagnostics.depthInsensitivity[modeNo], NaN,
                           aboveSubstrateST, relPeakHeight, method, score, 20 * score,
                           strong ? :strong : :weak)
end

makeLimitingEntry(ritz, modeNo::Integer, vLim::Real) =
    SignificantMode(:limiting, ritz.vels[modeNo], modeNo, ritz.vels[modeNo], NaN, NaN, NaN,
                    ritz.diagnostics.depthInsensitivity[modeNo], vLim, missing, NaN, :none,
                    NaN, NaN, :none)

# ---- calc ---------------------------------------------------------------------------------

"""
    calc(mc::ModeClassification, ritz::RitzSolverSAW, structure = nothing,
         christoffelResult = nothing, greenResult = nothing, ritzScore = nothing)
        -> ModeClassification

Classify significant modes (port of `ModeClassification.calc`). `ritz` is a calc()'d
RitzSolverSAW (duck-typed: `vels`, `diagnostics`, `depthInsensitivity`, and for the limiting
association `calcDisplacement`, so plain stand-ins work). Dispatch:
- no `christoffelResult` / `greenResult` -> the Ritz-only path; `ritzScore` must be a
  NamedTuple with `cands` (`RitzScore.direction` / `directionLayered` output), else an error;
- `greenResult.refined` -> recognizeRefined (half-space, `peaks` after `GreenSpotScore.identify`);
- `structure` layered -> recognizeLayered (`greenResult.greenScoreLayered.peaks`,
  `christoffelResult` with `subLongLimVel, subLongLimPol, subShearLimVels, subShearLimPols`);
- otherwise recognizeFull (`greenResult.greenScore.peaks`, `christoffelResult` with
  `longLimVel, longLimPol, shearLimVels, shearLimPols`).
"""
function calc(mc::ModeClassification, ritz, structure = nothing,
              christoffelResult = nothing, greenResult = nothing, ritzScore = nothing)
    isLayered = structure !== nothing && nlay(structure) > 1
    if christoffelResult === nothing || greenResult === nothing
        ritzScore === nothing && throw(ArgumentError(
            "ModeClassification:ritzOnly:noScore -- The Ritz-only path needs the ritzScore input."))
        hasproperty(ritzScore, :cands) || throw(ArgumentError(
            "ModeClassification:ritzOnly:noCandidates -- ritzScore needs the field cands " *
            "(RitzScore.direction / directionLayered output)."))
        sm, peaks = recognizeRitzOnlyScore(mc, ritz, ritzScore, isLayered)
        return _withOutputs(mc; mode = :ritzOnly, method = :ritz, sm = sm, ritzPeaks = peaks)
    elseif hasproperty(greenResult, :refined)
        sm = recognizeRefined(mc, ritz, christoffelResult, greenResult.refined)
        return _withOutputs(mc; mode = :full, method = :refined, sm = sm)
    elseif isLayered
        sm, pk, acc = recognizeLayered(mc, ritz, christoffelResult, greenResult)
        return _withOutputs(mc; mode = :layered, method = :green, sm = sm, greenPeaks = pk,
                            greenAcc = acc)
    else
        sm, pk, acc = recognizeFull(mc, ritz, christoffelResult, greenResult)
        return _withOutputs(mc; mode = :full, method = :green, sm = sm, greenPeaks = pk,
                            greenAcc = acc)
    end
end

# ---- Green-peak matching -------------------------------------------------------------------

"""
    matchGreenPeaks(mc, ritz, pk, acc, strong, isLayered) -> Vector{NamedTuple}

The accepted Green candidates in candidate order, each matched to its best Ritz mode:
`(vel, modeNo, score, confident, greenPeakVel, peakHeight, peakProminence, peakWidth,
relPeakHeight, existScore, strong)`.
"""
function matchGreenPeaks(mc::ModeClassification, ritz, pk::AbstractVector, acc, strong,
                         isLayered::Bool)
    out = NamedTuple{(:vel, :modeNo, :score, :confident, :greenPeakVel, :peakHeight,
                      :peakProminence, :peakWidth, :relPeakHeight, :existScore, :strong),
                     Tuple{Float64,Int,Float64,Bool,Float64,Float64,Float64,Float64,Float64,
                           Float64,Bool}}[]
    for j in eachindex(pk)
        acc[j] || continue
        p = pk[j]
        modeNo, score, confident = matchRitzToGreenPeak(ritz, p.vel, mc, isLayered)
        push!(out, (vel = ritz.vels[modeNo], modeNo = modeNo, score = score,
                    confident = confident, greenPeakVel = p.vel, peakHeight = p.height,
                    peakProminence = p.prominence, peakWidth = p.relWidth * p.vel,
                    relPeakHeight = p.relHeight, existScore = p.score, strong = strong[j]))
    end
    return out
end

"""
    matchRitzToGreenPeak(ritz, vGreen, mc, isLayered = false) -> (modeNo, score, confident)

The Ritz mode nearest a single-direction Green peak (the 'green' rule): among the modes within `GREEN_TIE_REL * vGreen` of the nearest distance, the one with the
largest [`identityWeight`](@ref) (first maximum, NaN skipped, as MATLAB `max`). `score` is
[`matchScoreOf`](@ref) that mode (how well it matches, not how it was chosen); `confident =
score > matchScoreThresh`, never a rejection. Always returns a real mode.
"""
function matchRitzToGreenPeak(ritz, vGreen::Real, mc::ModeClassification, isLayered::Bool = false)
    vels = ritz.vels
    dMin = Inf
    @inbounds for n in eachindex(vels)
        d = abs(vels[n] - vGreen)
        d < dMin && (dMin = d)
    end
    lim = dMin + GREEN_TIE_REL * vGreen
    best = -Inf; modeNo = 0; first = 0
    @inbounds for n in eachindex(vels)
        abs(vels[n] - vGreen) <= lim || continue
        first == 0 && (first = n)
        w = _identityWeight1(ritz, n, isLayered)
        if w > best          # strict: the first maximum, as MATLAB max
            best = w; modeNo = n
        end
    end
    modeNo == 0 && (modeNo = first)   # every weight NaN (MATLAB max -> index 1 of near)
    score = _matchScore1(ritz, vGreen, mc, isLayered, modeNo)
    return modeNo, score, score > mc.matchScoreThresh
end

"""
    matchSpotPeak(ritz, vSpot, mc) -> (modeNo, score, confident)

Match a 'refined' (spot-averaged) half-space peak with no single-direction pole in its window
(the fallback of `GreenSpotScore.identify`): argmax over all modes of [`matchScoreOf`](@ref)
(Gaussian velocity proximity x `oopSurfDispl * freqDisp`). Always returns a mode; NaN scores
are skipped (MATLAB `max`).
"""
function matchSpotPeak(ritz, vSpot::Real, mc::ModeClassification)
    best = -Inf; modeNo = 0
    @inbounds for n in eachindex(ritz.vels)
        s = _matchScore1(ritz, vSpot, mc, false, n)
        if s > best          # strict: the first maximum, as MATLAB max
            best = s; modeNo = n
        end
    end
    modeNo == 0 && (modeNo = 1; best = NaN)   # every score NaN (MATLAB max -> NaN, index 1)
    return modeNo, best, best > mc.matchScoreThresh
end

"""
    matchScoreOf(ritz, v, mc, isLayered, modeNos) -> Vector{Float64}

`exp(-(|v_n - v| / v / velMatchTol)^2) * identityWeight` of the given modes:
[`matchSpotPeak`](@ref)'s ranking score and the `matchScore` every 'green' / 'refined' entry
reports.
"""
matchScoreOf(ritz, v::Real, mc::ModeClassification, isLayered::Bool, modeNos) =
    [_matchScore1(ritz, v, mc, isLayered, n) for n in modeNos]

@inline function _identityWeight1(ritz, n::Integer, isLayered::Bool)
    w = ritz.diagnostics.oopSurfDispl[n]
    isLayered || (w *= ritz.diagnostics.freqDisp[n])
    return w
end

@inline function _matchScore1(ritz, v::Real, mc::ModeClassification, isLayered::Bool, n::Integer)
    relDist = abs(ritz.vels[n] - v) / v
    return exp(-(relDist / mc.velMatchTol)^2) * _identityWeight1(ritz, n, isLayered)
end

"""
    identityWeight(ritz, modeNos, isLayered) -> Vector{Float64}

`oopSurfDispl * freqDisp` (half-space) or `oopSurfDispl` (layered) of the given modes.
"""
identityWeight(ritz, modeNos, isLayered::Bool) =
    Float64[_identityWeight1(ritz, n, isLayered) for n in modeNos]

"""
    associateLimitingMode(ritz, vLim, polLim, mc) -> (modeNo::Union{Int,Nothing}, score)

Associate a Christoffel limiting velocity with a Ritz mode within +-5 modes of the nearest:
`sensN * projection * velProximity`, `sensN` the min-max normalised raw
`ritz.depthInsensitivity`, `projection` = |normalised summed |u| profile . polLim|. `nothing`
if the best score does not exceed `matchScoreThresh`.
"""
function associateLimitingMode(ritz, vLim::Real, polLim::AbstractVector{<:Real},
                               mc::ModeClassification, U = nothing)
    di = ritz.depthInsensitivity
    lo, hi = minimum(di), maximum(di)
    vels = ritz.vels
    w1, w2 = _limitingWindow(vels, vLim)
    # U: calcDisplacement of at least modes 1..w2 ([depthPts x n x 3]); appendLimiting
    # computes it once for all limits of the angle, as MATLAB
    (U === nothing || size(U, 2) < w2) && (U = calcDisplacement(ritz, w2))
    best = -Inf; kBest = 0
    for (k, m) in enumerate(w1:w2)
        s1 = 0.0; s2 = 0.0; s3 = 0.0
        @inbounds for i in axes(U, 1)
            s1 += abs(U[i, m, 1]); s2 += abs(U[i, m, 2]); s3 += abs(U[i, m, 3])
        end
        tot = max(s1 + s2 + s3, eps())
        proj = abs((s1 * polLim[1] + s2 * polLim[2] + s3 * polLim[3]) / tot)
        sensN = (di[m] - lo) / max(hi - lo, eps())
        relDist = abs(vels[m] - vLim) / vLim
        s = sensN * proj * exp(-(relDist / mc.velMatchTol)^2)
        if s > best
            best = s; kBest = k
        end
    end
    best <= mc.matchScoreThresh && return nothing, NaN
    return w1 + kBest - 1, best
end

"""
    appendLimiting(mc, ritz, lim, sm) -> Vector{SignificantMode}

Append the limiting-velocity entries (each shear limit, then the longitudinal one) -- or, with
`associateLimiting = false` (a reclassify), the carried `carryEntries` unchanged. `lim` has
`longLimVel, longLimPol, shearLimVels, shearLimPols`.
"""
function _limitingWindow(vels, vLim)
    closestIdx = argmin(abs.(vels .- vLim))
    return max(1, closestIdx - 5), min(length(vels), closestIdx + 5)
end

function appendLimiting(mc::ModeClassification, ritz, lim, sm::Vector{SignificantMode})
    if !mc.associateLimiting
        return vcat(sm, mc.carryEntries)
    end
    vLims = vcat(lim.shearLimVels, lim.longLimVel)
    wMax = maximum(_limitingWindow(ritz.vels, v)[2] for v in vLims)
    U = calcDisplacement(ritz, wMax)                  # one displacement solve per angle
    pols = hcat(reshape(lim.shearLimPols, 3, :), vec(lim.longLimPol))
    # each shear limit, then the longitudinal one; a degenerate group (one velocity,
    # orthogonal polarisations: the shear pair at an acoustic axis) jointly
    grpId = degenerateLimitGroups(vLims, pols)
    modeOf = Vector{Union{Int,Nothing}}(nothing, length(vLims))
    for g in unique(grpId)
        members = findall(==(g), grpId)
        if length(members) == 1
            m = members[1]
            modeOf[m], _ = associateLimitingMode(ritz, vLims[m], view(pols, :, m), mc, U)
        else
            modeNos, _ = associateLimitingGroup(ritz, vLims[members], pols[:, members], mc, U)
            for (k, n) in enumerate(modeNos)
                modeOf[members[k]] = n
            end
        end
    end
    for i in eachindex(vLims)
        modeOf[i] === nothing || push!(sm, makeLimitingEntry(ritz, modeOf[i], vLims[i]))
    end
    return sm
end

"""
    degenerateLimitGroups(vLims, pols) -> Vector{Int}

Group index of each limiting velocity of one angle (shear limits, then `longLimVel`; `pols`
`[3 x n]`); a group of one is a non-degenerate limit. Same group: velocities within
`DEGENERATE_REL_TOL` (relative) AND polarisations orthogonal to [`LIMIT_ORTH_TOL`](@ref).
Velocity alone is not enough: a mirror-symmetric sagittal plane gives every root a twin at
-theta* with the same velocity and a mirrored, NOT orthogonal, polarisation (one physical
root, associated on its own). `findRoots` gives a degenerate group an orthonormal canonical
basis, so the test separates the two. MATLAB `ModeClassification.degenerateLimitGroups`.
"""
function degenerateLimitGroups(vLims::AbstractVector{<:Real}, pols::AbstractMatrix{<:Real})
    n = length(vLims)
    grpId = zeros(Int, n)
    nrm = [norm(view(pols, :, j)) for j in 1:n]
    g = 0
    for a in 1:n
        grpId[a] > 0 && continue
        g += 1
        grpId[a] = g
        members = [a]
        for b in a+1:n
            if grpId[b] == 0 &&
               abs(vLims[b] - vLims[a]) <= DEGENERATE_REL_TOL * abs(vLims[a]) &&
               all(abs(dot(view(pols, :, m), view(pols, :, b))) <= LIMIT_ORTH_TOL * nrm[m] * nrm[b]
                   for m in members)
                grpId[b] = g
                push!(members, b)
            end
        end
    end
    return grpId
end

"""
    associateLimitingGroup(ritz, vGroup, polGroup, mc, U = nothing) -> (modeNos, scores)

Associate a degenerate group of limiting velocities with DISTINCT Ritz modes (MATLAB
`ModeClassification.associateLimitingGroup`): [`associateLimitingMode`](@ref)'s score with
the projection onto the group's span, `||Q' displNorm||` for an orthonormal basis `Q` of
`polGroup` (independent of which basis of the eigenspace the root finder returned), window
and velocity proximity from `vGroup[1]`. The `g` best modes above `matchScoreThresh` (ties to
the lower mode), returned in ascending mode order: member `k` takes `modeNos[k]`.
"""
function associateLimitingGroup(ritz, vGroup::AbstractVector{<:Real},
                                polGroup::AbstractMatrix{<:Real}, mc::ModeClassification,
                                U = nothing)
    di = ritz.depthInsensitivity
    lo, hi = minimum(di), maximum(di)
    vels = ritz.vels
    vLim = vGroup[1]
    w1, w2 = _limitingWindow(vels, vLim)
    (U === nothing || size(U, 2) < w2) && (U = calcDisplacement(ritz, w2))
    Q = Matrix(qr(polGroup).Q)[:, 1:size(polGroup, 2)]
    window = collect(w1:w2)
    s = zeros(length(window))
    for (k, m) in enumerate(window)
        d = zeros(3)
        @inbounds for c in 1:3, i in axes(U, 1)
            d[c] += abs(U[i, m, c])
        end
        d ./= max(sum(d), eps())
        proj = norm(Q' * d)
        sensN = (di[m] - lo) / max(hi - lo, eps())
        relDist = abs(vels[m] - vLim) / vLim
        s[k] = sensN * proj * exp(-(relDist / mc.velMatchTol)^2)
    end
    ord = sortperm(s; rev = true, alg = Base.Sort.DEFAULT_STABLE)   # ties: lower mode first
    ord = filter(k -> s[k] > mc.matchScoreThresh, ord)
    ord = ord[1:min(length(vGroup), length(ord))]
    iSort = sortperm(window[ord])
    return window[ord][iSort], s[ord][iSort]
end

# ---- the four paths ------------------------------------------------------------------------

"""
    recognizeFull(mc, ritz, christoffelResult, greenResult) -> (sm, peaks, accepted)

Half-space Green path: every GreenScore candidate clearing the threshold, bound to its
nearest Ritz mode ([`matchRitzToGreenPeak`](@ref)), labelled Rayleigh below `minShear (1 + RAYLEIGH_LABEL_TOL)` (by the Green peak
velocity), else pseudoRayleigh; then the limiting entries.
"""
function recognizeFull(mc::ModeClassification, ritz, cr, greenResult)
    hasproperty(greenResult, :greenScore) || throw(ArgumentError(
        "ModeClassification:recognizeFull:noScore -- greenResult carries no greenScore input: " *
        "the Green path needs the scored candidates."))
    pk = greenResult.greenScore.peaks
    acc = GreenScore.accepted(pk, mc.threshold)
    strong = GreenScore.accepted(pk, :low)
    matched = matchGreenPeaks(mc, ritz, pk, acc, strong, false)
    vLabel = minimum(cr.shearLimVels) * (1 + RAYLEIGH_LABEL_TOL)
    sm = SignificantMode[]
    for m in matched
        type = m.greenPeakVel < vLabel ? :Rayleigh : :pseudoRayleigh
        push!(sm, makeSAWEntry(type, ritz, m.modeNo, m.score, m.greenPeakVel, m.relPeakHeight,
                               missing, m.existScore, :green, m.strong))
    end
    sm = appendLimiting(mc, ritz, cr, sm)
    return sm, collect(pk), acc
end

"""
    recognizeRefined(mc, ritz, christoffelResult, rf) -> Vector{SignificantMode}

Refined path: one entry per accepted GreenSpotScore candidate (identity and label already
fixed by `GreenSpotScore.identify`), `relPeakHeight = height / max height`; then the limiting
entries.
"""
function recognizeRefined(mc::ModeClassification, ritz, cr, rf)
    sm = SignificantMode[]
    pk = rf.peaks
    if !isempty(pk)
        acc = GreenSpotScore.accepted(pk, mc.threshold)
        strong = GreenSpotScore.accepted(pk, :low)
        hMax = maximum(p.height for p in pk)
        for j in eachindex(pk)
            acc[j] || continue
            p = pk[j]
            push!(sm, makeSAWEntry(Symbol(p.type), ritz, p.modeNo, p.matchScore, p.vel,
                                   p.height / hMax, missing, p.scoreHigh, :refined, strong[j]))
        end
    end
    return appendLimiting(mc, ritz, cr, sm)
end

"""
    recognizeLayered(mc, ritz, subChristoffelResult, greenResult) -> (sm, peaks, accepted)

Layered Green path: every GreenScoreLayered candidate clearing the threshold, matched with
the 'green' rule (nearest mode; layered weight for ties and matchScore), sorted by the Green
peak velocity (= the entries' `vel`), one 'SAW' type with
`aboveSubstrateST = greenPeakVel > min(subShearLimVels)`; then the substrate's limiting
entries.
"""
function recognizeLayered(mc::ModeClassification, ritz, sc, greenResult)
    hasproperty(greenResult, :greenScoreLayered) || throw(ArgumentError(
        "ModeClassification:recognizeLayered:noScore -- greenResult carries no " *
        "greenScoreLayered input: the layered Green path needs the scored candidates."))
    pk = greenResult.greenScoreLayered.peaks
    acc = GreenScoreLayered.accepted(pk, mc.threshold)
    strong = GreenScoreLayered.accepted(pk, :low)
    matched = matchGreenPeaks(mc, ritz, pk, acc, strong, true)
    matched = matched[sortperm([m.greenPeakVel for m in matched])]  # = the entries' vel; stable, as MATLAB sort
    minSubShear = minimum(sc.subShearLimVels)
    sm = SignificantMode[]
    for m in matched
        push!(sm, makeSAWEntry(:SAW, ritz, m.modeNo, m.score, m.greenPeakVel, m.relPeakHeight,
                               m.greenPeakVel > minSubShear, m.existScore, :green, m.strong))
    end
    lim = (longLimVel = sc.subLongLimVel, longLimPol = sc.subLongLimPol,
           shearLimVels = sc.subShearLimVels, shearLimPols = sc.subShearLimPols)
    sm = appendLimiting(mc, ritz, lim, sm)
    return sm, collect(pk), acc
end

"""
    recognizeRitzOnlyScore(mc, ritz, rs, isLayered) -> (sm, peaks)

Ritz-only path: the stored RitzScore candidates (`rs.cands`), thresholded
(`RitzScore.acceptedHS` / `acceptedLayered`) and labelled. Half-space: Rayleigh if
`belowShear`, else pseudoRayleigh, `matchScore = freqDisp`, `score = dScore`. Layered: 'SAW',
`aboveSubstrateST = !belowShear`, `matchScore = share`, `score = dScore + polTerm`. No
limiting entries; `greenPeakVel` NaN.
"""
function recognizeRitzOnlyScore(mc::ModeClassification, ritz, rs, isLayered::Bool)
    cands = rs.cands
    if isLayered
        acc = RitzScore.acceptedLayered(cands, mc.threshold)
        strong = RitzScore.acceptedLayered(cands, :low)
    else
        acc = RitzScore.acceptedHS(cands, mc.threshold)
        strong = RitzScore.acceptedHS(cands, :low)
    end
    peaks = RitzScore.Peak[_with(cands[j]; accepted = acc[j]) for j in eachindex(cands)]
    sm = SignificantMode[]
    for j in eachindex(peaks)
        acc[j] || continue
        p = peaks[j]
        if isLayered
            push!(sm, makeSAWEntry(:SAW, ritz, p.modeNo, p.share, NaN, p.relHeight,
                                   !p.belowShear, p.dScore + p.polTerm, :ritz, strong[j]))
        else
            push!(sm, makeSAWEntry(p.belowShear ? :Rayleigh : :pseudoRayleigh, ritz, p.modeNo,
                                   p.freqDisp, NaN, p.relHeight, missing, p.dScore, :ritz,
                                   strong[j]))
        end
    end
    return sm, peaks
end
