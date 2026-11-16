# =========================================================================================
# GreenSpotScore.jl -- port of guliwerr/forwardCalcs/GreenSpotScore.m: the
# refiner of the half-space 'full' classification ('refined'). G13_eff summed coherently over
# a Gaussian spread of directions (the probe spot), late-time spectra of that sum at gates 4
# (detection) and 2 (position), 'high' against the direction's own un-gated level fU1, 'low'
# against the structure's fT12 reference. See the MATLAB header and
# docs/reference/architecture/forwardCalcs_modeClassification.md.
#
# EXPERIMENTAL (as in MATLAB): checked against the measured maps published in Grabec et al.
# 2024 and the UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
#
# The direction cache (`DirectionCache`) is the MATLAB `dc` struct; `solveDirections` extends
# a cache (returns a new one, the input is not modified). `identify` takes the two Ritz
# matchers as functions (the classification's matchRitzToGreenPeak and matchSpotPeak), so this
# module does not depend on ModeClassification.
# =========================================================================================

module GreenSpotScore

using Statistics: median
using ...CoreTypes: Structure, nlay
using ..ForwardCalcs: rotateStructure, limitingAndBulk, forEachElement, resolveParallel,
                      thresholdValue, surfaceVRef, shearRange, scoreGreenSolve, mround, _with,
                      ScoreMemo, memoGet, memoPut!, memoClear!, memoSize, structureKey
using ..SawSpectrum: SawSpectrum, latticeAxis
using ..GreenScore: GreenScore

const ETA            = 7e-3
const VEL_STEP       = 0.008
const FWHM           = 2.5
const FWHM_TO_SIGMA  = 2.3548
const TRUNC_SIGMA    = 2.5
const DIR_STEP       = 0.5
const GATE_DETECT    = 4
const GATE_POS       = 2
const BAND_LO        = GreenScore.BAND_LO
const BAND_HI        = 1.00
const SEARCH_HI_TOL  = GreenScore.SEARCH_HI_TOL
const THRESHOLD_HIGH = -0.125
const THRESHOLD_LOW  = -0.550
const REF_ANGLES     = collect(0.0:15.0:165.0)
const REF_CACHE_SIZE = 32
const LABEL_TOL      = 7e-3
const POLE_WIN_REL   = 0.012   # identify: a spot peak's single-direction pole lies within
const POLE_WIN_ABS   = 0.030   # max(POLE_WIN_REL v, POLE_WIN_ABS) km/s (= refinedSawAgreement)

"""
    Peak

One refined candidate (MATLAB `GreenSpotScore.emptyPeaks()` field set): `vel` (reported:
gate-2-linked), `velDet` (gate 4), `posLinked`, `prominence`, `relWidth`, `height` (gate 4),
`fU1`, `scoreHigh`, `scoreLow`; identity fields `modeNo` (0 = none), `matchScore`,
`confident`, `type` ("" | "Rayleigh" | "pseudoRayleigh") after `identify`; `accepted`.
"""
struct Peak
    vel::Float64
    velDet::Float64
    posLinked::Bool
    prominence::Float64
    relWidth::Float64
    height::Float64
    fU1::Float64
    scoreHigh::Float64
    scoreLow::Float64
    modeNo::Int
    matchScore::Float64
    confident::Bool
    type::String
    accepted::Bool
end

"""
    Reference

`reference` result: `vRef`, `fT12` (median over REF_ANGLES of each direction's largest gate-4
prominence), `pmaxPerAngle`, `angles`.
"""
struct Reference
    vRef::Float64
    fT12::Float64
    pmaxPerAngle::Vector{Float64}
    angles::Vector{Float64}
end

"""
    DirectionCache

MATLAB `dc`: `key` (Int64 round(dir 1e6), folded to [0, 180e6)), `dir`, `vL` (longLimVel),
`vS` (min shear), `n` (axis length), `G` ([nMax x nDir] complex G13_eff at ETA, NaN-padded),
`v` (the common axis, VEL_MIN + k VEL_STEP).
"""
struct DirectionCache
    key::Vector{Int64}
    dir::Vector{Float64}
    vL::Vector{Float64}
    vS::Vector{Float64}
    n::Vector{Int}
    G::Matrix{ComplexF64}
    v::Vector{Float64}
end
DirectionCache() = DirectionCache(Int64[], Float64[], Float64[], Float64[], Int[],
                                  zeros(ComplexF64, 0, 0), Float64[])

"""
    kernel() -> (o, w)

Direction offsets [deg] and Gaussian weights (summing to 1) of the probe spot.
"""
function kernel()
    sig = FWHM / FWHM_TO_SIGMA
    m = floor(Int, TRUNC_SIGMA * sig / DIR_STEP + 1e-9)
    o = [k * DIR_STEP for k in -m:m]
    w = exp.(-o .^ 2 ./ (2 * sig^2)); w ./= sum(w)
    return o, w
end

"""
    foldDirection(x) -> (xFolded, key)

Directions folded into [0, 180) deg, rounded to 1e-6 deg, and their Int64 keys.
"""
function foldDirection(x::Real)
    key = mod(Int64(mround(x * 1e6)), Int64(180_000_000))
    return Float64(key) / 1e6, key
end
function foldDirection(x::AbstractArray{<:Real})
    keys = [mod(Int64(mround(xi * 1e6)), Int64(180_000_000)) for xi in x]
    return Float64.(keys) ./ 1e6, keys
end

"""
    directionSet(angles) -> Vector{Float64}

Every folded direction the kernels of `angles` need, unique and sorted.
"""
function directionSet(angles::AbstractVector{<:Real})
    o, _ = kernel()
    _, key = foldDirection([a + oi for a in angles, oi in o])
    return Float64.(sort!(unique(vec(key)))) ./ 1e6
end

"""
    solveOne(structureBase, x) -> NamedTuple (G, vL, vS, n)

One direction: its own axis `GreenScore.gridAxis(AXIS_HI_X_VL vL, VEL_STEP)` and the complex
G13_eff on it at ETA.
"""
function solveOne(structureBase::Structure, x::Real)
    sr = rotateStructure(structureBase, x)
    L = limitingAndBulk(sr.mater[nlay(sr)])
    v = GreenScore.gridAxis(SawSpectrum.AXIS_HI_X_VL * L.longLimVel, VEL_STEP)
    g = scoreGreenSolve(sr, v; eta = ETA)
    return (G = Vector{ComplexF64}(g.G_eff), vL = L.longLimVel, vS = minimum(L.shearLimVels),
            n = length(v))
end

"""
    solveDirections(structureBase, D; useParallel = :auto, dc = DirectionCache()) -> DirectionCache

A cache holding every direction of `D` (folded): `dc`'s entries plus the missing ones, solved
(each alone, on its own axis), in parallel over the directions when allowed.
"""
function solveDirections(structureBase::Structure, D::AbstractVector{<:Real};
                         useParallel = :auto, dc::DirectionCache = DirectionCache())
    x, key = foldDirection(collect(Float64, D))
    perm = sortperm(key); key = key[perm]; x = x[perm]
    keepU = [i == 1 || key[i] != key[i-1] for i in eachindex(key)]
    key = key[keepU]; x = x[keepU]
    new = [!(k in dc.key) for k in key]
    x = x[new]; key = key[new]
    isempty(key) && return dc
    useP = resolveParallel(Symbol(useParallel), length(x))
    sol = forEachElement(Any[xi for xi in x], useP, (xi, _) -> solveOne(structureBase, xi))
    nNew = [s.n for s in sol]
    nMax = max(size(dc.G, 1), maximum(nNew))
    nOld = length(dc.key)
    G = fill(complex(NaN, 0.0), nMax, nOld + length(key))
    G[1:size(dc.G, 1), 1:nOld] = dc.G
    for (i, s) in enumerate(sol)
        G[1:s.n, nOld+i] = s.G
    end
    return DirectionCache([dc.key; key], [dc.dir; x], [dc.vL; [s.vL for s in sol]],
                          [dc.vS; [s.vS for s in sol]], [dc.n; nNew], G,
                          latticeAxis(0, nMax - 1, VEL_STEP))
end

"""
    average(dc, angle) -> (X, v, minShearSpot)

The spot-averaged complex G13_eff of one direction (summed term by term in kernel order, on
the shortest kernel axis) and the spot-weighted minShear.
"""
function average(dc::DirectionCache, angle::Real)
    o, w = kernel()
    _, key = foldDirection(angle .+ o)
    col = [findfirst(==(k), dc.key) for k in key]
    any(isnothing, col) && throw(ArgumentError("GreenSpotScore:average:missingDirection -- " *
        "the direction cache lacks $(count(isnothing, col)) of the $(length(col)) kernel " *
        "directions of $angle deg."))
    n = minimum(dc.n[c] for c in col)
    X = zeros(ComplexF64, n)
    for (j, c) in enumerate(col)
        wj = w[j]
        @inbounds for i in 1:n
            X[i] += wj * dc.G[i, c]
        end
    end
    mss = 0.0
    for (j, c) in enumerate(col)
        mss += w[j] * dc.vS[c]
    end
    return X, dc.v[1:n], mss
end

"""
    direction(X, v, vS, vTf, vL, vRef, ref = nothing) -> Vector{Peak}

Scored candidates of one direction from its spot-averaged spectrum: gate-4 peaks below vL,
velocity linked to the gate-2 peak whose half-prominence half-width contains it (nearest),
kept at a reported velocity <= vTf (1 + SEARCH_HI_TOL). `ref` (a `Reference`) fills scoreLow.
"""
function direction(X::AbstractVector{<:Number}, v::AbstractVector{<:Real}, vS::Real, vTf::Real,
                   vL::Real, vRef::Real, ref = nothing)
    xb = Float64[]
    for i in eachindex(v)
        (v[i] >= BAND_LO * vS && v[i] <= BAND_HI * vTf) || continue
        a = abs(X[i]); isnan(a) || push!(xb, a)
    end
    fU1 = isempty(xb) ? NaN : median(xb)
    A4, vO = SawSpectrum.lateTime(X, v, vRef; gatePeriods = GATE_DETECT)
    P4 = SawSpectrum.peaksOf(A4, vO, vL)
    isempty(P4) && return Peak[]
    A2, vO2 = SawSpectrum.lateTime(X, v, vRef; gatePeriods = GATE_POS)
    P2 = SawSpectrum.peaksOf(A2, vO2, vL)
    hw = P2.relWidth .* P2.vel ./ 2
    fT12 = ref === nothing ? NaN : ref.fT12
    out = Peak[]
    for j in eachindex(P4.vel)
        vel = P4.vel[j]; linked = false
        best = Inf
        for k in eachindex(P2.vel)
            d = abs(P4.vel[j] - P2.vel[k])
            if d <= hw[k] && d < best
                best = d; vel = P2.vel[k]; linked = true
            end
        end
        vel <= vTf * (1 + SEARCH_HI_TOL) || continue
        push!(out, Peak(vel, P4.vel[j], linked, P4.prom[j], P4.relWidth[j], P4.height[j], fU1,
                        log10(P4.prom[j] / fU1), log10(P4.prom[j] / fT12), 0, NaN, false, "", false))
    end
    return out
end

"""
    identify(pk, matchGreen, matchSpot, minShearSpot, poleVels = Float64[]) -> Vector{Peak}

Ritz identity and label of every candidate (MATLAB `GreenSpotScore.identify(pk, ritz,
minShearSpot, poleVels)`). Mode: the spot peak is mapped to the
nearest single-direction pole (`poleVels`: the vel of every GreenScore candidate of this angle,
accepted or not) within `max(POLE_WIN_REL v, POLE_WIN_ABS)`, and that pole bound by
`matchGreen(vPole) -> (modeNo, score, confident)` (`ModeClassification.matchRitzToGreenPeak(ritz,
vPole, mc, false)`); with no pole in the window `matchSpot(vSpot)`
(`ModeClassification.matchSpotPeak(ritz, vSpot, mc)`) decides. Label 'Rayleigh' below
`minShearSpot (1 + LABEL_TOL)`, else 'pseudoRayleigh' (from the spot velocity alone).
"""
function identify(pk::AbstractVector{Peak}, matchGreen, matchSpot, minShearSpot::Real,
                  poleVels::AbstractVector{<:Real} = Float64[])
    return Peak[_identifyOne(p, matchGreen, matchSpot, minShearSpot, poleVels) for p in pk]
end
function _identifyOne(p::Peak, matchGreen, matchSpot, minShearSpot::Real, poleVels)
    v = p.vel
    d = Inf; iPole = 0
    for i in eachindex(poleVels)            # first minimum, as MATLAB min
        di = abs(poleVels[i] - v)
        di < d && (d = di; iPole = i)
    end
    if iPole > 0 && d <= max(POLE_WIN_REL * v, POLE_WIN_ABS)
        modeNo, sc, conf = matchGreen(poleVels[iPole])
    else
        modeNo, sc, conf = matchSpot(v)
    end
    t = v < minShearSpot * (1 + LABEL_TOL) ? "Rayleigh" : "pseudoRayleigh"
    return _with(p; modeNo = modeNo, matchScore = sc, confident = conf, type = t)
end

"""
    resolveThreshold(t) -> (value, scale)

`:high` -> scoreHigh, `:low` -> scoreLow.
"""
resolveThreshold(t) = thresholdValue(t, THRESHOLD_HIGH, THRESHOLD_LOW,
                                     "GreenSpotScore:resolveThreshold:badName")

"""
    accepted(pk, threshold) -> Vector{Bool}

'low' is NOT a subset of 'high' (different normalisers). NaN on the compared scale errors.
"""
function accepted(pk::AbstractVector{Peak}, threshold)
    value, scale = resolveThreshold(threshold)
    isempty(pk) && return Bool[]
    x = scale === :low ? [p.scoreLow for p in pk] : [p.scoreHigh for p in pk]
    any(isnan, x) && throw(ArgumentError("GreenSpotScore:accepted:noScore -- a candidate has " *
        "no $scale-scale score (the 'low' scale needs direction() with the reference)."))
    return x .>= value
end

const _REF_MEMO = ScoreMemo{Reference}(REF_CACHE_SIZE)

"""
    referenceFromCache(structure, dc) -> Reference

fT12 from a direction cache that holds every REF_ANGLES kernel direction (no solve).
"""
function referenceFromCache(structure::Structure, dc::DirectionCache)
    vRef = surfaceVRef(structure)
    pmax = fill(NaN, length(REF_ANGLES))
    for (i, a) in enumerate(REF_ANGLES)
        sr = rotateStructure(structure, a)
        L = limitingAndBulk(sr.mater[nlay(sr)])
        X, v, _ = average(dc, a)
        A, vO = SawSpectrum.lateTime(X, v, vRef; gatePeriods = GATE_DETECT)
        P = SawSpectrum.peaksOf(A, vO, L.longLimVel)
        isempty(P.prom) || (pmax[i] = maximum(P.prom))
    end
    fin = filter(isfinite, pmax)
    return Reference(vRef, isempty(fin) ? NaN : median(fin), pmax, copy(REF_ANGLES))
end

"""
    reference(structure; dc = DirectionCache(), useParallel = :auto) -> Reference

The 'low' reference of a structure (memoised per session by content;
`clearReferenceCache()`). Directions already in `dc` are not solved again.
"""
function reference(structure::Structure; dc::DirectionCache = DirectionCache(), useParallel = :auto)
    key = structureKey(structure)
    cached = memoGet(_REF_MEMO, key)
    cached === nothing || return cached
    dc = solveDirections(structure, directionSet(REF_ANGLES); useParallel, dc)
    ref = referenceFromCache(structure, dc)
    memoPut!(_REF_MEMO, key, ref)
    return ref
end

clearReferenceCache() = memoClear!(_REF_MEMO)
referenceCacheSize() = memoSize(_REF_MEMO)

end # module GreenSpotScore
