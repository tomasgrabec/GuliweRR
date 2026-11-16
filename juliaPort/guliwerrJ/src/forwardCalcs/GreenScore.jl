# =========================================================================================
# GreenScore.jl -- port of guliwerr/forwardCalcs/GreenScore.m: the default
# existence score of the half-space 'full' path. Each G13_eff peak is fitted as one damped
# pole; its height is corrected analytically for the late-time gate and the angular-aperture
# dephasing. See the MATLAB header and docs/reference/architecture/forwardCalcs_modeClassification.md.
#
# EXPERIMENTAL (as in MATLAB): checked against the measured maps published in Grabec et al.
# 2024 and the UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
#
# Structure of the port: the pure maths (bandPeaks, poleFit3, fitPeaks, candidates,
# probeSlopeFromSolve, combineSlopes, score, lowScore, dephaseFactor) take arrays and return
# immutable `Peak` records; the three solver-touching steps (probeSlope's two probe solves,
# reference's band solves) go through `scoreGreenSolve` (CalcCommonScores.jl), so the maths is
# golden-tested on MATLAB-dumped solver outputs.
# =========================================================================================

module GreenScore

using LinearAlgebra: lu!, ldiv!, issuccess
using Statistics: median
using SpecialFunctions: erfcx
using ...CoreTypes: Structure, nlay, rcond_matlab
using ..ForwardCalcs: rotateStructure, limitingAndBulk, forEachElement, resolveParallel,
                      thresholdValue, surfaceVRef, shearRange, scoreGreenSolve, mround,
                      _with, ScoreMemo, memoGet, memoPut!, memoClear!, memoSize, structureKey
using ..SawSpectrum: SawSpectrum, latticeAxis, findpeaksProm

const GATE_PERIODS   = SawSpectrum.GATE_PERIODS
const APERTURE_SIGMA = 0.4
const MAX_REL_WIDTH  = SawSpectrum.MAX_REL_WIDTH
const BAND_LO        = 0.80
const BAND_SOLVE_HI  = 1.03
const SEARCH_HI_TOL  = 0.012
const FIT_SPAN       = 0.016
const PROBE_OFFSET   = 0.5
const PROBE_STEP     = 0.016
const PROBE_HALF     = 3
const DEFAULT_STEP   = 0.006
const MAX_STEP       = 0.016
const THRESHOLD_HIGH = 0.16
const THRESHOLD_LOW  = -0.79
const N_REF_ANGLES   = 12
const REF_CACHE_SIZE = 32

"""
    Peak

One scored GreenScore candidate (MATLAB `GreenScore.emptyPeaks()` field set): `vel`
(log-parabola refined |G| peak), `polVel` (fitted pole), `gamma` (pole half-width [km/s],
>= ETA vel/2), `poleOk`, `height`, `prominence`, `base`, `poleHeight`, `relWidth`,
`relHeight`, `fBand`, `rawScore`; `slope`/`slopeFwd`/`slopeBwd` [km/s per deg] after
`probeSlope`; `gateTerm`, `dephaseTerm`, `score`, `logStrength` after `score`; `lowScore`
after `lowScore` (NaN until filled).
"""
struct Peak
    vel::Float64
    polVel::Float64
    gamma::Float64
    poleOk::Bool
    height::Float64
    prominence::Float64
    base::Float64
    poleHeight::Float64
    relWidth::Float64
    relHeight::Float64
    fBand::Float64
    rawScore::Float64
    slope::Float64
    slopeFwd::Float64
    slopeBwd::Float64
    gateTerm::Float64
    dephaseTerm::Float64
    score::Float64
    logStrength::Float64
    lowScore::Float64
end

"""
    Reference

`GreenScore.reference` result: `vRef` [km/s], `logStrengthTop` (median over the directions of
each direction's largest logStrength), `logStrengthPerAngle`, `angles`.
"""
struct Reference
    vRef::Float64
    logStrengthTop::Float64
    logStrengthPerAngle::Vector{Float64}
    angles::Vector{Float64}
end

# ---- thresholds ------------------------------------------------------------------------

"""
    resolveThreshold(t) -> (value, scale)

`:high` scale compares `pk.score`, `:low` compares `pk.lowScore`; a number is on `:high`.
"""
resolveThreshold(t) = thresholdValue(t, THRESHOLD_HIGH, THRESHOLD_LOW,
                                     "GreenScore:resolveThreshold:badName")

"""
    accepted(pk, threshold) -> Vector{Bool}

Which candidates pass `threshold` (`"high"` | `"low"` | number). A NaN on the compared scale
is an error (`GreenScore:accepted:noScore`): `lowScore` is NaN without the reference.
"""
function accepted(pk::AbstractVector{Peak}, threshold)
    value, scale = resolveThreshold(threshold)
    isempty(pk) && return Bool[]
    x = scale === :low ? [p.lowScore for p in pk] : [p.score for p in pk]
    any(isnan, x) && throw(ArgumentError("GreenScore:accepted:noScore -- a candidate has " *
        "no $scale-scale score (the 'low' scale needs direction() with the reference)."))
    return x .>= value
end

const _STEP_WARNED = Ref(false)

"""
    resolveStep(step) -> Float64

A Green axis step the score is valid for: positive and finite, capped at MAX_STEP with a
once-per-session warning (`GreenScore._STEP_WARNED[] = false` re-enables it).
"""
function resolveStep(step::Real)
    (isfinite(step) && step > 0) ||
        throw(ArgumentError("GreenScore.resolveStep: step must be positive and finite; got $step"))
    if step > MAX_STEP + 1e-12
        if !_STEP_WARNED[]
            @warn "Green axis step $step km/s is above the largest step the existence score " *
                  "(GreenScore) is validated for; using $MAX_STEP km/s. (Shown once per session.)"
            _STEP_WARNED[] = true
        end
        return Float64(MAX_STEP)
    end
    return Float64(step)
end

# ---- axes ------------------------------------------------------------------------------

"""
    gridAxis(vMax, step) -> Vector{Float64}

Calc Green axis `VEL_MIN + n*step` up to `vMax + EDGE_TAPER` (lattice, so a short axis is the
leading part of a longer one).
"""
function gridAxis(vMax::Real, step::Real)
    n = floor(Int, (vMax + SawSpectrum.EDGE_TAPER - SawSpectrum.VEL_MIN) / step + 1e-9)
    return latticeAxis(0, n, step)
end

"""
    bandAxis(vS, vTf, step) -> Vector{Float64}

The band `[BAND_LO vS, BAND_SOLVE_HI vTf]` of one direction on the Calc lattice.
"""
function bandAxis(vS::Real, vTf::Real, step::Real)
    v0 = SawSpectrum.VEL_MIN
    n1 = floor(Int, (BAND_LO * vS - v0) / step)
    n2 = ceil(Int, (BAND_SOLVE_HI * vTf - v0) / step)
    return latticeAxis(max(n1, 0), n2, step)
end

# ---- shared with GreenScoreLayered: band peaks, pole fit ---------------------------------

"""
    bandPeaks(G, velAxis, lo, hi, errId) -> (found, v, G, s, h, loc, wid, prom)

The band `[lo, hi]` (+-1.5 steps) of a uniform axis and the local maxima of |G| on it
(findpeaks, half-prominence widths, prominence floor 1e-4 x the largest value), unfiltered.
`found = false` when the band holds fewer than 5 (finite) samples. Errors `errId` if the axis
is not uniform over the band.
"""
function bandPeaks(G::AbstractVector{<:Number}, velAxis::AbstractVector{<:Real},
                   lo::Real, hi::Real, errId::AbstractString)
    E = Float64[]
    vAll = velAxis
    if length(vAll) < 5
        return (false, E, ComplexF64[], E, E, E, E, E)
    end
    dv0 = vAll[2] - vAll[1]
    b = findall(x -> x >= lo - 1.5 * dv0 && x <= hi + 1.5 * dv0, vAll)
    v = Float64.(vAll[b]); Gb = ComplexF64.(G[b]); s = abs.(Gb)
    length(v) < 5 && return (false, v, Gb, s, E, E, E, E)
    dv = v[2] - v[1]
    for i in 1:length(v)-1
        abs((v[i+1] - v[i]) - dv) > 1e-9 &&
            throw(ArgumentError("$errId -- velAxis must be uniform over the band."))
    end
    m = findall(isfinite, s)
    length(m) < 5 && return (false, v, Gb, s, E, E, E, E)
    sm = s[m]
    h, loc, wid, prom = findpeaksProm(sm, v[m], 1e-4 * maximum(sm))
    return (true, v, Gb, s, h, loc, wid, prom)
end

"""
    poleFit3(v, g) -> (pV, pG, hP, r)

3-point pole model `G(x) ~ (a0 + a1 x)/(x - q)`, `x = v - v[2]`, through 3 complex samples:
`pV = Re pole`, `pG = -Im pole` (half-width, > 0 for a damped resonance), `hP = max |model|`
over the span (81 points), `r` the residue `(a0 + a1 q) sc` of the pole. NaN if
ill-conditioned (`rcond(M'M) < 1e-14`) or a sample is non-finite.
"""
function poleFit3(v::AbstractVector{<:Real}, g::AbstractVector{<:Number})
    nanOut = (NaN, NaN, NaN, complex(NaN, NaN))
    (length(v) == 3 && length(g) == 3) || return nanOut
    all(isfinite, g) || return nanOut
    x1 = v[1] - v[2]; x3 = v[3] - v[2]
    sc = max(abs(x1), 0.0, abs(x3))
    xs1 = x1 / sc; xs2 = 0.0 / sc; xs3 = x3 / sc
    M = ComplexF64[1 xs1 g[1]; 1 xs2 g[2]; 1 xs3 g[3]]
    rhs = ComplexF64[g[1] * xs1, g[2] * xs2, g[3] * xs3]
    rcond_matlab(M' * M) < 1e-14 && return nanOut
    F = lu!(copy(M); check = false)
    issuccess(F) || return nanOut
    c = ldiv!(F, rhs)
    q = c[3]
    p = v[2] + q * sc
    pV = real(p); pG = -imag(p)
    hP = 0.0
    for k in 0:80
        xf = -1 + k * (2 / 80)
        hP = max(hP, abs((c[1] + c[2] * xf) / (xf - q)))
    end
    r = (c[1] + c[2] * q) * sc
    return (pV, pG, hP, r)
end

"""
    fitPeaks(v, G, h, loc, wid) -> (vel, polVel, gamma, poleOk, poleHeight)

Pole fit of every grid peak (columns): reported velocity (log-parabola refined |G| peak),
pole, half-width (floored at the damping width ETA vel/2; half the half-prominence width when
the fit failed), poleOk (a damped pole within the FIT_SPAN stencil), poleHeight (the model's
maximum, never below the grid peak; the grid peak when the fit failed).
"""
function fitPeaks(v::AbstractVector{Float64}, G::AbstractVector{ComplexF64},
                  h::AbstractVector{Float64}, loc::AbstractVector{Float64},
                  wid::AbstractVector{Float64})
    n = length(loc); dv = v[2] - v[1]; s = abs.(G)
    hs = max(1, Int(mround(FIT_SPAN / dv)))
    vel = copy(loc); polVel = fill(NaN, n); gamma = fill(NaN, n)
    poleOk = falses(n); poleHeight = copy(h)
    nv = length(v)
    for j in 1:n
        i0 = findfirst(x -> abs(x - loc[j]) < dv / 4, v)
        if i0 > 1 && i0 < nv && s[i0-1] > 0 && s[i0] > 0 && s[i0+1] > 0
            y1 = log(s[i0-1]); y2 = log(s[i0]); y3 = log(s[i0+1])
            den = y1 - 2 * y2 + y3
            if den < 0
                dl = 0.5 * (y1 - y3) / den
                abs(dl) <= 1 && (vel[j] = loc[j] + dl * dv)
            end
        end
        pV = NaN; pG = NaN; hP = NaN
        if i0 > hs && i0 + hs <= nv
            ix = (i0 - hs, i0, i0 + hs)
            pV, pG, hP, _ = poleFit3([v[ix[1]], v[ix[2]], v[ix[3]]], [G[ix[1]], G[ix[2]], G[ix[3]]])
        end
        if !isfinite(hP) || hP < h[j] || abs(pV - loc[j]) > hs * dv
            hP = h[j]
        end
        ok = isfinite(pV) && abs(pV - loc[j]) <= hs * dv && pG > 0
        ok || (hP = h[j])
        gam = pG
        if !ok || !isfinite(gam) || gam <= 0
            gam = wid[j] / 2
        end
        polVel[j] = pV; poleOk[j] = ok; poleHeight[j] = hP
        gamma[j] = max(gam, SawSpectrum.ETA * vel[j] / 2)
    end
    return (vel = vel, polVel = polVel, gamma = gamma, poleOk = poleOk, poleHeight = poleHeight)
end

# ---- the score steps ---------------------------------------------------------------------

"""
    candidates(G, velAxis, vS, vTf) -> Vector{Peak}

Pole-fitted candidates of one direction's complex G13_eff at SawSpectrum.ETA on a uniform
axis covering the band `[BAND_LO vS, BAND_SOLVE_HI vTf]`. Kept: `loc <= vTf (1 +
SEARCH_HI_TOL)` and half-prominence width <= MAX_REL_WIDTH x velocity. fBand = median |G| over
`[BAND_LO vS, vTf]` (the whole band if no finite sample lies there).
"""
function candidates(G::AbstractVector{<:Number}, velAxis::AbstractVector{<:Real},
                    vS::Real, vTf::Real)
    found, v, Gb, s, h, loc, wid, prom = bandPeaks(G, velAxis, BAND_LO * vS,
        BAND_SOLVE_HI * vTf, "GreenScore:candidates:nonUniform")
    found || return Peak[]
    nb = [isfinite(s[i]) && v[i] >= BAND_LO * vS - 1e-9 && v[i] <= vTf + 1e-9 for i in eachindex(v)]
    any(nb) || (nb = isfinite.(s))
    fBand = median(s[nb])
    keep = [loc[k] <= vTf * (1 + SEARCH_HI_TOL) && wid[k] / loc[k] <= MAX_REL_WIDTH for k in eachindex(loc)]
    h = h[keep]; loc = loc[keep]; wid = wid[keep]; prom = prom[keep]
    f = fitPeaks(v, Gb, h, loc, wid)
    hMax = max(isempty(h) ? -Inf : maximum(h), eps())
    out = Vector{Peak}(undef, length(loc))
    for j in eachindex(loc)
        base = h[j] - prom[j]
        out[j] = Peak(f.vel[j], f.polVel[j], f.gamma[j], f.poleOk[j], h[j], prom[j], base,
            f.poleHeight[j], wid[j] / f.vel[j], h[j] / hMax, fBand,
            log10(max(f.poleHeight[j] - base, floatmin()) / fBand),
            NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    end
    return out
end

"""
    probeStencil(pk) -> (V, vu, iu)

The probe stencil of every candidate: `V[k, j] = vp_j + k PROBE_STEP`, k = -PROBE_HALF ..
PROBE_HALF, `vp` the fitted pole (the reported vel where the fit failed); `vu` the unique
points (rounded to 1e-9 km/s, sorted) one probe solve evaluates, `iu` the index of every
`V` entry in `vu`. Works on any candidate type with `polVel`, `poleOk`, `vel` fields.
"""
function probeStencil(pk::AbstractVector)
    n = length(pk)
    vp = [(p.poleOk && isfinite(p.polVel)) ? p.polVel : p.vel for p in pk]
    nk = 2 * PROBE_HALF + 1
    V = Matrix{Float64}(undef, nk, n)
    for j in 1:n, (ik, k) in enumerate(-PROBE_HALF:PROBE_HALF)
        V[ik, j] = vp[j] + k * PROBE_STEP
    end
    r = [mround(x * 1e9) / 1e9 for x in vec(V)]
    vu = sort!(unique(r))
    iu = [searchsortedfirst(vu, x) for x in r]
    return V, vu, iu, vp
end

"""
    probeSlopeFromSolve(V, vp, Gs, offset) -> Vector{Float64}

One side of the probe: slope [km/s per deg] of every candidate from the probe solve's
G13_eff `Gs` on the stencil `V` (`[nk x nPk]`). The probe pole is the interior local maximum
of |G| nearest the stencil centre, refined by the 3-point pole fit (or a log-parabola
vertex); NaN where the stencil has no interior local maximum.
"""
function probeSlopeFromSolve(V::AbstractMatrix{Float64}, vp::AbstractVector{Float64},
                             Gs::AbstractMatrix{<:Number}, offset::Real)
    nk, n = size(V)
    sl = fill(NaN, n)
    centre = PROBE_HALF + 1
    for j in 1:n
        s = [abs(Gs[i, j]) for i in 1:nk]
        best = 0; bestD = typemax(Int)
        for i in 2:nk-1
            if s[i] > s[i-1] && s[i] >= s[i+1]
                d = abs(i - centre)
                if d < bestD
                    bestD = d; best = i
                end
            end
        end
        best == 0 && continue
        i = best
        pV, pG, _, _ = poleFit3([V[i-1, j], V[i, j], V[i+1, j]], [Gs[i-1, j], Gs[i, j], Gs[i+1, j]])
        if isfinite(pV) && abs(pV - V[i, j]) <= PROBE_STEP && pG > 0
            v2 = pV
        else
            y1 = log(s[i-1]); y2 = log(s[i]); y3 = log(s[i+1])
            den = y1 - 2 * y2 + y3; v2 = V[i, j]
            if den < 0
                dl = 0.5 * (y1 - y3) / den
                abs(dl) <= 1 && (v2 = V[i, j] + dl * PROBE_STEP)
            end
        end
        sl[j] = (v2 - vp[j]) / offset
    end
    return sl
end

"""
    probeOneSide(structureBase, angle, pk, offset) -> Vector{Float64}

One probe solve at `angle + offset` (the stencil points of every candidate, shared points
once) and the slope of every candidate from it.
"""
function probeOneSide(structureBase::Structure, angle::Real, pk::AbstractVector, offset::Real)
    V, vu, iu, vp = probeStencil(pk)
    sr = rotateStructure(structureBase, angle + offset)
    g = scoreGreenSolve(sr, vu; eta = SawSpectrum.ETA)
    Gs = reshape(g.G_eff[iu], size(V))
    return probeSlopeFromSolve(V, vp, Gs, offset)
end

"""
    combineSlopes(sF, sB) -> Vector{Float64}

Same sign and within a factor 2 -> mean; otherwise the gentler one; one side NaN -> the
other; both NaN -> NaN. Equal and opposite to 1e-9 relative (a mirror direction, where the
branch is stationary) -> the mean |slope|: "the gentler side" would pick the sign by
round-off. Only |slope| enters the score.
"""
function combineSlopes(sF::AbstractVector{<:Real}, sB::AbstractVector{<:Real})
    s2 = Float64.(sF)
    for i in eachindex(s2)
        f = sF[i]; b = sB[i]
        both = isfinite(f) && isfinite(b)
        agree = both && sign(f) == sign(b) && max(abs(f), abs(b)) <= 2 * min(abs(f), abs(b))
        if agree
            s2[i] = (f + b) / 2
        elseif both && abs(abs(f) - abs(b)) <= 1e-9 * max(abs(f), abs(b))
            s2[i] = (abs(f) + abs(b)) / 2
        elseif both && abs(b) < abs(f)
            s2[i] = b
        elseif !isfinite(f) && isfinite(b)
            s2[i] = b
        end
    end
    return s2
end

"""
    withSlopes(pk, sF, sB) -> Vector

The candidates with `slopeFwd`, `slopeBwd` and `slope = combineSlopes(sF, sB)` filled (any
candidate type with those fields: GreenScore.Peak, GreenScoreLayered.Peak).
"""
function withSlopes(pk::AbstractVector{T}, sF, sB) where {T}
    s2 = combineSlopes(sF, sB)
    return T[_with(pk[j]; slopeFwd = sF[j], slopeBwd = sB[j], slope = s2[j]) for j in eachindex(pk)]
end

"""
    probeSlope(structureBase, angle, pk) -> Vector

Branch slope dv/dtheta of every candidate from the two probes at `angle +- PROBE_OFFSET`
(`structureBase` UNROTATED). Shared with GreenScoreLayered.
"""
function probeSlope(structureBase::Structure, angle::Real, pk::AbstractVector)
    isempty(pk) && return pk
    sF = probeOneSide(structureBase, angle, pk, PROBE_OFFSET)
    sB = probeOneSide(structureBase, angle, pk, -PROBE_OFFSET)
    return withSlopes(pk, sF, sB)
end

"""
    dephaseFactor(gamma, dv, vRef)

Gated height ratio of a pole of half-width `gamma` with a Gaussian velocity spread `dv`
(std) [km/s]: `F = rho sqrt(pi/2) erfcx((rho + q)/sqrt(2)) exp(-q^2/2)`, `rho = gamma/dv`,
`q = 2 pi GATE_PERIODS dv / vRef`; F = 1 for dv = 0 (or non-finite). Scalar or broadcast.
"""
function dephaseFactor(gamma::Real, dv::Real, vRef::Real)
    dv = abs(dv)
    (isfinite(dv) && dv > 0) || return 1.0
    rho = gamma / dv; q = 2 * GATE_PERIODS * pi * dv / vRef
    return rho * sqrt(pi / 2) * erfcx((rho + q) / sqrt(2)) * exp(-q^2 / 2)
end
dephaseFactor(gamma::AbstractArray, dv::AbstractArray, vRef::Real) = dephaseFactor.(gamma, dv, vRef)

"""
    score(pk, vRef) -> Vector{Peak}

Gate and dephasing terms, `score = rawScore + gateTerm + dephaseTerm` and
`logStrength = score + log10(fBand)`. An undetermined slope is set to the stencil bound
PROBE_HALF PROBE_STEP / PROBE_OFFSET.
"""
function score(pk::AbstractVector{Peak}, vRef::Real)
    sLB = PROBE_HALF * PROBE_STEP / PROBE_OFFSET
    out = Vector{Peak}(undef, length(pk))
    for (j, p) in pairs(pk)
        sl = isfinite(p.slope) ? p.slope : sLB
        g = p.gamma
        gateTerm = -2 * GATE_PERIODS * pi * g / vRef / log(10)
        F = dephaseFactor(g, abs(sl) * APERTURE_SIGMA, vRef)
        dT = log10(max(F, floatmin()))
        sc = p.rawScore + gateTerm + dT
        out[j] = _with(p; gateTerm = gateTerm, dephaseTerm = dT, score = sc,
                       logStrength = sc + log10(p.fBand))
    end
    return out
end

"""
    lowScore(pk, ref) -> Vector{Peak}

`lowScore = logStrength - ref.logStrengthTop` (the 'low' scale). `ref` is a `Reference` or
anything with a `logStrengthTop` field.
"""
lowScore(pk::AbstractVector{Peak}, ref) =
    Peak[_with(p; lowScore = p.logStrength - ref.logStrengthTop) for p in pk]

"""
    direction(structureBase, angle, G, velAxis, vS, vTf, vRef, ref = nothing) -> Vector{Peak}

All steps for one direction: candidates -> probeSlope -> score (-> lowScore with `ref`).
`structureBase` UNROTATED; `G` the complex G13_eff at SawSpectrum.ETA of the rotated
structure on a uniform axis covering the band.
"""
function direction(structureBase::Structure, angle::Real, G::AbstractVector{<:Number},
                   velAxis::AbstractVector{<:Real}, vS::Real, vTf::Real, vRef::Real,
                   ref = nothing)
    pk = candidates(G, velAxis, vS, vTf)
    pk = probeSlope(structureBase, angle, pk)
    pk = score(pk, vRef)
    ref === nothing || (pk = lowScore(pk, ref))
    return pk
end

# ---- the 'low' reference ------------------------------------------------------------------

const _REF_MEMO = ScoreMemo{Reference}(REF_CACHE_SIZE)

"""
    strongest(structure, angle, vRef) -> Float64

Largest logStrength of one reference direction (-Inf without candidates).
"""
function strongest(structure::Structure, angle::Real, vRef::Real)
    sr = rotateStructure(structure, angle)
    L = limitingAndBulk(sr.mater[nlay(sr)])
    vS, vTf = shearRange(L)
    v = bandAxis(vS, vTf, DEFAULT_STEP)
    g = scoreGreenSolve(sr, v; eta = SawSpectrum.ETA)
    pk = direction(structure, angle, g.G_eff, v, vS, vTf, vRef)
    return isempty(pk) ? -Inf : max(-Inf, maximum(p.logStrength for p in pk))
end

"""
    reference(structure; nAngles = N_REF_ANGLES, useParallel = :auto) -> Reference

The 'low' reference of a structure: median over `nAngles` equally spaced directions over
[0, 180) of each direction's largest logStrength. Memoised per session by content
(`clearReferenceCache()` empties it).
"""
function reference(structure::Structure; nAngles::Integer = N_REF_ANGLES, useParallel = :auto)
    String(useParallel) in ("auto", "always", "never") || throw(ArgumentError(
        "GreenScore:reference:badParallel -- useParallel must be 'auto', 'always' or 'never'; got '$useParallel'."))
    key = (structureKey(structure), Int(nAngles))
    cached = memoGet(_REF_MEMO, key)
    cached === nothing || return cached
    angles = [(k * 180) / nAngles for k in 0:nAngles-1]
    vRef = surfaceVRef(structure)
    useP = String(useParallel) != "never" && resolveParallel(Symbol(useParallel), nAngles)
    best = forEachElement(Any[a for a in angles], useP, (a, _) -> strongest(structure, a, vRef))
    perAngle = Float64[b for b in best]
    fin = filter(isfinite, perAngle)
    ref = Reference(vRef, isempty(fin) ? NaN : median(fin), perAngle, angles)
    memoPut!(_REF_MEMO, key, ref)
    return ref
end

clearReferenceCache() = memoClear!(_REF_MEMO)
referenceCacheSize() = memoSize(_REF_MEMO)

end # module GreenScore
