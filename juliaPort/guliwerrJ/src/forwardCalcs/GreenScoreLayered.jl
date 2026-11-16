# =========================================================================================
# GreenScoreLayered.jl -- port of guliwerr/forwardCalcs/GreenScoreLayered.m:
# the existence score of the layered 'full' path. GreenScore's candidates, pole fit and
# two-sided probe on the all-layer band, an own-period late-time gate, aperture dephasing, a
# direction-own normaliser and the surface-polarisation term K_POL (1 - 2 SH). See the
# MATLAB header and docs/reference/architecture/forwardCalcs_modeClassification.md ("Layered
# stacks").
#
# EXPERIMENTAL (as in MATLAB): checked against the measured maps published in Grabec et al.
# 2024 and the UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
# =========================================================================================

module GreenScoreLayered

using Statistics: median
using ...CoreTypes: Structure, nlay
using ..ForwardCalcs: rotateStructure, limitingAndBulk, thresholdValue, scoreGreenSolve, _with
using ..SawSpectrum: SawSpectrum, latticeAxis
using ..GreenScore: GreenScore

const THRESHOLD_LOW    = 0.70
const THRESHOLD_HIGH   = 0.20
const K_POL            = 0.75
const GATE_OWN_PERIODS = 10
const BAND_LO_X_VSMIN  = 0.3
const BAND_HI_X_VLMAX  = 1.05

"""
    Peak

One layered candidate (MATLAB `GreenScoreLayered.emptyPeaks()` field set): `vel`, `polVel`,
`gamma`, `poleOk`, `height`, `prominence`, `base`, `poleHeight`, `relWidth`, `relHeight`,
`slope`/`slopeFwd`/`slopeBwd`, `SH`, `norm`, `gateTerm`, `dephaseTerm`, `polTerm`,
`logStrength`, `score`.
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
    slope::Float64
    slopeFwd::Float64
    slopeBwd::Float64
    SH::Float64
    norm::Float64
    gateTerm::Float64
    dephaseTerm::Float64
    polTerm::Float64
    logStrength::Float64
    score::Float64
end

"""
    resolveThreshold(t) -> Float64

`'low'` | `'high'` | numeric, all on the one score scale.
"""
resolveThreshold(t) = first(thresholdValue(t, THRESHOLD_HIGH, THRESHOLD_LOW,
                                           "GreenScoreLayered:resolveThreshold:badName"))

"""
    accepted(pk, threshold) -> Vector{Bool}

`score >= threshold`; a NaN score (empty normalising band) is an error
(`GreenScoreLayered:accepted:noScore`).
"""
function accepted(pk::AbstractVector{Peak}, threshold)
    thr = resolveThreshold(threshold)
    isempty(pk) && return Bool[]
    x = [p.score for p in pk]
    any(isnan, x) && throw(ArgumentError("GreenScoreLayered:accepted:noScore -- a candidate " *
        "has no score (normaliser undefined: empty band)."))
    return x .>= thr
end

"""
    bandLimits(structureRotated) -> (vSmin, vLmax)

Slowest shear and fastest longitudinal limiting velocity over ALL layers [km/s] (the band's
edges only, never a score input).
"""
bandLimits(sr::Structure) = bandLimits([limitingAndBulk(sr.mater[k]) for k in 1:nlay(sr)])

"`bandLimits` from every layer's `limitingAndBulk`, already computed (layer order)."
function bandLimits(Ls::AbstractVector)
    vSmin = Inf; vLmax = -Inf
    for L in Ls
        vSmin = min(vSmin, minimum(L.shearLimVels))
        vLmax = max(vLmax, L.longLimVel)
    end
    return vSmin, vLmax
end

"""
    band(vSmin, vLmax) -> (vLo, vHi)
"""
band(vSmin::Real, vLmax::Real) = (BAND_LO_X_VSMIN * vSmin, BAND_HI_X_VLMAX * vLmax)

"""
    gridAxis(vLo, vMax) -> Vector{Float64}

Uniform 4 m/s axis on the lattice VEL_MIN + n VEL_STEP from min(VEL_MIN, the lattice point at
or below vLo) to vMax + EDGE_TAPER.
"""
function gridAxis(vLo::Real, vMax::Real)
    v0 = SawSpectrum.VEL_MIN; dv = SawSpectrum.VEL_STEP
    n1 = min(0, floor(Int, (vLo - v0) / dv + 1e-9))
    n2 = floor(Int, (vMax + SawSpectrum.EDGE_TAPER - v0) / dv + 1e-9)
    return latticeAxis(n1, n2, dv)
end

"""
    candidates(G, velAxis, vLo, vHi) -> Vector{Peak}

GreenScore's peak finder, width filter and pole fit on the whole band [vLo, vHi] (peaks must
be interior to it).
"""
function candidates(G::AbstractVector{<:Number}, velAxis::AbstractVector{<:Real},
                    vLo::Real, vHi::Real)
    found, v, Gb, _, h, loc, wid, prom = GreenScore.bandPeaks(G, velAxis, vLo, vHi,
        "GreenScoreLayered:candidates:nonUniform")
    found || return Peak[]
    keep = [loc[k] >= vLo && loc[k] <= vHi && wid[k] / loc[k] <= GreenScore.MAX_REL_WIDTH
            for k in eachindex(loc)]
    h = h[keep]; loc = loc[keep]; wid = wid[keep]; prom = prom[keep]
    f = GreenScore.fitPeaks(v, Gb, h, loc, wid)
    hMax = max(isempty(h) ? -Inf : maximum(h), eps())
    return [Peak(f.vel[j], f.polVel[j], f.gamma[j], f.poleOk[j], h[j], prom[j], h[j] - prom[j],
                 f.poleHeight[j], wid[j] / f.vel[j], h[j] / hMax,
                 NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN) for j in eachindex(loc)]
end

"""
    normaliser(G, velAxis, vLo, vHi) -> Float64

log10 of the median |G13_eff| over the band (finite, positive samples); NaN if none.
"""
function normaliser(G::AbstractVector{<:Number}, velAxis::AbstractVector{<:Real},
                    vLo::Real, vHi::Real)
    x = Float64[]
    for i in eachindex(velAxis)
        s = abs(G[i]); vi = velAxis[i]
        (vi >= vLo && vi <= vHi && isfinite(s) && s > 0) && push!(x, s)
    end
    isempty(x) && return NaN
    return log10(median(x))
end

"""
    surfaceSH(v, g11, g22, g33) -> Float64

Shear-horizontal share `|r22| / (|r11| + |r22| + |r33|)` of a pole from 3 samples of the
diagonal G (residues of GreenScore.poleFit3's model); NaN if a fit is ill-conditioned.
"""
function surfaceSH(v, g11, g22, g33)
    r11 = GreenScore.poleFit3(v, g11)[4]
    r22 = GreenScore.poleFit3(v, g22)[4]
    r33 = GreenScore.poleFit3(v, g33)[4]
    return abs(r22) / max(abs(r11) + abs(r22) + abs(r33), floatmin())
end

"""
    polarisationStencil(pk) -> V

`[3 x nPk]` velocities (pole - FIT_SPAN, pole, pole + FIT_SPAN) of the full-tensor solve
(the reported vel where the pole fit failed).
"""
function polarisationStencil(pk::AbstractVector)
    n = length(pk)
    V = Matrix{Float64}(undef, 3, n)
    for (j, p) in pairs(pk)
        vp = (p.poleOk && isfinite(p.polVel)) ? p.polVel : p.vel
        V[1, j] = vp - GreenScore.FIT_SPAN; V[2, j] = vp; V[3, j] = vp + GreenScore.FIT_SPAN
    end
    return V
end

"""
    polarisationFromSolve(pk, V, g11, g22, g33) -> Vector{Peak}

SH of every candidate from the diagonal G11/G22/G33 (`[3 x nPk]` each) on the stencil V.
"""
polarisationFromSolve(pk::AbstractVector{Peak}, V, g11, g22, g33) =
    Peak[_with(pk[j]; SH = surfaceSH(view(V, :, j), view(g11, :, j), view(g22, :, j),
                                     view(g33, :, j))) for j in eachindex(pk)]

"""
    polarisation(structureRotated, pk) -> Vector{Peak}

One full-tensor Green solve (3 points per candidate, computeMinus = false) and SH of every
candidate.
"""
function polarisation(sr::Structure, pk::AbstractVector{Peak})
    isempty(pk) && return pk
    n = length(pk)
    V = polarisationStencil(pk)
    g = scoreGreenSolve(sr, vec(V); eta = SawSpectrum.ETA, computeMinus = false)
    GF = reshape(g.GFull, 3, 3, 3, n)
    return polarisationFromSolve(pk, V, GF[1, 1, :, :], GF[2, 2, :, :], GF[3, 3, :, :])
end

"""
    score(pk, N) -> Vector{Peak}

Gate (own periods, tau = GATE_OWN_PERIODS / vel), dephasing (GreenScore.dephaseFactor with
vRef = GATE_PERIODS / tau), polarisation and the score `logStrength - N + polTerm`.
"""
function score(pk::AbstractVector{Peak}, N::Real)
    sLB = GreenScore.PROBE_HALF * GreenScore.PROBE_STEP / GreenScore.PROBE_OFFSET
    out = Vector{Peak}(undef, length(pk))
    for (j, p) in pairs(pk)
        tau = GATE_OWN_PERIODS / p.vel
        sl = abs(p.slope); isfinite(sl) || (sl = sLB)
        g = p.gamma
        gateTerm = -2 * pi * g * tau / log(10)
        F = GreenScore.dephaseFactor(g, sl * GreenScore.APERTURE_SIGMA, GreenScore.GATE_PERIODS / tau)
        dT = log10(max(F, floatmin()))
        logS = log10(max(p.poleHeight - p.base, floatmin())) + gateTerm + dT
        pol = K_POL * (1 - 2 * p.SH)
        isfinite(pol) || (pol = -K_POL)
        out[j] = _with(p; gateTerm = gateTerm, dephaseTerm = dT, logStrength = logS,
                       polTerm = pol, norm = N, score = logS - N + pol)
    end
    return out
end

"""
    direction(structureBase, angle, G, velAxis, vSmin, vLmax) -> Vector{Peak}

All steps for one direction: candidates -> probes -> polarisation -> score
(threshold-independent). `structureBase` UNROTATED.
"""
function direction(structureBase::Structure, angle::Real, G::AbstractVector{<:Number},
                   velAxis::AbstractVector{<:Real}, vSmin::Real, vLmax::Real)
    vLo, vHi = band(vSmin, vLmax)
    pk = candidates(G, velAxis, vLo, vHi)
    isempty(pk) && return pk
    N = normaliser(G, velAxis, vLo, vHi)
    pk = GreenScore.probeSlope(structureBase, angle, pk)
    pk = polarisation(rotateStructure(structureBase, angle), pk)
    return score(pk, N)
end

end # module GreenScoreLayered
