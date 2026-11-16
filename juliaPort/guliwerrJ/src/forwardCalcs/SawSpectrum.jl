# =========================================================================================
# SawSpectrum.jl -- port of guliwerr/forwardCalcs/SawSpectrum.m and of
# GreenPeaks.m (the generic |G| peak finder, a sibling module below).
#
# EXPERIMENTAL (as in MATLAB): checked against the measured maps published in Grabec et al.
# 2024 and the UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
#
# SawSpectrum is the late-time spectrum machinery the existence scores share: the damping of
# the default existence solves (ETA), the axis lattice (latticeAxis, velocityAxis), the
# late-time gate (lateTime), the peak search (peaksOf) and the scored candidates of one
# spectrum (candidates). See the MATLAB header and
# docs/reference/architecture/forwardCalcs_modeClassification.md for the physics.
#
# A MATLAB static namespace -> a Julia submodule of ForwardCalcs, so the call sites read as
# in MATLAB: `SawSpectrum.lateTime(G, v, vRef)`, `SawSpectrum.ETA`.
#
# Efficiency (the late-time spectrum runs on every 'refined'/'ritz' angle, three times on a
# GreenSpotScore angle): FFTW plans are cached per (length) in a locked Dict and executed
# with `mul!` into preallocated buffers; the window, gate and output axis are built once per
# call (they depend on the axis and vRef only) and reused over the columns.
# =========================================================================================

module SawSpectrum

import FFTW
using FFTW: ESTIMATE
using LinearAlgebra: mul!
using ...CoreTypes: findpeaks_matlab

const ETA           = 3e-3
const GATE_PERIODS  = 3
const RAMP_PERIODS  = 0.6
const MAX_REL_WIDTH = 0.12
const VEL_MIN       = 0.3
const VEL_STEP      = 0.004
const PAD_FACTOR    = 4
const TAIL_TAPER    = 0.2
const EDGE_TAPER    = 0.1
const AXIS_HI_X_VL  = 1.05

"""
    velocityAxis(vMax, step = VEL_STEP) -> Vector{Float64}

`SawSpectrum.velocityAxis`: `(VEL_MIN:step:vMax + EDGE_TAPER)'`. Built with MATLAB's colon
semantics (`VEL_MIN + k*step`, k = 0..floor((hi-lo)/step + tol)), see `_colon`.
"""
velocityAxis(vMax::Real, step::Real = VEL_STEP) = _colon(VEL_MIN, step, vMax + EDGE_TAPER)

"""
MATLAB `lo:step:hi` for positive step: n = floor((hi-lo)/step + 4 eps-ish tolerance), the
points lo + k*step for k < n/2 and hi - (n-k)*step from the far end (MATLAB fills a colon
range from both ends). Reproduced exactly so the Ritz-only axis (velocityAxis) matches.
"""
function _colon(lo::Float64, step::Float64, hi::Float64)
    # MATLAB colon: n = floor((hi - lo)/step + tol), tol = 3*eps(max(|lo|,|hi|))/|step|-ish;
    # values a(k) = lo + k*step for k <= n/2, a(k) = hi_eff - (n-k)*step above, where
    # hi_eff = lo + n*step.
    tol = 3 * eps(max(abs(lo), abs(hi))) / abs(step)
    n = floor(Int, (hi - lo) / step + tol)
    out = Vector{Float64}(undef, n + 1)
    hiEff = lo + n * step
    half = n ÷ 2
    @inbounds for k in 0:n
        out[k+1] = k <= half ? lo + k * step : hiEff - (n - k) * step
    end
    return out
end
_colon(lo::Real, step::Real, hi::Real) = _colon(Float64(lo), Float64(step), Float64(hi))

"""
    latticeAxis(nLo, nHi, step) -> Vector{Float64}

`VEL_MIN + (nLo:nHi)' * step`: the one builder of every score axis, so a point has the same
value in a short axis and in a longer one with the same step.
"""
latticeAxis(nLo::Integer, nHi::Integer, step::Real) =
    [VEL_MIN + k * Float64(step) for k in nLo:nHi]

# ---- FFT plan cache ---------------------------------------------------------------------
const _PLAN_LOCK = ReentrantLock()
const _FWD = Dict{Int,Any}()
const _INV = Dict{Int,Any}()

# FFTW's own `plan_fft` method, selected with `invoke`: AppleAccelerate.jl's AbstractFFTs
# extension (if a user loads that package) adds a more specific `plan_fft(::Vector{ComplexF64})`
# that accepts power-of-2 lengths only and would throw here (see BlasBackend.jl). The inverse
# is FFTW's `inv` of that plan (backward plan scaled by 1/n), as `plan_ifft` builds it.
_fftwPlanFwd(n::Int) = invoke(FFTW.plan_fft, Tuple{StridedArray{ComplexF64,1},Any},
                              Vector{ComplexF64}(undef, n), (1,); flags = ESTIMATE)
function _fwdPlan(n::Int)
    lock(_PLAN_LOCK) do
        get!(() -> _fftwPlanFwd(n), _FWD, n)
    end
end
function _invPlan(n::Int)
    lock(_PLAN_LOCK) do
        get!(() -> inv(_fftwPlanFwd(n)), _INV, n)
    end
end

"""
    fillNaNLinear!(g::AbstractVector{ComplexF64}) -> g

MATLAB `fillmissing(G, 'linear', 1, 'EndValues', 'nearest')` on one column: NaN entries
(either part NaN) between finite ones are linearly interpolated over the sample index; leading
/ trailing NaNs take the nearest finite value. All-NaN columns are left as they are.
"""
function fillNaNLinear!(g::AbstractVector{ComplexF64})
    n = length(g)
    good = findall(z -> !isnan(z), g)
    isempty(good) && return g
    i1 = first(good); i2 = last(good)
    @inbounds for i in 1:i1-1; g[i] = g[i1]; end
    @inbounds for i in i2+1:n; g[i] = g[i2]; end
    @inbounds for k in 1:length(good)-1
        a = good[k]; b = good[k+1]
        b - a > 1 || continue
        ga = g[a]; gb = g[b]
        for i in a+1:b-1
            t = (i - a) / (b - a)
            g[i] = ga + (gb - ga) * t
        end
    end
    return g
end

"""
    lateTime(G, velAxis, vRef; pad = PAD_FACTOR, tailFrac = TAIL_TAPER,
             gatePeriods = GATE_PERIODS) -> (A, vOut)
    lateTime(G, velAxis, vRef, pad, tailFrac, gatePeriods)   # MATLAB positional order;
                                                             # `nothing` = default

`SawSpectrum.lateTime`: |late-time spectrum| with the first `gatePeriods` periods of `vRef`
removed. `G` is a complex vector (one direction) or matrix `[nVel x nAng]` on the uniform
`velAxis` [km/s]. Returns `A` (`[pad*nVel x nAng]` or a vector, NaN on the tapered axis ends)
and its axis `vOut`. RitzScore passes `pad = 1, tailFrac = 0`.
"""
lateTime(G, velAxis, vRef, pad, tailFrac = nothing, gatePeriods = nothing) =
    lateTime(G, velAxis, vRef; pad = something(pad, PAD_FACTOR),
             tailFrac = something(tailFrac, TAIL_TAPER),
             gatePeriods = something(gatePeriods, GATE_PERIODS))

function lateTime(G::AbstractVector{<:Number}, velAxis::AbstractVector{<:Real}, vRef::Real;
                  pad::Integer = PAD_FACTOR, tailFrac::Real = TAIL_TAPER,
                  gatePeriods::Real = GATE_PERIODS)
    A, vOut = lateTime(reshape(G, :, 1), velAxis, vRef; pad, tailFrac, gatePeriods)
    return vec(A), vOut
end

function lateTime(G::AbstractMatrix{<:Number}, velAxis::AbstractVector{<:Real}, vRef::Real;
                  pad::Integer = PAD_FACTOR, tailFrac::Real = TAIL_TAPER,
                  gatePeriods::Real = GATE_PERIODS)
    nVel = length(velAxis)
    size(G, 1) == nVel || throw(DimensionMismatch("SawSpectrum.lateTime: G has $(size(G,1)) rows, velAxis $nVel"))
    v = velAxis
    dv = Float64(v[2] - v[1])
    tAx = [(k) / (nVel * dv) for k in 0:nVel-1]
    edgeW = max(2, Int(mround_(EDGE_TAPER / dv)))
    w = ones(nVel)
    @inbounds for k in 0:edgeW-1
        r = 0.5 * (1 - cos(pi * k / edgeW))
        w[k+1] = r
        w[nVel-k] = r
    end
    # gate
    tCut = gatePeriods / vRef
    ramp = RAMP_PERIODS / vRef
    tHalf = tAx[end] / 2
    gate = zeros(nVel)
    @inbounds for i in 1:nVel
        t = tAx[i]
        (t >= tCut && t < tHalf) || continue
        gate[i] = t < tCut + ramp ? 0.5 * (1 - cos(pi * (t - tCut) / ramp)) : 1.0
    end
    if tailFrac > 0
        tEnd = tHalf; t0 = tEnd * (1 - tailFrac)
        @inbounds for i in 1:nVel
            if gate[i] > 0 && tAx[i] > t0
                gate[i] = gate[i] * 0.5 * (1 + cos(pi * (tAx[i] - t0) / (tEnd - t0)))
            end
        end
    end
    nOut = pad * nVel
    pf = _fwdPlan(nVel)
    pin = pad > 1 ? _invPlan(nOut) : _invPlan(nVel)
    buf = Vector{ComplexF64}(undef, nVel)
    S = Vector{ComplexF64}(undef, nVel)
    Sg = zeros(ComplexF64, nOut)
    out = Vector{ComplexF64}(undef, nOut)
    nAng = size(G, 2)
    A = Matrix{Float64}(undef, nOut, nAng)
    for ja in 1:nAng
        @inbounds for i in 1:nVel
            buf[i] = G[i, ja]
        end
        any(z -> !isfinite(z), buf) && fillNaNLinear!(buf)
        @inbounds for i in 1:nVel
            buf[i] *= w[i]
        end
        mul!(S, pf, buf)
        @inbounds for i in 1:nVel
            Sg[i] = S[i] * gate[i]
        end
        mul!(out, pin, Sg)
        @inbounds for i in 1:nOut
            A[i, ja] = abs(out[i]) * pad
        end
    end
    vOut = [Float64(v[1]) + k * dv / pad for k in 0:nOut-1]
    # wOut = interp1(velAxis, w, vOut, 'linear', 0); A(wOut < 1, :) = NaN
    j = 1
    @inbounds for i in 1:nOut
        x = vOut[i]
        wo = 0.0
        if x >= v[1] && x <= v[end]
            while j < nVel - 1 && x > v[j+1]
                j += 1
            end
            if x == v[j+1]
                wo = w[j+1]
            else
                wo = w[j] + (w[j+1] - w[j]) * ((x - v[j]) / (v[j+1] - v[j]))
            end
        end
        if wo < 1
            for ja in 1:nAng
                A[i, ja] = NaN
            end
        end
    end
    return A, vOut
end

mround_(x) = round(x, RoundNearestTiesAway)

"""
    PeakSet

`SawSpectrum.peaksOf`'s struct of column vectors: `vel`, `prom`, `height`, `relWidth`.
"""
struct PeakSet
    vel::Vector{Float64}
    prom::Vector{Float64}
    height::Vector{Float64}
    relWidth::Vector{Float64}
end
PeakSet() = PeakSet(Float64[], Float64[], Float64[], Float64[])
Base.isempty(P::PeakSet) = isempty(P.vel)
Base.length(P::PeakSet) = length(P.vel)

"""
    findpeaksProm(y, x, minProm) -> (height, loc, width, prom)

MATLAB `findpeaks(y, x, 'WidthReference', 'halfprom', 'MinPeakProminence', minProm)`:
[`findpeaks_matlab`](@ref) (the verified port: boundary exclusion, edge minima) with the
prominence floor applied after the prominence computation, as MATLAB orders it.
"""
function findpeaksProm(y::AbstractVector{<:Real}, x::AbstractVector{<:Real}, minProm::Real)
    length(y) < 3 && return (Float64[], Float64[], Float64[], Float64[])
    pk = findpeaks_matlab(y, x)
    keep = pk.prominence .>= minProm
    return (Float64.(pk.height[keep]), Float64.(pk.loc[keep]),
            Float64.(pk.width[keep]), Float64.(pk.prominence[keep]))
end

"""
    peaksOf(A, velAxis, vSearchMax) -> PeakSet

`SawSpectrum.peaksOf`: the local maxima of `A` below `vSearchMax` (search runs to
`vSearchMax + 0.05`, docs/pitfalls.md MTH4) with half-prominence width <= MAX_REL_WIDTH x
velocity; prominence floor 1e-4 x the largest value.
"""
function peaksOf(A::AbstractVector{<:Real}, velAxis::AbstractVector{<:Real}, vSearchMax::Real)
    m = findall(i -> isfinite(A[i]) && velAxis[i] < vSearchMax + 0.05, eachindex(A))
    length(m) < 3 && return PeakSet()
    y = A[m]; x = velAxis[m]
    h, loc, wid, prom = findpeaksProm(y, x, 1e-4 * maximum(y))
    keep = findall(k -> loc[k] < vSearchMax && wid[k] / loc[k] <= MAX_REL_WIDTH, eachindex(loc))
    return PeakSet(loc[keep], prom[keep], h[keep], wid[keep] ./ loc[keep])
end

"""
    Candidate

One scored candidate of `candidates`: `vel`, `score = log10(prominence / ref.floor)`,
`prominence`, `height`, `relWidth`, `relHeight` (height / tallest candidate here -- display
only).
"""
struct Candidate
    vel::Float64
    score::Float64
    prominence::Float64
    height::Float64
    relWidth::Float64
    relHeight::Float64
end

"""
    candidates(G, velAxis, vSearchMax, vRef, floor; pad = PAD_FACTOR, tailFrac = TAIL_TAPER)
        -> Vector{Candidate}

`SawSpectrum.candidates(G, velAxis, vSearchMax, ref, pad, tailFrac)` with `ref.vRef`,
`ref.floor` passed as numbers: the peaks of the late-time spectrum (peaksOf), scored as
`log10(prominence / floor)`. RitzScore passes `pad = 1, tailFrac = 0`.
"""
function candidates(G::AbstractVector{<:Number}, velAxis::AbstractVector{<:Real},
                    vSearchMax::Real, vRef::Real, floor::Real;
                    pad::Integer = PAD_FACTOR, tailFrac::Real = TAIL_TAPER)
    A, vO = lateTime(G, velAxis, vRef; pad, tailFrac)
    P = peaksOf(A, vO, vSearchMax)
    hMax = max(isempty(P.height) ? -Inf : maximum(P.height), eps())
    return [Candidate(P.vel[j], log10(P.prom[j] / floor), P.prom[j], P.height[j],
                      P.relWidth[j], P.height[j] / hMax) for j in eachindex(P.vel)]
end

end # module SawSpectrum

# =========================================================================================
# GreenPeaks -- port of guliwerr/forwardCalcs/GreenPeaks.m: the generic |G| peak finder
# (isoFreq spectra; not part of the surface-wave classification).
# =========================================================================================

module GreenPeaks

using ...CoreTypes: findpeaks_matlab

const HEIGHT_FRAC = 0.15

"""
    Peak

One `GreenPeaks.find` peak: `vel`, `height`, `prominence`, `width` (half-prominence, in the
units of the axis).
"""
struct Peak
    vel::Float64
    height::Float64
    prominence::Float64
    width::Float64
end

"""
    find(absG, velAxis, mask, heightFrac = HEIGHT_FRAC) -> Vector{Peak}

`GreenPeaks.find`: local maxima of `absG[mask]` (on `velAxis[mask]`, increasing) whose
height is at least `heightFrac x maximum(absG)` over the WHOLE array. Prominence and width
are reported, not filtered on.
"""
function find(absG::AbstractVector{<:Real}, velAxis::AbstractVector{<:Real},
              mask::AbstractVector{Bool}, heightFrac::Real = HEIGHT_FRAC)
    any(mask) || return Peak[]
    globalMax = maximum(absG)
    y = absG[mask]; x = velAxis[mask]
    length(y) < 3 && return Peak[]
    pk = findpeaks_matlab(y, x; minheight = heightFrac * globalMax)
    return [Peak(pk.loc[i], pk.height[i], pk.prominence[i], pk.width[i]) for i in eachindex(pk.loc)]
end

end # module GreenPeaks
