# =========================================================================================
# RitzScore.jl -- port of guliwerr/forwardCalcs/RitzScore.m: the Ritz-only
# existence score. A gated Ritz modal sum S(v) = sum_n u3_n conj(u1_n) / (v_n^2 - v^2 (1 +
# i ETA)) that mimics G13, scored with the shared late-time machinery (SawSpectrum) against
# the direction's own floor; half-space rule (direction / acceptedHS: dephased score,
# freqDisp veto, window mass, bottom-trapped gate) and layered rule (directionLayered /
# acceptedLayered: polarisation term, trapped gate). See the MATLAB header and
# docs/reference/architecture/forwardCalcs_modeClassification.md (the Ritz method).
#
# EXPERIMENTAL (as in MATLAB): checked against the measured maps published in Grabec et al.
# 2024 and the UTGS paper (arXiv:2510.10696); the scores, their scales and the interface may change.
#
# Structure of the port: `directionCore` / `directionLayeredCore` hold the whole rule on plain
# arrays (mode velocities, group velocities, modal numerators and profile measures, and a
# slope function for the supporting modes), so they are golden-tested on MATLAB's own Ritz
# data; `direction` / `directionLayered` / `modeProfiles` read them off a RitzSolverSAW
# (calcDisplacement, angleSlope).
# =========================================================================================

module RitzScore

using LinearAlgebra: Symmetric, eigvals
using Statistics: median
using ...CoreTypes: Structure, Material, nlay
using ...Solvers: Solvers, calcDisplacement
using ..ForwardCalcs: thresholdValue, surfaceVRef, standardElastTensor, _with
using ..SawSpectrum: SawSpectrum, latticeAxis
using ..GreenScore: GreenScore

const VEL_STEP        = 0.002
const ETA             = 1e-2
const MODE_VEL_FACTOR = 1.2
const HS_THRESHOLD_HIGH   = 0.16
const HS_THRESHOLD_LOW    = 0.70
const HS_MIN_WINDOW_MASS  = -1.55
const HS_WINDOW_DEPTH_REF = 30
const HS_BOTTOM10_MAX     = 0.15
const HS_BOTTOM_RATIO_MAX = 0.8
const HS_MIN_LEG_PER_DEPTH = 2.3
const HS_FLOOR_LO   = 0.8
const WINDOW_REL    = 0.012
const WINDOW_ABS    = 0.030
const WINDOW_TOP    = 0.98
const LOCAL_MIN_HALF = 0.004
const LAY_THRESHOLD_LOW  = 0.65
const LAY_THRESHOLD_HIGH = -0.7
const LAY_K_POL          = 0.75
const LAY_EBOT_MAX       = 1e-2
const MIN_FREQDISP       = 0.3

"""
    Peak

One Ritz candidate (MATLAB `RitzScore.emptyPeaks()` / ModeClassification.ritzScorePeaks
field set; both rules fill the same fields): `vel` (modal-sum peak), `modeNo` / `modeVel`
(local supporting mode), `score`, `dScore` (= score + aperture dephasing), `prominence`,
`relWidth`, `relHeight`, `share`, `freqDisp`, `slope` [km/s per deg], `windowMass`,
`belowShear`, `vetoed`, `bottom10`, `bottom25`, `SH`, `polTerm`, `accepted` (false here;
the classification sets it).
"""
struct Peak
    vel::Float64
    modeNo::Int
    modeVel::Float64
    score::Float64
    dScore::Float64
    prominence::Float64
    relWidth::Float64
    relHeight::Float64
    share::Float64
    freqDisp::Float64
    slope::Float64
    windowMass::Float64
    belowShear::Bool
    vetoed::Bool
    bottom10::Float64
    bottom25::Float64
    SH::Float64
    polTerm::Float64
    accepted::Bool
end

"""
    OwnReference

`ownReference` result: `vRef` [km/s] (surface layer sqrt(G_Hill/rho)), `floor` (median
late-time |S| over the floor band), `slabRes` (half-space: slabResolution [km/s], the
Rayleigh side's tolerance; NaN until set -- direction() requires it).
"""
struct OwnReference
    vRef::Float64
    floor::Float64
    slabRes::Float64
end
OwnReference(vRef, floor) = OwnReference(vRef, floor, NaN)

"""
    ModeProfiles

`modeProfiles` result for the kept modes 1..nKeep: `num` (u3 conj(u1) at the free surface),
`E1` (|u|^2 fraction within one wavelength of the surface), `r` (|u3(surface)|^2 relative to
the mode's spread), `Ebot` ([nKeep x 2], fraction in the deepest 10 % / 25 %), `SH`
(|u2|^2/|u|^2 at the surface).
"""
struct ModeProfiles
    num::Vector{ComplexF64}
    E1::Vector{Float64}
    r::Vector{Float64}
    Ebot::Matrix{Float64}
    SH::Vector{Float64}
end

# ---- thresholds / acceptance ------------------------------------------------------------

resolveThresholdHS(t) = first(thresholdValue(t, HS_THRESHOLD_HIGH, HS_THRESHOLD_LOW,
                                             "RitzScore:resolveThresholdHS:badName"))
resolveThresholdLayered(t) = first(thresholdValue(t, LAY_THRESHOLD_HIGH, LAY_THRESHOLD_LOW,
                                                  "RitzScore:resolveThresholdLayered:badName"))

"""
    candidateScore(pk, isLayered) -> Vector{Float64}

The score the rule thresholds: half-space `dScore`; layered `dScore + polTerm`.
"""
candidateScore(pk::AbstractVector{Peak}, isLayered::Bool) =
    isLayered ? [p.dScore + p.polTerm for p in pk] : [p.dScore for p in pk]

"""
    acceptedHS(pk, threshold) -> Vector{Bool}

`~vetoed & dScore >= thr & (belowShear | windowMass >= HS_MIN_WINDOW_MASS & ~bottomTrapped)`,
bottomTrapped = `bottom10 > HS_BOTTOM10_MAX & bottom10/bottom25 > HS_BOTTOM_RATIO_MAX`.
"""
function acceptedHS(pk::AbstractVector{Peak}, threshold)
    thr = resolveThresholdHS(threshold)
    return Bool[_acceptHS(p, thr) for p in pk]
end

function _acceptHS(p::Peak, thr::Real)
    trapped = p.bottom10 > HS_BOTTOM10_MAX && p.bottom10 / p.bottom25 > HS_BOTTOM_RATIO_MAX
    return !p.vetoed && p.dScore >= thr &&
           (p.belowShear || (p.windowMass >= HS_MIN_WINDOW_MASS && !trapped))
end

"""
    acceptedLayered(pk, threshold) -> Vector{Bool}

`bottom25 <= LAY_EBOT_MAX & dScore + polTerm >= threshold`.
"""
function acceptedLayered(pk::AbstractVector{Peak}, threshold)
    thr = resolveThresholdLayered(threshold)
    return Bool[p.bottom25 <= LAY_EBOT_MAX && p.dScore + p.polTerm >= thr for p in pk]
end

# ---- the modal sum and its pieces ----------------------------------------------------------

"""
    spectrum(vn, num, velAxis) -> Vector{ComplexF64}

Complex Ritz modal sum on a velocity axis: `S(v) = sum_n num_n / (vn_n^2 - v^2 (1 + i ETA))`.
"""
function spectrum(vn::AbstractVector{<:Real}, num::AbstractVector{<:Number},
                  velAxis::AbstractVector{<:Real})
    # Per point: d_n = (vn_n^2 - v^2) - i v^2 ETA, term = num_n conj(d_n) / |d_n|^2, in real
    # arithmetic so the mode loop vectorises (@simd: the sum is split over vector lanes, so
    # its order -- and MATLAB's complex division -- differ from MATLAB by rounding only;
    # measured ~1e-15 relative on the golden spectra, far inside the score tolerances).
    nV = length(velAxis); nK = length(num)
    vn2 = [Float64(vn[n])^2 for n in 1:nK]
    nr = [Float64(real(num[n])) for n in 1:nK]; ni = [Float64(imag(num[n])) for n in 1:nK]
    S = Vector{ComplexF64}(undef, nV)
    @inbounds for i in 1:nV
        x2 = Float64(velAxis[i])^2
        wi = x2 * ETA
        sr = 0.0; si = 0.0
        @simd for n in 1:nK
            dr = vn2[n] - x2
            inv = 1 / (dr * dr + wi * wi)
            # num * conj(d) / |d|^2, conj(d) = dr + i wi
            sr += (nr[n] * dr - ni[n] * wi) * inv
            si += (ni[n] * dr + nr[n] * wi) * inv
        end
        S[i] = complex(sr, si)
    end
    return S
end

"""
    decompose(vn, num, v) -> Vector{Float64}

Each mode's share |term_n| / sum |term| of the modal sum at velocity `v`.
"""
function decompose(vn, num, v::Real)
    t = [abs(num[n] / (vn[n]^2 - v^2 * complex(1.0, ETA))) for n in eachindex(num)]
    return t ./ max(sum(t), eps())
end

"""
    localMode(vn, num, vPeak, relWidth) -> (modeNo, share)

The supporting mode of a modal-sum peak: argmax |num_n / (v_n^2 - vPeak^2 (1 + i ETA))| over
the modes within max(relWidth vPeak / 2, LOCAL_MIN_HALF) of vPeak (the nearest mode if none);
`share` = that term's fraction of the sum over ALL modes.
"""
function localMode(vn::AbstractVector{<:Real}, num::AbstractVector{<:Number},
                   vPeak::Real, relWidth::Real)
    c = vPeak^2 * complex(1.0, ETA)
    t = [abs(num[n] / (vn[n]^2 - c)) for n in eachindex(num)]
    half = max(relWidth * vPeak / 2, LOCAL_MIN_HALF)
    modeNo = 0; best = -Inf
    for n in eachindex(vn)
        if abs(vn[n] - vPeak) <= half && t[n] > best
            best = t[n]; modeNo = n
        end
    end
    if modeNo == 0
        modeNo = argmin(abs.(vn .- vPeak))
    end
    return modeNo, t[modeNo] / max(sum(t), eps())
end

"""
    windowWeights(vn, fd, E1, r, vS, vL) -> Vector{Float64}

Per-mode near-surface weight: `r fd` below vS, `r E1 fd` at or above, 0 above
WINDOW_TOP x vL.
"""
function windowWeights(vn, fd, E1, r, vS::Real, vL::Real)
    w = r .* fd
    for n in eachindex(vn)
        vn[n] >= vS && (w[n] *= E1[n])
        vn[n] > WINDOW_TOP * vL && (w[n] = 0.0)
    end
    return w
end

"""
    windowMass(vn, w, vPeak, vS, above) -> Float64

log10 of the summed `w` over the modes within max(WINDOW_REL vPeak, WINDOW_ABS) of vPeak on
one side of vS (`above`: at or above), floored at 1e-12.
"""
function windowMass(vn, w, vPeak::Real, vS::Real, above::Bool)
    half = max(WINDOW_REL * vPeak, WINDOW_ABS)
    s = 0.0
    for n in eachindex(vn)
        (abs(vn[n] - vPeak) <= half && (vn[n] >= vS) == above) && (s += w[n])
    end
    return log10(max(s, 1e-12))
end

"""
    dephaseTerm(vMode, slope, vRef) -> Vector{Float64}

log10 GreenScore.dephaseFactor(ETA vMode / 2, |slope| APERTURE_SIGMA, vRef), floored at -6.
"""
dephaseTerm(vMode::AbstractVector, slope::AbstractVector, vRef::Real) =
    [log10(max(GreenScore.dephaseFactor(ETA * vMode[i] / 2,
                                        abs(slope[i]) * GreenScore.APERTURE_SIGMA, vRef), 1e-6))
     for i in eachindex(vMode)]

# ---- the half-space rule ------------------------------------------------------------------

"""
    directionCore(vn, vg, num, E1, r, Ebot, velAxis, ref, vS, vL, relDepth, slopeOf)
        -> Vector{Peak}

The half-space rule on plain arrays (RitzScore.direction without the solver): `vn`, `vg` the
phase / group velocities of the kept modes 1..nKeep, `num, E1, r, Ebot` modeProfiles
outputs, `ref::OwnReference` (with slabRes), `relDepth` the slab depth, `slopeOf(um)` the
angular slopes [km/s per deg] of the distinct supporting modes `um` (sorted).
"""
function directionCore(vn::AbstractVector{<:Real}, vg::AbstractVector{<:Real},
                       num::AbstractVector{<:Number}, E1, r, Ebot::AbstractMatrix,
                       velAxis::AbstractVector{<:Real}, ref::OwnReference, vS::Real, vL::Real,
                       relDepth::Real, slopeOf)
    nK = length(num)
    (isfinite(ref.slabRes) && size(Ebot) == (nK, 2)) || throw(ArgumentError(
        "RitzScore:direction:missingInput -- RitzScore.direction needs ref.slabRes " *
        "(slabResolution) and Ebot [nKeep x 2] (modeProfiles)."))
    vn = Float64.(vn[1:nK]); vg = Float64.(vg[1:nK])
    fd = [exp(-10 * ((vn[i] - vg[i]) / max(vn[i], eps()))^2) for i in 1:nK]
    c = SawSpectrum.candidates(spectrum(vn, num, velAxis), velAxis, vL, ref.vRef, ref.floor;
                               pad = 1, tailFrac = 0)
    nP = length(c)
    nP == 0 && return Peak[]
    w = windowWeights(vn, fd, E1, r, vS, vL)
    vB = vS + ref.slabRes
    wmShift = log10(relDepth / HS_WINDOW_DEPTH_REF)
    md = zeros(Int, nP); sh = zeros(nP); wm = zeros(nP); below = falses(nP)
    for j in 1:nP
        md[j], sh[j] = localMode(vn, num, c[j].vel, c[j].relWidth)
        below[j] = vn[md[j]] < vB
        wm[j] = windowMass(vn, w, c[j].vel, vS, !below[j]) + wmShift
    end
    um = sort(unique(md))
    slU = slopeOf(um)
    sl = [slU[searchsortedfirst(um, m)] for m in md]
    dT = dephaseTerm(vn[md], sl, ref.vRef)
    return [Peak(c[j].vel, md[j], vn[md[j]], c[j].score, c[j].score + dT[j], c[j].prominence,
                 c[j].relWidth, c[j].relHeight, sh[j], fd[md[j]], sl[j], wm[j], below[j],
                 fd[md[j]] < MIN_FREQDISP, Ebot[md[j], 1], Ebot[md[j], 2], NaN, NaN, false)
            for j in 1:nP]
end

"""
    angleSlope(ritz, modeNos, dKdC_cell) -> Vector{Float64}

The supporting modes' angular slopes dv/dtheta from the Ritz solver (MATLAB
`ritz.angleSlope`, Hellmann-Feynman), via `Solvers.angleSlope`.
"""
angleSlope(ritz, modeNos, dKdC_cell) = Solvers.angleSlope(ritz, modeNos, dKdC_cell)

"""
    direction(ritz, prof::ModeProfiles, velAxis, ref, vS, vL, dKdC_cell = nothing)
        -> Vector{Peak}

`RitzScore.direction`: scored half-space candidates of one direction (threshold-independent)
from a calc()'d RitzSolverSAW with eigenvectors present.
"""
function direction(ritz, prof::ModeProfiles, velAxis::AbstractVector{<:Real},
                   ref::OwnReference, vS::Real, vL::Real, dKdC_cell = nothing)
    nK = length(prof.num)
    relDepth = ritz.structure.calcParams[1].relDepth
    return directionCore(ritz.vels[1:nK], ritz.velsGroup[1:nK], prof.num, prof.E1, prof.r,
                         prof.Ebot, velAxis, ref, vS, vL, relDepth,
                         um -> angleSlope(ritz, um, dKdC_cell))
end

# ---- the layered rule ---------------------------------------------------------------------

"""
    directionLayeredCore(vn, num, Ebot, SH, velAxis, ref, vLo, vHi, vSsub, slopeOf)
        -> Vector{Peak}

The layered rule on plain arrays (RitzScore.directionLayered without the solver).
"""
function directionLayeredCore(vn::AbstractVector{<:Real}, num::AbstractVector{<:Number},
                              Ebot::AbstractMatrix, SH::AbstractVector,
                              velAxis::AbstractVector{<:Real}, ref::OwnReference,
                              vLo::Real, vHi::Real, vSsub::Real, slopeOf)
    nK = length(num)
    vn = Float64.(vn[1:nK])
    c = SawSpectrum.candidates(spectrum(vn, num, velAxis), velAxis, vHi, ref.vRef, ref.floor;
                               pad = 1, tailFrac = 0)
    c = filter(x -> x.vel >= vLo, c)
    nP = length(c)
    nP == 0 && return Peak[]
    md = zeros(Int, nP); sh = zeros(nP)
    for j in 1:nP
        md[j], sh[j] = localMode(vn, num, c[j].vel, c[j].relWidth)
    end
    um = sort(unique(md))
    slU = slopeOf(um)
    sl = [slU[searchsortedfirst(um, m)] for m in md]
    dT = dephaseTerm(vn[md], sl, ref.vRef)
    out = Vector{Peak}(undef, nP)
    for j in 1:nP
        m = md[j]
        out[j] = Peak(c[j].vel, m, vn[m], c[j].score, c[j].score + dT[j], c[j].prominence,
                      c[j].relWidth, c[j].relHeight, sh[j], NaN, sl[j], NaN, vn[m] < vSsub,
                      false, Ebot[m, 1], Ebot[m, 2], SH[m], LAY_K_POL * (1 - 2 * SH[m]), false)
    end
    return out
end

"""
    directionLayered(ritz, prof::ModeProfiles, velAxis, ref, vLo, vHi, vSsub, dKdC_cell = nothing)
        -> Vector{Peak}
"""
function directionLayered(ritz, prof::ModeProfiles, velAxis::AbstractVector{<:Real},
                          ref::OwnReference, vLo::Real, vHi::Real, vSsub::Real,
                          dKdC_cell = nothing)
    nK = length(prof.num)
    return directionLayeredCore(ritz.vels[1:nK], prof.num, prof.Ebot, prof.SH, velAxis, ref,
                                vLo, vHi, vSsub, um -> angleSlope(ritz, um, dKdC_cell))
end

"""
    layeredAxis(vLo, vHi) -> Vector{Float64}

VEL_STEP lattice from min(VEL_MIN, the lattice point at or below vLo) to vHi + 0.1 +
EDGE_TAPER.
"""
function layeredAxis(vLo::Real, vHi::Real)
    v0 = SawSpectrum.VEL_MIN; dv = VEL_STEP
    n1 = min(0, floor(Int, (vLo - v0) / dv + 1e-9))
    n2 = floor(Int, (vHi + 0.1 + SawSpectrum.EDGE_TAPER - v0) / dv + 1e-9)
    return latticeAxis(n1, n2, dv)
end

# ---- floors, profiles, slab resolution ---------------------------------------------------------

"""
    ownReference(structure, vn, num, velAxis, fLo, fHi) -> OwnReference

Floor of ONE direction from its modal sum: `vRef = surfaceVRef(structure)`, `floor` = median
late-time |S| (pad 1, no taper) over [fLo, fHi]. `slabRes` NaN (set it for the half-space).
"""
function ownReference(structure::Structure, vn, num, velAxis::AbstractVector{<:Real},
                      fLo::Real, fHi::Real)
    vRef = surfaceVRef(structure)
    A, _ = SawSpectrum.lateTime(spectrum(vn, num, velAxis), velAxis, vRef; pad = 1, tailFrac = 0)
    x = [A[i] for i in eachindex(A) if isfinite(A[i]) && velAxis[i] >= fLo && velAxis[i] <= fHi]
    return OwnReference(vRef, isempty(x) ? NaN : median(x))
end

withSlabRes(ref::OwnReference, slabRes::Real) = OwnReference(ref.vRef, ref.floor, slabRes)

"MATLAB linspace(d1, d2, n): d1 + (k (d2 - d1)) / (n - 1), last point d2."
function mlinspace(d1::Real, d2::Real, n::Integer)
    n == 1 && return [Float64(d2)]
    y = [d1 + (k * (d2 - d1)) / (n - 1) for k in 0:n-1]
    y[end] = d2
    return y
end

"""
    modeProfiles(ritz, vLong) -> ModeProfiles

Modal-sum weights and near-surface measures of the kept modes (vel <= MODE_VEL_FACTOR x
vLong), from one calcDisplacement call (surface = last row). Half-space: plain |u|^2 on the
uniform grid; layered: density- and trapezoid-weighted over the whole domain.
"""
function modeProfiles(ritz, vLong::Real)
    isempty(ritz.eigVecs_con) && throw(ArgumentError("RitzScore:modalNumerators:noEigVecs -- " *
        "Ritz eigenvectors are needed (compute before retentionLevel pruning)."))
    nKeep = findlast(<=(MODE_VEL_FACTOR * vLong), ritz.vels[1:ritz.nCompModes])
    nKeep === nothing && return ModeProfiles(ComplexF64[], Float64[], Float64[], zeros(0, 2), Float64[])
    D = calcDisplacement(ritz, nKeep)                 # [nPts x nKeep x 3], surface = last row
    return profilesFromDisplacement(ritz.structure, D)
end

"""
    profilesFromDisplacement(structure, D) -> ModeProfiles

modeProfiles' measures from a calcDisplacement array `D` ([nRows x nKeep x 3]).
"""
function profilesFromDisplacement(s::Structure, D::AbstractArray{<:Number,3})
    nRows, nKeep, _ = size(D)
    us = [D[nRows, n, d] for n in 1:nKeep, d in 1:3]
    num = us[:, 3] .* conj.(us[:, 1])
    u2 = abs2.(us)
    SH = [u2[n, 2] / max(u2[n, 1] + u2[n, 2] + u2[n, 3], floatmin()) for n in 1:nKeep]
    P = [abs2(D[i, n, 1]) + abs2(D[i, n, 2]) + abs2(D[i, n, 3]) for i in 1:nRows, n in 1:nKeep]
    if nlay(s) > 1
        E1, r, Ebot = layeredProfiles(s, P, us)
        return ModeProfiles(num, E1, r, Ebot, SH)
    end
    relDepth = s.calcParams[1].relDepth
    z = reverse(mlinspace(0.0, relDepth, nRows))
    tot = vec(sum(P; dims = 1))
    E1 = [sum(P[i, n] for i in 1:nRows if z[i] <= 1; init = 0.0) for n in 1:nKeep] ./ max.(tot, eps())
    r = abs2.(us[:, 3]) ./ max.(vec(sum(P; dims = 1)) ./ nRows, eps()) ./ relDepth
    Et = max.(tot, eps())
    Ebot = hcat([sum(P[i, n] for i in 1:nRows if z[i] >= 0.90 * relDepth; init = 0.0) for n in 1:nKeep] ./ Et,
                [sum(P[i, n] for i in 1:nRows if z[i] >= 0.75 * relDepth; init = 0.0) for n in 1:nKeep] ./ Et)
    return ModeProfiles(num, E1, r, Ebot, SH)
end

function layeredProfiles(s::Structure, Pw::AbstractMatrix, us::AbstractMatrix)
    nL = nlay(s)
    nP = size(Pw, 1) ÷ nL
    z = zeros(nP, nL); rho = zeros(nP, nL); wz = zeros(nP, nL)
    zTop = 0.0
    for il in nL:-1:1
        h = s.calcParams[il].relDepth
        z[:, il] = reverse(mlinspace(zTop, zTop + h, nP))
        rho[:, il] .= s.mater[il].rho
        dz = h / (nP - 1); w = fill(dz, nP); w[1] = dz / 2; w[end] = dz / 2
        wz[:, il] = w
        zTop += h
    end
    z = vec(z); wz = vec(wz) .* vec(rho); dTot = zTop
    Etot = max.(vec(wz' * Pw), floatmin())
    E1 = vec((wz .* (z .<= 1))' * Pw) ./ Etot
    r = abs2.(us[:, 3]) ./ Etot
    Ebot = hcat(vec((wz .* (z .>= 0.90 * dTot))' * Pw), vec((wz .* (z .>= 0.75 * dTot))' * Pw)) ./ Etot
    return E1, r, Ebot
end

"""
    modalNumerators(ritz, vLong) -> Vector{ComplexF64}
"""
modalNumerators(ritz, vLong::Real) = modeProfiles(ritz, vLong).num

"""
    traceVel(C, rho, q) -> Float64

Trace velocity along x1 of the slowest bulk wave with k_z/k_x = q (slabResolution):
`sqrt(min eig sym(Gamma(n)) / rho) sqrt(1 + q^2)`, n = [1, 0, q]/sqrt(1 + q^2).
"""
function traceVel(C::AbstractArray{<:Real,4}, rho::Real, q::Real)
    s = sqrt(1 + q^2)
    n = (1 / s, 0.0, q / s)
    G = zeros(3, 3)
    for i in 1:3, k in 1:3
        acc = 0.0
        for j in 1:3, l in 1:3
            acc += C[i, j, k, l] * n[j] * n[l]
        end
        G[i, k] = acc
    end
    return sqrt(minimum(eigvals(Symmetric((G + G') / 2))) / rho) * s
end

"""
    slabResolution(mat, vS, thetaStar, relDepth) -> Float64

How far above the true min shear velocity the Ritz slab can place it [km/s] (>= 0): the
minimum over 11 window positions [a, a + w], w = 1/(2 relDepth), around q0 of
max(v_tr(a), v_tr(a + w)), minus vS. `mat` the half-space Material of the ROTATED structure,
`thetaStar` the sagittal root angle of the min-shear entry (q0 from tan, -tan or cot,
whichever reproduces vS).
"""
function slabResolution(mat::Material, vS::Real, thetaStar::Real, relDepth::Real)
    C = standardElastTensor(mat)
    rho = mat.rho
    qs = filter(isfinite, [tan(thetaStar), -tan(thetaStar), cot(thetaStar)])
    q0 = qs[argmin([abs(traceVel(C, rho, q) - vS) for q in qs])]
    w = 1 / (2 * relDepth)
    best = Inf
    for k in 0:10
        a = q0 - w + (k * w) / 10
        best = min(best, max(traceVel(C, rho, a), traceVel(C, rho, a + w)))
    end
    return max(best - vS, 0.0)
end

const _BASIS_WARNED = Ref(false)

"""
    checkBasisHS(calcParams) -> Bool

`legOrder >= HS_MIN_LEG_PER_DEPTH x relDepth`; warns once per session if not.
"""
function checkBasisHS(cp)
    ok = cp.legOrder >= HS_MIN_LEG_PER_DEPTH * cp.relDepth
    if !ok && !_BASIS_WARNED[]
        @warn "Ritz-only half-space classification: legOrder $(cp.legOrder) < " *
              "$HS_MIN_LEG_PER_DEPTH x relDepth ($(cp.relDepth)). The slab does not resolve the " *
              "band up to the longitudinal limit and the calibrated thresholds do not hold. " *
              "Use legOrder >= $(ceil(Int, HS_MIN_LEG_PER_DEPTH * cp.relDepth)). (Shown once per session.)"
        _BASIS_WARNED[] = true
    end
    return ok
end

end # module RitzScore
