# =========================================================================================
# VisResults.jl -- plot-data helpers that read classification results
#
# The Vis helpers that read the Calc RESULTS of the current classification system, plus the
# polarisation-RGB view. Included into `module Visual` after `VisData.jl`.
#
#   greenComponentData           <- @Vis/greenComponentData.m
#   visValidateGreenComponent    <- @Vis/private/visValidateGreenComponent.m
#   visIsLayered                 <- @Vis/private/visIsLayered.m
#   visHasRitzSpectrum           <- @Vis/private/visHasRitzSpectrum.m
#   visRitzSpectrumData          <- @Vis/private/visRitzSpectrumData.m
#   ritzModeScores               <- @Vis/ritzModeScores.m
#   classificationScoreScale     <- @Vis/classificationScoreScale.m
#   visRitzCriterion(+Names)     <- @Vis/private/visRitzCriterion.m
#   visClassificationSizes       <- @Vis/private/visClassificationSizes.m
#   calcPolarisationRgb          <- @Vis/calcPolarisationRgb.m
#   visPolarisationColours       <- @Vis/private/visPolarisationColours.m
#   visIsRgbDirs                 <- @Vis/private/visIsRgbDirs.m
#
# The 'spot' component of greenComponentData (the stored spot-averaged Green curve / map,
# CalcSingleHS.calcSpotGreen / CalcAngSweepHS.calcSpotGreenMap, EXPERIMENTAL as in MATLAB) is
# ported: `greenComponentData(result, "spot")` returns the
# stored |X| on its own `spotGreenVelAxis`, `Vis:greenComponentData:noSpotGreenData` when
# nothing is stored and `:noSpotSupport` on a layered result. The plot entry points
# (plotGreenMap / plotGreenCurve / plotClassificationDiagnostic / plotClassificationCurve) take
# `component = "spot"` unchanged.
# =========================================================================================

const _SingleResult = Union{CalcSingleHS,CalcSingleLayered}
const _SweepResult = Union{CalcAngSweepHS,CalcAngSweepLayered}

"""
    visValidateGreenComponent(component, fnName) -> String

Reject a Green-component argument that is not `"GIJ"` / `"GIJeff"` (I, J in 1..3) or
`"spot"`; returns the component as a `String`. Port of
`@Vis/private/visValidateGreenComponent.m` (error `Vis:<fnName>:badComponent`).
"""
function visValidateGreenComponent(component, fnName::AbstractString)
    c = component isa Union{Symbol,AbstractString} ? String(component) : ""
    occursin(r"^(G[123][123](eff)?|spot)$", c) ||
        throw(ArgumentError("Vis:$fnName:badComponent -- component must be 'GIJ' or 'GIJeff' " *
            "with I, J in 1..3 -- e.g. 'G13eff' (stored, default), 'G13' (one-sided G13_+), " *
            "'G33', or 'spot' (got '$component')."))
    return c
end

"""
    greenComponentData(result, component = "G13eff") -> (absG, componentLabel, velAxis)

|G| data, TeX label and velocity axis of one Green component of a Calc result. Port of
`@Vis/greenComponentData.m`.

`result` -- `CalcSingleHS`/`CalcSingleLayered` (a curve, `absG` a `Vector{Float32}`) or
`CalcAngSweepHS`/`CalcAngSweepLayered` (a map, `[nVel x nAngles]`), calculated with
`calcMode = :full`. `component` -- `"G13eff"` (default: the only component a result stores,
the one the existence score reads), `"GIJ"` (one-sided `|G_IJ(+k)|`) or `"GIJeff"`
(`+-k` symmetrised); anything but `"G13eff"` is RECOMPUTED on the result's own axis
(`greenComponent` / `greenComponentMap`: one Green solve per angle) and not stored.
Labels: `"G_{13}^{eff}"`, `"G_{13}^{+}"`, `"G_{33}"` (diagonal: no side).

Errors: `Vis:greenComponentData:noGreenData` (empty `greenVelAxis`: not `:full`),
`:noG13EffData` (empty stored G13_eff).

`component = "spot"` returns the stored spot-averaged |G| of a half-space result
(`calcSpotGreen` / `calcSpotGreenMap` with a `ProbeKernel`, EXPERIMENTAL) on its own
`spotGreenVelAxis` (works under `:ritzOnly` too), labelled `"G_{13}^{spot}"` (the kernel's
component); errors `:noSpotGreenData` (nothing stored yet) and `:noSpotSupport` (a layered
result).
"""
function greenComponentData(result::Union{_SingleResult,_SweepResult}, component = "G13eff")
    c = visValidateGreenComponent(component === nothing ? "G13eff" : component,
                                  "greenComponentData")
    c == "spot" && return _spotGreenData(result)
    I, J, isEff = parseGreenComponent(c)
    label = I == J ? "G_{$I$J}" : "G_{$I$J}^{$(isEff ? "eff" : "+")}"
    velAxis = result.greenVelAxis
    isempty(velAxis) && throw(ArgumentError("Vis:greenComponentData:noGreenData -- " *
        "greenVelAxis is empty -- the result must be calculated with calcMode=:full to " *
        "populate Green's function data."))
    isSweep = result isa _SweepResult
    if c == "G13eff"
        absG = isSweep ? result.greenAbsG13EffMap : result.greenAbsG13Eff
        isempty(absG) && throw(ArgumentError("Vis:greenComponentData:noG13EffData -- " *
            "$(isSweep ? "greenAbsG13EffMap" : "greenAbsG13Eff") is empty -- the result has " *
            "no Green data (calcMode :ritzOnly or not calculated)."))
    else
        absG = isSweep ? greenComponentMap(result, c) : greenComponent(result, c)
    end
    return absG, label, velAxis
end

# The stored spot-averaged curve / map (calcSpotGreen / calcSpotGreenMap); MATLAB's local
# `spotData`. Checked before the Green-data checks: it does not need calcMode :full.
_spotGreenData(result::CalcAngSweepHS) = _spotGreenData(result, result.spotGreenMap,
    "spotGreenMap", "calcSpotGreenMap(sweepResult, kernel)")
_spotGreenData(result::CalcSingleHS) = _spotGreenData(result, result.spotGreen,
    "spotGreen", "calcSpotGreen(singleResult, kernel)")
_spotGreenData(result) = throw(ArgumentError("Vis:greenComponentData:noSpotSupport -- 'spot' " *
    "needs a half-space result (CalcAngSweepHS / CalcSingleHS); got $(nameof(typeof(result)))."))
function _spotGreenData(result, absG, field, how)
    isempty(absG) && throw(ArgumentError("Vis:greenComponentData:noSpotGreenData -- $field " *
        "is empty -- run $how first (kernel: a ProbeKernel)."))
    label = "G_{$(result.spotGreenInfo.kernel.component[2:3])}^{spot}"
    return absG, label, result.spotGreenVelAxis
end

"""
    visIsLayered(result) -> Bool

`true` for `CalcSingleLayered` / `CalcAngSweepLayered`. Port of `@Vis/private/visIsLayered.m`.
"""
visIsLayered(result) = result isa Union{CalcSingleLayered,CalcAngSweepLayered}

"""
    visHasRitzSpectrum(result) -> Bool

Whether a result stores the Ritz modal numerators |S(v)| is rebuilt from: single angle --
`ritzScoreNum` and `ritzScoreVelAxis` non-empty; sweep -- true for every angle (and at least
one). Port of `@Vis/private/visHasRitzSpectrum.m`.
"""
visHasRitzSpectrum(result::_SingleResult) =
    !isempty(result.ritzScoreNum) && !isempty(result.ritzScoreVelAxis)
visHasRitzSpectrum(result::_SweepResult) =
    !isempty(result.cs) && all(visHasRitzSpectrum, result.cs)

"""
    visRitzSpectrumData(singleResult, velAxis = nothing) -> (absS, velAxis)

The Ritz spectrum `|S(v)| = |RitzScore.spectrum(vels[1:nNum], ritzScoreNum, velAxis)|` of one
angle, rebuilt from the stored modal numerators (no solve); `velAxis` defaults to the
result's own `ritzScoreVelAxis`. Port of `@Vis/private/visRitzSpectrumData.m` -- the one
evaluation shared by the spectrum curve, the spectrum map and the classification curve.
"""
function visRitzSpectrumData(singleResult::_SingleResult, velAxis = nothing)
    num = singleResult.ritzScoreNum
    v = (velAxis === nothing || isempty(velAxis)) ? singleResult.ritzScoreVelAxis :
        collect(Float64, vec(velAxis))
    absS = abs.(RitzScore.spectrum(view(singleResult.ritzSolver.vels, 1:length(num)), num, v))
    return absS, v
end

"""
    ritzModeScores(singleResult) -> Vector{Float64}

Per-mode Ritz score `[nCompModes]`: each Ritz mode gets the score (`RitzScore.candidateScore`,
log10; half-space `dScore`, layered `dScore + polTerm`) of the Ritz candidate it supports
(its `modeNo`), the highest where it supports several, `NaN` where none. Every scored
candidate counts, accepted or not. Port of `@Vis/ritzModeScores.m`.
"""
function ritzModeScores(singleResult::_SingleResult)
    nM = singleResult.ritzSolver.nCompModes
    scores = fill(NaN, nM)
    cands = singleResult.ritzScoreCands
    isempty(cands) && return scores
    sc = RitzScore.candidateScore(cands, visIsLayered(singleResult))
    for j in eachindex(cands)
        m = cands[j].modeNo
        (m < 1 || m > nM) && continue
        if isnan(scores[m]) || sc[j] > scores[m]
            scores[m] = sc[j]
        end
    end
    return scores
end

"""
    classificationScoreScale(method, isLayered) -> NamedTuple

The score scale of a classification method: `(high, low, lowScale, sizeRange, name)`.
`method` -- `:green` | `:refined` | `:ritz` (layered: `:green` | `:ritz`). `high`/`low` are the
classes' own threshold constants; `lowScale` -- `:score` (the `'low'` standard is compared with
the entry score) or `:structureRelative` (half-space `'green'`/`'refined'`: a separate
structure-relative score); `sizeRange` -- the score range the marker size spans,
`[high, max(high, low) + 1]` on the `:score` scale, else `[high, high + 1.5]`; `name` --
`"Green"` | `"Refined"` | `"Ritz"`. Port of `@Vis/classificationScoreScale.m`
(MATLAB's `lowScale` `'structure-relative'` is `:structureRelative`).
"""
function classificationScoreScale(method, isLayered::Bool)
    m = Symbol(method)
    if m === :green
        high, low, lowScale = isLayered ?
            (GreenScoreLayered.THRESHOLD_HIGH, GreenScoreLayered.THRESHOLD_LOW, :score) :
            (GreenScore.THRESHOLD_HIGH, GreenScore.THRESHOLD_LOW, :structureRelative)
    elseif m === :refined
        isLayered && throw(ArgumentError("Vis:classificationScoreScale:noRefinedLayered -- " *
            "Method 'refined' does not exist for layered stacks."))
        high, low, lowScale = GreenSpotScore.THRESHOLD_HIGH, GreenSpotScore.THRESHOLD_LOW,
                              :structureRelative
    elseif m === :ritz
        high, low, lowScale = isLayered ?
            (RitzScore.LAY_THRESHOLD_HIGH, RitzScore.LAY_THRESHOLD_LOW, :score) :
            (RitzScore.HS_THRESHOLD_HIGH, RitzScore.HS_THRESHOLD_LOW, :score)
    else
        throw(ArgumentError("Vis:classificationScoreScale:badMethod -- method must be " *
            "'green', 'refined' or 'ritz'; got '$method'."))
    end
    ids, names = classificationMethods(false)
    name = names[findfirst(==(m), ids)]
    sizeRange = lowScale === :score ? (high, max(high, low) + 1.0) : (high, high + 1.5)
    return (high = Float64(high), low = Float64(low), lowScale = lowScale,
            sizeRange = Float64.(sizeRange), name = name)
end

"""
    visRitzCriterionNames() -> (:ritzScore, :freqDisp, :depthInsensitivity)

The Ritz-mode criteria `plotRitzModeMap` / `plotRitzModeCurve` validate against.
"""
visRitzCriterionNames() = (:ritzScore, :freqDisp, :depthInsensitivity)

"""
    visRitzCriterion(singleResult, criterion) -> (vals, info)

Per-mode values `[nCompModes]` of a Ritz-mode criterion and how to show them. `:ritzScore`
-- `20 * ritzModeScores` in dB (`NaN` = the mode supports no candidate); `:freqDisp` /
`:depthInsensitivity` -- straight from `ritzSolver.diagnostics`. `info = (label, vetoLevel)`:
`vetoLevel` = `RitzScore.MIN_FREQDISP` for `:freqDisp` on a half-space, `NaN` otherwise.
Port of `@Vis/private/visRitzCriterion.m`.
"""
function visRitzCriterion(singleResult::_SingleResult, criterion)
    c = Symbol(criterion)
    ritz = singleResult.ritzSolver
    nM = ritz.nCompModes
    label = String(c); veto = NaN
    if c === :ritzScore
        vals = 20 .* ritzModeScores(singleResult)
        label = "Ritz score [dB]"
    elseif c === :freqDisp
        vals = Float64.(ritz.diagnostics.freqDisp[1:nM])
        visIsLayered(singleResult) || (veto = RitzScore.MIN_FREQDISP)
    elseif c === :depthInsensitivity
        vals = Float64.(ritz.diagnostics.depthInsensitivity[1:nM])
    else
        throw(ArgumentError("Vis:visRitzCriterion:badCriterion -- Unknown criterion '$criterion'."))
    end
    return vals, (label = label, vetoLevel = Float64(veto))
end

"""
    visClassificationSizes(pts, isLayered, minSz, maxSz, defaultSz, limitingSz) -> Vector{Float64}

Marker size of each classified point by its score: per method, the score range
`classificationScoreScale(method, isLayered).sizeRange` maps linearly (clipped) onto
`[minSz, maxSz]`; `defaultSz` where a point has no score (or no method); `limitingSz` for every
`:limiting` point. `pts` is [`visFlattenSignificantModes`](@ref)'s output. Port of
`@Vis/private/visClassificationSizes.m`.
"""
function visClassificationSizes(pts, isLayered::Bool, minSz::Real, maxSz::Real,
                                defaultSz::Real, limitingSz::Real)
    sizeVals = fill(Float64(defaultSz), length(pts.vel))
    for m in unique(pts.method)
        m === :none && continue
        lo, hi = classificationScoreScale(m, isLayered).sizeRange
        for i in eachindex(sizeVals)
            (pts.method[i] === m && !isnan(pts.score[i])) || continue
            frac = clamp((pts.score[i] - lo) / (hi - lo), 0.0, 1.0)
            sizeVals[i] = minSz + (maxSz - minSz) * frac
        end
    end
    sizeVals[pts.type .=== :limiting] .= limitingSz
    return sizeVals
end

# -----------------------------------------------------------------------------------------
# Polarisation-RGB view (dirs = 'rgb')
# -----------------------------------------------------------------------------------------

"""
    visPolarisationColours() -> Matrix{Float64} [3 x 3]

Row i = the RGB colour of a displacement purely along direction i: u1 (propagation) red,
u2 (in-plane transverse) blue, u3 (normal) green darkened to 0.65 for balanced luminance.
Port of `@Vis/private/visPolarisationColours.m` (shared by the panels and the key).
"""
visPolarisationColours() = [1.0 0.0 0.0;
                            0.0 0.0 1.0;
                            0.0 0.65 0.0]

"""
    visIsRgbDirs(dirs, fnName) -> Bool

A modal-shape plot's `dirs`: a vector of directions 1-3 (`false`) or `:rgb`/`"rgb"` (`true`,
the merged polarisation view); anything else throws `Vis:<fnName>:badDirs`. Port of
`@Vis/private/visIsRgbDirs.m`.
"""
function visIsRgbDirs(dirs, fnName::AbstractString)
    if dirs isa Union{Symbol,AbstractString}
        lowercase(String(dirs)) == "rgb" ||
            throw(ArgumentError("Vis:$fnName:badDirs -- dirs must be a vector of directions " *
                "1-3 or 'rgb' (got '$dirs')."))
        return true
    end
    (dirs isa Union{Integer,AbstractVector{<:Integer}} && !isempty(dirs) &&
     all(d -> 1 <= d <= 3, dirs)) ||
        throw(ArgumentError("Vis:$fnName:badDirs -- dirs must be a vector of directions 1-3 or 'rgb'."))
    return false
end

"""
    calcPolarisationRgb(dispField) -> (rgbCell, ampCell)

Merge a displacement-field triplet into one true-colour polarisation image per layer. Port
of `@Vis/calcPolarisationRgb.m` (the UTGS paper's Figs 3(g)/4(b) scheme).

`dispField` -- `(3, nLay)` matrix of signed real 2D fields `[nMeshDepth x nMeshLateral]`
(as from [`calcDispField2D`](@ref) / `calcGreenDispField2D(Layered)`; layers may differ in
row count). Returns `rgbCell::Vector{Array{Float64,3}}` (`[nR x nC x 3]` in [0,1]) and
`ampCell::Vector{Matrix{Float64}}` (`|u| / max|u|`, one normaliser across all layers).
Hue = `p1*red + p2*blue + p3*green` with `p = (|u1|,|u2|,|u3|)/|u|`; brightness = amplitude:
`pixel = (1 - a) + a * hue`, white at zero.
"""
function calcPolarisationRgb(dispField::AbstractMatrix)
    size(dispField, 1) == 3 || throw(ArgumentError("calcPolarisationRgb: dispField must have " *
        "3 rows (one per displacement direction); got $(size(dispField, 1))."))
    nLay = size(dispField, 2)
    base = visPolarisationColours()
    amps = Vector{Matrix{Float64}}(undef, nLay)
    ampMax = 0.0
    for il in 1:nLay
        u1, u2, u3 = dispField[1, il], dispField[2, il], dispField[3, il]
        amps[il] = sqrt.(abs.(u1) .^ 2 .+ abs.(u2) .^ 2 .+ abs.(u3) .^ 2)
        ampMax = max(ampMax, maximum(amps[il]))
    end
    ampMax = max(ampMax, eps())
    rgbCell = Vector{Array{Float64,3}}(undef, nLay)
    ampCell = Vector{Matrix{Float64}}(undef, nLay)
    for il in 1:nLay
        amp = amps[il]
        nR, nC = size(amp)
        rgb = Array{Float64,3}(undef, nR, nC, 3)
        a = amp ./ ampMax
        @inbounds for j in 1:nC, i in 1:nR
            d = max(amp[i, j], eps())
            p1 = abs(dispField[1, il][i, j]) / d
            p2 = abs(dispField[2, il][i, j]) / d
            p3 = abs(dispField[3, il][i, j]) / d
            ai = a[i, j]
            for c in 1:3
                hue = p1 * base[1, c] + p2 * base[2, c] + p3 * base[3, c]
                rgb[i, j, c] = clamp(1 - ai + ai * hue, 0.0, 1.0)
            end
        end
        rgbCell[il] = rgb
        ampCell[il] = a
    end
    return rgbCell, ampCell
end
