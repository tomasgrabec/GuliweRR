# =========================================================================================
# PlotsDiagnostic.jl -- Makie plotting extension, "classification/dispersion/
# inversion diagnostic" group (7 functions): plotClassificationCurve,
# plotClassificationDiagnostic, plotDispersionCurves, plotDispersionCurvesPlate,
# plotAssociationDiagnostic, plotInversionFit, plotSensitivitySpectrum.
#
# `include`d directly into `module GuliwerrJMakieExt` (see GuliwerrJMakieExt.jl) -- `using
# GuliwerrJ` and `using Makie` are already in scope from that file, so every exported
# GuliwerrJ name (visValidateChoice, visFlattenSignificantModes, dispersionCurvesData,
# ColormapRYG, ...) is reachable unqualified here. Every method below is added to the
# existing generic function via the fully-qualified `function GuliwerrJ.plotXxx(...)` form,
# per PlotStubs.jl's own instruction -- never a bare `function plotXxx(...)`.
#
# DEVIATIONS FROM THE MATLAB SOURCE (ported for semantics, not pixel-for-pixel style; see
# PlotStubs.jl -- exact MATLAB styling is not required):
#
#   1. GREEN BACKGROUND DRAWN INLINE, NOT VIA plotGreenCurve/plotGreenMap. MATLAB's
#      plotClassificationCurve.m/plotClassificationDiagnostic.m call `Vis.plotGreenCurve`/
#      `Vis.plotGreenMap` for their background. Those two stubs are owned by a DIFFERENT
#      group file (PlotsSurfaceCurve.jl/PlotsMaps.jl) -- calling them
#      here would make this file depend on another group's method being loaded
#      (a `MethodError` otherwise), purely for a background element. Instead the
#      background is drawn directly here from `scaleGreenData` (VisData.jl), using the SAME data and the SAME log-normalised semantics MATLAB's own
#      plotGreenCurve/plotGreenMap use -- functionally identical, just not delegated.
#   2. NO PRIVATE-HELPER NAME COLLISION. `visClassificationScatterByType`,
#      `visDispersionCurveAxes` and `visGreenColourAxis` are MATLAB private
#      plotting helpers with no Julia port of their own. Their logic is
#      reproduced below as `_diag`-prefixed local functions/constants (`_diagScatterByType!`,
#      `_diagDispersionCurveAxes`, `_diagGreenColourAxis!`) rather than under their MATLAB
#      names, since this file is `include`d into a MODULE SHARED with the other three group
#      files -- a same-named top-level function/constant defined independently in more than
#      one of those files would silently overwrite the other (Julia has no per-file
#      namespacing inside one `include`d module), so every helper here is deliberately
#      prefixed to be collision-proof regardless of what the other groups choose to name
#      theirs.
#   3. COLORBARS ONLY WHEN THIS FUNCTION OWNS THE FIGURE. Makie's `Colorbar` needs a
#      `GridPosition` to place itself into; there is no general way to find "the empty cell
#      next to an arbitrary caller-supplied `Axis`" the way MATLAB's `colorbar(hAx)` manages
#      automatically. So a colorbar is added only when `ax`/`fig` was `nothing` (this
#      function created its own `Figure`) -- when the caller supplies an existing axis for
#      overlay/composition, no colorbar is added (the caller already owns that figure's
#      layout). Every other visual element (points, lines, labels) is drawn either way.
#   4. RETURN TYPE matches PlotStubs.jl exactly (`Axis`, or `(Axis, Axis)` for
#      plotSensitivitySpectrum) -- MATLAB's richer multi-output signatures
#      (`[hPoints, hAx, hBg, hLim]`, ...) are not reproduced; a caller wanting the individual
#      plot objects can query `ax.scene.plots` (as the structural checks of
#      these functions do).
#   5. plotDispersionCurves'/plotDispersionCurvesPlate's Green/heatmap x-axis is sorted
#      ascending before calling `heatmap!` -- Makie's `heatmap!(xs, ys, zs)` assumes
#      monotonically increasing coordinate vectors, unlike `pcolorCenteredGrid`'s explicit
#      contract supporting a decreasing axis (`:k_vphase`'s `2*pi ./ wavelengths` is
#      decreasing for increasing wavelengths). Sorting `xs` (and reindexing `zs`'s matching
#      dimension) changes nothing physically -- a heatmap's cell order carries no information
#      beyond position.
#   6. plotSensitivitySpectrum's heat-map: Makie's default y-axis increases upward (row 1 at
#      the bottom); MATLAB's `imagesc` puts row 1 at the top. The cell CONTENTS/labels are
#      identical either way (this function's `text!` annotations are placed at the same
#      (column, row) coordinate the heatmap cell occupies), so this is a vertical mirror of
#      MATLAB's own layout, not a data error.
#   7. LEGEND CALLS mirror MATLAB's own choice per function, not blanket-applied:
#      plotClassificationCurve/Diagnostic (colorBy=:type) and plotAssociationDiagnostic/
#      plotInversionFit show a legend (matching `legend(hAx, ...)` in the MATLAB source);
#      plotDispersionCurves/plotDispersionCurvesPlate do NOT (MATLAB never calls `legend`
#      there either, only tags each line's `DisplayName`/`UserData` for a caller's own use).
# =========================================================================================

# -----------------------------------------------------------------------------------------
# Shared private helpers (collision-proof `_diag` prefix -- see header point 2)
# -----------------------------------------------------------------------------------------

const _DIAG_TYPE_COLORS = Dict(
    :Rayleigh       => Makie.RGBf(0.0, 0.4470, 0.7410),
    :pseudoRayleigh => Makie.RGBf(0.8500, 0.3250, 0.0980),
    :SAW            => Makie.RGBf(0.8500, 0.3250, 0.0980),
    :limiting       => Makie.RGBf(0.4660, 0.6740, 0.1880),
    :other          => Makie.RGBf(0.5, 0.5, 0.5),
)
const _DIAG_TYPE_MARKERS = Dict(
    :Rayleigh       => :circle,
    :pseudoRayleigh => :diamond,
    :SAW            => :utriangle,
    :limiting       => :xcross,
    :other          => :circle,
)
const _DIAG_TYPE_DRAW_ORDER = (:other, :limiting, :pseudoRayleigh, :Rayleigh, :SAW)

# MATLAB's default `lines(7)` categorical palette, cycled beyond 7 -- used by
# plotDispersionCurvesPlate's :branch colour mode.
const _DIAG_LINES_CMAP = Makie.RGBf[
    Makie.RGBf(0.0000, 0.4470, 0.7410), Makie.RGBf(0.8500, 0.3250, 0.0980),
    Makie.RGBf(0.9290, 0.6940, 0.1250), Makie.RGBf(0.4940, 0.1840, 0.5560),
    Makie.RGBf(0.4660, 0.6740, 0.1880), Makie.RGBf(0.3010, 0.7450, 0.9330),
    Makie.RGBf(0.6350, 0.0780, 0.1840)]

"""
    _diagRgbVector(cmap::AbstractMatrix{<:Real}) -> Vector{Makie.RGBf}

Convert one of `Colormaps.jl`'s `N x 3` RGB matrices (`ColMapGrayOrange`/`ColormapRWB`/
`ColormapRYG`) into a Makie colour vector usable as a `colormap =` value.
"""
_diagRgbVector(cmap::AbstractMatrix{<:Real}) = [Makie.RGBf(cmap[i, 1], cmap[i, 2], cmap[i, 3])
                                                  for i in axes(cmap, 1)]

"""
    _diagInterp1(x, y, xq) -> Vector{Float64}

Piecewise-linear interpolation of `y(x)` at `xq`, linearly EXTRAPOLATED beyond `x`'s own
range (MATLAB's `interp1(x, y, xq, 'linear', 'extrap')`). `x` assumed sorted ascending
(every `greenVelAxis` in this codebase is built that way).
"""
function _diagInterp1(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                       xq::AbstractVector{<:Real})
    n = length(x)
    out = zeros(Float64, length(xq))
    for (k, xv) in pairs(xq)
        i = if xv <= x[1]
            1
        elseif xv >= x[end]
            n - 1
        else
            clamp(searchsortedlast(x, xv), 1, n - 1)
        end
        t = (xv - x[i]) / (x[i + 1] - x[i])
        out[k] = y[i] + t * (y[i + 1] - y[i])
    end
    return out
end

"""
    _diagScatterByType!(ax, xVals, yVals, sizeVals, allTypes) -> Vector{Symbol}

Port of `guliwerr/visual/@Vis/private/visClassificationScatterByType.m`: one scatter series
per classification `type` actually present, in the FIXED draw order (`:other`/`:limiting`
first so a genuine mode point drawn later is never hidden underneath a limiting-velocity
reference marker), with a legend. Returns the list of types actually drawn (draw order).
"""
function _diagScatterByType!(ax, xVals::AbstractVector{<:Real}, yVals::AbstractVector{<:Real},
                              sizeVals::AbstractVector{<:Real}, allTypes::AbstractVector{Symbol},
                              isWeak::AbstractVector{Bool} = falses(length(allTypes)))
    uTypesPresent = unique(allTypes)
    orderedKnown = [t for t in _DIAG_TYPE_DRAW_ORDER if t in uTypesPresent]
    orderedRest = [t for t in uTypesPresent if !(t in _DIAG_TYPE_DRAW_ORDER)]
    uTypes = vcat(orderedKnown, orderedRest)

    for t in uTypes
        mask = allTypes .== t
        col = get(_DIAG_TYPE_COLORS, t, Makie.RGBf(0.3, 0.3, 0.3))
        # 'weak' entries (accepted at the looser standard only) are hollow, as MATLAB's
        # zero face alpha with an opaque edge
        faces = [w ? Makie.RGBAf(0, 0, 0, 0) : Makie.RGBAf(col.r, col.g, col.b, 1) for w in isWeak[mask]]
        scatter!(ax, xVals[mask], yVals[mask]; markersize = sizeVals[mask], color = faces,
                 strokecolor = col, strokewidth = 1,
                 marker = get(_DIAG_TYPE_MARKERS, t, :circle), label = String(t))
    end
    if !isempty(uTypes)
        try
            extAxisLegend!(ax)
        catch
        end
    end
    return uTypes
end

"""
    _diagDispersionCurveAxes(curveType, wavelength, k, vel, velGroup, freq,
                              wavelengthLabel, kLabel, freqLabel) -> (xAll, yAll, xLabel, yLabel)

Port of `guliwerr/visual/@Vis/private/visDispersionCurveAxes.m`: map a `curveType` onto the
x/y data and axis labels of a dispersion-curve plot. Shared by [`plotDispersionCurves`](@ref)
(Layered) and [`plotDispersionCurvesPlate`](@ref) (Plate) -- each caller does its own physics
and hands the finished `[nWav x nCompModes]` arrays plus its OWN axis-label strings here
(the Layered/Plate sides mean slightly different things by "k"/"wavelength"/"frequency" --
see each caller).
"""
function _diagDispersionCurveAxes(curveType::Symbol, wavelength::AbstractVector{<:Real},
                                   k::AbstractVector{<:Real}, vel::AbstractMatrix{<:Real},
                                   velGroup::AbstractMatrix{<:Real}, freq::AbstractMatrix{<:Real},
                                   wavelengthLabel::AbstractString, kLabel::AbstractString,
                                   freqLabel::AbstractString)
    nCompModes = size(vel, 2)
    VPHASE_LABEL = "v_φ [km/s]"
    VGROUP_LABEL = "v_g [km/s]"

    if curveType === :wavelength_freq
        return repeat(wavelength, 1, nCompModes), freq, wavelengthLabel, freqLabel
    elseif curveType === :wavelength_vphase
        return repeat(wavelength, 1, nCompModes), vel, wavelengthLabel, VPHASE_LABEL
    elseif curveType === :k_freq
        return repeat(k, 1, nCompModes), freq, kLabel, freqLabel
    elseif curveType === :k_vphase
        return repeat(k, 1, nCompModes), vel, kLabel, VPHASE_LABEL
    elseif curveType === :freq_vphase
        return freq, vel, freqLabel, VPHASE_LABEL
    elseif curveType === :freq_vgroup
        return freq, velGroup, freqLabel, VGROUP_LABEL
    else
        throw(ArgumentError("_diagDispersionCurveAxes:badCurveType -- unknown curveType " *
            ":$curveType"))
    end
end

"""
    _diagConfidenceToRGB(vals, cmap) -> Vector{Makie.RGBf}

Port of `plotAssociationDiagnostic.m`'s local `localConfidenceToRGB`: map confidence values
in `[0,1]` to explicit RGB rows via a Colormaps.jl `N x 3` matrix, for use as a per-point
scatter colour spec (overlay mode). `NaN` -> neutral mid-grey (never an colormap extreme).
"""
function _diagConfidenceToRGB(vals::AbstractVector{<:Real}, cmap::AbstractMatrix{<:Real})
    n = size(cmap, 1)
    out = Vector{Makie.RGBf}(undef, length(vals))
    for (i, v) in pairs(vals)
        if isnan(v)
            out[i] = Makie.RGBf(0.6, 0.6, 0.6)
        else
            idx = clamp(round(Int, 1 + clamp(v, 0.0, 1.0) * (n - 1)), 1, n)
            out[i] = Makie.RGBf(cmap[idx, 1], cmap[idx, 2], cmap[idx, 3])
        end
    end
    return out
end

"""
    _diagMatchCalculatedVelocity(pts, angleList, velMeasList) -> Vector{Float64}

Port of `plotInversionFit.m`'s local `localMatchCalculatedVelocity`: for each
`(angle, measuredVelocity)` pair, find the matching `AssocPoint` in a DIFFERENT
`ModeAssociation` snapshot of the same points and return its `calculatedVelocity` (`NaN` if
no match).
"""
function _diagMatchCalculatedVelocity(pts, angleList::AbstractVector{<:Real},
                                       velMeasList::AbstractVector{<:Real})
    velOut = fill(NaN, length(angleList))
    ptsAngle = [p.angleDeg for p in pts]
    ptsVelMeas = [p.measuredVelocity for p in pts]
    ptsVelCalc = [p.calculatedVelocity for p in pts]
    for k in eachindex(angleList)
        onAngle = abs.(ptsAngle .- angleList[k]) .< 1e-6
        match = onAngle .& (abs.(ptsVelMeas .- velMeasList[k]) .< 1e-9)
        if any(match)
            velOut[k] = ptsVelCalc[findfirst(match)]
        end
    end
    return velOut
end

"""
    _diagCheckAngleIdx(angleIdx, nAngles, fname)

Shared bounds check for `angleIdx` (mirrors MATLAB's `validateattributes(angleIdx, ...,
{'scalar','integer','positive','<=',numel(sweepResult.angles)}, ...)`).
"""
function _diagCheckAngleIdx(angleIdx::Integer, nAngles::Integer, fname::AbstractString)
    (1 <= angleIdx <= nAngles) || throw(ArgumentError(
        "$fname:badAngleIdx -- angleIdx must be an integer in 1:$nAngles, got $angleIdx."))
end

# -----------------------------------------------------------------------------------------
# plotClassificationCurve -- Vis.plotClassificationCurve.m
# -----------------------------------------------------------------------------------------

"""Marker diameter (Makie px) from a MATLAB scatter area (pt^2)."""
_diagMarkerDiam(area) = 1.5 .* sqrt.(area)

"""Default classification background: G13_eff, or the Ritz spectrum for a `:ritzOnly` result."""
_diagDefaultComponent(result, component) = component !== nothing ? String(component) :
    (isempty(result.greenVelAxis) && visHasRitzSpectrum(result)) ? "ritzSpectrum" : "G13eff"

function _diagMissingBackground(component, isSweep)
    component == "spot" && return "No spot-averaged Green $(isSweep ? "map" : "curve") " *
        "(run $(isSweep ? "calcSpotGreenMap" : "calcSpotGreen") first)"
    component == "ritzSpectrum" &&
        return "No Ritz spectrum (this result stores no Ritz modal numerators)"
    return "No $component $(isSweep ? "map" : "curve") (calcMode was not :full)"
end

"""Whether the stored spot-averaged curve / map is there (`calcSpotGreen` / `calcSpotGreenMap`;
a layered result has none)."""
_diagHasSpot(result::CalcSingleHS) = !isempty(result.spotGreen)
_diagHasSpot(result::CalcAngSweepHS) = !isempty(result.spotGreenMap)
_diagHasSpot(result) = false

"""Score-coloured scatter (black = low, white = high). Weak entries: MATLAB draws them hollow
(face alpha 0); a colormapped Makie scatter has no per-point face alpha, so they get a thick
grey edge instead -- one series either way."""
function _diagScatterByScore!(ax, fig, col, x, y, sizes, pts, colorBy, isWeak)
    cVals, cLabel = visClassificationColorVals(colorBy, pts.modeVel, pts.scoreDb, pts.greenPeakVel)   # modeVel: vel IS greenPeakVel
    finiteVals = filter(isfinite, cVals)
    crange = isempty(finiteVals) ? (0.0, 1.0) : (minimum(finiteVals), maximum(finiteVals))
    crange = crange[1] == crange[2] ? (crange[1], crange[1] + 1.0) : crange
    sc = scatter!(ax, x, y; markersize = sizes, color = cVals, colormap = :grays,
                  colorrange = crange, strokewidth = [w ? 2.0 : 0.5 for w in isWeak],
                  strokecolor = [w ? Makie.RGBf(0.55, 0.55, 0.55) : Makie.RGBf(0, 0, 0) for w in isWeak])
    fig !== nothing && Colorbar(fig[1, col], sc; label = cLabel)
    return sc
end

"""
    GuliwerrJ.plotClassificationCurve(singleResult, colorBy = :type, showGreenBackground = true,
                                       refLevel = :full; ax = nothing, component = nothing)
        -> Axis

The classified entries of one angle on its response curve. Port of
`Vis.plotClassificationCurve`. Background `component`: `"G13eff"` (default
with Green data), any `"GIJ"`/`"GIJeff"`, or `"ritzSpectrum"` (default for a `:ritzOnly`
result). Marker size by score on the entry's method scale ([`visClassificationSizes`](@ref));
`'weak'` entries hollow; `colorBy` -- `:type` | `:score` (`scoreDb`) | `:confidence`.
"""
function GuliwerrJ.plotClassificationCurve(singleResult, colorBy::Symbol = :type,
                                            showGreenBackground::Bool = true,
                                            refLevel::Symbol = :full; ax = nothing,
                                            component = nothing)
    visValidateChoice(colorBy, (:type, :score, :confidence), "plotClassificationCurve", "colorBy")
    visValidateChoice(refLevel, (:sawOnly, :sawLimit, :full), "plotClassificationCurve", "refLevel")
    comp = _diagDefaultComponent(singleResult, component)
    comp == "ritzSpectrum" || visValidateGreenComponent(comp, "plotClassificationCurve")

    fig = nothing
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1])
    end

    hasBg = comp == "ritzSpectrum" ? visHasRitzSpectrum(singleResult) :
            comp == "spot" ? _diagHasSpot(singleResult) : !isempty(singleResult.greenVelAxis)
    if hasBg
        absCurve, label, velAxis = comp == "ritzSpectrum" ?
            (visRitzSpectrumData(singleResult)[1], "S", singleResult.ritzScoreVelAxis) :
            greenComponentData(singleResult, comp)
        curveYData, yLabelStr = scaleGreenData(absCurve, :log, label)
    end
    if showGreenBackground
        if !hasBg
            @warn "Vis:plotClassificationCurve:noGreenBackground -- " *
                  "$(_diagMissingBackground(comp, false)) -- skipping the background."
        else
            lines!(ax, velAxis, curveYData; color = (:gray20, 0.8), linewidth = 1.5,
                   label = "|$label|")
        end
    end

    pts = visFlattenSignificantModes([singleResult.classification.significantModes], [NaN],
                                      visClassificationTypesAllowed(refLevel))
    allYPos = hasBg ? _diagInterp1(velAxis, curveYData, pts.vel) : zeros(Float64, length(pts.vel))
    isLayered = visIsLayered(singleResult)
    sizeVals = _diagMarkerDiam(visClassificationSizes(pts, isLayered, 30, 140, 60, 26))
    isWeak = pts.grade .=== :weak

    if colorBy === :type
        _diagScatterByType!(ax, pts.vel, allYPos, sizeVals, pts.type, isWeak)
    else
        _diagScatterByScore!(ax, fig, 2, pts.vel, allYPos, sizeVals, pts, colorBy, isWeak)
    end

    if refLevel !== :sawOnly
        if !isLayered
            longLimVel = singleResult.longLimVel
            isnan(longLimVel) || vlines!(ax, [longLimVel]; color = Makie.RGBf(0.9290, 0.6940, 0.1250),
                                          linewidth = 1.5, label = "longitudinal limit")
        else
            shearVels = visPropOrEmpty(singleResult, (:subShearLimVels,))
            if shearVels !== nothing && !isempty(shearVels)
                vlines!(ax, shearVels; color = Makie.RGBf(0.4660, 0.6740, 0.1880),
                        linestyle = :dash, label = "shear limit")
            end
        end
    end

    ax.xlabel = "velocity [km/s]"
    hasBg && (ax.ylabel = yLabelStr)
    ax.title = "Classification curve (angle = $(round(singleResult.angle, digits = 1)) deg)"
    return ax
end

# -----------------------------------------------------------------------------------------
# plotClassificationDiagnostic -- Vis.plotClassificationDiagnostic.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotClassificationDiagnostic(sweepResult, colorBy = :type,
                                            showGreenBackground = true, refLevel = :full;
                                            ax = nothing, component = nothing) -> Axis

The classified entries of an angular sweep over its response map. Port of
`Vis.plotClassificationDiagnostic`: background [`plotGreenMap`](@ref)
(`component`, default `"G13eff"`) or [`plotRitzSpectrumMap`](@ref) (`"ritzSpectrum"`, default
for a `:ritzOnly` sweep); marker size by score on the method scale, `'weak'` entries hollow.
"""
function GuliwerrJ.plotClassificationDiagnostic(sweepResult, colorBy::Symbol = :type,
                                                 showGreenBackground::Bool = true,
                                                 refLevel::Symbol = :full; ax = nothing,
                                                 component = nothing)
    visValidateChoice(colorBy, (:type, :score, :confidence), "plotClassificationDiagnostic",
                       "colorBy")
    visValidateChoice(refLevel, (:sawOnly, :sawLimit, :full), "plotClassificationDiagnostic",
                       "refLevel")
    comp = _diagDefaultComponent(sweepResult, component)
    comp == "ritzSpectrum" || visValidateGreenComponent(comp, "plotClassificationDiagnostic")

    fig = nothing
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1])
    end

    if showGreenBackground
        noBg = comp == "ritzSpectrum" ? !visHasRitzSpectrum(sweepResult) :
               comp == "spot" ? !_diagHasSpot(sweepResult) : isempty(sweepResult.greenVelAxis)
        if noBg
            @warn "Vis:plotClassificationDiagnostic:noGreenBackground -- " *
                  "$(_diagMissingBackground(comp, true)) -- skipping the background."
        elseif comp == "ritzSpectrum"
            GuliwerrJ.plotRitzSpectrumMap(sweepResult; ax = ax)
        else
            GuliwerrJ.plotGreenMap(sweepResult; ax = ax, component = comp)
        end
    end

    angles = sweepResult.angles
    pts = visFlattenSignificantModes(sweepResult.significantModesMap, angles,
                                      visClassificationTypesAllowed(refLevel))
    isLayered = visIsLayered(sweepResult)
    sizeVals = _diagMarkerDiam(visClassificationSizes(pts, isLayered, 12, 100, 36, 16))
    isWeak = pts.grade .=== :weak

    if colorBy === :type
        _diagScatterByType!(ax, pts.angle, pts.vel, sizeVals, pts.type, isWeak)
    else
        _diagScatterByScore!(ax, fig, 3, pts.angle, pts.vel, sizeVals, pts, colorBy, isWeak)
    end

    if refLevel !== :sawOnly
        if !isLayered
            longVel = sweepResult.longVelMap
            if !isempty(longVel) && !all(isnan, longVel)
                lines!(ax, angles, longVel; color = Makie.RGBf(0.9290, 0.6940, 0.1250),
                       linewidth = 1.5, label = "longitudinal limit")
            end
        else
            shearAngles = Float64[]; shearVels = Float64[]
            for ia in eachindex(angles)
                sv = visPropOrEmpty(sweepResult.cs[ia], (:subShearLimVels,))
                (sv === nothing || isempty(sv)) && continue
                append!(shearAngles, fill(angles[ia], length(sv)))
                append!(shearVels, sv)
            end
            if !isempty(shearVels)
                scatter!(ax, shearAngles, shearVels; marker = :cross, markersize = 9,
                         color = Makie.RGBf(0.4660, 0.6740, 0.1880), label = "shear limit")
            end
        end
    end

    ax.xlabel = "angle [deg]"
    ax.ylabel = "velocity [km/s]"
    return ax
end

# -----------------------------------------------------------------------------------------
# plotDispersionCurves -- Vis.plotDispersionCurves.m (CalcWavevecAngSweep2Layered)
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotDispersionCurves(sweepResult::CalcWavevecAngSweep2Layered,
                                         curveType::Symbol = :wavelength_freq,
                                         angleIdx::Integer = 1,
                                         showRitzBranches::Bool = true,
                                         showGreenBackground::Bool = false,
                                         scaleMode::Symbol = :log,
                                         climRange = (-4.0, -1.0); ax = nothing)
    visValidateChoice(curveType, visDispersionCurveTypes(), "plotDispersionCurves", "curveType")
    _diagCheckAngleIdx(angleIdx, length(sweepResult.angles), "plotDispersionCurves")
    visValidateChoice(scaleMode, (:log, :linear), "plotDispersionCurves", "scaleMode")

    VPHASE_CURVE_TYPES = (:wavelength_vphase, :k_vphase)
    if showGreenBackground && !(curveType in VPHASE_CURVE_TYPES)
        throw(ArgumentError("plotDispersionCurves:backgroundNeedsVphase -- " *
            "showGreenBackground=true requires curveType :wavelength_vphase or :k_vphase " *
            "(velocity y-axis) -- got :$curveType."))
    end

    fig = nothing
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1])
    end

    if showGreenBackground
        xAxis, velCommon, plotData, cbLabel = wavevecGreenBackgroundData(sweepResult, angleIdx,
                                                                          curveType, scaleMode)
        ord = sortperm(xAxis)
        colorKw = scaleMode === :log ? (; colorrange = (Float64(climRange[1]), Float64(climRange[2]))) :
                                        NamedTuple()
        hm = heatmap!(ax, xAxis[ord], velCommon, permutedims(plotData[:, ord]);
                       colormap = _diagRgbVector(ColMapGrayOrange()), colorKw...)
        if fig !== nothing
            Colorbar(fig[1, 2], hm; label = cbLabel)
        end
    end

    branchModeNoMap, wavelengths, k, vel, velGroup, freq, highlightModeNos1 =
        dispersionCurvesData(sweepResult, angleIdx)
    nCompModes = size(branchModeNoMap, 2)

    xAll, yAll, xLabelStr, yLabelStr = _diagDispersionCurveAxes(
        curveType, wavelengths, k, vel, velGroup, freq,
        "wavelength", "k [2π/wavelength]", "frequency (= v_φ / wavelength)")

    for im in 1:nCompModes
        highlighted = branchModeNoMap[1, im] in highlightModeNos1
        col = highlighted ? Makie.RGBf(0.80, 0.20, 0.02) : Makie.RGBf(0.6, 0.6, 0.6)
        lw = highlighted ? 2.0 : 0.5
        lines!(ax, xAll[:, im], yAll[:, im]; color = col, linewidth = lw,
               visible = showRitzBranches)
    end

    ax.xlabel = xLabelStr
    ax.ylabel = yLabelStr
    return ax
end

# -----------------------------------------------------------------------------------------
# plotDispersionCurvesPlate -- Vis.plotDispersionCurvesPlate.m (CalcWavevecAngSweep2Plate)
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotDispersionCurvesPlate(sweepResult::CalcWavevecAngSweep2Plate,
                                              curveType::Symbol = :wavelength_freq,
                                              angleIdx::Integer = 1,
                                              colorMode::Symbol = :type,
                                              confidenceThresh::Real = 0.3,
                                              showChristoffelBounds::Union{Bool,Nothing} = nothing;
                                              ax = nothing)
    visValidateChoice(curveType, visDispersionCurveTypes(), "plotDispersionCurvesPlate", "curveType")
    _diagCheckAngleIdx(angleIdx, length(sweepResult.angles), "plotDispersionCurvesPlate")
    visValidateChoice(colorMode, (:type, :branch), "plotDispersionCurvesPlate", "colorMode")

    VPHASE_CURVE_TYPES = (:wavelength_vphase, :k_vphase, :freq_vphase)
    showBounds = showChristoffelBounds === nothing ? (curveType in VPHASE_CURVE_TYPES) :
                                                       showChristoffelBounds

    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1])
    end

    classification = classifyBranches(sweepResult, angleIdx, confidenceThresh)
    branchModeNoMap = classification.branchModeNoMap
    nWav, nCompModes = size(branchModeNoMap)
    wavVecs = sweepResult.wavVecs

    vel = zeros(Float64, nWav, nCompModes)
    velGroup = zeros(Float64, nWav, nCompModes)
    freq = zeros(Float64, nWav, nCompModes)
    for iw in 1:nWav, im in 1:nCompModes
        rawIdx = branchModeNoMap[iw, im]
        vel[iw, im] = sweepResult.velsMap3D[rawIdx, angleIdx, iw]
        velGroup[iw, im] = sweepResult.velsGroupMap3D[rawIdx, angleIdx, iw]
        freq[iw, im] = sweepResult.freqsPhysicalMap3D[rawIdx, angleIdx, iw]
    end
    wavelength = 2 * pi ./ wavVecs

    xAll, yAll, xLabelStr, yLabelStr = _diagDispersionCurveAxes(
        curveType, wavelength, wavVecs, vel, velGroup, freq, "wavelength", "k", "frequency")

    typeNames, typeColorsMat, parityNames, parityColorsMat = plateModeTypeStyle()
    typeColors = _diagRgbVector(typeColorsMat)
    parityColors = _diagRgbVector(parityColorsMat)

    for im in 1:nCompModes
        isTypeMode = colorMode === :type
        typeNo = classification.modesTypeByBranch[1, im]
        lineStyle = :solid
        if isTypeMode && classification.trustLabelByBranch[1, im]
            col = typeColors[typeNo]
            dispName = typeNames[typeNo]
        elseif isTypeMode && classification.isExactMirrorSymmetric
            # parity exact, Lamb/TH not: the family colour, dashed
            col = parityColors[typeNo]
            dispName = parityNames[typeNo]
            lineStyle = :dash
        else
            col = _DIAG_LINES_CMAP[mod1(im, length(_DIAG_LINES_CMAP))]
            dispName = "branch $im"
        end
        lines!(ax, xAll[:, im], yAll[:, im]; color = col, linewidth = 1.5, linestyle = lineStyle,
               label = dispName)
    end

    ax.xlabel = xLabelStr
    ax.ylabel = yLabelStr

    if showBounds
        structure = sweepResult.structure
        angleDeg = sweepResult.angles[angleIdx]
        longLimBot, shearLimBot = visPlateLimVels(structure, 1, angleDeg)
        hlines!(ax, [shearLimBot]; color = (:gray30, 0.9), linestyle = :dash)
        hlines!(ax, [longLimBot]; color = (:gray30, 0.9), linestyle = :dot)
        if nlay(structure) > 1
            longLimTop, shearLimTop = visPlateLimVels(structure, nlay(structure), angleDeg)
            sameSurface = abs(shearLimTop - shearLimBot) < 1e-9 && abs(longLimTop - longLimBot) < 1e-9
            if !sameSurface
                hlines!(ax, [shearLimTop]; color = (:gray60, 0.9), linestyle = :dash)
                hlines!(ax, [longLimTop]; color = (:gray60, 0.9), linestyle = :dot)
            end
        end
    end
    return ax
end

# -----------------------------------------------------------------------------------------
# plotAssociationDiagnostic -- Vis.plotAssociationDiagnostic.m
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotAssociationDiagnostic(association::ModeAssociation,
                                              overlayMode::Bool = false; ax = nothing)
    fig = nothing
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1])
    end

    pts = association.points
    macTracked = [p.bindingClass in ("ritzMode", "pseudoSAW") for p in pts]

    if any(!, macTracked)
        pNA = pts[.!macTracked]
        scatter!(ax, [p.angleDeg for p in pNA], [p.calculatedVelocity for p in pNA];
                 markersize = 11, color = Makie.RGBf(0.55, 0.55, 0.55), marker = :diamond,
                 strokewidth = 1.25, strokecolor = :black,
                 label = "limiting (Christoffel-routed, no MAC) ($(length(pNA)))")
    end

    pMAC = pts[macTracked]
    angleDeg = [p.angleDeg for p in pMAC]
    velCalc = [p.calculatedVelocity for p in pMAC]
    macValue = [p.macValue for p in pMAC]
    collided = [p.modeCollision for p in pMAC]

    ryg = ColormapRYG()
    sc = nothing
    if overlayMode
        macColour = _diagConfidenceToRGB(macValue, ryg)
        scatter!(ax, angleDeg, velCalc; markersize = 9, color = macColour,
                 strokewidth = 0.5, strokecolor = :black, label = "Ritz-tracked, by MAC")
    else
        sc = scatter!(ax, angleDeg, velCalc; markersize = 10, color = macValue,
                       colormap = _diagRgbVector(ryg), colorrange = (0.0, 1.0),
                       label = "Ritz-tracked, coloured by MAC")
    end
    if any(collided)
        scatter!(ax, angleDeg[collided], velCalc[collided]; markersize = 16,
                 color = Makie.RGBAf(0, 0, 0, 0), marker = :circle, strokewidth = 1.5,
                 strokecolor = :red, label = "collision / no match")
    end
    nCollided = count(collided)

    if overlayMode
        return ax
    end

    if fig !== nothing && sc !== nothing
        Colorbar(fig[1, 2], sc; label = "MAC value (Ritz-tracked points only)")
    end

    ax.xlabel = "propagation angle [deg]"
    ax.ylabel = "calculated velocity [km/s]"
    ax.title = nCollided > 0 ?
        "Association diagnostic ($nCollided point(s) ringed = collision/no match)" :
        "Association diagnostic (no collisions)"
    try
        extAxisLegend!(ax; position = :rb)
    catch
    end
    return ax
end

# -----------------------------------------------------------------------------------------
# plotInversionFit -- Vis.plotInversionFit.m
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotInversionFit(measured::MeasuredVelocities, associationFinal::ModeAssociation,
                                     associationInitial::Union{ModeAssociation,Nothing} = nothing;
                                     ax = nothing)
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1])
    end

    pts = associationFinal.points
    assocAngle = [p.angleDeg for p in pts]
    assocVelMeas = [p.measuredVelocity for p in pts]

    isAssociatedMask = falses(measured.nPoints)
    for i in 1:measured.nPoints
        onAngle = abs.(assocAngle .- measured.angleDeg[i]) .< 1e-6
        isAssociatedMask[i] = any(onAngle) &&
            any(abs.(assocVelMeas[onAngle] .- measured.velocity[i]) .< 1e-9)
    end

    markerOf = Dict("ritzMode" => :circle, "limiting" => :rect, "pseudoSAW" => :utriangle)

    if any(.!isAssociatedMask)
        nDropped = count(.!isAssociatedMask)
        scatter!(ax, measured.angleDeg[.!isAssociatedMask], measured.velocity[.!isAssociatedMask];
                 markersize = 11, color = Makie.RGBf(0.6, 0.6, 0.6), marker = :xcross,
                 strokewidth = 1.5, label = "measured, dropped ($nDropped)")
    end

    bindingClass = [p.bindingClass for p in pts]
    classes = unique(bindingClass)
    for cls in classes
        mask = bindingClass .== cls
        marker = get(markerOf, cls, :circle)
        angleK = assocAngle[mask]
        velMeasK = assocVelMeas[mask]

        scatter!(ax, angleK, velMeasK; markersize = 9, color = :black, marker = marker,
                 label = "measured ($cls)")

        if associationInitial !== nothing
            velGuess = _diagMatchCalculatedVelocity(associationInitial.points, angleK, velMeasK)
            scatter!(ax, angleK, velGuess; markersize = 8, color = Makie.RGBf(0.8500, 0.3250, 0.0980),
                     marker = marker, label = "initial guess ($cls)")
        end

        velRecovered = [p.calculatedVelocity for p in pts[mask]]
        scatter!(ax, angleK, velRecovered; markersize = 8, color = :red, marker = marker,
                 label = "recovered ($cls)")
    end

    ax.xlabel = "propagation angle [deg]"
    ax.ylabel = "velocity [km/s]"
    ax.title = "Inversion fit -- $(count(isAssociatedMask)) of $(measured.nPoints) measured " *
               "point(s) associated"
    try
        extAxisLegend!(ax; position = :rb, merge = true)
    catch
    end
    return ax
end

# -----------------------------------------------------------------------------------------
# plotSensitivitySpectrum -- Vis.plotSensitivitySpectrum.m
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotSensitivitySpectrum(sensitivityAnalysis::SensitivityAnalysis; fig = nothing)
    f = fig === nothing ? Figure() : fig

    L = sensitivityAnalysis.nIndepConsts
    roots = sensitivityAnalysis.normalisedSensitivityRoots
    P = sensitivityAnalysis.nOptimisedConstants
    T = sensitivityAnalysis.transformMatrix

    posList = IndepVoigtPositions(sensitivityAnalysis.symmClassNo)
    constLabels = ["c$(rc[1])$(rc[2])" for rc in posList]

    axBar = Axis(f[1, 1]; yscale = log10, xlabel = "transformed combination p",
                 ylabel = "(√ψ_p)_norm = √(ψ_p/ψ_1)",
                 title = "Sensitivity spectrum ($P of $L optimised)")
    barColors = [i <= P ? Makie.RGBf(0.0, 0.4470, 0.7410) : Makie.RGBf(0.7, 0.7, 0.7) for i in 1:L]
    barplot!(axBar, 1:L, max.(roots, eps()); color = barColors)
    hlines!(axBar, [sensitivityAnalysis.sensitivityCutoff]; color = :red, linestyle = :dash)
    axBar.xticks = 1:L

    axHeat = Axis(f[1, 2]; xlabel = "transformed combination p",
                  ylabel = "independent Voigt constant",
                  title = "T -- signed composition of each combination", aspect = 1)
    maxAbsT = max(maximum(abs, T; init = 0.0), eps())
    hm = heatmap!(axHeat, 1:L, 1:L, permutedims(T); colormap = _diagRgbVector(ColormapRWB()),
                  colorrange = (-maxAbsT, maxAbsT))
    Colorbar(f[1, 3], hm)
    axHeat.xticks = 1:L
    axHeat.yticks = (1:L, constLabels)

    for iRow in 1:L, jCol in 1:L
        val = T[iRow, jCol]
        textColor = abs(val) / maxAbsT > 0.6 ? :white : :black
        text!(axHeat, jCol, iRow; text = string(round(val; digits = 2)),
              align = (:center, :center), color = textColor, fontsize = 10)
    end

    return axBar, axHeat
end
