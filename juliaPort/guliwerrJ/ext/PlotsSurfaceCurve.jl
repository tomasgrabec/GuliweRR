# =========================================================================================
# PlotsSurfaceCurve.jl -- Makie plotting extension, "surface/curve" group
#
# Real implementations for the 9 stub functions declared in `src/visual/PlotStubs.jl` that
# draw a 3D directional surface over the whole sphere (`plotVelSurface`,
# `plotGroupVelSurface`, `plotYoungsSurface`) or a 2D line/curve over angle or velocity
# (`plotVelCurve`, `plotYoungsCurve`, `plotGreenCurve`, `plotRitzSpectrumCurve`,
# `plotRitzModeCurve`, `plotRitzModeCurvePlate`). Ports of
# `guliwerr/visual/@Vis/plotVelSurface.m` etc. -- see each function's own docstring below for
# the exact MATLAB file it mirrors.
#
# This file is `include`d directly into `module GuliwerrJMakieExt` (see
# `ext/GuliwerrJMakieExt.jl`), which already has `using GuliwerrJ` and `using Makie` in scope
# -- neither is repeated here. Every plotting method below is added to the REAL generic
# function declared in `PlotStubs.jl` via the fully-qualified `function GuliwerrJ.plotXxx(...)`
# form, never a bare `function plotXxx(...)`, which would silently create an unrelated second
# binding instead of extending the stub (see `PlotStubs.jl`'s own header for why this matters).
#
# Local helpers (`_psc_*` prefix, "Plots Surface Curve") are plain top-level functions in this
# same module namespace -- shared with three sibling group files (`PlotsMaps.jl`,
# `PlotsModalShape.jl`, `PlotsDiagnostic.jl`), so names are prefixed to avoid an accidental
# collision; the MATLAB-private helpers this file needs
# (`visWaveSurfaceFigure`, `visWaveSurfacePanels`) have no port in `VisData.jl`: they touch
# graphics handles, so this extension reimplements them in Makie's idiom, not MATLAB's.
#
# Return types match each stub's docstring exactly: `Figure` for the three whole-sphere
# `*Surface` functions (each builds/owns its own multi-panel figure), `Axis` for the
# curve functions (draw into a supplied `ax` or create+return a fresh one via `_psc_new_axis`).
# `plotRitzSpectrumCurve` additionally returns its spectrum data (`(Axis, spec)`), matching
# MATLAB's multi-output signature minus the line handle `hOut`, which none of these Julia
# methods return -- see `PlotStubs.jl`'s documented signatures, none of which carries a
# MATLAB-style graphics-handle output.
#
# Deliberately NOT reproduced from MATLAB (documented deviations, not oversights):
#   - No colour-handle/line-handle outputs (`hOut`) -- Makie plot objects are not part of any
#     of the 9 stub signatures in `PlotStubs.jl`; a caller wanting the underlying plot object
#     can walk `ax.scene.plots`.
#   - `plotVelSurface`/`plotGroupVelSurface`/`plotYoungsSurface` colour with Makie's default
#     `:viridis` colormap rather than MATLAB's default (`parula`, implicit -- the MATLAB
#     source never calls `colormap(...)` for these three); colours have no numeric reference
#     to match (see this extension module's own docstring), and `:viridis` is
#     Makie's own idiomatic sequential default.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# Shared private helpers
# -----------------------------------------------------------------------------------------

"""
    _psc_new_axis(ax) -> Axis

Return `ax` unchanged if it is already an `Axis`; otherwise create a fresh `Figure` and
return a new `Axis` inside it. Shared by every curve-drawing function below -- the Makie
idiom for the MATLAB pattern `if nargin < N || isempty(hAx); figure(); hAx = axes(); end`.
"""
function _psc_new_axis(ax)
    ax !== nothing && return ax
    fig = Makie.Figure()
    return Makie.Axis(fig[1, 1])
end

"""
    _psc_wave_surface_figure() -> Figure

Fresh figure sized for a 1x3 row of 3D wave-surface panels. Port of
`guliwerr/visual/@Vis/private/visWaveSurfaceFigure.m` -- a default `Figure()` is roughly
4:3, which would squeeze a 1x3 row of `Axis3` panels far too narrow next to their colorbars;
widening fixes it, matching MATLAB's own `3 * 480` width choice. Used only when creating a
fresh figure, never to override a caller-supplied `fig`'s size.
"""
_psc_wave_surface_figure() = Makie.Figure(size = (3 * 480, 450))

"""
    _psc_wave_surface_panels(fig, X, Y, Z, C) -> (axs, plts)

Render the qL/qT1/qT2 row of 3D wave-surface panels shared by [`plotVelSurface`](@ref) and
[`plotGroupVelSurface`](@ref). Port of `guliwerr/visual/@Vis/private/visWaveSurfacePanels.m`:
only the surface GEOMETRY differs between those two callers (a scalar phase
velocity/slowness scaling the wave normal vs. the genuine group-velocity vector), so the
panel layout/labels/titles/colour treatment live here once.

`X`,`Y`,`Z`,`C` are `[nPhi x nTheta x 3]` (third dimension = wave type, qL/qT1/qT2 in that
order) -- the direct output shape of `PlotData.velSurfaceData`/`PlotData.groupVelSurfaceData`.
Axis labels name the project's standard frame (`x_3` = surface normal), per
`visStandardFrameTensor`'s own docblock.
"""
function _psc_wave_surface_panels(fig, X::AbstractArray{<:Real,3}, Y::AbstractArray{<:Real,3},
                                   Z::AbstractArray{<:Real,3}, C::AbstractArray{<:Real,3})
    titles = ("qL", "qT1", "qT2")
    axs = Vector{Any}(undef, 3)
    plts = Vector{Any}(undef, 3)
    for k in 1:3
        ax = Makie.Axis3(fig[1, k]; aspect = :data,
                          xlabel = "x_1 (in-plane ref)", ylabel = "x_2 (in-plane)",
                          zlabel = "x_3 (surface normal)", title = titles[k])
        plt = Makie.surface!(ax, X[:, :, k], Y[:, :, k], Z[:, :, k];
                              color = C[:, :, k], colormap = :viridis)
        Makie.Colorbar(fig[2, k], plt; vertical = false)
        axs[k] = ax
        plts[k] = plt
    end
    return axs, plts
end

# -----------------------------------------------------------------------------------------
# plotVelSurface -- port of plotVelSurface.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotVelSurface(mater::Material, surfaceType::Symbol = :velocity,
                              nTheta::Integer = 61, nPhi::Integer = 121; fig = nothing)
        -> Makie.Figure

3D directional bulk-wave velocity (or slowness) surfaces -- one panel per wave type
(qL, qT1, qT2). Port of `guliwerr/visual/@Vis/plotVelSurface.m`; numeric half is
`PlotData.velSurfaceData`. `fig`, if given, is drawn into directly (a fresh 1x3-sized
`Figure` is created otherwise, via [`_psc_wave_surface_figure`](@ref)).
"""
function GuliwerrJ.plotVelSurface(mater::Material, surfaceType::Symbol = :velocity,
                                   nTheta::Integer = 61, nPhi::Integer = 121; fig = nothing)
    _, _, Xs, Ys, Zs, V = velSurfaceData(mater, surfaceType, nTheta, nPhi)
    fig = fig === nothing ? _psc_wave_surface_figure() : fig
    _psc_wave_surface_panels(fig, Xs, Ys, Zs, V)
    return fig
end

# -----------------------------------------------------------------------------------------
# plotGroupVelSurface -- port of plotGroupVelSurface.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotGroupVelSurface(mater::Material, nTheta::Integer = 61, nPhi::Integer = 121;
                                   fig = nothing) -> Makie.Figure

3D group-velocity ("ray") surface -- one panel per wave type. Port of
`guliwerr/visual/@Vis/plotGroupVelSurface.m`; numeric half is `PlotData.groupVelSurfaceData`.
Coloured by group-velocity MAGNITUDE (`mag`), not the surface radius, since the plotted point
need not lie along the wave normal (see that function's own docstring -- a self-intersecting
surface is real physics, not a bug).
"""
function GuliwerrJ.plotGroupVelSurface(mater::Material, nTheta::Integer = 61,
                                        nPhi::Integer = 121; fig = nothing)
    _, _, X, Y, Z, mag = groupVelSurfaceData(mater, nTheta, nPhi)
    fig = fig === nothing ? _psc_wave_surface_figure() : fig
    _psc_wave_surface_panels(fig, X, Y, Z, mag)
    return fig
end

# -----------------------------------------------------------------------------------------
# plotYoungsSurface -- port of plotYoungsSurface.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotYoungsSurface(mater::Material, nTheta::Integer = 61, nPhi::Integer = 121;
                                 fig = nothing) -> Makie.Figure

3D directional Young's-modulus surface -- single panel (one scalar field, not three
wave-type sheets). Port of `guliwerr/visual/@Vis/plotYoungsSurface.m`; numeric half is
`PlotData.youngsSurfaceData`.
"""
function GuliwerrJ.plotYoungsSurface(mater::Material, nTheta::Integer = 61,
                                      nPhi::Integer = 121; fig = nothing)
    _, _, Xs, Ys, Zs, E = youngsSurfaceData(mater, nTheta, nPhi)
    fig = fig === nothing ? Makie.Figure() : fig
    ax = Makie.Axis3(fig[1, 1]; aspect = :data,
                      xlabel = "x_1 (in-plane ref)", ylabel = "x_2 (in-plane)",
                      zlabel = "x_3 (surface normal)", title = "Young's modulus [GPa]")
    plt = Makie.surface!(ax, Xs, Ys, Zs; color = E, colormap = :viridis)
    Makie.Colorbar(fig[1, 2], plt)
    return fig
end

# -----------------------------------------------------------------------------------------
# plotVelCurve -- port of plotVelCurve.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotVelCurve(mater::Material, angleDeg = range(0.0, 360.0, length = 361);
                            ax = nothing) -> Makie.Axis

Bulk-wave phase velocity vs. in-plane angle, at the fixed surface normal in `mater.dirMat` --
one line per wave type (qL, qT1, qT2). Port of `guliwerr/visual/@Vis/plotVelCurve.m`; numeric
half is `PlotData.velCurveData`.
"""
function GuliwerrJ.plotVelCurve(mater::Material, angleDeg = range(0.0, 360.0, length = 361);
                                 ax = nothing)
    angleDeg, vels = velCurveData(mater, angleDeg)
    ax = _psc_new_axis(ax)
    labels = ("qL", "qT1", "qT2")
    for k in 1:3
        Makie.lines!(ax, angleDeg, vels[k, :]; label = labels[k])
    end
    ax.xlabel = "angle [deg]"
    ax.ylabel = "phase velocity [km/s]"
    extAxisLegend!(ax)
    return ax
end

# -----------------------------------------------------------------------------------------
# plotYoungsCurve -- port of plotYoungsCurve.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotYoungsCurve(mater::Material, angleDeg = range(0.0, 360.0, length = 361);
                               ax = nothing) -> Makie.Axis

Directional Young's modulus vs. in-plane angle, at the fixed surface normal in
`mater.dirMat`. Port of `guliwerr/visual/@Vis/plotYoungsCurve.m`; numeric half is
`PlotData.youngsCurveData`.
"""
function GuliwerrJ.plotYoungsCurve(mater::Material, angleDeg = range(0.0, 360.0, length = 361);
                                    ax = nothing)
    angleDeg, E = youngsCurveData(mater, angleDeg)
    ax = _psc_new_axis(ax)
    Makie.lines!(ax, angleDeg, E)
    ax.xlabel = "angle [deg]"
    ax.ylabel = "Young's modulus [GPa]"
    return ax
end

# -----------------------------------------------------------------------------------------
# plotGreenCurve -- port of plotGreenCurve.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotGreenCurve(singleResult::Union{CalcSingleHS,CalcSingleLayered},
                              scaleMode::Symbol = :log; ax = nothing, component = "G13eff")
        -> Makie.Axis

`|G|` vs. velocity for one propagation angle. Port of `@Vis/plotGreenCurve.m`. `component`
-- `"G13eff"` (default, stored) or any `"GIJ"`/`"GIJeff"` (recomputed, see
[`greenComponentData`](@ref)). Throws `Vis:greenComponentData:noGreenData` when the result
has no Green data (`calcMode = :ritzOnly`).
"""
function GuliwerrJ.plotGreenCurve(singleResult::Union{CalcSingleHS,CalcSingleLayered},
                                   scaleMode::Symbol = :log; ax = nothing, component = "G13eff")
    absG, label, velAxis = greenComponentData(singleResult, component)
    visValidateChoice(scaleMode, (:log, :linear), "plotGreenCurve", "scaleMode")
    plotData, yLabel = scaleGreenData(absG, scaleMode, label)
    ax = _psc_new_axis(ax)
    Makie.lines!(ax, velAxis, plotData; label = "|$label|")
    ax.xlabel = "velocity [km/s]"
    ax.ylabel = yLabel
    ax.title = "Green's function curve ($label, angle = $(round(singleResult.angle, digits = 1)) deg)"
    return ax
end

# -----------------------------------------------------------------------------------------
# plotRitzSpectrumCurve -- port of plotRitzSpectrumCurve.m (replaced plotSurfaceResponseCurve)
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotRitzSpectrumCurve(singleResult, scaleMode::Symbol = :log; ax = nothing,
                                     showCands = false) -> (Makie.Axis, spec)

The Ritz spectrum `|S(v)|` of one angle from the stored modal numerators (no solve; works
under `calcMode = :ritzOnly`). `showCands` marks the Ritz candidates (red down-triangles at
the candidate velocity, black circles at the supporting mode's velocity). Returns the axis
and `spec = (velAxis, absS, cands)`. Port of `@Vis/plotRitzSpectrumCurve.m`; numeric half
`ritzSpectrumCurveData`.
"""
function GuliwerrJ.plotRitzSpectrumCurve(singleResult::Union{CalcSingleHS,CalcSingleLayered},
                                          scaleMode::Symbol = :log; ax = nothing,
                                          showCands::Bool = false)
    velAxis, plotData, yLabel, absS, cands = ritzSpectrumCurveData(singleResult, scaleMode)
    ax = _psc_new_axis(ax)
    Makie.lines!(ax, velAxis, plotData; label = "|S|")
    if showCands && !isempty(cands)
        yAt(v) = _psc_interp(velAxis, plotData, v)
        cv = [c.vel for c in cands]; mv = [c.modeVel for c in cands]
        Makie.scatter!(ax, cv, yAt.(cv); color = :red, marker = :dtriangle, markersize = 12)
        Makie.scatter!(ax, mv, yAt.(mv); color = :black, marker = :circle, markersize = 9)
    end
    ax.xlabel = "velocity [km/s]"
    ax.ylabel = yLabel
    ax.title = "Ritz spectrum (angle = $(round(singleResult.angle, digits = 1)) deg)"
    return ax, (velAxis = velAxis, absS = absS, cands = cands)
end

"""Linear interpolation of `y(x)` at `v` (`x` ascending), `NaN` outside -- MATLAB
`interp1(x, y, v, 'linear', NaN)`."""
function _psc_interp(x::AbstractVector, y::AbstractVector, v::Real)
    (isnan(v) || v < x[1] || v > x[end]) && return NaN
    i = clamp(searchsortedlast(x, v), 1, length(x) - 1)
    t = (v - x[i]) / (x[i + 1] - x[i])
    return (1 - t) * y[i] + t * y[i + 1]
end

# -----------------------------------------------------------------------------------------
# plotRitzModeCurve -- port of plotRitzModeCurve.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotRitzModeCurve(singleResult::Union{CalcSingleHS,CalcSingleLayered},
                                 criterion::Symbol = :ritzScore,
                                 velCutoffFactor::Real = 1.10; ax = nothing) -> Makie.Axis

Every computed Ritz mode as a velocity x criterion curve for one angle. Port of
`@Vis/plotRitzModeCurve.m`. `criterion` -- `:ritzScore` (default: the per-mode Ritz score in
dB, [`ritzModeScores`](@ref); modes supporting no candidate drawn as faint grey dots 6 dB below
the lowest score) | `:freqDisp` (with the half-space veto line `RitzScore.MIN_FREQDISP`) |
`:depthInsensitivity` ([`visRitzCriterion`](@ref)). `velCutoffFactor` -- upper velocity
cutoff as a multiple of the longitudinal limiting velocity (`Inf` disables it).
"""
function GuliwerrJ.plotRitzModeCurve(singleResult::Union{CalcSingleHS,CalcSingleLayered},
                                      criterion::Symbol = :ritzScore,
                                      velCutoffFactor::Real = 1.10; ax = nothing)
    visValidateChoice(criterion, visRitzCriterionNames(), "plotRitzModeCurve", "criterion")
    if isfinite(velCutoffFactor)
        longLimVel = visPropOrEmpty(singleResult, (:longLimVel, :subLongLimVel))
        longLimVel === nothing &&
            throw(ArgumentError("Vis:plotRitzModeCurve:noLongVel -- longitudinal limiting " *
                "velocity is unavailable (calcMode was not :full) -- cannot apply " *
                "velCutoffFactor; pass Inf to disable it."))
        velCutoff = velCutoffFactor * longLimVel
    else
        velCutoff = Inf
    end
    ritz = singleResult.ritzSolver
    nM = ritz.nCompModes
    vels = ritz.vels[1:nM]
    critVal, info = visRitzCriterion(singleResult, criterion)
    keep = vels .<= velCutoff
    velKeep = vels[keep]; critKeep = critVal[keep]
    order = sortperm(velKeep)
    velSort = velKeep[order]; critSort = critKeep[order]

    ax = _psc_new_axis(ax)
    if criterion === :ritzScore
        faint = isnan.(critSort)
        yBase = all(faint) ? 0.0 : minimum(critSort[.!faint]) - 6
        blue = Makie.RGBf(0, 0.4470, 0.7410)
        Makie.lines!(ax, velSort[.!faint], critSort[.!faint]; color = blue)
        Makie.scatter!(ax, velSort[.!faint], critSort[.!faint]; color = blue, markersize = 6)
        any(faint) && Makie.scatter!(ax, velSort[faint], fill(yBase, count(faint));
                                      color = Makie.RGBf(0.7, 0.7, 0.7), markersize = 4)
    else
        Makie.lines!(ax, velSort, critSort)
        isnan(info.vetoLevel) || Makie.hlines!(ax, [info.vetoLevel];
                                                color = Makie.RGBf(0.8, 0.1, 0.1),
                                                linestyle = :dash)
    end
    ax.xlabel = "velocity [km/s]"
    ax.ylabel = info.label
    ax.title = "Ritz mode curve (angle = $(round(singleResult.angle, digits = 1)) deg) -- " *
               info.label
    return ax
end

# -----------------------------------------------------------------------------------------
# plotRitzModeCurvePlate -- port of plotRitzModeCurvePlate.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotRitzModeCurvePlate(singleResult::CalcSinglePlate,
                                      criterion::Symbol = :modeType,
                                      velCutoffFactor::Real = 1.10; ax = nothing) -> Makie.Axis

Plate-geometry counterpart of [`plotRitzModeCurve`](@ref): every computed Ritz mode for a
`CalcSinglePlate` result as a velocity x criterion scatter, coloured by mode type. Port of
`guliwerr/visual/@Vis/plotRitzModeCurvePlate.m`; numeric half is `visPlateRitzModeData`.

`criterion` -- `:modeType` (default, categorical y-axis: Asymm/Symm/TH-Asymm/TH-Symm, ticks
labelled via [`plateModeTypeStyle`](@ref)) | `:confidence` (`modeTypeConfidence` in `[0,1]`).
Marker colour always encodes mode type either way. `velCutoffFactor` -- upper velocity cutoff
as a multiple of the structure's own longitudinal limiting velocity (computed fresh via
[`visPlateLongLimVel`](@ref), since `CalcSinglePlate` stores no `longLimVel`); pass `Inf` to
disable.
"""
function GuliwerrJ.plotRitzModeCurvePlate(singleResult::CalcSinglePlate,
                                           criterion::Symbol = :modeType,
                                           velCutoffFactor::Real = 1.10; ax = nothing)
    visValidateChoice(criterion, (:modeType, :confidence), "plotRitzModeCurvePlate", "criterion")

    velCutoff = isfinite(velCutoffFactor) ?
        velCutoffFactor * visPlateLongLimVel(singleResult.structure, nothing) : Inf

    vels, types, conf = visPlateRitzModeData(singleResult.ritzSolver)
    keep = vels .<= velCutoff

    velKeep = vels[keep]
    typeKeep = types[keep]
    confKeep = conf[keep]
    order = sortperm(velKeep)
    velSort = velKeep[order]
    typeSort = typeKeep[order]
    confSort = confKeep[order]

    typeNames, typeColors = plateModeTypeStyle()
    colSort = [Makie.RGBf(typeColors[t, 1], typeColors[t, 2], typeColors[t, 3]) for t in typeSort]
    ySort = criterion === :modeType ? Float64.(typeSort) : confSort

    ax = _psc_new_axis(ax)
    Makie.scatter!(ax, velSort, ySort; color = colSort, markersize = 10)
    ax.xlabel = "velocity [km/s]"
    if criterion === :modeType
        ax.ylabel = "mode type"
        ax.yticks = (1:4, collect(typeNames))
        Makie.ylims!(ax, 0.5, 4.5)
    else
        ax.ylabel = "mode-type confidence"
        Makie.ylims!(ax, 0, 1)
    end
    ax.title = "Ritz mode curve (angle = $(round(singleResult.angle, digits = 1)) deg) -- " *
               string(criterion)
    return ax
end
