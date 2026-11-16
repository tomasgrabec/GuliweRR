# =========================================================================================
# PlotsMaps.jl -- Makie plotting extension, "maps/heatmaps" group
#
# Real Makie implementations for 7 of the 26 `src/visual/PlotStubs.jl` stubs: plotGreenMap,
# plotDispersionMap, plotPolarHeatmap, plotAngularMap, plotRitzModeMap, plotRitzModeMapPlate,
# plotRitzSpectrumMap -- ports of `guliwerr/visual/@Vis/plot*.m` of the same names. This
# file is `include`d directly into `module GuliwerrJMakieExt` (see that file's own header),
# which already has `using GuliwerrJ`/`using Makie` in scope -- no repeated `using` here.
# Every method below is added to the REAL stub function via the fully-qualified
# `function GuliwerrJ.plotXxx(...)`, per the stub file's own instruction, never a bare
# `function plotXxx(...)` (which would silently create an unrelated second binding).
#
# All private helpers below are prefixed `_maps` (not just `_`): this file shares ONE Julia
# module namespace with the other three plotting groups
# (PlotsSurfaceCurve.jl/PlotsModalShape.jl/PlotsDiagnostic.jl), so a generically-named helper
# (`_newax`, `_finite_extrema`, ...) would be a real collision risk even though the FILES
# themselves never overlap.
#
# ---------------------------------------------------------------------------------------
# pcolorCenteredGrid vs. Makie's own heatmap! centring convention
# ---------------------------------------------------------------------------------------
# `pcolorCenteredGrid` exists because MATLAB's `pcolor(x,y,C)` treats `C(:,j)` as the cell
# spanning `x(j):x(j+1)` (a lower-left-CORNER convention), not centred on `x(j)` -- see the
# `pcolorCenteredGrid` docstring in VisData.jl and docs/pitfalls.md NUM3.
#
# A tiny 3x3 check (x=[1,2,3], y=[10,20,30], heatmap!(ax,x,y,C) then reading back
# the plot object's OWN converted x/y via `hm[1][]`/`hm[2][]`) shows Makie's `heatmap!`
# ALREADY performs exactly `pcolorCenteredGrid`'s own edge construction internally, whenever
# `length(x) == size(C,1)` (sample CENTRES, not edges): interior edges at consecutive
# midpoints, endpoints linearly extrapolated by the nearest interior half-gap -- e.g.
# x=[1,2,3] -> edges [0.5,1.5,2.5,3.5], bit-for-bit what `_local_edges` computes. The same holds
# on a DECREASING axis ([3,2,1] -> edges [3.5,2.5,1.5,0.5]) and a DEGENERATE
# single-sample axis ([5.0] -> edges [4.5,5.5]) -- both match `_local_edges` exactly, so
# `pcolorCenteredGrid`'s own "monotonic, not necessarily increasing" contract carries over
# for free. Passing `length(x) == size(C,1)+1` values instead makes Makie treat them as
# edges VERBATIM (no further shift) -- i.e. Makie's `heatmap!` cleanly supports BOTH
# conventions and auto-detects which one it was given from the length alone.
#
# CONCLUSION: **`pcolorCenteredGrid` must NOT be called before a Cartesian `heatmap!`** --
# passing the raw sample coordinates (angles, velAxis, kParRange, ...) straight to
# `heatmap!` already reproduces the exact centring `pcolorCenteredGrid` was built to force
# onto MATLAB's `pcolor`. The one real adaptation needed is a TRANSPOSE: this port's
# established data convention (`PlotData.jl`/`VisData.jl`, matching MATLAB's own
# `pcolor(x,y,C)`/`meshgrid(x,y)`) is `C[i,j] <-> (x[j], y[i])` (row = y-index, column =
# x-index; e.g. `ritzSpectrumMapData`'s `plotData` is `[nVel x nAngles]`), whereas Makie's
# `heatmap!(ax,x,y,z)` convention is `z[i,j] <-> (x[i], y[j])` (row = x-index, column =
# y-index) -- the OPPOSITE orientation. Every Cartesian heatmap call below therefore passes
# `permutedims(plotData)`, not `plotData` itself, with no other correction needed.
#
# `pcolorCenteredGrid` is NOT dead code for this group, though: `plotPolarHeatmap` still
# calls it directly (see below) -- it is exactly the right tool there, just not for the
# reason MATLAB used it.
#
# ---------------------------------------------------------------------------------------
# DEVIATION FROM THE MATLAB SOURCE -- polar maps use a hand-built quad mesh, NOT PolarAxis,
# for the same underlying reason plotGreenMap.m/plotPolarHeatmap.m give for their own choice
# ---------------------------------------------------------------------------------------
# MATLAB's `plotPolarHeatmap.m` docstring is explicit: hAx is "ALWAYS a plain axes(), never
# polaraxes(), for both plotTypes" -- because MATLAB's own `pcolor` cannot draw on a
# `polaraxes()` in any current release, so MATLAB manually converts (angle,radius) EDGES to
# Cartesian (X,Y) via `pol2cart` and draws a normal `pcolor(X,Y,CPadded)` on a plain axes,
# then hand-draws a substitute polar grid (`drawPolarGrid`, rings+spokes+labels).
#
# Makie has a native `PolarAxis`, and a check found the SAME split MATLAB's own
# docstring already flags: `scatter!`/`lines!` on a `PolarAxis` render correctly (verified),
# but `heatmap!`/`surface!` on a `PolarAxis` do NOT -- in this Makie/CairoMakie build
# (Makie 0.24.15 / CairoMakie 0.15.15), `heatmap!(polarAxis, theta, r, C)` silently renders a
# garbled/near-empty plot (constructs without error, but the drawn shape does not match the
# data at all -- reproduced on two independent test cases), and `surface!(ax, X, Y, zeros(...);
# color=C, shading=NoShading)` (the natural curvilinear-mesh alternative) throws inside
# CairoMakie's own mesh rasteriser (`MethodError: no method matching iterate(::RGBA{Float32})`
# in `CairoMakie.draw_pattern`). Both are almost certainly genuine upstream limitations/bugs,
# not a one-off misuse (reproduced on independent test cases).
#
# So: `plotPolarHeatmap` reproduces MATLAB's OWN workaround one-for-one, for the SAME reason
# (a native polar heatmap primitive is unreliable in both toolchains) -- `pcolorCenteredGrid`
# converts (angle,radius) to EDGE coordinates, each resulting quadrilateral cell is built as
# an explicit `Makie.Polygon` in Cartesian (x,y) = (r*cosd(a), r*sind(a)) space, and the whole
# batch is drawn with ONE `poly!(ax, polys; color = cellvalues, colormap, colorrange)` call
# (verified correct on a half-disc test case, and fast enough at ~91,000 cells: ~9s total
# including JIT warmup, dominated by first-call compilation). `plotAngularMap`'s polar branch,
# by contrast, is a plain SCATTER of discrete points (mirroring MATLAB's own `polarplot`, not
# `pcolor`) -- exactly the primitive that DOES work natively on `PolarAxis`, verified
# separately -- so it uses a real `PolarAxis`, matching MATLAB's `polaraxes()` choice there.
# This means `plotPolarHeatmap`/`plotGreenMap`/`plotDispersionMap`'s polar branches always
# return a plain `Makie.Axis` (never `PolarAxis`), while `plotAngularMap`'s polar branch
# returns a `PolarAxis` -- exactly mirroring MATLAB's own dual axes()/polaraxes() split
# (plotAngularMap.m's own docstring: "axes() for cartesian, polaraxes() for polar").
# =========================================================================================

const _MAPS_GREEN_CMAP = [RGBf(row[1], row[2], row[3]) for row in eachrow(ColMapGrayOrange())]

"""
    _maps_finite_extrema(A) -> (Float64, Float64)

`extrema` over the finite entries of `A` only (skips `NaN`/`Inf`, matching
[`scaleGreenData`](@ref)'s own `NaN`-means-"blank" contract used throughout this group).
Falls back to `(0.0, 1.0)` if nothing finite remains, and widens a degenerate `lo == hi` by
`0.5` either way so `Colorbar`/`colorrange` never receive a zero-width range.
"""
function _maps_finite_extrema(A)
    vals = filter(isfinite, vec(A))
    isempty(vals) && return (0.0, 1.0)
    lo, hi = extrema(vals)
    lo == hi && return (Float64(lo) - 0.5, Float64(hi) + 0.5)
    return (Float64(lo), Float64(hi))
end

"""
    _maps_newax(ax; xlabel="", ylabel="", title="") -> Axis

Shared "use `ax` if given, else make a fresh `Figure`+`Axis`" helper for this group's
Cartesian plots, setting axis labels/title either way (matching MATLAB's own
`xlabel(hAx,...)`/`ylabel(hAx,...)` calls, which run regardless of whether `hAx` was
freshly created or supplied by the caller).
"""
function _maps_newax(ax; xlabel = "", ylabel = "", title = "")
    if ax === nothing
        fig = Figure()
        ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, title = title)
    else
        ax.xlabel = xlabel
        ax.ylabel = ylabel
        ax.title = title
    end
    return ax
end

"""
    _maps_key_cells(ax) -> (right, below)

The layout cells immediately to the right of and below `ax` (found as in
[`_maps_colorbar_beside`](@ref)), where a map's colour bar and legend go, after deleting a
`Colorbar` left in the right cell or a `Legend` left in the one below by an earlier draw -- so
a viewer switching between a colour-mapped and a legend-keyed map never keeps both.
"""
function _maps_key_cells(ax)
    function cells()
        gc = ax.layoutobservables.gridcontent[]
        gc === nothing && return ax.parent[1, end + 1], ax.parent[end + 1, 1]
        return (gc.parent[gc.span.rows, gc.span.cols.stop + 1],
                gc.parent[gc.span.rows.stop + 1, gc.span.cols])
    end
    right, below = cells()
    oldKeys = [filter(c -> c isa Colorbar, contents(right)); filter(c -> c isa Legend, contents(below))]
    isempty(oldKeys) && return right, below
    foreach(delete!, oldKeys)
    # the emptied row/column would otherwise keep its share of the space (the map shrank to a
    # third of the viewer after Confidence -> Mode type); trimming can renumber the cells
    gc = ax.layoutobservables.gridcontent[]
    gc === nothing || trim!(gc.parent)
    return cells()
end

"""
    _maps_legend_below(ax) -> Legend

Horizontal `Legend` of the labelled plots on `ax`, in the cell below it. For maps whose data
reach every corner, where an `axislegend` covers points (MATLAB's `'Location','best'` has no
Makie equivalent). Below rather than beside: in a side-by-side figure the cell to the right
is the next panel.
"""
function _maps_legend_below(ax)
    _, below = _maps_key_cells(ax)
    return Legend(below, ax; orientation = :horizontal, tellheight = true, tellwidth = false)
end

"""
    _maps_colorbar_beside(ax; colormap, colorrange, label="") -> Colorbar

Places a standalone `Colorbar` in the column immediately to the right of `ax` inside its
OWN parent `GridLayout`, whether `ax` was just created by [`_maps_newax`](@ref) (column 2 of
a fresh 1x1 figure) or supplied by the caller already sitting somewhere inside a larger
multi-panel `Figure` (found via `ax.layoutobservables.gridcontent[]`, which every `Axis`
placed via `fig[i,j]` populates) -- verified empirically on both a freshly-created axis and
one nested at `fig[2,3]` of a bigger layout.

A `Colorbar` already in that cell (a viewer redrawing into the same axis) is deleted first,
so redraws never stack bars. The bar follows later changes to the colour range of the first
colour-mapped `Heatmap`/`Scatter` on `ax` (a viewer's colour fields, `setAxisRange!(…, :c, …)`):
built with fixed `limits`, it would otherwise keep showing the range it was drawn with.
"""
function _maps_colorbar_beside(ax; colormap, colorrange, label = "")
    pos, _ = _maps_key_cells(ax)
    cb = Colorbar(pos; colormap = colormap, limits = colorrange, label = label)
    for p in ax.scene.plots
        (p isa Makie.Heatmap || p isa Makie.Scatter) || continue
        p.colorrange[] isa Makie.Automatic && continue
        on(p.colorrange) do r
            cb.limits[] = (Float64(r[1]), Float64(r[2]))
        end
        break
    end
    return cb
end

"""
    _maps_polar_grid!(ax, angles, radius) -> ax

Hand-drawn radial rings + angular spokes + labels, the decoration a real `PolarAxis` gives
for free -- needed here because [`plotPolarHeatmap`](@ref) draws on a plain `Axis` (see this
file's header for why). Port of `guliwerr/visual/@Vis/private/visPcolorCentered.m`'s sibling
local function `drawPolarGrid` (`plotPolarHeatmap.m`): 4 rings (skipping r=0), spokes every
30 degrees across the DATA's own angular extent (so a half-disc sweep draws a half-disc grid,
not a full circle implying data that isn't there).
"""
function _maps_polar_grid!(ax, angles, radius)
    rMin, rMax = minimum(radius), maximum(radius)
    rMax <= rMin && return ax   # degenerate single-radius input -- nothing sensible to draw

    angMin, angMax = minimum(angles), maximum(angles)
    thArc = range(angMin, angMax, length = 200)
    gridColor = RGBf(0.55, 0.55, 0.55)

    nRings = 4
    ringRadii = collect(range(rMin, rMax, length = nRings + 1))[2:end]   # skip r=0
    for r in ringRadii
        lines!(ax, r .* cosd.(thArc), r .* sind.(thArc); color = gridColor, linewidth = 0.75)
        text!(ax, r * cosd(thArc[end]), r * sind(thArc[end]); text = string(round(r, sigdigits = 3)),
              color = gridColor, fontsize = 10, align = (:left, :top))
    end

    angStep = 30.0
    angTicks = collect(angMin:angStep:angMax)
    (isempty(angTicks) || angTicks[end] != angMax) && push!(angTicks, angMax)
    for a in angTicks
        lines!(ax, [0.0, rMax * cosd(a)], [0.0, rMax * sind(a)]; color = gridColor, linewidth = 0.75)
        text!(ax, 1.06 * rMax * cosd(a), 1.06 * rMax * sind(a); text = string(round(a, digits = 3), "°"),
              color = gridColor, fontsize = 10, align = (:center, :center))
    end
    return ax
end

# -----------------------------------------------------------------------------------------
# plotPolarHeatmap -- port of guliwerr/visual/@Vis/plotPolarHeatmap.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotPolarHeatmap(angles, radius, C; ax = nothing,
                                colormap = :viridis, colorrange = nothing) -> Axis

Generic (angle, radius, colour) map drawn in polar form on a plain Cartesian `Axis` -- see
this file's header for why `PolarAxis` is not used here. `angles` [deg], `radius`
(nonnegative), `C` is `[length(radius) x length(angles)]`, already scaled by the caller
(colour scaling is this function's caller's job, matching MATLAB's own docstring: "Colour
scaling..., colormap, clim and colorbar are left to the CALLER").

`colormap`/`colorrange` are extra optional keywords beyond the stub's documented
`(angles, radius, C; ax=nothing)` -- an additive, backward-compatible extension (not part of
the ported MATLAB signature, which has no notion of a Makie colormap object) that lets a
direct standalone call render sensibly on its own (`colorrange` defaults to the data's own
finite extrema) AND lets [`plotGreenMap`](@ref)/[`plotDispersionMap`](@ref) pass through the
EXACT colour treatment their own Cartesian branch uses, for a consistent look between the two
`plotType`s.
"""
function GuliwerrJ.plotPolarHeatmap(angles::AbstractVector{<:Real}, radius::AbstractVector{<:Real},
                                     C::AbstractMatrix{<:Real}; ax = nothing,
                                     colormap = :viridis, colorrange = nothing)
    all(isfinite, angles) || throw(ArgumentError("plotPolarHeatmap: angles must be finite"))
    all(isfinite, radius) || throw(ArgumentError("plotPolarHeatmap: radius must be finite"))
    all(>=(0), radius) || throw(ArgumentError("plotPolarHeatmap: radius must be nonnegative"))
    size(C) == (length(radius), length(angles)) ||
        throw(ArgumentError("plotPolarHeatmap:incorrectSize -- C must be " *
            "$(length(radius))x$(length(angles)), got $(size(C))"))

    angEdges, radEdges, _ = pcolorCenteredGrid(angles, radius, C)
    radEdges = max.(radEdges, 0.0)   # an extrapolated inner edge can dip below the innermost
                                       # sample; negative radius is meaningless here

    nR = length(radius)
    nA = length(angles)
    polys = Vector{Makie.Polygon}(undef, nR * nA)
    cols = Vector{Float64}(undef, nR * nA)
    k = 0
    for ia in 1:nA, ir in 1:nR
        k += 1
        a0, a1 = angEdges[ia], angEdges[ia + 1]
        r0, r1 = radEdges[ir], radEdges[ir + 1]
        polys[k] = Makie.Polygon(Point2f[
            (r0 * cosd(a0), r0 * sind(a0)),
            (r1 * cosd(a0), r1 * sind(a0)),
            (r1 * cosd(a1), r1 * sind(a1)),
            (r0 * cosd(a1), r0 * sind(a1)),
        ])
        cols[k] = C[ir, ia]
    end

    cr = colorrange === nothing ? _maps_finite_extrema(cols) : colorrange
    axOut = _maps_newax(ax)
    poly!(axOut, polys; color = cols, colormap = colormap, colorrange = cr, strokewidth = 0)
    axOut.aspect = DataAspect()
    hidedecorations!(axOut)
    hidespines!(axOut)
    _maps_polar_grid!(axOut, angles, radius)
    return axOut
end

# -----------------------------------------------------------------------------------------
# plotGreenMap -- port of guliwerr/visual/@Vis/plotGreenMap.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotGreenMap(sweepResult, scaleMode = :log, climRange = (-4.0, -1.0),
                            plotType = :cartesian; ax = nothing, component = "G13eff")
        -> (Axis, Colorbar)

Green's-function angular map -- `|G13_eff|` (default; any `"GIJ"`/`"GIJeff"` is recomputed,
see [`greenComponentData`](@ref)) heatmap, angle x velocity. Port of `Vis.plotGreenMap`
(as MATLAB: standard clim `[-4, -1]`). Throws `Vis:greenComponentData:noGreenData` for a
`:ritzOnly` sweep. `plotType` -- `:cartesian` (default) | `:polar`
([`plotPolarHeatmap`](@ref), same colour treatment).
"""
function GuliwerrJ.plotGreenMap(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered},
                                 scaleMode::Symbol = :log,
                                 climRange = (-4.0, -1.0),
                                 plotType::Symbol = :cartesian;
                                 ax = nothing, component = "G13eff")
    absMap, label, velAxis = greenComponentData(sweepResult, component)
    visValidateChoice(scaleMode, (:log, :linear), "plotGreenMap", "scaleMode")
    visValidateChoice(plotType, (:cartesian, :polar), "plotGreenMap", "plotType")

    plotData, cbLabel = scaleGreenData(absMap, scaleMode, label)
    colorrange = scaleMode === :log ? (Float64(climRange[1]), Float64(climRange[2])) :
                 _maps_finite_extrema(plotData)

    axOut = if plotType === :cartesian
        a = _maps_newax(ax; xlabel = "angle [deg]", ylabel = "velocity [km/s]")
        heatmap!(a, sweepResult.angles, velAxis, permutedims(plotData);
                 colormap = _MAPS_GREEN_CMAP, colorrange = colorrange)
        a
    else
        GuliwerrJ.plotPolarHeatmap(sweepResult.angles, velAxis, plotData;
                                    ax = ax, colormap = _MAPS_GREEN_CMAP, colorrange = colorrange)
    end

    cb = _maps_colorbar_beside(axOut; colormap = _MAPS_GREEN_CMAP, colorrange = colorrange,
                                label = cbLabel)
    return axOut, cb
end

# -----------------------------------------------------------------------------------------
# plotDispersionMap -- port of guliwerr/visual/@Vis/plotDispersionMap.m
# -----------------------------------------------------------------------------------------

"""
    _maps_isofreq_heatmap(slice, scaleMode, climRange, plotType, ax) -> Axis

Shared angle x kParRange Green `|G|` heatmap, used both directly (`CalcIsofreqAngSweepLayered`)
and for one selected omega0 slice of `CalcFreqAngSweep2Layered`. Port of
`guliwerr/visual/@Vis/plotDispersionMap.m`'s local function `plotIsofreqAngSweepHeatmap`.
"""
function _maps_isofreq_heatmap(slice, scaleMode::Symbol, climRange, plotType::Symbol, ax)
    visValidateChoice(scaleMode, (:log, :linear), "plotDispersionMap", "scaleMode")
    visValidateChoice(plotType, (:cartesian, :polar), "plotDispersionMap", "plotType")

    plotData, cbLabel = scaleGreenData(slice.absGMap, scaleMode, "G_{33}")
    colorrange = scaleMode === :log ? (Float64(climRange[1]), Float64(climRange[2])) :
                 _maps_finite_extrema(plotData)

    axOut = if plotType === :cartesian
        a = _maps_newax(ax; xlabel = "angle [deg]", ylabel = "kParRange")
        heatmap!(a, slice.angles, slice.kParRange, permutedims(plotData);
                 colormap = _MAPS_GREEN_CMAP, colorrange = colorrange)
        a
    else
        GuliwerrJ.plotPolarHeatmap(slice.angles, slice.kParRange, plotData;
                                    ax = ax, colormap = _MAPS_GREEN_CMAP, colorrange = colorrange)
    end

    _maps_colorbar_beside(axOut; colormap = _MAPS_GREEN_CMAP, colorrange = colorrange,
                           label = cbLabel)
    return axOut
end

"""
    GuliwerrJ.plotDispersionMap(sweepResult::CalcWavevecAngSweep2Layered, omegaIdx = 1,
                                 scaleMode = :log, climRange = (-5.0, -2.0),
                                 plotType = :cartesian; ax = nothing) -> Axis

isoWavevec 2D case: scatter of classified SAW points (angle, wavelength), coloured by LOCAL
RANK within each point's own cluster (`sawVelMap2D[iw,ia]` is already velocity-sorted
ascending, so rank is position in that array). These are ranks, not tracked branches -- for
connected branch lines use `plotDispersionCurves`. No Green magnitude here, so
`scaleMode`/`climRange`/`plotType` are ignored: always a Cartesian scatter, matching MATLAB.
"""
function GuliwerrJ.plotDispersionMap(sweepResult::CalcWavevecAngSweep2Layered,
                                      omegaIdx::Integer = 1, scaleMode::Symbol = :log,
                                      climRange = (-5.0, -2.0), plotType::Symbol = :cartesian;
                                      ax = nothing)
    nWav, nAngles = size(sweepResult.sawVelMap2D)
    allAng = Float64[]
    allWav = Float64[]
    allRank = Float64[]
    for iw in 1:nWav, ia in 1:nAngles
        v = sweepResult.sawVelMap2D[iw, ia]
        n = length(v)
        n == 0 && continue
        append!(allAng, fill(sweepResult.angles[ia], n))
        append!(allWav, fill(sweepResult.wavelengths[iw], n))
        append!(allRank, 1:n)
    end

    axOut = _maps_newax(ax; xlabel = "angle [deg]", ylabel = "wavelength",
                         title = "colour = local branch rank (no branch tracking yet)")
    scatter!(axOut, allAng, allWav; color = allRank, colormap = :viridis, markersize = 8)
    _maps_colorbar_beside(axOut; colormap = :viridis, colorrange = _maps_finite_extrema(allRank))
    return axOut
end

"""
    GuliwerrJ.plotDispersionMap(sweepResult::CalcFreqAngSweep2Layered, omegaIdx = 1,
                                 scaleMode = :log, climRange = (-5.0, -2.0),
                                 plotType = :cartesian; ax = nothing) -> Axis

isoFreq 2D case (angle x omega0, each omega0 its own angle x kParRange Green `|G|` heatmap
slice): plots ONE selected `omegaIdx` slice.
"""
function GuliwerrJ.plotDispersionMap(sweepResult::CalcFreqAngSweep2Layered,
                                      omegaIdx::Integer = 1, scaleMode::Symbol = :log,
                                      climRange = (-5.0, -2.0), plotType::Symbol = :cartesian;
                                      ax = nothing)
    (1 <= omegaIdx <= length(sweepResult.omega0List)) ||
        throw(ArgumentError("plotDispersionMap: omegaIdx must be an integer in " *
            "1:$(length(sweepResult.omega0List)), got $omegaIdx"))
    slice = sweepResult.sweeps[omegaIdx]
    axOut = _maps_isofreq_heatmap(slice, scaleMode, climRange, plotType, ax)
    axOut.title = "omega0 = $(round(sweepResult.omega0List[omegaIdx], sigdigits = 4))"
    return axOut
end

"""
    GuliwerrJ.plotDispersionMap(sweepResult::CalcIsofreqAngSweepLayered, omegaIdx = 1,
                                 scaleMode = :log, climRange = (-5.0, -2.0),
                                 plotType = :cartesian; ax = nothing) -> Axis

isoFreq, one omega0 (angle x kParRange): Green `|G|` heatmap, the same `ColMapGrayOrange`
treatment [`plotGreenMap`](@ref) uses. `omegaIdx` is ignored (no outer omega0 loop here).
"""
function GuliwerrJ.plotDispersionMap(sweepResult::CalcIsofreqAngSweepLayered,
                                      omegaIdx::Integer = 1, scaleMode::Symbol = :log,
                                      climRange = (-5.0, -2.0), plotType::Symbol = :cartesian;
                                      ax = nothing)
    return _maps_isofreq_heatmap(sweepResult, scaleMode, climRange, plotType, ax)
end

"""
    GuliwerrJ.plotDispersionMap(sweepResult::CalcIsofreqSingleLayered, omegaIdx = 1,
                                 scaleMode = :log, climRange = nothing,
                                 plotType = :cartesian; ax = nothing) -> Axis

isoFreq, one angle, one omega0 (kParRange only): 1D line, `|G_{33}|` (scaled) vs. kParRange.
No angle axis at all, so `omegaIdx`/`plotType` are ignored: always a Cartesian line.
`climRange` is the y-range for `:log`, default `(-5, 0)` (the normalised peak sits at 0);
the x-range always spans the whole kParRange.
"""
function GuliwerrJ.plotDispersionMap(sweepResult::CalcIsofreqSingleLayered,
                                      omegaIdx::Integer = 1, scaleMode::Symbol = :log,
                                      climRange = nothing, plotType::Symbol = :cartesian;
                                      ax = nothing)
    visValidateChoice(scaleMode, (:log, :linear), "plotDispersionMap", "scaleMode")
    plotData, yLabelStr = scaleGreenData(sweepResult.absG, scaleMode, "G_{33}")
    axOut = _maps_newax(ax; xlabel = "kParRange", ylabel = yLabelStr)
    lines!(axOut, sweepResult.kParRange, plotData)
    if scaleMode === :log
        yr = climRange === nothing ? (-5.0, 0.0) : climRange
        ylims!(axOut, Float64(yr[1]), Float64(yr[2]))
    end
    xlims!(axOut, extrema(sweepResult.kParRange)...)
    return axOut
end

# -----------------------------------------------------------------------------------------
# plotAngularMap -- port of guliwerr/visual/@Vis/plotAngularMap.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotAngularMap(sweepResult, plotType = :cartesian; ax = nothing) -> Axis

SAW velocity vs. angle -- Cartesian (velocity) or polar (SLOWNESS = 1/velocity, the standard
slowness-surface representation -- a real behavioural difference from `:cartesian`, not just
a coordinate-system swap, matching MATLAB). Every `sawVelMap[ia]` entry may be empty or hold
several velocities; flattened into one flat series across every angle. Limiting-velocity
reference lines are deliberately NOT drawn here -- `plotClassificationDiagnostic` is where
they belong (matches MATLAB's own docstring).

`plotType = :polar` returns a `PolarAxis` (`:cartesian` returns a plain `Axis`) -- mirroring
MATLAB's own `axes()`/`polaraxes()` split for this function. Unlike
[`plotPolarHeatmap`](@ref)'s pcolor-style maps, this is a plain scatter of discrete points --
exactly the primitive Makie's native `PolarAxis` renders correctly (verified separately from
the heatmap-on-`PolarAxis` finding documented at the top of this file), so no manual
pol2cart/quad-mesh workaround is needed here.
"""
function GuliwerrJ.plotAngularMap(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered},
                                   plotType::Symbol = :cartesian; ax = nothing)
    visValidateChoice(plotType, (:cartesian, :polar), "plotAngularMap", "plotType")

    angles = sweepResult.angles
    allAngles = Float64[]
    allVels = Float64[]
    for ia in eachindex(angles)
        v = sweepResult.sawVelMap[ia]
        append!(allAngles, fill(angles[ia], length(v)))
        append!(allVels, v)
    end

    if plotType === :cartesian
        axOut = _maps_newax(ax; xlabel = "angle [deg]", ylabel = "SAW velocity [km/s]")
        scatter!(axOut, allAngles, allVels; marker = :xcross, color = :red, markersize = 10)
        return axOut
    else
        axOut = ax === nothing ? PolarAxis(Figure()[1, 1]) : ax
        scatter!(axOut, deg2rad.(allAngles), 1 ./ allVels; marker = :xcross, color = :red,
                 markersize = 10)
        return axOut
    end
end

# -----------------------------------------------------------------------------------------
# plotRitzModeMap -- port of guliwerr/visual/@Vis/plotRitzModeMap.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotRitzModeMap(sweepResult, criterion = :ritzScore,
                               velCutoffFactor = 1.10; ax = nothing) -> Axis

Every computed Ritz mode as an angle x velocity scatter, coloured by a Ritz-mode criterion
([`visRitzCriterion`](@ref)): `:ritzScore` (default, dB; modes supporting no candidate as
small dark-grey dots) | `:freqDisp` (half-space: red ramp below the veto
`RitzScore.MIN_FREQDISP`, grey above, colour range [0, 1]) | `:depthInsensitivity`. Port of
`Vis.plotRitzModeMap`. `velCutoffFactor` -- upper velocity cutoff as a
multiple of the sweep's max longitudinal limiting velocity (`nothing`/`Inf` disables it).
Sorted ascending so the highest values draw on top; black background.
"""
function GuliwerrJ.plotRitzModeMap(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered},
                                    criterion::Symbol = :ritzScore,
                                    velCutoffFactor::Union{Real,Nothing} = 1.10;
                                    ax = nothing)
    visValidateChoice(criterion, visRitzCriterionNames(), "plotRitzModeMap", "criterion")

    velCutoff = if velCutoffFactor === nothing || !isfinite(velCutoffFactor)
        Inf
    else
        longVel = visPropOrEmpty(sweepResult, (:longVelMap, :subLongVelMap))
        if longVel === nothing || isempty(longVel) || all(isnan, longVel)
            longVel = [something(visPropOrEmpty(c, (:longLimVel, :subLongLimVel)), NaN)
                       for c in sweepResult.cs]
        end
        (isempty(longVel) || all(isnan, longVel)) &&
            throw(ArgumentError("Vis:plotRitzModeMap:noLongVel -- longitudinal limiting velocity " *
                "is unavailable (calcMode was not :full) -- cannot apply velCutoffFactor; " *
                "pass nothing or Inf to disable it."))
        velCutoffFactor * maximum(filter(!isnan, longVel))
    end

    angles = sweepResult.angles
    allAng = Float64[]; allVel = Float64[]; allCrit = Float64[]
    info = (label = String(criterion), vetoLevel = NaN)
    for iAngle in eachindex(angles)
        ritz = sweepResult.cs[iAngle].ritzSolver
        vels = ritz.vels[1:ritz.nCompModes]
        critVal, info = visRitzCriterion(sweepResult.cs[iAngle], criterion)
        keep = vels .<= velCutoff
        append!(allAng, fill(Float64(angles[iAngle]), count(keep)))
        append!(allVel, vels[keep])
        append!(allCrit, critVal[keep])
    end

    axOut = _maps_newax(ax; xlabel = "angle [deg]", ylabel = "velocity [km/s]",
                         title = "Ritz mode map -- $(info.label)")
    axOut.backgroundcolor = :black
    faint = isnan.(allCrit)
    if criterion === :ritzScore && any(faint)
        scatter!(axOut, allAng[faint], allVel[faint]; color = Makie.RGBf(0.35, 0.35, 0.35),
                 markersize = 3)
    end
    keepC = .!faint
    order = sortperm(allCrit[keepC])   # highest plotted last (on top)
    cAng = allAng[keepC][order]; cVel = allVel[keepC][order]; cVal = allCrit[keepC][order]
    cmap, critRange = _maps_ritz_criterion_colour(info, cVal)
    scatter!(axOut, cAng, cVel; color = cVal, colormap = cmap, colorrange = critRange,
             markersize = 6)
    _maps_colorbar_beside(axOut; colormap = cmap, colorrange = critRange, label = info.label)
    return axOut
end

"""Colormap and range of a Ritz-mode criterion -- port of private/visRitzCriterionColour.m:
plain grey, auto range; with a veto level, a red ramp below the veto then grey, range [0, 1]."""
function _maps_ritz_criterion_colour(info, vals)
    N = 256
    isnan(info.vetoLevel) && return (Makie.to_colormap(:grays), _maps_finite_extrema(vals))
    nVeto = round(Int, N * info.vetoLevel)
    ramp = range(0.25, 1, length = nVeto)
    vetoRows = [Makie.RGBf(r, 0.25r, 0.25r) for r in ramp]
    keepRows = [Makie.RGBf(g, g, g) for g in range(0, 1, length = N - nVeto)]
    return (vcat(vetoRows, keepRows), (0.0, 1.0))
end

# -----------------------------------------------------------------------------------------
# plotRitzModeMapPlate -- port of guliwerr/visual/@Vis/plotRitzModeMapPlate.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotRitzModeMapPlate(sweepResult::CalcAngSweepPlate, criterion = :modeType,
                                    velCutoffFactor = 1.10; ax = nothing) -> Axis

Plate analogue of [`plotRitzModeMap`](@ref): every computed Ritz mode as an angle x velocity
scatter, shaded by mode-TYPE classification (`RitzSolverLamb.modesType`/
`modeTypeConfidence`) rather than SAW's `modeRecog` criteria -- and there is no Green's
function map option at all. `criterion = :modeType` (default) colours by the 4 discrete
Lamb/plate mode types ([`plateModeTypeStyle`](@ref)'s names/colours, one labelled `scatter!`
series per type plus a legend below the axis, a discrete classification having no natural colorbar);
`criterion = :confidence` uses a continuous `:grays` colormap + colorbar in `[0,1]`, as
[`plotRitzModeMap`](@ref) does. `velCutoffFactor` -- upper velocity cutoff as a multiple of
the MAXIMUM over every angle of the structure's bottom-layer (and top-layer, if different)
Christoffel longitudinal limiting velocity, computed fresh via [`visPlateLongLimVel`](@ref)
since `CalcAngSweepPlate` stores no `longVelMap` of its own (pass `nothing`/`Inf` to disable).
"""
function GuliwerrJ.plotRitzModeMapPlate(sweepResult::CalcAngSweepPlate,
                                         criterion::Symbol = :modeType,
                                         velCutoffFactor::Union{Real,Nothing} = 1.10;
                                         ax = nothing)
    visValidateChoice(criterion, (:modeType, :confidence), "plotRitzModeMapPlate", "criterion")
    angles = sweepResult.angles

    velCutoff = if velCutoffFactor === nothing || !isfinite(velCutoffFactor)
        Inf
    else
        longLimByAngle = [visPlateLongLimVel(sweepResult.structure, a) for a in angles]
        velCutoffFactor * maximum(longLimByAngle)
    end

    allAng = Float64[]
    allVel = Float64[]
    allType = Int[]
    allConf = Float64[]
    for iAngle in eachindex(angles)
        vels, types, conf = visPlateRitzModeData(sweepResult.cs[iAngle].ritzSolver)
        keep = vels .<= velCutoff
        append!(allAng, fill(angles[iAngle], count(keep)))
        append!(allVel, vels[keep])
        append!(allType, types[keep])
        append!(allConf, conf[keep])
    end

    axOut = _maps_newax(ax; xlabel = "angle [deg]", ylabel = "velocity [km/s]",
                         title = "Ritz mode map (Plate) -- $criterion")

    if criterion === :modeType
        typeNames, typeColors = plateModeTypeStyle()
        for it in 1:4
            mask = allType .== it
            any(mask) || continue
            col = RGBf(typeColors[it, 1], typeColors[it, 2], typeColors[it, 3])
            scatter!(axOut, allAng[mask], allVel[mask]; color = col, label = typeNames[it],
                     markersize = 6)
        end
        axOut.backgroundcolor = :white
        _maps_legend_below(axOut)   # the branches reach every corner of the map
    else
        order = sortperm(allConf)   # highest confidence plotted last (on top)
        axOut.backgroundcolor = :black
        scatter!(axOut, allAng[order], allVel[order]; color = allConf[order], colormap = :grays,
                 colorrange = (0.0, 1.0), markersize = 6)
        _maps_colorbar_beside(axOut; colormap = :grays, colorrange = (0.0, 1.0),
                               label = "mode-type confidence")
    end
    return axOut
end

# -----------------------------------------------------------------------------------------
# plotRitzSpectrumMap -- port of guliwerr/visual/@Vis/plotRitzSpectrumMap.m
# -----------------------------------------------------------------------------------------

"""
    GuliwerrJ.plotRitzSpectrumMap(sweepResult, scaleMode = :log, climRange = (-4.0, -1.0),
                                   plotType = :cartesian; ax = nothing, showCands = false)
        -> (Axis, Colorbar)

Angle x velocity map of the Ritz spectrum `|S(v)|`, rebuilt from the stored modal numerators
(no solve; the default background of a `:ritzOnly` result). Same scale conventions as
[`plotGreenMap`](@ref) (clim `[-4, -1]`); `plotType` only `:cartesian` (the slot keeps the
argument order of `plotGreenMap`). `showCands` marks every angle's Ritz candidates at their
supporting modes' velocities (red dots). Port of `Vis.plotRitzSpectrumMap`; numeric half
[`ritzSpectrumMapData`](@ref).
"""
function GuliwerrJ.plotRitzSpectrumMap(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered},
                                        scaleMode::Symbol = :log, climRange = (-4.0, -1.0),
                                        plotType::Symbol = :cartesian; ax = nothing,
                                        showCands::Bool = false)
    visValidateChoice(plotType, (:cartesian,), "plotRitzSpectrumMap", "plotType")
    angles, velAxis, plotData, cbLabel, _ = ritzSpectrumMapData(sweepResult, scaleMode)
    colorrange = scaleMode === :log ? (Float64(climRange[1]), Float64(climRange[2])) :
                 _maps_finite_extrema(plotData)
    axOut = _maps_newax(ax; xlabel = "angle [deg]", ylabel = "velocity [km/s]",
                         title = "Ritz spectrum map")
    heatmap!(axOut, angles, velAxis, permutedims(plotData); colormap = _MAPS_GREEN_CMAP,
             colorrange = colorrange)
    if showCands
        for (ia, c) in enumerate(sweepResult.cs)
            isempty(c.ritzScoreCands) && continue
            mv = [p.modeVel for p in c.ritzScoreCands]
            scatter!(axOut, fill(Float64(angles[ia]), length(mv)), mv; color = :red,
                     markersize = 5)
        end
    end
    cb = _maps_colorbar_beside(axOut; colormap = _MAPS_GREEN_CMAP, colorrange = colorrange,
                                label = cbLabel)
    return axOut, cb
end
