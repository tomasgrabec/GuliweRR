# =========================================================================================
# PlotsModalShape.jl -- Makie extension, "modal shape / displacement field" group
#
# Implements the 3 stubs declared in `src/visual/PlotStubs.jl` for the signed 2D
# (depth x lateral) displacement-field family: `plotModalShape` (Ritz eigenmode),
# `plotGreenModalShape` (Green's response, half-space), `plotGreenModalShapeLayered`
# (Green's response, layered). Ports of `guliwerr/visual/@Vis/plotModalShape.m`,
# `plotGreenModalShape.m`, `plotGreenModalShapeLayered.m` and their shared private helper
# `visModalShapePanel.m`. The numeric halves (`calcDispField2D`, `calcGreenDispField2D`,
# `calcGreenDispField2DLayered`, `calcSubstrateCutoffDepth`, `calcModalShapeRowHeights`,
# `visDepthMagnitude`, `visDecayCutoffIdx`) live in `src/visual/VisData.jl` and
# are called here unqualified -- `GuliwerrJ` re-exports all of them (see `src/GuliwerrJ.jl`),
# and this extension does `using GuliwerrJ` at its top (`GuliwerrJMakieExt.jl`).
#
# METHODS ADDED TO THE STUBS, fully qualified per the extension's own house rule (never a
# bare `function plotXxx(...)`, which would create an unrelated second binding instead of
# extending `GuliwerrJ.plotXxx`).
#
# ---------------------------------------------------------------------------------------
# GRID LAYOUT -- ONE Axis PER DIRECTION, not one per (direction, layer) pair
# ---------------------------------------------------------------------------------------
# MATLAB builds an `{numel(dirs) x nLay}` grid of SEPARATE `axes()` objects, each manually
# `Position`-ed as a figure-fraction rectangle, to stack `nLay` layers vertically per
# direction column (layer 1 = substrate at the bottom, layer nLay = outermost film at the
# top -- see plotModalShape.m's own docblock). The stub this file implements is documented
# to return `Vector{Axis}` of length `numel(dirs)` (PlotStubs.jl), not an `nDirs x nLay`
# grid -- so the nLay-layer stack for one direction is consolidated into a SINGLE Makie
# `Axis` per direction, with each layer drawn as its own `heatmap!` call occupying its own
# vertical band `[layBottom[ilay], layBottom[ilay] + layHeightFrac[ilay]]` of that one
# axis's `y in [0,1]` range -- the exact same row-height fractions
# `calcModalShapeRowHeights`/proportional-to-`relDepth` allocation MATLAB computes for its
# per-panel `Position` heights, just realised as multiple heatmaps sharing one Axis's data
# coordinates instead of multiple sibling axes objects. A thin black rectangle is drawn
# around each layer's own band (`lines!`) to reproduce MATLAB's per-panel frame
# (`rectangle(...,'EdgeColor','k')`) and to make the layer boundaries visible when they
# are stacked inside one axis. This is a DELIBERATE STRUCTURAL DEVIATION from the MATLAB
# source, dictated by the stub's own documented return type -- not an oversight.
#
# ---------------------------------------------------------------------------------------
# COLOUR RANGE -- data-driven and symmetric about zero, NOT a hardcoded clim([-1,1])
# ---------------------------------------------------------------------------------------
# `visDispField2D` (VisData.jl) normalises by the GLOBAL MAX OVER POSITIVE VALUES ONLY, so
# a negative-dominant mode's most negative value can fall below -1 even though the positive
# side never exceeds +1. MATLAB's `visModalShapePanel.m` then hardcodes `clim(ax,[-1,1])`,
# which its own comment admits "clips a negative-dominant mode by design". This port does
# NOT reproduce that clip: `_modalShapeColLimit` reads the ACTUAL max-abs value across the
# full (uncropped, all-3-direction, all-layer) field -- the same population
# `visDispField2D` normalised over -- and uses `(-climLimit, climLimit)` as the Makie
# `colorrange` for every heatmap in the figure. The range is kept SYMMETRIC ABOUT ZERO
# (not simply `[min,max]`) deliberately: `ColormapRWB`'s centre row is exactly white at the
# colormap's midpoint (Colormaps.jl), and only a zero-centred range keeps that white
# midpoint aligned with physical zero displacement -- an asymmetric `[dataMin,dataMax]`
# range would shift white away from zero and misrepresent the field's sign structure, which
# matters more here than matching MATLAB's exact (clipping) numeric convention. A
# `Colorbar` is added next to each figure (not part of the stub's returned `Vector{Axis}`)
# specifically so this non-[-1,1] range is visible to the caller, not silently implied.
#
# ---------------------------------------------------------------------------------------
# PHASE CONVENTION
# ---------------------------------------------------------------------------------------
# `calcDispField2D`'s underlying `calcDisplacement` output has its phase fixed by
# `phaseFix` (largest in-plane Legendre coefficient real and positive), as in MATLAB's
# `RitzSolver.phaseFix`. The SIGN of every field plotted here is therefore well-defined and
# matches MATLAB's (up to ties in magnitude); nothing in this file needs to compensate for it.
# =========================================================================================

const _MODALSHAPE_DIR_NAMES = ("u_x (propag.dir.)", "u_y (in-plane)", "u_z (out-of-plane)")

"""
    _modalShapeValidateDirs(dirs)

Every requested direction must be a valid displacement-component index (1-3). Not validated
by the MATLAB source (which trusts its caller); added here since an out-of-range `dirs`
entry would otherwise surface as an opaque `BoundsError` deep inside plotting code.
"""
function _modalShapeValidateDirs(dirs)
    all(d -> 1 <= d <= 3, dirs) ||
        throw(ArgumentError("plotModalShape family: dirs entries must be integers in 1:3, " *
            "got $dirs"))
    return dirs
end

"""
    _modalShapeColLimit(dispField) -> Float64

Symmetric colour-range half-width: the actual max-abs value across every direction/layer of
`dispField` (a `(3, nLay)` matrix of real 2D fields, `visDispField2D`'s own output
convention), floored at `eps()`. See this file's header, "COLOUR RANGE".
"""
function _modalShapeColLimit(dispField::AbstractMatrix{<:AbstractMatrix{<:Real}})
    m = 0.0
    for field2D in dispField
        isempty(field2D) && continue
        m = max(m, maximum(abs, field2D))
    end
    return max(m, eps())
end

"""
    _modalShapeColormap() -> Vector{<:Makie.Colorant}

`ColormapRWB`'s `N x 3` RGB matrix (Colormaps.jl), converted to the row-of-colours form
Makie's `colormap`/`colorrange` keywords expect.
"""
function _modalShapeColormap()
    cm = ColormapRWB()
    return [Makie.RGBf(cm[i, 1], cm[i, 2], cm[i, 3]) for i in axes(cm, 1)]
end

"""
    _modalShapeCropSubstrate!(dispField, keepFrac) -> dispField

Crop layer 1 (the substrate)'s field to its shallow end (the LAST `keepN` rows -- row order
is deepest-first, see `calcDispField2D`'s own docblock) for every one of the 3 directions,
in place. Port of the identical block shared by `plotModalShape.m`/
`plotGreenModalShapeLayered.m`.
"""
function _modalShapeCropSubstrate!(dispField::AbstractMatrix{<:AbstractMatrix{<:Real}},
                                    keepFrac::Real)
    for iDir in 1:3
        nRows = size(dispField[iDir, 1], 1)
        keepN = max(1, round(Int, keepFrac * nRows))
        dispField[iDir, 1] = dispField[iDir, 1][(end - keepN + 1):end, :]
    end
    return dispField
end

"""
    _modalShapeRowHeightFracs(substrateCutoff, dispRelDepths) -> Vector{Float64}

Per-layer axis-height fraction (sums to 1): `calcModalShapeRowHeights` for `:boosted`
(non-proportional, every layer gets a guaranteed floor share -- VisData.jl), otherwise
plain-proportional to `dispRelDepths`. Port of the identical `layHeightFrac` block shared by
`plotModalShape.m`/`plotGreenModalShapeLayered.m` (MATLAB's `topMargin`/`bottomMargin`
figure-chrome reservation is dropped -- see this file's header; the PROPORTIONS, which are
the physically meaningful part, are preserved exactly).
"""
function _modalShapeRowHeightFracs(substrateCutoff::Symbol, dispRelDepths::AbstractVector{<:Real})
    if substrateCutoff === :boosted
        frac = calcModalShapeRowHeights(dispRelDepths)
        return frac isa Real ? [frac] : frac
    else
        return dispRelDepths ./ sum(dispRelDepths)
    end
end

"""
    _modalShapeLayerBands(layHeightFrac) -> Vector{Float64}

Bottom-edge `y` coordinate of every layer's band, bottom (substrate, layer 1) to top
(outermost film), given each layer's height fraction. `layBottom[ilay] +
layHeightFrac[ilay]` is that layer's top edge; the full stack spans `y in [0,1]`.
"""
function _modalShapeLayerBands(layHeightFrac::AbstractVector{<:Real})
    bottoms = Vector{Float64}(undef, length(layHeightFrac))
    acc = 0.0
    for ilay in eachindex(layHeightFrac)
        bottoms[ilay] = acc
        acc += layHeightFrac[ilay]
    end
    return bottoms
end

"""
    _modalShapeDrawColumn!(ax, bands, cMap, climLimit) -> ax

Render one direction's full layer stack into a single `Axis`. `bands` is a vector of
`(field2D, bandLo, bandHi)` triples, bottom (substrate) to top (outermost film); each is
drawn as its own `heatmap!` occupying `x in [0,1]` (lateral / wavelength) x `y in
[bandLo,bandHi]`, with a thin black frame around it (`visModalShapePanel.m`'s
`rectangle(...)`, one per layer here since layers share one axis). Axis decorations/spines
are hidden to match MATLAB's `axis(ax,'off')`.
"""
function _modalShapeDrawColumn!(ax::Axis, bands, cMap, climLimit::Real)
    for (field2D, bandLo, bandHi) in bands
        nMeshDepth, nMeshLateral = size(field2D)
        xNorm = range(0.0, 1.0; length = nMeshLateral)     # x / wavelength
        yNorm = range(bandLo, bandHi; length = nMeshDepth)  # this layer's own band
        heatmap!(ax, xNorm, yNorm, permutedims(field2D);
                 colormap = cMap, colorrange = (-climLimit, climLimit))
        lines!(ax, [0.0, 1.0, 1.0, 0.0, 0.0], [bandLo, bandLo, bandHi, bandHi, bandLo];
               color = :black, linewidth = 1)
    end
    hidedecorations!(ax)
    hidespines!(ax)
    xlims!(ax, 0.0, 1.0)
    ylims!(ax, 0.0, 1.0)
    return ax
end

"""
    _modalShapeFigureAndAxes(fig, nDirs) -> (Figure, Vector{Axis})

Fresh `Figure()` if `fig === nothing`, else `fig` itself cleared via `empty!` -- the Makie
analogue of MATLAB's own `hFig = figure()` / `clf(hFig)` pair. Builds one `GridLayout` with
`nDirs` `Axis` columns plus one trailing `Colorbar` column (see this file's header, "COLOUR
RANGE") and returns `(figure, axesVector)`; `axesVector` (length `nDirs`, one per requested
direction, in `dirs` order) is what every public function here returns, matching each
stub's documented `-> Vector{Axis}`.
"""
function _modalShapeFigureAndAxes(fig::Union{Figure,Nothing}, nDirs::Integer, climLimit::Real, cMap)
    f = fig === nothing ? Figure() : fig
    fig === nothing || empty!(f)
    gl = f[1, 1] = GridLayout()
    axes_ = Vector{Axis}(undef, nDirs)
    for ic in 1:nDirs
        axes_[ic] = Axis(gl[1, ic])
    end
    Colorbar(gl[1, nDirs + 1]; colorrange = (-climLimit, climLimit), colormap = cMap,
             label = "displacement (normalised)")
    return f, axes_
end

"""
    _modalShapeLayeredPlot(dispField, relDepths, substrateCutoff, cutoffParam, dirs, fig,
                            errPrefix) -> Vector{Axis}

Shared substrate-cutoff / row-height / draw logic behind [`GuliwerrJ.plotModalShape`](@ref)
and [`GuliwerrJ.plotGreenModalShapeLayered`](@ref) -- both mirror `plotModalShape.m`'s own
structure (per that MATLAB file's docblock: "The axis-stacking problem here is the SAME one
`plotModalShape` already solves ... this function's structure mirrors `plotModalShape.m`
directly"). `dispField` is the ALREADY-COMPUTED `(3, nLay)` field (mutated in place by the
substrate crop); `relDepths` is Structure order, layer 1 = substrate.
"""
function _modalShapeLayeredPlot(dispField::AbstractMatrix{<:AbstractMatrix{<:Real}},
                                 relDepths::AbstractVector{<:Real},
                                 substrateCutoff::Symbol,
                                 cutoffParam::Union{Real,Nothing},
                                 dirs,
                                 fig::Union{Figure,Nothing},
                                 errPrefix::AbstractString)
    isRgb = visIsRgbDirs(dirs, errPrefix)
    isRgb || _modalShapeValidateDirs(dirs)
    nLay = length(relDepths)
    climLimit = _modalShapeColLimit(dispField)

    dispRelDepths = copy(relDepths)
    if substrateCutoff === :decay
        substrateDispMag = visDepthMagnitude([dispField[iDir, 1] for iDir in 1:3])
        dispRelDepths[1] = calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths,
                                                     substrateDispMag)
    elseif substrateCutoff === :fixedMultiple
        dispRelDepths[1] = calcSubstrateCutoffDepth(substrateCutoff, cutoffParam, relDepths)
    elseif substrateCutoff !== :none && substrateCutoff !== :boosted
        throw(ArgumentError("$errPrefix:badSubstrateCutoff -- substrateCutoff must be :none, " *
            ":decay, :fixedMultiple, or :boosted (got '$substrateCutoff')."))
    end

    if dispRelDepths[1] < relDepths[1]
        keepFrac = dispRelDepths[1] / relDepths[1]
        _modalShapeCropSubstrate!(dispField, keepFrac)
    end

    layHeightFrac = _modalShapeRowHeightFracs(substrateCutoff, dispRelDepths)
    layBottom = _modalShapeLayerBands(layHeightFrac)
    isRgb && return _modalShapeRgbPlot(dispField, layBottom, layHeightFrac, fig)

    nDirs = length(dirs)
    cMap = _modalShapeColormap()
    _, axes_ = _modalShapeFigureAndAxes(fig, nDirs, climLimit, cMap)

    for (ic, iDir) in enumerate(dirs)
        bands = [(dispField[iDir, ilay], layBottom[ilay], layBottom[ilay] + layHeightFrac[ilay])
                 for ilay in 1:nLay]
        _modalShapeDrawColumn!(axes_[ic], bands, cMap, climLimit)
        axes_[ic].title = _MODALSHAPE_DIR_NAMES[iDir]
    end

    return axes_
end

# -----------------------------------------------------------------------------------------
# Polarisation-RGB view (dirs = :rgb) -- private/visRgbShapeStack.m + plotPolarisationKey.m
# -----------------------------------------------------------------------------------------

"""
    _modalShapeRgbPlot(dispField, layBottom, layHeightFrac, fig) -> Vector{Axis}

The `dirs = :rgb` branch of the modal-shape plots: one panel with every layer's
[`calcPolarisationRgb`](@ref) image in its band (thin grey amplitude contours at 0.2:0.2:0.8,
black frame), plus the colour key ([`plotPolarisationKey`](@ref)) beside it. Returns
`[panelAxis, keyAxis]`.
"""
function _modalShapeRgbPlot(dispField, layBottom, layHeightFrac, fig)
    rgbCell, ampCell = calcPolarisationRgb(dispField)
    f = fig === nothing ? Figure() : fig
    fig === nothing || empty!(f)
    gl = f[1, 1] = GridLayout()
    ax = Axis(gl[1, 1]; title = "u (polarisation RGB)")
    for ilay in eachindex(rgbCell)
        lo = layBottom[ilay]; hi = lo + layHeightFrac[ilay]
        rgb = rgbCell[ilay]; amp = ampCell[ilay]
        nR, nC = size(amp)
        img = [Makie.RGBf(rgb[i, j, 1], rgb[i, j, 2], rgb[i, j, 3]) for j in 1:nC, i in 1:nR]
        image!(ax, (0.0, 1.0), (lo, hi), img)
        if maximum(amp) - minimum(amp) > eps()
            contour!(ax, range(0.0, 1.0; length = nC), range(lo, hi; length = nR),
                     permutedims(amp); levels = 0.2:0.2:0.8, color = Makie.RGBf(0.5, 0.5, 0.5),
                     linewidth = 0.3)
        end
        lines!(ax, [0.0, 1.0, 1.0, 0.0, 0.0], [lo, lo, hi, hi, lo]; color = :black, linewidth = 1)
    end
    hidedecorations!(ax); hidespines!(ax)
    xlims!(ax, 0.0, 1.0); ylims!(ax, 0.0, 1.0)
    keyAx = Axis(gl[1, 2])
    GuliwerrJ.plotPolarisationKey(; ax = keyAx)
    colsize!(gl, 2, Relative(0.4))
    return [ax, keyAx]
end

"""
    GuliwerrJ.plotPolarisationKey(; ax = nothing, nSub = 40) -> Axis

Colour key of the polarisation-RGB view: the positive octant of the unit sphere of
directions `(|u1|, |u2|, |u3|)/|u|`, triangulated (`nSub` subdivisions), coloured exactly as
[`calcPolarisationRgb`](@ref) colours a saturated pixel, viewed down the (1,1,1) axis, with
the three boundary arcs and corner labels. Port of `Vis.plotPolarisationKey`.
"""
function GuliwerrJ.plotPolarisationKey(; ax = nothing, nSub::Integer = 40)
    nSub >= 2 || throw(ArgumentError("plotPolarisationKey: nSub must be an integer >= 2."))
    if ax === nothing
        ax = Axis(Figure()[1, 1])
    else
        empty!(ax)
    end
    idx = zeros(Int, nSub + 1, nSub + 1)
    P = Vector{NTuple{3,Float64}}()
    for j in 0:nSub, i in 0:nSub     # MATLAB meshgrid column order
        i + j <= nSub || continue
        k = nSub - i - j
        v = (i, j, k) ./ nSub
        push!(P, v ./ sqrt(sum(abs2, v)))
        idx[i + 1, j + 1] = length(P)
    end
    faces = Vector{NTuple{3,Int}}()
    for i in 0:(nSub - 1), j in 0:(nSub - 1 - i)
        push!(faces, (idx[i + 1, j + 1], idx[i + 2, j + 1], idx[i + 1, j + 2]))
        i + j < nSub - 1 && push!(faces, (idx[i + 2, j + 1], idx[i + 2, j + 2], idx[i + 1, j + 2]))
    end
    eR = (1.0, -1.0, 0.0) ./ sqrt(2); eU = (-1.0, -1.0, 2.0) ./ sqrt(6)
    toScreen(q) = Point2f(sum(q .* eR), sum(q .* eU))
    base = visPolarisationColours()
    cols = [Makie.RGBf((p[1] * base[1, c] + p[2] * base[2, c] + p[3] * base[3, c] for c in 1:3)...)
            for p in P]
    mesh!(ax, toScreen.(P), [f[k] for f in faces, k in 1:3]; color = cols, shading = Makie.NoShading)
    t = range(0, pi / 2; length = 60)
    for arc in ([(cos(x), sin(x), 0.0) for x in t], [(0.0, cos(x), sin(x)) for x in t],
                [(sin(x), 0.0, cos(x)) for x in t])
        lines!(ax, toScreen.(arc); color = :black, linewidth = 0.75)
    end
    labels = ("u_1 (propag.)", "u_2 (in-plane)", "u_3 (out-of-plane)")
    offs = ((0.0, -0.12), (0.0, -0.12), (0.0, 0.06))
    for k in 1:3
        c = toScreen(ntuple(m -> m == k ? 1.0 : 0.0, 3))
        text!(ax, c[1] + offs[k][1], c[2] + offs[k][2]; text = labels[k], fontsize = 11,
              align = (:center, k == 3 ? :bottom : :top),
              color = Makie.RGBf((0.8 .* base[k, :])...))
    end
    hidedecorations!(ax); hidespines!(ax)
    ax.aspect = DataAspect()
    xlims!(ax, -1.45, 1.45); ylims!(ax, -0.85, 1.15)
    ax.title = "polarisation"
    return ax
end

# -----------------------------------------------------------------------------------------
# plotModalShape -- plotModalShape.m
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotModalShape(ritzSolver::RitzSolver, modeNo::Integer,
                                   substrateCutoff::Symbol = :none,
                                   cutoffParam::Union{Real,Nothing} = nothing,
                                   dirs = 1:3,
                                   meshSpacing::Real = 0.01;
                                   fig::Union{Figure,Nothing} = nothing)
    nColsAvail = size(ritzSolver.eigVecs_con, 2)
    if isempty(ritzSolver.eigVecs_con) || nColsAvail < modeNo
        throw(ArgumentError("plotModalShape:modeNotAvailable -- mode $modeNo is not " *
            "available on this RitzSolver (eigVecs_con has $nColsAvail column(s)) -- " *
            "retention-aware recompute-on-demand is not built yet; construct/calc this " *
            "RitzSolver with retentionLevel = :full."))
    end

    relDepths = [cp.relDepth for cp in ritzSolver.structure.calcParams]
    dispField = calcDispField2D(ritzSolver, modeNo, meshSpacing)   # (3, nLay)

    return _modalShapeLayeredPlot(dispField, relDepths, substrateCutoff, cutoffParam, dirs,
                                   fig, "plotModalShape")
end

# -----------------------------------------------------------------------------------------
# plotGreenModalShape -- plotGreenModalShape.m (half-space, always exactly one "layer")
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotGreenModalShape(greenObj::GreenSolver, forceDir::Integer,
                                        cutoffParam::Real = 0.02,
                                        dirs = 1:3,
                                        meshSpacing::Real = 0.01;
                                        fig::Union{Figure,Nothing} = nothing)
    (0 < cutoffParam < 1) ||
        throw(ArgumentError("plotGreenModalShape: cutoffParam must be in (0,1), got " *
            "$cutoffParam"))
    isRgb = visIsRgbDirs(dirs, "plotGreenModalShape")
    isRgb || _modalShapeValidateDirs(dirs)

    dispField = calcGreenDispField2D(greenObj, forceDir, meshSpacing)   # (3, 1)
    climLimit = _modalShapeColLimit(dispField)

    keepIdx = visDecayCutoffIdx(visDepthMagnitude([dispField[iDir, 1] for iDir in 1:3]),
                                 cutoffParam)
    for iDir in 1:3
        dispField[iDir, 1] = dispField[iDir, 1][keepIdx:end, :]
    end
    isRgb && return _modalShapeRgbPlot(dispField, [0.0], [1.0], fig)

    nDirs = length(dirs)
    cMap = _modalShapeColormap()
    _, axes_ = _modalShapeFigureAndAxes(fig, nDirs, climLimit, cMap)

    for (ic, iDir) in enumerate(dirs)
        _modalShapeDrawColumn!(axes_[ic], [(dispField[iDir, 1], 0.0, 1.0)], cMap, climLimit)
        axes_[ic].title = _MODALSHAPE_DIR_NAMES[iDir]
    end

    return axes_
end

# -----------------------------------------------------------------------------------------
# plotGreenModalShapeLayered -- plotGreenModalShapeLayered.m
# -----------------------------------------------------------------------------------------

function GuliwerrJ.plotGreenModalShapeLayered(greenObj::GreenSolver, structure::Structure,
                                               forceDir::Integer,
                                               substrateCutoff::Symbol = :decay,
                                               cutoffParam::Union{Real,Nothing} = nothing,
                                               dirs = 1:3,
                                               meshSpacing::Real = 0.01;
                                               fig::Union{Figure,Nothing} = nothing)
    relDepths = [cp.relDepth for cp in structure.calcParams]
    dispField = calcGreenDispField2DLayered(greenObj, structure, forceDir, meshSpacing)

    return _modalShapeLayeredPlot(dispField, relDepths, substrateCutoff, cutoffParam, dirs,
                                   fig, "plotGreenModalShapeLayered")
end
