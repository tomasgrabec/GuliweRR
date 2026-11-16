# =========================================================================================
# GuiDispersionViewer.jl -- viewer for 2D wavelength/wavevector x angle sweeps ("map + shape
# grid in one window, with branch tracking")
#
# Julia/Makie port of `guliwerr/gui/DispersionViewerGUI.m` (17 methods) and
# `DispersionViewerGUIPlate.m`. The pair is a near-duplicate in MATLAB with one real feature
# delta: identical 17-method surface; `wavelengths` -> `wavVecs`; the Ritz-branch/Green-
# background checkbox pair is replaced by a colour-mode dropdown + Christoffel-bounds checkbox
# for Plate; frequency read from `ritz.freqsPhysical`. The viewer is a genuine hybrid of the
# angle-sweep map viewer and the modal-shape viewer in one class, not reducible to
# composition: ONE window shows a dispersion-curve MAP (the `plotDispersionCurves`/`...Plate`
# machinery, click-to-select, colour-range control) stacked with a depth-resolved modal-SHAPE
# display (the `plotModalShape` machinery, stepper-driven), tied together by tracking which raw
# Ritz mode index continues each physical dispersion BRANCH across the wavelength/wavevector
# sweep.
#
# `include`d directly into `module GuliwerrJMakieExt` (see that file's own header), which
# already has `using GuliwerrJ`/`using Makie` in scope -- neither repeated here. Built on
# `GuiWidgets.jl` (`guiHostLayout`, `guiDividedTopBottom`, `guiPlotWithSidePanel`,
# `guiDropdownStepperRow`, `guiIndexStepperRow`, `guiRangeFieldPair`, `guiBindRangeToAxis!`,
# `guiOnAxisClick`, `guiNearestPoint`, `guiCrosshair!`, `guiStatusLabel`/`guiSetStatusOk!`/
# `guiSetStatusError!`) and on `GuiSingleAngleViewer.jl`'s `guiSignificantModeType`/
# `guiPlateModeTypeText` -- this file MUST be `include`d AFTER `GuiSingleAngleViewer.jl` AND
# after `PlotsDiagnostic.jl` (`plotDispersionCurves`/`plotDispersionCurvesPlate`), as
# `GuliwerrJMakieExt.jl` does.
#
# ---------------------------------------------------------------------------------------
# CONSTRUCTOR INPUTS
# ---------------------------------------------------------------------------------------
# `DispersionViewerGUI(sweepResult, windowMode, angleIdx, parentContainer)` requires
# `sweepResult` to be a `CalcWavevecAngSweep2Layered`
# (`validateattributes(sweepResult, {'CalcWavevecAngSweep2Layered'}, ...)`), already `calc()`'d
# with `retentionLevel='full'` (needed by `ModeClassification.trackDispersionBranches`, called
# directly in the constructor: `obj.branchModeNoMap = ModeClassification.trackDispersionBranches
# (sweepResult, angleIdx)`). `DispersionViewerGUIPlate` is the same shape for
# `CalcWavevecAngSweep2Plate`. Both are 2D (angle x wavelength-or-wavevector) sweeps
# (`src/forwardCalcs/CalcWavevecAngSweep2Layered.jl`/`CalcWavevecAngSweep2Plate.jl`).
# `angleIdx` (default 1) just slices which one of possibly
# several swept angles this fixed-angle view is pinned to -- it adds no angle stepper of its
# own (`GuiDoubleSweepViewer.jl` is what re-embeds a fresh viewer at a new `angleIdx` for
# angle-switching).
#
# ---------------------------------------------------------------------------------------
# "BRANCH TRACKING" -- what it means, and where it is computed
# ---------------------------------------------------------------------------------------
# `ModeClassification.trackDispersionBranches` (`guliwerr/forwardCalcs/ModeClassification.m`)
# tracks each Ritz mode's PHYSICAL IDENTITY across
# a fixed-angle wavelength/wavevector sweep by maximum eigenvector overlap (a Modal Assurance
# Criterion, MAC = |unit-normalised eigenvector inner product|) with the PREVIOUS step: seed
# every branch at its own raw column index at step 1 (`branchModeNoMap(1,:) = 1:nCompModes`),
# then at each later step compute one `[nCompModes x nCompModes]` MAC matrix between the current
# REFERENCE eigenvectors and that step's full eigenvector matrix, take each row's argmax as that
# branch's continuation, and CHAIN the reference forward to the matched column. This is needed
# because the RAW Ritz mode index at a given step is just a numerical ordering by eigenvalue,
# not a physical identity, and it can shift between steps as branches approach/cross -- the
# whole point of `branchModeNoMap[iw, ib]` is "which raw mode index continues tracked branch
# `ib`'s identity at wavelength/wavevector step `iw`".
#
# `trackDispersionBranches` is ported in `src/forwardCalcs/ModeClassification.jl` and is
# also consumed by `CalcWavevecAngSweep2Plate.jl`'s `classifyBranches`/
# `reindexByBranch`/`locateZGVs`/`ZGVPoint` and by `PlotsDiagnostic.jl`'s
# `plotDispersionCurves`/`plotDispersionCurvesPlate`
# (`dispersionCurvesData`, `PlotData.jl`, calls it directly). So "branch tracking" here is
# the SAME `trackDispersionBranches` call every one of those pieces makes, called ONCE at
# construction (`DispersionViewer.branchModeNoMap`, exactly mirroring MATLAB's own
# `obj.branchModeNoMap` property) and consumed two ways: (1) the "Branch:" stepper below steps
# through COLUMNS of it (a physical branch, not a raw mode number) --
# `branchModeNoMap[wavNo, branchIdx]` is the raw Ritz mode index the shape grid renders;
# (2) `PlotsDiagnostic.jl`'s `plotDispersionCurves`/`Plate` draw ONE LINE PER TRACKED BRANCH
# on the map (calling `trackDispersionBranches`/`classifyBranches` a SECOND time internally --
# a small, deliberate duplication, see "CLICK-TO-SELECT" below for why this file does not try
# to share that internal computation).
#
# ---------------------------------------------------------------------------------------
# THE SHAPE GRID OWNS ITS OWN AXES -- this viewer needs only ONE Figure
# ---------------------------------------------------------------------------------------
# `ext/GuiModalShapeViewer.jl`'s own header documents the obstacle: `plotModalShape`/
# `plotGreenModalShape`/`plotGreenModalShapeLayered` (`PlotsModalShape.jl`) each own their OWN
# `Figure`/`Axis` construction via `_modalShapeFigureAndAxes`, which calls `empty!(fig)` on
# WHATEVER `Figure` it is handed before rebuilding its own `GridLayout` at `fig[1,1]` -- so
# that `Figure` can never safely be shared with anything else a caller might also want to draw
# into it (a controls panel, or here, a dispersion-curve MAP). The modal-shape viewer's
# resolution is "TWO FIGURES, ALWAYS": a controls-host `fig` plus an ALWAYS-STANDALONE
# `shapeFig`, in effect making MATLAB's own 'split' windowMode the only case that survives.
#
# This viewer chooses DIFFERENTLY, and deliberately: both MATLAB sources' own
# `redrawShape` methods show the shape display here is ALWAYS the plain Ritz modal shape at
# `'boosted'` cutoff -- `Vis.plotModalShape(ritz, modeNo, 'boosted', [], 1:3, [], [],
# obj.shapeAxesCell)` -- literally, unconditionally, for BOTH the Layered and Plate variants.
# UNLIKE the modal-shape viewer's base variant, there is NO Source dropdown here, no
# Green's-response branch, no call that would ever need `plotGreenModalShape`'s/
# `plotGreenModalShapeLayered`'s own Figure-owning machinery (both class docblocks state: "No
# Green's-response toggle ... no SAW/limiting override panel"). Given that, the numeric
# building blocks this ONE display mode needs -- `calcDispField2D`, `calcModalShapeRowHeights`,
# `ColormapRWB` -- are all PUBLIC, Figure-agnostic functions (`src/visual/VisData.jl`/
# `Colormaps.jl`, exported from `GuliwerrJ`), not entangled with `PlotsModalShape.jl`'s own
# private Figure-clearing internals (`_modalShapeFigureAndAxes`, `_modalShapeDrawColumn!`,
# etc.). Following the convention explained in `GuiSingleAngleViewer.jl`'s deviation 5, a
# viewer file does not reach into another file's `_`-prefixed private helpers.
#
# So this file writes its OWN small, file-local drawing helper (`_dvBuildShapeAxes!`/
# `_dvRedrawShapeGrid!` below) that re-implements the SAME "one Axis per direction, multiple
# stacked `heatmap!` layer-bands, boosted row-height fractions, symmetric-about-zero colour
# range" technique `PlotsModalShape.jl`'s own `_modalShapeDrawColumn!`/`_modalShapeLayeredPlot`
# use for `plotModalShape` (a DELIBERATE DUPLICATION of ~30 lines of drawing
# logic, in line with the "each viewer file keeps its own private copy" convention:
# `_aswReapplyRange!`/`_savReapplyRange!`,
# `_msvWireMenuStepper!`/`_aswWireStepper!`/`_savWireStepper!`), but draws directly into `Axis`
# objects THIS file creates and owns inside a `GridLayout` region it shares with the map --
# never touching `plotModalShape`'s own Figure-clearing call at all. The 3 direction axes are
# built ONCE at construction and reused via `empty!(ax)` on every redraw (mirroring MATLAB's own
# pre-built `{3 x nLay}` `shapeAxesCell`, reused across redraws via
# `Vis.plotModalShape(...,obj.shapeAxesCell)` -- MATLAB's own "build once, redraw into" idiom);
# only the small colourbar cell is
# cleared+rebuilt each redraw (via the SAME `Makie.GridLayoutBase.remove_from_gridlayout!`
# technique `ext/GuiAngSweepViewer.jl`'s `guiClearGridLayout!` uses, which works for a
# `Colorbar` `Block` exactly as it does for an `Axis`/`Menu`), since
# each redraw's colour half-width (`climLimit`) genuinely changes per mode and rebuilding one
# small `Block` is simpler and safer than relying on live-mutable `Colorbar.colorrange`
# semantics.
#
# NET RESULT: this viewer needs only ONE Makie `Figure`, not two -- a genuine improvement over
# the modal-shape viewer's compromise, possible ONLY because this viewer's shape source has no
# Source-toggle branch that would need a Figure-owning plotting call. `windowMode`
# ('combined'/'split') therefore collapses entirely, as in the angle-sweep and single-angle
# viewers (no `fig`-splitting concept ported at all): there is no reason a caller would ever
# want two separate windows here. This is a deliberate divergence from MATLAB.
#
# ---------------------------------------------------------------------------------------
# CLICK-TO-SELECT -- a JOINT 2D search across every branch line, `guiNearestPoint`'s own
# purpose-built case (unlike the angle-sweep viewer's two INDEPENDENT 1D searches)
# ---------------------------------------------------------------------------------------
# `DispersionViewerGUI.onMapClicked` resolves a click (via `GuiCommon.nearestLinePoint`) as ONE
# joint nearest-point search across
# EVERY branch's plotted line, in AXIS-RANGE-normalised distance (`xRange = max(diff(xlim(ax)),
# eps)`) -- exactly `guiNearestPoint`'s own documented purpose ("Nearest plotted point ... across
# several line/point series ... the map viewer's use case"), reused DIRECTLY here rather than
# re-derived the way the angle-sweep viewer's two-stage 1D search needed a bespoke helper
# (`_aswNearestIndex`) instead. `nearestLinePoint` also "Skips HIDDEN lines, not just missing
# ones" -- reproduced via `guiNearestPoint`'s own `visible` mask, read from the CURRENT "Show
# Ritz branches" toggle state for the Layered variant (every branch is hidden/shown together,
# since `Vis.plotDispersionCurves`'s own `showRitzBranches` is a single all-or-nothing flag --
# see `plotDispersionCurves`'s own `visible = showRitzBranches` kwarg, `PlotsDiagnostic.jl`);
# always fully visible for the Plate variant (no such toggle exists there at all).
#
# `guiNearestPoint` needs the plotted (x,y) DATA for every branch, and `plotDispersionCurves`/
# `plotDispersionCurvesPlate` return only the `Axis`, not those
# arrays -- and the private `_diagDispersionCurveAxes` helper that reshuffles
# `(wavelength,k,vel,velGroup,freq)` into "whichever goes on x, whichever on y, for this
# `curveType`" is `_diag`-prefixed, i.e. PRIVATE to `PlotsDiagnostic.jl` (same "do not reach into
# another file's `_`-prefixed helpers" rule as above). So this file duplicates that ~15-line
# reshuffle as its own `_dvCurveAxes` (verbatim logic, `_diag`-prefix dropped for `_dv`), fed by
# the PUBLIC `dispersionCurvesData` (Layered) / a small duplicate of `plotDispersionCurvesPlate`'s
# own public-API-only branch-value assembly (`_dvPlateCurveData`, built from the PUBLIC
# `classifyBranches`) -- both exported, safe to call directly. The resulting
# `[nWav x nCompModes]` `xAllCache`/`yAllCache` pair is recomputed and cached on every map
# redraw, used for BOTH click resolution and the crosshair's own (wavNo, branchIdx) -> (x,y)
# lookup (mirroring MATLAB's own `syncCrosshair`, which reads the already-plotted line's
# `XData`/`YData` directly rather than recomputing physics a second time).
#
# ---------------------------------------------------------------------------------------
# REENTRANCY GUARD -- the pattern documented in `ext/GuiModalShapeViewer.jl`
# ---------------------------------------------------------------------------------------
# `goTo!` (this file's `DispersionViewerGUI.goTo` analogue) syncs the wavelength dropdown AND
# the branch field to the new state on every call -- including when `goTo!` was ITSELF invoked
# by one of those two widgets' own change listener. A Makie `Observable` fires on EVERY
# assignment unconditionally (unlike MATLAB's `.Value`, which never re-fires its own
# `ValueChangedFcn` from a programmatic set), so without a guard this would recurse. The
# `updatingWidgets::Base.RefValue{Bool}` pattern documented in `ext/GuiModalShapeViewer.jl`'s
# header is reused here unchanged.
#
# ---------------------------------------------------------------------------------------
# DELIBERATE DEVIATIONS FROM THE MATLAB SOURCE (documented, not oversights)
# ---------------------------------------------------------------------------------------
#   1. NO `windowMode`/`parentContainer` SPLIT -- see "THE SHAPE GRID OWNS ITS OWN AXES"
#      above. `fig` keyword: `nothing` (default) for a standalone window, or a
#      `GridPosition` to embed into (matching every other viewer's own `fig` convention).
#   2. `<<`/`>>` STEPPER BUTTONS on the curve-type/colour-mode dropdowns -- same added
#      convenience as every other viewer file's own deviation 1, no MATLAB counterpart to
#      contradict (the Wavelength/Wavevector and Branch rows already had `<<`/`>>` in MATLAB;
#      only the curve-type/colour-mode dropdowns gain NEW stepper buttons here).
#   3. MATLAB's `uicheckbox` ("Show Ritz branches"/"Show Green background"/"Show Christoffel
#      bounds") becomes a Makie `Toggle` (the same widget `guiRebuildTypeToggles!` uses
#      for per-type visibility) -- functionally identical, no MATLAB behaviour lost. Never
#      "greyed out"/disabled the way MATLAB's `showGreenBgCheckbox.Enable` toggles: same
#      Makie limitation as the angle-sweep viewer's deviation 2 (`Toggle` exposes no
#      `enabled` attribute in this Makie version) -- the redraw recomputes the EFFECTIVE
#      showGreenBackground value from scratch every time (`hasGreenData && toggle.active[] &&
#      curveType in VPHASE_TYPES`), exactly mirroring MATLAB's own comment on
#      `showGreenBgCheckbox`'s redraw-time recompute ("Recomputed here from scratch ... so this
#      call is correct regardless of whether Enable has been synced").
#   4. AXIS-RANGE-FIELD LABEL TEXT STAYS STATIC ("X min/max"/"Y min/max"), unlike MATLAB's own
#      `obj.xRangeLabel.Text = obj.mapAxes.XLabel.String` re-labelling on every redraw (X/Y
#      meaning genuinely changes with `curveType` here, unlike every other viewer).
#      MATLAB needed this because its own side panel sits visually apart from the axis inside a
#      uifigure; in this Makie layout the axis's OWN `xlabel`/`ylabel` are drawn directly next to
#      the fields already (both live in the same `GridLayout` region), so duplicating that text
#      into the field's label adds little. `guiRangeFieldPair` (`GuiWidgets.jl`) has no hook to
#      mutate a mounted label's text after construction, and adding one would change a shared
#      API used by three other viewer files. The NUMERIC min/max VALUES, functionally important,
#      ARE still refreshed every redraw (see below) -- only the descriptive label text is static.
#   5. NUMERIC RANGE-FIELD DEFAULTS ARE STILL RE-DERIVED EVERY REDRAW, faithfully porting
#      MATLAB's own `xMinField.Value = num2str(xl(1))` etc. -- done here by writing the freshly
#      computed `(min,max)` of `xAllCache`/`yAllCache` directly into the MOUNTED `Textbox`'s own
#      `.stored_string[]` (not `GuiRangeField.loText[]`/`hiText[]` directly), the
#      technique `GuiSingleAngleViewer.jl`'s isoFreq variant uses for its
#      analogous "re-default the Value range on every scale-mode switch" case
#      (`yLoBox.stored_string[] = "-5"`) -- this both updates the VISIBLE field text and cascades
#      through the SAME already-wired `Textbox.stored_string -> GuiRangeField.loText/hiText ->
#      .range -> guiBindRangeToAxis!'s xlims!/ylims!` listener chain bound once at construction,
#      with no extra manual `xlims!`/`ylims!` call needed in the redraw body itself.
#   6. NO `OnModeSelected`/`OnPointSelected`-STYLE PUBLIC CALLBACK -- absent from BOTH
#      MATLAB sources (no such property exists on either class): this viewer is
#      self-contained (map and shape are shown together, in the SAME window), unlike the
#      angle-sweep, single-angle, and modal-shape viewers, which exist partly to notify an
#      external drill-down/composite parent. Nothing here is dropped; there was nothing to port.
#   7. NO SPLIT-BY-VARIANT SHARED INTERNAL BUILDER -- like `ext/GuiAngSweepViewer.jl`'s
#      3-variant design (not `ext/GuiSingleAngleViewer.jl`'s/`ext/GuiModalShapeViewer.jl`'s
#      shared-`_xxxBuildBase`-plus-thin-wrapper style), the Layered/Plate methods below are each
#      fully self-contained rather than factored through one shared internal builder function --
#      as in `GuiAngSweepViewer.jl`, for a pair/trio whose top-bar widget SET
#      genuinely differs (branches+Green-background vs colour-mode+Christoffel-bounds), even
#      though the surrounding shell (shape grid, wavelength/branch steppers, map content, status
#      lines) is identical -- easier to read as two complete pictures than as one heavily
#      parameterised function. Small pure/drawing helpers ARE shared (`_dvCurveAxes`,
#      `_dvBuildShapeAxes!`, `_dvRedrawShapeGrid!`, `_dvNearestBranchPoint`), matching every
#      other viewer's own "share pure helpers, duplicate the shell" convention.
# =========================================================================================

const _DV_DIR_NAMES = ("u_x (propag.dir.)", "u_y (in-plane)", "u_z (out-of-plane)")
const _DV_CURVE_TYPE_ITEMS = ["wavelength-freq", "wavelength-vphase", "k-freq", "k-vphase",
                               "freq-vphase", "freq-vgroup"]
const _DV_CURVE_TYPE_SYMBOLS = (:wavelength_freq, :wavelength_vphase, :k_freq, :k_vphase,
                                 :freq_vphase, :freq_vgroup)
const _DV_VPHASE_CURVE_TYPES = (:wavelength_vphase, :k_vphase)

# -----------------------------------------------------------------------------------------
# Shared pure helpers (`_dv`-prefixed -- private to this file, see this file's header,
# deviation 7, and the general cross-file `_xxx`-privacy convention every viewer file
# follows)
# -----------------------------------------------------------------------------------------

"""
    _dvCurveAxes(curveType, wavelength, k, vel, velGroup, freq,
                 wavelengthLabel, kLabel, freqLabel) -> (xAll, yAll, xLabel, yLabel)

Map a `curveType` onto the x/y data and axis labels of a dispersion-curve plot -- a deliberate
FILE-LOCAL DUPLICATE of `PlotsDiagnostic.jl`'s own private `_diagDispersionCurveAxes` (see this
file's header, "CLICK-TO-SELECT", for why calling that helper directly is not done). Used here
only to CACHE the plotted branch-line data for click-resolution/crosshair lookups -- the actual
on-screen drawing is always done by the PUBLIC `plotDispersionCurves`/
`plotDispersionCurvesPlate`.
"""
function _dvCurveAxes(curveType::Symbol, wavelength::AbstractVector{<:Real},
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
        throw(ArgumentError("_dvCurveAxes: unknown curveType :$curveType"))
    end
end

"""
    _dvPlateCurveData(sweepResult::CalcWavevecAngSweep2Plate, angleIdx::Integer)
        -> (branchModeNoMap, wavelength, wavVecs, vel, velGroup, freq)

FILE-LOCAL DUPLICATE of `plotDispersionCurvesPlate`'s own internal branch-value assembly
(`PlotsDiagnostic.jl`), built from the PUBLIC `classifyBranches`/`DispersionBranches.track` only
-- see this file's header, "CLICK-TO-SELECT".
"""
function _dvPlateCurveData(sweepResult::CalcWavevecAngSweep2Plate, angleIdx::Integer)
    branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx)
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
    return branchModeNoMap, wavelength, wavVecs, vel, velGroup, freq
end

"""
    _dvItemStrs(xs) -> Vector{String}

`"%.4g"`-equivalent labels for the wavelength/wavevector dropdown (MATLAB's
`sprintf('%.4g', w)`) -- display text only, the real selection is by `i_selected` index, same
convention `GuiModalShapeViewer.jl`'s own `_msvAngleItems` already established.
"""
_dvItemStrs(xs::AbstractVector{<:Real}) = [string(round(x; sigdigits = 4)) for x in xs]

"""
    _dvColLimit(dispField) -> Float64

Symmetric colour-range half-width: the actual max-abs value across every direction/layer of
`dispField` (a `(3, nLay)` matrix of real 2D fields) -- FILE-LOCAL DUPLICATE of
`PlotsModalShape.jl`'s own private `_modalShapeColLimit` (see this file's header, "THE
SHAPE GRID OWNS ITS OWN AXES").
"""
function _dvColLimit(dispField::AbstractMatrix{<:AbstractMatrix{<:Real}})
    m = 0.0
    for field2D in dispField
        isempty(field2D) && continue
        m = max(m, maximum(abs, field2D))
    end
    return max(m, eps())
end

"""
    _dvShapeColormap() -> Vector{<:Makie.Colorant}

`ColormapRWB`'s `N x 3` RGB matrix, converted to the row-of-colours form Makie's `colormap`
keyword expects -- same technique `PlotsModalShape.jl`'s own `_modalShapeColormap` uses.
"""
function _dvShapeColormap()
    cm = ColormapRWB()
    return [Makie.RGBf(cm[i, 1], cm[i, 2], cm[i, 3]) for i in axes(cm, 1)]
end

"""
    _dvBuildShapeAxes!(gl::GridLayout) -> (axes::Vector{Axis}, cbarGL::GridLayout)

Build the 3 direction-column `Axis` objects ONCE (decorations/spines hidden, titles set) plus a
small clearable `GridLayout` cell for the colourbar -- see this file's header, "THE SHAPE GRID
OWNS ITS OWN AXES". The 3 axes are reused via `empty!(ax)` on every subsequent
redraw ([`_dvRedrawShapeGrid!`](@ref)); only `cbarGL` is cleared+rebuilt each time (a fresh
`Colorbar` for the new `climLimit`).
"""
function _dvBuildShapeAxes!(gl::GridLayout)
    axes_ = Vector{Axis}(undef, 3)
    for i in 1:3
        ax = Axis(gl[1, i]; title = _DV_DIR_NAMES[i])
        hidedecorations!(ax)
        hidespines!(ax)
        axes_[i] = ax
    end
    cbarGL = GridLayout(gl[1, 4])   # natural (colour-bar) width, not guiGridLayout
    return axes_, cbarGL
end

"""
    _dvRedrawShapeGrid!(axes_::Vector{Axis}, cbarGL::GridLayout, ritz, modeNo::Integer) -> Nothing

Redraw the depth-stacked modal-shape grid for `ritz`'s mode `modeNo` at `'boosted'` substrate
cutoff (the ONE mode this viewer ever needs, per both MATLAB sources' own literal
`Vis.plotModalShape(ritz, modeNo, 'boosted', [], 1:3, ...)` call) into the 3 pre-built `axes_`
(one per direction; cleared via `empty!` first) plus a freshly-rebuilt colourbar in `cbarGL`.
Re-implements `PlotsModalShape.jl`'s own `_modalShapeDrawColumn!`/`_modalShapeLayeredPlot`
drawing technique directly (one `Axis` per direction, each layer its own `heatmap!` band with a
thin frame, row heights from `calcModalShapeRowHeights`, symmetric-about-zero colour range) --
see this file's header for why this is a deliberate duplication, not a call into that file's own
private helpers. `'boosted'` never crops the substrate band (unlike `:decay`/`:fixedMultiple`),
so no substrate-cropping step is needed here at all.
"""
function _dvRedrawShapeGrid!(axes_::Vector{Axis}, cbarGL::GridLayout, ritz, modeNo::Integer)
    relDepths = [cp.relDepth for cp in ritz.structure.calcParams]
    dispField = calcDispField2D(ritz, modeNo, 0.01)   # (3, nLay)
    climLimit = _dvColLimit(dispField)

    layHeightFrac = calcModalShapeRowHeights(relDepths)
    layHeightFrac = layHeightFrac isa Real ? [layHeightFrac] : layHeightFrac
    nLay = length(relDepths)
    layBottom = Vector{Float64}(undef, nLay)
    acc = 0.0
    for ilay in 1:nLay
        layBottom[ilay] = acc
        acc += layHeightFrac[ilay]
    end

    cMap = _dvShapeColormap()
    for iDir in 1:3
        ax = axes_[iDir]
        guiEmptyAxis!(ax)
        for ilay in 1:nLay
            field2D = dispField[iDir, ilay]
            bandLo, bandHi = layBottom[ilay], layBottom[ilay] + layHeightFrac[ilay]
            nMeshDepth, nMeshLateral = size(field2D)
            xNorm = range(0.0, 1.0; length = nMeshLateral)
            yNorm = range(bandLo, bandHi; length = nMeshDepth)
            heatmap!(ax, xNorm, yNorm, permutedims(field2D); colormap = cMap,
                     colorrange = (-climLimit, climLimit))
            lines!(ax, [0.0, 1.0, 1.0, 0.0, 0.0], [bandLo, bandLo, bandHi, bandHi, bandLo];
                   color = :black, linewidth = 1)
        end
        hidedecorations!(ax)
        hidespines!(ax)
        xlims!(ax, 0.0, 1.0)
        ylims!(ax, 0.0, 1.0)
    end

    guiClearGridLayout!(cbarGL)
    Colorbar(cbarGL[1, 1]; colorrange = (-climLimit, climLimit), colormap = cMap,
             label = "displacement (normalised)")
    return nothing
end

"""
    _dvNearestBranchPoint(xAll, yAll, x, y, xSpan, ySpan, visible) -> Union{Nothing,Tuple{Int,Int}}

Pure core of `onMapClicked!`: nearest plotted (wavelength/wavevector-step, branch) point across
EVERY branch's line, in axis-range-normalised distance -- routes through
[`guiNearestPoint`](@ref) directly (see this file's header, "CLICK-TO-SELECT", for why this
viewer's click resolution is a JOINT 2D search across several series, `guiNearestPoint`'s own
purpose-built case, unlike the angle-sweep viewer's two independent 1D searches). Returns `(wavNo,
branchIdx)`, or `nothing` if `xAll` is empty or every branch is hidden.
"""
function _dvNearestBranchPoint(xAll::AbstractMatrix{<:Real}, yAll::AbstractMatrix{<:Real},
                                x::Real, y::Real, xSpan::Real, ySpan::Real,
                                visible::AbstractVector{Bool})
    isempty(xAll) && return nothing
    nBranch = size(xAll, 2)
    seriesXs = [view(xAll, :, im) for im in 1:nBranch]
    seriesYs = [view(yAll, :, im) for im in 1:nBranch]
    hit = guiNearestPoint(seriesXs, seriesYs, (Float64(x), Float64(y));
                           xrange = xSpan, yrange = ySpan, visible = visible)
    hit === nothing && return nothing
    return (hit.index, hit.series)   # (wavNo, branchIdx)
end

"""
    _dvAxisSpan(rf::GuiRangeField) -> Float64

Current axis span tracked by `rf` (the same value [`guiBindRangeToAxis!`](@ref) applies to the
axis), or `1.0` if not yet set -- used to normalise [`_dvNearestBranchPoint`](@ref)'s distance
metric exactly the way MATLAB's own `nearestLinePoint` normalises by `diff(xlim(ax))`/
`diff(ylim(ax))`. Reading the tracked `GuiRangeField` (rather than re-deriving the axis's own
CURRENT rendered limits, e.g. via `ax.finallimits[]`) is deliberate: this file always keeps the
axis's explicit limits and `GuiRangeField.range` in lock-step (see this file's header, deviation
5), so the two are equivalent by construction, and using the field avoids depending on Makie's
own auto-limit computation having already run a render pass (unreliable
immediately after a headless plot call with no explicit `xlims!`/`ylims!`).
"""
function _dvAxisSpan(rf::GuiRangeField)
    r = rf.range[]
    return r === nothing ? 1.0 : max(r[2] - r[1], eps())
end

"""
    _dvToggleRow(parent, row::Integer, labelText::AbstractString; active::Bool = true) -> Toggle

Label + `Toggle` (Makie's checkbox-equivalent widget), one row -- the replacement for MATLAB's
`uicheckbox` (see this file's header, deviation 3). Mirrors `GuiWidgets.jl`'s own
`guiDropdownStepperRow`'s "own nested `GridLayout` at `parent[row,1]`" construction exactly, just
with a `Toggle` instead of a `Menu`+stepper buttons.
"""
function _dvToggleRow(parent, row::Integer, labelText::AbstractString; active::Bool = true)
    rowGL = guiGridLayout(parent[row, 1])
    Label(rowGL[1, 1], labelText; halign = :left, tellwidth = false)
    return Toggle(rowGL[1, 2], active = active)
end

"""
    _dvWireMenuStepper!(menu::Menu, prevBtn::Button, nextBtn::Button, n::Integer) -> Nothing

Wire a dropdown-stepper row's `<<`/`>>` buttons to cycle `menu.i_selected` through `1:n`,
wrapping via `mod1` -- private near-duplicate of every other viewer file's own
`_xxxWireStepper!`/`_xxxWireMenuStepper!` helper (`_aswWireStepper!`, `_savWireStepper!`,
`_msvWireMenuStepper!`), kept as this file's own copy per the same "never depend on another
viewer file's private internals" convention.
"""
function _dvWireMenuStepper!(menu, prevBtn, nextBtn, n::Integer)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] - 1, n)), prevBtn.clicks)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] + 1, n)), nextBtn.clicks)
    return nothing
end

# -----------------------------------------------------------------------------------------
# DispersionViewer -- the one handle type returned by both variants
# -----------------------------------------------------------------------------------------

"""
    DispersionViewer

Handles for one "map + shape grid in one window, with branch tracking" viewer,
returned by both [`dispersionViewer`](@ref) methods (see this file's header for the dispatch
design).

Fields:
- `fig` -- the ONE owning `Figure` (embeddable via the `fig` keyword; see this file's header,
  "THE SHAPE GRID OWNS ITS OWN AXES", for why this viewer needs only one,
  unlike the modal-shape viewer's two).
- `mapAxes` -- the dispersion-curve map `Axis`, built once, reused via `empty!` every redraw.
- `shapeAxes` -- the 3 depth-stacked modal-shape `Axis` objects (one per direction), likewise
  built once and reused.
- `shapeCbarGL` -- the shape grid's colourbar cell, cleared+rebuilt every redraw.
- `result` -- the `CalcWavevecAngSweep2Layered`/`CalcWavevecAngSweep2Plate` this viewer displays
  (untyped field, kept generic so this ONE struct serves both variants, matching every other
  viewer's own `result`/`cs` precedent).
- `angleIdx` -- the fixed angle index this whole view is pinned at (constructor argument,
  default 1; see this file's header).
- `branchModeNoMap` -- `[nWav x nCompModes]`, [`DispersionBranches.track`](@ref)'s own output,
  computed ONCE at construction (mirrors MATLAB's own `SetAccess = private` property of the same
  name).
- `wavMenu`, `wavStepperBtns` -- the wavelength/wavevector dropdown-stepper row.
- `branchBox`, `branchStepperBtns`, `velLabel` -- the Branch index-stepper row (`velLabel` is
  the trailing "vel: ... km/s, freq: ..." label).
- `curveTypeMenu`, `curveTypeStepperBtns` -- the curve-type dropdown-stepper row (every variant).
- `showBranchesToggle`, `showGreenBgToggle` -- the Layered variant's two checkboxes (Makie
  `Toggle`); `nothing` for Plate, which has neither concept (absent, not disabled, matching
  every other viewer's own precedent for a variant-specific widget).
- `colorModeMenu`, `colorModeStepperBtns`, `showChristoffelToggle` -- the Plate variant's own
  colour-mode dropdown + Christoffel-bounds checkbox; `nothing` for Layered.
- `xRangeField`/`yRangeField` -- persistent axis-range state (`GuiRangeField`); `xRangeBoxes`/
  `yRangeBoxes` -- the mounted `(loBox, hiBox)` `Textbox` pairs (re-defaulted every redraw, see
  this file's header, deviation 5).
- `infoLabel` -- the branch/mode info panel (map side panel).
- `shapeStatusLabel`/`mapStatusLabel` -- two independent ok/error status lines, mirroring
  MATLAB's own `shapeStatusLine`/`mapStatusLine` split.
- `crosshair` -- the map's selection-crosshair `Observable{Point2f}` (built once, moved via
  [`_dvSyncCrosshair!`](@ref); no drill-down/content-teardown in this viewer, so unlike
  the angle-sweep viewer's `Union{Observable,Nothing}` field this is never destroyed).
- `currentWavNo`/`currentBranchIdx` -- current selection (`Int`, mutated only through `goTo!`).
- `xAllCache`/`yAllCache` -- the current curve-type's `[nWav x nCompModes]` plotted branch data,
  recomputed every map redraw ([`_dvCurveAxes`](@ref)) and used for both click resolution and
  the crosshair lookup (see this file's header, "CLICK-TO-SELECT").
- `hasGreenData` -- Layered variant only (`sweeps[1].greenAbsG13EffMap` non-empty, i.e. built with
  `calcMode = :full`); always `false` (unused) for Plate.
- `updatingWidgets` -- reentrancy guard (see this file's header); not meant to be read/set by
  callers, exposed only because every field on every viewer struct is.
- `redrawMap!` -- `() -> Nothing`, redraws the dispersion-curve map for the current curve-type/
  toggle selection and refreshes `xAllCache`/`yAllCache`/the range-field defaults; also callable
  directly (e.g. from a test).
- `redrawShape!` -- `() -> Nothing`, redraws the shape grid for the current (wavNo, branchIdx)
  selection.
- `goTo!` -- `(wavNo, branchIdx) -> Nothing`, MATLAB's `goTo` -- jumps directly to a selection,
  clamping both indices, syncing every widget, and redrawing the shape grid + crosshair (NOT the
  map -- matches MATLAB's own `goTo`/`syncControlsAndRedraw` exactly).
- `onMapClicked!` -- `(x::Float64, y::Float64) -> Nothing`, the click handler
  [`guiOnAxisClick`](@ref) is wired to on `mapAxes`; also callable directly with a known
  data-space point (matching every other viewer's own convention).
"""
mutable struct DispersionViewer
    fig::Figure
    mapAxes::Axis
    shapeAxes::Vector{Axis}
    shapeCbarGL::GridLayout

    result::Any
    angleIdx::Int
    branchModeNoMap::Matrix{Int}

    wavMenu::Menu
    wavStepperBtns::Tuple{Button,Button}
    branchBox::Textbox
    branchStepperBtns::Tuple{Button,Button}
    velLabel::Label

    curveTypeMenu::Menu
    curveTypeStepperBtns::Tuple{Button,Button}

    showBranchesToggle::Union{Toggle,Nothing}
    showGreenBgToggle::Union{Toggle,Nothing}

    colorModeMenu::Union{Menu,Nothing}
    colorModeStepperBtns::Union{Tuple{Button,Button},Nothing}
    showChristoffelToggle::Union{Toggle,Nothing}

    xRangeField::GuiRangeField
    yRangeField::GuiRangeField
    xRangeBoxes::Tuple{Textbox,Textbox}
    yRangeBoxes::Tuple{Textbox,Textbox}

    infoLabel::Label
    shapeStatusLabel::Label
    mapStatusLabel::Label

    crosshair::Observable{Point2f}

    currentWavNo::Int
    currentBranchIdx::Int
    xAllCache::Matrix{Float64}
    yAllCache::Matrix{Float64}

    hasGreenData::Bool

    updatingWidgets::Base.RefValue{Bool}

    redrawMap!::Function
    redrawShape!::Function
    goTo!::Function
    onMapClicked!::Function
end

# -----------------------------------------------------------------------------------------
# dispersionViewer -- Layered variant (CalcWavevecAngSweep2Layered)
# -----------------------------------------------------------------------------------------

"""
    dispersionViewer(sweepResult::CalcWavevecAngSweep2Layered, angleIdx::Integer = 1;
                      fig = nothing) -> DispersionViewer

Map + shape-grid dispersion viewer for an already-`calc()`'d `CalcWavevecAngSweep2Layered`
result -- port of `guliwerr/gui/DispersionViewerGUI.m`. `fig` -- `nothing` (default) for a
standalone window, or a `GridPosition` to embed into (see this file's header, deviation 1, for
why there is no `windowMode`/`parentContainer` split).

Requires `retentionLevel = :full` ([`DispersionBranches.track`](@ref)'s own requirement, called
directly here at construction -- its `ArgumentError` propagates unchanged if not satisfied).
Throws `ArgumentError` if `sweepResult.sweeps` is empty (not yet `calc()`'d) or `angleIdx` is
out of range.

Top bar: curve-type dropdown (6 types, matching [`visDispersionCurveTypes`](@ref)'s own
vocabulary), "Show Ritz branches" (default on) and "Show Green background" (default off,
effective only for the 2 v_phase curve types AND only when `sweepResult` was built with
`calcMode = :full`) toggles, and a static "Angle: X deg" label.
"""
function GuliwerrJ.dispersionViewer(sweepResult::CalcWavevecAngSweep2Layered, angleIdx::Integer = 1;
                           fig = nothing)
    isempty(sweepResult.sweeps) && throw(ArgumentError(
        "dispersionViewer: sweepResult has not been calc()'d yet (sweeps is empty) -- " *
        "nothing to view."))
    (1 <= angleIdx <= length(sweepResult.angles)) || throw(ArgumentError(
        "dispersionViewer: angleIdx must be an integer in 1:$(length(sweepResult.angles)), " *
        "got $angleIdx."))

    branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx)
    hasGreenData = !isempty(sweepResult.sweeps[1].greenAbsG13EffMap)

    fig, mainLayout = guiHostLayout(fig; size = (1200, 900))
    shapeGL, shapeControlsGL, bottomGL, _ = guiDividedTopBottom(mainLayout[1, 1])

    shapeContentGL = guiGridLayout(shapeGL[1, 1])
    shapeAxes, shapeCbarGL = _dvBuildShapeAxes!(shapeContentGL)

    wavMenu, wavPrev, wavNext = guiDropdownStepperRow(shapeControlsGL, 1, "Wavelength:",
                                                       _dvItemStrs(sweepResult.wavelengths))
    branchBox, branchPrev, branchNext, velLabel =
        guiIndexStepperRow(shapeControlsGL, 2, "Branch:")
    shapeStatusLabel = guiStatusLabel(shapeControlsGL[3, 1])

    mapMainLayout = guiGridLayout(bottomGL[1, 1])
    topBarGL = guiGridLayout(mapMainLayout[1, 1])
    curveTypeMenu, ctPrev, ctNext =
        guiDropdownStepperRow(topBarGL, 1, "Curve type:", _DV_CURVE_TYPE_ITEMS)
    showBranchesToggle = _dvToggleRow(topBarGL, 2, "Show Ritz branches"; active = true)
    showGreenBgToggle = _dvToggleRow(topBarGL, 3, "Show Green background"; active = false)
    Label(topBarGL[4, 1], "Angle: $(round(sweepResult.angles[angleIdx]; sigdigits = 4)) deg";
          halign = :left, tellwidth = false)

    mapAxes, sideGL, _ = guiPlotWithSidePanel(mapMainLayout[2, 1])
    xRangeField, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "X min/max")
    yRangeField, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "Y min/max")
    infoLabel = Label(sideGL[5, 1], "Click the map to select a branch.";
                       valign = :top, halign = :left)
    mapStatusLabel = guiStatusLabel(mapMainLayout[3, 1])

    state = DispersionViewer(fig, mapAxes, shapeAxes, shapeCbarGL,
        sweepResult, angleIdx, branchModeNoMap,
        wavMenu, (wavPrev, wavNext), branchBox, (branchPrev, branchNext), velLabel,
        curveTypeMenu, (ctPrev, ctNext),
        showBranchesToggle, showGreenBgToggle,
        nothing, nothing, nothing,
        xRangeField, yRangeField, (xLoBox, xHiBox), (yLoBox, yHiBox),
        infoLabel, shapeStatusLabel, mapStatusLabel,
        Observable(Point2f(NaN, NaN)),
        1, 1, zeros(Float64, 0, 0), zeros(Float64, 0, 0),
        hasGreenData, Ref(false),
        () -> nothing, () -> nothing, (_, _) -> nothing, (_, _) -> nothing)

    guiBindRangeToAxis!(mapAxes, state.xRangeField, :x)
    guiBindRangeToAxis!(mapAxes, state.yRangeField, :y)
    drawnCurveType = Ref(:none)   # range fields refill from the data only when curveType changes

    function _redrawMap!()
        ax = state.mapAxes
        guiEmptyAxis!(ax)
        curveType = _DV_CURVE_TYPE_SYMBOLS[state.curveTypeMenu.i_selected[]]
        showBranches = state.showBranchesToggle.active[]
        showBg = state.hasGreenData && state.showGreenBgToggle.active[] &&
                 curveType in _DV_VPHASE_CURVE_TYPES
        ok, errMsg = true, ""
        try
            plotDispersionCurves(state.result, curveType, state.angleIdx, showBranches, showBg;
                                  ax = ax)
            _, wavelengths, k, vel, velGroup, freq, _ =
                dispersionCurvesData(state.result, state.angleIdx)
            xAll, yAll, _, _ = _dvCurveAxes(curveType, wavelengths, k, vel, velGroup, freq,
                                             "wavelength", "k [2π/wavelength]",
                                             "frequency (= v_φ / wavelength)")
            state.xAllCache = xAll
            state.yAllCache = yAll
            if curveType !== drawnCurveType[]
                xLoBox2, xHiBox2 = state.xRangeBoxes
                yLoBox2, yHiBox2 = state.yRangeBoxes
                xLoBox2.stored_string[] = string(minimum(xAll))
                xHiBox2.stored_string[] = string(maximum(xAll))
                yLoBox2.stored_string[] = string(minimum(yAll))
                yHiBox2.stored_string[] = string(maximum(yAll))
                drawnCurveType[] = curveType
            else   # display toggle only -- keep the set range
                r = state.xRangeField.range[]; r === nothing || xlims!(ax, r)
                r = state.yRangeField.range[]; r === nothing || ylims!(ax, r)
            end
        catch err
            ok, errMsg = false, sprint(showerror, err)
            state.xAllCache = zeros(Float64, 0, 0)
            state.yAllCache = zeros(Float64, 0, 0)
        end
        ok ? guiSetStatusOk!(state.mapStatusLabel, "") :
             guiSetStatusError!(state.mapStatusLabel, errMsg)
        return nothing
    end
    state.redrawMap! = _redrawMap!

    function _redrawShape!()
        iw, ib = state.currentWavNo, state.currentBranchIdx
        ok, errMsg = true, ""
        try
            ritz = state.result.sweeps[iw].cs[state.angleIdx].ritzSolver
            modeNo = state.branchModeNoMap[iw, ib]
            _dvRedrawShapeGrid!(state.shapeAxes, state.shapeCbarGL, ritz, modeNo)

            wav = state.result.wavelengths[iw]
            vel = ritz.vels[modeNo]
            state.velLabel.text[] =
                "vel: $(round(vel, digits = 4)) km/s, freq: $(round(vel / wav, sigdigits = 4))"

            typeStr = string(guiSignificantModeType(
                state.result.sweeps[iw].significantModesMap[state.angleIdx], modeNo))
            state.infoLabel.text[] = join([
                "wavelength = $(round(wav; sigdigits = 4))",
                "branch # = $ib (mode $modeNo)",
                "velocity = $(round(vel, digits = 4)) km/s",
                "type = $typeStr"], "\n")
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end
        ok ? guiSetStatusOk!(state.shapeStatusLabel, "") :
             guiSetStatusError!(state.shapeStatusLabel, errMsg)
        return nothing
    end
    state.redrawShape! = _redrawShape!

    function _syncWidgets!()
        state.updatingWidgets[] = true
        try
            state.wavMenu.i_selected[] = state.currentWavNo
            state.branchBox.stored_string[] = string(state.currentBranchIdx)
        finally
            state.updatingWidgets[] = false
        end
        return nothing
    end

    function _syncCrosshair!()
        im = state.currentBranchIdx
        iw = state.currentWavNo
        if !isempty(state.xAllCache) && im <= size(state.xAllCache, 2) && iw <= size(state.xAllCache, 1)
            state.crosshair[] = Point2f(state.xAllCache[iw, im], state.yAllCache[iw, im])
        end
        return nothing
    end

    function _goTo!(wavNo::Real, branchIdx::Real)
        nWav, nBranch = size(state.branchModeNoMap)
        state.currentWavNo = clamp(round(Int, wavNo), 1, nWav)
        state.currentBranchIdx = clamp(round(Int, branchIdx), 1, nBranch)
        _syncWidgets!()
        state.redrawShape!()
        _syncCrosshair!()
        return nothing
    end
    state.goTo! = _goTo!

    function _onMapClicked!(x::Float64, y::Float64)
        nBranch = size(state.branchModeNoMap, 2)
        visMask = fill(state.showBranchesToggle.active[], nBranch)
        hit = _dvNearestBranchPoint(state.xAllCache, state.yAllCache, x, y,
                                     _dvAxisSpan(state.xRangeField), _dvAxisSpan(state.yRangeField),
                                     visMask)
        hit === nothing && return nothing
        wavNo, branchIdx = hit
        state.goTo!(wavNo, branchIdx)
        return nothing
    end
    state.onMapClicked! = _onMapClicked!

    on(idx -> (state.updatingWidgets[] || state.goTo!(idx, state.currentBranchIdx)),
       wavMenu.i_selected)
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo - 1, state.currentBranchIdx)),
       wavPrev.clicks)
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo + 1, state.currentBranchIdx)),
       wavNext.clicks)
    on(branchBox.stored_string) do s
        state.updatingWidgets[] && return nothing
        v = tryparse(Int, s)
        v === nothing || state.goTo!(state.currentWavNo, v)
        return nothing
    end
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo, state.currentBranchIdx - 1)),
       branchPrev.clicks)
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo, state.currentBranchIdx + 1)),
       branchNext.clicks)

    function _onCurveTypeChanged!()
        state.redrawMap!()
        _syncCrosshair!()
        return nothing
    end
    on(_ -> _onCurveTypeChanged!(), curveTypeMenu.i_selected)
    on(_ -> _onCurveTypeChanged!(), showBranchesToggle.active)
    on(_ -> _onCurveTypeChanged!(), showGreenBgToggle.active)
    _dvWireMenuStepper!(curveTypeMenu, ctPrev, ctNext, length(_DV_CURVE_TYPE_ITEMS))

    guiOnAxisClick(state.onMapClicked!, mapAxes)

    state.redrawMap!()
    state.goTo!(1, 1)

    return state
end

# -----------------------------------------------------------------------------------------
# dispersionViewer -- Plate variant (CalcWavevecAngSweep2Plate)
# -----------------------------------------------------------------------------------------

"""
    dispersionViewer(sweepResult::CalcWavevecAngSweep2Plate, angleIdx::Integer = 1;
                      fig = nothing) -> DispersionViewer

Map + shape-grid dispersion viewer for an already-`calc()`'d `CalcWavevecAngSweep2Plate`
result -- port of `guliwerr/gui/DispersionViewerGUIPlate.m`. Same shell as the Layered method
above; differs only in: `wavVecs` (not `wavelengths`) drives the top stepper, no Green's-
response/significant-modes concept at all (`plotDispersionCurvesPlate`'s own `colorMode`
dropdown + Christoffel-bounds toggle replace the branches/Green-background pair), and
`ritz.freqsPhysical` is read directly for frequency rather than re-derived from `vel/wavelength`
(`RitzSolverLamb`'s own convention).

Requires `retentionLevel = :full`. Throws `ArgumentError` if `sweepResult.sweeps` is empty or
`angleIdx` is out of range.
"""
function GuliwerrJ.dispersionViewer(sweepResult::CalcWavevecAngSweep2Plate, angleIdx::Integer = 1;
                           fig = nothing)
    isempty(sweepResult.sweeps) && throw(ArgumentError(
        "dispersionViewer: sweepResult has not been calc()'d yet (sweeps is empty) -- " *
        "nothing to view."))
    (1 <= angleIdx <= length(sweepResult.angles)) || throw(ArgumentError(
        "dispersionViewer: angleIdx must be an integer in 1:$(length(sweepResult.angles)), " *
        "got $angleIdx."))

    branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx)

    fig, mainLayout = guiHostLayout(fig; size = (1200, 900))
    shapeGL, shapeControlsGL, bottomGL, _ = guiDividedTopBottom(mainLayout[1, 1])

    shapeContentGL = guiGridLayout(shapeGL[1, 1])
    shapeAxes, shapeCbarGL = _dvBuildShapeAxes!(shapeContentGL)

    wavMenu, wavPrev, wavNext = guiDropdownStepperRow(shapeControlsGL, 1, "Wavevector:",
                                                       _dvItemStrs(sweepResult.wavVecs))
    branchBox, branchPrev, branchNext, velLabel =
        guiIndexStepperRow(shapeControlsGL, 2, "Branch:")
    shapeStatusLabel = guiStatusLabel(shapeControlsGL[3, 1])

    mapMainLayout = guiGridLayout(bottomGL[1, 1])
    topBarGL = guiGridLayout(mapMainLayout[1, 1])
    curveTypeMenu, ctPrev, ctNext =
        guiDropdownStepperRow(topBarGL, 1, "Curve type:", _DV_CURVE_TYPE_ITEMS)
    colorModeMenu, cmPrev, cmNext = guiDropdownStepperRow(topBarGL, 2, "Colour by:", ["type", "branch"])
    showChristoffelToggle = _dvToggleRow(topBarGL, 3, "Show Christoffel bounds"; active = true)
    Label(topBarGL[4, 1], "Angle: $(round(sweepResult.angles[angleIdx]; sigdigits = 4)) deg";
          halign = :left, tellwidth = false)

    mapAxes, sideGL, _ = guiPlotWithSidePanel(mapMainLayout[2, 1])
    xRangeField, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "X min/max")
    yRangeField, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "Y min/max")
    infoLabel = Label(sideGL[5, 1], "Click the map to select a branch.";
                       valign = :top, halign = :left)
    mapStatusLabel = guiStatusLabel(mapMainLayout[3, 1])

    state = DispersionViewer(fig, mapAxes, shapeAxes, shapeCbarGL,
        sweepResult, angleIdx, branchModeNoMap,
        wavMenu, (wavPrev, wavNext), branchBox, (branchPrev, branchNext), velLabel,
        curveTypeMenu, (ctPrev, ctNext),
        nothing, nothing,
        colorModeMenu, (cmPrev, cmNext), showChristoffelToggle,
        xRangeField, yRangeField, (xLoBox, xHiBox), (yLoBox, yHiBox),
        infoLabel, shapeStatusLabel, mapStatusLabel,
        Observable(Point2f(NaN, NaN)),
        1, 1, zeros(Float64, 0, 0), zeros(Float64, 0, 0),
        false, Ref(false),
        () -> nothing, () -> nothing, (_, _) -> nothing, (_, _) -> nothing)

    guiBindRangeToAxis!(mapAxes, state.xRangeField, :x)
    guiBindRangeToAxis!(mapAxes, state.yRangeField, :y)
    drawnCurveType = Ref(:none)   # range fields refill from the data only when curveType changes

    function _redrawMap!()
        ax = state.mapAxes
        guiEmptyAxis!(ax)
        curveType = _DV_CURVE_TYPE_SYMBOLS[state.curveTypeMenu.i_selected[]]
        colorMode = state.colorModeMenu.i_selected[] == 1 ? :type : :branch
        showChrist = state.showChristoffelToggle.active[]
        ok, errMsg = true, ""
        try
            plotDispersionCurvesPlate(state.result, curveType, state.angleIdx, colorMode, 0.3,
                                       showChrist; ax = ax)
            _, wavelength, wavVecs, vel, velGroup, freq =
                _dvPlateCurveData(state.result, state.angleIdx)
            xAll, yAll, _, _ = _dvCurveAxes(curveType, wavelength, wavVecs, vel, velGroup, freq,
                                             "wavelength", "k", "frequency")
            state.xAllCache = xAll
            state.yAllCache = yAll
            if curveType !== drawnCurveType[]
                xLoBox2, xHiBox2 = state.xRangeBoxes
                yLoBox2, yHiBox2 = state.yRangeBoxes
                xLoBox2.stored_string[] = string(minimum(xAll))
                xHiBox2.stored_string[] = string(maximum(xAll))
                yLoBox2.stored_string[] = string(minimum(yAll))
                yHiBox2.stored_string[] = string(maximum(yAll))
                drawnCurveType[] = curveType
            else   # display toggle only -- keep the set range
                r = state.xRangeField.range[]; r === nothing || xlims!(ax, r)
                r = state.yRangeField.range[]; r === nothing || ylims!(ax, r)
            end
        catch err
            ok, errMsg = false, sprint(showerror, err)
            state.xAllCache = zeros(Float64, 0, 0)
            state.yAllCache = zeros(Float64, 0, 0)
        end
        ok ? guiSetStatusOk!(state.mapStatusLabel, "") :
             guiSetStatusError!(state.mapStatusLabel, errMsg)
        return nothing
    end
    state.redrawMap! = _redrawMap!

    function _redrawShape!()
        iw, ib = state.currentWavNo, state.currentBranchIdx
        ok, errMsg = true, ""
        try
            ritz = state.result.sweeps[iw].cs[state.angleIdx].ritzSolver
            modeNo = state.branchModeNoMap[iw, ib]
            _dvRedrawShapeGrid!(state.shapeAxes, state.shapeCbarGL, ritz, modeNo)

            wavVec = state.result.wavVecs[iw]
            vel = ritz.vels[modeNo]
            freq = ritz.freqsPhysical[modeNo]   # RitzSolverLamb's own convention -- NOT
                                                  # re-derived here, matching MATLAB exactly
            state.velLabel.text[] =
                "vel: $(round(vel, digits = 4)) km/s, freq: $(round(freq, sigdigits = 4))"

            typeStr = guiPlateModeTypeText(ritz, modeNo)
            state.infoLabel.text[] = join([
                "wavevector = $(round(wavVec; sigdigits = 4))",
                "branch # = $ib (mode $modeNo)",
                "velocity = $(round(vel, digits = 4)) km/s",
                "type = $typeStr"], "\n")
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end
        ok ? guiSetStatusOk!(state.shapeStatusLabel, "") :
             guiSetStatusError!(state.shapeStatusLabel, errMsg)
        return nothing
    end
    state.redrawShape! = _redrawShape!

    function _syncWidgets!()
        state.updatingWidgets[] = true
        try
            state.wavMenu.i_selected[] = state.currentWavNo
            state.branchBox.stored_string[] = string(state.currentBranchIdx)
        finally
            state.updatingWidgets[] = false
        end
        return nothing
    end

    function _syncCrosshair!()
        im = state.currentBranchIdx
        iw = state.currentWavNo
        if !isempty(state.xAllCache) && im <= size(state.xAllCache, 2) && iw <= size(state.xAllCache, 1)
            state.crosshair[] = Point2f(state.xAllCache[iw, im], state.yAllCache[iw, im])
        end
        return nothing
    end

    function _goTo!(wavNo::Real, branchIdx::Real)
        nWav, nBranch = size(state.branchModeNoMap)
        state.currentWavNo = clamp(round(Int, wavNo), 1, nWav)
        state.currentBranchIdx = clamp(round(Int, branchIdx), 1, nBranch)
        _syncWidgets!()
        state.redrawShape!()
        _syncCrosshair!()
        return nothing
    end
    state.goTo! = _goTo!

    function _onMapClicked!(x::Float64, y::Float64)
        nBranch = size(state.branchModeNoMap, 2)
        visMask = trues(nBranch)   # no branches-visibility toggle for Plate -- always clickable
        hit = _dvNearestBranchPoint(state.xAllCache, state.yAllCache, x, y,
                                     _dvAxisSpan(state.xRangeField), _dvAxisSpan(state.yRangeField),
                                     visMask)
        hit === nothing && return nothing
        wavNo, branchIdx = hit
        state.goTo!(wavNo, branchIdx)
        return nothing
    end
    state.onMapClicked! = _onMapClicked!

    on(idx -> (state.updatingWidgets[] || state.goTo!(idx, state.currentBranchIdx)),
       wavMenu.i_selected)
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo - 1, state.currentBranchIdx)),
       wavPrev.clicks)
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo + 1, state.currentBranchIdx)),
       wavNext.clicks)
    on(branchBox.stored_string) do s
        state.updatingWidgets[] && return nothing
        v = tryparse(Int, s)
        v === nothing || state.goTo!(state.currentWavNo, v)
        return nothing
    end
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo, state.currentBranchIdx - 1)),
       branchPrev.clicks)
    on(_ -> (state.updatingWidgets[] || state.goTo!(state.currentWavNo, state.currentBranchIdx + 1)),
       branchNext.clicks)

    function _onCurveTypeChanged!()
        state.redrawMap!()
        _syncCrosshair!()
        return nothing
    end
    on(_ -> _onCurveTypeChanged!(), curveTypeMenu.i_selected)
    on(_ -> _onCurveTypeChanged!(), colorModeMenu.i_selected)
    on(_ -> _onCurveTypeChanged!(), showChristoffelToggle.active)
    _dvWireMenuStepper!(curveTypeMenu, ctPrev, ctNext, length(_DV_CURVE_TYPE_ITEMS))
    _dvWireMenuStepper!(colorModeMenu, cmPrev, cmNext, 2)

    guiOnAxisClick(state.onMapClicked!, mapAxes)

    state.redrawMap!()
    state.goTo!(1, 1)

    return state
end

"""
    setAxisRange!(viewer::DispersionViewer, axis::Symbol, lims) -> Nothing

Set the curve map's `:x` or `:y` range from code, as if typed into its min/max fields (no
colour axis here). The fields start blank, so a NaN limit keeps the axis' current one -- see
[`guiSetRangeFromCode!`](@ref).
"""
function GuliwerrJ.setAxisRange!(viewer::DispersionViewer, axis::Symbol, lims)
    owner = "dispersionViewer"
    guiCheckRangeArgs(owner, axis, lims; hasColour = false) || return nothing
    rf, boxes = axis === :x ? (viewer.xRangeField, viewer.xRangeBoxes) : (viewer.yRangeField, viewer.yRangeBoxes)
    guiSetRangeFromCode!(owner, rf, boxes, lims, guiAxisLimits(viewer.mapAxes, axis))
    return nothing
end
