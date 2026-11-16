# =========================================================================================
# GuiSingleAngleViewer.jl -- viewer for single-angle results ("1D curve + click-to-select")
#
# Julia/Makie port of `guliwerr/gui/SingleAngleViewerGUI.m`, `SingleAngleViewerGUIPlate.m` and
# `IsofreqSingleLayeredViewerGUI.m`. The base/Plate pair is a near-duplicate in MATLAB (the
# angle-sweep viewer's deltas one dimension down). `include`d directly into
# `module GuliwerrJMakieExt` (see that file's own header), which already has `using
# GuliwerrJ`/`using Makie` in scope -- neither repeated here. Built on `GuiWidgets.jl`
# (`guiHostLayout`, `guiPlotWithSidePanel`, `guiDropdownStepperRow`, `guiRangeFieldPair`,
# `guiBindRangeToAxis!`, `guiOnAxisClick`, `guiNearestPoint`, `guiCrosshair!`,
# `guiRebuildTypeToggles!`, `guiStatusLabel`/`guiSetStatusOk!`/`guiSetStatusError!`), included
# earlier in `GuliwerrJMakieExt.jl`.
#
# ---------------------------------------------------------------------------------------
# THIS IS A VIEWER, NOT A CALCULATOR
# ---------------------------------------------------------------------------------------
# All three MATLAB constructors take an ALREADY-`calcAll()`'d / `calc()`'d result and only
# ERROR if it looks unsolved (`isempty(singleResult.classification)` /
# `ritzSolver.nCompModes == 0` / `isempty(singleResult.absG)`) -- no solver call anywhere in
# `buildFigure`/the constructor itself. `singleAngleViewer` below does exactly the same: it
# never calls `calcAll`/`calc`/any `Calc.*` entry point. The Ritz spectrum is rebuilt from
# the stored modal numerators, no solve.
#
# ---------------------------------------------------------------------------------------
# DISPATCH DESIGN -- ONE struct, multiple-dispatch `singleAngleViewer` methods
# ---------------------------------------------------------------------------------------
# `SingleAngleViewer` (below) is the ONE handle type returned by all three variants --
# `singleAngleViewer(result::Union{CalcSingleHS,CalcSingleLayered}; fig)`,
# `singleAngleViewer(result::CalcSinglePlate; fig)`,
# `singleAngleViewer(result::CalcIsofreqSingleLayered; fig)` -- rather than three near-copies
# of a bespoke type each. Each method builds the SAME widget shell (host layout, plot +
# side panel, status line, axis-range fields, crosshair, click handler) via ordinary shared
# code (no dispatch needed there -- it is identical across all three), and differs only at
# these few genuine points:
#   - which plotting function draws the curve (`plotClassificationCurve`/
#     `plotRitzModeCurve`/`plotGreenCurve`/`plotRitzSpectrumCurve` vs
#     `plotRitzModeCurvePlate` vs `plotDispersionMap`),
#   - the criterion/graph-type VOCABULARY (three graph types x three criteria for the base
#     variant; one graph type x two criteria for Plate; one graph x two scale modes for
#     isoFreq),
#   - the default x-range SOURCE (`greenVelAxis` falling back to `ritz.vels` vs `ritz.vels`
#     alone vs `kParRange`),
#   - the info-panel TEXT (`guiSignificantModeType` vs `guiPlateModeTypeText` vs a raw
#     kParRange/velocity/|G_33| readout) and the crosshair y-VALUE lookup
#     (`_savValueAtVelocity` methods vs `_savValueAtKPar`).
# Multiple dispatch on `result`'s concrete type selects each of these at its own call site
# (`_savDefaultXRange`, `_savModeInfoLines`, `_savValueAtVelocity`) -- there is no separate
# `SingleAngleViewer`/`SingleAngleViewerPlate`/`IsofreqSingleAngleViewer` type family the way
# MATLAB keeps three separate classdefs.
#
# ---------------------------------------------------------------------------------------
# CLICK-TO-SELECT -- what it resolves against, and why `guiNearestPoint` is used 1D
# ---------------------------------------------------------------------------------------
# MATLAB's `onCurveClicked` in all three source files snaps a click to the nearest sample by
# VELOCITY (or kParRange) ALONE -- `[~, modeNo] = min(abs(vels - velPos))` -- ignoring the
# click's y-coordinate and whichever curve is currently displayed entirely: the click always
# resolves against the FULL computed Ritz mode list (`ritz.vels(1:nCompModes)`) for the base/
# Plate variants, or the full `kParRange` sample list for isoFreq, regardless of graph type.
# `guiNearestPoint` is built for a 2D axis-range-normalised search across several visible/
# hidden SERIES (the map viewer's use case); here there is exactly one series and the
# y-coordinate must never participate. Both are reproduced by calling `guiNearestPoint` with a
# SINGLE series and the y-value forced to `0.0` for both the query and every series point --
# `dy` is then always `0`, so the normalised-distance argmin reduces exactly to
# nearest-by-velocity(or-kParRange), while still routing through the one shared pure function
# (`_savNearestModeByVelocity`/`_savNearestKPar` below), rather than hand-rolling a second,
# parallel nearest-search primitive. `xrange` only needs to be a positive finite constant for
# this to hold (dividing every candidate's `dx` by the same positive constant never changes
# which is smallest) -- the series' own `max - min` (or `eps()` if degenerate) is used,
# matching `guiNearestPoint`'s own normalisation intent without depending on the CURRENTLY
# DISPLAYED axis limits (which change per graph type/criterion and would make an identical
# click resolve to a different mode depending on what happens to be on screen -- MATLAB's own
# click target does not do that either).
#
# ---------------------------------------------------------------------------------------
# DELIBERATE DEVIATIONS FROM THE MATLAB SOURCE (documented, not oversights)
# ---------------------------------------------------------------------------------------
#   1. `<<`/`>>` STEPPER BUTTONS on every dropdown. `guiDropdownStepperRow`
#      always builds a label + dropdown + prev/next button row --
#      SingleAngleViewerGUI's plain `uidropdown` has no such buttons in MATLAB. Wired here to
#      cycle to the previous/next item (wrapping via `mod1`), a small added convenience with no
#      MATLAB counterpart to contradict, not a removed behaviour.
#   2. The base variant's CRITERION dropdown is always visible/enabled, not shown/hidden per
#      graph type the way MATLAB's `obj.criterionLabel.Visible`/`obj.criterionDropdown.Visible`
#      toggle. Makie's `GridLayout` has no per-row "hide this whole row of widgets" primitive
#      the way a MATLAB `uigridlayout` child's `Visible` property provides (`GuiWidgets.jl`
#      exposes no such helper either). The dropdown simply has no effect unless "Ritz-mode
#      curve" is the active graph type, which is discoverable (its own row is right there,
#      always) rather than actively misleading -- a documented simplification, not a missing
#      feature.
#   3. TYPE-TOGGLE SERIES IDENTIFICATION reads back `ax.scene.plots`, filtering for
#      `p isa Makie.Scatter && haskey(p, :label)` and using each match's own `p.label[]` as the
#      type name -- because `plotClassificationCurve`'s return type is `Axis` alone
#      (`PlotStubs.jl`), not MATLAB's one-handle-per-type array. `_diagScatterByType!`
#      (`PlotsDiagnostic.jl`) draws each type's scatter with `label = String(t)`, and Makie
#      genuinely stores an arbitrary `label` kwarg as a live plot attribute
#      (`haskey(plotobj, :label)` is `false` for a scatter/line drawn with no `label`, so the
#      crosshair's own un-labelled scatter is never mistaken for a classification-type series).
#   4. EXPLICIT AXIS-RANGE RE-APPLICATION after every redraw (`_savReapplyRange!`). Once
#      `guiBindRangeToAxis!` sets a NON-`nothing` `ax.limits` entry, Makie keeps honouring it
#      across `empty!(ax)` + a fresh plot call -- UNLESS the newly-called
#      plotting function itself calls `xlims!`/`ylims!` again with its OWN default
#      (`plotDispersionMap`'s `climRange` default, applied whenever `scaleMode === :log`),
#      which would silently override a user-typed range on the next redraw. Reapplying the
#      current (possibly user-edited) `GuiRangeField` value at the end of every redraw mirrors
#      MATLAB's own explicit `obj.applyAxisRange()` call at the end of
#      `onGraphTypeChanged`/`onScaleModeChanged` -- "the fields win", every redraw, exactly as
#      in the MATLAB source, not a new behaviour.
#   5. `_savInterp1Extrap` duplicates `PlotsDiagnostic.jl`'s private `_diagInterp1` (scalar
#      query, not vector) rather than calling it directly. `GuiWidgets.jl`'s own header
#      explains the convention this follows: only `gui`-prefixed names are meant to be reused
#      ACROSS files sharing this module; every `_xxx`-prefixed helper (`_psc_`, `_maps_`,
#      `_diag_`, and this file's own `_sav_`) is private to the file that defines it, so that
#      the files never depend on each other's internals surviving a later refactor.
#   6. `guiSignificantModeType`/`guiPlateModeTypeText` (below) ARE `gui`-prefixed and meant to
#      be reused by the other viewer files that build a mode-info panel: they are
#      physics-result formatters, not generic widget mechanics (see `GuiWidgets.jl`'s own
#      header), so they live here rather than there.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# Shared, reusable result-reading helpers (`gui`-prefixed -- see deviation 6 above)
# -----------------------------------------------------------------------------------------

"""
    guiSignificantModeType(significantModes::AbstractVector{SignificantMode}, modeNo::Integer)
        -> Symbol

Classification type of Ritz mode `modeNo`, if it appears in `significantModes` -- port of
`guliwerr/gui/GuiCommon.m`'s `significantModeType`. Returns `:unclassified` (not an error) when
`modeNo` names a genuinely computed Ritz mode that simply isn't a significant one -- every
viewer lets the user click/step onto ANY computed mode, most of which are
legitimately not significant.
"""
function guiSignificantModeType(significantModes::AbstractVector{SignificantMode}, modeNo::Integer)
    idx = findfirst(sm -> sm.modeNo == modeNo, significantModes)
    return idx === nothing ? :unclassified : significantModes[idx].type
end

"""
    guiPlateModeTypeText(ritz::RitzSolverLamb, modeNo::Integer) -> String

`"<type name> (confidence <x.xx>)"` for one Lamb/plate mode -- port of `GuiCommon.m`'s
`plateModeTypeText`. Names come from [`plateModeTypeStyle`](@ref), the one source for the
plate mode-type visual vocabulary, so a viewer's own label always agrees with the plot it sits
next to.
"""
function guiPlateModeTypeText(ritz::RitzSolverLamb, modeNo::Integer)
    typeNames, _ = plateModeTypeStyle()
    conf = round(ritz.modeTypeConfidence[modeNo]; digits = 2)
    return "$(typeNames[ritz.modesType[modeNo]]) (confidence $conf)"
end

# -----------------------------------------------------------------------------------------
# Private pure helpers (`_sav`-prefixed -- private to this file, see deviation 5 above)
# -----------------------------------------------------------------------------------------

"""
    _savInterp1Extrap(x, y, xq) -> Float64

Linear interpolation of `y(x)` at scalar `xq`, linearly extrapolated beyond `x`'s own range
(MATLAB's `interp1(x, y, xq, 'linear', 'extrap')`). `x` assumed sorted ascending (every
`greenVelAxis`/`kParRange` in this codebase is built that way).
"""
function _savInterp1Extrap(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::Real)
    n = length(x)
    i = if xq <= x[1]
        1
    elseif xq >= x[end]
        n - 1
    else
        clamp(searchsortedlast(x, xq), 1, n - 1)
    end
    t = (xq - x[i]) / (x[i + 1] - x[i])
    return y[i] + t * (y[i + 1] - y[i])
end

"""
    _savNearestModeByVelocity(vels::AbstractVector{<:Real}, velQuery::Real) -> Union{Int,Nothing}

Index into `vels` nearest `velQuery`, by velocity alone -- the pure core of every variant's
`onCurveClicked`. Routes through [`guiNearestPoint`](@ref) with a single series and `y` forced
to `0` on both sides (see this file's header for why that reduces to a 1D search); `nothing`
only when `vels` is empty.
"""
function _savNearestModeByVelocity(vels::AbstractVector{<:Real}, velQuery::Real)
    isempty(vels) && return nothing
    xr = max(maximum(vels) - minimum(vels), eps())
    zeros_ = zeros(Float64, length(vels))
    hit = guiNearestPoint([vels], [zeros_], (Float64(velQuery), 0.0); xrange = xr, yrange = 1.0)
    return hit === nothing ? nothing : hit.index
end

"""
    _savNearestKPar(kParRange::AbstractVector{<:Real}, kQuery::Real) -> Union{Int,Nothing}

`kParRange` analogue of [`_savNearestModeByVelocity`](@ref), for the isoFreq variant (no Ritz
modes -- the continuous `kParRange` axis is the only thing to click onto).
"""
function _savNearestKPar(kParRange::AbstractVector{<:Real}, kQuery::Real)
    isempty(kParRange) && return nothing
    xr = max(maximum(kParRange) - minimum(kParRange), eps())
    zeros_ = zeros(Float64, length(kParRange))
    hit = guiNearestPoint([kParRange], [zeros_], (Float64(kQuery), 0.0); xrange = xr, yrange = 1.0)
    return hit === nothing ? nothing : hit.index
end

"""
    _savReapplyRange!(ax::Axis, rf::GuiRangeField, kind::Symbol) -> Nothing

Re-apply `rf`'s CURRENT value to `ax` right now, unconditionally (unlike
[`guiBindRangeToAxis!`](@ref), which only reacts to future changes) -- see deviation 4 in this
file's header for why every variant's redraw calls this after drawing.
"""
function _savReapplyRange!(ax, rf::GuiRangeField, kind::Symbol)
    r = rf.range[]
    r === nothing && return nothing
    kind === :x ? xlims!(ax, r) : ylims!(ax, r)
    return nothing
end

"""
    _savWireStepper!(menu::Menu, prevBtn::Button, nextBtn::Button, n::Integer) -> Nothing

Wire a dropdown-stepper row's `<<`/`>>` buttons to cycle `menu.i_selected` through `1:n`,
wrapping (see deviation 1 in this file's header).
"""
function _savWireStepper!(menu, prevBtn, nextBtn, n::Integer)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] - 1, n)), prevBtn.clicks)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] + 1, n)), nextBtn.clicks)
    return nothing
end

# -- default x-range source, dispatched per variant --

"""
    _savDefaultXRange(result::Union{CalcSingleHS,CalcSingleLayered}) -> (String, String)

Default velocity-axis text-field values: `greenVelAxis`'s own min/max under `calcMode = :full`,
falling back to the Ritz spectrum's own velocity range under `:ritzOnly` (no `greenVelAxis` at
all then) -- port of `SingleAngleViewerGUI.buildFigure`'s own fallback logic.
"""
function _savDefaultXRange(result::Union{CalcSingleHS,CalcSingleLayered})
    if !isempty(result.greenVelAxis)
        return string(minimum(result.greenVelAxis)), string(maximum(result.greenVelAxis))
    end
    ritz = result.ritzSolver
    vels = ritz.vels[1:ritz.nCompModes]
    return string(minimum(vels)), string(maximum(vels))
end

"""
    _savDefaultXRange(result::CalcSinglePlate) -> (String, String)

Default velocity-axis text-field values: always the Ritz spectrum's own min/max (no Green's
function exists for a free plate) -- port of `SingleAngleViewerGUIPlate.buildFigure`.
"""
function _savDefaultXRange(result::CalcSinglePlate)
    ritz = result.ritzSolver
    vels = ritz.vels[1:ritz.nCompModes]
    return string(minimum(vels)), string(maximum(vels))
end

# -- info-panel text, dispatched per variant --

"""
    _savModeInfoLines(result::Union{CalcSingleHS,CalcSingleLayered}, modeNo::Integer,
                       velSel::Real) -> Vector{String}

Four-line info-panel text for one selected Ritz mode -- port of `SingleAngleViewerGUI.
onCurveClicked`'s `obj.infoPanel.Text` assignment. `type` comes from
[`guiSignificantModeType`](@ref) against `result.classification.significantModes`.
"""
function _savModeInfoLines(result::Union{CalcSingleHS,CalcSingleLayered}, modeNo::Integer,
                            velSel::Real)
    typeStr = result.classification === nothing ? "unclassified" :
              string(guiSignificantModeType(result.classification.significantModes, modeNo))
    return ["angle = $(round(result.angle, digits = 1)) deg",
            "mode # = $modeNo",
            "velocity = $(round(velSel, digits = 4)) km/s",
            "type = $typeStr"]
end

"""
    _savModeInfoLines(result::CalcSinglePlate, modeNo::Integer, velSel::Real) -> Vector{String}

Plate counterpart of the method above -- port of `SingleAngleViewerGUIPlate.onCurveClicked`.
`type` comes from [`guiPlateModeTypeText`](@ref) (mode type + confidence), not a
`classification` lookup -- a `CalcSinglePlate` has no significant-vs-insignificant distinction
to classify into (every Ritz mode already IS a genuine physical branch).
"""
function _savModeInfoLines(result::CalcSinglePlate, modeNo::Integer, velSel::Real)
    return ["angle = $(round(result.angle, digits = 1)) deg",
            "mode # = $modeNo",
            "velocity = $(round(velSel, digits = 4)) km/s",
            "type = $(guiPlateModeTypeText(result.ritzSolver, modeNo))"]
end

"""
    _savPointInfoLines(result::CalcIsofreqSingleLayered, idx::Integer, kParVal::Real)
        -> Vector{String}

Four-line info-panel text for one selected `kParRange` sample -- port of
`IsofreqSingleLayeredViewerGUI.onCurveClicked`'s `obj.infoPanel.Text` assignment. `velocity`
is `omega0 / kParVal` (this class never runs a Ritz solve; velocity is derived directly from
the fixed reference frequency), matching the MATLAB source exactly.
"""
function _savPointInfoLines(result::CalcIsofreqSingleLayered, idx::Integer, kParVal::Real)
    velVal = result.omega0 / kParVal
    return ["angle = $(round(result.angle, digits = 1)) deg",
            "kParRange = $(round(kParVal, sigdigits = 4))",
            "velocity = $(round(velVal, digits = 4)) km/s",
            "|G_{33}| = $(round(result.absG[idx], sigdigits = 4))"]
end

# -- crosshair y-value lookup, dispatched per variant/graph-type --

"""
    _savValueAtVelocity(result::Union{CalcSingleHS,CalcSingleLayered}, graphType::Symbol,
                         criterion::Symbol, velVal::Real) -> Float64

Y-value the CURRENTLY DISPLAYED curve would show at `velVal` -- port of
`SingleAngleViewerGUI.valueAtVelocity`. `graphType` -- `:classification` |
`:green` | `:ritzMode` | `:ritzSpectrum`; `criterion` (`:ritzScore`/`:freqDisp`/
`:depthInsensitivity`) only for `:ritzMode`.

- `:ritzMode`: the nearest mode's criterion value ([`visRitzCriterion`](@ref); 0 where `NaN`).
- `:ritzSpectrum`, or `:classification` on a `:ritzOnly` result: linear interpolation onto the
  log-normalised Ritz spectrum.
- `:green` / `:classification`: linear interpolation (extrapolated) onto the log-normalised
  `|G13_eff|` curve ([`scaleGreenData`](@ref)); `0.0` without Green data.
"""
function _savValueAtVelocity(result::Union{CalcSingleHS,CalcSingleLayered}, graphType::Symbol,
                              criterion::Symbol, velVal::Real)
    if graphType === :ritzMode
        ritz = result.ritzSolver
        vels = ritz.vels[1:ritz.nCompModes]
        critVal, _ = visRitzCriterion(result, criterion)
        v = critVal[argmin(abs.(vels .- velVal))]
        return isnan(v) ? 0.0 : v
    end
    useSpectrum = graphType === :ritzSpectrum ||
                  (graphType === :classification && isempty(result.greenVelAxis))
    if useSpectrum
        visHasRitzSpectrum(result) || return 0.0
        absS, v = visRitzSpectrumData(result)
        return _savInterp1Extrap(v, scaleGreenData(absS, :log, "S")[1], velVal)
    end
    isempty(result.greenAbsG13Eff) && return 0.0
    yData, _ = scaleGreenData(result.greenAbsG13Eff, :log, "G_{13}^{eff}")
    return _savInterp1Extrap(result.greenVelAxis, yData, velVal)
end

"""
    _savValueAtVelocity(result::CalcSinglePlate, criterion::Symbol, velVal::Real) -> Float64

Plate counterpart: nearest mode's own `modesType`/`modeTypeConfidence` value -- port of
`SingleAngleViewerGUIPlate.valueAtVelocity`. `criterion` -- `:modeType` | `:confidence`.
"""
function _savValueAtVelocity(result::CalcSinglePlate, criterion::Symbol, velVal::Real)
    ritz = result.ritzSolver
    nM = ritz.nCompModes
    vels = ritz.vels[1:nM]
    idx = argmin(abs.(vels .- velVal))
    return criterion === :modeType ? Float64(ritz.modesType[idx]) : ritz.modeTypeConfidence[idx]
end

"""
    _savValueAtKPar(result::CalcIsofreqSingleLayered, scaleMode::Symbol, kParVal::Real) -> Float64

Y-value the isoFreq curve shows at `kParVal` -- linear interpolation (extrapolated) onto the
same (log-normalised or raw) `|G_33|` curve `plotDispersionMap` itself draws for this result
type, via [`scaleGreenData`](@ref). Port of `IsofreqSingleLayeredViewerGUI.valueAtKPar`.
"""
function _savValueAtKPar(result::CalcIsofreqSingleLayered, scaleMode::Symbol, kParVal::Real)
    yData, _ = scaleGreenData(result.absG, scaleMode, "G_{33}")
    return _savInterp1Extrap(result.kParRange, yData, kParVal)
end

# -----------------------------------------------------------------------------------------
# SingleAngleViewer -- the one handle type returned by every variant
# -----------------------------------------------------------------------------------------

"""
    SingleAngleViewer

Handles for one "1D curve + click-to-select" viewer window, returned by every
[`singleAngleViewer`](@ref) method regardless of which `CalcSingleHS`/`CalcSingleLayered`/
`CalcSinglePlate`/`CalcIsofreqSingleLayered` variant built it (see this file's header for the
dispatch design).

Fields:
- `fig`, `ax` -- the owning `Figure` and the one curve `Axis`.
- `result` -- the `Calc*` result this viewer displays (untyped field; concretely one of the
  four supported types -- kept generic so this ONE struct serves every variant).
- `graphTypeMenu` -- the graph-type dropdown (`Menu`) for the base (HS/Layered) variant, or
  `nothing` for Plate/isoFreq, which have exactly one graph type each (absent, not disabled,
  matching MATLAB's own precedent for those two classes).
- `modeMenu` -- the variant's OTHER dropdown: the Ritz-mode criterion selector
  (`ritzScore`/`freqDisp`/`depthInsensitivity`) for the base variant, the criterion selector
  (`Mode type`/`Confidence`) for Plate, or the scale-mode selector (`log`/`linear`) for
  isoFreq. Never `nothing` -- every variant has exactly one of these.
- `xRangeField`/`yRangeField` -- persistent axis-range state (`GuiRangeField`); `xRangeBoxes`/
  `yRangeBoxes` -- the mounted `(loBox, hiBox)` `Textbox` pairs.
- `statusLabel` -- ok/error status line (`guiStatusLabel`/`guiSetStatusOk!`/
  `guiSetStatusError!`).
- `infoLabel` -- the mode/point info panel, updated on each click.
- `typeTogglePanel` -- the "Show types:" per-type visibility panel (`GridLayout`), populated
  only when the base variant's Classification curve is displayed (see deviation 3 in this
  file's header); `nothing` for Plate/isoFreq, which have no per-type legend to filter.
- `crosshair` -- the current crosshair position `Observable{Point2f}` (`guiCrosshair!`);
  reassigned on every redraw (the previous one is destroyed by `empty!(ax)`, exactly mirroring
  MATLAB's own `cla` + recreate-on-demand `hCrosshair` handling).
- `onSelect` -- caller-settable callback (`Ref{Union{Nothing,Function}}`, `nothing` initially),
  MATLAB's `OnModeSelected`/`OnPointSelected` property. Called after every successful click:
  `(angle, modeNo, velSnapped, velRaw)` for the base/Plate variants, `(angle, kParVal)` for
  isoFreq -- set it with `viewer.onSelect[] = (args...) -> ...`.
- `redraw!` -- `() -> Nothing`, redraws the curve for whichever dropdown selection is current;
  also callable directly (e.g. from a test) without going through a real widget event.
- `onClicked!` -- `(x::Float64, y::Float64) -> Nothing`, the click handler `guiOnAxisClick` is
  wired to; also callable directly with a known data-space point, in place of a simulated
  mouse event (e.g. under headless CairoMakie).
"""
mutable struct SingleAngleViewer
    fig::Figure
    ax::Axis
    result::Any
    graphTypeMenu::Union{Menu,Nothing}
    modeMenu::Menu
    xRangeField::GuiRangeField
    yRangeField::GuiRangeField
    xRangeBoxes::Tuple{Textbox,Textbox}
    yRangeBoxes::Tuple{Textbox,Textbox}
    statusLabel::Label
    infoLabel::Label
    typeTogglePanel::Union{GridLayout,Nothing}
    crosshair::Observable{Point2f}
    onSelect::Base.RefValue{Union{Nothing,Function}}
    redraw!::Function
    onClicked!::Function
end

# -----------------------------------------------------------------------------------------
# singleAngleViewer -- base variant (CalcSingleHS / CalcSingleLayered)
# -----------------------------------------------------------------------------------------

"""
    singleAngleViewer(result::Union{CalcSingleHS,CalcSingleLayered}; fig = nothing)
        -> SingleAngleViewer

Single-angle curve viewer for an already-`calcAll()`'d `CalcSingleHS`/`CalcSingleLayered`
result -- port of `guliwerr/gui/SingleAngleViewerGUI.m`. `fig` -- `nothing` (default) for a
standalone window, or anything [`guiHostLayout`](@ref) accepts as `parent` (e.g. a
`GridPosition` such as `otherFig[2, 1]`) to embed this viewer's content into an existing
figure -- forwarded there directly.

Four graph types in MATLAB's order (`graphTypeMenu`, default "Classification curve"):
"Classification curve" ([`plotClassificationCurve`](@ref): G13_eff background, or the Ritz
spectrum for a `:ritzOnly` result), "Green's-function curve" ([`plotGreenCurve`](@ref),
G13_eff), "Ritz-mode curve" ([`plotRitzModeCurve`](@ref)), "Ritz spectrum"
([`plotRitzSpectrumCurve`](@ref)). `modeMenu` is the Ritz-mode criterion dropdown
(`ritzScore`/`freqDisp`/`depthInsensitivity`), consulted only while "Ritz-mode curve" is
selected (see deviation 2 in this file's header).

Throws `ArgumentError` if `result.classification === nothing` (not yet `calcAll()`'d).
"""
function GuliwerrJ.singleAngleViewer(result::Union{CalcSingleHS,CalcSingleLayered}; fig = nothing)
    result.classification === nothing && throw(ArgumentError(
        "singleAngleViewer: result has not been calcAll()'d yet (classification is nothing) " *
        "-- nothing to view."))

    fig, mainLayout = guiHostLayout(fig; size = (1000, 700))
    topBarGL = guiGridLayout(mainLayout[1, 1])
    ax, sideGL, _ = guiPlotWithSidePanel(mainLayout[2, 1])
    statusLabel = guiStatusLabel(mainLayout[3, 1])

    graphItems = ["Classification curve", "Green's-function curve", "Ritz-mode curve",
                  "Ritz spectrum"]
    graphTypes = (:classification, :green, :ritzMode, :ritzSpectrum)
    critSymbols = visRitzCriterionNames()
    critItems = [String(c) for c in critSymbols]

    graphMenu, graphPrev, graphNext = guiDropdownStepperRow(topBarGL, 1, "Graph type:", graphItems)
    criterionMenu, critPrev, critNext = guiDropdownStepperRow(topBarGL, 2, "Criterion:", critItems)

    xLoDefault, xHiDefault = _savDefaultXRange(result)
    xRangeField, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "Velocity [km/s] min/max";
                                                     loInit = xLoDefault, hiInit = xHiDefault)
    # Value (y) range left blank/auto -- the y-axis meaning changes per graph type/criterion,
    # so no single default range covers all three (SingleAngleViewerGUI.buildFigure's own
    # reasoning).
    yRangeField, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "Value [Y] min/max")
    guiBindRangeToAxis!(ax, xRangeField, :x)
    guiBindRangeToAxis!(ax, yRangeField, :y)

    typeTogglePanel = guiGridLayout(sideGL[5, 1])
    infoLabel = Label(sideGL[6, 1], "Click the curve to select a mode.";
                       valign = :top, halign = :left)

    state = SingleAngleViewer(fig, ax, result, graphMenu, criterionMenu,
                               xRangeField, yRangeField, (xLoBox, xHiBox), (yLoBox, yHiBox),
                               statusLabel, infoLabel, typeTogglePanel,
                               Observable(Point2f(NaN, NaN)),
                               Ref{Union{Nothing,Function}}(nothing),
                               () -> nothing, (x, y) -> nothing)

    function _redraw!()
        guiEmptyAxis!(ax)
        graphType = graphTypes[state.graphTypeMenu.i_selected[]]
        criterion = critSymbols[state.modeMenu.i_selected[]]
        ok, errMsg = true, ""
        try
            if graphType === :ritzSpectrum
                plotRitzSpectrumCurve(result, :log; ax = ax)
            elseif graphType === :ritzMode
                plotRitzModeCurve(result, criterion; ax = ax)
            elseif graphType === :green
                plotGreenCurve(result, :log; ax = ax)
            else # :classification
                plotClassificationCurve(result; ax = ax)
            end
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end

        # "Show types:" toggle panel only for the classification curve's per-type series --
        # see deviation 3 in this file's header for why `ax.scene.plots` is walked instead of
        # a returned handle list.
        if ok && graphType === :classification
            series = [(string(p.label[]), p) for p in ax.scene.plots
                      if p isa Makie.Scatter && haskey(p, :label)]
            guiRebuildTypeToggles!(state.typeTogglePanel, series)
        else
            guiRebuildTypeToggles!(state.typeTogglePanel, Tuple{String,Any}[])
        end

        state.crosshair = guiCrosshair!(ax)
        _savReapplyRange!(ax, state.xRangeField, :x)
        _savReapplyRange!(ax, state.yRangeField, :y)

        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _onClicked!(x::Float64, y::Float64)
        ritz = result.ritzSolver
        vels = ritz.vels[1:ritz.nCompModes]
        idx = _savNearestModeByVelocity(vels, x)
        idx === nothing && return nothing
        velSel = vels[idx]

        state.infoLabel.text[] = join(_savModeInfoLines(result, idx, velSel), "\n")
        graphType = graphTypes[state.graphTypeMenu.i_selected[]]
        criterion = critSymbols[state.modeMenu.i_selected[]]
        yPos = _savValueAtVelocity(result, graphType, criterion, velSel)
        state.crosshair[] = Point2f(velSel, yPos)

        cb = state.onSelect[]
        cb === nothing || cb(result.angle, idx, velSel, x)
        return nothing
    end
    state.onClicked! = _onClicked!

    on(_ -> state.redraw!(), graphMenu.i_selected)
    on(_ -> state.redraw!(), criterionMenu.i_selected)
    _savWireStepper!(graphMenu, graphPrev, graphNext, length(graphItems))
    _savWireStepper!(criterionMenu, critPrev, critNext, length(critItems))
    guiOnAxisClick(state.onClicked!, ax)

    state.redraw!()
    return state
end

# -----------------------------------------------------------------------------------------
# singleAngleViewer -- Plate variant (CalcSinglePlate)
# -----------------------------------------------------------------------------------------

"""
    singleAngleViewer(result::CalcSinglePlate; fig = nothing) -> SingleAngleViewer

Single-angle curve viewer for an already-`calcAll()`'d `CalcSinglePlate` result -- port of
`guliwerr/gui/SingleAngleViewerGUIPlate.m`. No graph-type dropdown at all -- ABSENT, not
disabled (`graphTypeMenu` is `nothing`): there is exactly one graph type here
([`plotRitzModeCurvePlate`](@ref)), no Green's-function curve (no analogue for a free plate)
and no classification curve (every Ritz mode already IS a genuine physical branch, nothing to
classify into). `modeMenu` offers "Mode type" | "Confidence"
(`plotRitzModeCurvePlate`'s own `:modeType`/`:confidence`).

Throws `ArgumentError` if `result.ritzSolver.nCompModes == 0` (not yet `calcAll()`'d).
"""
function GuliwerrJ.singleAngleViewer(result::CalcSinglePlate; fig = nothing)
    result.ritzSolver.nCompModes == 0 && throw(ArgumentError(
        "singleAngleViewer: result has not been calcAll()'d yet (ritzSolver.nCompModes == 0) " *
        "-- nothing to view."))

    fig, mainLayout = guiHostLayout(fig; size = (1000, 700))
    topBarGL = guiGridLayout(mainLayout[1, 1])
    ax, sideGL, _ = guiPlotWithSidePanel(mainLayout[2, 1])
    statusLabel = guiStatusLabel(mainLayout[3, 1])

    critItems = ["Mode type", "Confidence"]
    critSymbols = (:modeType, :confidence)
    criterionMenu, critPrev, critNext = guiDropdownStepperRow(topBarGL, 1, "Criterion:", critItems)

    xLoDefault, xHiDefault = _savDefaultXRange(result)
    xRangeField, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "Velocity [km/s] min/max";
                                                     loInit = xLoDefault, hiInit = xHiDefault)
    yRangeField, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "Value [Y] min/max")
    guiBindRangeToAxis!(ax, xRangeField, :x)
    guiBindRangeToAxis!(ax, yRangeField, :y)

    infoLabel = Label(sideGL[5, 1], "Click the curve to select a mode.";
                       valign = :top, halign = :left)

    state = SingleAngleViewer(fig, ax, result, nothing, criterionMenu,
                               xRangeField, yRangeField, (xLoBox, xHiBox), (yLoBox, yHiBox),
                               statusLabel, infoLabel, nothing,
                               Observable(Point2f(NaN, NaN)),
                               Ref{Union{Nothing,Function}}(nothing),
                               () -> nothing, (x, y) -> nothing)

    function _redraw!()
        guiEmptyAxis!(ax)
        criterion = critSymbols[state.modeMenu.i_selected[]]
        ok, errMsg = true, ""
        try
            plotRitzModeCurvePlate(result, criterion; ax = ax)
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end

        state.crosshair = guiCrosshair!(ax)
        _savReapplyRange!(ax, state.xRangeField, :x)
        _savReapplyRange!(ax, state.yRangeField, :y)

        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _onClicked!(x::Float64, y::Float64)
        ritz = result.ritzSolver
        vels = ritz.vels[1:ritz.nCompModes]
        idx = _savNearestModeByVelocity(vels, x)
        idx === nothing && return nothing
        velSel = vels[idx]

        state.infoLabel.text[] = join(_savModeInfoLines(result, idx, velSel), "\n")
        criterion = critSymbols[state.modeMenu.i_selected[]]
        yPos = _savValueAtVelocity(result, criterion, velSel)
        state.crosshair[] = Point2f(velSel, yPos)

        cb = state.onSelect[]
        cb === nothing || cb(result.angle, idx, velSel, x)
        return nothing
    end
    state.onClicked! = _onClicked!

    on(_ -> state.redraw!(), criterionMenu.i_selected)
    _savWireStepper!(criterionMenu, critPrev, critNext, length(critItems))
    guiOnAxisClick(state.onClicked!, ax)

    state.redraw!()
    return state
end

# -----------------------------------------------------------------------------------------
# singleAngleViewer -- isoFreq variant (CalcIsofreqSingleLayered)
# -----------------------------------------------------------------------------------------

"""
    singleAngleViewer(result::CalcIsofreqSingleLayered; fig = nothing) -> SingleAngleViewer

Single-angle, single-`omega0` isoFreq curve viewer -- port of
`guliwerr/gui/IsofreqSingleLayeredViewerGUI.m`. Exactly one possible graph -- the Green
`|G_33|` field vs. `kParRange` ([`plotDispersionMap`](@ref)'s `CalcIsofreqSingleLayered`
method) -- since this class never runs a Ritz solve at all (fixed `omega0`, k-swept
Green's-function scan only). No `graphTypeMenu`/criterion dropdown at all (`nothing`);
`modeMenu` is the scale-mode dropdown ("log" | "linear"), which re-defaults the Value range
fields to `-5`/`-2` (log) or blank (linear) on every change, matching
`onScaleModeChanged`'s own re-defaulting.

`onSelect`, if set, is called `(angle, kParVal)` after a click -- NOT the base/Plate variants'
`(angle, modeNo, velSnapped, velRaw)`: `kParRange` is a continuous, physically-meaningful axis
here, not a discrete Ritz mode index, matching `IsofreqSingleLayeredViewerGUI`'s own
`OnPointSelected` (vs. `OnModeSelected`) naming.

Throws `ArgumentError` if `result.absG` is empty (not yet `calc()`'d).
"""
function GuliwerrJ.singleAngleViewer(result::CalcIsofreqSingleLayered; fig = nothing)
    isempty(result.absG) && throw(ArgumentError(
        "singleAngleViewer: result has not been calc()'d yet (absG is empty) -- nothing to view."))

    fig, mainLayout = guiHostLayout(fig; size = (1000, 700))
    topBarGL = guiGridLayout(mainLayout[1, 1])
    ax, sideGL, _ = guiPlotWithSidePanel(mainLayout[2, 1])
    statusLabel = guiStatusLabel(mainLayout[3, 1])

    scaleItems = ["log", "linear"]
    scaleSymbols = (:log, :linear)
    scaleMenu, scalePrev, scaleNext = guiDropdownStepperRow(topBarGL, 1, "Scale:", scaleItems)

    xRangeField, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "kParRange min/max";
        loInit = string(minimum(result.kParRange)), hiInit = string(maximum(result.kParRange)))
    # Value (y) range defaults to plotDispersionMap's own 1D-line default ([-5,0], top at the normalised peak) since
    # 'log' is the initial scale mode -- re-defaulted on every scale switch below, same
    # "re-default per active mode" mechanic AngSweepViewerGUI's own colour-range fields use.
    yRangeField, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "Value [Y] min/max";
                                                     loInit = "-5", hiInit = "0")
    guiBindRangeToAxis!(ax, xRangeField, :x)
    guiBindRangeToAxis!(ax, yRangeField, :y)

    infoLabel = Label(sideGL[5, 1], "Click the curve to select a point.";
                       valign = :top, halign = :left)

    state = SingleAngleViewer(fig, ax, result, nothing, scaleMenu,
                               xRangeField, yRangeField, (xLoBox, xHiBox), (yLoBox, yHiBox),
                               statusLabel, infoLabel, nothing,
                               Observable(Point2f(NaN, NaN)),
                               Ref{Union{Nothing,Function}}(nothing),
                               () -> nothing, (x, y) -> nothing)

    function _redraw!()
        guiEmptyAxis!(ax)
        scaleMode = scaleSymbols[state.modeMenu.i_selected[]]
        if scaleMode === :log
            yLoBox.stored_string[] = "-5"
            yHiBox.stored_string[] = "0"
        else
            yLoBox.stored_string[] = ""
            yHiBox.stored_string[] = ""
        end

        ok, errMsg = true, ""
        try
            plotDispersionMap(result, 1, scaleMode; ax = ax)
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end

        state.crosshair = guiCrosshair!(ax)
        _savReapplyRange!(ax, state.xRangeField, :x)
        _savReapplyRange!(ax, state.yRangeField, :y)

        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _onClicked!(x::Float64, y::Float64)
        idx = _savNearestKPar(result.kParRange, x)
        idx === nothing && return nothing
        kParVal = result.kParRange[idx]

        state.infoLabel.text[] = join(_savPointInfoLines(result, idx, kParVal), "\n")
        scaleMode = scaleSymbols[state.modeMenu.i_selected[]]
        yPos = _savValueAtKPar(result, scaleMode, kParVal)
        state.crosshair[] = Point2f(kParVal, yPos)

        cb = state.onSelect[]
        cb === nothing || cb(result.angle, kParVal)
        return nothing
    end
    state.onClicked! = _onClicked!

    on(_ -> state.redraw!(), scaleMenu.i_selected)
    _savWireStepper!(scaleMenu, scalePrev, scaleNext, length(scaleItems))
    guiOnAxisClick(state.onClicked!, ax)

    state.redraw!()
    return state
end

"""
    setAxisRange!(viewer::SingleAngleViewer, axis::Symbol, lims) -> Nothing

Set the curve's `:x` or `:y` range from code, as if typed into its min/max fields (no colour
axis here) -- see [`guiSetRangeFromCode!`](@ref).
"""
function GuliwerrJ.setAxisRange!(viewer::SingleAngleViewer, axis::Symbol, lims)
    owner = "singleAngleViewer"
    guiCheckRangeArgs(owner, axis, lims; hasColour = false) || return nothing
    rf, boxes = axis === :x ? (viewer.xRangeField, viewer.xRangeBoxes) : (viewer.yRangeField, viewer.yRangeBoxes)
    guiSetRangeFromCode!(owner, rf, boxes, lims, guiAxisLimits(viewer.ax, axis))
    return nothing
end
