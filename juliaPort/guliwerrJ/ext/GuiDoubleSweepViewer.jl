# =========================================================================================
# GuiDoubleSweepViewer.jl -- viewer for double sweeps ("outer slice-stepper shell over one
# extra sweep axis")
#
# Julia/Makie port of `guliwerr/gui/DoubleSweepViewerGUI.m`, `DoubleSweepViewerGUIPlate.m` and
# `FreqAngSweep2LayeredViewerGUI.m`: an outer slice-stepper shell over one extra sweep axis,
# a no-plotting composite host owning layout and child lifecycle. `include`d directly into the
# EXISTING `module GuliwerrJMakieExt`, after `GuiCombinedViewer.jl` -- this file's whole job is
# to rebuild that file's `combinedViewer`/`CombinedViewer` against a different outer-sweep
# slice on every step.
#
# ---------------------------------------------------------------------------------------
# THE WRAPPED CHILD IS A COMBINED VIEWER, NOT A BARE ANGLE-SWEEP OR SINGLE-ANGLE VIEWER
# ---------------------------------------------------------------------------------------
# `DoubleSweepViewerGUI.buildAngularContent` builds `ModalShapeViewerGUI` + `AngSweepViewerGUI`
# TOGETHER, wired exactly as `CombinedViewerGUI` wires them (its own docblock says so
# explicitly: "wired exactly as CombinedViewerGUI wires them");
# `DoubleSweepViewerGUIPlate.buildAngularContent` does the same with the Plate pair;
# `FreqAngSweep2LayeredViewerGUI.buildContent` does the same with the isoFreq pair ("wired
# exactly as IsofreqCombinedViewerGUI wires them"). So the child this file wraps, for every
# one of its three variants, is a composite (`combinedViewer`/`CombinedViewer`,
# `ext/GuiCombinedViewer.jl`), never a bare angle-sweep or single-angle viewer on its own.
#
# ---------------------------------------------------------------------------------------
# TWO MODES -- `DoubleSweepViewerGUI`/`...Plate` ARE BIMODAL IN MATLAB;
# `FreqAngSweep2LayeredViewerGUI` IS NOT
# ---------------------------------------------------------------------------------------
# Both base and Plate `DoubleSweepViewerGUI` classes toggle between TWO qualitatively different
# whole-pair views via a "Switch to dispersion view ->" button: "angular" mode (the combined
# viewer above, at a fixed wavelength/wavevector slice) and "dispersion" mode
# (`DispersionViewerGUI`/`...Plate` embedded via their own `parentContainer` support, at a fixed
# angle slice -- `ext/GuiDispersionViewer.jl`). `FreqAngSweep2LayeredViewerGUI` has NO second
# mode at all (its own docblock:
# "UNLIKE DoubleSweepViewerGUI, there is no second mode to switch to... every omega0 step
# rebuilds the same one view, always resetting to the map" -- no Ritz, no branch tracking, no
# dispersion-curve concept exists for isoFreq; `FreqAngSweep2LayeredViewerGUI.m` contains no
# `switchModeBtn`/`SwitchMode`). So the "Switch to dispersion view ->" button is wired to REAL
# dispersion mode for the two variants whose MATLAB source has one
# (`CalcWavevecAngSweep2Layered`/`CalcWavevecAngSweep2Plate`); the isoFreq method
# (`CalcFreqAngSweep2Layered`) builds NO switch button at all (`switchModeBtn = nothing`).
#
# Following `DoubleSweepViewerGUI.m`'s (`onSwitchMode`/`buildAngularContent`/
# `buildDispersionContent`) and `DoubleSweepViewerGUIPlate.m`'s identical-shape methods:
#   - ANGULAR -> DISPERSION: if `angSweepViewer` is currently drilled into its own single-angle
#     child, undrill first (`toggleDrillDown!()` once -- `AngSweepViewer.isDrilledIn[]` is
#     always readable, matching MATLAB's own `obj.angularSweepViewer.isDrilledIn` check; no lost
#     angle -- `currentAngNo` is this file's OWN tracked value, updated on every click/step
#     regardless of drill state, see below), then build [`dispersionViewer`](@ref) at
#     `state.currentAngNo`.
#   - DISPERSION -> ANGULAR: read the dispersion viewer's OWN `currentWavNo` back into
#     `state.currentIdx[]` (MATLAB's `obj.currentWavNo = obj.dispersionViewer.currentWavNo`) --
#     the two viewers track the SAME physical wavelength/wavevector independently, and only the
#     wavelength side needs carrying over (branch index has no cross-mode equivalent, exactly
#     MATLAB's own comment: "Branch index has no cross-mode equivalent to restore, so it stays
#     at 1").
#   - `currentAngNo` TRACKING, while in angular mode: MATLAB's `buildAngularContent` wraps
#     `angularSweepViewer.OnModeSelected`/`angularModalShapeViewer.OnSelectionChanged` in ITS OWN
#     handlers (`onAngularModeSelected`/`onAngularSelectionChanged`) that update
#     `obj.currentAngNo` THEN forward to the SAME calls `CombinedViewerGUI` would make anyway.
#     Since this port's angular mode reuses [`combinedViewer`](@ref)
#     WHOLESALE -- which already performs that exact forwarding internally via
#     `_cvWireModeFamily!` -- reproducing MATLAB's own extra wrapping layer here means CHAINING
#     onto whatever callback `combinedViewer`'s own construction already installed, not
#     replacing it: this file reads `cv.angSweepViewer.onSelect[]`/
#     `cv.modalShapeViewer.onSelectionChanged[]` right after `combinedViewer(...)` returns, then
#     overwrites each with a new closure that updates `state.currentAngNo` and THEN calls the
#     saved original -- so the combined viewer's own bidirectional modal-shape<->map link keeps
#     working exactly as it does standalone, with angle-tracking layered on top, not duplicated.
#   - AXIS-RANGE SAVE/RESTORE ACROSS A REBUILD (MATLAB's `GuiCommon.captureAxisRange`/
#     `restoreAxisRange`, `obj.savedAxisRangeAngular`/`savedAxisRangeDispersion`) is reproduced
#     by `_dswCaptureRange`/`_dswRestoreRange!`: the min/max field TEXT of the outgoing map (or
#     dispersion) viewer is kept in `savedRangeAngular`/`savedRangeDispersion` and written into
#     the next viewer of the SAME kind via `Makie.set!` on its mounted boxes, so the new
#     `GuiRangeField`s and axes follow through their own listeners. As in MATLAB, each mode keeps
#     its own snapshot -- nothing is carried between the two kinds. Taken on every slice step
#     (MATLAB's `onWavChanged`/`onOmegaChanged`, skipped while drilled in) and on every mode
#     switch.
#
# ---------------------------------------------------------------------------------------
# `currentAngNo` TRACKING NEEDS NO NEW STATE -- ALREADY EXPOSED BY THE EMBEDDED CombinedViewer
# ---------------------------------------------------------------------------------------
# MATLAB's `onAngularModeSelected`/`onAngularSelectionChanged` track `obj.currentAngNo` purely
# so a LATER switch to dispersion mode lands on the right angle. Every outer-sweep type this
# file wraps (`CalcWavevecAngSweep2Layered.sweeps::Vector{CalcAngSweepLayered}`,
# `CalcWavevecAngSweep2Plate.sweeps::Vector{CalcAngSweepPlate}`,
# `CalcFreqAngSweep2Layered.sweeps::Vector{CalcIsofreqAngSweepLayered}`) holds ONLY genuine
# multi-angle sweep elements (never a single-point result) -- per each struct's
# own field type -- so the embedded `CombinedViewer`'s own `.angSweepViewer` field is ALWAYS a
# real `AngSweepViewer` (never a `SingleAngleViewer`), which already exposes its own
# `currentAngNo::Base.RefValue{Int}`, readable as
# `doubleSweepViewer.combinedViewer[].angSweepViewer.currentAngNo[]`; this file adds no
# redundant tracking field of its own.
#
# ---------------------------------------------------------------------------------------
# GENUINE-DOUBLE-SWEEP GUARD
# ---------------------------------------------------------------------------------------
# `DoubleSweepViewerGUI`/`...Plate`'s own constructors reject `numel(angles) < 2 ||
# numel(<outer axis>) < 2` with a message pointing the caller at `CombinedViewerGUI`
# (angle-only) or `DispersionViewerGUI` (outer-axis-only) directly. This port keeps the SAME
# guard against the angle axis (steering toward [`combinedViewer`](@ref)) -- a useful "you
# probably want a simpler viewer" check.
# `FreqAngSweep2LayeredViewerGUI` never had this guard in MATLAB either (no second mode to steer
# away from), so its Julia method has none.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# DoubleSweepViewer -- the one handle type returned by every variant
# -----------------------------------------------------------------------------------------

"""
    DoubleSweepViewer

Handles for one "outer slice-stepper shell" viewer, returned by every
[`doubleSweepViewer`](@ref) method. Owns the outer stepper (wavelength/wavevector/omega0) and a
stable content slot into which either a fresh [`combinedViewer`](@ref) ("angular"
mode) or a fresh [`dispersionViewer`](@ref) ("dispersion" mode -- Layered/Plate
only) is rebuilt on every step/mode switch -- mirrors `DoubleSweepViewerGUI`'s own
`contentParent` + `buildAngularContent`/`buildDispersionContent` pattern (see this file's
header for the isoFreq exception, which has no dispersion mode at all).

Fields:
- `fig::Figure` -- the (possibly ancestor) `Figure`.
- `sweepResult` -- the outer 2D-sweep result (`CalcWavevecAngSweep2Layered`/
  `CalcWavevecAngSweep2Plate`/`CalcFreqAngSweep2Layered`; untyped so this one struct serves all
  three, matching every other viewer file's own "kept generic" precedent).
- `contentGL::GridLayout` -- the stable content slot the embedded child is rebuilt into on
  every step/mode switch ([`guiClearGridLayout!`](@ref), reused verbatim from
  `ext/GuiAngSweepViewer.jl`).
- `axisMenu::Menu` -- the outer-axis dropdown-stepper (wavelength/wavevector/omega0); stays
  mounted (and visible) across a mode switch, but its `on(...)` listener is a no-op while
  `mode !== :angular` -- Makie's `Menu` exposes no `enabled` attribute to disable it outright
  (the same limitation every other viewer file's header records).
- `switchModeBtn::Union{Button,Nothing}` -- the "Switch to dispersion view ->"/"<< Switch to
  angular view" toggle button for the Layered/Plate methods; `nothing` for the isoFreq method,
  which has no second mode at all (ABSENT, not disabled -- matching every other viewer's own
  precedent for a variant with no such concept).
- `statusLabel::Label` -- ok + "current slice"/"current angle" status line.
- `currentIdx::Base.RefValue{Int}` -- index into the outer sweep axis; tracked independently of
  `currentAngNo` and shared by both modes.
- `currentAngNo::Int` -- index into `sweepResult.angles`; tracked independently of `currentIdx`
  and shared by both modes (MATLAB's own `obj.currentAngNo`) -- only ever updated while in
  angular mode (a click/step on the embedded `CombinedViewer`) or read back from the dispersion
  viewer's own selection; unused (stays at its initial value `1`) for the isoFreq method, which
  never builds dispersion mode.
- `mode::Symbol` -- `:angular` (default, every variant) or `:dispersion` (Layered/Plate only,
  after a switch).
- `combinedViewer::Base.RefValue{Union{Nothing,CombinedViewer}}` -- the currently-embedded
  combined-viewer child while `mode === :angular`; `nothing` while `mode === :dispersion`.
- `dispersionViewer::Base.RefValue{Union{Nothing,DispersionViewer}}` -- the currently-embedded
  dispersion-viewer child while `mode === :dispersion`; `nothing` otherwise (always `nothing` for the
  isoFreq method).
- `updatingWidgets::Base.RefValue{Bool}` -- reentrancy guard (the same pattern
  `ext/GuiModalShapeViewer.jl`'s/`ext/GuiDispersionViewer.jl`'s own headers document as
  reusable) -- guards `axisMenu`'s own listener against
  the one genuine re-entrant case here: syncing it back to `currentIdx[]` when returning FROM
  dispersion mode must not itself trigger a second, redundant angular rebuild.
- `rebuild!::Function` -- `() -> Nothing`, rebuilds the embedded `CombinedViewer` for
  `currentIdx[]` (ALWAYS the angular-mode rebuild, regardless of the CURRENT `mode` -- matching
  this field's own meaning); also callable directly (e.g. from a
  test), matching every other viewer file's own `redraw!`/`rebuild!` convention.
- `onSwitchMode!::Union{Function,Nothing}` -- `() -> Nothing`, MATLAB's `onSwitchMode` --
  toggles between angular and dispersion mode, carrying over `currentIdx[]`/`currentAngNo`
  exactly as MATLAB's own version does and saving/restoring each mode's axis ranges;
  `nothing` for the isoFreq method, which has no second mode to switch to at all.
- `savedRangeAngular`/`savedRangeDispersion` -- `Ref`, the last captured range-field text per
  mode (`_dswCaptureRange`), or `nothing` before any capture.
"""
mutable struct DoubleSweepViewer
    fig::Figure
    sweepResult::Any
    contentGL::GridLayout
    axisMenu::Menu
    switchModeBtn::Union{Button,Nothing}
    statusLabel::Label
    currentIdx::Base.RefValue{Int}
    currentAngNo::Int
    mode::Symbol
    combinedViewer::Base.RefValue{Union{Nothing,CombinedViewer}}
    dispersionViewer::Base.RefValue{Union{Nothing,DispersionViewer}}
    updatingWidgets::Base.RefValue{Bool}
    savedRangeAngular::Base.RefValue{Any}
    savedRangeDispersion::Base.RefValue{Any}
    rebuild!::Function
    onSwitchMode!::Union{Function,Nothing}
end

"""
    _dswNearestAngle(angles, angleVal) -> Int

Index into `angles` nearest `angleVal` -- MATLAB's own literal `[~, obj.currentAngNo] =
min(abs(obj.sweepResult.angles - angle))`, used by both `onAngularModeSelected` and
`onAngularSelectionChanged` in the MATLAB source. Private near-duplicate of every other
viewer file's own `_xxxNearestAngle` helper (`_msvNearestAngle`, `_gcvNearestAngle`), kept as
this file's own copy per the "never depend on another viewer file's
private internals" convention.
"""
_dswNearestAngle(angles::AbstractVector{<:Real}, angleVal::Real) = argmin(abs.(angles .- angleVal))

# -----------------------------------------------------------------------------------------
# Private helpers (`_dsw`-prefixed -- private to this file, same per-file-duplication
# convention every other viewer file's header explains for its own `_xxx`-prefixed stepper
# wiring)
# -----------------------------------------------------------------------------------------

"""
    _dswWireStepper!(menu::Menu, prevBtn::Button, nextBtn::Button, n::Integer) -> Nothing

Wire a dropdown-stepper row's `<<`/`>>` buttons to cycle `menu.i_selected` through `1:n`,
wrapping -- duplicates every other viewer file's own `_xxxWireStepper!` verbatim.
"""
function _dswWireStepper!(menu, prevBtn, nextBtn, n::Integer)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] - 1, n)), prevBtn.clicks)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] + 1, n)), nextBtn.clicks)
    return nothing
end

"""
    _dswCaptureRange(viewer) -> Union{Nothing, Dict{Symbol,Tuple{String,String}}}

The min/max field TEXT of every `x`/`y`/`c` range a viewer has -- MATLAB's
`GuiCommon.captureAxisRange`. `nothing` for an `AngSweepViewer` drilled into a single angle
(its map fields are not mounted; the caller then keeps its previous snapshot).
"""
function _dswCaptureRange(viewer)
    viewer isa AngSweepViewer && viewer.isDrilledIn[] && return nothing
    s = Dict{Symbol,Tuple{String,String}}()
    for a in (:x, :y, :c)
        f = Symbol(a, :RangeField)
        hasproperty(viewer, f) || continue
        rf = getproperty(viewer, f)
        s[a] = (rf.loText[], rf.hiText[])
    end
    return s
end

"""
    _dswRestoreRange!(viewer, s) -> Nothing

Write a [`_dswCaptureRange`](@ref) snapshot into `viewer`'s mounted range boxes (`Makie.set!`,
the same path as typing, so its `GuiRangeField`s and axes follow) -- MATLAB's
`GuiCommon.restoreAxisRange`. A no-op for `nothing`.
"""
function _dswRestoreRange!(viewer, s)
    s === nothing && return nothing
    for (a, (lo, hi)) in s
        b = Symbol(a, :RangeBoxes)
        hasproperty(viewer, b) || continue
        boxes = getproperty(viewer, b)
        boxes === nothing && continue
        Makie.set!(boxes[1], lo)
        Makie.set!(boxes[2], hi)
    end
    return nothing
end

# -----------------------------------------------------------------------------------------
# _dswBuild -- shared shell for every variant. Differs only in which outer-axis vector/label
# is stepped over -- everything else (layout, rebuild mechanics, switch button) is
# identical across the three MATLAB sources.
# -----------------------------------------------------------------------------------------

function _dswBuild(sweepResult, axisVals::AbstractVector{<:Real}, axisLabel::AbstractString,
                    sweeps::AbstractVector; supportsDispersion::Bool = false, fig = nothing)
    figOut, mainLayout = guiHostLayout(fig; size = (1200, 900))
    topBarGL = guiGridLayout(mainLayout[1, 1])
    contentGL = guiGridLayout(mainLayout[2, 1])
    statusLabel = guiStatusLabel(mainLayout[3, 1])

    axisItems = [string(round(w; sigdigits = 4)) for w in axisVals]
    axisMenu, axisPrev, axisNext = guiDropdownStepperRow(topBarGL, 1, axisLabel, axisItems)
    switchBtn = supportsDispersion ?
        Button(topBarGL[2, 1], label = "Switch to dispersion view ->", halign = :left, tellwidth = false) : nothing

    state = DoubleSweepViewer(figOut, sweepResult, contentGL, axisMenu, switchBtn, statusLabel,
                               Ref(1), 1, :angular,
                               Ref{Union{Nothing,CombinedViewer}}(nothing),
                               Ref{Union{Nothing,DispersionViewer}}(nothing),
                               Ref(false), Ref{Any}(nothing), Ref{Any}(nothing),
                               () -> nothing, nothing)

    function _buildAngular!()
        cvOld = state.combinedViewer[]   # slice step: keep the outgoing map's ranges
        if cvOld !== nothing
            s = _dswCaptureRange(cvOld.angSweepViewer)
            s === nothing || (state.savedRangeAngular[] = s)
        end
        guiClearGridLayout!(state.contentGL)
        i = state.currentIdx[]
        cv = combinedViewer(sweeps[i]; windowMode = :combined, fig = state.contentGL[1, 1])
        _dswRestoreRange!(cv.angSweepViewer, state.savedRangeAngular[])
        if supportsDispersion
            # layer currentAngNo tracking ON TOP OF combinedViewer's own already-complete
            # bidirectional wiring -- see this file's header, "TWO MODES", for why this
            # chains onto the existing callback rather than replacing it.
            angles = sweeps[i].angles
            existingOnSelect = cv.angSweepViewer.onSelect[]
            cv.angSweepViewer.onSelect[] = (angle, modeNo, velSnap, velRaw) -> begin
                state.currentAngNo = _dswNearestAngle(angles, angle)
                existingOnSelect === nothing || existingOnSelect(angle, modeNo, velSnap, velRaw)
                return nothing
            end
            existingOnSelectionChanged = cv.modalShapeViewer.onSelectionChanged[]
            cv.modalShapeViewer.onSelectionChanged[] = (angle, vel) -> begin
                state.currentAngNo = _dswNearestAngle(angles, angle)
                existingOnSelectionChanged === nothing ||
                    existingOnSelectionChanged(angle, vel)
                return nothing
            end
        end
        state.combinedViewer[] = cv
        state.dispersionViewer[] = nothing
        state.mode = :angular
        state.switchModeBtn === nothing ||
            (state.switchModeBtn.label[] = "Switch to dispersion view ->")
        guiSetStatusOk!(state.statusLabel, "$(axisLabel) $(axisItems[i])")
        return nothing
    end
    state.rebuild! = _buildAngular!

    on(idx -> begin
        (state.updatingWidgets[] || state.mode !== :angular) && return nothing
        state.currentIdx[] = idx
        state.rebuild!()
    end, axisMenu.i_selected)
    _dswWireStepper!(axisMenu, axisPrev, axisNext, length(axisVals))

    if supportsDispersion
        function _buildDispersion!()
            guiClearGridLayout!(state.contentGL)
            dv = dispersionViewer(sweepResult, state.currentAngNo; fig = state.contentGL[1, 1])
            dv.goTo!(state.currentIdx[], 1)   # restore the wavelength/wavevector position --
                                                # branch has no cross-mode equivalent, stays at 1
            _dswRestoreRange!(dv, state.savedRangeDispersion[])
            state.dispersionViewer[] = dv
            state.combinedViewer[] = nothing
            state.mode = :dispersion
            state.switchModeBtn.label[] = "<- Switch to angular view"
            guiSetStatusOk!(state.statusLabel, "Dispersion curve at angle = " *
                "$(round(sweepResult.angles[state.currentAngNo]; sigdigits = 4)) deg")
            return nothing
        end

        function _onSwitchMode!()
            if state.mode === :angular
                cv = state.combinedViewer[]
                if cv !== nothing && cv.angSweepViewer isa AngSweepViewer &&
                   cv.angSweepViewer.isDrilledIn[]
                    cv.angSweepViewer.toggleDrillDown!()   # undrill first -- MATLAB's own
                                                             # returnToSweepView() guard
                end
                cv === nothing || (state.savedRangeAngular[] = _dswCaptureRange(cv.angSweepViewer))
                _buildDispersion!()
            else
                dv = state.dispersionViewer[]
                state.savedRangeDispersion[] = _dswCaptureRange(dv)
                state.currentIdx[] = dv.currentWavNo
                state.updatingWidgets[] = true
                try
                    state.axisMenu.i_selected[] = state.currentIdx[]
                finally
                    state.updatingWidgets[] = false
                end
                _buildAngular!()
            end
            return nothing
        end
        state.onSwitchMode! = _onSwitchMode!
        on(_ -> state.onSwitchMode!(), switchBtn.clicks)
    end

    _buildAngular!()
    return state
end

# -----------------------------------------------------------------------------------------
# doubleSweepViewer -- one method per outer-sweep family
# -----------------------------------------------------------------------------------------

"""
    doubleSweepViewer(sweepResult::CalcWavevecAngSweep2Layered; fig = nothing)
        -> DoubleSweepViewer

Outer wavelength-stepper shell over a [`combinedViewer`](@ref) ("angular" mode) or a
[`dispersionViewer`](@ref) ("dispersion" mode, toggled via the returned viewer's own
`switchModeBtn`/`onSwitchMode!`) -- port of `guliwerr/gui/DoubleSweepViewerGUI.m`, both modes.
`fig` -- `nothing` (default) for a standalone window, or a `GridPosition` to embed into.

Throws `ArgumentError` if `sweepResult` is not a genuine double sweep
(`numel(angles) < 2 || numel(wavelengths) < 2`) or has not been `calc()`'d yet
(`isempty(sweepResult.sweeps)`). Switching to dispersion mode additionally requires
`retentionLevel = :full` ([`dispersionViewer`](@ref)'s own requirement, via
[`trackDispersionBranches`](@ref)) -- not checked here, since angular mode alone does not need
it (matches `DoubleSweepViewerGUI`'s own class docblock: "requiring it unconditionally would
block viewing a large double sweep computed cheaply for its angular side").
"""
function GuliwerrJ.doubleSweepViewer(sweepResult::CalcWavevecAngSweep2Layered; fig = nothing)
    (length(sweepResult.angles) < 2 || length(sweepResult.wavelengths) < 2) && throw(ArgumentError(
        "doubleSweepViewer: needs a genuine double sweep (numel(angles)>1 AND " *
        "numel(wavelengths)>1); got $(length(sweepResult.angles)) angle(s), " *
        "$(length(sweepResult.wavelengths)) wavelength(s) -- use combinedViewer directly for " *
        "an angle-only sweep."))
    isempty(sweepResult.sweeps) && throw(ArgumentError(
        "doubleSweepViewer: sweepResult has not been calc()'d yet (sweeps is empty) -- " *
        "nothing to view."))
    return _dswBuild(sweepResult, sweepResult.wavelengths, "Wavelength:", sweepResult.sweeps;
                      supportsDispersion = true, fig = fig)
end

"""
    doubleSweepViewer(sweepResult::CalcWavevecAngSweep2Plate; fig = nothing) -> DoubleSweepViewer

Plate counterpart -- port of `guliwerr/gui/DoubleSweepViewerGUIPlate.m`, both modes.

Throws `ArgumentError` under the same two conditions as the Layered method, against `wavVecs`.
"""
function GuliwerrJ.doubleSweepViewer(sweepResult::CalcWavevecAngSweep2Plate; fig = nothing)
    (length(sweepResult.angles) < 2 || length(sweepResult.wavVecs) < 2) && throw(ArgumentError(
        "doubleSweepViewer: needs a genuine double sweep (numel(angles)>1 AND " *
        "numel(wavVecs)>1); got $(length(sweepResult.angles)) angle(s), " *
        "$(length(sweepResult.wavVecs)) wavevector(s) -- use combinedViewer directly for an " *
        "angle-only sweep."))
    isempty(sweepResult.sweeps) && throw(ArgumentError(
        "doubleSweepViewer: sweepResult has not been calc()'d yet (sweeps is empty) -- " *
        "nothing to view."))
    return _dswBuild(sweepResult, sweepResult.wavVecs, "Wavevector:", sweepResult.sweeps;
                      supportsDispersion = true, fig = fig)
end

"""
    doubleSweepViewer(sweepResult::CalcFreqAngSweep2Layered; fig = nothing) -> DoubleSweepViewer

isoFreq counterpart -- port of `guliwerr/gui/FreqAngSweep2LayeredViewerGUI.m`. No genuine-
double-sweep guard, and NO dispersion mode at all (`switchModeBtn`/`onSwitchMode!` are both
`nothing` on the returned viewer) -- MATLAB's own source has neither; see this file's header.

Throws `ArgumentError` if `sweepResult` has not been `calc()`'d yet
(`isempty(sweepResult.sweeps)`).
"""
function GuliwerrJ.doubleSweepViewer(sweepResult::CalcFreqAngSweep2Layered; fig = nothing)
    isempty(sweepResult.sweeps) && throw(ArgumentError(
        "doubleSweepViewer: sweepResult has not been calc()'d yet (sweeps is empty) -- " *
        "nothing to view."))
    return _dswBuild(sweepResult, sweepResult.omega0List, "omega0:", sweepResult.sweeps;
                      supportsDispersion = false, fig = fig)
end

"""
    setAxisRange!(viewer::DoubleSweepViewer, axis::Symbol, lims) -> Nothing

Set one range of the CURRENT mode's view from code -- forwarded to the angular combined
viewer, or in dispersion mode to the dispersion viewer (no colour axis there: `:c` throws).
As in MATLAB, the range is kept across a slice step and a mode switch and back.
"""
function GuliwerrJ.setAxisRange!(viewer::DoubleSweepViewer, axis::Symbol, lims)
    target = viewer.mode === :angular ? viewer.combinedViewer[] : viewer.dispersionViewer[]
    return GuliwerrJ.setAxisRange!(target, axis, lims)
end
