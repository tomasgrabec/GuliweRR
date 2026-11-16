# =========================================================================================
# GuiModalShapeViewer.jl -- viewer for modal shapes ("stacked modal-shape depth panels +
# steppers")
#
# Julia/Makie port of `guliwerr/gui/ModalShapeViewerGUI.m` (23 methods),
# `ModalShapeViewerGUIPlate.m` (14 methods) and `IsofreqModalShapeViewerGUI.m`.
# `include`d directly into
# `module GuliwerrJMakieExt` (see that file's own header), which already has `using
# GuliwerrJ`/`using Makie` in scope -- neither repeated here. Built on `GuiWidgets.jl`
# (`guiHostLayout`, `guiDropdownStepperRow`, `guiIndexStepperRow`, `guiStatusLabel`/
# `guiSetStatusOk!`/`guiSetStatusError!`) and `GuiSingleAngleViewer.jl`'s
# [`guiPlateModeTypeText`](@ref) (a `gui`-prefixed helper written there for reuse by the
# other viewer files). The actual depth-stacked
# panel drawing is `PlotsModalShape.jl`'s trio --
# `plotModalShape`, `plotGreenModalShape`, `plotGreenModalShapeLayered` -- called here, never
# reimplemented.
#
# ---------------------------------------------------------------------------------------
# THE BASE AND PLATE MATLAB CLASSES GENUINELY DIVERGED: WHAT IS AND IS NOT PORTED
# ---------------------------------------------------------------------------------------
# Unlike the other viewer families, base `ModalShapeViewerGUI` and its Plate sibling are not
# near-duplicates. The base class has 8 methods with NO Plate/isoFreq equivalent:
# `goToVelocity`, `onSourceChanged`, `onForceDirChanged`, `onMarkAsShown`,
# `classifyCurrentMode`, `findOverride`, `followPhysicalWave`, `refreshLimModePanel`.
# They split into two genuinely separate clusters:
#
#   CLUSTER 1 -- the Ritz/Green's-response SOURCE TOGGLE (`goToVelocity`, `onSourceChanged`,
#   `onForceDirChanged`): this is the SAME depth-stacked-panel display, just fed by
#   `plotGreenModalShape`/`plotGreenModalShapeLayered` instead of `plotModalShape` at an arbitrary
#   continuous velocity instead of a discrete Ritz mode index. It reuses the existing
#   plotting functions and slots onto the shared stepper/redraw machinery with no new
#   widget concepts (`guiIndexStepperRow`'s validator-swap trick exists in
#   `GuiWidgets.jl` FOR EXACTLY THIS FEATURE -- see that file's own `guiIndexStepperRow`
#   docstring, which names `ModalShapeViewerGUI.onSourceChanged` directly). IMPLEMENTED IN
#   FULL for the base (HS/Layered) variant.
#
#   CLUSTER 2 -- the session-local classification-override QC panel (`onMarkAsShown`,
#   `classifyCurrentMode`, `findOverride`, `refreshLimModePanel`) plus the angle-stepping
#   heuristic built on top of it (`followPhysicalWave`, which reads the override list FIRST
#   before falling back to `significantModesMap`). This is a distinct bookkeeping/QC layer,
#   not part of the depth-stacked-panel display itself or the mode/velocity stepper: it
#   neither draws a shape panel nor changes what value is being
#   plotted, only which mode a later angle-step lands on and a manual annotation of "this is
#   the SAW I was showing." Porting it faithfully requires a new mutable overrides list, a
#   dynamically-many-rows scrollable sub-panel with per-entry "Go To" buttons (Makie has no
#   scrollable-panel primitive at all, unlike MATLAB's `uipanel(...,'Scrollable','on')`), and
#   `significantModesMap` plumbing this file otherwise never
#   needs. NOT PORTED:
#     - `stepToAngle` for the Ritz-mode branch uses Plate's OWN
#       fallback -- "keep the same raw mode index" (`ModalShapeViewerGUIPlate.m`'s own
#       docblock: "a plain default rather than an unrequested mode-type-tracking heuristic").
#       This is a behaviour the MATLAB codebase itself already ships and documents as
#       acceptable for exactly this situation. The Green's-response branch of `stepToAngle`
#       needs NO change at all -- MATLAB's own Green's-response branch ALREADY keeps the same
#       velocity across an angle step (no mode/type to follow for a continuous value), which
#       is exactly what this port does too.
#     - The SAW/limiting mode list + "Mark as..." QC panel is not built. A user who needs to
#       see which modes are classified SAW/limiting at the current angle can open
#       `GuiSingleAngleViewer.jl`'s `singleAngleViewer` on the SAME result,
#       whose classification curve/info panel shows exactly that -- this file's own job is the
#       depth-resolved SHAPE, not a second copy of the mode browser.
#   `sweepResult.significantModesMap`/`SignificantMode` are
#   available in Julia, so nothing blocks adding Cluster 2 later.
#
# NO CLICK-TO-SELECT ANYWHERE IN THIS VIEWER -- confirmed against all three MATLAB sources.
# Unlike the angle-sweep and single-angle viewers, none of `ModalShapeViewerGUI`/`...Plate`/
# `IsofreqModalShapeViewerGUI`
# ever attaches a click handler to their own shape axes; every state change is
# stepper/dropdown/field-driven, or arrives externally via `goToMode`/`goToVelocity`/
# `goToSelection`/`goToPoint` (called by another viewer's own `OnModeSelected`/
# `OnPointSelected` wiring, e.g. `win2.OnModeSelected = @(angle, modeNo, ~, rawVel)
# win1.goToSelection(angle, modeNo, rawVel)`). This file reproduces exactly that: no
# `guiOnAxisClick` anywhere, and `goTo*!`/`stepToAngle!` are the public entry points the
# combined viewer (`GuiCombinedViewer.jl`) wiring calls directly, exactly mirroring MATLAB's own
# "public, not Access=private ... so tests can call these directly, and so
# AngSweepViewerGUI.OnModeSelected can be wired straight to goToSelection" comment.
#
# ---------------------------------------------------------------------------------------
# TWO FIGURES, ALWAYS -- `windowMode` ('combined'/'split') COLLAPSES ENTIRELY
# ---------------------------------------------------------------------------------------
# MATLAB's `windowMode` chooses between ONE uifigure (shapes-left/controls-right inside a
# single window, via a pre-built `{3 x nLay}` uiaxes grid) and TWO windows (a uifigure control
# panel plus a separate classic `figure()` for the shapes) -- needed only because a classic
# MATLAB `axes()` cannot be parented inside a `uifigure` at all.
#
# `GuiWidgets.jl`'s own header flags half of why this collapses in Julia:
# `shapeAxesGrid`/the `hAxCellIn` pre-built-grid mechanism dissolves entirely, because
# `plotModalShape`/`plotGreenModalShape`/`plotGreenModalShapeLayered`
# each own their OWN `Figure`/`Axis` construction via `_modalShapeFigureAndAxes`, which --
# this is the OTHER half -- calls `empty!(fig)` on whatever `Figure` it
# is given before rebuilding its own `GridLayout` inside it (see `PlotsModalShape.jl`).
# That means a plotting call can NEVER safely share a `Figure` with anything else this
# file draws (a controls panel, say): the NEXT redraw would wipe it out. So every
# `ModalShapeViewer` here owns exactly TWO Makie `Figure`s unconditionally:
#   - `fig` -- the controls host, embeddable exactly like every other viewer's `fig`
#     keyword (`nothing` for a standalone window, or a `GridPosition` to build controls into a
#     caller-owned region, via `guiHostLayout`).
#   - `shapeFig` -- ALWAYS its own standalone `Figure()`, NEVER embeddable, passed directly to
#     `plotModalShape`/`plotGreenModalShape`/`plotGreenModalShapeLayered` as their `fig`
#     keyword on every redraw.
# This is, in effect, MATLAB's OWN 'split' mode as the only case that now exists -- a real but
# narrow UI-convenience loss (one visual window instead of two arranged side by side), not a
# functionality loss: every angle/mode/velocity/source control and every value it can show is
# identical either way. A composite viewer wanting shape panels genuinely co-mingled with
# other content in ONE window needs its own solution to this same `empty!(fig)`-ownership
# constraint (e.g. not reusing this file's `shapeFig` wrapping at all, and calling the
# plotting functions directly against a `Figure` it owns and clears itself).
#
# `show(obj)` (MATLAB's `Visible = 'on'` toggle for both windows) is not ported, for the same
# reason `guiHostLayout`'s own docstring gives: a Makie `Figure` is never constructed
# hidden, so there is nothing to "reveal" -- this applies to `shapeFig` exactly as it
# applies to `fig`.
#
# ---------------------------------------------------------------------------------------
# DISPATCH DESIGN -- ONE struct, multiple-dispatch `modalShapeViewer` methods (mirrors
# `SingleAngleViewer`'s own precedent in `GuiSingleAngleViewer.jl`)
# ---------------------------------------------------------------------------------------
# `ModalShapeViewer` (below) is the ONE handle type returned by every variant --
# `modalShapeViewer(result::Union{CalcSingleHS,CalcSingleLayered,CalcAngSweepHS,
# CalcAngSweepLayered}; fig)`, `modalShapeViewer(result::Union{CalcSinglePlate,
# CalcAngSweepPlate}; fig)`, `modalShapeViewer(result::Union{CalcIsofreqSingleLayered,
# CalcIsofreqAngSweepLayered}; fig)` -- rather than three (or six) near-copies of a bespoke
# type each. Fields that only apply to some variants (`sourceMenu`/`forceDirMenu`/`modeBox`/
# `kParBox`/`modeTypeLabel`/`goToMode!`/`goToVelocity!`/`goToSelection!`/`goToPoint!`) are
# `nothing` -- ABSENT, not disabled -- for the variants that have no such concept, exactly
# `SingleAngleViewer.graphTypeMenu`'s own precedent for Plate/isoFreq.
#
# SINGLE-POINT-AS-DEGENERATE-1-ANGLE-SWEEP -- ported verbatim, not reinvented. All three
# MATLAB constructors accept EITHER a genuine multi-angle sweep OR a single-point result,
# adapting the latter into "the exact shape one element of a real sweep would have" (the
# MATLAB source's own words). This file does the same: `modalShapeViewer(::CalcSingleHS)`
# builds `angles = [result.angle]`, `cs = [result]` and hands them to the SAME internal
# `_msvBuildBase` a genuine `CalcAngSweepHS` also calls -- one redraw/stepping implementation
# per family, not a SingleAngle-specific rendering path (MATLAB's own stated reason). Unlike
# MATLAB's adapter, `structure`/`significantModesMap` are NOT threaded through: this file never
# pre-builds a `shapeAxesGrid` from `structure.calcParams` (dissolved, see `GuiWidgets.jl`'s
# header) and never reads `significantModesMap` at all (Cluster 2, not ported, above) -- so the
# adapter here only needs `angles`/`cs` (base/Plate) or `angles`/`structure`/`omega0`/
# `kParRange`/`structureAlreadyRotated` (isoFreq), each built inline at the six outer
# `modalShapeViewer` methods below, not via a separate NamedTuple factory.
#
# ---------------------------------------------------------------------------------------
# REENTRANCY GUARD -- Makie `Observable`s fire on EVERY assignment, MATLAB `.Value` does not
# ---------------------------------------------------------------------------------------
# `goToMode!`/`goToVelocity!`/`goToPoint!`/`stepToAngle!` are public entry points meant to be
# called externally (by a test, or by cross-viewer wiring) as well as internally (by
# this file's own stepper-button/dropdown/field handlers) -- exactly MATLAB's `goToMode`/
# `goToVelocity`/`goToPoint`. Each syncs EVERY widget's displayed value to the new state
# (`angleMenu`, `sourceMenu`, `modeBox`/`kParBox`, ...). In MATLAB this is safe by construction:
# assigning a `uidropdown`/`uieditfield`'s `.Value` NEVER fires that widget's own
# `ValueChangedFcn` unless the value genuinely changes from the user's perspective, and even
# then only once. Makie's `Observable` has no such guard -- `menu.i_selected[] = v` calls every
# attached `on(...)` listener unconditionally, every time, including listeners THIS very sync
# step itself attached to drive redraws on genuine user interaction. Without a guard, a single
# `goToMode!` call could cascade into several redundant (though not incorrect -- every
# cascaded redraw reads the SAME already-updated state) nested redraws, or worse, tangle two
# widgets' handlers into firing each other.
#
# `updatingWidgets::Base.RefValue{Bool}` reproduces MATLAB's "programmatic assignment never
# re-triggers its own handler" contract, generically, in ONE place: every widget-sync helper
# below sets it `true` for the duration of its own assignments, and every USER-interaction
# handler checks it first and bails out (`state.updatingWidgets[] && return nothing`) if a
# programmatic sync is already in progress. The same pattern suits any viewer that exposes a
# "jump to this exact state" public API syncing several widgets at once.
#
# ---------------------------------------------------------------------------------------
# `guiIndexStepperRow`'s `onChange` FORWARDER IS NOT USED HERE -- it hardcodes `tryparse(Int,
# ...)`
# ---------------------------------------------------------------------------------------
# `GuiWidgets.jl`'s `guiIndexStepperRow` builds the label+field+steppers+trailing-label row
# this file uses for the base variant's Mode/velocity field, and its own docstring explicitly
# anticipates relaxing `box.validator[]` to `Float64` for the continuous-velocity
# Green's-response repurposing. Its OPTIONAL `onChange` callback
# argument, however, parses the field's text via a hardcoded `tryparse(Int, s)` regardless of
# the box's OWN (possibly relaxed) validator, so it silently never fires for a fractional
# string like `"3.45"` once `onSourceChanged` has switched the field over to velocity entry.
# This file does not pass `onChange` to `guiIndexStepperRow` at all; it wires
# `modeBox.stored_string` directly with its own listener that always parses via
# `tryparse(Float64, ...)` (which also accepts plain-integer strings like `"5"` -- one parser
# covers both of this field's two meanings).
#
# ---------------------------------------------------------------------------------------
# `<<`/`>>` BUTTONS ON `sourceMenu`/`forceDirMenu` -- carrying `GuiSingleAngleViewer.jl`'s
# deviation 1 forward, not introducing a new one
# ---------------------------------------------------------------------------------------
# MATLAB's `sourceDropdown`/`forceDirDropdown` (base) and `forceDirDropdown` (isoFreq) are
# plain `uidropdown`s built by hand, with no stepper buttons at all -- ONLY the angle row
# (`GuiCommon.dropdownStepperRow`, wired to `onPrevAngle`/`onNextAngle`) and the mode row
# (`GuiCommon.indexStepperRow`, wired to `onPrevMode`/`onNextMode`) get real `<<`/`>>` buttons
# in the MATLAB source. This file builds EVERY dropdown here (angle, source, force-direction)
# via [`guiDropdownStepperRow`](@ref), which always attaches a `<<`/`>>` pair -- so
# `sourceMenu`/`forceDirMenu` gain stepper buttons MATLAB's own versions of this class do not
# have. This is the exact same convenience `GuiSingleAngleViewer.jl`'s header names as
# its own deviation 1 ("a small added convenience with no MATLAB counterpart to contradict,
# not a removed behaviour") -- carried forward here for consistency across viewer files,
# not a second, independently-invented deviation. The angle and mode rows' OWN `<<`/`>>`
# buttons ARE a faithful port (MATLAB already has them there), wired to `max(...,1)`/
# `min(..., length(angles))` CLAMPING (not the wrapping `mod1` `_msvWireMenuStepper!` uses for
# source/force-dir) -- exactly `onPrevAngle`/`onNextAngle`'s own clamp, since angle has a real
# start/end, unlike a 2- or 3-item mode-vocabulary dropdown where wrapping reads more naturally
# (again matching `GuiSingleAngleViewer.jl`'s own precedent for ITS graph-type/criterion
# dropdowns).
# =========================================================================================

const MSV_VEL_STEP = 0.05   # km/s per << >> step in Green's-response mode -- MATLAB's own
                              # ModalShapeViewerGUI.VEL_STEP constant, unchanged

# -----------------------------------------------------------------------------------------
# Private pure/near-pure helpers (`_msv`-prefixed -- private to this file, same convention
# `GuiSingleAngleViewer.jl`'s own header explains for its `_sav`-prefixed helpers: the
# viewer files must never depend on each other's internals)
# -----------------------------------------------------------------------------------------

"""
    _msvNearestAngle(angles, angleVal) -> Int

Index into `angles` nearest `angleVal` -- the pure core of every variant's `goToMode`/
`goToVelocity`/`goToPoint` (MATLAB's `[~, angNo] = min(abs(obj.sweepResult.angles -
angleVal))`).
"""
_msvNearestAngle(angles::AbstractVector{<:Real}, angleVal::Real) = argmin(abs.(angles .- angleVal))

"""
    _msvClampMode(nCompModes, modeNo) -> Int

Round-and-clamp a requested mode number into `1:nCompModes` -- MATLAB's `min(max(round(modeNo),
1), nM)`. `nCompModes` is the STORED `RitzSolver`'s own count, which stays full-size even under
partial retention (only `eigVecs_con` shrinks -- see [`_msvResolveDisplayRitz`](@ref)'s own
docstring), so this clamp is correct regardless of `retentionLevel`.
"""
_msvClampMode(nCompModes::Integer, modeNo::Real) = clamp(round(Int, modeNo), 1, nCompModes)

"""
    _msvResolveDisplayRitz(cs_i::Union{CalcSingleHS,CalcSingleLayered}, cache, angNo)
        -> RitzSolverSAW
    _msvResolveDisplayRitz(cs_i::CalcSinglePlate, cache, angNo) -> RitzSolverLamb

Return a `RitzSolver` with a FULL `eigVecs_con` for this angle, recomputing on demand when
`cs_i.retentionLevel != :full` pruned it away -- direct port of `ModalShapeViewerGUI.
resolveDisplayRitz`/`ModalShapeViewerGUIPlate.resolveDisplayRitz`. `plotModalShape` itself
deliberately errors rather than guessing when `eigVecs_con` has too few columns (see
`PlotsModalShape.jl`'s own `plotModalShape` docstring); this is the one layer up that keeps
that contract while making the viewer transparent to `retentionLevel`.

Recomputed ONCE per angle then cached in `cache[angNo]` (a per-viewer `Vector{Any}`,
`nothing` until filled) -- a genuine Ritz solve, so repeated mode/velocity stepping at one
angle must not re-solve on every redraw. NEVER written back into `cs_i` itself, which would
defeat `retentionLevel`'s whole purpose (matches MATLAB exactly).

Angle `0.0` is passed to the solver constructor, not the real propagation angle: `cs_i.
structure` is ALREADY ROTATED for this angle (the same convention `CalcSingleHS.calcRitz`
itself uses) -- passing the real angle here would double-rotate it.
"""
function _msvResolveDisplayRitz(cs_i::Union{CalcSingleHS,CalcSingleLayered},
                                 cache::Vector{Any}, angNo::Integer)
    cs_i.retentionLevel === :full && return cs_i.ritzSolver
    cache[angNo] !== nothing && return cache[angNo]
    rs = RitzSolverSAW(; structure = cs_i.structure, angleDeg = 0.0)
    ritz = calc(rs, true, cs_i.precompStructForRitz)
    cache[angNo] = ritz
    return ritz
end

function _msvResolveDisplayRitz(cs_i::CalcSinglePlate, cache::Vector{Any}, angNo::Integer)
    cs_i.retentionLevel === :full && return cs_i.ritzSolver
    cache[angNo] !== nothing && return cache[angNo]
    rs = RitzSolverLamb(; structure = cs_i.structure, angleDeg = 0.0,
                           realThickness = cs_i.realThickness, realWavVec = cs_i.realWavVec)
    ritz = calc(rs, true, cs_i.precompStructForRitz)
    cache[angNo] = ritz
    return ritz
end

"""
    _msvAngleItems(angles) -> Vector{String}

`"%.1f"`-equivalent labels for the angle dropdown (MATLAB's `arrayfun(@(a) sprintf('%.1f', a),
angles)`).
"""
_msvAngleItems(angles::AbstractVector{<:Real}) = [string(round(a; digits = 1)) for a in angles]

"""
    _msvWireMenuStepper!(menu::Menu, prevBtn::Button, nextBtn::Button, n::Integer) -> Nothing

Wire a dropdown-stepper row's `<<`/`>>` buttons to cycle `menu.i_selected` through `1:n`,
wrapping via `mod1` -- private near-duplicate of `GuiSingleAngleViewer.jl`'s own
`_savWireStepper!` (see this file's header for why every viewer file keeps its own copy of
such a helper rather than reaching into another file's `_xxx`-prefixed internals). Used for
this file's `angleMenu`/`sourceMenu`/`forceDirMenu` steppers -- MATLAB's own `<<`/`>>` buttons
on every one of those dropdown rows.
"""
function _msvWireMenuStepper!(menu, prevBtn, nextBtn, n::Integer)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] - 1, n)), prevBtn.clicks)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] + 1, n)), nextBtn.clicks)
    return nothing
end

# -----------------------------------------------------------------------------------------
# ModalShapeViewer -- the one handle type returned by every variant
# -----------------------------------------------------------------------------------------

"""
    ModalShapeViewer

Handles for one "stacked modal-shape depth panels + steppers" viewer, returned by
every [`modalShapeViewer`](@ref) method regardless of which family built it (see this file's
header for the dispatch design and the Cluster-1/Cluster-2 scoping decision).

Fields:
- `fig` -- the controls-host `Figure` (embeddable; see this file's header, "TWO FIGURES,
  ALWAYS").
- `shapeFig` -- the ALWAYS-standalone `Figure` the depth-stacked panels are drawn into on
  every redraw (never embeddable -- see this file's header).
- `angles` -- propagation angles, one per element of `cs` (base/Plate) or the sweep itself
  (isoFreq).
- `cs` -- `Vector{CalcSingleHS}`/`Vector{CalcSingleLayered}`/`Vector{CalcSinglePlate}` for the
  base/Plate variants; `nothing` for isoFreq (untyped field so this ONE struct serves every
  variant, matching `SingleAngleViewer.result`'s own precedent).
- `isofreqStructure`/`isofreqOmega0`/`isofreqKParRange`/`isofreqAlreadyRotated` -- isoFreq-only
  state (`nothing`/`NaN`/empty/`false` otherwise); see `IsofreqModalShapeViewerGUI`'s own
  "IMPORTANT DIFFERENCE" docblock note for what `isofreqAlreadyRotated` guards against (the
  rotation convention: a genuine sweep's `structure` is unrotated, a wrapped single-point
  result's is already rotated to its own angle).
- `angleMenu`, `angleStepperBtns` -- the angle dropdown-stepper row, every variant.
- `modeBox`, `modeStepperBtns`, `velLabel` -- the Mode/velocity stepper row (base/Plate).
  For isoFreq `modeBox` is `nothing`, while `modeStepperBtns`/`velLabel` hold the kParRange
  row's `<<`/`>>` buttons (one kParRange sample per step) and its `v = omega0/kParRange` label.
- `kParBox` -- the kParRange field (isoFreq only); `nothing` for base/Plate.
- `sourceMenu` -- the Ritz-modal-shape/Green's-response source dropdown (base only, Cluster 1);
  `nothing` for Plate/isoFreq, which have exactly one source each (absent, not disabled,
  matching MATLAB's own precedent for those two classes).
- `forceDirMenu` -- the force-direction dropdown (base + isoFreq, both of which can display a
  Green's response); `nothing` for Plate (no Green's response at all for a free plate).
- `modeTypeLabel` -- the current mode's type/confidence label (Plate only, replacing base's
  SAW/limiting list -- Cluster 2, not ported); `nothing` for base/isoFreq.
- `statusLabel` -- ok/error status line.
- `currentAngNo`/`currentModeNo`/`currentVel`/`currentKParVal` -- current selection state,
  mirroring MATLAB's own `SetAccess = private` properties of the same names (`currentModeNo`/
  `currentVel` unused -- left at their initial values -- for isoFreq; `currentKParVal` unused
  for base/Plate).
- `fullRitzCache` -- per-angle `Vector{Any}` cache for [`_msvResolveDisplayRitz`](@ref)
  (base/Plate only; empty for isoFreq, which never runs a Ritz solve at all).
- `onSelectionChanged` -- caller-settable callback (`Ref{Union{Nothing,Function}}`, `nothing`
  initially), MATLAB's `OnSelectionChanged` property. Called after every successful state
  change: `(angle, vel)` for base/Plate, `(angle, kParVal)` for isoFreq -- set it with
  `viewer.onSelectionChanged[] = (args...) -> ...`.
- `updatingWidgets` -- reentrancy guard (see this file's header); not meant to be read/set by
  callers, exposed only because every field on this struct is.
- `redraw!` -- `() -> Nothing`, redraws `shapeFig` for the current selection; also callable
  directly (e.g. from a test) without going through a real widget event, matching
  `SingleAngleViewer.redraw!`'s own precedent.
- `goToMode!` -- `(angleVal, modeNo) -> Nothing`, MATLAB's `goToMode` (base/Plate only;
  `nothing` for isoFreq, which has no discrete Ritz mode).
- `goToVelocity!` -- `(angleVal, velVal) -> Nothing`, MATLAB's `goToVelocity` (base only,
  Cluster 1; `nothing` for Plate/isoFreq).
- `goToSelection!` -- `(angleVal, modeNo, rawVel) -> Nothing`, MATLAB's `goToSelection` (base/
  Plate only -- the external entry point an angle-sweep/single-angle viewer's `OnModeSelected` wiring calls;
  `nothing` for isoFreq, which uses `goToPoint!` instead, matching `IsofreqModalShapeViewerGUI`
  having no `goToSelection` method at all).
- `goToPoint!` -- `(angleVal, kParVal) -> Nothing`, MATLAB's `goToPoint` (isoFreq only;
  `nothing` for base/Plate).
- `stepToAngle!` -- `(newAngNo) -> Nothing`, MATLAB's `stepToAngle` -- present for every
  variant (every variant steps through angles, even though WHAT it keeps fixed while doing so
  differs -- see this file's header, Cluster 2).
"""
mutable struct ModalShapeViewer
    fig::Figure
    shapeFig::Figure

    angles::Vector{Float64}
    cs::Any

    isofreqStructure::Any
    isofreqOmega0::Float64
    isofreqKParRange::Vector{Float64}
    isofreqAlreadyRotated::Bool

    angleMenu::Menu
    angleStepperBtns::Tuple{Button,Button}
    modeBox::Union{Textbox,Nothing}
    modeStepperBtns::Union{Tuple{Button,Button},Nothing}
    velLabel::Union{Label,Nothing}
    kParBox::Union{Textbox,Nothing}
    sourceMenu::Union{Menu,Nothing}
    forceDirMenu::Union{Menu,Nothing}
    modeTypeLabel::Union{Label,Nothing}
    statusLabel::Label

    currentAngNo::Int
    currentModeNo::Int
    currentVel::Float64
    currentKParVal::Float64

    fullRitzCache::Vector{Any}

    onSelectionChanged::Base.RefValue{Union{Nothing,Function}}
    updatingWidgets::Base.RefValue{Bool}

    redraw!::Function
    goToMode!::Union{Function,Nothing}
    goToVelocity!::Union{Function,Nothing}
    goToSelection!::Union{Function,Nothing}
    goToPoint!::Union{Function,Nothing}
    stepToAngle!::Function
end

# -----------------------------------------------------------------------------------------
# _msvBuildBase -- shared body for the base variant (CalcSingleHS/CalcSingleLayered/
# CalcAngSweepHS/CalcAngSweepLayered), Cluster 1 (Ritz/Green's-response source toggle) in full
# -----------------------------------------------------------------------------------------

function _msvBuildBase(angles::Vector{Float64}, cs::AbstractVector; fig = nothing)
    fig, layout = guiHostLayout(fig; size = (620, 380))
    mainLayout = guiGridLayout(layout[1, 1])

    angleMenu, anglePrev, angleNext =
        guiDropdownStepperRow(mainLayout, 1, "Angle:", _msvAngleItems(angles))
    modeBox, modePrev, modeNext, velLabel = guiIndexStepperRow(mainLayout, 2, "Mode:")
    sourceMenu, sourcePrev, sourceNext =
        guiDropdownStepperRow(mainLayout, 3, "Source:", ["Ritz modal shape", "Green's response"])
    forceDirMenu, forceDirPrev, forceDirNext =
        guiDropdownStepperRow(mainLayout, 4, "Force dir:",
                               ["x (propag. dir.)", "y (in-plane)", "z (out-of-plane)"])
    statusLabel = guiStatusLabel(mainLayout[5, 1])

    shapeFig = Figure(size = (900, 500))

    state = ModalShapeViewer(fig, shapeFig, angles, cs,
                              nothing, NaN, Float64[], false,
                              angleMenu, (anglePrev, angleNext), modeBox,
                              (modePrev, modeNext), velLabel, nothing,
                              sourceMenu, forceDirMenu, nothing, statusLabel,
                              1, 1, NaN, NaN,
                              Any[nothing for _ in 1:length(angles)],
                              Ref{Union{Nothing,Function}}(nothing), Ref(false),
                              () -> nothing, nothing, nothing, nothing, nothing,
                              (_...) -> nothing)

    function _redraw!()
        angNo = state.currentAngNo
        cs_i = state.cs[angNo]
        ok, errMsg = true, ""
        try
            if state.sourceMenu.i_selected[] == 1
                ritz = _msvResolveDisplayRitz(cs_i, state.fullRitzCache, angNo)
                plotModalShape(ritz, state.currentModeNo, :boosted; fig = state.shapeFig)
            else
                forceDir = state.forceDirMenu.i_selected[]
                greenObj = calc(GreenSolver, cs_i.structure, 0.0, [state.currentVel];
                                 forceDir = forceDir, displDir = 3, fullOutput = true)
                if length(cs_i.structure.calcParams) > 1
                    plotGreenModalShapeLayered(greenObj, cs_i.structure, forceDir, :boosted;
                                               fig = state.shapeFig)
                else
                    plotGreenModalShape(greenObj, forceDir; fig = state.shapeFig)
                end
            end
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end
        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    # sets EVERY widget's displayed value from current state, under the reentrancy guard --
    # see this file's header, "REENTRANCY GUARD"
    function _syncWidgets!(desiredSourceIdx::Integer)
        state.updatingWidgets[] = true
        try
            state.sourceMenu.i_selected[] = desiredSourceIdx
            state.angleMenu.i_selected[] = state.currentAngNo
            cs_i = state.cs[state.currentAngNo]
            if desiredSourceIdx == 1
                state.modeBox.validator[] = Int
                state.modeBox.stored_string[] = string(state.currentModeNo)
                vel = cs_i.ritzSolver.vels[state.currentModeNo]
                state.velLabel.text[] = "velocity: $(round(vel, digits = 4)) km/s"
            else
                state.modeBox.validator[] = Float64
                state.modeBox.stored_string[] = string(round(state.currentVel, digits = 4))
                state.velLabel.text[] = ""
            end
        finally
            state.updatingWidgets[] = false
        end
        return nothing
    end

    function _syncControlsAndRedraw!(desiredSourceIdx::Integer)
        _syncWidgets!(desiredSourceIdx)
        _redraw!()
        cb = state.onSelectionChanged[]
        cb === nothing || cb(state.angles[state.currentAngNo], state.currentVel)
        return nothing
    end

    function _goToMode!(angleVal::Real, modeNo::Real)
        angNo = _msvNearestAngle(state.angles, angleVal)
        cs_i = state.cs[angNo]
        modeNo_ = _msvClampMode(cs_i.ritzSolver.nCompModes, modeNo)
        state.currentAngNo = angNo
        state.currentModeNo = modeNo_
        state.currentVel = cs_i.ritzSolver.vels[modeNo_]
        _syncControlsAndRedraw!(1)
        return nothing
    end
    state.goToMode! = _goToMode!

    function _goToVelocity!(angleVal::Real, velVal::Real)
        angNo = _msvNearestAngle(state.angles, angleVal)
        cs_i = state.cs[angNo]
        state.currentAngNo = angNo
        state.currentVel = Float64(velVal)
        # keep currentModeNo pointed at the nearest Ritz mode purely so a later switch BACK to
        # Ritz modal shape starts somewhere sensible (matches MATLAB's goToVelocity exactly)
        ritz = cs_i.ritzSolver
        nM = ritz.nCompModes
        state.currentModeNo = argmin(abs.(view(ritz.vels, 1:nM) .- velVal))
        _syncControlsAndRedraw!(2)
        return nothing
    end
    state.goToVelocity! = _goToVelocity!

    function _goToSelection!(angleVal::Real, modeNo::Real, rawVel::Real)
        if state.sourceMenu.i_selected[] == 1
            state.goToMode!(angleVal, modeNo)
        else
            state.goToVelocity!(angleVal, rawVel)
        end
        return nothing
    end
    state.goToSelection! = _goToSelection!

    function _stepToAngle!(newAngNo::Integer)
        if state.sourceMenu.i_selected[] == 1
            # Cluster 2 (followPhysicalWave) is not ported -- keep the SAME RAW MODE INDEX, exactly
            # ModalShapeViewerGUIPlate's own documented fallback (see this file's
            # header, "WHAT IS AND IS NOT PORTED")
            state.goToMode!(state.angles[newAngNo], state.currentModeNo)
        else
            # no discrete mode/type to follow for an arbitrary velocity -- same velocity,
            # different angle (matches MATLAB's OWN Green's-response stepToAngle exactly,
            # not a simplification introduced here)
            state.goToVelocity!(state.angles[newAngNo], state.currentVel)
        end
        return nothing
    end
    state.stepToAngle! = _stepToAngle!

    on(idx -> (state.updatingWidgets[] || state.stepToAngle!(idx)), angleMenu.i_selected)
    on(idx -> (state.updatingWidgets[] || _syncControlsAndRedraw!(idx)), sourceMenu.i_selected)
    on(idx -> (state.updatingWidgets[] || _redraw!()), forceDirMenu.i_selected)
    on(modeBox.stored_string) do s
        state.updatingWidgets[] && return nothing
        v = tryparse(Float64, s)
        v === nothing && return nothing
        if state.sourceMenu.i_selected[] == 1
            state.goToMode!(state.angles[state.currentAngNo], v)
        else
            state.goToVelocity!(state.angles[state.currentAngNo], v)
        end
        return nothing
    end

    on(_ -> state.stepToAngle!(max(state.currentAngNo - 1, 1)), anglePrev.clicks)
    on(_ -> state.stepToAngle!(min(state.currentAngNo + 1, length(state.angles))), angleNext.clicks)
    on(_ -> begin
        if state.sourceMenu.i_selected[] == 1
            state.goToMode!(state.angles[state.currentAngNo], state.currentModeNo - 1)
        else
            state.goToVelocity!(state.angles[state.currentAngNo], state.currentVel - MSV_VEL_STEP)
        end
    end, modePrev.clicks)
    on(_ -> begin
        if state.sourceMenu.i_selected[] == 1
            state.goToMode!(state.angles[state.currentAngNo], state.currentModeNo + 1)
        else
            state.goToVelocity!(state.angles[state.currentAngNo], state.currentVel + MSV_VEL_STEP)
        end
    end, modeNext.clicks)
    _msvWireMenuStepper!(sourceMenu, sourcePrev, sourceNext, 2)
    _msvWireMenuStepper!(forceDirMenu, forceDirPrev, forceDirNext, 3)

    state.goToMode!(angles[1], 1)
    return state
end

"""
    modalShapeViewer(result::Union{CalcSingleHS,CalcSingleLayered}; fig = nothing)
        -> ModalShapeViewer

Stacked modal-shape depth-panel viewer for an already-`calcAll()`'d `CalcSingleHS`/
`CalcSingleLayered` result -- port of `guliwerr/gui/ModalShapeViewerGUI.m`, packaged as a
degenerate 1-angle "sweep" (see this file's header). `fig` -- `nothing` (default) for a
standalone controls window, or anything [`guiHostLayout`](@ref) accepts as `parent` to embed
the CONTROLS into an existing figure (`shapeFig` is always its own standalone `Figure` --
see this file's header, "TWO FIGURES, ALWAYS").

Throws `ArgumentError` if `result.classification === nothing` (not yet `calcAll()`'d).
"""
function GuliwerrJ.modalShapeViewer(result::Union{CalcSingleHS,CalcSingleLayered}; fig = nothing)
    result.classification === nothing && throw(ArgumentError(
        "modalShapeViewer: result has not been calcAll()'d yet (classification is nothing) " *
        "-- nothing to view."))
    return _msvBuildBase([result.angle], [result]; fig = fig)
end

"""
    modalShapeViewer(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}; fig = nothing)
        -> ModalShapeViewer

Stacked modal-shape depth-panel viewer for an already-`calc()`'d `CalcAngSweepHS`/
`CalcAngSweepLayered` sweep -- port of `guliwerr/gui/ModalShapeViewerGUI.m`.

Throws `ArgumentError` if `sweepResult.significantModesMap` is empty (not yet `calc()`'d) --
this file never reads `significantModesMap` itself (Cluster 2, not ported), but its emptiness is
still the correct "has this been calculated" signal, exactly as in MATLAB.
"""
function GuliwerrJ.modalShapeViewer(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}; fig = nothing)
    isempty(sweepResult.significantModesMap) && throw(ArgumentError(
        "modalShapeViewer: sweepResult has not been calc()'d yet (significantModesMap is " *
        "empty) -- nothing to view."))
    return _msvBuildBase(sweepResult.angles, sweepResult.cs; fig = fig)
end

# -----------------------------------------------------------------------------------------
# _msvBuildPlate -- shared body for the Plate variant (CalcSinglePlate/CalcAngSweepPlate).
# No source toggle, no force-direction control (no Green's response for a free plate) -- ABSENT
# not disabled, matching ModalShapeViewerGUIPlate.m's own precedent exactly.
# -----------------------------------------------------------------------------------------

function _msvBuildPlate(angles::Vector{Float64}, cs::AbstractVector; fig = nothing)
    fig, layout = guiHostLayout(fig; size = (560, 260))
    mainLayout = guiGridLayout(layout[1, 1])

    angleMenu, anglePrev, angleNext =
        guiDropdownStepperRow(mainLayout, 1, "Angle:", _msvAngleItems(angles))
    modeBox, modePrev, modeNext, velLabel = guiIndexStepperRow(mainLayout, 2, "Mode:")
    modeTypeLabel = Label(mainLayout[3, 1], ""; halign = :left, tellwidth = false)
    statusLabel = guiStatusLabel(mainLayout[4, 1])

    shapeFig = Figure(size = (900, 500))

    state = ModalShapeViewer(fig, shapeFig, angles, cs,
                              nothing, NaN, Float64[], false,
                              angleMenu, (anglePrev, angleNext), modeBox,
                              (modePrev, modeNext), velLabel, nothing,
                              nothing, nothing, modeTypeLabel, statusLabel,
                              1, 1, NaN, NaN,
                              Any[nothing for _ in 1:length(angles)],
                              Ref{Union{Nothing,Function}}(nothing), Ref(false),
                              () -> nothing, nothing, nothing, nothing, nothing,
                              (_...) -> nothing)

    function _redraw!()
        angNo = state.currentAngNo
        cs_i = state.cs[angNo]
        ok, errMsg = true, ""
        try
            ritz = _msvResolveDisplayRitz(cs_i, state.fullRitzCache, angNo)
            plotModalShape(ritz, state.currentModeNo, :boosted; fig = state.shapeFig)
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end
        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _syncWidgets!()
        state.updatingWidgets[] = true
        try
            state.angleMenu.i_selected[] = state.currentAngNo
            state.modeBox.stored_string[] = string(state.currentModeNo)
            ritz = state.cs[state.currentAngNo].ritzSolver
            vel = ritz.vels[state.currentModeNo]
            state.velLabel.text[] = "velocity: $(round(vel, digits = 4)) km/s"
            state.modeTypeLabel.text[] =
                "Type: $(guiPlateModeTypeText(ritz, state.currentModeNo))"
        finally
            state.updatingWidgets[] = false
        end
        return nothing
    end

    function _syncControlsAndRedraw!()
        _syncWidgets!()
        _redraw!()
        cb = state.onSelectionChanged[]
        cb === nothing || cb(state.angles[state.currentAngNo], state.currentVel)
        return nothing
    end

    function _goToMode!(angleVal::Real, modeNo::Real)
        angNo = _msvNearestAngle(state.angles, angleVal)
        cs_i = state.cs[angNo]
        modeNo_ = _msvClampMode(cs_i.ritzSolver.nCompModes, modeNo)
        state.currentAngNo = angNo
        state.currentModeNo = modeNo_
        state.currentVel = cs_i.ritzSolver.vels[modeNo_]
        _syncControlsAndRedraw!()
        return nothing
    end
    state.goToMode! = _goToMode!

    # MATLAB's own 4-arg duck-typed signature (4th arg unused) -- kept so an angle-sweep/single-angle
    # viewer's `OnModeSelected` wiring can call every family's `goToSelection!` identically
    function _goToSelection!(angleVal::Real, modeNo::Real, ::Real)
        state.goToMode!(angleVal, modeNo)
        return nothing
    end
    state.goToSelection! = _goToSelection!

    function _stepToAngle!(newAngNo::Integer)
        # keeps the SAME RAW MODE INDEX -- ModalShapeViewerGUIPlate.m's own documented choice
        # (no followPhysicalWave analogue for Plate; see its class docblock), not a
        # simplification introduced by this port
        state.goToMode!(state.angles[newAngNo], state.currentModeNo)
        return nothing
    end
    state.stepToAngle! = _stepToAngle!

    on(idx -> (state.updatingWidgets[] || state.stepToAngle!(idx)), angleMenu.i_selected)
    on(modeBox.stored_string) do s
        state.updatingWidgets[] && return nothing
        v = tryparse(Float64, s)
        v === nothing && return nothing
        state.goToMode!(state.angles[state.currentAngNo], v)
        return nothing
    end
    on(_ -> state.stepToAngle!(max(state.currentAngNo - 1, 1)), anglePrev.clicks)
    on(_ -> state.stepToAngle!(min(state.currentAngNo + 1, length(state.angles))), angleNext.clicks)
    on(_ -> state.goToMode!(state.angles[state.currentAngNo], state.currentModeNo - 1),
       modePrev.clicks)
    on(_ -> state.goToMode!(state.angles[state.currentAngNo], state.currentModeNo + 1),
       modeNext.clicks)

    state.goToMode!(angles[1], 1)
    return state
end

"""
    modalShapeViewer(result::CalcSinglePlate; fig = nothing) -> ModalShapeViewer

Stacked modal-shape depth-panel viewer for an already-`calc()`'d `CalcSinglePlate` result --
port of `guliwerr/gui/ModalShapeViewerGUIPlate.m`, packaged as a degenerate 1-angle "sweep".

Throws `ArgumentError` if `result.ritzSolver.nCompModes == 0` (not yet `calc()`'d).
"""
function GuliwerrJ.modalShapeViewer(result::CalcSinglePlate; fig = nothing)
    result.ritzSolver.nCompModes == 0 && throw(ArgumentError(
        "modalShapeViewer: result has not been calc()'d yet (ritzSolver.nCompModes == 0) " *
        "-- nothing to view."))
    return _msvBuildPlate([result.angle], [result]; fig = fig)
end

"""
    modalShapeViewer(sweepResult::CalcAngSweepPlate; fig = nothing) -> ModalShapeViewer

Stacked modal-shape depth-panel viewer for an already-`calc()`'d `CalcAngSweepPlate` sweep --
port of `guliwerr/gui/ModalShapeViewerGUIPlate.m`.

Throws `ArgumentError` if `sweepResult.velsMap` is empty (not yet `calc()`'d -- `cs` is
populated by the CONSTRUCTOR itself with unsolved `CalcSinglePlate` objects, so `isempty(cs)`
would never be a valid "has calc() run" signal; matches MATLAB's own identical pitfall-avoidance
note).
"""
function GuliwerrJ.modalShapeViewer(sweepResult::CalcAngSweepPlate; fig = nothing)
    isempty(sweepResult.velsMap) && throw(ArgumentError(
        "modalShapeViewer: sweepResult has not been calc()'d yet (velsMap is empty) -- " *
        "nothing to view."))
    return _msvBuildPlate(sweepResult.angles, sweepResult.cs; fig = fig)
end

# -----------------------------------------------------------------------------------------
# _msvBuildIsofreq -- shared body for the isoFreq variant (CalcIsofreqSingleLayered/
# CalcIsofreqAngSweepLayered). Green's response is the ONLY source (no Ritz solve exists for
# this family at all) -- no sourceMenu, no Ritz cache, no modeBox; `kParRange` is a continuous,
# physically-meaningful axis, not a discrete Ritz mode index, so state is (currentAngNo,
# currentKParVal) rather than the base variant's (currentAngNo, currentModeNo/currentVel)
# dual-meaning pair -- port of `guliwerr/gui/IsofreqModalShapeViewerGUI.m`.
# -----------------------------------------------------------------------------------------

function _msvBuildIsofreq(angles::Vector{Float64}, structure::Structure, omega0::Real,
                           kParRange::Vector{Float64}, alreadyRotated::Bool; fig = nothing)
    fig, layout = guiHostLayout(fig; size = (560, 220))
    mainLayout = guiGridLayout(layout[1, 1])

    angleMenu, anglePrev, angleNext =
        guiDropdownStepperRow(mainLayout, 1, "Angle:", _msvAngleItems(angles))

    # << / >> step to the neighbouring kParRange sample; the trailing label shows the matching
    # phase velocity v = omega0/kParRange -- IsofreqModalShapeViewerGUI's kParRange row
    kParBox, kParPrev, kParNext, velLabel = guiIndexStepperRow(mainLayout, 2, "kParRange:")
    kParBox.validator[] = Float64
    kParBox.stored_string[] = string(kParRange[1])

    forceDirMenu, forceDirPrev, forceDirNext =
        guiDropdownStepperRow(mainLayout, 3, "Force dir:",
                               ["x (propag. dir.)", "y (in-plane)", "z (out-of-plane)"])
    statusLabel = guiStatusLabel(mainLayout[4, 1])

    shapeFig = Figure(size = (900, 500))

    state = ModalShapeViewer(fig, shapeFig, angles, nothing,
                              structure, Float64(omega0), kParRange, alreadyRotated,
                              angleMenu, (anglePrev, angleNext), nothing,
                              (kParPrev, kParNext), velLabel, kParBox,
                              nothing, forceDirMenu, nothing, statusLabel,
                              1, 1, NaN, kParRange[1],
                              Any[],
                              Ref{Union{Nothing,Function}}(nothing), Ref(false),
                              () -> nothing, nothing, nothing, nothing, nothing,
                              (_...) -> nothing)

    function _redraw!()
        angleVal = state.angles[state.currentAngNo]
        forceDir = state.forceDirMenu.i_selected[]
        velForGreen = state.isofreqOmega0 / state.currentKParVal
        angleForGreen = state.isofreqAlreadyRotated ? 0.0 : angleVal
        ok, errMsg = true, ""
        try
            greenObj = calc(GreenSolver, state.isofreqStructure, angleForGreen, [velForGreen];
                             forceDir = forceDir, displDir = 3, scanMode = :isoFreq,
                             omega0 = state.isofreqOmega0, computeMinus = false,
                             fullOutput = true)
            str = state.isofreqStructure
            if length(str.calcParams) > 1
                # layered (the usual isoFreq case): one row per layer, 'boosted' rows as the
                # Ritz viewers -- matches IsofreqModalShapeViewerGUI.redrawShape
                plotGreenModalShapeLayered(greenObj, str, forceDir, :boosted; fig = state.shapeFig)
            else
                plotGreenModalShape(greenObj, forceDir; fig = state.shapeFig)
            end
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end
        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _syncWidgets!()
        state.updatingWidgets[] = true
        try
            state.angleMenu.i_selected[] = state.currentAngNo
            state.kParBox.stored_string[] = string(state.currentKParVal)
            vel = state.isofreqOmega0 / state.currentKParVal
            state.velLabel.text[] = "$(round(vel, digits = 3)) km/s"
        finally
            state.updatingWidgets[] = false
        end
        return nothing
    end

    function _syncControlsAndRedraw!()
        _syncWidgets!()
        _redraw!()
        cb = state.onSelectionChanged[]
        cb === nothing || cb(state.angles[state.currentAngNo], state.currentKParVal)
        return nothing
    end

    function _goToPoint!(angleVal::Real, kParVal::Real)
        angNo = _msvNearestAngle(state.angles, angleVal)
        state.currentAngNo = angNo
        state.currentKParVal = Float64(kParVal)
        _syncControlsAndRedraw!()
        return nothing
    end
    state.goToPoint! = _goToPoint!

    function _stepToAngle!(newAngNo::Integer)
        # keeps the SAME kParRange across the step -- direct analogue of the base variant's
        # own Green's-response angle-stepping behaviour (no discrete mode/type to follow for a
        # continuous value), matching MATLAB exactly
        state.goToPoint!(state.angles[newAngNo], state.currentKParVal)
        return nothing
    end
    state.stepToAngle! = _stepToAngle!

    # neighbouring kParRange sample, clamped at the ends -- IsofreqModalShapeViewerGUI.stepKPar
    function _stepKPar!(direction::Integer)
        grid = sort(state.isofreqKParRange)
        tol = 1e-9 * maximum(abs, grid)
        k = direction > 0 ? something(findfirst(>(state.currentKParVal + tol), grid), length(grid)) :
                            something(findlast(<(state.currentKParVal - tol), grid), 1)
        state.goToPoint!(state.angles[state.currentAngNo], grid[k])
        return nothing
    end

    on(idx -> (state.updatingWidgets[] || state.stepToAngle!(idx)), angleMenu.i_selected)
    on(idx -> (state.updatingWidgets[] || _redraw!()), forceDirMenu.i_selected)
    on(_ -> _stepKPar!(-1), kParPrev.clicks)
    on(_ -> _stepKPar!(+1), kParNext.clicks)
    on(kParBox.stored_string) do s
        state.updatingWidgets[] && return nothing
        v = tryparse(Float64, s)
        v === nothing && return nothing
        state.goToPoint!(state.angles[state.currentAngNo], v)
        return nothing
    end
    on(_ -> state.stepToAngle!(max(state.currentAngNo - 1, 1)), anglePrev.clicks)
    on(_ -> state.stepToAngle!(min(state.currentAngNo + 1, length(state.angles))), angleNext.clicks)
    _msvWireMenuStepper!(forceDirMenu, forceDirPrev, forceDirNext, 3)

    state.goToPoint!(angles[1], kParRange[1])
    return state
end

"""
    modalShapeViewer(result::CalcIsofreqSingleLayered; fig = nothing) -> ModalShapeViewer

On-demand Green's-response depth-resolved displacement-field viewer for an already-`calc()`'d
`CalcIsofreqSingleLayered` result -- port of `guliwerr/gui/IsofreqModalShapeViewerGUI.m`,
packaged as a degenerate 1-angle "sweep" with `isofreqAlreadyRotated = true` (`result.structure`
is already rotated to `result.angle` -- see this file's header note on
`ModalShapeViewer.isofreqStructure`/the rotation convention).

A layered structure (the usual isoFreq case) is drawn with `plotGreenModalShapeLayered`, one
row per layer; a half-space with `plotGreenModalShape`. Depths are already in units of the
selected wave's own wavelength (`GreenSolver`'s `hScaleOut = omega0/v`), so nothing
isoFreq-specific is needed (see `IsofreqModalShapeViewerGUI.m`'s own docblock).

Throws `ArgumentError` if `result.absG` is empty (not yet `calc()`'d).
"""
function GuliwerrJ.modalShapeViewer(result::CalcIsofreqSingleLayered; fig = nothing)
    isempty(result.absG) && throw(ArgumentError(
        "modalShapeViewer: result has not been calc()'d yet (absG is empty) -- nothing to view."))
    return _msvBuildIsofreq([result.angle], result.structure, result.omega0, result.kParRange,
                            true; fig = fig)
end

"""
    modalShapeViewer(sweepResult::CalcIsofreqAngSweepLayered; fig = nothing) -> ModalShapeViewer

On-demand Green's-response depth-resolved displacement-field viewer for an already-`calc()`'d
`CalcIsofreqAngSweepLayered` sweep -- port of `guliwerr/gui/IsofreqModalShapeViewerGUI.m`, with
`isofreqAlreadyRotated = false` (`sweepResult.structure` is the UNROTATED base; `redraw!`
rotates it per angle itself via `GreenSolver.calc`'s own `angle` argument).

Throws `ArgumentError` if `sweepResult.absGMap` is empty (not yet `calc()`'d).
"""
function GuliwerrJ.modalShapeViewer(sweepResult::CalcIsofreqAngSweepLayered; fig = nothing)
    isempty(sweepResult.absGMap) && throw(ArgumentError(
        "modalShapeViewer: sweepResult has not been calc()'d yet (absGMap is empty) -- " *
        "nothing to view."))
    return _msvBuildIsofreq(sweepResult.angles, sweepResult.structure, sweepResult.omega0,
                            sweepResult.kParRange, false; fig = fig)
end
