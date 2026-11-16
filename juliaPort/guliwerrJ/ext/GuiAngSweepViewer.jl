# =========================================================================================
# GuiAngSweepViewer.jl -- viewer for angle-sweep maps ("2D map + click-to-select + drill-down")
#
# Julia/Makie port of `guliwerr/gui/AngSweepViewerGUI.m`, `AngSweepViewerGUIPlate.m` and
# `IsofreqAngSweepViewerGUI.m`. The three MATLAB classes are near-duplicates (identical
# ~11-method surface; Plate drops the 3-way graph-type dropdown for one plot function and
# swaps the criteria vocabulary), the same shape as the base/Plate/isoFreq split of
# `ext/GuiSingleAngleViewer.jl`. `include`d directly into `module GuliwerrJMakieExt` (see
# that file's own header), which already has `using GuliwerrJ`/`using Makie` in scope --
# neither repeated here. Built on `GuiWidgets.jl` (`guiHostLayout`, `guiPlotWithSidePanel`,
# `guiDropdownStepperRow`, `guiRangeFieldPair`, `guiBindRangeToAxis!`, `guiOnAxisClick`,
# `guiNearestPoint`, `guiCrosshair!`, `guiRebuildTypeToggles!`, `guiStatusLabel`/
# `guiSetStatusOk!`/`guiSetStatusError!`) AND on `GuiSingleAngleViewer.jl`'s
# `singleAngleViewer`/`SingleAngleViewer`/`guiSignificantModeType`/`guiPlateModeTypeText` --
# this file MUST be `include`d AFTER `GuiSingleAngleViewer.jl` (as `GuliwerrJMakieExt.jl` does).
#
# ---------------------------------------------------------------------------------------
# THIS IS A VIEWER, NOT A CALCULATOR
# ---------------------------------------------------------------------------------------
# All three MATLAB constructors take an ALREADY-`calc()`'d sweep result and only ERROR if it
# looks unsolved (`isempty(sweepResult.significantModesMap)` / `isempty(sweepResult.velsMap)` /
# `isempty(sweepResult.absGMap)`) -- no solver call anywhere in `buildFigure`/the constructor
# itself. `angSweepViewer` below does exactly the same. The Ritz spectrum map is rebuilt
# from the stored modal numerators, and the map crosshair only needs the clicked point's own
# coordinates, so nothing is cached for it.
#
# ---------------------------------------------------------------------------------------
# DISPATCH DESIGN -- ONE struct, multiple-dispatch `angSweepViewer` methods
# ---------------------------------------------------------------------------------------
# `AngSweepViewer` (below) is the ONE handle type returned by all three variants --
# `angSweepViewer(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}; fig)`,
# `angSweepViewer(sweepResult::CalcAngSweepPlate; fig)`,
# `angSweepViewer(sweepResult::CalcIsofreqAngSweepLayered; fig)` -- the same
# "one struct, several dispatch methods sharing the shell" design as `GuiSingleAngleViewer.jl`.
# Each method builds the SAME shell (host layout, top bar, stable content slot, status line)
# and differs only at the genuine points of difference: which plotting function draws the
# map, the criterion/scale/plot-type VOCABULARY, the default y-range SOURCE, and the
# info-panel TEXT/click-resolution semantics (nearest mode by velocity vs. nearest `kParRange`
# sample).
#
# ---------------------------------------------------------------------------------------
# CLICK-TO-SELECT -- a NESTED two-stage 1D search, not `guiNearestPoint`'s joint 2D search
# ---------------------------------------------------------------------------------------
# `GuiWidgets.jl`'s own header calls a 2D axis-range-normalised search across several series
# "the map viewer's use case" for `guiNearestPoint` -- but every one of this trio's
# `onMapClicked` methods resolves a click SEQUENTIALLY, not jointly: first the nearest ANGLE
# alone (`[~, angNo] = min(abs(angles - angPos))`, ignoring velocity/kParRange entirely), THEN
# the nearest mode-velocity (or `kParRange` sample) AT THAT ONE ANGLE alone. A joint 2D
# nearest-neighbour search over the whole flattened (angle, velocity) point cloud would give a
# DIFFERENT answer whenever a nearby angle happens to have a mode closer in y than the
# genuinely nearest angle's own mode -- not what MATLAB computes. So this file reproduces the
# two independent 1D `[~, idx] = min(abs(x - q))` searches literally (`_aswNearestIndex`
# below). `GuiSingleAngleViewer.jl` gets a 1D search out of `guiNearestPoint` by fixing
# `y = 0`, but here TWO separate axes are searched, angle and velocity/kParRange, never
# together, so `guiNearestPoint`'s 2D form is the wrong tool for BOTH stages.
#
# ---------------------------------------------------------------------------------------
# DRILL-DOWN -- a REAL embedded child viewer, via `singleAngleViewer`'s own `fig` keyword
# ---------------------------------------------------------------------------------------
# All three MATLAB classes DO drill down into a genuine child viewer window: a "View single
# angle ->" button (`onToggleSingleAngleView`/`viewSingleAngleBtn`) swaps the map content for
# an embedded `SingleAngleViewerGUI`/`...Plate`/`IsofreqSingleLayeredViewerGUI` at
# `obj.currentAngNo` (the last-clicked angle, or the first if none was clicked) -- NOT
# triggered directly by a map click itself, which only sets `currentAngNo` and fires
# `OnModeSelected`/`OnPointSelected`. `singleAngleViewer` accepts `fig` for embedding into a
# caller-supplied `GridPosition`, so `angSweepViewer`'s own "View single angle ->" button
# embeds an actual `SingleAngleViewer` at the content slot, dispatching to the matching
# `singleAngleViewer` method automatically via `sweepResult.cs[angNo]`'s own concrete type.
# This gives full behavioural parity with MATLAB.
#
# `GuiWidgets.jl`'s own `guiRebuildTypeToggles!` clears a panel via `delete!(gc.content)` for
# each `gc in panel.content`, which works ONLY because a toggle panel holds `Block`s directly
# (`Toggle`/`Label`, and `Makie.Block` has its own `delete!` method). This viewer's content
# slot instead holds NESTED `GridLayout`s (`guiPlotWithSidePanel`'s own `outer`, or
# `singleAngleViewer`'s `guiHostLayout`-built `layout`), and `delete!` has NO method for
# `GridLayout` (`methods(delete!, (GridLayout,))` returns zero methods), so reusing that idiom
# here throws `MethodError`. The lower-level
# `Makie.GridLayoutBase.remove_from_gridlayout!(gc::GridContent)` -- the very function
# `delete!(::Block)` itself calls internally -- removes ANY grid content, `Block` or nested
# `GridLayout` alike, with no stacking or leftover content across repeated map -> drill-down
# -> map cycles. `guiClearGridLayout!` below wraps this; it is `gui`-prefixed (not
# `_asw`-private) so any other viewer that swaps a region's content wholesale and repeatedly
# can reuse it, as every viewer file already reuses `guiNearestPoint`/`guiCrosshair!`.
#
# NOT REPRODUCED: MATLAB's public `showCrosshairAt` (the hand-off point a linked
# `ModalShapeViewerGUI` attaches to, so a selection made in THAT window can mark this map's
# crosshair without re-triggering click/mode-lookup logic, re-drilling first if the linked
# angle differs from `currentAngNo`). The OUTGOING half of the same hand-off (`onSelect`,
# this viewer telling a listener what was clicked) IS implemented below. The INCOMING half
# is not a method of `AngSweepViewer`: `GuiCombinedViewer.jl` implements it from outside,
# through the public `crosshair`, `currentAngNo` and `toggleDrillDown!` fields.
#
# ---------------------------------------------------------------------------------------
# COLOUR-RANGE CONTROL
# ---------------------------------------------------------------------------------------
# `guiBindRangeToAxis!`'s own docstring says the colour-range case is "deliberately NOT
# covered here... wire that one directly at the call site instead" -- this file is that call
# site (a 1D curve viewer has no colour axis to control). Makie has no MATLAB-style single
# axes-wide `CLim`: colour range lives on the PLOT OBJECT (`Heatmap`/`Scatter`), and a map
# plotting function returns only the `Axis`, never that object -- so `guiColorPlotObject`
# (below) walks `ax.scene.plots` (the same "read back an unreturned handle via its live
# attributes" technique `GuiSingleAngleViewer.jl`'s own header documents for `label`)
# for the first `Heatmap`/`Scatter` carrying a GENUINELY data-driven `colorrange`. A plot
# drawn with a FIXED colour (e.g. `_diagScatterByType!`'s one-RGB-per-type scatters, or
# `plotRitzModeMapPlate`'s `criterion = :modeType` scatters) reports
# `colorrange[] isa Makie.Automatic`, while one drawn with a genuine `color = <data array>`
# (every Heatmap; `plotRitzModeMap`'s/`plotRitzModeMapPlate`'s `:confidence` scatter) reports
# a concrete 2-element value -- so `Makie.Automatic` reliably distinguishes "nothing to
# control" from "the one object Colour min/max should drive," with no per-graph-type
# special-casing needed. Every graph type this viewer draws puts its ONE genuinely
# colour-mapped object FIRST when more than one plot object exists (the green background
# heatmap, when `plotClassificationDiagnostic` shows one, is drawn before its own
# fixed-colour per-type scatters), so "first match" is reliable.
#
# `guiRebindColorRange!` (below) re-derives the Colour min/max fields from that object's own
# CURRENT `colorrange` on every redraw (MATLAB's own "re-default from the just-drawn CLim"
# mechanic, this class's own docblock: one shared field pair stays safe across graph types
# with fundamentally different colour scales by re-defaulting on every switch, never by
# picking one fixed range) and (re)binds a LIVE listener so a later user edit writes straight
# into `plotobj.colorrange[]`. Because the underlying plot object is destroyed and rebuilt by
# every redraw (`empty!(ax)`), the previous listener is explicitly `off()`'d first (the
# handle `on(f, obs)` returns stops future calls), so listeners never accumulate across
# repeated redraws pointing at dead plot objects. Both functions are `gui`-prefixed (not
# `_asw`-private): any other viewer drawing its own colour-mapped map (the dispersion maps,
# say) can reuse this exact pair.
#
# ---------------------------------------------------------------------------------------
# DELIBERATE DEVIATIONS FROM THE MATLAB SOURCE (documented, not oversights)
# ---------------------------------------------------------------------------------------
#   1. `<<`/`>>` STEPPER BUTTONS on every dropdown -- same added convenience as
#      `GuiSingleAngleViewer.jl`'s deviation 1, no MATLAB counterpart to contradict.
#   2. GRAPH-TYPE/CRITERION/SCALE/PLOT-TYPE DROPDOWNS ARE NEVER DISABLED WHILE DRILLED IN,
#      unlike MATLAB's `obj.graphTypeDropdown.Enable = 'off'`/`obj.criterionDropdown.Enable =
#      'off'`. `Makie.Menu`/`Textbox`/`Button` expose no `enabled`/`visible` attribute in
#      this Makie version at all, only `Label` and plot objects do (same finding
#      `GuiStructureControl.jl`'s own header records). Harmless rather than
#      misleading: `redraw!` returns immediately whenever `ax === nothing` (drilled in), so
#      switching a dropdown while drilled in silently has no effect until the user returns to
#      the map view -- discoverable (the map is right there when you go back), not actively
#      wrong.
#   3. IsoFreq's `plotType = :polar`: the map is drawn on a plain `Axis` in pol2cart X/Y
#      (`plotPolarHeatmap`), so -- as in MATLAB -- a click is converted back to (angle, kPar),
#      the crosshair is placed at pol2cart, and the angle/kParRange range fields do not set
#      the axis limits (the colour range still applies). MATLAB also disables those fields;
#      here they stay editable and take effect back in the Cartesian view.
#   4. NO `showCrosshairAt` METHOD -- see the DRILL-DOWN section above.
#   5. BASE VARIANT'S "0"/"10" HARDCODED VELOCITY-RANGE FALLBACK IS PORTED FAITHFULLY, even
#      though it differs from `GuiSingleAngleViewer.jl`'s more principled Ritz-spectrum-derived
#      fallback for the analogous single-angle case (`_savDefaultXRange`). It follows
#      `AngSweepViewerGUI.buildMapContent`'s own comment ("Falls back to 0/10 (the old
#      hardcoded default) only when [greenVelAxis is] empty") -- a minor MATLAB
#      inconsistency between two sibling viewer classes, not a Julia translation choice;
#      ported as-is (`_aswDefaultYRange` below) rather than silently harmonised.
#   6. TYPE-TOGGLE PANEL reuses `GuiSingleAngleViewer.jl`'s technique (`_diagScatterByType!`
#      draws each classification type with `label = String(t)`; `ax.scene.plots` is filtered
#      for `p isa Makie.Scatter && haskey(p, :label)`) verbatim -- only the base variant's
#      Classification-diagnostic graph type has one (Ritz-mode/Green's-function maps colour
#      continuously; Plate/isoFreq have no per-type legend at all, matching each MATLAB
#      class's own docblock).
#   7. EXPLICIT AXIS-RANGE RE-APPLICATION after every redraw (`_aswReapplyRange!`) -- the same
#      technique and the same underlying Makie behaviour (`ax.limits` persists across
#      `empty!(ax)` + a fresh plot UNLESS the new plot call sets it again itself)
#      `GuiSingleAngleViewer.jl`'s own header documents;
#      `plotDispersionMap`'s own `climRange`-driven `ylims!` under the isoFreq variant's
#      `:log` scale mode is the one plotting call in THIS file's scope that can trigger it.
#   8. PRIVATE HELPER DUPLICATION (`_asw`-prefixed private to this file; `_savXxx`-prefixed
#      ones in `GuiSingleAngleViewer.jl` are NOT called from here even where the underlying
#      arithmetic is similar) -- same reasoning as `GuiSingleAngleViewer.jl`'s deviation 5:
#      the two viewer files share the substrate but do not depend on each other's private
#      internals. `guiClearGridLayout!`/`guiColorPlotObject`/`guiRebindColorRange!` are the
#      exception, `gui`-prefixed and genuinely meant for reuse -- see the sections above.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# Shared, reusable helpers (`gui`-prefixed -- see this file's header)
# -----------------------------------------------------------------------------------------

"""
    guiClearGridLayout!(gl::GridLayout) -> Nothing

Wholesale-remove every child of `gl` -- `Block`s (`Axis`, `Menu`, `Textbox`, `Label`, ...) AND
nested `GridLayout`s alike -- so a fresh call building either the map content or an embedded
drill-down viewer at `gl[1, 1]` never stacks on top of the previous content. See this file's
header ("DRILL-DOWN") for why `GuiWidgets.jl`'s own `guiRebuildTypeToggles!` idiom
(`delete!(gc.content)`) does not generalise to a slot holding nested `GridLayout`s on its
own. Removing a child from the layout (`Makie.GridLayoutBase.remove_from_gridlayout!`) is not
enough either: the `Block`s inside stay in the figure and its scene, still drawn (and still
clickable) under the new content -- 5 slice steps of a frequency-sweep viewer left 7 `Axis`,
7 `Colorbar` and 65 `Button`s alive. So nested layouts are cleared
recursively and every `Block` is `delete!`d, which removes it from layout, figure and scene.
An `Axis`' `axislegend`s are not layout content, so they are deleted with it
([`guiEmptyAxis!`](@ref)).
"""
function guiClearGridLayout!(gl::GridLayout)
    for gc in copy(gl.content)
        c = gc.content
        if c isa Makie.Block
            c isa Axis && guiEmptyAxis!(c)   # its axislegends sit outside the layout
            delete!(c)
        else
            c isa GridLayout && guiClearGridLayout!(c)
            Makie.GridLayoutBase.remove_from_gridlayout!(gc)
        end
    end
    return nothing
end

"""
    guiColorPlotObject(ax::Axis) -> Union{Nothing,Any}

The first `Heatmap`/`Scatter` plot object on `ax` carrying a genuinely DATA-DRIVEN
`colorrange` -- see this file's header ("COLOUR-RANGE CONTROL") for the
`Makie.Automatic()` behaviour this relies on. `nothing` when nothing colour-mapped was drawn at
all (e.g. `plotRitzModeMapPlate` under `criterion = :modeType`, or a classification diagnostic
drawn with no Green background to show).
"""
function guiColorPlotObject(ax::Axis)
    for p in ax.scene.plots
        (p isa Makie.Heatmap || p isa Makie.Scatter) || continue
        cr = p.colorrange[]
        cr isa Makie.Automatic && continue
        return p
    end
    return nothing
end

"""
    guiRebindColorRange!(ax::Axis, cRangeField::GuiRangeField,
                          cRangeBoxes::Tuple{Textbox,Textbox},
                          colorListener::Base.RefValue) -> Nothing

Re-derive `cRangeBoxes`' text from `ax`'s current colour-mapped object
([`guiColorPlotObject`](@ref)) and (re)bind a live listener so a later edit of `cRangeField`
writes straight into that object's own `colorrange`. Call once per redraw, AFTER the redraw's
own plotting call and BEFORE the crosshair is (re)created (the crosshair's own fixed-colour
marker would never match `guiColorPlotObject` anyway, but calling before keeps the "first
colour-mapped object on this axis" search scoped to what the redraw itself just drew).

`colorListener[]`, if non-`nothing`, is `off()`'d first -- the previous colour-mapped object
was destroyed by the redraw's own `empty!(ax)`, so its listener must not be left dangling
against a rebuild's worth of accumulated dead objects. When nothing colour-mapped is found,
`cRangeBoxes`' text is cleared (blank) and no listener is installed -- editing the fields then
has no effect until a redraw draws something colour-mapped again, matching MATLAB's own
`cMinField.Enable = false` intent without an actual disable primitive to use (see this file's
header, deviation 2).
"""
function guiRebindColorRange!(ax::Axis, cRangeField::GuiRangeField,
                               cRangeBoxes::Tuple{Textbox,Textbox}, colorListener::Base.RefValue)
    colorListener[] === nothing || off(colorListener[])
    colorListener[] = nothing
    loBox, hiBox = cRangeBoxes
    p = guiColorPlotObject(ax)
    if p === nothing
        loBox.stored_string[] = ""
        hiBox.stored_string[] = ""
        return nothing
    end
    cr = p.colorrange[]
    loBox.stored_string[] = string(round(Float64(cr[1]); sigdigits = 4))
    hiBox.stored_string[] = string(round(Float64(cr[2]); sigdigits = 4))
    colorListener[] = on(cRangeField.range) do r
        r === nothing && return
        p.colorrange[] = (Float32(r[1]), Float32(r[2]))
    end
    return nothing
end

# -----------------------------------------------------------------------------------------
# Private pure helpers (`_asw`-prefixed -- private to this file, see deviation 8 above)
# -----------------------------------------------------------------------------------------

"""
    _aswNearestIndex(xs::AbstractVector{<:Real}, xQuery::Real) -> Union{Int,Nothing}

Index into `xs` nearest `xQuery` by plain absolute difference -- MATLAB's own literal
`[~, idx] = min(abs(xs - xQuery))`. See this file's header ("CLICK-TO-SELECT") for why this is
NOT routed through `guiNearestPoint`: every `onMapClicked` in this trio resolves a click as
TWO INDEPENDENT sequential 1D searches (angle, then mode-velocity/kParRange AT that angle),
never a single joint 2D nearest-point search. `nothing` only when `xs` is empty.
"""
function _aswNearestIndex(xs::AbstractVector{<:Real}, xQuery::Real)
    isempty(xs) && return nothing
    return argmin(abs.(xs .- xQuery))
end

"""
    _aswReapplyRange!(ax::Axis, rf::GuiRangeField, kind::Symbol) -> Nothing

Re-apply `rf`'s CURRENT value to `ax` right now, unconditionally -- duplicates
`GuiSingleAngleViewer.jl`'s own `_savReapplyRange!` verbatim (see deviation 7/8 above for why
this is a deliberate duplication, not a missed opportunity to share).
"""
function _aswReapplyRange!(ax, rf::GuiRangeField, kind::Symbol)
    r = rf.range[]
    r === nothing && return nothing
    kind === :x ? xlims!(ax, r) : ylims!(ax, r)
    return nothing
end

"""
    _aswWireStepper!(menu::Menu, prevBtn::Button, nextBtn::Button, n::Integer) -> Nothing

Wire a dropdown-stepper row's `<<`/`>>` buttons to cycle `menu.i_selected` through `1:n`,
wrapping -- duplicates `GuiSingleAngleViewer.jl`'s own `_savWireStepper!` verbatim.
"""
function _aswWireStepper!(menu, prevBtn, nextBtn, n::Integer)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] - 1, n)), prevBtn.clicks)
    on(_ -> (menu.i_selected[] = mod1(menu.i_selected[] + 1, n)), nextBtn.clicks)
    return nothing
end

# -- default y-range source, dispatched per variant --

"""
    _aswDefaultYRange(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}) -> (String,String)

Default velocity-axis text-field values: `greenVelAxis`'s own min/max under `calcMode = :full`,
falling back to the LITERAL "0"/"10" MATLAB hardcodes -- port of
`AngSweepViewerGUI.buildMapContent`'s own fallback. See this file's header, deviation 5, for
why this is NOT the same fallback `GuiSingleAngleViewer.jl`'s `_savDefaultXRange` uses for the
single-angle case (a genuine, minor MATLAB inconsistency, ported faithfully).
"""
function _aswDefaultYRange(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered})
    if !isempty(sweepResult.greenVelAxis)
        return string(minimum(sweepResult.greenVelAxis)), string(maximum(sweepResult.greenVelAxis))
    end
    return "0", "10"
end

"""
    _aswDefaultYRange(sweepResult::CalcAngSweepPlate) -> (String,String)

Blank: the velocity fields are filled from the drawn map on the first redraw -- port of
`AngSweepViewerGUIPlate.buildMapContent`/`onGraphTypeChanged`. No Green's function/
`greenVelAxis` exists for a free plate, and the whole Ritz spectrum reaches hundreds of km/s,
far above [`plotRitzModeMapPlate`](@ref)'s own velocity cutoff.
"""
_aswDefaultYRange(::CalcAngSweepPlate) = ("", "")

"""
    _aswDefaultYRange(sweepResult::CalcIsofreqAngSweepLayered) -> (String,String)

Default `kParRange`-axis text-field values: the sweep's own shared `kParRange` min/max (the
SAME axis at every angle) -- port of `IsofreqAngSweepViewerGUI.buildMapContent`.
"""
function _aswDefaultYRange(sweepResult::CalcIsofreqAngSweepLayered)
    return string(minimum(sweepResult.kParRange)), string(maximum(sweepResult.kParRange))
end

# -- info-panel text, dispatched per variant --

"""
    _aswModeInfoLines(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}, angNo::Integer,
                       modeNo::Integer, angleVal::Real, velSel::Real) -> Vector{String}

Four-line info-panel text for one selected (angle, Ritz mode) -- port of
`AngSweepViewerGUI.onMapClicked`'s `obj.infoPanel.Text` assignment. `type` comes from
[`guiSignificantModeType`](@ref) against `sweepResult.significantModesMap[angNo]`.
"""
function _aswModeInfoLines(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}, angNo::Integer,
                            modeNo::Integer, angleVal::Real, velSel::Real)
    typeStr = string(guiSignificantModeType(sweepResult.significantModesMap[angNo], modeNo))
    return ["angle = $(round(angleVal, digits = 1)) deg",
            "mode # = $modeNo",
            "velocity = $(round(velSel, digits = 4)) km/s",
            "type = $typeStr"]
end

"""
    _aswModeInfoLines(sweepResult::CalcAngSweepPlate, angNo::Integer, modeNo::Integer,
                       angleVal::Real, velSel::Real) -> Vector{String}

Plate counterpart -- port of `AngSweepViewerGUIPlate.onMapClicked`. `type` comes from
[`guiPlateModeTypeText`](@ref) (mode type + confidence), not a significant-modes lookup.
"""
function _aswModeInfoLines(sweepResult::CalcAngSweepPlate, angNo::Integer, modeNo::Integer,
                            angleVal::Real, velSel::Real)
    ritz = sweepResult.cs[angNo].ritzSolver
    return ["angle = $(round(angleVal, digits = 1)) deg",
            "mode # = $modeNo",
            "velocity = $(round(velSel, digits = 4)) km/s",
            "type = $(guiPlateModeTypeText(ritz, modeNo))"]
end

"""
    _aswPointInfoLines(sweepResult::CalcIsofreqAngSweepLayered, angNo::Integer, kIdx::Integer,
                        angleVal::Real, kParVal::Real) -> Vector{String}

THREE-line info-panel text for one selected (angle, `kParRange` sample) -- port of
`IsofreqAngSweepViewerGUI.onMapClicked`'s `obj.infoPanel.Text` assignment. Deliberately only
THREE lines (angle, kParRange, |G_33|), NOT the four-line `angle/kParRange/velocity/|G_33|`
form `GuiSingleAngleViewer.jl`'s `_savPointInfoLines` uses for the analogous single-angle isoFreq
viewer -- MATLAB's own `IsofreqAngSweepViewerGUI.onMapClicked` genuinely omits the derived
velocity line here; ported as written, not harmonised with its single-angle sibling.
"""
function _aswPointInfoLines(sweepResult::CalcIsofreqAngSweepLayered, angNo::Integer, kIdx::Integer,
                             angleVal::Real, kParVal::Real)
    absGVal = sweepResult.absGMap[kIdx, angNo]
    return ["angle = $(round(angleVal, digits = 1)) deg",
            "kParRange = $(round(kParVal, sigdigits = 4))",
            "|G_{33}| = $(round(absGVal, sigdigits = 4))"]
end

# -----------------------------------------------------------------------------------------
# AngSweepViewer -- the one handle type returned by every variant
# -----------------------------------------------------------------------------------------

"""
    AngSweepViewer

Handles for one "2D map + click-to-select + drill-down" viewer window, returned by
every [`angSweepViewer`](@ref) method regardless of which `CalcAngSweepHS`/`CalcAngSweepLayered`/
`CalcAngSweepPlate`/`CalcIsofreqAngSweepLayered` variant built it (see this file's header for
the dispatch design).

Fields:
- `fig` -- the owning `Figure`.
- `ax` -- the current map `Axis`, or `nothing` while drilled into a single-angle view (the map
  content, including `ax`, is destroyed by [`guiClearGridLayout!`](@ref) during drill-down --
  there is no MATLAB-style stale-but-still-referenced handle to guard against here).
- `result` -- the sweep result this viewer displays (untyped field, kept generic so this ONE
  struct serves every variant -- same reasoning as `SingleAngleViewer.result`).
- `contentGL` -- the stable content slot (`GridLayout`) the map/drill-down content is (re)built
  into -- itself never destroyed, mirroring MATLAB's `obj.contentParent`.
- `graphTypeMenu` -- the graph-type dropdown for the base (HS/Layered) variant, or `nothing`
  for Plate/isoFreq (absent, not disabled -- matching MATLAB's own precedent for those two
  classes, same as `SingleAngleViewer`).
- `modeMenu` -- the variant's criterion/scale-mode dropdown: criterion (base/Plate) or scale
  mode `log`/`linear` (isoFreq). Never `nothing`.
- `plotTypeMenu` -- the isoFreq variant's `cartesian`/`polar` dropdown, or `nothing` for
  base/Plate (neither exposes a plot-type choice at all).
- `xRangeField`/`yRangeField`/`cRangeField` -- persistent axis-range state (`GuiRangeField`);
  `xRangeBoxes`/`yRangeBoxes`/`cRangeBoxes` -- the CURRENTLY MOUNTED `(loBox, hiBox)` `Textbox`
  pairs, `nothing` while drilled in (destroyed along with the map content; the underlying
  `GuiRangeField`s themselves are never destroyed -- see `guiRangeFieldPair`'s own docstring:
  MATLAB's `GuiCommon.m` captureAxisRange/restoreAxisRange needed no Julia port at all --
  passing the SAME field back in as
  `state` on the next `buildMapContent!` reproduces MATLAB's `captureAxisRange`/
  `restoreAxisRange` round trip for free).
- `statusLabel` -- ok/error status line; persists across drill-down (part of the outer
  layout, not `contentGL`).
- `infoLabel` -- the mode/point info panel, `nothing` while drilled in (the embedded child
  viewer has its own).
- `typeTogglePanel` -- the "Show types:" panel, populated only for the base variant's
  Classification-diagnostic graph type; `nothing` for Plate/isoFreq (no per-type legend) AND
  while drilled in.
- `viewSingleAngleBtn` -- the drill-down toggle button; its own label text flips between
  `"View single angle ->"`/`"<< Back to angular sweep"`, matching MATLAB.
- `crosshair` -- the current crosshair position `Observable{Point2f}`, or `nothing` while
  drilled in.
- `colorListener` -- the currently-installed colour-range listener handle (`Ref{Any}`, see
  [`guiRebindColorRange!`](@ref)), or `nothing`.
- `currentAngNo` -- last-clicked (or default first) angle index (`Ref{Int}`) -- the
  drill-down target `viewSingleAngleBtn` uses next.
- `isDrilledIn` -- `Ref{Bool}`, readable externally the same reason MATLAB's own
  `isDrilledIn` is `SetAccess = private` but still PUBLICLY READABLE (a composite host could
  need to know before touching map-only fields).
- `childViewer` -- the embedded `SingleAngleViewer` while drilled in, else `nothing`
  (`Ref{Any}`, since its concrete type varies by variant exactly the way `result` does).
- `onSelect` -- caller-settable callback (`Ref{Union{Nothing,Function}}`, `nothing` initially),
  MATLAB's `OnModeSelected`/`OnPointSelected`. Called after every successful click OR after a
  click inside the drilled-down child viewer is forwarded through (MATLAB's own
  `dispatchModeSelected`/`dispatchPointSelected`): `(angle, modeNo, velSnapped, velRaw)` for
  the base/Plate variants, `(angle, kParVal)` for isoFreq.
- `redraw!` -- `() -> Nothing`, redraws the map for whichever dropdown selection is current;
  a no-op while drilled in (`ax === nothing`). Also callable directly (e.g. from a test).
- `onClicked!` -- `(x::Float64, y::Float64) -> Nothing`, the click handler `guiOnAxisClick` is
  wired to on the CURRENT map `ax` (re-wired fresh every time `buildMapContent!` runs); also
  callable directly with a known data-space point, the same pattern as `SingleAngleViewer.onClicked!`.
- `toggleDrillDown!` -- `() -> Nothing`, the "View single angle ->"/"<< Back" button's own
  handler; also callable directly (e.g. from a test) without a simulated click.
"""
mutable struct AngSweepViewer
    fig::Figure
    ax::Union{Axis,Nothing}
    result::Any
    contentGL::GridLayout
    graphTypeMenu::Union{Menu,Nothing}
    modeMenu::Menu
    plotTypeMenu::Union{Menu,Nothing}
    xRangeField::GuiRangeField
    yRangeField::GuiRangeField
    cRangeField::GuiRangeField
    xRangeBoxes::Union{Tuple{Textbox,Textbox},Nothing}
    yRangeBoxes::Union{Tuple{Textbox,Textbox},Nothing}
    cRangeBoxes::Union{Tuple{Textbox,Textbox},Nothing}
    statusLabel::Label
    infoLabel::Union{Label,Nothing}
    typeTogglePanel::Union{GridLayout,Nothing}
    viewSingleAngleBtn::Button
    crosshair::Union{Observable{Point2f},Nothing}
    colorListener::Base.RefValue{Any}
    currentAngNo::Base.RefValue{Int}
    isDrilledIn::Base.RefValue{Bool}
    childViewer::Base.RefValue{Any}
    onSelect::Base.RefValue{Union{Nothing,Function}}
    redraw!::Function
    onClicked!::Function
    toggleDrillDown!::Function
end

# -----------------------------------------------------------------------------------------
# angSweepViewer -- base variant (CalcAngSweepHS / CalcAngSweepLayered)
# -----------------------------------------------------------------------------------------

"""
    angSweepViewer(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}; fig = nothing)
        -> AngSweepViewer

Angular-sweep map viewer for an already-`calc()`'d `CalcAngSweepHS`/`CalcAngSweepLayered`
result -- port of `guliwerr/gui/AngSweepViewerGUI.m`. `fig` -- `nothing` (default) for a
standalone window, or a `GridPosition` to embed into.

Four graph types in MATLAB's order (`graphTypeMenu`, default "Classification diagnostic"):
"Classification diagnostic" ([`plotClassificationDiagnostic`](@ref): G13_eff background, or
the Ritz spectrum for a `:ritzOnly` sweep), "Green's-function map" ([`plotGreenMap`](@ref)),
"Ritz-mode map" ([`plotRitzModeMap`](@ref)), "Ritz spectrum map"
([`plotRitzSpectrumMap`](@ref), stored numerators, no solve). `modeMenu` is the Ritz-mode
criterion dropdown (`ritzScore`/`freqDisp`/`depthInsensitivity`), consulted only while
"Ritz-mode map" is selected.

Throws `ArgumentError` if `sweepResult.significantModesMap` is empty (not yet `calc()`'d).
"""
function GuliwerrJ.angSweepViewer(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered}; fig = nothing)
    isempty(sweepResult.significantModesMap) && throw(ArgumentError(
        "angSweepViewer: sweepResult has not been calc()'d yet (significantModesMap is empty) " *
        "-- nothing to view."))

    fig, mainLayout = guiHostLayout(fig; size = (1000, 700))
    topBarGL = guiGridLayout(mainLayout[1, 1])
    contentGL = guiGridLayout(mainLayout[2, 1])
    statusLabel = guiStatusLabel(mainLayout[3, 1])

    graphItems = ["Classification diagnostic", "Green's-function map", "Ritz-mode map",
                  "Ritz spectrum map"]
    graphTypes = (:classification, :green, :ritzMode, :ritzSpectrum)
    critSymbols = visRitzCriterionNames()
    critItems = [String(c) for c in critSymbols]

    graphMenu, graphPrev, graphNext = guiDropdownStepperRow(topBarGL, 1, "Graph type:", graphItems)
    criterionMenu, critPrev, critNext = guiDropdownStepperRow(topBarGL, 2, "Criterion:", critItems)
    viewBtn = Button(topBarGL[3, 1], label = "View single angle ->", halign = :left, tellwidth = false)

    xLoDefault = string(minimum(sweepResult.angles))
    xHiDefault = string(maximum(sweepResult.angles))
    yLoDefault, yHiDefault = _aswDefaultYRange(sweepResult)

    state = AngSweepViewer(fig, nothing, sweepResult, contentGL, graphMenu, criterionMenu, nothing,
                            GuiRangeField(xLoDefault, xHiDefault), GuiRangeField(yLoDefault, yHiDefault),
                            GuiRangeField(), nothing, nothing, nothing, statusLabel, nothing, nothing,
                            viewBtn, nothing, Ref{Any}(nothing), Ref(1), Ref(false), Ref{Any}(nothing),
                            Ref{Union{Nothing,Function}}(nothing), () -> nothing, (x, y) -> nothing,
                            () -> nothing)

    function _buildMapContent!()
        ax, sideGL, _ = guiPlotWithSidePanel(state.contentGL[1, 1])
        state.ax = ax
        _, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "Angle [deg] min/max"; state = state.xRangeField)
        _, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "Velocity [km/s] min/max"; state = state.yRangeField)
        _, cLoBox, cHiBox = guiRangeFieldPair(sideGL, 5, "Colour min/max"; state = state.cRangeField)
        guiBindRangeToAxis!(ax, state.xRangeField, :x)
        guiBindRangeToAxis!(ax, state.yRangeField, :y)
        state.xRangeBoxes = (xLoBox, xHiBox)
        state.yRangeBoxes = (yLoBox, yHiBox)
        state.cRangeBoxes = (cLoBox, cHiBox)
        state.typeTogglePanel = guiGridLayout(sideGL[7, 1])
        state.infoLabel = Label(sideGL[8, 1], "Click the map to select a mode.";
                                 valign = :top, halign = :left)
        guiOnAxisClick(state.onClicked!, ax)
        return nothing
    end

    function _redraw!()
        ax = state.ax
        ax === nothing && return nothing
        guiEmptyAxis!(ax)
        graphType = graphTypes[state.graphTypeMenu.i_selected[]]
        criterion = critSymbols[state.modeMenu.i_selected[]]
        ok, errMsg = true, ""
        try
            if graphType === :ritzSpectrum
                plotRitzSpectrumMap(state.result; ax = ax)
            elseif graphType === :ritzMode
                plotRitzModeMap(state.result, criterion; ax = ax)
            elseif graphType === :green
                plotGreenMap(state.result; ax = ax)
            else   # :classification
                plotClassificationDiagnostic(state.result; ax = ax)
            end
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end

        if ok && graphType === :classification
            series = [(string(p.label[]), p) for p in ax.scene.plots
                      if p isa Makie.Scatter && haskey(p, :label)]
            guiRebuildTypeToggles!(state.typeTogglePanel, series)
        else
            guiRebuildTypeToggles!(state.typeTogglePanel, Tuple{String,Any}[])
        end

        guiRebindColorRange!(ax, state.cRangeField, state.cRangeBoxes, state.colorListener)
        state.crosshair = guiCrosshair!(ax)
        _aswReapplyRange!(ax, state.xRangeField, :x)
        _aswReapplyRange!(ax, state.yRangeField, :y)

        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _onClicked!(x::Float64, y::Float64)
        state.ax === nothing && return nothing
        angNo = _aswNearestIndex(state.result.angles, x)
        angNo === nothing && return nothing
        state.currentAngNo[] = angNo
        ritz = state.result.cs[angNo].ritzSolver
        vels = ritz.vels[1:ritz.nCompModes]
        modeNo = _aswNearestIndex(vels, y)
        modeNo === nothing && return nothing
        velSel = vels[modeNo]
        angleVal = state.result.angles[angNo]

        state.infoLabel.text[] = join(_aswModeInfoLines(state.result, angNo, modeNo, angleVal, velSel), "\n")
        state.crosshair === nothing || (state.crosshair[] = Point2f(angleVal, velSel))

        cb = state.onSelect[]
        cb === nothing || cb(angleVal, modeNo, velSel, y)
        return nothing
    end
    state.onClicked! = _onClicked!

    function _drillDownToAngle!(angNo::Integer)
        guiClearGridLayout!(state.contentGL)
        state.ax = nothing
        state.crosshair = nothing
        singleResult = state.result.cs[angNo]
        child = singleAngleViewer(singleResult; fig = state.contentGL[1, 1])
        child.onSelect[] = (angle, modeNo, velSnap, velRaw) -> begin
            cb = state.onSelect[]
            cb === nothing || cb(angle, modeNo, velSnap, velRaw)
            return nothing
        end
        state.childViewer[] = child
        state.currentAngNo[] = angNo
        state.viewSingleAngleBtn.label[] = "<< Back to angular sweep"
        state.isDrilledIn[] = true
        return nothing
    end

    function _returnToSweepView!()
        state.childViewer[] = nothing
        guiClearGridLayout!(state.contentGL)
        _buildMapContent!()
        state.redraw!()
        state.viewSingleAngleBtn.label[] = "View single angle ->"
        state.isDrilledIn[] = false
        return nothing
    end

    function _toggleDrillDown!()
        state.isDrilledIn[] ? _returnToSweepView!() : _drillDownToAngle!(state.currentAngNo[])
        return nothing
    end
    state.toggleDrillDown! = _toggleDrillDown!

    on(_ -> state.toggleDrillDown!(), viewBtn.clicks)
    on(_ -> state.redraw!(), graphMenu.i_selected)
    on(_ -> state.redraw!(), criterionMenu.i_selected)
    _aswWireStepper!(graphMenu, graphPrev, graphNext, length(graphItems))
    _aswWireStepper!(criterionMenu, critPrev, critNext, length(critItems))

    _buildMapContent!()
    state.redraw!()

    return state
end

# -----------------------------------------------------------------------------------------
# angSweepViewer -- Plate variant (CalcAngSweepPlate)
# -----------------------------------------------------------------------------------------

"""
    angSweepViewer(sweepResult::CalcAngSweepPlate; fig = nothing) -> AngSweepViewer

Angular-sweep map viewer for an already-`calc()`'d `CalcAngSweepPlate` result -- port of
`guliwerr/gui/AngSweepViewerGUIPlate.m`. No graph-type dropdown at all -- ABSENT, not disabled
(`graphTypeMenu` is `nothing`): there is exactly one map here
([`plotRitzModeMapPlate`](@ref)). `modeMenu` offers "Mode type" | "Confidence"
(`plotRitzModeMapPlate`'s own `:modeType`/`:confidence`); colour-range control is only
functionally meaningful for `:confidence` (a genuine continuous scale) -- see
[`guiRebindColorRange!`](@ref) for how `:modeType`'s fixed-per-type-RGB scatters make the
Colour min/max fields simply inert rather than visually disabled.

Throws `ArgumentError` if `sweepResult.velsMap` is empty (not yet `calc()`'d).
"""
function GuliwerrJ.angSweepViewer(sweepResult::CalcAngSweepPlate; fig = nothing)
    isempty(sweepResult.velsMap) && throw(ArgumentError(
        "angSweepViewer: sweepResult has not been calc()'d yet (velsMap is empty) -- nothing to view."))

    fig, mainLayout = guiHostLayout(fig; size = (1000, 700))
    topBarGL = guiGridLayout(mainLayout[1, 1])
    contentGL = guiGridLayout(mainLayout[2, 1])
    statusLabel = guiStatusLabel(mainLayout[3, 1])

    critItems = ["Mode type", "Confidence"]
    critSymbols = (:modeType, :confidence)
    criterionMenu, critPrev, critNext = guiDropdownStepperRow(topBarGL, 1, "Criterion:", critItems)
    viewBtn = Button(topBarGL[2, 1], label = "View single angle ->", halign = :left, tellwidth = false)

    xLoDefault = string(minimum(sweepResult.angles))
    xHiDefault = string(maximum(sweepResult.angles))
    yLoDefault, yHiDefault = _aswDefaultYRange(sweepResult)

    state = AngSweepViewer(fig, nothing, sweepResult, contentGL, nothing, criterionMenu, nothing,
                            GuiRangeField(xLoDefault, xHiDefault), GuiRangeField(yLoDefault, yHiDefault),
                            GuiRangeField(), nothing, nothing, nothing, statusLabel, nothing, nothing,
                            viewBtn, nothing, Ref{Any}(nothing), Ref(1), Ref(false), Ref{Any}(nothing),
                            Ref{Union{Nothing,Function}}(nothing), () -> nothing, (x, y) -> nothing,
                            () -> nothing)

    function _buildMapContent!()
        ax, sideGL, _ = guiPlotWithSidePanel(state.contentGL[1, 1])
        state.ax = ax
        _, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "Angle [deg] min/max"; state = state.xRangeField)
        _, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "Velocity [km/s] min/max"; state = state.yRangeField)
        _, cLoBox, cHiBox = guiRangeFieldPair(sideGL, 5, "Colour min/max"; state = state.cRangeField)
        guiBindRangeToAxis!(ax, state.xRangeField, :x)
        guiBindRangeToAxis!(ax, state.yRangeField, :y)
        state.xRangeBoxes = (xLoBox, xHiBox)
        state.yRangeBoxes = (yLoBox, yHiBox)
        state.cRangeBoxes = (cLoBox, cHiBox)
        state.infoLabel = Label(sideGL[7, 1], "Click the map to select a mode.";
                                 valign = :top, halign = :left)
        guiOnAxisClick(state.onClicked!, ax)
        return nothing
    end

    function _redraw!()
        ax = state.ax
        ax === nothing && return nothing
        guiEmptyAxis!(ax)
        criterion = critSymbols[state.modeMenu.i_selected[]]
        ok, errMsg = true, ""
        try
            plotRitzModeMapPlate(state.result, criterion; ax = ax)
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end
        if ok && state.yRangeField.range[] === nothing   # first draw: the map's own range
            reset_limits!(ax)
            yLo, yHi = ax.finallimits[].origin[2], ax.finallimits[].origin[2] + ax.finallimits[].widths[2]
            state.yRangeBoxes[1].stored_string[] = string(round(yLo; sigdigits = 4))
            state.yRangeBoxes[2].stored_string[] = string(round(yHi; sigdigits = 4))
        end

        guiRebindColorRange!(ax, state.cRangeField, state.cRangeBoxes, state.colorListener)
        state.crosshair = guiCrosshair!(ax)
        _aswReapplyRange!(ax, state.xRangeField, :x)
        _aswReapplyRange!(ax, state.yRangeField, :y)

        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _onClicked!(x::Float64, y::Float64)
        state.ax === nothing && return nothing
        angNo = _aswNearestIndex(state.result.angles, x)
        angNo === nothing && return nothing
        state.currentAngNo[] = angNo
        ritz = state.result.cs[angNo].ritzSolver
        vels = ritz.vels[1:ritz.nCompModes]
        modeNo = _aswNearestIndex(vels, y)
        modeNo === nothing && return nothing
        velSel = vels[modeNo]
        angleVal = state.result.angles[angNo]

        state.infoLabel.text[] = join(_aswModeInfoLines(state.result, angNo, modeNo, angleVal, velSel), "\n")
        state.crosshair === nothing || (state.crosshair[] = Point2f(angleVal, velSel))

        cb = state.onSelect[]
        cb === nothing || cb(angleVal, modeNo, velSel, y)
        return nothing
    end
    state.onClicked! = _onClicked!

    function _drillDownToAngle!(angNo::Integer)
        guiClearGridLayout!(state.contentGL)
        state.ax = nothing
        state.crosshair = nothing
        singleResult = state.result.cs[angNo]
        child = singleAngleViewer(singleResult; fig = state.contentGL[1, 1])
        child.onSelect[] = (angle, modeNo, velSnap, velRaw) -> begin
            cb = state.onSelect[]
            cb === nothing || cb(angle, modeNo, velSnap, velRaw)
            return nothing
        end
        state.childViewer[] = child
        state.currentAngNo[] = angNo
        state.viewSingleAngleBtn.label[] = "<< Back to angular sweep"
        state.isDrilledIn[] = true
        return nothing
    end

    function _returnToSweepView!()
        state.childViewer[] = nothing
        guiClearGridLayout!(state.contentGL)
        _buildMapContent!()
        state.redraw!()
        state.viewSingleAngleBtn.label[] = "View single angle ->"
        state.isDrilledIn[] = false
        return nothing
    end

    function _toggleDrillDown!()
        state.isDrilledIn[] ? _returnToSweepView!() : _drillDownToAngle!(state.currentAngNo[])
        return nothing
    end
    state.toggleDrillDown! = _toggleDrillDown!

    on(_ -> state.toggleDrillDown!(), viewBtn.clicks)
    on(_ -> state.redraw!(), criterionMenu.i_selected)
    _aswWireStepper!(criterionMenu, critPrev, critNext, length(critItems))

    _buildMapContent!()
    state.redraw!()

    return state
end

# -----------------------------------------------------------------------------------------
# angSweepViewer -- isoFreq variant (CalcIsofreqAngSweepLayered)
# -----------------------------------------------------------------------------------------

"""
    angSweepViewer(sweepResult::CalcIsofreqAngSweepLayered; fig = nothing) -> AngSweepViewer

isoFreq angular-sweep map viewer -- port of `guliwerr/gui/IsofreqAngSweepViewerGUI.m`. Exactly
one possible map -- the Green `|G_33|` field vs. (angle, `kParRange`)
([`plotDispersionMap`](@ref)'s `CalcIsofreqAngSweepLayered` method) -- since this class never
runs a Ritz solve at all. No `graphTypeMenu`/criterion dropdown (`nothing`); `modeMenu` is the
scale-mode dropdown (`log`/`linear`); `plotTypeMenu` is the `cartesian`/`polar` dropdown --
see this file's header, deviation 3, for clicks and ranges in the polar view.

`onSelect`, if set, is called `(angle, kParVal)` -- NOT the base/Plate variants'
`(angle, modeNo, velSnapped, velRaw)` -- matching `IsofreqAngSweepViewerGUI`'s own
`OnPointSelected` (vs. `OnModeSelected`) naming.

Throws `ArgumentError` if `sweepResult.absGMap` is empty (not yet `calc()`'d).
"""
function GuliwerrJ.angSweepViewer(sweepResult::CalcIsofreqAngSweepLayered; fig = nothing)
    isempty(sweepResult.absGMap) && throw(ArgumentError(
        "angSweepViewer: sweepResult has not been calc()'d yet (absGMap is empty) -- nothing to view."))

    fig, mainLayout = guiHostLayout(fig; size = (1000, 700))
    topBarGL = guiGridLayout(mainLayout[1, 1])
    contentGL = guiGridLayout(mainLayout[2, 1])
    statusLabel = guiStatusLabel(mainLayout[3, 1])

    scaleItems = ["log", "linear"]
    scaleSymbols = (:log, :linear)
    plotTypeItems = ["cartesian", "polar"]
    plotTypeSymbols = (:cartesian, :polar)
    _isPolar() = plotTypeSymbols[state.plotTypeMenu.i_selected[]] === :polar

    scaleMenu, scalePrev, scaleNext = guiDropdownStepperRow(topBarGL, 1, "Scale:", scaleItems)
    plotTypeMenu, ptPrev, ptNext = guiDropdownStepperRow(topBarGL, 2, "Plot type:", plotTypeItems)
    viewBtn = Button(topBarGL[3, 1], label = "View single angle ->", halign = :left, tellwidth = false)

    xLoDefault = string(minimum(sweepResult.angles))
    xHiDefault = string(maximum(sweepResult.angles))
    yLoDefault, yHiDefault = _aswDefaultYRange(sweepResult)

    state = AngSweepViewer(fig, nothing, sweepResult, contentGL, nothing, scaleMenu, plotTypeMenu,
                            GuiRangeField(xLoDefault, xHiDefault), GuiRangeField(yLoDefault, yHiDefault),
                            GuiRangeField(), nothing, nothing, nothing, statusLabel, nothing, nothing,
                            viewBtn, nothing, Ref{Any}(nothing), Ref(1), Ref(false), Ref{Any}(nothing),
                            Ref{Union{Nothing,Function}}(nothing), () -> nothing, (x, y) -> nothing,
                            () -> nothing)

    function _buildMapContent!()
        ax, sideGL, _ = guiPlotWithSidePanel(state.contentGL[1, 1])
        state.ax = ax
        _, xLoBox, xHiBox = guiRangeFieldPair(sideGL, 1, "Angle [deg] min/max"; state = state.xRangeField)
        _, yLoBox, yHiBox = guiRangeFieldPair(sideGL, 3, "kParRange min/max"; state = state.yRangeField)
        _, cLoBox, cHiBox = guiRangeFieldPair(sideGL, 5, "Colour min/max"; state = state.cRangeField)
        # the angle/kParRange fields set the limits in the Cartesian view only (deviation 3)
        on(r -> (r === nothing || _isPolar() || xlims!(ax, r); nothing), state.xRangeField.range)
        on(r -> (r === nothing || _isPolar() || ylims!(ax, r); nothing), state.yRangeField.range)
        state.xRangeBoxes = (xLoBox, xHiBox)
        state.yRangeBoxes = (yLoBox, yHiBox)
        state.cRangeBoxes = (cLoBox, cHiBox)
        state.infoLabel = Label(sideGL[7, 1], "Click the map to select a point.";
                                 valign = :top, halign = :left)
        guiOnAxisClick(state.onClicked!, ax)
        return nothing
    end

    function _redraw!()
        ax = state.ax
        ax === nothing && return nothing
        guiEmptyAxis!(ax)
        scaleMode = scaleSymbols[state.modeMenu.i_selected[]]
        plotType = plotTypeSymbols[state.plotTypeMenu.i_selected[]]
        # plotPolarHeatmap's aspect and hidden decorations outlive empty!(ax)
        ax.aspect = nothing
        for a in (:xlabelvisible, :ylabelvisible, :xticklabelsvisible, :yticklabelsvisible,
                  :xticksvisible, :yticksvisible, :xgridvisible, :ygridvisible,
                  :leftspinevisible, :rightspinevisible, :topspinevisible, :bottomspinevisible)
            getproperty(ax, a)[] = true
        end
        ok, errMsg = true, ""
        try
            plotDispersionMap(state.result, 1, scaleMode, (-5.0, -2.0), plotType; ax = ax)
        catch err
            ok, errMsg = false, sprint(showerror, err)
        end

        guiRebindColorRange!(ax, state.cRangeField, state.cRangeBoxes, state.colorListener)
        state.crosshair = guiCrosshair!(ax)
        if plotType === :polar
            autolimits!(ax)   # drop limits a Cartesian range left behind
        else
            _aswReapplyRange!(ax, state.xRangeField, :x)
            _aswReapplyRange!(ax, state.yRangeField, :y)
        end

        ok ? guiSetStatusOk!(state.statusLabel, "") : guiSetStatusError!(state.statusLabel, errMsg)
        return nothing
    end
    state.redraw! = _redraw!

    function _onClicked!(x::Float64, y::Float64)
        state.ax === nothing && return nothing
        if _isPolar()   # the axes hold pol2cart X/Y; nearest angle on the circle
            th, kPos = atand(y, x), hypot(x, y)
            angs = state.result.angles
            angNo = isempty(angs) ? nothing : argmin(abs.(mod.(angs .- th .+ 180, 360) .- 180))
        else
            kPos = y
            angNo = _aswNearestIndex(state.result.angles, x)
        end
        angNo === nothing && return nothing
        state.currentAngNo[] = angNo
        kIdx = _aswNearestIndex(state.result.kParRange, kPos)
        kIdx === nothing && return nothing
        kParVal = state.result.kParRange[kIdx]
        angleVal = state.result.angles[angNo]

        state.infoLabel.text[] = join(_aswPointInfoLines(state.result, angNo, kIdx, angleVal, kParVal), "\n")
        state.crosshair === nothing || (state.crosshair[] = _isPolar() ?
            Point2f(kParVal * cosd(angleVal), kParVal * sind(angleVal)) : Point2f(angleVal, kParVal))

        cb = state.onSelect[]
        cb === nothing || cb(angleVal, kParVal)
        return nothing
    end
    state.onClicked! = _onClicked!

    function _drillDownToAngle!(angNo::Integer)
        guiClearGridLayout!(state.contentGL)
        state.ax = nothing
        state.crosshair = nothing
        singleResult = state.result.cs[angNo]
        child = singleAngleViewer(singleResult; fig = state.contentGL[1, 1])
        child.onSelect[] = (angle, kParVal) -> begin
            cb = state.onSelect[]
            cb === nothing || cb(angle, kParVal)
            return nothing
        end
        state.childViewer[] = child
        state.currentAngNo[] = angNo
        state.viewSingleAngleBtn.label[] = "<< Back to angular sweep"
        state.isDrilledIn[] = true
        return nothing
    end

    function _returnToSweepView!()
        state.childViewer[] = nothing
        guiClearGridLayout!(state.contentGL)
        _buildMapContent!()
        state.redraw!()
        state.viewSingleAngleBtn.label[] = "View single angle ->"
        state.isDrilledIn[] = false
        return nothing
    end

    function _toggleDrillDown!()
        state.isDrilledIn[] ? _returnToSweepView!() : _drillDownToAngle!(state.currentAngNo[])
        return nothing
    end
    state.toggleDrillDown! = _toggleDrillDown!

    on(_ -> state.toggleDrillDown!(), viewBtn.clicks)
    on(_ -> state.redraw!(), scaleMenu.i_selected)
    on(_ -> state.redraw!(), plotTypeMenu.i_selected)
    _aswWireStepper!(scaleMenu, scalePrev, scaleNext, length(scaleItems))
    _aswWireStepper!(plotTypeMenu, ptPrev, ptNext, length(plotTypeItems))

    _buildMapContent!()
    state.redraw!()

    return state
end

# -----------------------------------------------------------------------------------------
# setAxisRange! -- MATLAB's AngSweepViewerGUI.setAxisRange (all three variants share this struct)
# -----------------------------------------------------------------------------------------

"""
    setAxisRange!(viewer::AngSweepViewer, axis::Symbol, lims) -> Nothing

Set one map range (`:x`, `:y` or `:c`) from code, as if typed into its min/max fields -- see
[`guiSetRangeFromCode!`](@ref). Throws while drilled into a single angle (the map's fields are
not mounted then) and for `:c` when the current plot has no data-driven colour range.
"""
function GuliwerrJ.setAxisRange!(viewer::AngSweepViewer, axis::Symbol, lims)
    owner = "angSweepViewer"
    # MATLAB's order: the arguments (lims = [] returns), then the drilled-in state, then the colour axis
    guiCheckRangeArgs(owner, axis, lims; hasColour = true) || return nothing
    viewer.isDrilledIn[] && throw(ArgumentError("$owner.setAxisRange!: the map's range fields " *
        "do not exist while drilled into a single angle -- return to the sweep view first"))
    p = axis === :c ? guiColorPlotObject(viewer.ax) : nothing
    guiCheckRangeArgs(owner, axis, lims; hasColour = p !== nothing)
    if axis === :c
        guiSetRangeFromCode!(owner, viewer.cRangeField, viewer.cRangeBoxes, lims, Float64.(p.colorrange[]))
    elseif axis === :x
        guiSetRangeFromCode!(owner, viewer.xRangeField, viewer.xRangeBoxes, lims, guiAxisLimits(viewer.ax, :x))
    else
        guiSetRangeFromCode!(owner, viewer.yRangeField, viewer.yRangeBoxes, lims, guiAxisLimits(viewer.ax, :y))
    end
    return nothing
end
