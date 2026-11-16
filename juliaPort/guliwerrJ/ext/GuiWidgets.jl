# =========================================================================================
# GuiWidgets.jl -- shared widget/layout substrate of the GLMakie viewers
#
# Julia/Makie port of `guliwerr/gui/GuiCommon.m` (25 static helpers) -- the reuse substrate
# every viewer bottoms out in. `include`d directly into `module GuliwerrJMakieExt` (see that
# file's own header), which already has `using GuliwerrJ`/`using Makie` in scope -- no
# repeated `using` here.
#
# This file adds NO methods to any `GuliwerrJ`-owned stub (unlike PlotsMaps.jl/
# PlotsSurfaceCurve.jl/PlotsModalShape.jl/PlotsDiagnostic.jl): every name below is a genuinely
# NEW function with no MATLAB-side "stub now, extend later" pattern to honour, since nothing
# outside this extension calls them. They are named with a `gui` prefix (not the `_group`
# underscore-prefix the four plotting-group files use for THEIR private helpers) precisely
# because they ARE meant to be found and reused across files: this is the foundation the
# viewer files (Gui*Viewer.jl, GuiStructureControl.jl) are built on, not a one-file
# private helper set.
#
# ---------------------------------------------------------------------------------------
# WHAT DISSOLVED, AND WHY (most of `GuiCommon.m` dissolves into Makie's GridLayout + Observable)
# ---------------------------------------------------------------------------------------
# - `hostLayout`'s `gridSize` argument: a Makie `GridLayout` grows to whatever index its
#   children are placed at; there is no upfront row/column count to declare, unlike
#   `uigridlayout(parent, gridSize)`. `Padding`: `GridLayout`'s own default gap is already
#   sane. What DOES survive as `guiHostLayout` below is the "own Figure vs. build into a
#   caller-supplied region, either way hand back the ancestor Figure" logic -- genuinely
#   non-trivial once you need `Makie.get_figure` (see `guiHostLayout`'s docstring) to
#   recover the Figure from a caller-supplied `GridPosition`.
# - `shapeAxesGrid` is DROPPED ENTIRELY, not ported even as a thin wrapper. It doesn't
#   dissolve into "call GridLayout directly" -- it dissolves because `ext/PlotsModalShape.jl`
#   makes a documented STRUCTURAL DEVIATION from the MATLAB source:
#   `plotModalShape`/`plotGreenModalShape`/`plotGreenModalShapeLayered` each own their
#   OWN Figure/Axis construction internally (`_modalShapeFigureAndAxes`, one `Axis` per
#   DIRECTION with multiple stacked `heatmap!` bands, not one `axes()` per (direction, layer)
#   pair) and take an optional `fig` to draw into, not a pre-built axes grid. There is no
#   `hAxCellIn`-shaped parameter left on the Julia side for a GUI to pre-build a grid for;
#   a viewer calls `GuliwerrJ.plotModalShape(...; fig = someFig)` directly.
# - `clickPoint` dissolves into `mouseposition(ax.scene)` inside `guiOnAxisClick` below.
#   MATLAB needed it because `event.IntersectionPoint` only exists on an OBJECT-level hit
#   (clicking a plotted line/scatter/pcolor), while an axes-level hit on empty background
#   carries no such field -- two cases, two code paths. Makie's `mouseposition` resolves
#   data-space coordinates the same way regardless of what, if anything, is drawn under the
#   cursor: one case, one code path.
# - `setSeriesVisible` dissolves into `plotobj.visible[] = v` directly (every Makie plot
#   object already exposes `.visible` as a settable `Observable`) -- inlined into
#   `guiRebuildTypeToggles!`'s own per-row listener below, nothing left to wrap.
# - `captureAxisRange`/`restoreAxisRange` dissolve by construction, not by a ported function:
#   MATLAB captured VALUES (never handles) because the field widgets themselves are destroyed
#   by the teardown/rebuild a drill-down (e.g. `AngSweepViewerGUI`'s single-angle view) causes.
#   `GuiRangeField` (below) puts the persistent state in `Observable`s owned by a struct that
#   outlives any number of widget teardown/rebuild cycles; passing that SAME struct back in as
#   `guiRangeFieldPair`'s `state` argument after a rebuild reproduces "restore" for free, and
#   there is nothing left to "capture" -- see `guiRangeFieldPair`'s own docstring.
#
# ---------------------------------------------------------------------------------------
# NOT IN THIS FILE
# ---------------------------------------------------------------------------------------
# - `significantModeType`/`plateModeTypeText` are NOT ported here. Both are physics-RESULT
#   lookups/formatters (reading `ModeClassification` entries / `RitzSolverLamb.modesType`
#   `.modeTypeConfidence` and calling `Vis.plateModeTypeStyle`), not generic widget/layout
#   mechanics -- they live next to the viewer files that build a mode-info panel
#   (GuiSingleAngleViewer.jl, as `guiSignificantModeType`/`guiPlateModeTypeText`).
#
# ---------------------------------------------------------------------------------------
# NAMING / STYLE
# ---------------------------------------------------------------------------------------
# Every genuinely GUI-generic name below is `gui`-prefixed camelCase (`guiHostLayout`,
# `guiRangeFieldPair`, ...) and documented with a `@ref`-style pointer back to the
# `GuiCommon.m` method it replaces, for exactly the reason PlotsMaps.jl's own header gives
# for ITS `_maps`-prefixed private helpers: several files share this one
# module namespace, and an unprefixed generic name (`newax`, `parseRange`, ...) would be a
# real collision risk. The difference is these names are meant to be PUBLIC within the
# extension (every viewer file reuses them), not private to one file, hence no
# leading underscore.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# Shared colours (GuiCommon.m's `properties (Constant)` block)
# -----------------------------------------------------------------------------------------

const GUI_STATUS_OK_COLOR    = RGBf(0.0, 0.45, 0.0)    # dark green -- confirmations/info
const GUI_STATUS_ERROR_COLOR = RGBf(0.6, 0.0, 0.0)     # dark red -- errors/blocked actions
const GUI_CROSSHAIR_COLOR    = RGBf(0.80, 0.40, 0.02)  # selection crosshair
const GUI_DIVIDER_COLOR      = RGBf(0.6, 0.6, 0.6)     # thin between-halves divider panel

# ===========================================================================================
# Status line
# ===========================================================================================

"""
    guiStatusLabel(parent) -> Label

Build the standard status label: empty text, [`GUI_STATUS_OK_COLOR`](@ref). Port of
`GuiCommon.statusLabel`. Place it with ordinary `GridLayout` indexing
(`parent[row, col] = guiStatusLabel(parent[row, col])`) -- there is no separate `rowNo`
argument to thread through the way MATLAB's optional second argument was, since a Makie
`GridPosition` already says where.
"""
guiStatusLabel(parent) = Label(parent, ""; color = GUI_STATUS_OK_COLOR, halign = :left, tellwidth = false)

"""
    guiSetStatusOk!(label::Label, msg) -> Nothing

Show a non-error status message (confirmation/info) in [`GUI_STATUS_OK_COLOR`](@ref). Port of
`GuiCommon.setStatusOk`.
"""
function guiSetStatusOk!(label::Label, msg::AbstractString)
    label.color[] = GUI_STATUS_OK_COLOR
    label.text[] = msg
    return nothing
end

"""
    guiSetStatusError!(label::Label, msg) -> Nothing

Show an error/blocked-action status message in [`GUI_STATUS_ERROR_COLOR`](@ref). Port of
`GuiCommon.setStatusError`.
"""
function guiSetStatusError!(label::Label, msg::AbstractString)
    label.color[] = GUI_STATUS_ERROR_COLOR
    label.text[] = msg
    return nothing
end

# Legends drawn on an axis by `extAxisLegend!`, so a redraw can delete them (`guiEmptyAxis!`).
const _AXIS_LEGENDS = WeakKeyDict{Axis,Vector{Legend}}()

"""
    extAxisLegend!(ax::Axis; kwargs...) -> Legend

`axislegend(ax; kwargs...)`, recording the legend against `ax`. Every plotting function of
this extension draws its legends through this, so [`guiEmptyAxis!`](@ref) can find them: a
`Legend` keeps no reference to its axis, and matching it by position proved unreliable (the
axis can move after a re-layout while the legend's stored box does not).
"""
function extAxisLegend!(ax::Axis; kwargs...)
    l = axislegend(ax; kwargs...)
    push!(get!(Vector{Legend}, _AXIS_LEGENDS, ax), l)
    return l
end

"""
    guiEmptyAxis!(ax::Axis) -> Nothing

`empty!(ax)` for a viewer redraw, plus deleting the legends earlier draws put on `ax`
([`extAxisLegend!`](@ref)). A legend is a separate `Block` of the figure, not a plot of `ax`,
so `empty!` leaves it behind and every redraw added one more on top (9 per
pass in `singleAngleViewer`).
"""
function guiEmptyAxis!(ax::Axis)
    for l in pop!(_AXIS_LEGENDS, ax, Legend[])
        l in ax.parent.content && delete!(l)
    end
    empty!(ax)
    return nothing
end

# ===========================================================================================
# Figure / layout scaffolding
# ===========================================================================================

"""
    guiGridLayout(args...; kwargs...) -> GridLayout

`GridLayout(args...; kwargs...)` that fills the space its parent cell gives it (`width =
nothing`), the way MATLAB's `uigridlayout` fills its container. Every GUI layout is built with
it. A plain `GridLayout` has `width = Auto()`, which shrinks it to the natural width of its
content whenever that width is determinable -- in a viewer made only of fixed-width widgets
(buttons, menus, fields; the labels do not report a width) that is a few hundred pixels, and
the whole window collapsed into a narrow centre column. Only a holder that must keep its
natural width (a colour bar beside its plot) uses plain `GridLayout`.
"""
guiGridLayout(args...; kwargs...) = GridLayout(args...; width = nothing, kwargs...)

"""
    guiHostLayout(parent = nothing; size = (1000, 700)) -> (fig::Figure, layout::GridLayout)

Host a viewer's content either in its own `Figure` or inside a caller-supplied region -- the
part of `GuiCommon.hostLayout` that survives the port (see this file's header for what does
not: `gridSize`, `Padding`).

`parent === nothing` (standalone): a fresh `Figure(size = size)` is created and `layout` is
its own root `GridLayout` (`fig.layout`).

`parent` given (embedded -- e.g. `someOtherFig[2, 1]`, a `GridPosition`): `layout` is a NEW
nested `GridLayout` built at that position, and `fig` is `parent`'s own ancestor `Figure`
(via `Makie.get_figure`, which is exactly the `ancestor(parentContainer, 'figure')` MATLAB
called for the same reason) -- so `display(fig)`/further layout at `fig[...]` works
identically regardless of which branch built it, the same duck-typed contract the MATLAB
`[fig, mainLayout]` return offered.

Not constructed `Visible = 'off'` the way `uifigure(...)` was for headless/test safety --
Makie has no such state to opt out of: a `Figure` is simply never shown until the caller
`display`s it (GLMakie) or renders it to a file (CairoMakie), so there is nothing here to
suppress.
"""
function guiHostLayout(parent = nothing; size::Tuple{<:Real,<:Real} = (1000, 700))
    if parent === nothing
        fig = Figure(size = size)
        fig.layout.width[] = nothing   # fill the window -- see guiGridLayout
        return fig, fig.layout
    else
        fig = Makie.get_figure(parent)
        layout = guiGridLayout(parent)
        return fig, layout
    end
end

"""
    guiContentControlsSplit(parent; ratio = (3, 1)) -> (contentGL, controlsGL, outer)

The shared "content left, controls right" two-column split (`GuiCommon.shapeControlsSplit`),
`ratio`-weighted column widths (MATLAB's `{'3x', '1x'}` default). `parent` is a `GridPosition`
(or anything `GridLayout(parent)` accepts) to build the split into. `outer` is the 1x2
`GridLayout` itself, returned so a caller inside a multi-row parent can place/resize it
further.
"""
function guiContentControlsSplit(parent; ratio::Tuple{<:Real,<:Real} = (3, 1))
    outer = guiGridLayout(parent)
    contentGL = guiGridLayout(outer[1, 1])
    controlsGL = guiGridLayout(outer[1, 2])
    colsize!(outer, 1, Relative(ratio[1] / sum(ratio)))
    colsize!(outer, 2, Relative(ratio[2] / sum(ratio)))
    return contentGL, controlsGL, outer
end

"""
    guiDividedTopBottom(parent, rowHeights = (Relative(0.44), Fixed(3), Relative(0.56));
                         ratio = (3, 1))
        -> (shapeGL, controlsGL, bottomGL, outer)

The shared combined-viewer body (`GuiCommon.dividedTopBottom`): a content/controls split (via
[`guiContentControlsSplit`](@ref)) on top, a thin divider `Box`, and a `bottomGL` region below
-- e.g. modal-shape panels + their controls above a map/curve viewer's content.

`rowHeights` is the 3-row sizing (top, divider, bottom) as `GridLayoutBase` row-size specs
(`Relative`/`Fixed`/`Auto`) -- e.g. the default reproduces MATLAB's `{'44x', 3, '56x'}` ratio
with a 3-pixel-thick divider row. The divider's own thickness is genuinely a per-family
choice, the same reasoning the MATLAB docstring gives for taking it as a parameter.
"""
function guiDividedTopBottom(parent,
                              rowHeights::Tuple{Any,Any,Any} = (Relative(0.44), Fixed(3), Relative(0.56));
                              ratio::Tuple{<:Real,<:Real} = (3, 1))
    outer = guiGridLayout(parent)
    shapeGL, controlsGL, _ = guiContentControlsSplit(outer[1, 1]; ratio = ratio)
    Box(outer[2, 1]; color = GUI_DIVIDER_COLOR, strokevisible = false)
    bottomGL = guiGridLayout(outer[3, 1])
    rowsize!(outer, 1, rowHeights[1])
    rowsize!(outer, 2, rowHeights[2])
    rowsize!(outer, 3, rowHeights[3])
    return shapeGL, controlsGL, bottomGL, outer
end

"""
    guiPlotWithSidePanel(parent; ratio = (3, 1)) -> (ax::Axis, sideGL::GridLayout, outer)

The shared map/curve viewer body: one plot `Axis` beside a control side panel, in a
`ratio`-weighted 2-column split (`GuiCommon.plotWithSidePanel`'s "3x/1x" default).

Unlike the MATLAB original, `sideGL` needs no pre-declared row count/heights: a Makie
`GridLayout` grows to whatever rows its caller places into it (`sideGL[k, 1] = ...`), so the
"row budget" bookkeeping `plotWithSidePanel`'s `sideRowHeight` argument existed for has
nothing left to do -- just place controls at increasing row indices.

Click handling is deliberately NOT wired here (unlike the MATLAB original's `onClick`
argument) -- see [`guiOnAxisClick`](@ref), which replaces both the axes-level `ButtonDownFcn`
and the per-plotted-object `ButtonDownFcn` re-attachment MATLAB needed after every redraw.
"""
function guiPlotWithSidePanel(parent; ratio::Tuple{<:Real,<:Real} = (3, 1))
    outer = guiGridLayout(parent)
    # the Axis gets its own sub-layout, so a colour bar drawn "beside the plot"
    # (PlotsMaps._maps_colorbar_beside: the next column of the Axis's own layout) lands next
    # to it instead of on top of the side panel
    ax = Axis(guiGridLayout(outer[1, 1])[1, 1])
    sideGL = guiGridLayout(outer[1, 2])
    colsize!(outer, 1, Relative(ratio[1] / sum(ratio)))
    colsize!(outer, 2, Relative(ratio[2] / sum(ratio)))
    return ax, sideGL, outer
end

# ===========================================================================================
# Stepper rows (label + selector + << >> buttons)
# ===========================================================================================

"""
    guiDropdownStepperRow(parent, row, labelText, itemStrs; onSelect = nothing)
        -> (menu::Menu, prevBtn::Button, nextBtn::Button)

Label + index-valued dropdown + `<<`/`>>` buttons, on one row (`GuiCommon.dropdownStepperRow`).
Built as its OWN nested `GridLayout` at `parent[row, 1]` so this row's 4-way column split does
not force every OTHER row placed into the same `parent` to share it -- MATLAB's identical
reason for building a fresh `uigridlayout(parent, [1 4])` per row.

No `ItemsData = 1:numel(itemStrs)` trick is needed: a Makie `Menu`'s `i_selected` `Observable`
already IS the 1-based index into `itemStrs` directly -- attach a callback to THAT (not
`.selection`, which carries the displayed string) to reproduce MATLAB's ".Value is the index,
never the displayed string" contract. `onSelect`, if given, is attached to `menu.i_selected`
(called with the new integer index).

Wiring `prevBtn`/`nextBtn` (what "previous"/"next" means -- clamped? wrapped? which axis?) is
left to the caller, exactly as in MATLAB.
"""
function guiDropdownStepperRow(parent, row::Integer, labelText::AbstractString,
                                itemStrs::AbstractVector{<:AbstractString};
                                onSelect::Union{Function,Nothing} = nothing)
    rowGL = guiGridLayout(parent[row, 1])
    Label(rowGL[1, 1], labelText; halign = :left, tellwidth = false)
    menu = Menu(rowGL[1, 2], options = itemStrs)
    prevBtn = Button(rowGL[1, 3], label = "<<")
    nextBtn = Button(rowGL[1, 4], label = ">>")
    onSelect === nothing || on(onSelect, menu.i_selected)
    return menu, prevBtn, nextBtn
end

"""
    guiIndexStepperRow(parent, row, labelText; init = 1, onChange = nothing)
        -> (box::Textbox, prevBtn::Button, nextBtn::Button, trailLabel::Label)

Label + positive-integer field + `<<`/`>>` buttons + a trailing free-text label, on one row
(`GuiCommon.indexStepperRow`). `validator = Int` is the Makie `Textbox` equivalent of
MATLAB's `Limits [1 Inf]` / `RoundFractionalValues 'on'` -- both express "this is a 1-based
index" at construction. `onChange`, if given, is called with the parsed `Int` (not the raw
string) whenever the box's `stored_string` changes, mirroring MATLAB's "ValueChangedFcn
always sees an already-validated numeric `.Value`" contract; a value that fails to parse
(should not happen given `validator = Int`, but defensively) is simply not forwarded.

`box.validator` is itself an `Observable`, so it can be relaxed at run time exactly the way
`ModalShapeViewerGUI.onSourceChanged` repurposes the same MATLAB field as a continuous
velocity entry (`box.validator[] = Float64`) -- a per-caller override, not a second shape
this helper builds, the same division of responsibility the MATLAB docstring itself draws.
"""
function guiIndexStepperRow(parent, row::Integer, labelText::AbstractString;
                             init::Integer = 1, onChange::Union{Function,Nothing} = nothing)
    rowGL = guiGridLayout(parent[row, 1])
    Label(rowGL[1, 1], labelText; halign = :left, tellwidth = false)
    box = Textbox(rowGL[1, 2], stored_string = string(init), validator = Int)
    prevBtn = Button(rowGL[1, 3], label = "<<")
    nextBtn = Button(rowGL[1, 4], label = ">>")
    trailLabel = Label(rowGL[1, 5], ""; halign = :left, tellwidth = false)
    if onChange !== nothing
        on(box.stored_string) do s
            v = tryparse(Int, s)
            v === nothing || onChange(v)
        end
    end
    return box, prevBtn, nextBtn, trailLabel
end

# ===========================================================================================
# Axis-range controls
# ===========================================================================================

"""
    guiParseRange(loStr, hiStr) -> Union{Nothing, Tuple{Float64,Float64}}

Parse a min/max text-field pair into `(lo, hi)`, or `nothing` if either side is blank,
unparseable, or `hi <= lo`. This is `GuiCommon.applyFieldLimit`'s validity check, lifted out
as a PURE function so it can be unit-tested without any GUI/Makie object at all. The "silently
do nothing on bad input" contract itself survives unchanged from the MATLAB original: a
half-typed value must not disturb the plot, so this never throws or clamps, only reports
"not a valid range right now" via `nothing`.
"""
function guiParseRange(loStr::AbstractString, hiStr::AbstractString)
    lo = tryparse(Float64, loStr)
    hi = tryparse(Float64, hiStr)
    (lo === nothing || hi === nothing) && return nothing
    hi <= lo && return nothing
    return (lo, hi)
end

"""
    GuiRangeField

Persistent state behind one "min/max" axis-range control: `loText`/`hiText`
(`Observable{String}`, the raw field contents) and `range` (`Observable{Union{Nothing,
Tuple{Float64,Float64}}}`, kept in sync via [`guiParseRange`](@ref) whenever either text
observable changes).

This is the Julia replacement for MATLAB's `savedAxisRange` struct / `captureAxisRange` /
`restoreAxisRange`: because the state lives here, in a struct that is never destroyed,
rather than being read off `Textbox` widgets that ARE destroyed on every teardown/rebuild,
nothing needs to "capture" or "restore" it -- see [`guiRangeFieldPair`](@ref)'s `state`
argument.

Construct via `GuiRangeField(loInit, hiInit)`, not the bare 3-argument struct constructor
directly (the outer constructor is what wires the `loText`/`hiText` -> `range` listeners).
"""
struct GuiRangeField
    loText::Observable{String}
    hiText::Observable{String}
    range::Observable{Union{Nothing,Tuple{Float64,Float64}}}
end

function GuiRangeField(loInit::AbstractString = "", hiInit::AbstractString = "")
    loText = Observable(String(loInit))
    hiText = Observable(String(hiInit))
    range = Observable{Union{Nothing,Tuple{Float64,Float64}}}(guiParseRange(loInit, hiInit))
    rf = GuiRangeField(loText, hiText, range)
    on(_ -> (rf.range[] = guiParseRange(rf.loText[], rf.hiText[])), loText)
    on(_ -> (rf.range[] = guiParseRange(rf.loText[], rf.hiText[])), hiText)
    return rf
end

"""
    guiRangeFieldPair(parent, row, labelText; loInit = "", hiInit = "", state = nothing)
        -> (rf::GuiRangeField, loBox::Textbox, hiBox::Textbox)

Mount a `labelText` label (row `row`) over a side-by-side min/max `Textbox` pair (row
`row + 1`) into `parent` -- `GuiCommon.rangeFieldPair`'s layout, consuming two rows exactly as
the MATLAB original did. Text, not numeric, fields, deliberately: a blank field has to mean
"leave this axis alone", which a numeric field cannot express -- see [`guiParseRange`](@ref).

Pass an EXISTING `state::GuiRangeField` (returned by a previous call) to re-mount fresh
widgets against state that already survived a teardown/rebuild cycle -- e.g. drilling down
into a single-angle view and back (`AngSweepViewerGUI`'s `savedAxisRange`/`captureAxisRange`/
`restoreAxisRange` in the MATLAB source). Because `rf.loText`/`rf.hiText` are `Observable`s
owned by the RETURNED struct, not by the (destroyable) `Textbox` widgets themselves, they
outlive any number of teardown/rebuild cycles on their own: keep the SAME `GuiRangeField`
across the rebuild and pass it back in as `state`, and the new `Textbox`es simply start from
wherever the old ones left off -- `rf.range` already reflects the persisted value the instant
the new widgets exist, no extra "restore" call needed.
"""
function guiRangeFieldPair(parent, row::Integer, labelText::AbstractString;
                            loInit::AbstractString = "", hiInit::AbstractString = "",
                            state::Union{GuiRangeField,Nothing} = nothing)
    rf = state === nothing ? GuiRangeField(loInit, hiInit) : state
    Label(parent[row, 1], labelText; halign = :left, tellwidth = false)
    fieldsGL = guiGridLayout(parent[row + 1, 1])
    loBox = Textbox(fieldsGL[1, 1], stored_string = rf.loText[], placeholder = "min")
    hiBox = Textbox(fieldsGL[1, 2], stored_string = rf.hiText[], placeholder = "max")
    on(s -> rf.loText[] = s, loBox.stored_string)
    on(s -> rf.hiText[] = s, hiBox.stored_string)
    return rf, loBox, hiBox
end

"""
    guiBindRangeToAxis!(ax::Axis, rf::GuiRangeField, kind::Symbol) -> Nothing

Apply `rf.range` to `ax`'s x- or y-limits (`kind === :x` / `:y`) whenever it changes, and once
immediately with the CURRENT value. This is `GuiCommon.applyFieldLimit`'s `case 'x'`/`case
'y'` branches, dissolved into a direct `xlims!`/`ylims!` call -- Makie axis limits are already
just an assignable property, so there is nothing left to wrap beyond "when does this run".

The colour-range case (`case 'c'`) is deliberately NOT covered here: which object owns
`colorrange` varies by plot type (a `Heatmap`, a `Scatter`, a `Colorbar`, ...), so wire that
one directly at the call site instead, e.g.
`on(rf.range) do r; r === nothing || (hm.colorrange[] = r); end`.
"""
function guiBindRangeToAxis!(ax::Axis, rf::GuiRangeField, kind::Symbol)
    kind in (:x, :y) ||
        throw(ArgumentError("guiBindRangeToAxis!: kind must be :x or :y (got $(repr(kind)))"))
    apply = r -> begin
        r === nothing && return
        kind === :x ? xlims!(ax, r) : ylims!(ax, r)
    end
    on(apply, rf.range)
    apply(rf.range[])
    return nothing
end

"""
    guiCheckRangeArgs(owner, axis::Symbol, lims; hasColour::Bool) -> Bool

Argument checks shared by every `setAxisRange!` method (`GuiCommon.setAxisRange`'s opening
block): `axis` is `:x`/`:y`/`:c`, `lims` is a 2-element real vector/tuple, and `:c` needs a
colour axis. Returns `false` for an empty `lims` (keep both limits -- the caller then does
nothing), `true` otherwise. Throws `ArgumentError("<owner>.setAxisRange!: ...")`.
"""
function guiCheckRangeArgs(owner::AbstractString, axis::Symbol, lims; hasColour::Bool)
    axis in (:x, :y, :c) ||
        throw(ArgumentError("$owner.setAxisRange!: axis must be :x, :y or :c (got $(repr(axis)))"))
    isempty(lims) && return false
    (length(lims) == 2 && all(l -> l isa Real, lims)) ||
        throw(ArgumentError("$owner.setAxisRange!: lims must be (lo, hi) (NaN keeps that limit) or []"))
    axis === :c && !hasColour &&
        throw(ArgumentError("$owner.setAxisRange!: this view has no colour range to set"))
    return true
end

"""
    guiAxisLimits(ax::Axis, axis::Symbol) -> Tuple{Float64,Float64}

`ax`'s current x or y limits as shown (`ax.finallimits`) -- MATLAB's `ax.XLim`/`ax.YLim`.
"""
function guiAxisLimits(ax::Axis, axis::Symbol)
    lim = ax.finallimits[]
    k = axis === :x ? 1 : 2
    return (Float64(lim.origin[k]), Float64(lim.origin[k] + lim.widths[k]))
end

"""
    guiSetRangeFromCode!(owner, rf::GuiRangeField, boxes::Tuple{Textbox,Textbox}, lims, curLim)
        -> Nothing

`GuiCommon.setAxisRange`'s field logic: merge `lims` with the current field text (a NaN entry
keeps that limit; a blank field keeps `curLim`'s, the axis' current limit) and write the
result into the mounted `boxes` via `Makie.set!` -- the same path as typing, so `rf` and every
axis/colour range bound to it follow through their existing listeners. Throws
`ArgumentError` for an empty/inverted result, where typing would silently leave the plot alone.
"""
function guiSetRangeFromCode!(owner::AbstractString, rf::GuiRangeField,
                              boxes::Tuple{Textbox,Textbox}, lims, curLim)
    cur = [something(tryparse(Float64, rf.loText[]), NaN), something(tryparse(Float64, rf.hiText[]), NaN)]
    cur = ifelse.(isnan.(cur), Float64.(collect(curLim)), cur)
    new = Float64.(collect(lims))
    new = ifelse.(isnan.(new), cur, new)
    new[2] > new[1] ||
        throw(ArgumentError("$owner.setAxisRange!: resulting range ($(new[1]), $(new[2])) is empty or inverted"))
    Makie.set!(boxes[1], string(new[1]))
    Makie.set!(boxes[2], string(new[2]))
    return nothing
end

# ===========================================================================================
# Click handling / crosshair
# ===========================================================================================

"""
    guiOnAxisClick(f, ax::Axis; button = Mouse.left) -> Observable

Call `f(x, y)` with the DATA-SPACE coordinates of a mouse click released while the pointer is
over `ax`.

Together, this replaces THREE things `GuiCommon` needed separately: `clickPoint`
(data-space conversion of a click event), the axes-level `ButtonDownFcn` (background clicks),
and the per-plotted-object `ButtonDownFcn` re-attachment every viewer's redraw needed (an
object-level click never reached a MATLAB axes' own callback, so every `set(hPlot,
'ButtonDownFcn', ...)` had to be repeated after every draw). Makie's `mouseposition(scene)`
already returns data coordinates the same way regardless of what, if anything, is drawn under
the cursor, so there is exactly one case to handle here, not three, and nothing to
re-attach after a redraw.

Returns the listener handle from `on(...)` (useful only if the caller wants to `off(...)` it
later, e.g. while a viewer is torn down).
"""
function guiOnAxisClick(f::Function, ax::Axis; button = Mouse.left)
    return on(events(ax).mousebutton) do event
        if event.button == button && event.action == Mouse.release &&
           Makie.is_mouseinside(ax.scene)
            pos = mouseposition(ax.scene)
            f(Float64(pos[1]), Float64(pos[2]))
        end
    end
end

"""
    guiCrosshair!(ax::Axis; color = GUI_CROSSHAIR_COLOR, markersize = 14) -> Observable{Point2f}

Create the selection crosshair marker on `ax` and return the `Observable{Point2f}` that
controls its position (`GuiCommon.markCrosshair`). Move (or "create", the first time) it with
`pos[] = Point2f(xVal, yVal)`.

Starts at `Point2f(NaN, NaN)`, which Makie simply does not draw -- there is no separate
create-vs-move branch to write the way MATLAB's `isempty(h) || ~isgraphics(h)` guard needed (a
stale handle from a destroyed axes was a real MATLAB concern there); the marker here is
created exactly once, alongside `ax`, and never needs re-creating.

No `HitTest = 'off'` equivalent is needed either: this crosshair is meant to be read alongside
[`guiOnAxisClick`](@ref)-style axis-level click handling (`mouseposition`, not per-object hit
testing), which never resolves a click against the crosshair marker itself in the first place.
"""
function guiCrosshair!(ax::Axis; color = GUI_CROSSHAIR_COLOR, markersize::Real = 14)
    pos = Observable(Point2f(NaN, NaN))
    scatter!(ax, pos; marker = :cross, markersize = markersize, color = color,
             strokewidth = 2, strokecolor = color)
    return pos
end

"""
    guiNearestPoint(seriesXs, seriesYs, query; xrange, yrange, visible = trues(length(seriesXs)))
        -> Union{Nothing, NamedTuple}

Nearest plotted point to `query = (qx, qy)` across several line/point series, in
axis-range-normalised distance -- the PURE-function core of `GuiCommon.nearestLinePoint`,
with the graphics-handle bookkeeping stripped out: a `visible` mask (one `Bool` per series)
takes the place of MATLAB's `isgraphics`/`Visible == 'off'` checks, and a miss returns
`nothing` instead of a caller-supplied default index (the caller decides what "keep the
previous selection" means at the call site, rather than this function needing to know it).

`seriesXs`/`seriesYs`: one coordinate-vector pair per series (e.g. one per mode branch/type
actually drawn). `xrange`/`yrange`: the caller's OWN currently-displayed axis spans (e.g. the
width/height of `ax.limits[]`) -- normalising by these, not raw data units, is why a click
matches sensibly even when x and y are on wildly different physical scales (wavelength in
microns vs. velocity in km/s), exactly as in MATLAB.

Returns `(series, index, dist2)` -- the winning series' position in `seriesXs`/`seriesYs`, the
index ALONG that series' own data, and the winning squared normalised distance -- or `nothing`
if every series is invisible or empty (nothing to hit).
"""
function guiNearestPoint(seriesXs::AbstractVector, seriesYs::AbstractVector,
                          query::Tuple{<:Real,<:Real};
                          xrange::Real, yrange::Real,
                          visible::AbstractVector{Bool} = trues(length(seriesXs)))
    (length(seriesXs) == length(seriesYs) == length(visible)) ||
        throw(ArgumentError("guiNearestPoint: seriesXs, seriesYs and visible must have equal length"))
    xr = max(Float64(xrange), eps())
    yr = max(Float64(yrange), eps())
    qx, qy = query

    bestDist = Inf
    bestSeries = 0
    bestIdx = 0
    for (im, (xs, ys, vis)) in enumerate(zip(seriesXs, seriesYs, visible))
        (vis && !isempty(xs)) || continue
        for (j, (x, y)) in enumerate(zip(xs, ys))
            dx = (Float64(x) - qx) / xr
            dy = (Float64(y) - qy) / yr
            d = dx^2 + dy^2
            if d < bestDist
                bestDist = d
                bestSeries = im
                bestIdx = j
            end
        end
    end
    bestSeries == 0 && return nothing
    return (series = bestSeries, index = bestIdx, dist2 = bestDist)
end

# ===========================================================================================
# Per-type visibility toggles
# ===========================================================================================

"""
    guiRebuildTypeToggles!(panel::GridLayout, series; header = "Show types:") -> Vector{Toggle}

(Re)build a "Show types:" toggle list, one row per `(label, plotobj)` pair in `series`
(`GuiCommon.rebuildTypeCheckboxes`) -- `plotobj` is any Makie plot object; each row's `Toggle`
directly drives that object's own `.visible` Observable.

Wholesale-replaces `panel`'s current content first (every `Toggle`/`Label` a PREVIOUS call
placed there), the same "delete everything, rebuild from the handles actually on the plot"
contract the MATLAB original used, so the toggle list never drifts from what is really drawn.
Pass an empty `series` to clear the panel (MATLAB's `gobjects(0)` case: any other graph type,
or after a failed draw).

`GuiCommon.setSeriesVisible` gets no separate Julia counterpart: toggling a series is just
`plotobj.visible[] = v`, which is exactly what each row's own listener does inline below.
"""
function guiRebuildTypeToggles!(panel::GridLayout, series; header::AbstractString = "Show types:")
    for gc in copy(panel.content)
        delete!(gc.content)
    end
    isempty(series) && return Toggle[]

    Label(panel[1, 1:2], header; fontsize = 13, halign = :left)
    toggles = Vector{Toggle}(undef, length(series))
    for (k, (label, plotobj)) in enumerate(series)
        t = Toggle(panel[k + 1, 1], active = true)
        Label(panel[k + 1, 2], label; halign = :left, tellwidth = false)
        on(v -> plotobj.visible[] = v, t.active)
        toggles[k] = t
    end
    return toggles
end

# ===========================================================================================
# Constructor-argument validation
# ===========================================================================================

"""
    guiValidateChoice(value, allowed; name = "value") -> value

Throw `ArgumentError` if `value ∉ allowed`, else return `value` unchanged
(`GuiCommon.validateWindowMode`).

Julia has no equivalent of MATLAB's `<className>:tag` error-identifier convention -- there is
no per-class identifier namespace to keep greppable, since Julia exceptions are matched by
TYPE and already carry the throwing call site in their stacktrace -- so no `className`
argument is threaded through here; pass a more specific `name` if the default "value" would
read too generically in context.
"""
function guiValidateChoice(value, allowed; name::AbstractString = "value")
    value in allowed ||
        throw(ArgumentError("$name must be one of $(allowed) (got $(repr(value)))."))
    return value
end

"""
    guiRequireStandaloneFor(splitValue, mode, hasParent::Bool; name = "windowMode") -> Nothing

Throw `ArgumentError` if `mode == splitValue` while `hasParent` is `true`
(`GuiCommon.requireStandaloneForSplit`) -- e.g. a `:split` window mode only makes sense when a
viewer owns its own window(s); once a caller hosts it, there is no second window to split
into.
"""
function guiRequireStandaloneFor(splitValue, mode, hasParent::Bool; name::AbstractString = "windowMode")
    hasParent && mode == splitValue && throw(ArgumentError(
        "$name = $(repr(splitValue)) only applies to standalone construction (omit the parent argument)."))
    return nothing
end

"""
    guiValidatePairedParents(a, b; names = ("controlsParent", "shapeParent")) -> Nothing

Throw `ArgumentError` if exactly one of `a`, `b` is `nothing`
(`GuiCommon.validateEmbedParents`'s both-or-neither contract: a control panel with nowhere to
draw shapes, or the reverse, is never meaningful).

Compose with [`guiRequireStandaloneFor`](@ref) to reproduce `validateEmbedParents`'s full
contract -- the two checks were only bundled together in MATLAB because both happen to apply
to the SAME two constructor arguments in the modal-shape viewer family; nothing requires them
to travel together in the Julia port, so they are kept as two independently-usable checks:

```julia
guiValidatePairedParents(controlsParent, shapeParent)
guiRequireStandaloneFor(:split, windowMode, controlsParent !== nothing)
```
"""
function guiValidatePairedParents(a, b; names::Tuple{AbstractString,AbstractString} = ("controlsParent", "shapeParent"))
    xor(a === nothing, b === nothing) && throw(ArgumentError(
        "$(names[1]) and $(names[2]) must be given together, or both omitted."))
    return nothing
end
