# =========================================================================================
# GuiCombinedViewer.jl -- composite viewer ("two-pane composite shell: a modal-shape view
# above, an angular-sweep-or-single-angle view below")
#
# Julia/Makie port of `guliwerr/gui/CombinedViewerGUI.m`, `CombinedViewerGUIPlate.m` and
# `IsofreqCombinedViewerGUI.m`. In MATLAB, the base class wires
# `ModalShapeViewerGUI` (top) to `AngSweepViewerGUI`-or-`SingleAngleViewerGUI` (bottom,
# dispatched on the result's own class in `buildBottomViewer` -- `CalcSingleHS`/
# `CalcSingleLayered` get the 1D curve, `CalcAngSweepHS`/`CalcAngSweepLayered` get the 2D map);
# Plate wires the Plate pair identically; isoFreq wires the isoFreq pair with the renamed
# `OnPointSelected`/`goToPoint`/`showCrosshairAt` interface (`kParRange` is continuous, not a
# discrete Ritz mode index) in place of `OnModeSelected`/`goToSelection`. `include`d directly
# into the EXISTING `module GuliwerrJMakieExt` (see that file's own header), after
# `GuiModalShapeViewer.jl` -- this file needs `modalShapeViewer`/`ModalShapeViewer`
# (`GuiModalShapeViewer.jl`), `angSweepViewer`/`AngSweepViewer` (`GuiAngSweepViewer.jl`),
# `singleAngleViewer`/`SingleAngleViewer` (`GuiSingleAngleViewer.jl`) and
# `guiHostLayout`/`guiDividedTopBottom`/`guiValidateChoice`/`guiRequireStandaloneFor`
# (`GuiWidgets.jl`).
#
# THIS FILE BUILDS NOTHING NEW GRAPHICALLY -- it only instantiates, lays out, and wires
# together three child viewers, as a "no-plotting composite host owning layout and child
# lifecycle."
#
# ---------------------------------------------------------------------------------------
# THE TWO-FIGURES OBSTACLE -- THE ALWAYS-SEPARATE SHAPE WINDOW IS ACCEPTED, NOT RE-IMPLEMENTED
# ---------------------------------------------------------------------------------------
# `GuiModalShapeViewer.jl`'s own header ("TWO FIGURES, ALWAYS") describes the problem:
# `plotModalShape`/`plotGreenModalShape`/
# `plotGreenModalShapeLayered` (`ext/PlotsModalShape.jl`) each call `empty!(fig)` on whatever
# `Figure` they are handed, via `_modalShapeFigureAndAxes`, before rebuilding their OWN
# `GridLayout` at `fig[1,1]` -- they require ownership of an ENTIRE `Figure`, never a
# caller-supplied `GridPosition` sub-region, so a `ModalShapeViewer`'s `shapeFig` can never
# safely share a window with a controls panel or another viewer's content.
#
# Two options exist: (a) accept the always-separate shape window, or (b) rebuild the drawing
# logic against a caller-owned `GridPosition`.
# Option (b) would mean re-implementing, in THIS file, everything currently living in
# `PlotsModalShape.jl`'s private
# `_modalShape`-prefixed helpers: the data-driven, symmetric-about-zero colour limit
# (`_modalShapeColLimit`), the substrate-cutoff depth calculation feeding three different
# `substrateCutoff` modes (`_modalShapeLayeredPlot`'s `:decay`/`:fixedMultiple`/`:boosted`
# branches, calling into `calcSubstrateCutoffDepth`/`calcModalShapeRowHeights`), the per-layer
# row-height-fraction bookkeeping (`_modalShapeRowHeightFracs`/`_modalShapeLayerBands`), and the
# actual per-layer-banded heatmap+frame drawing (`_modalShapeDrawColumn!`) -- roughly 150 lines
# of GENUINELY PHYSICS/DISPLAY-BEARING logic, not boilerplate. Every one of those helpers is
# `_modalShape`-prefixed, i.e. deliberately private to that file (see the headers of
# `GuiSingleAngleViewer.jl`/`GuiAngSweepViewer.jl`: the viewer files never depend on each
# other's internals) -- reusing them directly from here would violate that convention, and
# duplicating them here instead would create a SECOND, independently-maintained copy of
# colour-range/cutoff/row-height physics that WILL silently drift from the canonical version
# the next time either file is touched.
#
# Weighed against that, option (a)'s cost is purely a UI-convenience loss, the SAME
# category `GuiModalShapeViewer.jl`'s own header accepts ("one visual
# window instead of two arranged side by side ... every angle/mode/velocity/source control and
# every value it can show is identical either way"). THIS FILE THEREFORE CHOOSES (a):
# reuse `modalShapeViewer` wholesale, unmodified, as the top pane. A `CombinedViewer` built in
# `:combined` windowMode (see below) therefore owns TWO real Makie `Figure`s -- the
# controls+bottom-viewer window this file itself lays out, and `ModalShapeViewer`'s own
# perpetually-separate `shapeFig`, popped up alongside it (never embedded) -- and THREE Figures
# total under `:split` mode (see below). The `:combined` window's own top strip carries a
# one-line hint label pointing at the separate shape window, precisely because that region can
# never hold real shape content (see `_cvBuild` below).
#
# ---------------------------------------------------------------------------------------
# `windowMode` COLLAPSES FROM THREE CHOICES TO TWO -- 'fullSplit' IS INDISTINGUISHABLE FROM
# 'split' ONCE `ModalShapeViewer` ITSELF ALREADY ALWAYS SPLITS
# ---------------------------------------------------------------------------------------
# MATLAB's three-way `windowMode` ('combined'/'split'/'fullSplit') exists because
# `ModalShapeViewerGUI` ITSELF has a two-way `windowMode` (shapes+controls together, or
# separate) that `CombinedViewerGUI` forwards a choice into: 'combined' = one window holding
# everything; 'split' = `ModalShapeViewerGUI` combined (its own shapes+controls together, one
# window) plus the bottom viewer standalone (one window) = 2 windows; 'fullSplit' =
# `ModalShapeViewerGUI` SPLIT (shapes/controls in two separate windows of their own) plus the
# bottom viewer standalone = 3 windows.
#
# `GuiModalShapeViewer.jl` collapses `ModalShapeViewerGUI`'s OWN two-way choice to ONE
# case (`modalShapeViewer`'s `fig` keyword controls only whether its CONTROLS are embedded;
# `shapeFig` is unconditionally, always separate -- see that file's header). That means the one
# axis MATLAB's 'split' vs 'fullSplit' distinguished (is the modal-shape viewer's OWN
# shapes/controls split further, or not) no longer exists as a choice on the Julia side AT ALL
# -- it is always "yes, split." So `combinedViewer`'s own `windowMode` keyword here has only TWO
# meaningful values:
#   - `:combined` (default) -- one shared `Figure` holds this file's OWN controls-strip-plus-
#     bottom-viewer layout (`ModalShapeViewer`'s controls embedded at the top via its own `fig`
#     keyword, the bottom viewer's content embedded below via ITS OWN `fig` keyword);
#     `shapeFig` pops up separately regardless (see above). Two windows total.
#   - `:split` -- covers what were MATLAB's DISTINCT 'split' and 'fullSplit' (now identical in
#     this port): every child is built fully standalone (`fig = nothing` throughout). Three
#     windows total (`ModalShapeViewer`'s controls window + its `shapeFig` + the bottom
#     viewer's own window) -- matching MATLAB's 'fullSplit' window COUNT, under the
#     'split' NAME, since that is now the only non-`:combined` case left.
# Validated via [`guiValidateChoice`](@ref); embedding (`fig` keyword) is rejected under
# `:split` via [`guiRequireStandaloneFor`](@ref) -- the helper for "an
# already-fully-split viewer has no single window to embed into".
#
# ---------------------------------------------------------------------------------------
# BIDIRECTIONAL WIRING -- OUTGOING HALF IS A DIRECT, EXACT CALL; INCOMING HALF
# (`showCrosshairAt`) IS A DOCUMENTED, NARROW APPROXIMATION FOR CURVE CHILDREN
# ---------------------------------------------------------------------------------------
# OUTGOING (bottom-viewer click -> shape viewer follows): exact, no approximation --
# `bottomViewer.onSelect[]` is wired straight to `modalViewer.goToSelection!`/`goToPoint!`,
# precisely MATLAB's own `OnModeSelected`/`OnPointSelected` wiring.
#
# INCOMING (shape-viewer state change -> mark the bottom viewer's own crosshair, WITHOUT
# triggering its click/mode-lookup logic -- MATLAB's `showCrosshairAt`, see
# `AngSweepViewerGUI.showCrosshairAt` and `SingleAngleViewerGUI.showCrosshairAt`).
# Neither `AngSweepViewer` nor `SingleAngleViewer` has such a method; this file implements
# the incoming half from outside, through their public surface:
#   - `AngSweepViewer` (2D map), not drilled in: EXACT -- `crosshair[] = Point2f(angle, val)`.
#     The map's own two axes ARE (angle, velocity-or-kParRange), so no y-value lookup is needed
#     at all, unlike MATLAB's curve case below.
#   - `AngSweepViewer`, drilled in: reproduces MATLAB's OWN re-drill-on-mismatch behaviour
#     (`showCrosshairAt`'s `if angNo ~= obj.currentAngNo; obj.drillDownToAngle(angNo); end`)
#     using ONLY `AngSweepViewer`'s PUBLIC surface -- `currentAngNo[]` and
#     `toggleDrillDown!()` called twice (once to return to the sweep view, once to re-drill at
#     the new angle) -- then forwards to the resulting child.
#   - `SingleAngleViewer` (1D curve, standalone OR as an `AngSweepViewer`'s drilled-down child):
#     APPROXIMATED. MATLAB's own `showCrosshairAt` snaps the marker's y-coordinate onto
#     whichever curve is CURRENTLY displayed (`valueAtVelocity`, ported as
#     `GuiSingleAngleViewer.jl`'s private, four-result-type x three-graph-type x
#     three-criterion `_savValueAtVelocity`/`_savValueAtKPar`). Reusing that private, multi-
#     dispatch helper from here would break the "no cross-file private dependency"
#     convention the TWO-FIGURES decision above also follows; duplicating its full dispatch
#     table here would be a second copy of real display logic with the same drift risk. Given
#     the OUTGOING half already carries the real physical link (a click always drives the
#     shape viewer correctly), and this incoming half is a purely cosmetic convenience in the
#     OTHER direction, this file sets the marker's x (velocity/kParRange) exactly but fixes its
#     y at `0.0` -- the same "no precise lookup available" fallback value
#     `_savValueAtVelocity` itself uses for its own Green-data-unavailable case. The
#     crosshair still lands at the correct velocity/kParRange; only its vertical snap onto the
#     exact curve height is approximate, for this one, purely-cosmetic direction of the link.
#
# NO REENTRANCY GUARD IS NEEDED: the incoming half above NEVER calls
# `onSelect`/`onClicked!` on the bottom viewer -- it only ever assigns the leaf `crosshair`
# Observable (or re-drills via the already-guarded `toggleDrillDown!`), so the outgoing ->
# incoming -> outgoing cycle terminates after exactly one hop in each direction (see
# `_cvWireModeFamily!` and `_cvWirePointFamily!` below).
# =========================================================================================

_gcvNearestAngle(angles::AbstractVector{<:Real}, angleVal::Real) = argmin(abs.(angles .- angleVal))

# -----------------------------------------------------------------------------------------
# CombinedViewer -- the one handle type returned by every variant
# -----------------------------------------------------------------------------------------

"""
    CombinedViewer

Handles for one "two-pane composite shell" viewer, returned by every
[`combinedViewer`](@ref) method. Mirrors `CombinedViewerGUI`'s own property list, including
keeping the field name `angSweepViewer` even when it holds a `SingleAngleViewer` -- MATLAB's own
class docblock gives the reason directly: "The property keeps one name across both cases,
external code needing only the shared duck-typed interface... never the concrete class."

Fields:
- `fig` -- the shared controls+bottom-viewer `Figure` under `:combined` windowMode; `nothing`
  under `:split` (no single shared window -- see this file's header).
- `modalShapeViewer::ModalShapeViewer` -- the top-pane child (`ext/GuiModalShapeViewer.jl`).
- `angSweepViewer` -- the bottom-pane child: an `AngSweepViewer` (angular-sweep `result`) or a
  `SingleAngleViewer` (single-angle `result`), untyped so this one field serves both -- matching
  MATLAB's own "one name, two possible concrete classes" precedent (see above).
- `windowMode::Symbol` -- `:combined` (default) or `:split` -- see this file's header for why
  MATLAB's third choice ('fullSplit') collapsed into `:split` in this port.
"""
mutable struct CombinedViewer
    fig::Union{Figure,Nothing}
    modalShapeViewer::ModalShapeViewer
    angSweepViewer::Any
    windowMode::Symbol
end

# -----------------------------------------------------------------------------------------
# Incoming crosshair sync ("showCrosshairAt") -- private to this file (`_gcv`-prefixed),
# dispatched on the bottom viewer's own concrete type. See this file's header for the exact
# (map) vs. approximated (curve) split.
# -----------------------------------------------------------------------------------------

"""
    _gcvShowCrosshairAt!(viewer::AngSweepViewer, angleVal::Real, valAxisVal::Real) -> Nothing
    _gcvShowCrosshairAt!(viewer::SingleAngleViewer, angleVal::Real, valAxisVal::Real) -> Nothing

Move (or create) `viewer`'s own crosshair marker to reflect a selection made in the LINKED
modal-shape viewer, WITHOUT triggering any click/mode-lookup logic on `viewer` itself -- see
this file's header ("BIDIRECTIONAL WIRING") for the exact/approximated split and for why no
reentrancy guard is needed here.

Both methods take the SAME `(angleVal, valAxisVal)` pair (the modal-shape viewer's own current
angle and velocity/kParVal) -- mirroring MATLAB's own dual signature exactly:
`AngSweepViewerGUI.showCrosshairAt(obj, angleVal, velVal)` uses `angleVal` as the map's x-axis;
`SingleAngleViewerGUI.showCrosshairAt(obj, ~, velVal)` IGNORES the angle argument (a curve
viewer has only one, fixed angle) and uses `velVal` as the curve's own x-axis instead -- see
that MATLAB method's own docblock: "First argument (angle) accepted-but-unused -- kept only so
this method's call signature matches AngSweepViewerGUI's own... letting CombinedViewerGUI wire
either class identically."
"""
function _gcvShowCrosshairAt!(viewer::AngSweepViewer, angleVal::Real, valAxisVal::Real)
    if viewer.isDrilledIn[]
        angNo = _gcvNearestAngle(viewer.result.angles, angleVal)
        if angNo != viewer.currentAngNo[]
            viewer.currentAngNo[] = angNo
            viewer.toggleDrillDown!()   # back to the sweep view
            viewer.toggleDrillDown!()   # re-drill at the new angle
        end
        child = viewer.childViewer[]
        child === nothing || _gcvShowCrosshairAt!(child, angleVal, valAxisVal)
    else
        viewer.crosshair === nothing || (viewer.crosshair[] = Point2f(angleVal, valAxisVal))
    end
    return nothing
end

function _gcvShowCrosshairAt!(viewer::SingleAngleViewer, ::Real, valAxisVal::Real)
    # y fixed at 0.0 -- see this file's header for why an exact per-curve lookup is not
    # reproduced here.
    viewer.crosshair[] = Point2f(valAxisVal, 0.0)
    return nothing
end

# -----------------------------------------------------------------------------------------
# _cvBuild -- shared shell for every variant; differs only in `wireFn!` (goToSelection! vs.
# goToPoint! call signature) and `rowHeights` (44/56 base+Plate vs. 34/66 isoFreq, matching
# each MATLAB variant's own `dividedTopBottom` ratio).
# -----------------------------------------------------------------------------------------

function _cvBuild(result, modalBuilderFn::Function, bottomBuilderFn::Function, wireFn!::Function;
                   windowMode::Symbol = :combined, fig = nothing,
                   rowHeights::Tuple{Any,Any,Any} = (Relative(0.44), Fixed(3), Relative(0.56)))
    windowMode = guiValidateChoice(windowMode, (:combined, :split); name = "windowMode")
    guiRequireStandaloneFor(:split, windowMode, fig !== nothing)

    local figOut, modalViewer, bottomViewer
    if windowMode === :combined
        figOut, layout = guiHostLayout(fig; size = (1100, 900))
        # ratio flipped to (1, 3) relative to guiDividedTopBottom's own (3, 1) default: the
        # "content" half of the top strip holds only a hint label (see this file's header, the
        # two-Figures decision), while the "controls" half holds ModalShapeViewer's real,
        # functional widgets -- the reverse of MATLAB's own weighting, where the content half
        # held the actual shape heatmaps and so genuinely needed the extra width.
        shapeGL, controlsGL, bottomGL, _ = guiDividedTopBottom(layout[1, 1], rowHeights;
                                                                 ratio = (1, 3))
        Label(shapeGL[1, 1], "Modal shape shown\nin a separate window ->";
              halign = :center, valign = :center, fontsize = 13,
              tellwidth = false, tellheight = false)
        modalViewer = modalBuilderFn(result; fig = controlsGL[1, 1])
        bottomViewer = bottomBuilderFn(result; fig = bottomGL[1, 1])
    else
        modalViewer = modalBuilderFn(result; fig = nothing)
        bottomViewer = bottomBuilderFn(result; fig = nothing)
        figOut = nothing
    end

    wireFn!(modalViewer, bottomViewer)

    return CombinedViewer(figOut, modalViewer, bottomViewer, windowMode)
end

"""
    _cvWireModeFamily!(modalViewer::ModalShapeViewer, bottomViewer) -> Nothing

Bidirectional wiring for the base/Plate families (`goToSelection!`/`(angle, modeNo, velSnap,
velRaw)`) -- MATLAB's `OnModeSelected`/`goToSelection` pair.
"""
function _cvWireModeFamily!(modalViewer::ModalShapeViewer, bottomViewer)
    bottomViewer.onSelect[] = (angle, modeNo, velSnap, velRaw) -> begin
        modalViewer.goToSelection!(angle, modeNo, velRaw)
        return nothing
    end
    modalViewer.onSelectionChanged[] = (angle, vel) -> begin
        _gcvShowCrosshairAt!(bottomViewer, angle, vel)
        return nothing
    end
    return nothing
end

"""
    _cvWirePointFamily!(modalViewer::ModalShapeViewer, bottomViewer) -> Nothing

Bidirectional wiring for the isoFreq family (`goToPoint!`/`(angle, kParVal)`) -- MATLAB's
`OnPointSelected`/`goToPoint` pair.
"""
function _cvWirePointFamily!(modalViewer::ModalShapeViewer, bottomViewer)
    bottomViewer.onSelect[] = (angle, kParVal) -> begin
        modalViewer.goToPoint!(angle, kParVal)
        return nothing
    end
    modalViewer.onSelectionChanged[] = (angle, kParVal) -> begin
        _gcvShowCrosshairAt!(bottomViewer, angle, kParVal)
        return nothing
    end
    return nothing
end

# -----------------------------------------------------------------------------------------
# combinedViewer -- one method per family
# -----------------------------------------------------------------------------------------

"""
    combinedViewer(result::Union{CalcSingleHS,CalcSingleLayered,CalcAngSweepHS,
                                  CalcAngSweepLayered};
                   windowMode::Symbol = :combined, fig = nothing) -> CombinedViewer

Two-pane composite viewer: a modal-shape depth-panel view (top, [`modalShapeViewer`](@ref))
wired to an angular-sweep map or single-angle curve view (bottom, [`angSweepViewer`](@ref) /
[`singleAngleViewer`](@ref), dispatched on `result`'s own concrete type -- a genuine sweep gets
the map, a single-point result gets the curve, exactly `CombinedViewerGUI.buildBottomViewer`'s
own dispatch). Port of `guliwerr/gui/CombinedViewerGUI.m`.

`windowMode` -- `:combined` (default) or `:split`; see this file's header for why MATLAB's
third choice ('fullSplit') collapsed into `:split` in this port, and for why `:combined` still
spawns TWO real `Figure`s (a `ModalShapeViewer`'s `shapeFig` is always separate). `fig` --
`nothing` (default), or a `GridPosition` to embed the SHARED `:combined`-mode window into
(rejected under `:split`, which has no single window to embed into).

Throws whatever [`modalShapeViewer`](@ref)/[`angSweepViewer`](@ref)/[`singleAngleViewer`](@ref)
themselves throw for a not-yet-solved `result` -- this function performs no solver call and no
extra validation of its own.
"""
function GuliwerrJ.combinedViewer(result::Union{CalcSingleHS,CalcSingleLayered,CalcAngSweepHS,CalcAngSweepLayered};
                         windowMode::Symbol = :combined, fig = nothing)
    bottomBuilder = result isa Union{CalcSingleHS,CalcSingleLayered} ? singleAngleViewer : angSweepViewer
    return _cvBuild(result, modalShapeViewer, bottomBuilder, _cvWireModeFamily!;
                     windowMode = windowMode, fig = fig)
end

"""
    combinedViewer(result::Union{CalcSinglePlate,CalcAngSweepPlate};
                   windowMode::Symbol = :combined, fig = nothing) -> CombinedViewer

Plate counterpart -- port of `guliwerr/gui/CombinedViewerGUIPlate.m`. Same dispatch/wiring as
the base method, against the Plate sibling viewers.
"""
function GuliwerrJ.combinedViewer(result::Union{CalcSinglePlate,CalcAngSweepPlate};
                         windowMode::Symbol = :combined, fig = nothing)
    bottomBuilder = result isa CalcSinglePlate ? singleAngleViewer : angSweepViewer
    return _cvBuild(result, modalShapeViewer, bottomBuilder, _cvWireModeFamily!;
                     windowMode = windowMode, fig = fig)
end

"""
    combinedViewer(result::Union{CalcIsofreqSingleLayered,CalcIsofreqAngSweepLayered};
                   windowMode::Symbol = :combined, fig = nothing) -> CombinedViewer

isoFreq counterpart -- port of `guliwerr/gui/IsofreqCombinedViewerGUI.m`. Wires the renamed
`goToPoint!`/`(angle, kParVal)` interface ([`_cvWirePointFamily!`](@ref)) rather than
`goToSelection!`/`(angle, modeNo, velRaw)` -- `kParRange` is a continuous axis here, not a
discrete Ritz mode. Row-height split is 34/66 (top/bottom), matching
`IsofreqCombinedViewerGUI.m`'s own comment: "the isoFreq shape half is a single Green's-response
row, not a {3 x nLay} stack, so it needs less of the window."
"""
function GuliwerrJ.combinedViewer(result::Union{CalcIsofreqSingleLayered,CalcIsofreqAngSweepLayered};
                         windowMode::Symbol = :combined, fig = nothing)
    bottomBuilder = result isa CalcIsofreqSingleLayered ? singleAngleViewer : angSweepViewer
    return _cvBuild(result, modalShapeViewer, bottomBuilder, _cvWirePointFamily!;
                     windowMode = windowMode, fig = fig,
                     rowHeights = (Relative(0.34), Fixed(3), Relative(0.66)))
end

"""
    setAxisRange!(viewer::CombinedViewer, axis::Symbol, lims) -> Nothing

Set one range of the bottom map/curve viewer from code -- forwarded to its own
`setAxisRange!`.
"""
GuliwerrJ.setAxisRange!(viewer::CombinedViewer, axis::Symbol, lims) =
    GuliwerrJ.setAxisRange!(viewer.angSweepViewer, axis, lims)
