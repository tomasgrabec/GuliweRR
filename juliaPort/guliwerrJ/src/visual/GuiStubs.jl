# =========================================================================================
# GuiStubs.jl -- launcher stubs for the GLMakie/WGLMakie GUI extension
#
# Empty function stubs for every archetype's top-level launcher (7 total). Each stub's real
# implementation lives in `ext/GuliwerrJMakieExt.jl`'s own included files (a package
# extension, per `Project.toml`'s `[extensions]` table, loaded automatically whenever the end
# user does `using GuliwerrJ, GLMakie` or `using GuliwerrJ, WGLMakie` -- the extension
# triggers on `Makie` being loaded, not a specific backend).
#
# Same "stub now, extend later" pattern as `PlotStubs.jl`, so the GUI layer has the calling
# convention every `plotXxx` function has: `using GuliwerrJ, GLMakie; structureControl()`
# works directly, with no need to know about `Base.get_extension`.
#
# The 7 launchers are `structureControl`, `singleAngleViewer`, `angSweepViewer`,
# `modalShapeViewer`, `combinedViewer`, `doubleSweepViewer` and `dispersionViewer`. The
# extension must add methods to these stubs (`function GuliwerrJ.name(...)`) rather than
# define fresh bindings of the same name inside the extension module: a fresh binding is not
# exported from `GuliwerrJ`, so `structureControl()` would fail with `UndefVarError`.
# =========================================================================================

"""
    structureControl(; fig = nothing) -> NamedTuple

Archetype G: structure editor / job-launcher form. Port of `guliwerr/gui/StructureControlGUI.m`.
Implemented in the Makie extension; see `ext/GuiStructureControl.jl`.
"""
function structureControl end

"""
    singleAngleViewer(result; fig = nothing) -> SingleAngleViewer

Archetype B: 1D curve + click-to-select viewer, for an already-solved single-angle result
(`CalcSingleHS`/`CalcSingleLayered`/`CalcSinglePlate`/`CalcIsofreqSingleLayered`). Port of
`guliwerr/gui/SingleAngleViewerGUI.m` and its Plate/isoFreq siblings. Implemented in the Makie
extension; see `ext/GuiSingleAngleViewer.jl`.
"""
function singleAngleViewer end

"""
    angSweepViewer(sweepResult; fig = nothing) -> AngSweepViewer

Archetype A: 2D map + click-to-select + drill-down viewer, for an already-solved angular
sweep (`CalcAngSweepHS`/`CalcAngSweepLayered`/`CalcAngSweepPlate`/
`CalcIsofreqAngSweepLayered`). Port of `guliwerr/gui/AngSweepViewerGUI.m` and its Plate/isoFreq
siblings. Implemented in the Makie extension; see `ext/GuiAngSweepViewer.jl`.
"""
function angSweepViewer end

"""
    modalShapeViewer(result; fig = nothing) -> ModalShapeViewer

Archetype C: stacked modal-shape depth panels + steppers, for an already-solved single-angle
result or angular sweep across every family (base/Plate/isoFreq). Port of
`guliwerr/gui/ModalShapeViewerGUI.m` and its Plate/isoFreq siblings. Implemented in the Makie
extension; see `ext/GuiModalShapeViewer.jl`.
"""
function modalShapeViewer end

"""
    combinedViewer(result; windowMode = :combined, fig = nothing) -> CombinedViewer

Archetype D: two-pane composite shell wiring a [`modalShapeViewer`](@ref) (top) to a
[`singleAngleViewer`](@ref)/[`angSweepViewer`](@ref) (bottom, dispatched on `result`'s own
concrete type). Port of `guliwerr/gui/CombinedViewerGUI.m` and its Plate/isoFreq siblings.
Implemented in the Makie extension; see `ext/GuiCombinedViewer.jl`.
"""
function combinedViewer end

"""
    doubleSweepViewer(sweepResult; fig = nothing) -> DoubleSweepViewer

Archetype F: outer slice-stepper shell over one extra sweep axis, embedding a fresh
[`combinedViewer`](@ref) per slice (and, for the Layered/Plate methods, a real
[`dispersionViewer`](@ref) mode switch). Port of `guliwerr/gui/DoubleSweepViewerGUI.m` and its
Plate/isoFreq siblings. Implemented in the Makie extension; see `ext/GuiDoubleSweepViewer.jl`.
"""
function doubleSweepViewer end

"""
    dispersionViewer(sweepResult, angleIdx = 1; fig = nothing) -> DispersionViewer

Archetype E: map + shape grid in one window, with branch tracking, for an already-solved
2D wavevector/angle sweep (`CalcWavevecAngSweep2Layered`/`CalcWavevecAngSweep2Plate`). Port of
`guliwerr/gui/DispersionViewerGUI.m` and its Plate sibling. Implemented in the Makie extension;
see `ext/GuiDispersionViewer.jl`.
"""
function dispersionViewer end

"""
    setAxisRange!(viewer, axis::Symbol, lims) -> Nothing

Set one axis range of a viewer from code, as if typed into its min/max fields: `axis` is `:x`,
`:y` or `:c` (colour, map viewers only); `lims` is `(lo, hi)`, a `NaN` entry keeps that limit and
`[]` keeps both. Composite viewers forward to their map/curve viewer. Port of MATLAB's
`setAxisRange` (`GuiCommon.setAxisRange`). Implemented in the Makie extension.
"""
function setAxisRange! end

export structureControl, singleAngleViewer, angSweepViewer, modalShapeViewer, combinedViewer,
       doubleSweepViewer, dispersionViewer, setAxisRange!
