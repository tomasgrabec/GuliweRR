"""
    GuliwerrJMakieExt

Package extension providing Makie-based plotting for GuliwerrJ. Loads automatically when
the end user does `using GuliwerrJ, GLMakie` (or `WGLMakie`) — see `Project.toml`'s
`[extensions]` table, which triggers on `Makie` alone (backend-agnostic; either concrete
backend loads `Makie` as a real dependency, which is enough to activate this extension).

Every plotting function here ADDS A METHOD to a generic function `GuliwerrJ.Visual`
already declared as an empty stub in `src/visual/PlotStubs.jl` (the "stub now, extend
later" pattern already used internally by `CalcCommon.jl`'s `calcGreen`/`identifyModes`/
`applyRetention`) — methods are added via the fully-qualified form
`function GuliwerrJ.plotXxx(...) ... end`, never a fresh `function plotXxx(...) ... end`,
which would create an unrelated second binding instead of extending the real one.

One module per top-level folder applies here too, in spirit: this file only `include`s
grouped sub-files as plain top-level code sharing this module's single `using Makie`/
`using GuliwerrJ` block — no per-file submodule.

Makie's rendering has no numeric MATLAB counterpart, so the plotting functions are
checked structurally rather than against reference values: they run on real solver output,
return the documented Makie object type, and the number of plotted series/points, axis
labels, and data range match the input arrays.
"""
module GuliwerrJMakieExt

using GuliwerrJ
using Makie

include("GuiWidgets.jl")
include("GuiStructureControl.jl")
include("PlotsSurfaceCurve.jl")
include("PlotsMaps.jl")
include("PlotsModalShape.jl")
include("PlotsDiagnostic.jl")
include("GuiSingleAngleViewer.jl")
include("GuiAngSweepViewer.jl")
include("GuiModalShapeViewer.jl")
include("GuiDispersionViewer.jl")
include("GuiCombinedViewer.jl")
include("GuiDoubleSweepViewer.jl")

end # module GuliwerrJMakieExt
