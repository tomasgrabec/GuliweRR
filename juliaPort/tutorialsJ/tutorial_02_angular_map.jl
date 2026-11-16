## Tutorial: anisotropic half-space -- full angular sweep + visualisation (Julia port)
#
# Same physics as examples/tutorials/tutorial_02_angular_map.m: Cu(001) surface, propagation
# swept [100] -> [010]. A :full sweep (Ritz + Christoffel + Green together), then the standard
# plots in the order you would reach for them: SAW velocity/slowness vs. angle, the |G_13,eff|
# Green's-function map (cartesian and polar), the classification diagnostic (classified
# SAW/limiting points over that map), every Ritz mode (exploratory), and a recomputed
# one-sided |G_13+| map for comparison.
#
# Copy this script as a starting point for a new angular-sweep calculation. See
# juliaPort/docsJ/quickStartJ.md for the call sequence, juliaPort/docsJ/pitfallsJ.md for the
# gotchas specific to the Julia port, and docs/pitfalls.md for the legO/relDepth rule (all
# still apply -- the physics and its numerical requirements are unchanged from the MATLAB
# original; only the calling language differs).
#
# How to run: see tutorial_01_isotropic_halfspace.jl's header for the two ways (terminal vs.
# a running REPL session via `include(...)`) -- not repeated in full for every tutorial after
# this point. With a plotting backend the REPL route is the one to use: a window opened by a
# script run from the terminal closes again as soon as the script ends.
#
# Julia concepts introduced in this tutorial (tutorial 1's concepts -- using, structs, keyword
# constructors, immutability, multiple dispatch, Symbols, 1-based indexing -- are assumed
# familiar from here on and not re-explained):
#   - matrix literals (`[1 0 0; 0 0 1]`) -- IDENTICAL syntax to MATLAB, nothing to relearn
#   - loading a plotting backend with its own `using` statement, and the GLMakie/CairoMakie
#     choice; showing a figure explicitly with `display`
#   - functions returning MULTIPLE values as a tuple, destructured at the call site
#     (`ax, cb = plotGreenMap(sweep)`)
#   - Ranges (`0:2:90`) as a compact, lazy stand-in for a MATLAB colon-expression vector

using GuliwerrJ

# Plotting needs a SECOND `using` statement, on top of `using GuliwerrJ` -- GuliwerrJ's plotting
# functions (plotAngularMap, plotGreenMap, ...) exist only as empty placeholders until a Makie
# backend is loaded (see src/visual/PlotStubs.jl's own docstring: "a package extension"). This
# keeps GuliwerrJ itself lightweight for someone who only wants the numerics -- the whole
# plotting machinery is opt-in.
#
#   using GLMakie     -- opens a real, interactive window; what this tutorial series uses,
#                        since you are meant to explore these plots by eye (zoom, rotate, ...).
#   using CairoMakie  -- produces the EXACT SAME Figure/Axis objects, but renders to a static
#                        image (PNG/PDF/SVG) instead of a window -- no display needed. Useful
#                        on a remote server, in a batch job, or simply to `save("fig.png", fig)`
#                        for a report. Swap this line for the one below to use it; every plotting
#                        call in this script works unchanged either way, since both backends
#                        implement the same Makie interface -- Figure/Axis/heatmap!/scatter!/...
#                        mean the same thing regardless of which one is loaded.
using GLMakie

## Material -- Cu(001), cubic. dirMat is Material's default here, set explicitly to show
# where an orientation goes (tutorial_S1 builds others).
dirMat   = [1 0 0;    # surface normal
            0 0 1]     # propagation reference (Material's default)
# This 2-row-by-3-column matrix literal is EXACTLY MATLAB syntax: spaces separate columns,
# `;` (here, a line break inside the brackets) separates rows. Nothing to relearn here.

# The Cu constants and density come from the materials library (`ListMaterialLibrary(; print = true)`
# shows values and sources). `MaterialFromLibrary` returns the Material and the library row; the
# row is not needed here, so it is discarded with `_`.
mat, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)

## Discretisation -- "easy" tier (docs/pitfalls.md MTH1): Cu's Zener ratio
# 2*C44/(C11-C12) = 3.2 < 5 on a principal plane. Use the "complex" tier
# CalcParams(; relDepth = 30, legOrder = 70) for a non-principal cut or a Zener ratio of 5-10.
relDepth = 20   # truncation depth of the computational domain, in wavelengths
legOrder = 50   # number of Legendre basis functions spanning that depth;
                # must be >= 2*relDepth (docs/pitfalls.md MTH1)
cp = CalcParams(; relDepth = relDepth, legOrder = legOrder)

strName   = "cu_001_100to010"
nMaxModes = 500
str = Structure(; materIn = mat, calcParamsIn = cp, name = strName, nMaxModes = nMaxModes)

## Angular sweep, 0:2:90 deg ([100] -> [010]), full mode (Ritz+Christoffel+Green)
angles         = 0.0:2.0:90.0
# `0.0:2.0:90.0` is a Julia RANGE, not an Array -- a compact (start, step, stop) triple that
# behaves like a vector wherever one is needed (it supports indexing, length, iteration, ...)
# without ever allocating all 46 numbers until something genuinely needs them materialised.
# GuliwerrJ's sweep constructors accept it directly (they call `collect` on it internally),
# exactly mirroring MATLAB's own `(0:2:90)'` colon-expression column vector.
calcMode       = :full       # :full (Ritz+Christoffel+Green) | :ritzOnly (SAW velocity only,
                             # no Green map -- much faster)
retentionLevel = :saw_only   # :full | :saw_only | :none -- :full only needed for modal-shape
                             # plots (tutorial 3; docs/pitfalls.md NUM4)
nVelGreen      = nothing     # number of Green-axis points (nothing = greenVelStep steps, 6 m/s)
useParallel    = :auto       # :auto (threads for >= 20 angles) | :always | :never (debugging) --
                             # threads only exist if Julia was started with `-t auto`
sweep = angSweepHS(; structure = str, angles = angles, calcMode = calcMode,
                      retentionLevel = retentionLevel, nVelGreen = nVelGreen,
                      useParallel = useParallel)
sweep = calc(sweep)
# `nothing` is Julia's "no value given" -- the counterpart of MATLAB's `[]` for an optional
# argument: the keyword falls back to its default.

## First look -- SAW velocity vs. angle (also available from a :ritzOnly sweep).
# Cartesian (velocity) and polar (slowness).
axCart = plotAngularMap(sweep)
axCart.title = "Cu (001), [100]→[010]: SAW velocity vs. angle"
display(axCart.parent)
# `plotAngularMap` returns a plain `Axis` (a Makie plotting pane). Julia has no `title(...)`
# free function for this the way MATLAB does -- axis properties are set by assigning directly
# to a field, `ax.title = "..."`, after the plotting call returns the axis.
#
# `display(axCart.parent)` is what puts the figure on screen: a plot function only BUILDS the
# figure (an Axis's `.parent` is its owning Figure). In a REPL a figure that is the value of a
# line is shown automatically; in a script nothing is, so each figure is `display`ed
# explicitly. With CairoMakie loaded `display` does not open a window (use `save("name.png",
# fig)` there).

plotType = :polar   # :cartesian (default, velocity) | :polar (slowness = 1/velocity)
axPolar = plotAngularMap(sweep, plotType)
axPolar.thetalimits = (0, π/2)   # cubic (001): 0-90 deg covers every direction (angles in rad)
axPolar.title = "Cu (001), [100]→[010]: SAW slowness (polar)"
display(axPolar.parent)

## |G13_eff| map (log scale by default) -- the ±k-symmetrised TGS response, the spectrum
# the classification reads and the only Green component stored
scaleMode = :log           # :log (default, log10(|G|/max|G|)) | :linear (raw |G|)
climRange = (-3.0, -1.0)   # default for :log is (-4, -1) (1e-4..1e-1 of the peak); automatic for :linear
axG13, cbG13 = plotGreenMap(sweep, scaleMode, climRange)
ylims!(axG13, 1, 5.5)   # the Green axis starts at 0.3 km/s, far below any surface wave
axG13.title = "Cu (001), [100]→[010]: |G_13^eff|"
display(axG13.parent)
# `plotGreenMap` returns a TUPLE of two objects -- the Axis AND its Colorbar -- unlike
# `plotAngularMap`'s single Axis above. Julia functions often return several values packed into
# one tuple, and the call site "destructures" them into separate names in one line, exactly as
# written here (`axG13, cbG13 = ...`) -- no cell array `{}` unpacking needed, and no fixed
# `nargout` count to worry about: `ax, cb = plotGreenMap(sweep)` always gives you both.
# `ylims!(ax, lo, hi)` is Makie's counterpart of MATLAB's `ylim([lo hi])`.

## Same map, polar plotType -- angle x velocity heatmaps (plotGreenMap, plotDispersionMap)
# all share this cartesian/polar toggle as their positional plotType argument.
plotType  = :polar         # :cartesian (default) | :polar
axG13Polar, cbG13Polar = plotGreenMap(sweep, scaleMode, climRange, plotType)
limits!(axG13Polar, 0, 5.5, 0, 5.5)   # plain Axis, velocity as radius: crop to 5.5 km/s
axG13Polar.title = "Cu (001), [100]→[010]: |G_13^eff| (polar)"
display(axG13Polar.parent)

## Classification diagnostic: classified SAW/limiting points and the longitudinal Christoffel
# reference line over the Green map, to check the classified points against the Green peaks.
colorBy              = :type    # :type (categorical) | :score (classification score, dB) |
                                # :confidence (|v_Ritz - v_GreenPeak|)
showGreenBackground  = false    # true draws |G13_eff| underneath at the default (-4, -1) range
                                # (a :ritzOnly result: its Ritz spectrum); here plotGreenMap
                                # draws it first with climRange, as in tutorial 05
refLevel             = :full    # :sawOnly (SAW points only) | :sawLimit (+ longitudinal reference
                                # line, HS) | :full (+ limiting-type points)
figDiag = Figure()
axDiag  = Axis(figDiag[1, 1])
plotGreenMap(sweep, scaleMode, climRange; ax = axDiag)
plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel; ax = axDiag)
ylims!(axDiag, 1, 5.5)
axDiag.title = "Cu (001), [100]→[010]: mode classification vs. |G_13^eff|"
display(figDiag)
# These are the default :green classification at its :high standard. Switch without
# recomputing: `reclassify(sweep, :low)` (another standard), `classify(sweep, :ritz)` or
# `classify(sweep, :refined)` (another method; :refined solves its spot directions once) --
# docs/quickStart.md, "Choosing the classification"; quickStartJ.md's "SAW classification:
# method and standard" has the Julia spelling. Both return the sweep (mutated in place).

## EXPLORATORY: every computed Ritz mode (not just the classified ones), shaded by a
# criterion -- shows what the classification left out.
criterion       = :depthInsensitivity   # :ritzScore (default) | :freqDisp | :depthInsensitivity
velCutoffFactor = 1.10   # hide modes above 1.1x the max longitudinal limiting velocity
                         # (nothing = default, Inf = none)
axRitz = plotRitzModeMap(sweep, criterion, velCutoffFactor)
ritzClim = (1.0, 8.0)              # colour limits of the criterion
for p in axRitz.scene.plots
    p.color[] isa AbstractVector{<:Real} && (p.colorrange[] = ritzClim)
end
for c in axRitz.parent.content
    c isa Colorbar && (c.limits[] = ritzClim)
end
axRitz.title = "Cu (001), [100]→[010]: all Ritz modes (depthInsensitivity)"
# Makie has no `clim`: the colour range sits on each colour-mapped plot (`colorrange`) and on
# the Colorbar (`limits`), so both are set.
display(axRitz.parent)

## Optional: any other Green component, e.g. the out-of-plane |G_33|.
# Not stored: recomputed here, one solve per angle.
component = "G33"   # "G13eff" (default, stored) | any "GIJ" / "GIJeff", e.g. "G13", "G33"
axG33, cbG33 = plotGreenMap(sweep, scaleMode, climRange, :cartesian; component = component)
ylims!(axG33, 1, 5.5)
axG33.title = "Cu (001), [100]→[010]: |G_33| (recomputed)"
display(axG33.parent)
# `component` is a keyword argument (after the `;` -- keywords are always named in Julia), where
# MATLAB takes it as its 6th positional argument.
