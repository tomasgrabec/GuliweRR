## Tutorial: layered full angular sweep -- Green's function on a film+substrate stack (Julia port)
#
# Same physics as examples/tutorials/tutorial_05_layered_green.m: the layered counterpart of
# tutorial 2's half-space sweep -- a film-on-substrate stack in :full calcMode (Ritz +
# Christoffel + Green). tutorial_06 (film on substrate) stays in :ritzOnly mode throughout (SAW
# velocity only). The layered launchers take one extra keyword, vMax; it is optional and only
# extends the Green velocity axis beyond its default end.
#
# Copy this script as a starting point for a new layered :full/Green calculation.
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - `AddLayer`, a plain function (not a method call) that returns a NEW `Structure` with one
#     more layer -- the free-function/immutability pattern from tutorial 1, now used to build a
#     multi-layer stack instead of a single half-space

using GuliwerrJ
using GLMakie   # opens an interactive window; `using CairoMakie` instead renders the same
                 # Figure/Axis objects headlessly (e.g. for save("fig.png", fig)) -- see
                 # tutorial_02_angular_map.jl for the full explanation

## Substrate -- Cu(001), cubic, from the materials library (`ListMaterialLibrary(; print = true)`
# shows values and sources)
dirMatSub = [1 0 0;    # surface normal
             0 0 1]     # propagation reference (Material's default)

## Film -- Ni(001), cubic (constants and density from the materials library)
dirMatFilm = [1 0 0;    # surface normal
              0 0 1]     # propagation reference (Material's default)

## Discretisation -- "easy" tier (docs/pitfalls.md MTH1): both layers are (001) cuts, a
# principal plane, with Zener ratios 2*C44/(C11-C12) of 3.2 (Cu) and 2.7 (Ni), both in the
# "easy" tier as a starting point. The film relDepth is set directly here;
# tutorial_06/07 derive it from a real thickness with ApplyWavelength. An entry's vel is its
# Green peak velocity. A leaky entry (above the substrate's shear velocity) is bound to the
# Ritz mode nearest that peak, but the slab has no mode exactly at a leaky resonance: the
# mode's own velocity (modeVel) can sit up to half a mode spacing (here a few m/s) off the
# peak, at any tier.
relDepthSub  = 20.0   # substrate depth, in wavelengths
relDepthFilm = 0.2    # film thickness, in wavelengths
legOrder     = 50     # Legendre basis order per layer; keep >= 2*relDepth
                      # (docs/pitfalls.md MTH1)

matSub, _ = MaterialFromLibrary("Cu"; dirMat = dirMatSub)
cpSub  = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)

matFilm, _ = MaterialFromLibrary("Ni"; dirMat = dirMatFilm)
cpFilm  = CalcParams(; relDepth = relDepthFilm, legOrder = legOrder)

strName   = "cu_sub_ni_film"
nMaxModes = 500
str = Structure(; materIn = matSub, calcParamsIn = cpSub, name = strName, nMaxModes = nMaxModes)
# layer 1 = substrate
str = AddLayer(str, matFilm, cpFilm)   # layer 2 = film, on top
# `AddLayer(str, matFilm, cpFilm)` -- a plain function, POSITIONAL (material, then calcParams),
# taking `str` as its first argument and returning a NEW `Structure` with the extra layer
# appended. This replaces MATLAB's `str = str.AddLayer(matFilm, cpFilm)` method-call syntax --
# same "reassign to the same name" idiom, since `str` is immutable and cannot be extended in
# place.

## Green velocity-axis upper end -- optional. nothing ends the axis at 1.05 x the fastest
# longitudinal limiting velocity of any layer (the end of the band the classification
# searches); a larger value only extends the displayed Green map.
vMax = nothing

## Angular sweep, 0:2:90 deg ([100] -> [010]), full mode (Ritz+Christoffel+Green)
angles         = 0.0:2.0:90.0
calcMode       = :full
retentionLevel = :full   # :full | :saw_only | :none -- :full keeps modal shapes for a viewer
sweep = angSweepLayered(; structure = str, angles = angles, vMax = vMax, calcMode = calcMode,
                           retentionLevel = retentionLevel)
sweep = calc(sweep)
# There is no `nVelGreen` keyword for layered stacks (MATLAB's slot for it is an unused empty
# argument): the Green axis has a fixed 4 m/s step.

## |G13_eff| map -- the spectrum the classification read; log scale, colour range
# (-3, -0.5) instead of the (-4, -1) default (a narrower, higher window that keeps the strong
# branches distinct; lower climRange[1] to see weaker ones). Pass e.g.
# `component = "G33"` to recompute |G33| instead.
scaleMode = :log          # :log (default) | :linear
climRange = (-3.0, -0.5)  # log10(|G|/max|G|) colour limits
axG13, cbG13 = plotGreenMap(sweep, scaleMode, climRange)
axG13.title = "Cu substrate / Ni film, [100]→[010]: |G_13^eff|"
display(axG13.parent)
# The colour range is passed straight to plotGreenMap as its own `climRange` argument, rather
# than MATLAB's post-hoc `clim([-3 -0.5])` on an already-drawn axis.

## Classification diagnostic: classified points and the substrate shear reference line over
# the Green map. On a stack every accepted wave has the single type :SAW -- no Rayleigh/Sezawa
# split (docs/pitfalls.md PHY2); marker size scales with the score.
#
# plotClassificationDiagnostic's OWN Green background always uses the default colour range, so
# to get the same (-3, -0.5) range the background is drawn first by plotGreenMap into an Axis we
# create, and the classified points are then drawn on top of it with showGreenBackground =
# false (`ax = ...` hands both functions the SAME axis).
colorBy             = :type    # :type | :score | :confidence
showGreenBackground = false    # background already drawn by plotGreenMap below
refLevel            = :full    # :sawOnly | :sawLimit (+ substrate shear line, layered) |
                               # :full (+ limiting-type points)
figDiag = Figure()
axDiag  = Axis(figDiag[1, 1])
plotGreenMap(sweep, scaleMode, climRange; ax = axDiag)
plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel; ax = axDiag)
axDiag.title = "Cu substrate / Ni film, [100]→[010]: mode classification vs. |G_13^eff|"
display(figDiag)
