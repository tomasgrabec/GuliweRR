## Tutorial 02 (code) -- anisotropic half-space: Cu(001) angular sweep + standard plots
# Narrated version: juliaPort/tutorialsJ/tutorial_02_angular_map.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_02_angular_map.jl")

using GuliwerrJ
using GLMakie                      # or CairoMakie for headless/static output

## Material -- Cu, cubic, (001) surface
dirMat      = [1 0 0; 0 0 1]       # row 1 = surface normal, row 2 = propagation at angle 0
material, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)   # constants + density from the materials library

## Discretisation -- "easy" tier (anisotropy 3.2, principal plane)
relDepth   = 20.0
legOrder   = 50
calcParams = CalcParams(; relDepth = relDepth, legOrder = legOrder)

## Structure
structureName = "cu_001_100to010"
nMaxModes     = 500
structure     = Structure(; materIn = material, calcParamsIn = calcParams,
                            name = structureName, nMaxModes = nMaxModes)

## Angular sweep -- Ritz + Christoffel + Green
anglesDeg      = 0.0:2.0:90.0      # [100] -> [010]
calcMode       = :full             # :full | :ritzOnly
retentionLevel = :saw_only         # :full | :saw_only | :none
nVelGreen      = nothing           # Green-axis points (nothing = greenVelStep steps, 6 m/s)
useParallel    = :auto             # :auto | :always | :never (threads need `julia -t auto`)
sweep = angSweepHS(; structure = structure, angles = anglesDeg, calcMode = calcMode,
                     retentionLevel = retentionLevel, nVelGreen = nVelGreen,
                     useParallel = useParallel)
sweep = calc(sweep)

## SAW velocity vs. angle -- cartesian (velocity) and polar (slowness)
# (a figure is only built by the plot call; `display(ax.parent)` shows it -- GLMakie)
plotTypeCartesian = :cartesian     # :cartesian | :polar
axVel = plotAngularMap(sweep, plotTypeCartesian)
axVel.title = "Cu (001): SAW velocity"
display(axVel.parent)
plotTypePolar = :polar
axSlow = plotAngularMap(sweep, plotTypePolar)
axSlow.thetalimits = (0, π/2)
axSlow.title = "Cu (001): SAW slowness"
display(axSlow.parent)

## |G_13,eff| Green map (the spectrum the classification read) -- cartesian and polar
scaleMode = :log                   # :log | :linear
climRange = (-3.0, -1.0)           # log10(|G|/max|G|) colour limits (default (-4, -1))
axGreen, _ = plotGreenMap(sweep, scaleMode, climRange, plotTypeCartesian)
ylims!(axGreen, 1, 5.5)
axGreen.title = "Cu (001): |G_13^eff|"
display(axGreen.parent)
axGreenPolar, _ = plotGreenMap(sweep, scaleMode, climRange, plotTypePolar)
limits!(axGreenPolar, 0, 5.5, 0, 5.5)
axGreenPolar.title = "Cu (001): |G_13^eff| (polar)"
display(axGreenPolar.parent)

## Classification diagnostic -- classified modes over the Green map
colorBy             = :type        # :type | :score | :confidence
showGreenBackground = false       # background drawn by plotGreenMap (same climRange) into the same axis
refLevel            = :full        # :sawOnly | :sawLimit | :full
figClass = Figure()
axClass  = Axis(figClass[1, 1]; title = "Cu (001): classification vs. |G_13^eff|")
plotGreenMap(sweep, scaleMode, climRange; ax = axClass)
plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel; ax = axClass)
ylims!(axClass, 1, 5.5)
display(figClass)
# switch without recomputing: reclassify(sweep, :low), classify(sweep, :ritz | :refined)

## All Ritz modes, shaded by a recognition criterion (diagnostic)
criterion       = :depthInsensitivity     # :ritzScore (default) | :freqDisp | :depthInsensitivity
velCutoffFactor = 1.10                    # x max longitudinal velocity; Inf = no cutoff
axRitz = plotRitzModeMap(sweep, criterion, velCutoffFactor)
ritzClim = (1.0, 8.0)              # colour limits of the criterion
for p in axRitz.scene.plots
    p.color[] isa AbstractVector{<:Real} && (p.colorrange[] = ritzClim)
end
for c in axRitz.parent.content
    c isa Colorbar && (c.limits[] = ritzClim)
end
axRitz.title = "Cu (001): all Ritz modes"
display(axRitz.parent)

## Any other component, recomputed on request (not stored) -- here |G_33|
component = "G33"                  # "G13eff" (default, stored) | "GIJ" | "GIJeff", e.g. "G13", "G33"
axG33, _ = plotGreenMap(sweep, scaleMode, climRange, plotTypeCartesian; component = component)
ylims!(axG33, 1, 5.5)
axG33.title = "Cu (001): |G_33| (recomputed)"
display(axG33.parent)
