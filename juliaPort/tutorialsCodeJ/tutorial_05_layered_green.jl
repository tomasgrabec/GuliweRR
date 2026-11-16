## Tutorial 05 (code) -- layered :full angular sweep: Ni film on Cu(001), Green map + classification
# Narrated version: juliaPort/tutorialsJ/tutorial_05_layered_green.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_05_layered_green.jl")

using GuliwerrJ
using GLMakie

## Substrate -- Cu(001), cubic (materials library, see below)
dirMatSub   = [1 0 0; 0 0 1]       # surface normal; propagation reference

## Film -- Ni(001), cubic (materials library, see below)
dirMatFilm = [1 0 0; 0 0 1]

## Discretisation -- "easy" tier ((001) principal plane); film relDepth set directly (no ApplyWavelength)
relDepthSub  = 20.0                # substrate depth, wavelengths
relDepthFilm = 0.2                 # film thickness, wavelengths
legOrder     = 50

matSub, _ = MaterialFromLibrary("Cu"; dirMat = dirMatSub)   # constants + density from the materials library
matFilm, _ = MaterialFromLibrary("Ni"; dirMat = dirMatFilm)
cpSub   = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)
cpFilm  = CalcParams(; relDepth = relDepthFilm, legOrder = legOrder)

## Structure -- layer 1 = substrate, AddLayer puts the film on top
structureName = "cu_sub_ni_film"
nMaxModes     = 500
structure = Structure(; materIn = matSub, calcParamsIn = cpSub, name = structureName,
                        nMaxModes = nMaxModes)
structure = AddLayer(structure, matFilm, cpFilm)

## Angular sweep -- vMax is optional (nothing = Green axis ends at 1.05 x the fastest vL of any layer)
vMax           = nothing           # km/s; a larger value only extends the displayed map
anglesDeg      = 0.0:2.0:90.0
calcMode       = :full             # :full | :ritzOnly
retentionLevel = :full             # :full | :saw_only | :none
sweep = angSweepLayered(; structure = structure, angles = anglesDeg, vMax = vMax,
                          calcMode = calcMode, retentionLevel = retentionLevel)
sweep = calc(sweep)                # no nVelGreen for layered stacks (fixed 4 m/s axis)

## |G_13,eff| map with a widened colour range (component = "G33" recomputes |G_33|)
scaleMode = :log                   # :log | :linear
climRange = (-3.0, -0.5)           # log10(|G|/max|G|) colour limits (default (-4, -1))
axGreen, _ = plotGreenMap(sweep, scaleMode, climRange)
axGreen.title = "Cu / Ni film: |G_13^eff|"
display(axGreen.parent)

## Classification diagnostic over a plotGreenMap background with the same climRange
colorBy             = :type        # :type | :score | :confidence
showGreenBackground = false        # background drawn by plotGreenMap into the same axis
refLevel            = :full        # :sawOnly | :sawLimit | :full
figClass = Figure()
axClass  = Axis(figClass[1, 1]; title = "Cu / Ni film: classification vs. |G_13^eff|")
plotGreenMap(sweep, scaleMode, climRange; ax = axClass)
plotClassificationDiagnostic(sweep, colorBy, showGreenBackground, refLevel; ax = axClass)
display(figClass)
