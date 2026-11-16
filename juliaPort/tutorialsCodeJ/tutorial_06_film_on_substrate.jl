## Tutorial 06 (code) -- Cu film on Al: SAW velocity vs. film relDepth, one tracked branch
# Narrated version: juliaPort/tutorialsJ/tutorial_06_film_on_substrate.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_06_film_on_substrate.jl")

using GuliwerrJ
using GLMakie

## Substrate -- isotropic Al (faster than the Cu film)
matAlCrystal, _ = MaterialFromLibrary("Al")   # single crystal; isotropic Al = its Hill average
alC11, alC44 = HillAverages(matAlCrystal.elConsts.voigt)
alConsts    = [alC11, alC44]       # [C11, C44] GPa
isoSymmNo   = 1
alDensity   = matAlCrystal.rho     # g/cm^3
relDepthSub = 20.0                 # substrate depth, wavelengths

## Film -- Cu(001), cubic, real thickness in nm
dirMatFilm      = [1 0 0; 0 0 1]
filmThicknessNm = 200.0            # ApplyWavelength divides this by each wavelength

## Discretisation -- legOrder >= 2x the largest relDepth reached (3 here)
legOrder            = 50
relDepthPlaceholder = 1.0          # film relDepth, overwritten per wavelength

ecSub   = ElastConsts(alConsts, isoSymmNo)
matSub  = Material(; elConsts = ecSub, rho = alDensity)
matFilm, _ = MaterialFromLibrary("Cu"; dirMat = dirMatFilm)   # constants + density from the materials library
cpSub   = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)
cpFilm  = CalcParams(; relDepth = relDepthPlaceholder, legOrder = legOrder,
                       realThickness = filmThicknessNm)

structureName = "al_sub_cu_film"
nMaxModes     = 500
structure = Structure(; materIn = matSub, calcParamsIn = cpSub, name = structureName,
                        nMaxModes = nMaxModes)
structure = AddLayer(structure, matFilm, cpFilm)

## Wavelength sweep at one angle -- dense steps help branch tracking
filmRelDepths  = vcat(0.01:0.01:1.0, 1.05:0.05:3.0)   # starts at relDepth 0.01
wavelengthsNm  = filmThicknessNm ./ filmRelDepths   # descending
angleDeg       = 0.0
calcMode       = :ritzOnly         # vMax not needed for :ritzOnly
retentionLevel = :full             # REQUIRED by DispersionBranches.track
sweep = wavevecAngSweep2Layered(; structure = structure, angles = [angleDeg],
                                  wavelengths = wavelengthsNm, calcMode = calcMode,
                                  retentionLevel = retentionLevel)
sweep = calc(sweep)

## Branch tracking by eigenvector overlap, seeded at the thinnest film
vSAWTracked, _, _ = DispersionBranches.trackSeededBranch(sweep)   # seed: slowest SAW at the thinnest film

## Half-space limits -- each material alone
calcModeHS  = :ritzOnly
cpHalfSpace = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)   # each material as its own half-space
structAlHS  = Structure(; materIn = matSub, calcParamsIn = cpHalfSpace, name = "al_hs",
                          nMaxModes = nMaxModes)
singleAl    = calcAll(singleHS(; structure = structAlHS, angle = angleDeg), calcModeHS)
smAl        = singleAl.classification.significantModes
vSAWAl      = minimum(m.vel for m in smAl if m.type === :Rayleigh)
structCuHS  = Structure(; materIn = matFilm, calcParamsIn = cpHalfSpace, name = "cu_hs",
                          nMaxModes = nMaxModes)
singleCu    = calcAll(singleHS(; structure = structCuHS, angle = angleDeg), calcModeHS)
smCu        = singleCu.classification.significantModes
vSAWCu      = minimum(m.vel for m in smCu if m.type === :Rayleigh)

## Plot 1 -- tracked branch between the two half-space limits
figTracked = Figure()
axTracked  = Axis(figTracked[1, 1]; xscale = log10, xlabel = "film relDepth (∝ kd)",
                  ylabel = "v_SAW [km/s]")
scatterlines!(axTracked, filmRelDepths, vSAWTracked)
hlines!(axTracked, [vSAWAl, vSAWCu]; linestyle = :dash)   # Al and Cu half-space values
display(figTracked)

## Plot 2 -- every tracked Ritz branch, SAW branches highlighted
curveType = :k_vphase              # :wavelength_freq | :wavelength_vphase | :k_freq | :k_vphase | :freq_vphase | :freq_vgroup
kAxisLimits   = (0.0, 0.03)
velAxisLimits = (1.0, 4.0)         # km/s
axCurves = plotDispersionCurves(sweep, curveType)
xlims!(axCurves, kAxisLimits...); ylims!(axCurves, velAxisLimits...)
display(axCurves.parent)

## Interactive dispersion viewer, preset to k-vphase
gui = dispersionViewer(sweep)
display(gui.fig)
gui.curveTypeMenu.i_selected[] = findfirst(==(curveType), visDispersionCurveTypes())
setAxisRange!(gui, :x, kAxisLimits); setAxisRange!(gui, :y, velAxisLimits)

## Optional: :full mode -- tracked branches over the layered |G13_eff| map
sweepFull = wavevecAngSweep2Layered(; structure = structure, angles = [angleDeg],
                                      wavelengths = wavelengthsNm, calcMode = :full,
                                      retentionLevel = retentionLevel)
sweepFull = calc(sweepFull)
angleIdx = 1; showRitzBranches = true; showGreenBackground = true
scaleMode = :log; climRange = (-3.0, -0.5)
axFull = plotDispersionCurves(sweepFull, curveType, angleIdx, showRitzBranches,
                              showGreenBackground, scaleMode, climRange)
axFull.title = "Al / Cu film: tracked branches over |G_13^eff|"
xlims!(axFull, kAxisLimits...); ylims!(axFull, velAxisLimits...)
display(axFull.parent)

# dispersion viewer = the combined viewer's counterpart for wavelength sweeps
guiFull = dispersionViewer(sweepFull)
display(guiFull.fig)
guiFull.curveTypeMenu.i_selected[] = findfirst(==(curveType), visDispersionCurveTypes())
guiFull.showGreenBgToggle.active[] = true
setAxisRange!(guiFull, :x, kAxisLimits); setAxisRange!(guiFull, :y, velAxisLimits)
