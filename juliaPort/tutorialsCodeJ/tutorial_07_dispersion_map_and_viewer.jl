## Tutorial 07 (code) -- Cu film on Si(100): angle sweep, wavelength sweep, double sweep + viewers
# Narrated version: juliaPort/tutorialsJ/tutorial_07_dispersion_map_and_viewer.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_07_dispersion_map_and_viewer.jl")

using GuliwerrJ
using GLMakie

## Substrate -- Si(100), cubic (materials library)
dirMatSi    = [1 0 0; 0 1 0]       # (100) surface, [010] reference
relDepthSub = 20.0                 # not swept by ApplyWavelength

## Film -- sputtered Cu, isotropic via the Hill average of cubic Cu (materials library)
matCuCrystal, _ = MaterialFromLibrary("Cu")   # single crystal
cuDensity       = matCuCrystal.rho            # g/cm^3
isoSymmNo       = 1
cuC11Hill, cuC44Hill = HillAverages(matCuCrystal.elConsts.voigt)   # isotropic (M, G) = (C11, C44)
filmThicknessUm = 3.0               # um, same unit as the wavelengths

## Discretisation -- substrate relDepth 20 (film at most 3/1.5 = 2) -> legOrder 50
legOrder            = 50
relDepthPlaceholder = 1.0          # film relDepth, overwritten per wavelength

ecFilm  = ElastConsts([cuC11Hill, cuC44Hill], isoSymmNo)
matSub, _ = MaterialFromLibrary("Si"; dirMat = dirMatSi)
matFilm = Material(; elConsts = ecFilm, rho = cuDensity)
cpSub   = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)
cpFilm  = CalcParams(; relDepth = relDepthPlaceholder, legOrder = legOrder,
                       realThickness = filmThicknessUm)

structureName = "si100_sub_cu_film"
nMaxModes     = 500
structure = Structure(; materIn = matSub, calcParamsIn = cpSub, name = structureName,
                        nMaxModes = nMaxModes)
structure = AddLayer(structure, matFilm, cpFilm)

## Shared axes and vMax
anglesDeg     = 0.0:2.0:90.0
wavelengthsUm = collect(1.5:0.25:12.0)   # um, evenly spaced; film relDepth 2 .. 0.25
vMax     = nothing                 # optional: nothing = Green axis ends at 1.05 x the fastest vL of any layer
calcMode = :full                   # :ritzOnly: no Green map, trapped waves only

## Part 1 -- one wavelength, many angles
wavelengthAngular = wavelengthsUm[end]                          # um, thinnest relative film
structureAtWav    = ApplyWavelength(structure, wavelengthAngular)   # returns a NEW Structure
retentionAngular  = :saw_only      # viewer recomputes pruned angles on demand
sweepAngular = angSweepLayered(; structure = structureAtWav, angles = anglesDeg, vMax = vMax,
                                 calcMode = calcMode, retentionLevel = retentionAngular)
sweepAngular = calc(sweepAngular)
velLimsAngular = (2.5, 3.0)         # km/s, figure and viewer alike
axAngular = plotAngularMap(sweepAngular)
ylims!(axAngular, velLimsAngular...)
axAngular.title = "Si(100)/$(filmThicknessUm) μm Cu: SAW velocity vs. angle, $(round(wavelengthAngular, digits = 1)) μm"
display(axAngular.parent)          # plot functions and viewers only build figures; display shows them

guiAngular = combinedViewer(sweepAngular)
display(guiAngular.fig); display(guiAngular.modalShapeViewer.shapeFig)
setAxisRange!(guiAngular, :y, velLimsAngular)
setAxisRange!(guiAngular, :c, (-2.0, 0.0))

## Part 2 -- one angle, many wavelengths, branch-tracked
angleDeg       = 0.0
retentionCurve = :full             # REQUIRED by DispersionBranches.track
sweepCurve = wavevecAngSweep2Layered(; structure = structure, angles = [angleDeg],
                                       wavelengths = wavelengthsUm, vMax = vMax,
                                       calcMode = calcMode, retentionLevel = retentionCurve)
sweepCurve = calc(sweepCurve)

curveType           = :k_vphase    # 6 options, see visDispersionCurveTypes()
angleIdx            = 1            # index into the sweep's angles
showRitzBranches    = true
showGreenBackground = true         # needs calcMode = :full and a velocity y-axis
kAxisLimits         = round.(2pi ./ (maximum(wavelengthsUm), minimum(wavelengthsUm)); sigdigits = 3)   # rad/um
velAxisLimits       = (1.0, 6.0)   # km/s

gui = dispersionViewer(sweepCurve)
display(gui.fig)
gui.curveTypeMenu.i_selected[]  = findfirst(==(curveType), visDispersionCurveTypes())
gui.showGreenBgToggle.active[]  = showGreenBackground
setAxisRange!(gui, :x, kAxisLimits); setAxisRange!(gui, :y, velAxisLimits)

## Static equivalent -- branches over the Green background
axCurve = plotDispersionCurves(sweepCurve, curveType, angleIdx, showRitzBranches, showGreenBackground)
xlims!(axCurve, kAxisLimits...); ylims!(axCurve, velAxisLimits...)
display(axCurve.parent)

## Part 3 -- genuine double sweep: angle x wavelength in one calculation
anglesDouble      = 0.0:10.0:90.0                                          # coarse: RAM scales with the grid
wavelengthsDouble = collect(1.5:0.5:12.0)   # um
retentionDouble   = :full          # :none suffices for angular maps only
sweepDouble = wavevecAngSweep2Layered(; structure = structure, angles = collect(anglesDouble),
                                        wavelengths = wavelengthsDouble, vMax = vMax,
                                        calcMode = calcMode, retentionLevel = retentionDouble)
sweepDouble = calc(sweepDouble)

## Angular |G_13,eff| map at one wavelength -- each sweeps[iWav] is a full angular sweep
scaleMode = :log
climRange = (-3.0, -1.0)
iWavShort = 1
axSlice, _ = plotGreenMap(sweepDouble.sweeps[iWavShort], scaleMode, climRange)
ylims!(axSlice, velAxisLimits...)
axSlice.title = "Si(100)/$(filmThicknessUm) μm Cu: |G_13^eff| at $(round(wavelengthsDouble[iWavShort], digits = 1)) μm"
display(axSlice.parent)

## Dispersion curve at one angle of the same grid
curveTypeDouble = :wavelength_vphase
angleIdxDouble  = findfirst(==(0.0), anglesDouble)
axDouble = plotDispersionCurves(sweepDouble, curveTypeDouble, angleIdxDouble, showRitzBranches,
                                showGreenBackground, scaleMode, climRange)
wavAxisLimits = (wavelengthsDouble[1], wavelengthsDouble[end])   # um, the swept span
xlims!(axDouble, wavAxisLimits...); ylims!(axDouble, velAxisLimits...)
display(axDouble.parent)

## Interactive double-sweep viewer -- toggle angular map <-> dispersion curve
guiDouble = doubleSweepViewer(sweepDouble)
display(guiDouble.fig)
setAxisRange!(guiDouble, :y, velAxisLimits)   # angular map at iWavShort, as its figure
setAxisRange!(guiDouble, :c, climRange)
guiDouble.onSwitchMode!()                     # open on the dispersion curve, same ranges as above
dvDouble = guiDouble.dispersionViewer[]
dvDouble.curveTypeMenu.i_selected[] = findfirst(==(curveTypeDouble), visDispersionCurveTypes())
dvDouble.showGreenBgToggle.active[] = showGreenBackground
setAxisRange!(guiDouble, :x, wavAxisLimits); setAxisRange!(guiDouble, :y, velAxisLimits)
for p in dvDouble.mapAxes.scene.plots         # no colour control in the viewer
    p isa Heatmap && (p.colorrange[] = climRange)
end
