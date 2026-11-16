## Tutorial: dispersion curves -- angle sweep and wavelength sweep, two viewers (Julia port)
#
# Same physics as examples/tutorials/tutorial_07_dispersion_map_and_viewer.m: two complementary
# real-velocity-units views of the SAME material pair (NOT the same underlying sweep type --
# see each Part's own header below):
#   Part 1 -- ONE fixed wavelength, MANY angles: SAW velocity vs. angle (plotAngularMap), plus
#             the combined viewer GUI.
#   Part 2 -- ONE fixed angle, MANY wavelengths: branch-tracked CURVES in the dispersion
#             viewer GUI, the wavelength-axis counterpart of Part 1.
#   Part 3 -- MANY angles x MANY wavelengths in ONE genuine double sweep: angular Green maps
#             sliced at fixed wavelengths, a dispersion curve sliced at a fixed angle, and the
#             double-sweep viewer switching between the two.
#
# Material pair -- Si(100) substrate / sputtered Cu film: the Cu is polycrystalline and assumed untextured, so
# isotropic (no crystal orientation -- see the Cu derivation below), and all angular variation
# comes from the cubic Si.
#
# Copy this script as a starting point for a new angle sweep (Part 1) or wavelength sweep
# (Part 2) on a film-on-substrate pair.
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - multiple return values used for real (not just destructured and ignored): `HillAverages`
#     returns a `(c11, c44)` tuple, and both are used immediately
#   - none else new -- this tutorial mainly exercises `AddLayer` (tutorial 5),
#     `ApplyWavelength`/`DispersionBranches.track`/`dispersionViewer` (tutorial 6), and adds
#     `combinedViewer` alongside them (same calling convention, no new concept)

using GuliwerrJ
using GLMakie   # opens an interactive window; `using CairoMakie` instead renders the same
                 # Figure/Axis objects headlessly (e.g. for save("fig.png", fig)) -- see
                 # tutorial_02_angular_map.jl for the full explanation

## Substrate -- Si(100), cubic, from the materials library (constants and density with
# their source: ListMaterialLibrary(print = true))
dirMatSi = [1 0 0;    # surface normal ([100])
            0 1 0]     # propagation reference ([010] -- arbitrary in-plane choice, cubic
                        # (100) has 4-fold symmetry
relDepthSub = 20.0   # substrate depth, in wavelengths; unlike the film's, it does not change with
                      # the swept wavelength

## Film -- Cu, isotropic (sputtered, so no crystal orientation and no dirMat): the Hill
# average of the single-crystal Cu constants (materials library, as in the other tutorials).
# Hill is the mean of the Voigt (uniform strain, upper) and Reuss (uniform stress, lower)
# bounds, the usual practical choice for a polycrystalline film; for Cu it gives a clearly
# softer shear constant (C44 47.0 vs 54.1 GPa Voigt), so it matters here. tutorial_S1
# compares the three averages.
matCuCrystal, _ = MaterialFromLibrary("Cu")
ecCuCubic     = matCuCrystal.elConsts
cuRho         = matCuCrystal.rho        # g/cm^3
cuC11Hill, cuC44Hill = HillAverages(ecCuCubic.voigt)   # -> (M, G), isotropic pair
# `HillAverages` returns TWO values as a tuple, destructured straight into two names here --
# the same pattern tutorial 2 used for `plotGreenMap`'s `(Axis, Colorbar)`. This replaces
# MATLAB's `[cuC11Hill, cuC44Hill] = ecCuCubic.GetHillAverages()`: the port keeps this
# computation as a plain function taking the Voigt matrix directly (`ecCuCubic.voigt`, a field
# access) rather than an object method, since `GetHillAverages` was itself only ever a
# one-line wrapper around the same free function in the MATLAB source.
filmThicknessUm = 3.0   # real film thickness [um] -- ApplyWavelength divides
                         # this by each swept wavelength (same unit) to derive relDepth

## Discretisation -- legOrder must stay >= 2x the LARGEST relDepth reached
# anywhere in the sweep. The largest is the substrate's 20, so legOrder = 50 clears the
# rule; the film reaches only relDepth 3/1.5 = 2 at the shortest wavelength (1.5 um).
legOrder = 50

matSub, _ = MaterialFromLibrary("Si"; dirMat = dirMatSi)
cpSub  = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)

ecFilm  = ElastConsts([cuC11Hill, cuC44Hill], 1)   # isotropic class -- no dirMat needed
matFilm = Material(; elConsts = ecFilm, rho = cuRho)
relDepthFilmPlaceholder = 1.0   # relDepth placeholder -- every point of the
                                  # sweep below overwrites it via ApplyWavelength
cpFilm = CalcParams(; relDepth = relDepthFilmPlaceholder, legOrder = legOrder,
                       realThickness = filmThicknessUm)

strName   = "si100_sub_cu_film"
nMaxModes = 500
str = Structure(; materIn = matSub, calcParamsIn = cpSub, name = strName, nMaxModes = nMaxModes)
str = AddLayer(str, matFilm, cpFilm)   # layer 2 = film, on top
# realThickness is set above via cpFilm, but ApplyWavelength itself is NOT
# called here -- wavevecAngSweep2Layered calls it internally, once per
# wavelength, on its own copy of str.

## =========================================================================================
## Part 1 -- fixed layer stack, SAW velocity vs. angle (angSweepLayered)
## =========================================================================================
#
# The layered counterpart of angSweepHS (tutorial 2): the film relDepth is fixed by one
# wavelength, and only the angle is swept. Part 1 uses the longest wavelength of the vector
# Part 2 sweeps.
angles      = 0.0:2.0:90.0                                       # deg
wavelengths = collect(1.5:0.25:12.0)   # um, evenly spaced; film relDepth 2 .. 0.25
# At relDepth 2 the film is thick enough that the substrate no longer shifts the Rayleigh
# wave (its velocity is within 1e-4 of a Cu half-space's).

## vMax = nothing (optional): the Green axis ends at 1.05 x the fastest longitudinal
# limiting velocity of any layer; a larger value only extends the displayed map. calcMode
# :full, not :ritzOnly: :ritzOnly has no Green map and reports trapped waves only (no leaky
# branch above the substrate's shear limit), while `classify(sweep, :ritz)` on a :full result
# still gives the Ritz-only entries at no extra cost.
vMax     = nothing
calcMode = :full

wavAngular = wavelengths[end]   # um; longest wavelength -> thinnest film (relDepth 0.25)
strAtWav   = ApplyWavelength(str, wavAngular)   # returns a NEW Structure, str itself untouched
retentionLevelAng = :saw_only   # enough: the viewer below recomputes a pruned angle's modal
                                  # shapes on demand
sweepAng = angSweepLayered(; structure = strAtWav, angles = angles, vMax = vMax,
                              calcMode = calcMode, retentionLevel = retentionLevelAng)
sweepAng = calc(sweepAng)

velLimsAng = (2.5, 3.0)    # km/s; the same range in the figure and the viewer below
axAng = plotAngularMap(sweepAng)   # creates its own figure/axes when none is passed
ylims!(axAng, velLimsAng...)
axAng.title = "Si(100)/$(filmThicknessUm) μm Cu film -- SAW velocity vs. angle, wavelength = $(round(wavAngular, digits = 1)) μm"
display(axAng.parent)

## Interactive companion to the static plot -- combinedViewer pairs modal shapes with the same
# angular map, click-to-inspect. Genuinely interactive only with a real GLMakie window (see
# tutorial 6's own note on this); the viewer object itself still builds under CairoMakie.
# `combinedViewer`, like `dispersionViewer` in tutorial 6, is a plain top-level GuliwerrJ
# function -- `using GuliwerrJ, GLMakie` is enough, no extra lookup needed.
guiAng = combinedViewer(sweepAng)
display(guiAng.fig)
display(guiAng.modalShapeViewer.shapeFig)   # the shape panels are always a window of their own
setAxisRange!(guiAng, :y, velLimsAng)   # sets the viewer's own range fields too
setAxisRange!(guiAng, :c, (-2.0, 0.0))   # colour: the two branches stand apart from the background
# Field names differ slightly from what you might guess by analogy with `DispersionViewer`
# (tutorial 6's `gui.mapAxes`) -- `AngSweepViewer` (the bottom pane here) calls its own axis
# just `ax`. When exploring a new viewer struct for the first time, `fieldnames(typeof(obj))`
# (or simply `dump(obj)`) lists every field name directly, the same job MATLAB's tab-completion
# on `gui.` does interactively.

## =========================================================================================
## Part 2 -- pick ONE angle from Part 1, follow it as branch-tracked curves
## =========================================================================================

angle          = 0.0   # deg
retentionLevel = :full   # REQUIRED: branch tracking compares every eigenvector of consecutive
                          # wavelength steps
calcModeCurve  = :full    # not :ritzOnly: each wavelength step keeps its Green map, needed
                           # for the Green background below
sweep1D = wavevecAngSweep2Layered(; structure = str, angles = [angle], wavelengths = wavelengths,
                                     vMax = vMax, calcMode = calcModeCurve,
                                     retentionLevel = retentionLevel)
sweep1D = calc(sweep1D)
# Expect two kinds of console warning from the Green solver here: a one-off "nearIsotropic"
# note (the Cu film is isotropic -- the solver switches to a closed-form shear basis, nothing
# to fix) and "nearSingular" messages for the thickest films (film relDepth up to 20): the
# boundary-condition matrix is badly conditioned there, so the Green map in the affected
# velocity band is less reliable. The Ritz branches are not affected.

branchModeNoMap = DispersionBranches.track(sweep1D)   # Ritz mode number of every branch at every step

## Launch the viewer -- one window (shapes+controls top, curve map bottom). Click a curve, or
# use the Wavelength/Branch steppers, to inspect that mode's own displacement field. Genuinely
# interactive only with a real GLMakie window.
gui = dispersionViewer(sweep1D)
display(gui.fig)

## Default the map to k-vphase with a fixed axis range (the viewer opens on wavelength-freq,
# autoscaled), over the Green background -- available for the two v_phase curve types on a
# :full sweep -- to see which tracked branch follows a Green peak. Setting these two
# Observables (curveTypeMenu.i_selected, showGreenBgToggle.active) already triggers the redraw
# -- unlike MATLAB's `gui.curveTypeDropdown.Value = ...; gui.onCurveTypeChanged();` two-step,
# there is no separate "changed" call to remember here (see tutorial 6's own note).
kVphaseIdx = findfirst(==(:k_vphase), visDispersionCurveTypes())
gui.curveTypeMenu.i_selected[] = kVphaseIdx
gui.showGreenBgToggle.active[] = true
kLims   = round.(2pi ./ (maximum(wavelengths), minimum(wavelengths)); sigdigits = 3)   # k [rad/um], the swept span
velLims = (1.0, 6.0)                            # km/s
setAxisRange!(gui, :x, kLims)   # as typed into the viewer's range fields
setAxisRange!(gui, :y, velLims)

## Static equivalent of the combined (branches + background) view above -- same data, one
# function call.
axStatic = plotDispersionCurves(sweep1D, :k_vphase, 1, true, true)   # creates its own figure/axes
xlims!(axStatic, kLims...)
ylims!(axStatic, velLims...)
axStatic.title = "Si(100)/$(filmThicknessUm) μm Cu film -- k-dispersion, angle = 0 deg, with Green's background"
display(axStatic.parent)

## "Show Ritz branches" is the paired toggle -- switch it off to see the Green background alone,
# with no line clutter over it (static equivalent: plotDispersionCurves(sweep1D, :k_vphase, 1,
# false, true)). The viewer's `mapAxes` limits are set once above and stay.
gui.showBranchesToggle.active[] = false

## ===================================================================
## Part 3 -- genuine double sweep: angle x wavelength in one calculation
## ===================================================================
#
# Parts 1 and 2 each fix one axis. wavevecAngSweep2Layered with BOTH several angles AND several
# wavelengths fills the whole grid at once: `sweep2D.sweeps[iw]` is a complete angular sweep
# (a CalcAngSweepLayered) at wavelength iw, so any single-wavelength plot from Part 1 works on a
# slice. RAM scales with nAngles x nWavelengths under retentionLevel = :full -- keep the grid
# coarse (docs/pitfalls.md NUM4); :full is needed here only for the dispersion-curve side.
anglesDouble      = 0.0:10.0:90.0                                        # deg, coarse
wavelengthsDouble = collect(1.5:0.5:12.0)   # um, coarse
retentionDouble   = :full   # :none is enough if only angular maps are wanted
sweep2D = wavevecAngSweep2Layered(; structure = str, angles = collect(anglesDouble),
                                     wavelengths = wavelengthsDouble, vMax = vMax,
                                     calcMode = calcMode, retentionLevel = retentionDouble)
sweep2D = calc(sweep2D)

## Angular |G_13,eff| maps at the shortest and longest wavelength (thickest vs. thinnest
# relative film) -- each slice is an ordinary angular sweep result. At 1.5 um (film 2
# wavelengths thick) the Rayleigh wave is that of a Cu half-space, flat in angle; the extra
# lines are standing waves across the film thickness. At 12 um the Si below shows through in
# the angle dependence.
iWavShort = 1
iWavLong  = length(wavelengthsDouble)
climRange = (-3.0, -1.0)   # one decade above the (-4, -1) default: the film branches stand out
velLimsGreen = (velLims, (1.0, nothing))   # km/s; the long wavelength reaches faster, Si-like
                                            # velocities (nothing = automatic upper limit)
for (iFig, iWav) in enumerate((iWavShort, iWavLong))
    axSlice, _ = plotGreenMap(sweep2D.sweeps[iWav], :log, climRange)
    ylims!(axSlice, velLimsGreen[iFig]...)
    axSlice.title = "Si(100)/$(filmThicknessUm) μm Cu film -- |G_13^eff|, wavelength = $(round(wavelengthsDouble[iWav], digits = 1)) μm"
    display(axSlice.parent)
end

## Dispersion curve at one fixed angle of the SAME grid -- angleIdx picks it
angleIdxDouble = findfirst(==(0.0), anglesDouble)   # index into anglesDouble, not an angle in deg
wavLimsDouble  = (wavelengthsDouble[1], wavelengthsDouble[end])   # um, the swept span
axDouble = plotDispersionCurves(sweep2D, :wavelength_vphase, angleIdxDouble, true, true, :log,
                                climRange)
xlims!(axDouble, wavLimsDouble...)
ylims!(axDouble, velLims...)
axDouble.title = "Si(100)/$(filmThicknessUm) μm Cu film -- wavelength dispersion at 0 deg, from the double sweep"
display(axDouble.parent)

## Interactive: the double-sweep viewer -- a button switches between the angular map (fixed
# wavelength, stepper) and the dispersion curve (fixed angle), remembering each side's axis
# ranges. Needs GLMakie (a window).
guiDouble = doubleSweepViewer(sweep2D)
display(guiDouble.fig)
# setAxisRange! acts on the mode currently shown, and each mode keeps its own ranges. The
# viewer opens on the angular map at the shortest wavelength: give it the ranges of that
# wavelength's Green-map figure above.
setAxisRange!(guiDouble, :y, velLimsGreen[1])
setAxisRange!(guiDouble, :c, climRange)
# Then open on the dispersion curve, as in the figure above (the button switches back to the
# angular map). The dispersion viewer has no colour control, so its Green background takes
# climRange directly on the heatmap (a display toggle there redraws it at the default range).
guiDouble.onSwitchMode!()
dvDouble = guiDouble.dispersionViewer[]
dvDouble.curveTypeMenu.i_selected[] = findfirst(==(:wavelength_vphase), visDispersionCurveTypes())
dvDouble.showGreenBgToggle.active[] = true
setAxisRange!(guiDouble, :x, wavLimsDouble)
setAxisRange!(guiDouble, :y, velLims)
for p in dvDouble.mapAxes.scene.plots
    p isa Heatmap && (p.colorrange[] = climRange)
end
