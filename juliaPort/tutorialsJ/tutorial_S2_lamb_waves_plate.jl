## Tutorial: free plate -- Lamb-wave dispersion, mode types and ZGV points (Julia port)
#
# Same physics as examples/tutorials/tutorial_S2_lamb_waves_plate.m. Every other tutorial in
# this series models a SURFACE wave on a semi-infinite substrate. This one models a
# free-standing PLATE: both faces traction-free, finite total thickness -- the Lamb-wave
# problem (Grabec, Sedlak, Seiner 2020, Wave Motion 96, 102567). Same Ritz method, different
# boundary conditions, so the family is deliberately leaner (docs/reference/architecture/plate.md):
#   - NO Green's function and NO calcMode argument -- only the Ritz solve runs;
#   - every computed mode is a genuine physical branch (no "which mode is the SAW" question) --
#     instead each mode gets a TYPE: Asymm (A0-like), Symm (S0-like), TH-Asymm, TH-Symm
#     (shear-horizontal-like);
#   - retentionLevel is :full | :none only (no :saw_only tier).
#
# Physical scaling: pass the plate's real thickness and a real wavevector k in CONSISTENT
# units -- here mm and rad/mm. Velocities stay in km/s (= mm/us), so frequency f = v*k/(2*pi)
# comes out in MHz.
#
#   Part 1 -- one angle, one wavevector: the mode list itself
#   Part 2 -- one angle, many wavevectors: dispersion curves, branch classification,
#             zero-group-velocity (ZGV) points
#   Part 3 -- many angles, one wavevector: angular map + combined viewer
#   Part 4 -- a coated (two-layer, NON-mirror-symmetric) plate: when mode-type labels stop
#             being rigorous
#   Part 5 -- a measured reference: the S1 ZGV resonances of a Si(001) wafer, compared with
#             published measurements
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series. Run
# with `-t auto` for the threaded sweeps to use more than one core.
#
# Julia concepts introduced in this tutorial:
#   - indexing a Tuple of names with an integer code (`modeTypeNames[code]`) -- identical to
#     MATLAB's cell-array `{code}` lookup, with round brackets replaced by square ones
#   - `mean` from the Statistics standard library (MATLAB has it built in)
#   - iterating over a Vector of structs with `for z in zgvList` -- no index variable needed

using GuliwerrJ
using GLMakie   # opens an interactive window; `using CairoMakie` instead renders the same
                 # Figure/Axis objects headlessly -- see tutorial_02_angular_map.jl
using Printf
using Statistics   # for `mean`

## Material -- Cu(001), cubic: an anisotropic plate with ZGV points for Part 2 to find.
# Constants and density come from the materials library (`ListMaterialLibrary(; print = true)`).
dirMat   = [1 0 0;    # plate normal
            0 0 1]     # in-plane propagation reference (Material's default)
mat, _   = MaterialFromLibrary("Cu"; dirMat = dirMat)

## Discretisation -- for a plate, relDepth is the plate thickness in wavelengths
# (k*d/(2*pi)), and it is OVERWRITTEN at every point from realThickness and realWavVec below --
# the value here is only a placeholder. legOrder must resolve the THICKEST relDepth reached:
# kMax*d/(2*pi) = 20/(2*pi) ~ 3.2 here; legOrder = 20 and 30 give the same branches below
# 6 MHz (Part 2), so 20 is ample there (a plate is thin in wavelengths, unlike a
# 20-wavelength-deep half-space).
relDepthPlaceholder = 0.2
legOrder = 20   # number of Legendre basis functions across the thickness
cp = CalcParams(; relDepth = relDepthPlaceholder, legOrder = legOrder)

strName   = "cu_001_plate"
nMaxModes = 500
str = Structure(; materIn = mat, calcParamsIn = cp, name = strName, nMaxModes = nMaxModes)

realThickness = 1.0   # mm -- plate thickness, FIXED throughout this script

## ===================================================================
## Part 1 -- one angle, one wavevector: the mode list
## ===================================================================
angle          = 0.0     # deg
realWavVec     = 2.0     # rad/mm -> wavelength 2*pi/k ~ 3.1 mm, k*d = 2
retentionLevel = :full   # :full (keep eigenvectors: modal shapes, tracking) | :none

singleResult = singlePlate(; structure = str, angle = angle, realThickness = realThickness,
                             realWavVec = realWavVec)
singleResult = calcAll(singleResult, retentionLevel)
ritz = singleResult.ritzSolver

modeTypeNames = ("Asymm", "Symm", "TH-Asymm", "TH-Symm")   # modesType codes 1..4
nShow = 8
@printf("\nCu(001) plate, d = %g mm, k = %g rad/mm, angle = %g deg\n", realThickness, realWavVec, angle)
@printf("%5s %10s %10s %10s %10s %11s\n", "mode", "v [km/s]", "vg [km/s]", "f [MHz]", "type", "confidence")
for iMode in 1:nShow
    @printf("%5d %10.4f %10.4f %10.4f %10s %11.2f\n", iMode, ritz.vels[iMode], ritz.velsGroup[iMode],
            ritz.freqsPhysical[iMode], modeTypeNames[ritz.modesType[iMode]],
            ritz.modeTypeConfidence[iMode])
end
# The type comes from mid-plane parity: each displacement component is split into its even and
# odd parts about the plate's mid-plane, and the mode goes to the class holding most of its
# energy (Symm: u1 even + u3 odd, Asymm: u1 odd + u3 even, TH-Symm: u2 even, TH-Asymm: u2 odd).
# modeTypeConfidence = largest fraction - runner-up, in [0,1]: near 1 = pure (every mode here,
# along the cubic axis), low = mixed: SH-Lamb coupling within one parity class off the axes,
# parity mixing on an asymmetric stack (Part 4).

## 1D view: every mode at this one point, coloured by type
criterion       = :modeType   # :modeType (categorical) | :confidence (continuous)
velCutoffFactor = 2.5         # x the longitudinal velocity: the 8 modes of the table above
                              # (default 1.1, a surface-wave habit that hides most plate modes)
axModes = plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor)
axModes.title = "Cu(001) plate, k = $realWavVec rad/mm: Ritz modes by type"
display(axModes.parent)   # a plot function only builds the figure; `display` shows it (tutorial 2)

## ===================================================================
## Part 2 -- one angle, many wavevectors: dispersion curves and ZGV
## ===================================================================
# wavevecAngSweep2Plate sweeps the REAL wavevector at one fixed realThickness (no
# ApplyWavelength step, unlike the layered family). Dense steps help branch tracking, which
# matches CONSECUTIVE steps by eigenvector overlap (MAC). The span matters as much as the
# density: a branch is drawn only over the k range swept, so the low branches need k up to
# ~20 rad/mm to reach 6 MHz in the group-velocity plot.
wavVecs        = collect(0.1:0.1:20.0)   # rad/mm
retentionSweep = :full    # REQUIRED for branch tracking, classifyBranches and locateZGVs
useParallel    = :auto    # :auto | :always | :never -- forwarded to each angle loop
sweepDisp = wavevecAngSweep2Plate(; structure = str, angles = [angle], wavVecs = wavVecs,
                                     realThickness = realThickness,
                                     retentionLevel = retentionSweep, useParallel = useParallel)
sweepDisp = calc(sweepDisp)

## Dispersion curves -- one tracked line per branch. curveType picks the axes:
# :wavelength_freq (default) | :wavelength_vphase | :k_freq | :k_vphase | :freq_vphase |
# :freq_vgroup ('k' is wavVecs directly). Note the underscores: MATLAB's 'k-vphase' becomes the
# Symbol :k_vphase, since a `:name` Symbol literal cannot contain a hyphen.
angleIdx              = 1       # which angle of the sweep (only one here)
colorMode             = :type   # :type (mode type where trusted) | :branch (branch identity)
confidenceThresh      = 0.3     # default; below it a label falls back to "branch N" (unvalidated)
showChristoffelBounds = true    # bulk shear/longitudinal reference lines (velocity y-axes only)
velLimsDisp           = (0.0, 8.0)   # km/s, phase-velocity axis (here and in the viewer below)
freqLimsDisp          = (0.0, 6.0)   # MHz; higher branches are less resolved and crowd the plot
axVphase = plotDispersionCurvesPlate(sweepDisp, :k_vphase, angleIdx, colorMode, confidenceThresh,
                                     showChristoffelBounds)
ylims!(axVphase, velLimsDisp...)
axVphase.xlabel = "k [rad/mm]"   # the plot cannot know the length unit: units follow realThickness and wavVecs
axVphase.title = "Cu(001) plate: phase velocity vs. k"
display(axVphase.parent)

axFreq = plotDispersionCurvesPlate(sweepDisp, :k_freq, angleIdx, colorMode, confidenceThresh)
ylims!(axFreq, freqLimsDisp...)
axFreq.xlabel = "k [rad/mm]"
axFreq.ylabel = "frequency [MHz]"   # km/s x rad/mm / (2*pi) = MHz
axFreq.title = "Cu(001) plate: frequency vs. k"
display(axFreq.parent)

axGroup = plotDispersionCurvesPlate(sweepDisp, :freq_vgroup, angleIdx, colorMode, confidenceThresh)
xlims!(axGroup, freqLimsDisp...)
axGroup.xlabel = "frequency [MHz]"
hlines!(axGroup, [0.0]; color = :black, linestyle = :dot)   # zero group velocity at finite k = ZGV
axGroup.title = "Cu(001) plate: group velocity vs. frequency"
display(axGroup.parent)

## Branch classification -- Cu(001) along [100] is symmetric under both the mid-plane mirror
# (exact parity: Symm/Asymm) and the sagittal mirror (exact Lamb/TH split), so every label is
# rigorous here. Off the symmetry axes (e.g. 30 deg) only the parity is exact, and a cut like
# (111) has neither; Part 4 shows a stack without the mid-plane mirror.
classResult = classifyBranches(sweepDisp, angleIdx, confidenceThresh)
@printf("\nMid-plane mirror: %s, sagittal mirror: %s (both: labels rigorous regardless of confidence)\n",
        classResult.isExactMirrorSymmetric, classResult.isExactSagittalSymmetric)
println("Branches 1-6 at the first k: ", join(classResult.displayLabel[1, 1:6], ", "))

## Zero-group-velocity points -- a scan along each tracked branch for a sign change of the
# group velocity, then a Newton refinement that solves the plate eigenproblem and the
# zero-group-velocity condition together (no mode number is held fixed, so a branch crossing
# nearby does not matter). A point is kept only where an independent solve confirms that the
# group velocity vanishes; MATLAB returns the same list.
maxCount = Inf   # cap on returned points (Inf = all)
zgvList  = locateZGVs(sweepDisp, angleIdx, maxCount)
zgvShow  = filter(z -> z.freqZGV <= freqLimsDisp[2], zgvList)   # the plotted frequency range
@printf("\nZGV points found: %d, of them below %g MHz:\n", length(zgvList), freqLimsDisp[2])
for z in zgvShow
    @printf("  branch %2d: k = %.4f rad/mm, f = %.4f MHz, true ZGV: %s\n", z.branchNo, z.wavVecZGV,
            z.freqZGV, z.isTrueZGV)
end
# Check higher-branch points for convergence (raise legOrder, refine wavVecs) before trusting
# them: high-order modes are the least resolved. Part 5 checks the same search against
# published measurements on a silicon wafer.
#
# The search runs along k at ONE fixed angle, so it finds where the group velocity along k
# vanishes. Along a symmetry direction of the plate ([100] here) the group velocity is parallel
# to k, and that is a true ZGV point. Along any other direction the group velocity also has a
# component perpendicular to k (velGroupTrans, from ritzSolver.velsGroupTrans), which need not
# vanish there; isTrueZGV says whether it does (Part 5).

## Interactive: dispersion viewer -- click a curve (or step wavevector/branch) to see that
# mode's displacement field through the thickness. The same `dispersionViewer` function as the
# layered tutorials: multiple dispatch picks the plate method from the argument's type.
guiDisp = dispersionViewer(sweepDisp)
display(guiDisp.fig)
# the viewer opens on :wavelength_freq, autoscaled: switch to the k-vphase view of the first
# figure (setting the Observable already redraws), then fix its range
guiDisp.curveTypeMenu.i_selected[] = findfirst(==(:k_vphase), visDispersionCurveTypes())
setAxisRange!(guiDisp, :y, velLimsDisp)   # as typed into the viewer's range fields

## ===================================================================
## Part 3 -- many angles, one wavevector: angular dependence
## ===================================================================
angles = 0.0:2.0:90.0   # deg, [100] -> [010] on the (001) plate
sweepAng = angSweepPlate(; structure = str, angles = angles, realThickness = realThickness,
                            realWavVec = realWavVec, retentionLevel = retentionLevel,
                            useParallel = useParallel)
sweepAng = calc(sweepAng)
# sweepAng.velsMap / .velsGroupMap / .modesTypeMap / .freqsPhysicalMap are dense
# [nModes x nAngles] matrices -- every angle has the same mode count.

# Mode type and its confidence side by side in one window: pass each plot an Axis of one
# Figure instead of letting it open its own (the confidence colour bar goes beside its map)
velLimsAng = (1.0, 5.5)   # km/s; the same range in both maps and the viewer below
figMaps   = Figure(; size = (1400, 520))
axMapType = Axis(figMaps[1, 1])
plotRitzModeMapPlate(sweepAng, :modeType; ax = axMapType)
ylims!(axMapType, velLimsAng...)
axMapType.title = "Cu(001) plate, k = $realWavVec rad/mm: by type"
axMapConf = Axis(figMaps[1, 2])
plotRitzModeMapPlate(sweepAng, :confidence; ax = axMapConf)
ylims!(axMapConf, velLimsAng...)
axMapConf.title = "mode-type confidence"
display(figMaps)

## Interactive: combined viewer -- angular map (bottom) linked to the modal shape (top).
# Accepts Part 1's single-angle `singleResult` too (dispatch again).
guiAng = combinedViewer(sweepAng)
display(guiAng.fig)
setAxisRange!(guiAng, :y, velLimsAng)   # sets the viewer's own range fields too
display(guiAng.modalShapeViewer.shapeFig)   # the shape panels are always a window of their own

## ===================================================================
## Part 4 -- coated plate: a thin isotropic Al layer on ONE face
## ===================================================================
# A different material on one face removes the mid-plane mirror symmetry, so symmetric/
# antisymmetric motion mixes and a mode-type label only names the dominant parity class,
# trusted where confidence >= confidenceThresh (a wave bound to one face is half even, half
# odd: confidence ~0). On four multi-layer test stacks the frequencies, group velocities and
# ZGV points agree with GEWtool (Kiefer, open source) to round-off; the confidenceThresh
# default is not validated (docs/reference/architecture/plate.md).
#
# realThickness is the TOTAL thickness; every layer's relDepth is scaled by the same factor, so
# the CalcParams relDepths below set the THICKNESS RATIO.
# Isotropic (polycrystalline) Al = the Hill average of the single-crystal constants from the
# materials library (tutorial_S1 section 9 shows the averages).
matAlCrystal, _ = MaterialFromLibrary("Al")
alC11, alC44 = HillAverages(matAlCrystal.elConsts.voigt)
alConsts = [alC11, alC44]   # [C11, C44] GPa, isotropic
alSymm   = 1
alRho    = matAlCrystal.rho   # g/cm^3
matAl    = Material(; elConsts = ElastConsts(alConsts, alSymm), rho = alRho)

coatingFraction = 0.1      # coating = 10% of the total thickness
cpCu = CalcParams(; relDepth = (1 - coatingFraction) * relDepthPlaceholder, legOrder = legOrder)
cpAl = CalcParams(; relDepth = coatingFraction * relDepthPlaceholder, legOrder = legOrder)
strCoated = Structure(; materIn = mat, calcParamsIn = cpCu, name = "cu_plate_al_coating",
                        nMaxModes = nMaxModes)          # layer 1 = bottom face
strCoated = AddLayer(strCoated, matAl, cpAl)            # layer 2 = top face

@printf("\nCoated plate mirror-symmetric: %s\n", isMirrorSymmetric(strCoated))

sweepCoated = wavevecAngSweep2Plate(; structure = strCoated, angles = [angle], wavVecs = wavVecs,
                                       realThickness = realThickness,
                                       retentionLevel = retentionSweep, useParallel = useParallel)
sweepCoated = calc(sweepCoated)

classCoated = classifyBranches(sweepCoated, angleIdx, confidenceThresh)
@printf("Coated plate: %.0f%% of (step, branch) labels trusted at confidenceThresh = %g\n",
        100 * mean(classCoated.trustLabelByBranch), confidenceThresh)

# :type colours by trusted label; :branch is the escape hatch when labels are not trusted
axCoatedType = plotDispersionCurvesPlate(sweepCoated, :k_vphase, angleIdx, :type, confidenceThresh)
ylims!(axCoatedType, velLimsDisp...)
axCoatedType.xlabel = "k [rad/mm]"
axCoatedType.title = "Cu plate + 10% Al coating: phase velocity vs. k (coloured by type)"
display(axCoatedType.parent)

axCoatedBranch = plotDispersionCurvesPlate(sweepCoated, :k_vphase, angleIdx, :branch)
ylims!(axCoatedBranch, velLimsDisp...)
axCoatedBranch.xlabel = "k [rad/mm]"
axCoatedBranch.title = "Cu plate + 10% Al coating: phase velocity vs. k (coloured by branch)"
display(axCoatedBranch.parent)

## ===================================================================
## Part 5 -- a measured reference: S1 ZGV resonances of a Si(001) wafer
## ===================================================================
# Prada, Clorennec, Murray & Royer (2009, J. Acoust. Soc. Am. 126, 620) and Kiefer, Mezil &
# Prada (2023, Sci. Adv. 9, eadk6846) excited a [001]-cut silicon wafer with a pulsed laser and
# measured the S1 ZGV resonances with an interferometer. Kiefer et al. give everything needed to
# repeat the calculation: thickness 524.6 um, C11 = 165.6, C12 = 63.9, C44 = 79.5 GPa, density
# 2330 kg/m^3. The ZGV points sit on the principal directions:
#   ZGV1 on <100>: computed 7.7079 MHz at k = 3.421 rad/mm, measured 7.7072 MHz
#   ZGV2 on <110>: computed 7.7551 MHz at k = 3.142 rad/mm, measured 7.7545 MHz
# (computed with their own semi-analytical code, GEWtool). Prada et al. give the same points as
# f*d = 4.044 and 4.069 MHz*mm (computed with DISPERSE).
# Kiefer et al.'s own Si constants, kept on purpose to repeat their calculation; the
# materials library's Si (McSkimin & Andreatch 1964) differs by about 0.1 %.
siConsts = [165.6, 63.9, 79.5]   # [C11, C12, C44] GPa
siRho    = 2.330                  # g/cm^3
siDirMat = [0 0 1;                # plate normal [001]
            1 0 0]                # angle 0 = [100], angle 45 = [110]
matSi = Material(; elConsts = ElastConsts(siConsts, 2), rho = siRho, dirMat = siDirMat)
strSi = Structure(; materIn = matSi, calcParamsIn = cp, name = "si_001_wafer", nMaxModes = nMaxModes)

siThickness = 0.5246                                  # mm
siAngles    = [0.0, 45.0, 18.4]                       # deg: [100], [110], 26.6 deg from [110]
siWavVecs   = collect(range(1.0, 6.0; length = 51))   # rad/mm, around the S1 minimum near 3.3 rad/mm
sweepSi = wavevecAngSweep2Plate(; structure = strSi, angles = siAngles, wavVecs = siWavVecs,
                                   realThickness = siThickness,
                                   retentionLevel = retentionSweep, useParallel = useParallel)
sweepSi = calc(sweepSi)

publishedComputed = [7.7079, 7.7551, 7.727]   # MHz, Kiefer et al. 2023
publishedMeasured = [7.7072, 7.7545, NaN]     # MHz, Kiefer et al. 2023, Fig. 8
@printf("\nSi(001) wafer, d = %.4f mm: lowest point of zero group velocity along k, per direction\n", siThickness)
@printf("%9s %12s %12s %12s %12s %6s %16s\n", "direction", "k [rad/mm]", "f [MHz]", "published",
        "measured", "type", "c_gtheta [mm/us]")
dirNames = ("[100]", "[110]", "18.4 deg")
for iAng in eachindex(siAngles)
    zgvSi = locateZGVs(sweepSi, iAng, Inf)
    zLow  = zgvSi[argmin([z.freqZGV for z in zgvSi])]   # the S1 point; the A2 point is near 15 MHz
    @printf("%9s %12.4f %12.4f %12.4f %12.4f %6s %16.4f\n", dirNames[iAng], zLow.wavVecZGV, zLow.freqZGV,
            publishedComputed[iAng], publishedMeasured[iAng], zLow.isTrueZGV ? "ZGV" : "TGV",
            zLow.velGroupTrans)
end
# GuliweRR reproduces the published ZGV frequencies to every digit given, the TGV frequency to 1
# kHz, and the measured ones to 0.01 %. Prada et al. used 2329 kg/m^3, which raises every
# frequency by sqrt(2330/2329), about 0.02 %; that is how their f*d = 4.044 and 4.069 MHz*mm
# come out.
#
# The two directions differ by about 50 kHz. ZGV1 is a minimum of the dispersion surface
# frequency(kx, ky); ZGV2 is a saddle point of it: a maximum across the angle, a minimum along
# k. Between the two axes, locateZGVs still finds a minimum along k at every angle, but there
# the group velocity has a component perpendicular to k, so the energy does not stay put:
# isTrueZGV is false. Kiefer et al. call these TGV (transverse-group-velocity) waves and
# measured them leaving a line source sideways. At 26.6 deg from [110] they computed 7.727 MHz
# and a transverse group velocity of -0.177 mm/us, and measured +-0.170 mm/us. Their angle runs
# from [110] towards [010], so that direction is 71.6 deg in GuliweRR's angle; 18.4 deg, used
# here, is its mirror image across [110]: an equivalent wave with the transverse component
# reversed. GuliweRR gives +0.178 mm/us here and -0.178 mm/us at 71.6 deg. Only the eight points
# on the principal directions are true ZGV resonances, and only they persist after about 40 us
# in the measurement.
