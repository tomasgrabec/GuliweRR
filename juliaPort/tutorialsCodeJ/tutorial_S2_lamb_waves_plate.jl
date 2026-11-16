## Tutorial S2 (code) -- free plate (Lamb waves): modes, dispersion, mode types, ZGV, coated plate, Si(001) ZGV check
# Narrated version: juliaPort/tutorialsJ/tutorial_S2_lamb_waves_plate.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_S2_lamb_waves_plate.jl")
# Units: thickness in mm, wavevector in rad/mm, velocity in km/s -> frequency in MHz.

using GuliwerrJ
using GLMakie
using Printf
using Statistics

## Material -- Cu(001) plate
dirMat      = [1 0 0; 0 0 1]       # plate normal; propagation reference
material, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)   # constants + density from the materials library

## Discretisation -- relDepth is rescaled from realThickness * k at every point
relDepthPlaceholder = 0.2          # overwritten; only the layer RATIO matters
legOrder   = 20                    # plate up to ~3.2 wavelengths thick; converged below 6 MHz
calcParams = CalcParams(; relDepth = relDepthPlaceholder, legOrder = legOrder)
nMaxModes  = 500
structure  = Structure(; materIn = material, calcParamsIn = calcParams, name = "cu_001_plate",
                         nMaxModes = nMaxModes)
realThickness = 1.0                # mm, total plate thickness

## Part 1 -- one angle, one wavevector: mode list
angleDeg       = 0.0
realWavVec     = 2.0               # rad/mm, k*d = 2
retentionLevel = :full             # :full | :none (no :saw_only for plates)
singleResult = singlePlate(; structure = structure, angle = angleDeg,
                             realThickness = realThickness, realWavVec = realWavVec)
singleResult = calcAll(singleResult, retentionLevel)
ritzSolver   = singleResult.ritzSolver

modeTypeNames = ("Asymm", "Symm", "TH-Asymm", "TH-Symm")   # modesType codes 1-4
nShow = 8
for iMode in 1:nShow
    @printf("%2d: v %.4f, vg %.4f km/s, f %.4f MHz, %-8s (confidence %.2f)\n", iMode,
            ritzSolver.vels[iMode], ritzSolver.velsGroup[iMode], ritzSolver.freqsPhysical[iMode],
            modeTypeNames[ritzSolver.modesType[iMode]], ritzSolver.modeTypeConfidence[iMode])
end

criterion       = :modeType        # :modeType | :confidence
velCutoffFactor = 2.5              # x longitudinal velocity (default 1.1 hides most plate modes)
axModes = plotRitzModeCurvePlate(singleResult, criterion, velCutoffFactor)
axModes.title = "Cu(001) plate, k = $realWavVec rad/mm: modes by type"
display(axModes.parent)            # plot functions and viewers only build figures; display shows them

## Part 2 -- one angle, many wavevectors: dispersion curves
wavVecs        = collect(0.1:0.1:20.0)   # rad/mm, dense for branch tracking, wide for whole branches
retentionSweep = :full             # REQUIRED for tracking, classifyBranches, locateZGVs
useParallel    = :auto             # :auto | :always | :never
sweepDisp = calc(wavevecAngSweep2Plate(; structure = structure, angles = [angleDeg],
                                         wavVecs = wavVecs, realThickness = realThickness,
                                         retentionLevel = retentionSweep, useParallel = useParallel))

angleIdx              = 1
colorMode             = :type      # :type | :branch
confidenceThresh      = 0.3        # default, unvalidated
showChristoffelBounds = true       # bulk-velocity reference lines
curveTypeVphase = :k_vphase        # :wavelength_freq | :wavelength_vphase | :k_freq | :k_vphase | :freq_vphase | :freq_vgroup
curveTypeFreq   = :k_freq
curveTypeGroup  = :freq_vgroup
velLimsDisp     = (0.0, 8.0)       # km/s
freqLimsDisp    = (0.0, 6.0)       # MHz
axVphase = plotDispersionCurvesPlate(sweepDisp, curveTypeVphase, angleIdx, colorMode, confidenceThresh, showChristoffelBounds)
ylims!(axVphase, velLimsDisp...)
axVphase.xlabel = "k [rad/mm]"     # units follow realThickness and wavVecs
axVphase.title = "Cu(001) plate: phase velocity"
display(axVphase.parent)
axFreq = plotDispersionCurvesPlate(sweepDisp, curveTypeFreq, angleIdx, colorMode, confidenceThresh)
ylims!(axFreq, freqLimsDisp...)
axFreq.xlabel = "k [rad/mm]"
axFreq.ylabel = "frequency [MHz]"
axFreq.title = "Cu(001) plate: frequency"
display(axFreq.parent)
axGroup = plotDispersionCurvesPlate(sweepDisp, curveTypeGroup, angleIdx, colorMode, confidenceThresh)
xlims!(axGroup, freqLimsDisp...)
axGroup.xlabel = "frequency [MHz]"
hlines!(axGroup, [0.0]; color = :black, linestyle = :dot)   # zero group velocity at finite k = ZGV
axGroup.title = "Cu(001) plate: group velocity"
display(axGroup.parent)

## Part 2b -- branch labels and zero-group-velocity points
classResult = classifyBranches(sweepDisp, angleIdx, confidenceThresh)
println("Mid-plane / sagittal mirror: ", classResult.isExactMirrorSymmetric, " / ",
        classResult.isExactSagittalSymmetric, "; first-k labels: ",
        join(classResult.displayLabel[1, 1:6], ", "))
maxCount = Inf                     # cap on returned ZGV points
zgvList  = locateZGVs(sweepDisp, angleIdx, maxCount)
zgvShow  = filter(z -> z.freqZGV <= freqLimsDisp[2], zgvList)   # the plotted frequency range
for zgv in zgvShow                 # fixed-angle search: isTrueZGV only on a symmetry axis (Part 5)
    @printf("ZGV branch %2d: k = %.4f rad/mm, f = %.4f MHz, true ZGV: %s\n", zgv.branchNo, zgv.wavVecZGV,
            zgv.freqZGV, zgv.isTrueZGV)
end

guiDisp = dispersionViewer(sweepDisp)
display(guiDisp.fig)
guiDisp.curveTypeMenu.i_selected[] = findfirst(==(curveTypeVphase), visDispersionCurveTypes())   # opens on :wavelength_freq
setAxisRange!(guiDisp, :y, velLimsDisp)

## Part 3 -- many angles, one wavevector
anglesDeg = 0.0:2.0:90.0
sweepAng  = calc(angSweepPlate(; structure = structure, angles = anglesDeg,
                                 realThickness = realThickness, realWavVec = realWavVec,
                                 retentionLevel = retentionLevel, useParallel = useParallel))
velLimsAng = (1.0, 5.5)            # km/s, both maps and the viewer
figMaps    = Figure(; size = (1400, 520))   # type and confidence side by side
axMapType  = Axis(figMaps[1, 1])
plotRitzModeMapPlate(sweepAng, criterion; ax = axMapType)   # .velsMap etc. are dense [nModes x nAngles]
ylims!(axMapType, velLimsAng...)
axMapType.title = "Cu(001) plate, k = $realWavVec rad/mm: by type"
axMapConf  = Axis(figMaps[1, 2])
plotRitzModeMapPlate(sweepAng, :confidence; ax = axMapConf)
ylims!(axMapConf, velLimsAng...)
axMapConf.title = "mode-type confidence"
display(figMaps)
guiAng = combinedViewer(sweepAng)
display(guiAng.fig); display(guiAng.modalShapeViewer.shapeFig)
setAxisRange!(guiAng, :y, velLimsAng)

## Part 4 -- coated plate (10% Al on one face): not mirror-symmetric, labels become heuristic
matAlCrystal, _ = MaterialFromLibrary("Al")   # isotropic Al = Hill average of the single crystal
alC11, alC44 = HillAverages(matAlCrystal.elConsts.voigt)
alConsts  = [alC11, alC44]         # [C11, C44] GPa
isoSymmNo = 1
alDensity = matAlCrystal.rho       # g/cm^3
ecAl      = ElastConsts(alConsts, isoSymmNo)
matAl     = Material(; elConsts = ecAl, rho = alDensity)
coatingFraction = 0.1              # relDepth ratio = thickness ratio
cpCu = CalcParams(; relDepth = (1 - coatingFraction) * relDepthPlaceholder, legOrder = legOrder)
cpAl = CalcParams(; relDepth = coatingFraction * relDepthPlaceholder, legOrder = legOrder)
structureCoated = Structure(; materIn = material, calcParamsIn = cpCu, name = "cu_plate_al_coating",
                              nMaxModes = nMaxModes)
structureCoated = AddLayer(structureCoated, matAl, cpAl)
println("Coated plate mirror-symmetric: ", isMirrorSymmetric(structureCoated))

sweepCoated = calc(wavevecAngSweep2Plate(; structure = structureCoated, angles = [angleDeg],
                                           wavVecs = wavVecs, realThickness = realThickness,
                                           retentionLevel = retentionSweep, useParallel = useParallel))
classCoated = classifyBranches(sweepCoated, angleIdx, confidenceThresh)
@printf("Trusted labels: %.0f%%\n", 100 * mean(classCoated.trustLabelByBranch))
colorModeBranch = :branch          # escape hatch when labels are not trusted
axCoatedType = plotDispersionCurvesPlate(sweepCoated, curveTypeVphase, angleIdx, colorMode, confidenceThresh)
ylims!(axCoatedType, velLimsDisp...)
axCoatedType.xlabel = "k [rad/mm]"
axCoatedType.title = "Cu + 10% Al: phase velocity, by type"
display(axCoatedType.parent)
axCoatedBranch = plotDispersionCurvesPlate(sweepCoated, curveTypeVphase, angleIdx, colorModeBranch, confidenceThresh)
ylims!(axCoatedBranch, velLimsDisp...)
axCoatedBranch.xlabel = "k [rad/mm]"
axCoatedBranch.title = "Cu + 10% Al: phase velocity, by branch"
display(axCoatedBranch.parent)

## Part 5 -- Si(001) wafer: S1 ZGV vs published measurements (Kiefer, Mezil & Prada 2023, Sci. Adv. 9, eadk6846)
siConsts  = [165.6, 63.9, 79.5]    # [C11, C12, C44] GPa, Kiefer et al. (library Si differs ~0.1 %)
siDensity = 2.330                  # g/cm^3
siDirMat  = [0 0 1; 1 0 0]         # normal [001]; angle 0 = [100], 45 = [110]
cubicSymmNo = 2
materialSi  = Material(; elConsts = ElastConsts(siConsts, cubicSymmNo), rho = siDensity, dirMat = siDirMat)
structureSi = Structure(; materIn = materialSi, calcParamsIn = calcParams, name = "si_001_wafer",
                          nMaxModes = nMaxModes)
siThickness = 0.5246               # mm
siAnglesDeg = [0.0, 45.0, 18.4]    # [100], [110], 26.6 deg from [110] (a TGV wave)
siWavVecs   = collect(range(1.0, 6.0; length = 51))   # rad/mm
sweepSi = calc(wavevecAngSweep2Plate(; structure = structureSi, angles = siAnglesDeg, wavVecs = siWavVecs,
                                       realThickness = siThickness,
                                       retentionLevel = retentionSweep, useParallel = useParallel))
publishedComputed = [7.7079, 7.7551, 7.727]   # MHz (GEWtool)
publishedMeasured = [7.7072, 7.7545, NaN]     # MHz
for iAng in eachindex(siAnglesDeg)
    zgvSi = locateZGVs(sweepSi, iAng, Inf)
    zLow  = zgvSi[argmin([z.freqZGV for z in zgvSi])]   # S1 point (A2 is near 15 MHz)
    @printf("%4.1f deg: k = %.4f rad/mm, f = %.4f MHz (published %.4f, measured %.4f), true ZGV %s, c_gtheta = %.4f mm/us\n",
            siAnglesDeg[iAng], zLow.wavVecZGV, zLow.freqZGV, publishedComputed[iAng], publishedMeasured[iAng],
            zLow.isTrueZGV, zLow.velGroupTrans)   # 18.4 deg: TGV, 0.178 (published -0.177 at the mirror direction, 71.6 deg)
end
