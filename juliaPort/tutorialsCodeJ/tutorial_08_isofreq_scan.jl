## Tutorial 08 (code) -- isoFreq scan: Ni film on Cu(001) at fixed real frequency, |G| vs k_parallel
# Narrated version: juliaPort/tutorialsJ/tutorial_08_isofreq_scan.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_08_isofreq_scan.jl")

using GuliwerrJ
using GLMakie

## Substrate -- Cu(001); treated analytically, its CalcParams are never read
dirMatSub   = [1 0 0; 0 0 1]
matSub, _   = MaterialFromLibrary("Cu"; dirMat = dirMatSub)   # constants + density from the materials library
cpSub       = CalcParams()         # defaults, kept only so the Structure is well-formed

## Film -- Ni(001), 100 nm
dirMatFilm          = [1 0 0; 0 0 1]
filmThicknessNm     = 100.0
matFilm, _          = MaterialFromLibrary("Ni"; dirMat = dirMatFilm)   # materials library
relDepthPlaceholder = 1.0          # overwritten by ApplyWavelength below
legOrder            = 50
cpFilm = CalcParams(; relDepth = relDepthPlaceholder, legOrder = legOrder,
                      realThickness = filmThicknessNm)

## Structure, normalised to a reference wavelength
structureName = "cu_sub_ni_film"
nMaxModes     = 500
lambdaRefNm   = 300.0              # free normalisation, ~ the expected SAW wavelength
structure = Structure(; materIn = matSub, calcParamsIn = cpSub, name = structureName,
                        nMaxModes = nMaxModes)
structure = AddLayer(structure, matFilm, cpFilm)
structure = ApplyWavelength(structure, lambdaRefNm)

## Fixed frequency -> normalised omega0, and a k_parallel axis from a velocity axis
freqGHz   = 5.0
omega0    = freqToOmega0(freqGHz * 1e9, lambdaRefNm)   # lambdaRefNm must match ApplyWavelength
velMax    = 6.0                    # km/s, brackets the expected bulk/SAW velocities
velMin    = 1.0                    # km/s
nVel      = 300
velAxis   = range(velMax, velMin; length = nVel)       # descending -> kParRange ascending
kParRange = omega0 ./ velAxis

## Step 1 -- single angle: 1D |G_33| vs k_parallel
angleDeg   = 0.0
singleFreq = isofreqSingleLayered(; structure = structure, angle = angleDeg, omega0 = omega0,
                                    kParRange = kParRange)
singleFreq = calc(singleFreq)
axSingle = plotDispersionMap(singleFreq, 1, :log, (-2.5, 0.0))   # normalised log: peak at 0
axSingle.title = "Cu(001)/$(round(Int, filmThicknessNm)) nm Ni, 0 deg, $freqGHz GHz: isoFreq spectrum"
display(axSingle.parent)           # plot functions and viewers only build figures; display shows them

## Step 2 -- angle sweep at the same frequency: 2D map
anglesDeg   = 0.0:2.0:90.0
useParallel = :auto                # :auto | :always | :never
sweep = isofreqAngSweepLayered(; structure = structure, angles = anglesDeg, omega0 = omega0,
                                 kParRange = kParRange, useParallel = useParallel)
sweep = calc(sweep)
axSweep = plotDispersionMap(sweep, 1, :log, (-3.5, -1.5))
tightlimits!(axSweep)
axSweep.title = "Cu(001)/$(round(Int, filmThicknessNm)) nm Ni, $freqGHz GHz: isoFreq |G_33| map"
display(axSweep.parent)

## Step 3 -- several frequencies: one shared kParRange must span all of them
freqListGHz    = [3.0, 5.0, 7.0]
omega0List     = freqToOmega0.(freqListGHz .* 1e9, lambdaRefNm)
nK             = 400
kParRangeMulti = collect(range(minimum(omega0List) / velMax, maximum(omega0List) / velMin; length = nK))
sweepMulti = freqAngSweep2Layered(; structure = structure, angles = anglesDeg,
                                    omega0List = omega0List, kParRange = kParRangeMulti,
                                    useParallel = useParallel)
sweepMulti = calc(sweepMulti)
for omegaIdx in eachindex(omega0List)
    axSlice = plotDispersionMap(sweepMulti, omegaIdx, :log, (-3.5, -1.5))   # one frequency slice per figure
    tightlimits!(axSlice)
    axSlice.title = "Cu(001)/$(round(Int, filmThicknessNm)) nm Ni, $(freqListGHz[omegaIdx]) GHz: isoFreq |G_33| map"
    display(axSlice.parent)
end

## Interactive frequency-slice viewer
gui = doubleSweepViewer(sweepMulti)
setAxisRange!(gui, :c, (-3.5, -1.5))
display(gui.fig)
