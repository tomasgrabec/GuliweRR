## Tutorial: isoFreq scan -- fixed real frequency, Green response vs k_parallel (Julia port)
#
# Same physics as examples/tutorials/tutorial_08_isofreq_scan.m: the isoFreq family answers a
# different question from every other tutorial so far: "at one fixed real frequency, which
# in-plane wavevectors k_parallel carry a guided mode?" (every other tutorial is isoWavevec: at
# fixed |k_parallel| = 1, which velocities are supported). The two coincide only for a
# non-dispersive medium (docs/background.md).
#
# This is the Green's-function-only path (no Ritz -- Ritz fixes |k_parallel| = 1 by
# construction, so it has no free parameter left to hold frequency fixed while k_parallel
# varies): calcMode/retentionLevel don't apply here.
#
# Copy this script as a starting point for a new isoFreq calculation.
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - broadcasting with the dot (`.`), Julia's general elementwise mechanism: `omega0 ./ velAxis`
#     is required here (no MATLAB-style implicit elementwise `/`); explained inline at the one
#     point it is genuinely new versus the plain scalar arithmetic used so far
#   - descending ranges (`range(6.0, 1.0; length = ...)`), which Julia handles exactly like an
#     ascending one -- nothing extra needed, noted only because it looks unusual at first glance

using GuliwerrJ
using GLMakie   # opens an interactive window; `using CairoMakie` instead renders the same
                 # Figure/Axis objects headlessly (e.g. for save("fig.png", fig)) -- see
                 # tutorial_02_angular_map.jl for the full explanation

## Substrate -- Cu(001), cubic. The Green solver treats it as a semi-infinite half-space, so
# its CalcParams are never read; the defaults just keep the Structure valid (and usable for a
# Ritz calculation later).
# The Cu constants and density come from the materials library (`ListMaterialLibrary(; print = true)`).
dirMatSub = [1 0 0;    # surface normal
             0 0 1]     # propagation reference (Material's default)

matSub, _ = MaterialFromLibrary("Cu"; dirMat = dirMatSub)
cpSub  = CalcParams()   # every keyword left at its default (relDepth = 20.0, legOrder = 50)

## Film -- Ni(001), cubic, thin, on top (constants and density from the materials library)
dirMatFilm = [1 0 0;    # surface normal
              0 0 1]     # propagation reference (Material's default)

filmThicknessNm = 100.0   # real film thickness [nm]
lambdaRefNm     = 300.0   # reference wavelength [nm] -- a free normalisation choice (any
                            # positive value gives the same physical answer), pick it on the
                            # order of the expected SAW wavelength so the derived relDepth
                            # stays well-conditioned (freqToOmega0's own derivation note)

matFilm, _ = MaterialFromLibrary("Ni"; dirMat = dirMatFilm)

relDepthFilmPlaceholder = 1.0   # overwritten by ApplyWavelength below
legOrder = 50   # not read by the Green solver; a valid value for later Ritz use
cpFilm = CalcParams(; relDepth = relDepthFilmPlaceholder, legOrder = legOrder,
                       realThickness = filmThicknessNm)

strName   = "cu_sub_ni_film"
nMaxModes = 500
str = Structure(; materIn = matSub, calcParamsIn = cpSub, name = strName, nMaxModes = nMaxModes)
str = AddLayer(str, matFilm, cpFilm)   # layer 2 = film, on top
str = ApplyWavelength(str, lambdaRefNm)   # film relDepth = thickness/lambdaRefNm

## Fixed real frequency -> GreenSolver's normalised omega0
freqGHz = 5.0
omega0  = freqToOmega0(freqGHz * 1e9, lambdaRefNm)

## k_parallel axis, built from a velocity axis for physical intuition.
# velAxis DESCENDING -> kParRange = omega0./velAxis comes out ASCENDING --
# range chosen to bracket the expected SAW/bulk velocities of this Cu/Ni pair.
nVel    = 300
velAxis = range(6.0, 1.0; length = nVel)   # km/s -- a descending Range; Julia has no trouble
                                             # with `step < 0` here, no special handling needed
kParRange = omega0 ./ velAxis
# The `.` in `omega0 ./ velAxis` is NOT optional here: `omega0` is one number, `velAxis` is a
# whole Range of 300 -- Julia's plain `/` only ever means "one honest division" (or a linear
# solve, for matrices), and refuses to guess that you meant "divide elementwise" the way
# MATLAB's `/` sometimes does. Putting a dot on ANY operator or function call means "apply this
# elementwise" (Julia calls this "broadcasting") -- `omega0 ./ velAxis` reads as "omega0, divided elementwise into velAxis".

## Step 1: single angle -- 1D spectrum, |G_33| vs k_parallel
angleSingle = 0.0   # deg
singleFreq = isofreqSingleLayered(; structure = str, angle = angleSingle, omega0 = omega0,
                                     kParRange = kParRange)
singleFreq = calc(singleFreq)

# 4th argument = y-range of the log10(|G|/max|G|) line; the peak sits at 0 by construction
# (default (-5, 0))
axSingle = plotDispersionMap(singleFreq, 1, :log, (-2.5, 0.0))
axSingle.title = "Cu(001)/$(round(Int, filmThicknessNm)) nm Ni film, 0 deg, f = $freqGHz GHz -- isoFreq spectrum"
display(axSingle.parent)   # a plot function only builds the figure; `display` shows it (tutorial 2)

## Step 2: angle sweep at the same fixed frequency -- 2D map, angle x k_parallel
angles = 0.0:2.0:90.0
sweep  = isofreqAngSweepLayered(; structure = str, angles = angles, omega0 = omega0,
                                   kParRange = kParRange)
sweep  = calc(sweep)

axSweep = plotDispersionMap(sweep, 1, :log, (-3.5, -1.5))   # 4th argument = colour range
tightlimits!(axSweep)   # crop the empty margin around the map
axSweep.title = "Cu(001)/$(round(Int, filmThicknessNm)) nm Ni film, f = $freqGHz GHz -- isoFreq |G_33| map"
display(axSweep.parent)

## Step 3: several frequencies at once -- freqAngSweep2Layered loops Step 2 over a whole
# omega0List (a 3rd axis). ONE kParRange is shared by every frequency, so it must span all of
# them: k = omega0/v for the fastest velocity at the lowest frequency up to the slowest at the
# highest.
freqListGHz    = [3.0, 5.0, 7.0]
omega0List     = freqToOmega0.(freqListGHz .* 1e9, lambdaRefNm)   # `.` broadcasts over the list
velMax         = 6.0   # km/s
velMin         = 1.0   # km/s
kParRangeMulti = collect(range(minimum(omega0List) / velMax, maximum(omega0List) / velMin;
                               length = 400))
useParallel    = :auto   # :auto | :always | :never -- forwarded to every angle loop
sweepMulti = freqAngSweep2Layered(; structure = str, angles = angles, omega0List = omega0List,
                                     kParRange = kParRangeMulti, useParallel = useParallel)
sweepMulti = calc(sweepMulti)

for omegaIdx in eachindex(omega0List)   # plotDispersionMap shows ONE frequency slice at a time
    axSlice = plotDispersionMap(sweepMulti, omegaIdx, :log, (-3.5, -1.5))
    tightlimits!(axSlice)
    axSlice.title = "Cu(001)/$(round(Int, filmThicknessNm)) nm Ni film, f = $(freqListGHz[omegaIdx]) GHz -- isoFreq |G_33| map"
    display(axSlice.parent)
end

## Interactive: step through the frequency slices with Green's-response fields, no
# recomputation. Needs GLMakie (a window). `doubleSweepViewer` is this port's counterpart of
# MATLAB's Gui.freqAngSweep2LayeredViewer.
gui = doubleSweepViewer(sweepMulti)
setAxisRange!(gui, :c, (-3.5, -1.5))   # colour limits, kept across frequency steps
display(gui.fig)

# Not shown: GreenPeaks.find extracts peak positions from (absG, axis, mask); calc() does not
# pick peaks itself.
