## Tutorial: thin film on substrate -- SAW velocity vs. kd (wavevector-thickness product) (Julia port)
#
# Same physics as examples/tutorials/tutorial_06_film_on_substrate.m: a film-on-substrate SAW
# velocity depends on kd (here, film relDepth = realThickness/wavelength): near relDepth = 0,
# v_SAW sits at the substrate's own value; for relDepth >> 1 it approaches the film's own
# half-space value.
#
# Material pair -- Al substrate / Cu film -- is a SLOW film on a FAST substrate (Cu's bulk
# velocities are lower than Al's despite higher elastic moduli, because it's 3.3x denser). For
# why a per-point pick can jump between branches in a layered sweep, and how the branch
# tracking below avoids it, see docs/pitfalls.md MTH11.
#
# The relDepth sweep runs at one fixed angle through wavevecAngSweep2Layered with a one-element
# angle vector rather than a hand-written loop: the sweep object holds every wavelength step,
# which DispersionBranches.trackSeededBranch needs to follow one physical branch by
# Ritz-eigenvector overlap (MAC).
#
# Copy this script as a starting point. See docs/pitfalls.md for the legO/relDepth rule.
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - `vcat` for building one flat numeric vector out of scalars and ranges mixed together --
#     genuinely different from MATLAB here, explained inline below at the point it matters

using GuliwerrJ
using GLMakie   # opens an interactive window; `using CairoMakie` instead renders the same
                 # Figure/Axis objects headlessly (e.g. for save("fig.png", fig)) -- see
                 # tutorial_02_angular_map.jl for the full explanation

## Substrate -- isotropic Al (faster than the Cu film below, both in
# v_L = sqrt(C11/rho) and v_T = sqrt(C44/rho), because Al's density is
# 3.3x lower than Cu's despite Cu's higher absolute elastic moduli). Isotropic
# (polycrystalline) Al = the Hill average of the single-crystal constants from the materials
# library (tutorial_S1 section 9 shows the averages).
matAlCrystal, _ = MaterialFromLibrary("Al")
alC11, alC44 = HillAverages(matAlCrystal.elConsts.voigt)
alConsts = [alC11, alC44]   # [C11, C44] GPa, isotropic class
alSymm   = 1
alRho    = matAlCrystal.rho   # g/cm^3
relDepthSub = 20.0   # substrate truncation depth, in wavelengths

## Film -- Cu(001), cubic, slower than the Al substrate
# Cu constants and density come from the materials library.
dirMatFilm = [1 0 0;    # surface normal
              0 0 1]     # propagation reference (Material's default)
filmThicknessNm = 200.0   # real film thickness [nm] -- ApplyWavelength (called internally, once
                           # per wavelength) divides this by each swept wavelength to derive
                           # relDepth

## Discretisation -- legOrder must stay >= 2x the LARGEST relDepth reached
# anywhere in the sweep (film relDepth sweeps up to 3, so legOrder=50 leaves a
# wide margin). Cu's Zener ratio 2*C44/(C11-C12) = 3.2 puts it in
# pitfalls.md MTH1's "easy" tier at this principal-plane angle.
legOrder = 50   # number of Legendre basis functions spanning the depth

ecSub  = ElastConsts(alConsts, alSymm)
matSub = Material(; elConsts = ecSub, rho = alRho)
cpSub  = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)

matFilm, _ = MaterialFromLibrary("Cu"; dirMat = dirMatFilm)
relDepthFilmPlaceholder = 1.0   # overwritten by ApplyWavelength once per wavelength below
cpFilm = CalcParams(; relDepth = relDepthFilmPlaceholder, legOrder = legOrder,
                       realThickness = filmThicknessNm)

strName   = "al_sub_cu_film"
nMaxModes = 500
str = Structure(; materIn = matSub, calcParamsIn = cpSub, name = strName, nMaxModes = nMaxModes)
str = AddLayer(str, matFilm, cpFilm)   # layer 2 = film, on top

## Sweep -- ONE fixed angle, many real wavelengths (relDepth = filmThicknessNm/wavelength).
# Descending wavelength = ascending relDepth, so tracking starts from the thinnest film, where
# the branches are clearest. Tracking matches CONSECUTIVE steps by eigenvector overlap, so the
# steps are dense (0.01) up to relDepth = 1 and coarser (0.05) beyond. The sweep starts at
# relDepth = 0.01 (a 200 nm film at a 20 um wavelength), where the branch is still close to the
# Al value but the film is no longer negligible.
filmRelDepths = vcat(0.01:0.01:1.0, 1.05:0.05:3.0)
# GENUINELY DIFFERENT FROM MATLAB: MATLAB's own `[0.01:0.01:1, 1.05:0.05:3]` concatenates
# several colon-expressions into ONE flat row vector automatically, because MATLAB's
# square-bracket literal always performs concatenation. Julia's square-bracket literal does NOT
# -- `[0.01:0.01:1.0, 1.05:0.05:3.0]` would build a 2-ELEMENT `Vector` holding the two Range
# objects side by side, not their flattened union. `vcat(...)` ("vertical concatenate", used
# here on 1-D vectors/ranges so the result is one flat vector) is the function that actually
# does what MATLAB's brackets did implicitly -- reach for it (or `collect` a range first)
# whenever you want to glue several ranges/scalars into one plain numeric vector.
wavelengths = filmThicknessNm ./ filmRelDepths   # nm, descending
# `./` is elementwise division: a scalar divided by a vector needs the dot, `filmThicknessNm /
# filmRelDepths` alone would be a `MethodError` (matrix "division" is not what this means, and
# Julia will not guess).

angle          = 0.0     # deg -- fixed propagation angle for this sweep
calcMode       = :ritzOnly   # vMax is not needed for :ritzOnly
retentionLevel = :full   # REQUIRED: branch tracking compares every eigenvector of consecutive
                          # wavelength steps
sweep = wavevecAngSweep2Layered(; structure = str, angles = [angle], wavelengths = wavelengths,
                                    calcMode = calcMode, retentionLevel = retentionLevel)
sweep = calc(sweep)

## Branch tracking -- one physical branch, followed by Ritz-eigenvector overlap (MAC). The
# seed is the slowest classified SAW at the thinnest film: the Rayleigh-type wave.
# (`DispersionBranches.track(sweep)` gives every branch.)
vSAW, _, _ = DispersionBranches.trackSeededBranch(sweep)
# `trackSeededBranch` returns THREE values (velocity of the branch at every step, its raw Ritz
# mode index at every step, and the seed's mode number); `_` is Julia's conventional name for
# a value you do not need.

println("film relDepth : ", join(round.(filmRelDepths, digits = 2), " "))
println("v_SAW [km/s]  : ", join(round.(vSAW, digits = 3), " "))
# `round.(x, digits = 3)` rounds every element (the dot broadcasts `round` over the vector);
# `join(..., " ")` glues the pieces into one space-separated string -- Julia's counterpart of
# MATLAB's `fprintf(' %6.3f', v)` trick for printing a whole vector on one line.

## Reference lines -- each material's OWN half-space SAW velocity (angle=0),
# the two limits the tracked branch moves between. Reuses matSub/matFilm
# built above, standing alone as a half-space instead of layered.
strAl = Structure(; materIn = matSub, calcParamsIn = CalcParams(; relDepth = relDepthSub,
                     legOrder = legOrder), name = "al_hs", nMaxModes = nMaxModes)
csAl  = singleHS(; structure = strAl, angle = angle)
csAl  = calcAll(csAl, :ritzOnly)
smAl    = csAl.classification.significantModes
vSAW_Al = minimum(m.vel for m in smAl if m.type === :Rayleigh)   # the Rayleigh wave

strCu = Structure(; materIn = matFilm, calcParamsIn = CalcParams(; relDepth = relDepthSub,
                     legOrder = legOrder), name = "cu_hs", nMaxModes = nMaxModes)
csCu  = singleHS(; structure = strCu, angle = angle)
csCu  = calcAll(csCu, :ritzOnly)
smCu    = csCu.classification.significantModes
vSAW_Cu = minimum(m.vel for m in smCu if m.type === :Rayleigh)

using Printf
@printf("Al half-space v_SAW: %.4f km/s\n", vSAW_Al)
@printf("Cu half-space v_SAW: %.4f km/s\n", vSAW_Cu)
# `significantModes` is the classification's list of accepted entries (one per classified wave:
# `type`, `vel`, `modeNo`, `method`, `score`, `scoreDb`, `grade`, ...); `type === :Rayleigh`
# picks the Rayleigh wave. The generator `(m.vel for m in smAl if ...)` is Julia's counterpart
# of MATLAB's `[smAl(strcmp({smAl.type}, 'Rayleigh')).vel]`.

## Plot 1 -- the tracked branch vs. relDepth, with the two half-space limits marked.
fig1 = Figure()
ax1 = Axis(fig1[1, 1]; xlabel = "film relDepth (∝ kd)", ylabel = "v_SAW [km/s]",
           xscale = log10, title = "Al substrate / Cu film -- tracked branch vs. film relDepth")
scatterlines!(ax1, filmRelDepths, vSAW; marker = :circle)
hlines!(ax1, [vSAW_Al]; color = (:green, 0.8), linestyle = :dash,
        label = "Al (substrate), $(round(vSAW_Al, digits = 3)) km/s")
hlines!(ax1, [vSAW_Cu]; color = (:orange, 0.8), linestyle = :dash,
        label = "Cu (film), $(round(vSAW_Cu, digits = 3)) km/s")
axislegend(ax1; position = :rb)
display(fig1)
# `scatterlines!`/`hlines!`/`axislegend` are Makie's own idiom for MATLAB's
# `semilogx(...,'o-')` + `yline(...)` + implicit legend -- a log-scale x-axis is set via the
# Axis keyword `xscale = log10` (a function, not a string), rather than a separate call.

## Plot 2 -- plotDispersionCurves: every tracked Ritz branch (thin grey), with SAW/limiting
# branches highlighted -- reproduces Plot 1's tracked branch plus the rest of the raw Ritz
# spectrum for context.
ax2 = plotDispersionCurves(sweep, :k_vphase)
ax2.title = "Al substrate / Cu film -- all tracked branches (plotDispersionCurves)"
xlims!(ax2, 0.0, 0.03)   # k [rad/nm]
ylims!(ax2, 1.0, 4.0)    # km/s
display(ax2.parent)

## Optional: interactive exploration -- click a curve (or step Wavelength/Branch) to see that
# mode's own displacement field. Genuinely interactive only with GLMakie's own window and mouse
# events -- run this part from a real terminal/REPL session with GLMakie loaded, exactly as
# MATLAB's own comment here says "Interactive MATLAB session only (not guarded-matlab.sh's
# -batch mode, no display)". The VIEWER OBJECT ITSELF still builds and can be driven
# programmatically under CairoMakie too (this tutorial's own automated verification does
# exactly that), it just won't respond to mouse clicks without a real window.
#
# `dispersionViewer`, like every plotting function used so far, is a plain top-level GuliwerrJ
# function -- it exists only once a Makie backend is loaded (its real body lives in the same
# package EXTENSION mechanism `using GLMakie` already triggered for the plots above), but
# needs no extra lookup to reach: `using GuliwerrJ, GLMakie` is enough for `dispersionViewer`
# to resolve directly, exactly like `plotDispersionCurves` did just above.
gui = dispersionViewer(sweep)
display(gui.fig)

# Default the map to k-vphase with a fixed axis range rather than the constructor's
# wavelength-freq/autoscale default. Unlike MATLAB's `gui.curveTypeDropdown.Value = 'k-vphase';
# gui.onCurveTypeChanged();` two-step, setting the underlying Observable (`i_selected`) below
# triggers the redraw AUTOMATICALLY -- there is no separate "changed" callback to invoke by
# hand, because the whole viewer is wired around Observables reacting to their own updates.
kVphaseIdx = findfirst(==(:k_vphase), visDispersionCurveTypes())
gui.curveTypeMenu.i_selected[] = kVphaseIdx
setAxisRange!(gui, :x, (0.0, 0.03))   # also fills the viewer's range fields
setAxisRange!(gui, :y, (1.0, 4.0))

## Optional: the same sweep in :full mode -- adds a layered Green map (|G13_eff|) at every
# wavelength, the response a measurement would see, behind the tracked Ritz branches. A few
# seconds extra here.
calcModeFull = :full   # Ritz + substrate Christoffel + Green
sweepFull = wavevecAngSweep2Layered(; structure = str, angles = [angle], wavelengths = wavelengths,
                                      calcMode = calcModeFull, retentionLevel = retentionLevel)
sweepFull = calc(sweepFull)

angleIdx            = 1        # the only angle
showRitzBranches    = true
showGreenBackground = true     # needs a v_phase curveType (:k_vphase | :wavelength_vphase)
scaleMode           = :log
climRange           = (-3.0, -0.5)   # log10(|G|/max|G|) colour limits
axFull = plotDispersionCurves(sweepFull, :k_vphase, angleIdx, showRitzBranches,
                              showGreenBackground, scaleMode, climRange)
axFull.title = "Al substrate / Cu film -- tracked branches over |G_13^eff|"
xlims!(axFull, 0.0, 0.03)   # k [rad/nm]
ylims!(axFull, 1.0, 4.0)    # km/s
display(axFull.parent)

# Interactive, Green background on. combinedViewer takes angular and single-angle results only;
# for a wavelength sweep the dispersion viewer is its counterpart (curves linked to the modal
# shape). Switching the toggle's Observable redraws the map, as `i_selected` did above.
guiFull = dispersionViewer(sweepFull)
display(guiFull.fig)
guiFull.curveTypeMenu.i_selected[] = kVphaseIdx
guiFull.showGreenBgToggle.active[] = true
setAxisRange!(guiFull, :x, (0.0, 0.03))
setAxisRange!(guiFull, :y, (1.0, 4.0))
