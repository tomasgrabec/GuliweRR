## Tutorial: modal-shape displacement fields -- Ritz mode + Green's response (Julia port)
#
# Same physics as examples/tutorials/tutorial_03_modal_shape.m: shows the SAW mode's
# displacement field itself (depth x lateral), two ways:
#   - plotModalShape      -- the Ritz eigenmode's own displacement field (needs
#                             retentionLevel = :full)
#   - plotGreenModalShape -- the Green's-function response at that same SAW velocity
#                             (half-space only)
# plus a |G| vs. velocity spectrum from a direct `calc(GreenSolver, ...)` call (own velocity
# axis and component), and the interactive combinedViewer on the same result.
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - dispatching on a TYPE itself, not a value of that type: `calc(GreenSolver, ...)`. This is
#     Julia's way of writing what MATLAB spelled as a "static"/class method
#     (`GreenSolver.calc(...)`) -- multiple dispatch works on the TYPES of all arguments, and
#     passing the type `GreenSolver` itself as a value (its own type is `Type{GreenSolver}`) is
#     dispatch like any other, just naming which type's behaviour you want
#   - `error(...)`, Julia's plain equivalent of MATLAB's `error(id, msg)` (Julia has no error
#     identifiers, only a message)

using GuliwerrJ
using GLMakie   # see tutorial_02_angular_map.jl for the GLMakie/CairoMakie choice; CairoMakie
                 # is the drop-in headless alternative (`using CairoMakie` instead)

## Material -- Cu(001), cubic. dirMat set explicitly below even though it's
# Material's own default value.
dirMat   = [1 0 0;    # surface normal
            0 0 1]     # propagation reference (Material's default)

# The Cu constants and density come from the materials library (`ListMaterialLibrary(; print = true)`
# shows values and sources). `MaterialFromLibrary` returns the Material and the library row; the
# row is not needed here, so it is discarded with `_`.
mat, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)

## Discretisation -- "easy" tier (Cu Zener ratio 2*C44/(C11-C12) = 3.2, principal (001) plane)
relDepth = 20   # truncation depth of the computational domain, in wavelengths
legOrder = 50   # number of Legendre-polynomial basis functions spanning that
                # depth -- must be >= 2*relDepth or the depth profile is
                # under-resolved (docs/pitfalls.md MTH1)
cp = CalcParams(; relDepth = relDepth, legOrder = legOrder)

strName   = "cu_001_20deg"
nMaxModes = 500
str = Structure(; materIn = mat, calcParamsIn = cp, name = strName, nMaxModes = nMaxModes)

## Single-angle full sweep
angle          = 20.0
calcMode       = :full     # Ritz+Christoffel+Green together
retentionLevel = :full     # needed by plotModalShape (keeps full eigenvector matrix)
nVelGreen      = nothing   # Green-axis points (nothing = greenVelStep steps, 6 m/s)
sweep = angSweepHS(; structure = str, angles = [angle], calcMode = calcMode,
                      retentionLevel = retentionLevel, nVelGreen = nVelGreen)
sweep = calc(sweep)

if isempty(sweep.sawModeNoMap[1])
    error("tutorial_03_modal_shape: no SAW mode identified at $angle deg for this " *
          "material/angle.")
end
# Julia's plain `error(msg)` is the direct analogue of MATLAB's `error(id, msg)` -- Julia has
# no error-IDENTIFIER system (only `MATLAB:someID`-style strings do), so the message alone is
# what gets raised, as an `ErrorException`. `*` between two strings is Julia's concatenation
# operator (not `+`, and not `[a b]` the way MATLAB sometimes concatenates char arrays);
# string interpolation (`"$angle deg"`, splicing a variable straight into a string) is usually
# the more idiomatic way to build a message like this one, used here alongside `*` to keep the
# two-line layout readable.
modeNo = sweep.sawModeNoMap[1][1]
vSAW   = sweep.sawVelMap[1][1]

## Ritz modal shape -- :none substrateCutoff (only option for a half-space)
ritzSolver      = sweep.cs[1].ritzSolver
substrateCutoff = :none   # :none (true-scale depth) | :decay | :fixedMultiple | :boosted --
                          # the last three matter for film-on-substrate stacks
cutoffParam     = nothing # decay fraction (:decay, default 0.05) or depth multiple
                          # (:fixedMultiple, default 4); ignored for :none/:boosted
dirs            = 1:3     # draw all three displacement directions -- a Julia Range again
meshSpacing     = 0.01    # depth/lateral mesh spacing passed to calcDispField2D
axesModal = plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing)
figModal = axesModal[1].parent   # every Axis in the returned Vector shares one owning Figure
sawType = sweep.sawTypeMap[1][1]
Label(figModal[0, :],
      "Cu (001), θ = $angle deg: Ritz $(sawType) mode, v = $(round(vSAW, digits = 3)) km/s";
      fontsize = 16)
display(figModal)   # a plot function only builds the figure; `display` shows it (tutorial 2)
# Makie has no single `sgtitle(...)` function the way MATLAB does -- the idiomatic way to add
# one caption spanning every panel of a multi-Axis Figure is to place a `Label` of its own in
# an extra row (`figModal[0, :]`, row 0 = above every existing row 1..N) of the SAME
# `GridLayout` the panels already live in.

## The same mode as ONE polarisation-coloured panel (dirs = :rgb)
# Colour = direction of the displacement (red u_x = propagation, blue u_y = in-plane
# transverse, green u_z = out-of-plane; key on the right), brightness = |u| (white where
# nothing moves). It is one snapshot: u_x and u_z of a Rayleigh wave are 90 deg out of phase,
# so they show as red and green stripes a quarter wavelength apart. :rgb also works for
# plotGreenModalShape / plotGreenModalShapeLayered; the Julia viewers do not offer it as a
# Display dropdown yet (MATLAB's do) -- call the plot function instead.
axesRgb = plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, :rgb, meshSpacing)
figRgb = axesRgb[1].parent
Label(figRgb[0, :],
      "Cu (001), θ = $angle deg: polarisation view, v = $(round(vSAW, digits = 3)) km/s";
      fontsize = 16)
display(figRgb)

## Green's-response field at that same SAW velocity
greenObj = calc(GreenSolver, str, angle, [vSAW]; fullOutput = true)
# `calc(GreenSolver, ...)` -- see this file's header. `velRange` must be a VECTOR even for one
# velocity (`[vSAW]`, not bare `vSAW`) since GreenSolver's own signature always sweeps a range;
# `fullOutput = true` is REQUIRED here (unlike MATLAB's own `nargout > 1` auto-detection,
# Julia has no way to tell how many outputs a caller "wants" -- there is only ever one return
# value, or one tuple of them -- so the extra partial-wave diagnostics plotGreenModalShape
# needs must be requested explicitly by name).
forceDir = 3   # 1 = propagation, 2 = in-plane transverse, 3 = surface normal
axesGreen = plotGreenModalShape(greenObj, forceDir)
figGreen = axesGreen[1].parent
Label(figGreen[0, :],
      "Cu (001), θ = $angle deg: Green response, v = $(round(vSAW, digits = 3)) km/s";
      fontsize = 16)
display(figGreen)

## Green's-function spectrum |G| vs. velocity at this one angle -- called directly (not
# through the sweep) to choose the velocity axis and the (force, displacement) component.
# Peaks mark surface-bound waves: a sharp pole = SAW, a sharp finite peak = pseudo-SAW.
longLimVel   = sweep.cs[1].longLimVel          # longitudinal limiting velocity, km/s
velMinGreen  = 1.5                             # km/s, lower end of the scan
velAxisGreen = collect(range(velMinGreen, 1.05 * longLimVel; length = 1000))
forceDirG    = 3   # force direction I: 1 | 2 | 3 (x1 = propagation, x3 = surface normal)
displDirG    = 3   # displacement direction J: 1 | 2 | 3 -- (1, 3) gives G_13 (TGS)
greenSpec    = calc(GreenSolver, str, angle, velAxisGreen; forceDir = forceDirG,
                    displDir = displDirG)
# forceDir/displDir are KEYWORD arguments here (after the `;`), unlike MATLAB's positional
# 4th/5th arguments -- they have defaults (1 and 3), so a keyword makes each one optional alone.
absG = abs.(greenSpec.G)   # `abs.` -- the dot broadcasts abs over every element

figSpec = Figure()
axSpec  = Axis(figSpec[1, 1]; xlabel = "velocity [km/s]", ylabel = "|G_$(forceDirG)$(displDirG)| [a.u.]",
               yscale = log10, title = "Cu (001), θ = $angle deg: Green's function vs. velocity")
lines!(axSpec, velAxisGreen, absG; linewidth = 1.5)
vlines!(axSpec, [vSAW]; linestyle = :dash, color = :black)   # the Ritz SAW velocity
display(figSpec)

## Interactive: the same result in combinedViewer -- map or curve (bottom) linked to the
# modal shape (top); click a point to see its shape. Needs GLMakie (a window).
gui = combinedViewer(sweep)
display(gui.fig)                         # the map/curve and the controls
display(gui.modalShapeViewer.shapeFig)   # the modal-shape panels: always a separate window
gui.angSweepViewer.toggleDrillDown!()   # open on the single-angle curve
# (same as the 'View single angle ->' button, top right; it toggles back to the map)
# A viewer, like a plot, only builds its figures; `display` opens them. The window that holds
# the shapes is a Figure of its own (`shapeFig`), so a combined viewer is two windows.

# Aside: the sign of the displacement field plotted above is reproducible run to run. The Ritz
# eigenvectors are put in a fixed real gauge right after the solve (the same convention as in
# MATLAB; `calcDisplacement` agrees with MATLAB's signed field to about 1e-8 of its scale), so
# the red/green stripes of the polarisation view keep their sense between runs.
