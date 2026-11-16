## Tutorial: 3D directional surfaces -- bulk waves and group velocity (Julia port)
#
# Same physics as examples/tutorials/tutorial_04_3d_surfaces.m: two self-contained 3D
# visualisations of the SAME anisotropic Material, each answering a different question about
# how bulk-wave behaviour depends on direction. Both sweep direction internally -- neither
# needs a Structure/CalcParams/sweep, only a Material. Exploratory/diagnostic use, not a
# pipeline step.
#
#   - plotVelSurface      -- bulk-wave (qL/qT1/qT2) phase velocity or slowness vs. WAVE-NORMAL
#                             direction (fast, pure Christoffel, no Ritz/Green)
#   - plotGroupVelSurface -- bulk-wave GROUP velocity (energy-propagation direction, which
#                             deviates from the wave normal) vs. wave-normal direction; can be
#                             genuinely "folded" for strongly anisotropic materials -- real
#                             physics, not a bug
#   - plotVelCurve / plotYoungsCurve -- the cheap 2D counterparts: one in-plane great circle
#                             at the dirMat surface normal instead of the whole sphere
#   - plotYoungsSurface   -- directional Young's modulus (one scalar field)
#
# SAW velocity has no analogous 3D-surface plot here: unlike a bulk wave, it isn't a function of
# direction alone -- it only exists once you fix a surface cut, and even then still varies with
# in-plane propagation azimuth on that cut, so it's a genuinely 2D-per-cut quantity, not a
# single-valued function of (theta,phi). See tutorial_02_angular_map.jl/plotAngularMap for the
# SAW velocity-vs-azimuth curve at one chosen surface cut (set via Material's dirMat).
#
# Copy this script as a starting point. Axis convention: plotVelSurface/plotGroupVelSurface use
# the SAME standard convention as Ritz/Christoffel/GreenSolver elsewhere in this project (axis 3
# = surface normal = dirMat row 1, axis 1 = dirMat row 2).
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - none new -- this script only reuses `Material`/`ElastConsts`/keyword construction from
#     tutorial 1 and plotting from tutorial 2. Worth noting explicitly: `plotVelSurface`/
#     `plotGroupVelSurface` return a whole `Figure` (several panels at once), not a single
#     `Axis` the way tutorial 2's angular-map plots did -- a function's return type in this
#     port always matches what it actually draws, one panel or several.

using GuliwerrJ
using GLMakie   # opens an interactive window; `using CairoMakie` instead renders the same
                 # Figure/Axis objects headlessly (e.g. for save("fig.png", fig)) -- see
                 # tutorial_02_angular_map.jl for the full explanation

## Material -- Cu(001), cubic (dirMat = Material's default, set explicitly)
dirMat   = [1 0 0;    # surface normal
            0 0 1]     # propagation reference (Material's default)

# The Cu constants and density come from the materials library (`ListMaterialLibrary(; print = true)`
# shows values and sources). `MaterialFromLibrary` returns the Material and the library row; the
# row is not needed here, so it is discarded with `_`.
mat, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)

## Bulk-wave phase-velocity surfaces (qL/qT1/qT2) -- cheap, pure Christoffel, no Ritz/Green
# involved. Mesh finer than the 61 x 121 default.
surfaceType = :velocity   # :velocity (default) | :slowness (= 1/velocity)
nTheta      = 91          # polar-angle mesh points (default 61)
nPhi        = 181         # azimuthal mesh points (default 121)
figVel = plotVelSurface(mat, surfaceType, nTheta, nPhi)
# `:velocity` is a Symbol, replacing MATLAB's `'velocity'` string argument (see tutorial 1's
# concept list).
Label(figVel[0, :], "Cu (001) -- bulk phase-velocity surfaces"; fontsize = 16)
display(figVel)   # a plot function only builds the figure; `display` shows it (tutorial 2)

## Bulk-wave group-velocity ("ray") surfaces -- energy-propagation direction,
# can genuinely differ in shape from the phase-velocity surfaces above.
figGroup = plotGroupVelSurface(mat, nTheta, nPhi)
Label(figGroup[0, :], "Cu (001) -- bulk group-velocity surfaces"; fontsize = 16)
display(figGroup)

## 2D cut: bulk phase velocities vs. in-plane angle on the (001) surface -- the Christoffel
# side of the SAW angular map (tutorial 2), sharing its angle axis.
inPlaneAngles = 0.0:1.0:360.0   # deg (default: the same full 1-deg circle)
axVelCurve = plotVelCurve(mat, inPlaneAngles)
axVelCurve.title = "Cu (001) -- bulk phase velocities vs. in-plane angle"
display(axVelCurve.parent)

## Directional Young's modulus -- 3D surface and the same 2D in-plane cut
figYoungs = plotYoungsSurface(mat, nTheta, nPhi)
content(figYoungs[1, 1]).title = "Cu (001) -- directional Young's modulus [GPa]"
# replaces the axis's own title (a Label above it would overlap); `content(fig[1, 1])` is the
# Axis3 in that grid cell.
display(figYoungs)

axYoungsCurve = plotYoungsCurve(mat, inPlaneAngles)
axYoungsCurve.title = "Cu (001) -- Young's modulus vs. in-plane angle"
display(axYoungsCurve.parent)
