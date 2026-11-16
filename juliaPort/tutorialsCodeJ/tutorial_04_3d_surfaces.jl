## Tutorial 04 (code) -- directional bulk quantities: velocity/group-velocity/Young's surfaces and curves
# Narrated version: juliaPort/tutorialsJ/tutorial_04_3d_surfaces.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_04_3d_surfaces.jl")

using GuliwerrJ
using GLMakie

## Material -- Cu(001); no Structure or sweep needed
dirMat      = [1 0 0; 0 0 1]       # surface normal; propagation reference
material, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)   # constants + density from the materials library

## Mesh -- finer than the 61 x 121 default
nTheta = 91                        # polar-angle points
nPhi   = 181                       # azimuthal points

## Bulk phase-velocity surfaces (qL, qT1, qT2)
surfaceType = :velocity            # :velocity | :slowness
figVel = plotVelSurface(material, surfaceType, nTheta, nPhi)
display(figVel)                    # plot functions only build figures; display shows them

## Bulk group-velocity (ray) surfaces -- folds are real physics
figGroup = plotGroupVelSurface(material, nTheta, nPhi)
display(figGroup)

## Young's modulus surface
figYoungs = plotYoungsSurface(material, nTheta, nPhi)
content(figYoungs[1, 1]).title = "Cu (001): Young's modulus [GPa]"   # replaces the axis title
display(figYoungs)

## 2D in-plane cuts at the dirMat surface normal
inPlaneAnglesDeg = 0.0:1.0:360.0
axVelCurve = plotVelCurve(material, inPlaneAnglesDeg)
axVelCurve.title = "Cu (001): bulk phase velocity vs. in-plane angle"
display(axVelCurve.parent)
axYoungsCurve = plotYoungsCurve(material, inPlaneAnglesDeg)
axYoungsCurve.title = "Cu (001): Young's modulus vs. in-plane angle"
display(axYoungsCurve.parent)
