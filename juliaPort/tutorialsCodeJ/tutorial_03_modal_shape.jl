## Tutorial 03 (code) -- modal shapes: Ritz mode, Green's response, Green spectrum, viewer
# Narrated version: juliaPort/tutorialsJ/tutorial_03_modal_shape.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_03_modal_shape.jl")

using GuliwerrJ
using GLMakie

## Material -- Cu(001)
dirMat      = [1 0 0; 0 0 1]       # surface normal; propagation reference
material, _ = MaterialFromLibrary("Cu"; dirMat = dirMat)   # constants + density from the materials library

## Discretisation + structure -- "easy" tier
relDepth      = 20.0
legOrder      = 50
calcParams    = CalcParams(; relDepth = relDepth, legOrder = legOrder)
structureName = "cu_001_20deg"
nMaxModes     = 500
structure     = Structure(; materIn = material, calcParamsIn = calcParams,
                            name = structureName, nMaxModes = nMaxModes)

## Single-angle full sweep, keeping every eigenvector
angleDeg       = 20.0
calcMode       = :full
retentionLevel = :full             # plotModalShape needs the full eigenvector matrix
nVelGreen      = nothing           # Green-axis points (nothing = 6 m/s steps)
sweep = angSweepHS(; structure = structure, angles = [angleDeg], calcMode = calcMode,
                     retentionLevel = retentionLevel, nVelGreen = nVelGreen)
sweep = calc(sweep)

modeNo = sweep.sawModeNoMap[1][1]  # Ritz index of the first classified SAW
vSAW   = sweep.sawVelMap[1][1]     # km/s

## Ritz modal shape
ritzSolver      = sweep.cs[1].ritzSolver
substrateCutoff = :none            # :none | :decay | :fixedMultiple | :boosted
cutoffParam     = nothing          # decay fraction or depth multiple; unused for :none
dispDirs        = 1:3              # displacement components to draw
meshSpacing     = 0.01             # depth/lateral mesh spacing
axesModal = plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dispDirs, meshSpacing)
display(axesModal[1].parent)       # plot functions only build the figure; display shows it

## Same mode, one polarisation-RGB panel (red u_x, blue u_y, green u_z; brightness = |u|)
rgbDirs = :rgb                     # merged view instead of one panel per component
axesRgb = plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, rgbDirs, meshSpacing)
display(axesRgb[1].parent)

## Green's-response displacement field at the same velocity (half-space only)
fullOutput = true                  # fills .dispField, needed by plotGreenModalShape
greenField = calc(GreenSolver, structure, angleDeg, [vSAW]; fullOutput = fullOutput)
forceDir   = 3                     # 1 = propagation, 2 = in-plane transverse, 3 = normal
axesGreen  = plotGreenModalShape(greenField, forceDir)
display(axesGreen[1].parent)

## Green spectrum |G_IJ| vs. velocity, own axis and component
longLimVel    = sweep.cs[1].longLimVel                    # km/s
velMin        = 1.5                                       # km/s
nVel          = 1000
velAxisGreen  = collect(range(velMin, 1.05 * longLimVel; length = nVel))
forceDirI     = 3                  # force direction I (1-3)
displDirJ     = 3                  # displacement direction J (1-3); (1, 3) = G_13 (TGS)
greenSpectrum = calc(GreenSolver, structure, angleDeg, velAxisGreen;
                     forceDir = forceDirI, displDir = displDirJ)
figSpec = Figure()
axSpec  = Axis(figSpec[1, 1]; xlabel = "velocity [km/s]", ylabel = "|G_$(forceDirI)$(displDirJ)|",
               yscale = log10)
lines!(axSpec, velAxisGreen, abs.(greenSpectrum.G))
vlines!(axSpec, [vSAW]; linestyle = :dash, color = :black)   # Ritz SAW velocity
display(figSpec)

## Interactive viewer -- Green/Ritz curve linked to the modal shape
gui = combinedViewer(sweep)
display(gui.fig)                         # map/curve + controls
display(gui.modalShapeViewer.shapeFig)   # modal-shape panels (always a separate window)
gui.angSweepViewer.toggleDrillDown!()   # open on the single-angle curve
# (same as the 'View single angle ->' button, top right; it toggles back to the map)
