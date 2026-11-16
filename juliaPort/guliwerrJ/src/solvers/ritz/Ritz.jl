"""
    Ritz

Aggregator for GULIWERR's Ritz-Rayleigh solver family: `GammaMatrix` (the static
Legendre-basis assembly pipeline), and `RitzSolver`/`RitzSolverSAW`/`RitzSolverLamb` (the
eigensolver family itself).

Nested exactly two levels below `GuliwerrJ` (`GuliwerrJ.Solvers.Ritz`), the same depth as
`GuliwerrJ.Solvers.Christoffel` — every file `include`d here that needs `CoreTypes` climbs
to it via `using ...CoreTypes: ...` (three dots), which only resolves at this exact nesting
depth. See `GammaMatrix.jl`'s and `Christoffel.jl`'s own header comments before changing
this module's position in the tree.

None of the files below declare their own `module` — see `CoreTypes.jl`'s header for why
(one module per top-level folder, not per file); `include` them directly, never `using`.
"""
module Ritz

include("GammaMatrix.jl")

export GammaMatrix, SAWBoundary, LambBoundary, calc,
       rawByDepth, derivsByDepth, elastByDepth, dummyForLayer,
       RawCellFor, LayerRotMat, RedElastCellFor, LayerNormFac, ContractDepth

include("RitzSolver.jl")
include("RitzSolverSAW.jl")
include("ZgvNewton.jl")      # module ZgvNewton (Newton refinement of ZGV points), used by RitzSolverLamb.jl
include("RitzSolverLamb.jl")
include("RitzWorkspace.jl")
include("RitzSensitivity.jl")

export RitzSolver, solveConstrainedEig, realGauge, gaugeVector, isGaugeReal, phaseFix,
       REAL_GAUGE_TOL, precomp_field,
       buildConstraintMatrix, solveEigenvalueProblem,
       calcGroupVelocity, calcDisplacement, angleSlope, pruneEigVecs
export RitzSolverSAW, RitzDiagnostics, modeDiagnostics, calcWithGamma
export RitzWorkspace, ritzSensitivities, ritzSensitivities!
export RitzSolverLamb, modeRecognition, parityFractions, fastCalcGroupVel, locateZGV,
       calcTransGroupVelocity, kPolynomial, ZgvNewton

end # module Ritz
