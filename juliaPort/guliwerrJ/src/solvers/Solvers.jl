"""
    Solvers

Aggregator for GULIWERR's `solvers/` layer: `Christoffel` (bulk-wave phase/group velocities
and limiting velocities), `Ritz` (the Ritz-Rayleigh SAW/Lamb eigensolver family) and
`Green` (the Green's-function surface-response solver).

One module per top-level folder's own submodule here (`Christoffel`, `Ritz`, `Green`), each
nested exactly two levels below `GuliwerrJ` (`GuliwerrJ.Solvers.Christoffel`, etc.) — see
`Christoffel.jl`'s own header comment for why that exact depth is load-bearing: its
`using ...CoreTypes: ...` climbs by lexical nesting, not by the `using` graph, so it only
resolves if `Solvers` sits directly under `GuliwerrJ`. Do not add another layer of nesting
here without updating every submodule's dot-count to match.

Each submodule's own `export` list only makes names reachable one level up (into `Solvers`
itself) — `using .Christoffel` here does NOT automatically re-export those names further up
into `GuliwerrJ`; this module must re-export them explicitly below, and `GuliwerrJ.jl` must
do the same again.
"""
module Solvers

include("christoffel/Christoffel.jl")
using .Christoffel

export matrix, GammaFromKron, phaseVels, phaseVelScan
export phaseVelAngleDeriv, groupVelVectorsForDir
export calc_anisotropy_ratio
export LimitingRoot, findRoots, refineRoot

include("ritz/Ritz.jl")
using .Ritz

export GammaMatrix, SAWBoundary, LambBoundary, calc,
       rawByDepth, derivsByDepth, elastByDepth, dummyForLayer,
       RawCellFor, LayerRotMat, RedElastCellFor, LayerNormFac, ContractDepth
export RitzSolver, solveConstrainedEig, realGauge, gaugeVector, isGaugeReal, phaseFix,
       REAL_GAUGE_TOL, precomp_field,
       buildConstraintMatrix, solveEigenvalueProblem,
       calcGroupVelocity, calcDisplacement, angleSlope, pruneEigVecs, ritzStructurePrecompute,
       ritzGramBlocks, gaugeBackMap
export RitzSolverSAW, RitzDiagnostics, modeDiagnostics, calcWithGamma, ritzSensitivities,
       ritzSensitivities!, RitzWorkspace
export RitzSolverLamb, modeRecognition, parityFractions, fastCalcGroupVel, locateZGV,
       calcTransGroupVelocity, kPolynomial, ZgvNewton

include("green/Green.jl")
using .Green

export GreenSolver, nvel, has_minus, has_half_space_full_output, has_layered_full_output
export GreenSexticPrepass, green_sextic_prepass, green_layer_setup, green_layers, GreenLayer,
       calc_detw_by_christoffel, solve_k3_allvel, canonicalize_k3,
       check_identical_records_allvel,
       calc_u_allvel, substrate_wave_indices, kvec_cols, outer_cols,
       OMEGA_TOLERANCE, ZERO_TOLERANCE, LIMITING_ABSORPTION_ETA,
       NEAR_ISO_ANGLE_STEP, NEAR_ISO_WARN_THRESH
export calc_traction_cols, calc_gfull_from_partial_waves, build_bc_matrix,
       green_component, green_gfull, DISP_FIELD_NDEPTH, DISP_FIELD_MAX_DEPTH, GREEN_CHUNK

end # module Solvers
