# =========================================================================================
# InverseCalcs.jl -- the inverse module (MATLAB guliwerr/inverseCalcs/), one module for the
# folder: MeasuredVelocities, ModeAssociation, SensitivityAnalysis, ElasticityInversion and
# the Inv launchers are `include`d as plain top-level definitions (each file's main type has
# the file's name, so no per-file module). Sits one level below GuliwerrJ, like ForwardCalcs:
# siblings are reached with two dots (`..Solvers`). `SensitivityAnalysis`'s `calc` extends
# the shared `calc` generic (hence `import ..Solvers: calc`).
#
# EXPERIMENTAL (as in MATLAB): validated on synthetic data only; the interface may change.
# =========================================================================================

module InverseCalcs

using LinearAlgebra: svd, Diagonal, norm, I
using Statistics: median

using ..CoreTypes: Structure, Material, nlay, ElastConsts, SetElastConsts, SetDirMat,
                    GetRotDirMat, GetStandardDirMat, GetRotElastTensor, Any2Tensor, Voigt2Bulk,
                    with_blas_threads, with_task_blas_single, matchpairs0, ModifyLayer

import ..Solvers: calc
using ..Solvers: dummyForLayer, matrix, refineRoot, LimitingRoot, SAWBoundary

using ..ForwardCalcs: limitingAndBulk, resolveParallel, rotateStructure, angSweepHS,
                      angSweepLayered

include("MeasuredVelocities.jl")
include("ModeAssociation.jl")
include("SensitivityAnalysis.jl")
include("ElasticityInversion.jl")
include("Inv.jl")

export WAVE_TYPES, MeasuredVelocities, PhaseVelocityEntry, addPhaseVelocity, setBulkModulus,
       fromTable
export ModeAssociation, AssocPoint, LIM_VELOCITY_WINDOW, LIM_THETA_WINDOW, LIM_BRANCH_PENALTY,
       RITZ_CLASSES, DEFAULT_TRACK_OPTS, DEFAULT_REBIND_WINDOW, activeRitzClasses,
       nPoints, uniqueRitzAngles, uniqueLimitingAngles, pointsAtAngle, confidence,
       pinnedPoints, setPinned, setLimitingAnchor, setLimitingRefined, setLimitingLost,
       setGreenPeak, addPoint, removeDuplicatePoints, rebind, rebindLimiting,
       fromSweep, defaultTrackOpts, limitingRootScore
export SensitivityAnalysis, DEFAULT_SENSITIVITY_CUTOFF, calc, waveTypeIndex, phaseVelocityAt,
       bulkVelocityAt, limitingStep, phaseBulkStep, limitingRootsAt, limitingVelocityNear,
       limitingVelocityTracked, rotatedTensorFor, channelDerivative, channelValue,
       ritzVelocityDerivative, checkDatasets
export ElasticityInversion, calc!, recalcForwardAtResult, pseudoSAWGreenPeaks
export measured, association, sensitivity, elasticity

end # module InverseCalcs
