# =========================================================================================
# DispersionBranches.jl -- port of guliwerr/forwardCalcs/DispersionBranches.m: tracks each Ritz
# mode's identity across a fixed-angle wavelength sweep by MASS-WEIGHTED eigenvector overlap
# (MAC) with the previous step, assigned ONE-TO-ONE (a linear assignment maximising the total
# MAC), rather than an argmax per branch.
# `trackSeededBranch` seeds the branch at the first step's slowest accepted SAW.
#
# `trackUnits` holds the algorithm on the mass-weighted unit eigenvector matrices, golden-
# tested on MATLAB's own matrices; `track(sweepResult, angleIdx)` reads them off the sweep's
# per-step RitzSolvers.
# =========================================================================================

module DispersionBranches

using LinearAlgebra: norm
using Hungarian: hungarian
using ...Solvers: LayerNormFac
using ...CoreTypes: nlay

"""
    massWeightedUnit(ritz, nCompModes) -> Matrix{ComplexF64}

`X = sqrt(w) .* (N * eigVecs_con[:, 1:nCompModes])`, w = each Legendre coefficient's layer
factor (GammaMatrix.LayerNormFac, the diagonal of the unconstrained M), columns normalised.
`X_a' X_b` is then the overlap in the physical inner product integral(rho u^H u dz).
"""
function massWeightedUnit(ritz, nCompModes::Integer)
    s = ritz.structure
    w = reduce(vcat, [fill(LayerNormFac(s, il), s.calcParams[il].basisVecL) for il in 1:nlay(s)])
    X = sqrt.(w) .* (ritz.nullConstrMat * ritz.eigVecs_con[:, 1:nCompModes])
    for j in axes(X, 2)
        X[:, j] ./= norm(@view X[:, j])
    end
    return X
end

"""
    assignOneToOne(MAC) -> Vector{Int}

MATLAB `matchpairs(1 - MAC, n + 1)` on the square MAC matrix (rows = branches, columns =
candidates): the complete one-to-one assignment minimising the total `1 - MAC` (the
unmatched cost n + 1 exceeds any possible total, so every branch is matched). Returns the
candidate of every branch.
"""
function assignOneToOne(MAC::AbstractMatrix{<:Real})
    asg, _ = hungarian(1 .- MAC)
    return Vector{Int}(asg)
end

"""
    trackUnits(units) -> Matrix{Int}

`branchModeNoMap` ([nWav x nCompModes]) from the mass-weighted unit eigenvector matrices of
the steps: row 1 = 1:nCompModes, then each step's one-to-one MAC assignment, chaining the
reference forward to the matched columns.
"""
function trackUnits(units::AbstractVector{<:AbstractMatrix{<:Number}})
    nWav = length(units)
    nC = size(units[1], 2)
    map = zeros(Int, nWav, nC)
    map[1, :] = 1:nC
    ref = units[1]
    for iw in 2:nWav
        cand = units[iw]
        MAC = abs.(ref' * cand)
        matchIdx = assignOneToOne(MAC)
        map[iw, :] = matchIdx
        ref = cand[:, matchIdx]
    end
    return map
end

"""
    track(sweepResult, angleIdx = 1) -> Matrix{Int}

`DispersionBranches.track`: the branch map of one angle of a CalcWavevecAngSweep2Layered /
CalcWavevecAngSweep2Plate result computed with retentionLevel 'full' (every step keeps its
full eigenvector matrix). Column im is one tracked branch; entry (iw, im) is the raw Ritz
mode index at step iw; every row is a permutation of 1:nCompModes.
"""
function track(sweepResult, angleIdx::Integer = 1)
    nA = length(sweepResult.angles)
    1 <= angleIdx <= nA || throw(ArgumentError("DispersionBranches.track: angleIdx must be " *
        "an integer in 1:$nA, got $angleIdx"))
    nWav = length(sweepResult.sweeps)
    ritz1 = sweepResult.sweeps[1].cs[angleIdx].ritzSolver
    nC = ritz1.nCompModes
    units = Vector{Matrix{ComplexF64}}(undef, nWav)
    for iw in 1:nWav
        rz = sweepResult.sweeps[iw].cs[angleIdx].ritzSolver
        isempty(rz.eigVecs_con) && throw(ArgumentError("DispersionBranches:track:missingEigVecs" *
            " -- eigVecs_con is empty at wavelength step $iw; DispersionBranches.track needs " *
            "the full Ritz eigenvector matrix at every step (retentionLevel = 'full')."))
        units[iw] = massWeightedUnit(rz, nC)
    end
    return trackUnits(units)
end

"""
    trackSeededBranch(sweepResult, angleIdx = 1) -> (vels, modeNos, seedModeNo)

`DispersionBranches.trackSeededBranch`: the velocity of the branch seeded at
the accepted `:SAW` entry with the slowest bound Ritz mode (`modeVel`) of the first wavelength step
(`sweepResult.sweeps[1].significantModesMap[angleIdx]`, i.e. the classification under the
sweep's active method and standard), followed through the sweep by [`track`](@ref).
`vels` [nWav] Ritz velocity of that branch at every step [km/s], `modeNos` [nWav] its raw
Ritz mode index, `seedModeNo` the seed entry's modeNo. Needs a CalcWavevecAngSweep2Layered
computed with retentionLevel 'full'. Throws `"DispersionBranches:trackSeededBranch:noSeed"`
if the first step has no accepted 'SAW' entry at that angle.
"""
function trackSeededBranch(sweepResult, angleIdx::Integer = 1)
    branchModeNoMap = track(sweepResult, angleIdx)
    sm = [m for m in sweepResult.sweeps[1].significantModesMap[angleIdx] if m.type === :SAW]
    isempty(sm) && throw(ArgumentError("DispersionBranches:trackSeededBranch:noSeed -- The " *
        "first wavelength step has no accepted 'SAW' entry at angle index $angleIdx to seed " *
        "the branch from."))
    seedModeNo = sm[argmin([m.modeVel for m in sm])].modeNo   # bound Ritz mode, the tracked quantity (vel is the Green peak)
    modeNos = branchModeNoMap[:, seedModeNo]
    vels = [sweepResult.sweeps[iw].cs[angleIdx].ritzSolver.vels[modeNos[iw]]
            for iw in eachindex(modeNos)]
    return vels, modeNos, seedModeNo
end

end # module DispersionBranches
