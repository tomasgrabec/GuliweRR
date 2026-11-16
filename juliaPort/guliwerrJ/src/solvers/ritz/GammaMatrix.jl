# =========================================================================================
# GammaMatrix.jl
#
# Port of guliwerr/solvers/ritz/@GammaMatrix/ (5 files): GammaMatrix.m (the type + `calc`/
# `calcLamb` instance methods + six `methods (Static)` helpers), rawByDepth.m, elastByDepth.m,
# derivsByDepth.m, dummyForLayer.m.
#
# Three-step static pipeline (Stoklasova 2015, Grabec 2024):
#   Step 1  rawByDepth    -- [bVL x bVL x 15 x 3] Legendre-integral tensor, material-independent
#   Step 2  derivsByDepth -- dK/dc, d^2K/(dc dD) cells for the angle-sweep (SAW) path
#   Step 3  elastByDepth  -- K0/K1/K2-with-elasticity cells for the frequency-sweep (Lamb) path
# `calc` assembles K/dK-dD/M for the SAW path; MATLAB's `calcLamb` assembles the same for the
# Lamb path -- both become ONE Julia `calc`, dispatching on a boundary-type singleton (see
# TRANSLATION DECISIONS note 2 below).
#
# NO per-file module (see CoreTypes.jl's / Christoffel.jl's headers): this file is
# `include`d directly into `module Ritz` (see Ritz.jl). No `module GammaMatrix` wrapper -- see
# ARCHITECTURAL NOTE 1 below for why that would break exactly like `CalcParams`/`ElastConsts`/
# `Material`/`Structure` would.
#
# MODULE NESTING / CROSS-MODULE DEPENDENCY (same load-bearing assumption as in
# Christoffel.jl's header): this file needs `Structure`, `CalcParams`, `Material`,
# `ElastConsts` and several of their free functions from `CoreTypes`, which is a SIBLING of
# `Ritz` (this file's parent module), not an ancestor reachable by plain lexical nesting.
# `Ritz` sits exactly TWO levels below `GuliwerrJ` (`GuliwerrJ.Solvers.Ritz`, exactly like
# `GuliwerrJ.Solvers.Christoffel`): from inside `module Ritz`, `..` reaches `Solvers` (Ritz's
# direct parent) and `...` reaches `GuliwerrJ` (Ritz's grandparent), where `CoreTypes` is
# bound directly (via `GuliwerrJ.jl`'s own `include("core/CoreTypes.jl")`), so
# `using ...CoreTypes: ...` (three dots) resolves. **If `Ritz.jl` is nested at a different
# depth than exactly two levels below `GuliwerrJ`, this `using ...CoreTypes` line needs its
# dot-count adjusted to match.**
#
# ARCHITECTURAL NOTE 1 -- GammaMatrix IS a result-bearing value type, not a namespace:
# `classdef GammaMatrix` (GammaMatrix.m) has NO `< handle` -- a MATLAB *value* class, exactly
# like `CalcParams`/`ElastConsts`/`Material`/`Structure`, NOT a static-namespace class like
# `Christoffel`/`ParallelPolicy`. Its properties (`gammaMatAll`, `dGammaDMatAll`,
# `normMatAll`) are genuine per-instance data (the assembled K, dK/dD, M matrices for one
# structure+angle). So this file declares `struct GammaMatrix` directly in `Ritz`'s
# namespace, NOT `module GammaMatrix` -- a same-named module would collide with the struct
# exactly as `module CalcParams` containing `struct CalcParams` would (CoreTypes.jl's header
# has the full writeup): `using .GammaMatrix` one level up could never rebind the
# auto-bound module name to the struct, and bare `GammaMatrix(...)` would fail with `objects
# of type Module are not callable`.
#
# ARCHITECTURAL NOTE 2 -- MUTABLE vs IMMUTABLE:
# In the MATLAB source, every write to `gammaMatAll`/`dGammaDMatAll`/`normMatAll` sits inside
# GammaMatrix.m's own `calc`/`calcLamb` methods (the `obj.field = value; ...; end` sequence a
# MATLAB value-class method uses before implicitly returning the (by-value) `obj`). Every
# caller (`RitzSolverSAW.m`, `RitzSolverLamb.m`) uses the `obj.gammaMat = obj.gammaMat.calc(...)`
# rebuild-and-return idiom, never `calc` for a mutating side effect; `ElasticityInversion.m`
# only ever *constructs* fresh `GammaMatrix()` objects. No `GammaMatrix` instance is mutated
# anywhere else in the MATLAB code.
# Conclusion: **immutable struct** (with mutable ARRAY fields -- cheap to pass, arrays
# written in place; the struct wrapper itself is never mutated in place), the same choice as
# for every other value class in the port.
#
# TRANSLATION DECISIONS:
# 1. **No "self" argument.** MATLAB's `obj = calc(obj, structure, angleDeg, ...)` takes a
#    receiver `obj`, but `calc`/`calcLamb` never actually READ any field of `obj` before
#    unconditionally overwriting all of them -- `obj` is a stateless placeholder receiver
#    required only by MATLAB's method-call syntax (`gm.calc(...)`), carrying no information
#    into the computation. Since Julia dispatches on all arguments and needs no receiver,
#    and `GammaMatrix` is immutable (so there is nothing useful an "old" instance could
#    contribute to building a "new" one here), the Julia `calc` simply takes the real inputs
#    (`structure`, `angleDeg`, a boundary-type singleton, and the optional precomputed cells
#    as keywords) and returns a fresh `GammaMatrix` -- there is no dummy
#    `gm0 = GammaMatrix(); gm = calc(gm0, ...)` step.
# 2. **`calc`/`calcLamb` unify into one `calc`, dispatching on a boundary-type singleton**
#    rather than two names -- the same pattern as `RitzSolver.buildConstraintMatrix`'s
#    `bcType = 'saw'|'lamb'` string switch. `SAWBoundary`/`LambBoundary` are defined HERE so
#    that `RitzSolver.jl`/`RitzSolverSAW.jl`/`RitzSolverLamb.jl` (`include`d into the same
#    `module Ritz`) reuse the identical singleton types for their own `buildConstraintMatrix`
#    dispatch with no further definition or import. `calc`'s positional argument order is
#    `(structure, ::SAWBoundary/LambBoundary, angleDeg=0.0; kwargs...)` -- the dispatch tag
#    must precede `angleDeg` positionally, since Julia requires every positional argument
#    with a default to come after every one without a default, and the singleton has (and
#    needs) no default.
# 3. **`GetReducedElastConsts` is reconstructed from `ElastRotReductMat`, not ported as its own
#    function.** `ElastConsts.jl` leaves `GetReducedElastConsts` (an object-facing wrapper
#    method, `ElastConsts.m`) out; its replacement is
#    `ElastRotReductMat(rotMat, ec.symmClassNo) * ec.vector`, as in `ElastConsts.m`:
#    `indep2RedElastMat = ElastConsts.ElastRotReductMat(rotMat, obj.symmClassNo);
#    redElast = indep2RedElastMat * obj.vector';`. `RedElastCellFor` below uses exactly that
#    expression (transpose is a no-op here since `ec.vector` is already a plain Julia column
#    `Vector`).
# 4. **No `squeeze()` needed anywhere** -- a genuine language difference, not a style choice.
#    `elastByDepth.m`'s own comment calls `squeeze()` "load-bearing": MATLAB's
#    `rawByDepthCell{ilay}(:,:,iC,:)` keeps a singleton 3rd dimension (`[bVL x bVL x 1 x 3]`)
#    that `squeeze` must strip before adding to a `[bVL x bVL x 3]` accumulator, because MATLAB
#    silently broadcasts otherwise-mismatched trailing singleton dimensions (that very bug is
#    what the MATLAB comment warns a future editor away from). Julia's `getindex` with a
#    SCALAR index (not a range/colon) at a given dimension always DROPS that dimension from the
#    result -- `raw[:, :, iC, :]` is already `bVL x bVL x 3` with no singleton to strip. The
#    arithmetic is the same; the MATLAB hazard simply does not exist in the Julia port.
# 5. **`blkdiag(A, B)` -> `cat(A, B; dims = (1, 2))`.** This matches MATLAB's
#    `blkdiag([], A) == A` accumulator-seeding idiom:
#    `cat(zeros(ComplexF64, 0, 0), A; dims = (1, 2)) == A`. `nLay` is 1-3 in practice, so the
#    O(nLay) reallocation this causes (a fresh matrix on every layer) is not worth trading
#    away for a hand-rolled offset/preallocate scheme at the cost of an extra code path.
# 6. **In-place `for iC = 1:15` accumulation:** `gammaMat_lay`/`dGammaDMat_lay`/`gByD`/`dG_lay`
#    are preallocated once per layer and accumulated into with `.+=`/`@views`, never rebuilt
#    from a fresh `zeros(...)` inside the `for iC in 1:15` loop body the way a literal
#    `A = A + B` translation would.
# =========================================================================================

using LinearAlgebra

using ...CoreTypes: Structure, CalcParams, Material, ElastConsts, nlay, str_basis_vec_l,
                     GetIndOf2DReduct_15of21, Any2Tensor, ElastRotReductMat,
                     GetStandardDirMat, CalcRotMat

export GammaMatrix, SAWBoundary, LambBoundary, calc,
       rawByDepth, derivsByDepth, elastByDepth, dummyForLayer,
       RawCellFor, LayerRotMat, RedElastCellFor, LayerNormFac, ContractDepth

# -----------------------------------------------------------------------------------------
# Boundary-type singletons -- see TRANSLATION DECISIONS note 2 above. Zero fields; existence
# alone selects a `calc` method.
# -----------------------------------------------------------------------------------------

"""
    SAWBoundary

Dispatch tag for the angle-sweep / SAW-and-layered-stack boundary condition (MATLAB's
`GammaMatrix.calc`, and `RitzSolver.buildConstraintMatrix`'s `bcType = 'saw'`).
"""
struct SAWBoundary end

"""
    LambBoundary

Dispatch tag for the frequency-sweep / free-plate (Lamb) boundary condition (MATLAB's
`GammaMatrix.calcLamb`, and `RitzSolver.buildConstraintMatrix`'s `bcType = 'lamb'`).
"""
struct LambBoundary end

# -----------------------------------------------------------------------------------------
# The type
# -----------------------------------------------------------------------------------------

"""
    GammaMatrix

Stiffness (K), mass (M) and dK/dD matrices for the Ritz-Rayleigh solver, assembled
block-diagonally across a `Structure`'s layers.

Port of `guliwerr/solvers/ritz/@GammaMatrix/GammaMatrix.m`. **Immutable struct** -- see the
file-level ARCHITECTURAL NOTE 2 comment in `GammaMatrix.jl` for why (the MATLAB source never
mutates an instance outside `calc`/`calcLamb`).

Not `Solvers.Christoffel`'s `matrix`/acoustic-tensor `Gamma(i,k)` -- same name in the MATLAB
source (`GammaMatrix.m`'s own disambiguating comment), different physics object.

# Fields
- `gammaMatAll::Matrix{ComplexF64}` -- K matrix (density-normalised, block-diagonal across
  layers). Complex because the `K1` depth-integral block is purely imaginary
  (`1i * (antisymmetric)`, see [`rawByDepth`](@ref)'s docstring) even though `K0`/`K2` are
  real; MATLAB's dynamically-typed array is promoted to `complex double` the first time a
  complex value is assigned into it, whole-array, hence `ComplexF64` here from construction
  rather than a real type that later needs promoting.
- `dGammaDMatAll::Matrix{ComplexF64}` -- dK/dD matrix (for group velocity), same block
  structure and complex-for-the-same-reason as `gammaMatAll`.
- `normMatAll::Matrix{Float64}` -- M matrix: block-diagonal, each block a real scaled
  identity (`eye(bVL) .* normFac`) -- genuinely real, never promoted, unlike the two fields
  above. `maximum(abs.(normMatAll - I))` is `< 1e-15` for a single layer (or a multi-layer
  stack with identical `relDepth*rho` per layer) and `> 1e-3` once layers differ in
  `relDepth*rho` -- the trigger for `RitzSolver`'s generalised-eigenproblem branch.

Construct via the zero-argument default (an empty placeholder, matching
`GammaMatrix()` in `GammaMatrix.m`) or via [`calc`](@ref). The Ritz solvers build one per
solve as a local and do not keep it (see [`RitzSolver`](@ref)).
"""
struct GammaMatrix
    gammaMatAll::Matrix{ComplexF64}
    dGammaDMatAll::Matrix{ComplexF64}
    normMatAll::Matrix{Float64}
end

GammaMatrix() = GammaMatrix(zeros(ComplexF64, 0, 0), zeros(ComplexF64, 0, 0),
                             zeros(Float64, 0, 0))

# -----------------------------------------------------------------------------------------
# ContractDepth.m -- shared by derivsByDepth (SAW path) and calc(..., ::LambBoundary, ...)
# (Lamb path); they differ in what the arrays hold (raw Legendre tensor vs
# elasticity-contracted K0/K1/K2), not in the depth convention. Generic over array shape
# (works on the 2-D per-layer matrices in the Lamb path and the 3-D per-constant stacks in
# the SAW path alike), matching MATLAB's own shape-polymorphic use of this function.
# -----------------------------------------------------------------------------------------

"""
    ContractDepth(K0, K1, K2, relDepth) -> (K, dKdD)

Apply the depth factors to a `(K0, K1, K2)` triple:

    K     = K0/dd1^2 + K1/(dd1*dd2/2) + K2/(dd2/2)^2
    dK/dD = -2*K1/(dd1*dd2^2) - 8*K2/dd2^3

with `dd1 = 1` permanently and `dd2 = relDepth` (Stoklasova 2015, Eq. (A1)-(A4)).

Port of `guliwerr/solvers/ritz/@GammaMatrix/GammaMatrix.m` (`GammaMatrix.ContractDepth`,
a `methods (Static)` helper, not a separate file).
"""
function ContractDepth(K0, K1, K2, relDepth::Real)
    K    = K0 .+ K1 ./ (relDepth / 2) .+ K2 ./ (relDepth / 2)^2
    dKdD = -2 .* K1 ./ relDepth^2 .- 8 .* K2 ./ relDepth^3
    return K, dKdD
end

# -----------------------------------------------------------------------------------------
# rawByDepth.m -- Step 1
# -----------------------------------------------------------------------------------------

"""
    rawByDepth(calcParams::CalcParams) -> Array{ComplexF64,4}

Step 1: material-independent Legendre-integral tensor, `bVL x bVL x 15 x 3`.
Ref: Stoklasova 2015, Eq. (A1)-(A4); Grabec 2024.

Dim 3 indexes the 15 reduced elastic constants ([`GetIndOf2DReduct_15of21`](@ref) order);
dim 4 is the depth-factor group: `(:,:,:,1) = dK0/dc` (factor `1/dd1^2 = 1`), `(:,:,:,2) =
dK1/dc` (factor `1/(dd1*dd2/2)`), `(:,:,:,3) = dK2/dc` (factor `1/(dd2/2)^2`).

`K0`/`K2` slices are real-symmetric; the `K1` slice is Hermitian (purely imaginary --
`1i * (antisymmetric real matrix)`). All three are stored in the SAME `ComplexF64` array
(MATLAB's own array, `zeros(bVL,bVL,nRedEl,3)`, gets promoted whole once the `K1` slice's
`1i * (...)` assignment lands on it -- ported here as `ComplexF64` from construction rather
than promoted after the fact, since Julia arrays do not change element type in place).

Port of `guliwerr/solvers/ritz/@GammaMatrix/rawByDepth.m`.
"""
function rawByDepth(calcParams::CalcParams)
    legO = calcParams.legOrder
    bVL  = calcParams.basisVecL
    legD = calcParams.legDiag
    legN = calcParams.legNondiag
    delta = Matrix{Float64}(I, legO, legO)

    indRedEl, nRedEl = GetIndOf2DReduct_15of21()
    rawTensor = zeros(ComplexF64, bVL, bVL, nRedEl, 3)

    for iconst in 1:nRedEl
        consts = zeros(Float64, 21)
        consts[indRedEl[iconst]] = 1.0
        cTen = Any2Tensor(consts, 10)   # triclinic dummy (symmNo = 10)

        triu_K0  = zeros(Float64, bVL, bVL)
        triu_K1a = zeros(Float64, bVL, bVL)
        triu_K1b = zeros(Float64, bVL, bVL)
        triu_K2  = zeros(Float64, bVL, bVL)

        for ind1 in 1:bVL
            iL1 = mod(ind1 - 1, legO)     # 0-indexed Legendre order of row
            iD1 = cld(ind1, legO)         # direction index 1-3 (ceil division)
            for ind2 in ind1:bVL          # upper triangle only
                iL2 = mod(ind2 - 1, legO)
                iD2 = cld(ind2, legO)

                triu_K0[ind1, ind2]  = cTen[iD1, 1, iD2, 1] * (2π)^2 * delta[iL1 + 1, iL2 + 1]
                triu_K1a[ind1, ind2] = cTen[iD1, 1, iD2, 3] *  (2π)  * legN[iL2 + 1, iL1 + 1]
                triu_K1b[ind1, ind2] = cTen[iD1, 3, iD2, 1] * (-2π)  * legN[iL1 + 1, iL2 + 1]
                triu_K2[ind1, ind2]  = cTen[iD1, 3, iD2, 3]          * legD[iL1 + 1, iL2 + 1]
            end
        end

        triu_K1 = triu_K1a .+ triu_K1b

        # symmetrise: K0 and K2 are real-symmetric; K1 is Hermitian (1i x antisymmetric).
        # `transpose`, not `'`, mirrors MATLAB's `.'` here -- these are still plain Float64
        # matrices at this point (the `im *` below is the only complex-promoting step).
        rawTensor[:, :, iconst, 1] = triu_K0 .+ transpose(triu(triu_K0, 1))
        rawTensor[:, :, iconst, 2] = im .* (triu_K1 .- transpose(triu(triu_K1, 1)))
        rawTensor[:, :, iconst, 3] = triu_K2 .+ transpose(triu(triu_K2, 1))
    end
    return rawTensor
end

# -----------------------------------------------------------------------------------------
# RawCellFor -- MATLAB `methods (Static)`, GammaMatrix.m
# -----------------------------------------------------------------------------------------

"""
    RawCellFor(structure::Structure) -> Vector{Array{ComplexF64,4}}

Step-1 Legendre-integral tensor ([`rawByDepth`](@ref)) for every layer of a structure. The
default-fill used by [`derivsByDepth`](@ref) and [`elastByDepth`](@ref) when no pre-built
cell is handed in. Public so a caller sweeping a structure at many wavelengths/angles can
build it once and pass it back down (as MATLAB's `CalcWavevecAngSweep2Layered` does).

Port of `guliwerr/solvers/ritz/@GammaMatrix/GammaMatrix.m`.
"""
function RawCellFor(structure::Structure)
    n = nlay(structure)
    rawCell = Vector{Array{ComplexF64,4}}(undef, n)
    for ilay in 1:n
        rawCell[ilay] = rawByDepth(structure.calcParams[ilay])
    end
    return rawCell
end

# -----------------------------------------------------------------------------------------
# derivsByDepth.m -- Step 2
# -----------------------------------------------------------------------------------------

"""
    derivsByDepth(structure::Structure, rawByDepthCell = nothing)
        -> (dKdC_cell, d2KdCdD_cell)

Step 2: contract depth factors ([`ContractDepth`](@ref)) to get dK/dc and d^2K/(dc dD) per
layer, each `bVL x bVL x 15` (`K = sum(dK/dc * c)`, `dK/dD = sum(d^2K/(dc dD) * c)`).

`rawByDepthCell` -- from [`RawCellFor`](@ref); computed if omitted.

Port of `guliwerr/solvers/ritz/@GammaMatrix/derivsByDepth.m`.
"""
function derivsByDepth(structure::Structure, rawByDepthCell = nothing)
    rawCell = rawByDepthCell === nothing ? RawCellFor(structure) : rawByDepthCell

    n = nlay(structure)
    dKdC_cell    = Vector{Array{ComplexF64,3}}(undef, n)
    d2KdCdD_cell = Vector{Array{ComplexF64,3}}(undef, n)

    for ilay in 1:n
        raw = rawCell[ilay]   # bVL x bVL x 15 x 3
        bVL = size(raw, 1)
        dKdC_cell[ilay]    = Array{ComplexF64,3}(undef, bVL, bVL, 15)
        d2KdCdD_cell[ilay] = Array{ComplexF64,3}(undef, bVL, bVL, 15)
        _contract_depth_rho!(dKdC_cell[ilay], d2KdCdD_cell[ilay], raw,
                             structure.calcParams[ilay].relDepth, structure.mater[ilay].rho)
    end
    return dKdC_cell, d2KdCdD_cell
end

"""
    _contract_depth_rho!(dKdC, d2KdCdD, raw, relDepth, rho)

[`ContractDepth`](@ref) of the raw tensor's three depth groups, divided by `rho`, written into
`dKdC`/`d2KdCdD` in one pass without the slice copies and intermediate arrays. Each element
goes through the same operations in the same order as `ContractDepth(...) ./ rho`, so the
result is bit-identical.
"""
function _contract_depth_rho!(dKdC, d2KdCdD, raw, relDepth::Real, rho::Real)
    h  = relDepth / 2
    h2 = (relDepth / 2)^2
    D2 = relDepth^2
    D3 = relDepth^3
    @inbounds for i in CartesianIndices(dKdC)
        k1 = raw[i, 2]
        k2 = raw[i, 3]
        dKdC[i]    = (raw[i, 1] + k1 / h + k2 / h2) / rho
        d2KdCdD[i] = (-2 * k1 / D2 - 8 * k2 / D3) / rho
    end
    return dKdC, d2KdCdD
end

# -----------------------------------------------------------------------------------------
# elastByDepth.m -- Step 3 (Lamb path)
# -----------------------------------------------------------------------------------------

"""
    elastByDepth(structure::Structure, angleDeg::Real = 0.0,
                 rawByDepthCell = nothing, redElastCell = nothing)
        -> Vector{Array{ComplexF64,3}}

Step 3 (Lamb path): contract [`rawByDepth`](@ref) with elasticity, keeping the depth factors
separate. Returns one `bVL x bVL x 3` array per layer (`K0`, `K1`, `K2` WITH elasticity,
WITHOUT rho/depth).

**No `squeeze()` needed here, unlike MATLAB** -- see this file's TRANSLATION DECISIONS note 4.

Port of `guliwerr/solvers/ritz/@GammaMatrix/elastByDepth.m`.
"""
function elastByDepth(structure::Structure, angleDeg::Real = 0.0,
                       rawByDepthCell = nothing, redElastCell = nothing)
    n = nlay(structure)
    if rawByDepthCell === nothing
        rawCell = RawCellFor(structure)
    else
        _check_cell_length(rawByDepthCell, n, "elastByDepth", "rawByDepthCell")
        rawCell = rawByDepthCell
    end
    if redElastCell === nothing
        redCell = RedElastCellFor(structure, angleDeg)
    else
        _check_cell_length(redElastCell, n, "elastByDepth", "redElastCell")
        redCell = redElastCell
    end

    gammaByDepthCell = Vector{Array{ComplexF64,3}}(undef, n)
    for ilay in 1:n
        bVL = structure.calcParams[ilay].basisVecL
        gByD = zeros(ComplexF64, bVL, bVL, 3)
        for iC in 1:15
            @views gByD .+= rawCell[ilay][:, :, iC, :] .* redCell[ilay][iC]
        end
        gammaByDepthCell[ilay] = gByD
    end
    return gammaByDepthCell
end

# -----------------------------------------------------------------------------------------
# LayerRotMat / RedElastCellFor / LayerNormFac -- MATLAB `methods (Static)`
# -----------------------------------------------------------------------------------------

"""
    LayerRotMat(structure::Structure, ilay::Integer, angleDeg::Real) -> Matrix{Float64}

Rotation matrix taking layer `ilay`'s crystal frame to the standard solver frame
(x1 = propagation, x2 = in-plane, x3 = depth) at `angleDeg`.

[`GetStandardDirMat`](@ref) applies the propagation angle within the layer's own surface
plane and reorders to the standard convention, so [`CalcRotMat`](@ref) takes angle `0` --
passing `angleDeg` to BOTH would double-rotate (`docs/pitfalls.md` OTH1).

Port of `guliwerr/solvers/ritz/@GammaMatrix/GammaMatrix.m`.
"""
function LayerRotMat(structure::Structure, ilay::Integer, angleDeg::Real)
    stdDirMat = GetStandardDirMat(structure.mater[ilay].dirMat, angleDeg)
    return CalcRotMat(stdDirMat, 0.0)
end

"""
    RedElastCellFor(structure::Structure, angleDeg::Real) -> Vector{Vector{Float64}}

Rotated reduced (15-of-21) elastic constants for every layer. The default-fill shared by
[`calc`](@ref)'s SAW-path method and [`elastByDepth`](@ref).

Port of `guliwerr/solvers/ritz/@GammaMatrix/GammaMatrix.m`, with
`ElastConsts.GetReducedElastConsts(rotMat)` expanded to
`ElastRotReductMat(rotMat, ec.symmClassNo) * ec.vector` -- see this file's TRANSLATION
DECISIONS note 3.
"""
function RedElastCellFor(structure::Structure, angleDeg::Real)
    n = nlay(structure)
    redElastCell = Vector{Vector{Float64}}(undef, n)
    for ilay in 1:n
        rotMat = LayerRotMat(structure, ilay, angleDeg)
        ec = structure.mater[ilay].elConsts
        redElastCell[ilay] = ElastRotReductMat(rotMat, ec.symmClassNo) * ec.vector
    end
    return redElastCell
end

"""
    LayerNormFac(structure::Structure, ilay::Integer) -> Float64

Per-layer normalisation factor of the block-diagonal K/M assembly:
`(relDepth_ilay * rho_ilay) / (relDepth_1 * rho_1)` -- layer 1 (the substrate / bottom layer)
is the reference, so its own factor is exactly `1`.

Port of `guliwerr/solvers/ritz/@GammaMatrix/GammaMatrix.m`.
"""
function LayerNormFac(structure::Structure, ilay::Integer)
    return (structure.calcParams[ilay].relDepth * structure.mater[ilay].rho) /
           (structure.calcParams[1].relDepth    * structure.mater[1].rho)
end

# -----------------------------------------------------------------------------------------
# calc / calcLamb -> ONE calc, dispatching on SAWBoundary/LambBoundary -- see TRANSLATION
# DECISIONS notes 1-2 above.
# -----------------------------------------------------------------------------------------

"""
    calc(structure::Structure, ::SAWBoundary, angleDeg::Real = 0.0;
         dKdC_cell = nothing, d2KdCdD_cell = nothing, redElastCell = nothing) -> GammaMatrix

Assemble K, dK/dD and M for the angle-sweep (SAW) path.

`dKdC_cell`/`d2KdCdD_cell` -- from [`derivsByDepth`](@ref); supplied TOGETHER or not at all
(both `nothing`: computed). A half-supplied pair, or a cell without one entry per layer, is an
`ArgumentError` `GammaMatrix:calc:badPrecompute` -- not a silent rebuild (as in MATLAB).
`redElastCell` -- rotated reduced elastic constants per layer ([`RedElastCellFor`](@ref));
computed if `nothing`, same length check.

Port of `GammaMatrix.calc`. The block-diagonal matrices are preallocated and filled block by
block, and each layer's `sum_c c * dK/dc` is one matrix-vector product over the 15 constants.
"""
function calc(structure::Structure, ::SAWBoundary, angleDeg::Real = 0.0;
              dKdC_cell = nothing, d2KdCdD_cell = nothing, redElastCell = nothing)
    n = nlay(structure)
    if dKdC_cell === nothing && d2KdCdD_cell === nothing
        dKdC_cell, d2KdCdD_cell = derivsByDepth(structure)
    elseif dKdC_cell === nothing || d2KdCdD_cell === nothing
        throw(ArgumentError("GammaMatrix:calc:badPrecompute -- dKdC_cell and d2KdCdD_cell " *
                            "must be supplied together (or both left empty)."))
    else
        _check_cell_length(dKdC_cell, n, "calc", "dKdC_cell")
        _check_cell_length(d2KdCdD_cell, n, "calc", "d2KdCdD_cell")
    end
    if redElastCell === nothing
        redElastCell = RedElastCellFor(structure, angleDeg)
    else
        _check_cell_length(redElastCell, n, "calc", "redElastCell")
    end

    strBVL = sum(structure.calcParams[i].basisVecL for i in 1:n; init = 0)
    gMatAll       = zeros(ComplexF64, strBVL, strBVL)
    gMatAll_dD    = zeros(ComplexF64, strBVL, strBVL)
    nMatAll       = zeros(Float64, strBVL, strBVL)

    r0 = 0
    for ilay in 1:n
        bVL = structure.calcParams[ilay].basisVecL
        c   = Vector{Float64}(redElastCell[ilay])
        normFac = LayerNormFac(structure, ilay)
        r = (r0 + 1):(r0 + bVL)
        # sum_c dK/dc * c as one (complex x real) matrix-vector product; one layer: straight
        # into the output (the block is the whole matrix), no temporaries. `cc`: `A * c`
        # converts `c` to complex first (zgemv); mul! with the real `c` takes another BLAS
        # route and differs in the last bit
        if n == 1
            cc = ComplexF64.(c)
            mul!(vec(gMatAll),    reshape(dKdC_cell[ilay],    bVL * bVL, 15), cc)
            mul!(vec(gMatAll_dD), reshape(d2KdCdD_cell[ilay], bVL * bVL, 15), cc)
            gMatAll    .*= normFac
            gMatAll_dD .*= normFac
        else
            gLay  = reshape(reshape(dKdC_cell[ilay],    bVL * bVL, 15) * c, bVL, bVL)
            dGLay = reshape(reshape(d2KdCdD_cell[ilay], bVL * bVL, 15) * c, bVL, bVL)
            @views gMatAll[r, r]    .= gLay  .* normFac
            @views gMatAll_dD[r, r] .= dGLay .* normFac
        end
        @inbounds for k in r
            nMatAll[k, k] = normFac
        end
        r0 += bVL
    end
    return GammaMatrix(gMatAll, gMatAll_dD, nMatAll)
end

"Error `GammaMatrix:<method>:badPrecompute` unless `c` has one entry per layer."
function _check_cell_length(c, n::Integer, methodName::AbstractString, argName::AbstractString)
    (c isa AbstractVector && length(c) == n) ||
        throw(ArgumentError("GammaMatrix:$(methodName):badPrecompute -- $(argName) must be a " *
                            "cell with one entry per layer ($(n)); got " *
                            "$(c isa AbstractArray ? length(c) : 1)."))
    return nothing
end

"""
    calc(structure::Structure, ::LambBoundary, angleDeg::Real = 0.0;
         gammaByDepthCell = nothing) -> GammaMatrix

Assemble K, dK/dD and M for the frequency-sweep (free-plate / Lamb) path.

`gammaByDepthCell` -- from [`elastByDepth`](@ref); computed if omitted.

Port of `guliwerr/solvers/ritz/@GammaMatrix/GammaMatrix.m` (`GammaMatrix.calcLamb`).
"""
function calc(structure::Structure, ::LambBoundary, angleDeg::Real = 0.0;
              gammaByDepthCell = nothing)
    gByDCell = gammaByDepthCell === nothing ? elastByDepth(structure, angleDeg) : gammaByDepthCell

    gMatAll    = zeros(ComplexF64, 0, 0)
    gMatAll_dD = zeros(ComplexF64, 0, 0)
    nMatAll    = zeros(Float64, 0, 0)

    n = nlay(structure)
    for ilay in 1:n
        gByD = gByDCell[ilay]
        gamma_woRho, dGamma_dD_woRho = ContractDepth(
            gByD[:, :, 1], gByD[:, :, 2], gByD[:, :, 3], structure.calcParams[ilay].relDepth)

        gamma_lay     = gamma_woRho     ./ structure.mater[ilay].rho
        dGamma_dD_lay = dGamma_dD_woRho ./ structure.mater[ilay].rho

        bVL     = structure.calcParams[ilay].basisVecL
        normFac = LayerNormFac(structure, ilay)

        gMatAll    = cat(gMatAll,    gamma_lay     .* normFac; dims = (1, 2))
        gMatAll_dD = cat(gMatAll_dD, dGamma_dD_lay .* normFac; dims = (1, 2))
        nMatAll    = cat(nMatAll,    Matrix{Float64}(I, bVL, bVL) .* normFac; dims = (1, 2))
    end

    return GammaMatrix(gMatAll, gMatAll_dD, nMatAll)
end

# -----------------------------------------------------------------------------------------
# dummyForLayer.m -- inverse-problem support
# -----------------------------------------------------------------------------------------

"""
    dummyForLayer(structure::Structure, noLay = nothing, dKdC_cell = nothing,
                   angleDeg::Real = 0.0) -> (dG_str_dummy, dG_lay)

Inverse support: dK/dc with zeros except for layer `noLay` (sensitivity-matrix input).

`noLay` -- layer index; `nothing` -> top layer (`nlay(structure)`). `dKdC_cell` -- from
[`derivsByDepth`](@ref); computed if omitted.

`dG_str_dummy` -- `strBasisVecL x strBasisVecL x nIndep`, zero except the `noLay` block.
`dG_lay` -- `bVL_lay x bVL_lay x nIndep`, dK/dc for layer `noLay` alone.

Port of `guliwerr/solvers/ritz/@GammaMatrix/dummyForLayer.m`.
"""
function dummyForLayer(structure::Structure, noLay::Union{Integer,Nothing} = nothing,
                        dKdC_cell = nothing, angleDeg::Real = 0.0)
    nl = noLay === nothing ? nlay(structure) : Int(noLay)
    if dKdC_cell === nothing
        dKdC_cell, _ = derivsByDepth(structure)
    end

    nIndep  = structure.mater[nl].elConsts.nIndepConsts
    bVL_lay = structure.calcParams[nl].basisVecL
    rotMat  = LayerRotMat(structure, nl, angleDeg)
    indep2RedMat = ElastRotReductMat(rotMat, structure.mater[nl].elConsts.symmClassNo)
    normFac = LayerNormFac(structure, nl)

    # dK/dc per independent constant for layer nl
    dG_lay = zeros(ComplexF64, bVL_lay, bVL_lay, nIndep)
    for iI in 1:nIndep, iR in 1:15
        @views dG_lay[:, :, iI] .+= dKdC_cell[nl][:, :, iR] .* indep2RedMat[iR, iI]
    end

    # build full-structure dummy: zeros except for the noLay block
    bVL_below = sum(structure.calcParams[i].basisVecL for i in 1:(nl - 1); init = 0)
    strBVL = str_basis_vec_l(structure)
    dG_str_dummy = zeros(ComplexF64, strBVL, strBVL, nIndep)
    rows = (bVL_below + 1):(bVL_below + bVL_lay)
    dG_str_dummy[rows, rows, :] .= dG_lay .* normFac

    return dG_str_dummy, dG_lay
end
