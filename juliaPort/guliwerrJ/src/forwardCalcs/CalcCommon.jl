# =========================================================================================
# CalcCommon.jl -- port of guliwerr/forwardCalcs/CalcCommon.m .
#
# Shared mechanical internals of the forwardCalcs Calc* family: nothing here decides WHICH
# combination of solves to run (each Calc* type's own job) and nothing here does wave
# physics. Plain top-level definitions in `module ForwardCalcs` (MATLAB's static namespace).
#
# One home per helper: every MATLAB `CalcCommon.*` function lives here (including
# `thresholdValue` / `surfaceVRef`); `CalcCommonScores.jl` keeps only the Julia-side
# infrastructure the score modules share (`scoreGreenSolve`, `ScoreMemo`, `structureKey`,
# `mround`, `_with`, `standardElastTensor`).
#
# Not in MATLAB's current CalcCommon, so not here: `validateVelAxisMode` (velAxisMode),
# `ritzVelMinAcrossAngles`, `greenMapsFrom` (|G33| / GFull maps; now `greenG13EffMapFrom`), the `sawRefEigVec_con`
# precompute field (`calcSAWRefEigVec`). Also not ported: `saveCalc` / `loadCalc` /
# `defaultResultPath` / `resultClasses` / `gitCommit` (MATLAB-file I/O),
# `flattenSweeps` / `scatterToSweeps` (the 2D sweeps' flat dispatch),
# `warnRitzPrecomputeBuiltOnWorker` (Julia threads share memory: a sweep's precompute is
# never lost on a worker, so nothing is rebuilt there).
#
# Load order: this file comes first in ForwardCalcs (the score modules import from it);
# functions that use `ModeClassification` or the Calc types refer to them only inside their
# bodies, which Julia resolves at call time.
# =========================================================================================

# ---- small generic helpers --------------------------------------------------------------

"""
    rotateStructure(structure, angleDeg) -> Structure

Every layer's Material rotated in-plane to the propagation angle (CalcParams preserved).
Callers rotate ONCE, in their constructor, and pass angle 0 to every solver afterwards (see
"Rotation convention: pass 0 to solvers after construction", docs/reference/architecture/forwardCalcs.md).
"""
function rotateStructure(structure::Structure, angleDeg::Real)
    structRot = structure
    for iLay in 1:nlay(structure)
        mat0 = structure.mater[iLay]
        mat0 = SetDirMat(mat0, GetRotDirMat(mat0, angleDeg))
        structRot = ModifyLayer(structRot, iLay, mat0)
    end
    return structRot
end

"""
    resolveParallel(useParallelSetting, nItems) -> Bool

`should_use_parallel(setting, nItems) && threads_available()` (MATLAB: the policy plus a
thread pool; Julia's threads are fixed at start-up).
"""
resolveParallel(useParallelSetting, nItems::Integer) =
    should_use_parallel(useParallelSetting, nItems) && threads_available()

"""
    nestedParallel(useParallelSetting, nOuter, nInner) -> (useOuter::Bool, innerSetting::Symbol)

Parallel decision for a 2D sweep whose outer loop runs inner sweeps that each loop again
(wavelengths x angles, frequencies x angles). Julia's tasks nest, so under `:auto` the
policy's threshold is applied to the TOTAL number of points `nOuter * nInner`: when that
reaches `AUTO_THRESHOLD`, both levels go parallel (`innerSetting = :always`), else neither
(`:never`). `:always` / `:never` pass through unchanged. Applying `:auto` per level instead
left e.g. 6 wavelengths x 10 angles (60 points) entirely serial, as neither level alone
reaches 20. Results do not depend on the decision (BLAS is pinned to one thread in
both branches of [`forEachElement`](@ref)).
"""
function nestedParallel(useParallelSetting::Symbol, nOuter::Integer, nInner::Integer)
    useParallelSetting === :auto ||
        return resolveParallel(useParallelSetting, nOuter), useParallelSetting
    useP = resolveParallel(:auto, nOuter * nInner)
    return useP, (useP ? :always : :never)
end

"""
    forEachElement(inCell, useP, stepFcn) -> Vector

`outCell[i] = stepFcn(inCell[i], i)`, over `Threads.@threads` when `useP`, serially
otherwise; BLAS pinned to one thread in BOTH cases (independent small solves per element,
measured serial, load ~100: one Fe3Al(316) angle's calcRitz 0.7-1.1 s
with 8 OpenBLAS threads vs 0.025-0.04 s with one; the per-element matrices (nCon ~150-600)
gain nothing from threaded BLAS. Pinning both branches also makes the serial and threaded
sweeps bitwise identical. Each task writes only its own slot. Under Apple
Accelerate, whose single-thread switch is per OS thread, every element is additionally run
inside `with_task_blas_single` (BlasBackend.jl).
"""
function forEachElement(inCell::AbstractVector, useP::Bool, stepFcn)
    outCell = copy(inCell)
    if useP
        with_blas_threads(1) do
            Threads.@threads :dynamic for iItem in eachindex(outCell)
                outCell[iItem] = with_task_blas_single(() -> stepFcn(outCell[iItem], iItem))
            end
        end
    else
        with_blas_threads(1) do
            for iItem in eachindex(outCell)
                outCell[iItem] = stepFcn(outCell[iItem], iItem)
            end
        end
    end
    return outCell
end

"""
    forEachElementPooled(inCell, useP, makeWorkspace, stepFcn) -> Vector

`out[i] = stepFcn(inCell[i], i, ws)` with a reusable workspace `ws` (e.g. a
[`RitzWorkspace`](@ref)) handed to every element: [`forEachElement`](@ref) for loops whose
elements would otherwise allocate the same scratch arrays again and again. Threaded
(`Threads.@threads :dynamic`) when `useP`, serial otherwise, BLAS pinned to one thread in both
cases exactly as `forEachElement` (so the two branches give bitwise-equal results whenever
`stepFcn`'s result does not depend on which workspace it got -- true for the `!` Ritz
functions, which overwrite everything they read).

Workspaces: `makeWorkspace()` is called once in the serial branch and
`min(Threads.threadpoolsize(), length(inCell))` times, up front, in the threaded one. They sit
in a `Channel` (a pool): each element `take!`s one, works, and `put!`s it back (also on an
error). A workspace is therefore used by one task at a time, whichever OS thread that task
runs on; indexing workspaces by `Threads.threadid()` would not guarantee this, because a task
can migrate between threads at any yield point (`pitfallsJ.md` jNUM1).

`stepFcn` must not let a workspace array escape into its result (copy what it keeps): the
next element on that workspace overwrites it. Unlike `forEachElement` (`copy(inCell)` updated
in place), the result is a new `Vector`, its element type narrowed to the results' common
type. Julia only (MATLAB's `parfor` gives each worker its own copies).
"""
function forEachElementPooled(inCell::AbstractVector, useP::Bool, makeWorkspace, stepFcn)
    Base.require_one_based_indexing(inCell)
    n   = length(inCell)
    out = Vector{Any}(undef, n)
    if useP && n > 0
        nPool = min(Threads.threadpoolsize(), n)
        ws1   = makeWorkspace()
        pool  = Channel{typeof(ws1)}(nPool)
        put!(pool, ws1)
        for _ in 2:nPool
            put!(pool, makeWorkspace())
        end
        with_blas_threads(1) do
            Threads.@threads :dynamic for i in 1:n
                # `local`: a name assigned anywhere else in this function would make the loop
                # body's variable that one, captured and shared by every task (a data race)
                local wsTask = take!(pool)
                try
                    out[i] = with_task_blas_single(() -> stepFcn(inCell[i], i, wsTask))
                finally
                    put!(pool, wsTask)
                end
            end
        end
        close(pool)
    else
        wsSerial = makeWorkspace()
        with_blas_threads(1) do
            for i in 1:n
                out[i] = stepFcn(inCell[i], i, wsSerial)
            end
        end
    end
    return map(identity, out)
end

"""
    orList(items) -> String

`"'a', 'b', or 'c'"` for error messages.
"""
function orList(items)
    quoted = ["'$item'" for item in items]
    isempty(quoted) && return ""
    length(quoted) == 1 && return quoted[1]
    length(quoted) == 2 && return "$(quoted[1]) or $(quoted[2])"
    return join(quoted[1:end-1], ", ") * ", or " * quoted[end]
end

# ---- validation -------------------------------------------------------------------------

"""
    validateCalcMode(mode, errPrefix)

`:full` | `:ritzOnly`; throws `"<errPrefix>:badMode"`.
"""
function validateCalcMode(mode, errPrefix::AbstractString)
    mode in (:full, :ritzOnly) ||
        throw(ArgumentError("$errPrefix:badMode -- mode must be $(orList((:full, :ritzOnly))); got '$mode'."))
    return nothing
end

"""
    validateRetentionLevel(retentionLevel, errPrefix, allowed)

`retentionLevel in allowed`; throws `"<errPrefix>:badRetentionLevel"`. The SAW family allows
`(:full, :saw_only, :none)`, the plate family only `(:full, :none)`.
"""
function validateRetentionLevel(retentionLevel, errPrefix::AbstractString, allowed)
    retentionLevel in allowed ||
        throw(ArgumentError("$errPrefix:badRetentionLevel -- retentionLevel must be " *
            "$(orList(allowed)); got '$retentionLevel'."))
    return nothing
end

"""
    validateConstructorFlags(calcMode, retentionLevel, errPrefix, allowedRetention)

Constructor-time check of the two flags (each skipped if `nothing`).
"""
function validateConstructorFlags(calcMode, retentionLevel, errPrefix::AbstractString, allowedRetention)
    calcMode === nothing || validateCalcMode(calcMode, errPrefix)
    retentionLevel === nothing || validateRetentionLevel(retentionLevel, errPrefix, allowedRetention)
    return nothing
end

"""
    normaliseThreshold(t) -> Union{Nothing,Symbol,Float64}

A sawThreshold in canonical form: `nothing` (the method's default), `:high` / `:low`
(Symbols or Strings, any case), or a Float64. Anything else is returned unchanged for
`validateSawThreshold` to reject.
"""
normaliseThreshold(t::Nothing) = nothing
normaliseThreshold(t::Union{Symbol,AbstractString}) = Symbol(lowercase(String(t)))
normaliseThreshold(t::Real) = Float64(t)
normaliseThreshold(t) = t

"""
    validateSawThreshold(t, errPrefix)

`nothing` | `:high` | `:low` | a finite real; throws `"<errPrefix>:badSawThreshold"`.
"""
function validateSawThreshold(t, errPrefix::AbstractString)
    t === nothing && return nothing
    t = normaliseThreshold(t)
    ok = (t isa Symbol && t in (:high, :low)) || (t isa Float64 && isfinite(t))
    ok || throw(ArgumentError("$errPrefix:badSawThreshold -- sawThreshold must be 'high', " *
        "'low', a finite real number or nothing (the method's default)."))
    return nothing
end

"""
    modeSawClassMethod(mode) -> Symbol

`:ritz` for `:ritzOnly`, `:green` otherwise.
"""
modeSawClassMethod(mode) = mode === :ritzOnly ? :ritz : :green

"""
    classificationMethods(isLayered) -> (ids, names)

`(:green, :refined, :ritz)` / `("Green", "Refined", "Ritz")`; layered without `:refined`.
"""
function classificationMethods(isLayered::Bool)
    isLayered && return ((:green, :ritz), ("Green", "Ritz"))
    return ((:green, :refined, :ritz), ("Green", "Refined", "Ritz"))
end

"""
    validateSawClassMethod(method, isLayered, errPrefix) -> Symbol

Throws `"<errPrefix>:noRefinedLayered"` for `:refined` on a layered stack,
`"<errPrefix>:badSawClassMethod"` for anything else not allowed.
"""
function validateSawClassMethod(method, isLayered::Bool, errPrefix::AbstractString)
    m = method isa Union{Symbol,AbstractString} ? Symbol(method) : method
    isLayered && m === :refined && throw(ArgumentError("$errPrefix:noRefinedLayered -- " *
        "Method 'refined' is not available for layered stacks (no spot score for layered stacks)."))
    ids, _ = classificationMethods(isLayered)
    (m isa Symbol && m in ids) ||
        throw(ArgumentError("$errPrefix:badSawClassMethod -- sawClassMethod must be $(orList(ids))."))
    return m
end

"""
    resolveSawClassMethod(method, threshold, mode, isLayered, errPrefix) -> Symbol

The classification a calculation runs: unset (`nothing`) follows calcMode (`:green` on
`:full`, `:ritz` on `:ritzOnly`). A numeric threshold needs the method set
(`"<errPrefix>:numericThresholdNeedsMethod"`); `:ritzOnly` allows only `:ritz`
(`"<errPrefix>:ritzOnlyNeedsRitz"`).
"""
function resolveSawClassMethod(method, threshold, mode, isLayered::Bool, errPrefix::AbstractString)
    validateSawThreshold(threshold, errPrefix)
    if method === nothing || method === Symbol("") || method == ""
        normaliseThreshold(threshold) isa Float64 && throw(ArgumentError(
            "$errPrefix:numericThresholdNeedsMethod -- A numeric sawThreshold needs " *
            "sawClassMethod set explicitly (each method has its own scale); or use 'high', " *
            "'low' or nothing (the default)."))
        method = modeSawClassMethod(mode)
    end
    m = validateSawClassMethod(method, isLayered, errPrefix)
    mode === :ritzOnly && m !== :ritz && throw(ArgumentError("$errPrefix:ritzOnlyNeedsRitz -- " *
        "calcMode 'ritzOnly' classifies with 'ritz' only; sawClassMethod is '$m' (set it to " *
        "'ritz', or leave it unset)."))
    return m
end

"""
    thresholdValue(t, high, low, errId) -> (value::Float64, scale::Symbol)

Value and scale of a standard: `"high"`/`:high` -> `(high, :high)`, `"low"` -> `(low, :low)`,
a finite number -> `(t, :high)`. The resolver behind every score's `resolveThreshold`.
"""
function thresholdValue(t, high::Real, low::Real, errId::AbstractString)
    if t isa AbstractString || t isa Symbol
        s = lowercase(String(t))
        s == "high" && return (Float64(high), :high)
        s == "low" && return (Float64(low), :low)
        throw(ArgumentError("$errId -- threshold must be 'high', 'low' or numeric; got '$s'"))
    elseif t isa Real
        isfinite(t) || throw(ArgumentError("$errId -- threshold must be finite; got $t"))
        return (Float64(t), :high)
    end
    throw(ArgumentError("$errId -- threshold must be 'high', 'low' or numeric; got $(typeof(t))"))
end

"""
    defaultThreshold(method, isLayered) -> Symbol

`:low` on a layered stack or for `:ritz`, `:high` otherwise.
"""
defaultThreshold(method, isLayered::Bool) = (isLayered || method === :ritz) ? :low : :high

# ---- Ritz precompute --------------------------------------------------------------------

"""
    ritzSAWPrecompute(structure, precomp = NamedTuple(), rawByDepthCell = nothing) -> NamedTuple

Fill any missing field of the angle-independent RitzSolverSAW precompute (`nullConstrMat`,
`dKdC_cell`, `d2KdCdD_cell`, `ritzStruct` = [`ritzStructurePrecompute`](@ref)): they depend
only on calcParams, depths and density, so a sweep builds them once from the unrotated
structure and shares them with every angle. Other fields of `precomp` are kept.
"""
function ritzSAWPrecompute(structure::Structure, precomp::NamedTuple = NamedTuple(),
                           rawByDepthCell = nothing)
    nullConstrMat = precomp_field(precomp, :nullConstrMat,
        () -> buildConstraintMatrix(structure.calcParams, SAWBoundary()))
    dKdC_cell, d2KdCdD_cell =
        if precomp_field(precomp, :dKdC_cell) === nothing ||
           precomp_field(precomp, :d2KdCdD_cell) === nothing
            rawByDepthCell === nothing ? derivsByDepth(structure) :
                                         derivsByDepth(structure, rawByDepthCell)
        else
            precomp_field(precomp, :dKdC_cell), precomp_field(precomp, :d2KdCdD_cell)
        end
    ritzStruct = precomp_field(precomp, :ritzStruct)
    if ritzStruct === nothing || ritzStruct.N !== nullConstrMat
        ritzStruct = ritzStructurePrecompute(structure, nullConstrMat)
    end
    return merge(precomp, (; nullConstrMat, dKdC_cell, d2KdCdD_cell, ritzStruct))
end

"""
    RitzWavelengthSweep

The wavelength-independent part of the RitzSolverSAW precompute of one layered structure,
built by [`ritzWavelengthSweep`](@ref) and turned into a full precompute at each wavelength by
`ritzSAWPrecompute(structureAtWavelength, sweep)`. Fields: `legOrders` (the Legendre order of
every layer, which is all that the other fields depend on), `rawByDepthCell`
([`RawCellFor`](@ref)), `nullConstrMat`, `gaugeBack` (the real-gauge back-map, see
[`ritzStructurePrecompute`](@ref)) and `gram` ([`ritzGramBlocks`](@ref)).
"""
struct RitzWavelengthSweep
    legOrders::Vector{Int}
    rawByDepthCell::Vector{Array{ComplexF64,4}}
    nullConstrMat::Matrix{Float64}
    gaugeBack::Any
    gram::Vector{Matrix{Float64}}
end

"""
    ritzWavelengthSweep(structure) -> RitzWavelengthSweep

Build once what a sweep over wavelengths (`ApplyWavelength`, i.e. over `|k|`) of one
structure shares: `ApplyWavelength` changes only the films' `relDepth`, while the raw
Legendre tensor, the constraint null space `N`, the gauge back-map and the per-layer Gram
blocks `N_b^T N_b` depend on the Legendre orders alone. At each wavelength,
`ritzSAWPrecompute(structureAtWavelength, sweep)` then only contracts the depth factors and
weights the Gram blocks (on the 4-layer, legOrder-50 stack about 10x cheaper than
`ritzSAWPrecompute(structure)`).

Elastic constants and densities may change between the structures handed to
`ritzSAWPrecompute` too (the raw tensor and `N` carry neither); the Legendre orders may not.
"""
function ritzWavelengthSweep(structure::Structure)
    cps = structure.calcParams
    N   = buildConstraintMatrix(cps, SAWBoundary())
    return RitzWavelengthSweep([cp.legOrder for cp in cps], RawCellFor(structure), N,
                               gaugeBackMap(cps, N), ritzGramBlocks(cps, N))
end

"""
    ritzSAWPrecompute(structure, sweep::RitzWavelengthSweep) -> NamedTuple

The full RitzSolverSAW precompute (same fields as `ritzSAWPrecompute(structure)`) of a
structure that differs from the one `sweep` was built from only in its layer depths
(wavelength), elastic constants or densities. `dKdC_cell`/`d2KdCdD_cell` are bit-identical to
the from-scratch build; `ritzStruct.M_con` equals it to round-off (summed from the Gram
blocks). `nullConstrMat` is the sweep's own matrix (shared, not copied), which is what lets
the solver use the cached pieces. `ArgumentError` when the Legendre orders differ.
"""
function ritzSAWPrecompute(structure::Structure, sweep::RitzWavelengthSweep)
    ws = RitzWorkspace(sweep.legOrders, size(sweep.nullConstrMat, 2); nModes = 1, solve = false)
    return ritzSAWPrecompute!(ws, structure, sweep)    # a temporary workspace: arrays owned here
end

"""
    ritzSAWPrecompute!(ws::RitzWorkspace, structure, sweep::RitzWavelengthSweep) -> NamedTuple

`ritzSAWPrecompute(structure, sweep)` written into the precompute part of `ws` (bit-identical;
the allocating method is this one on a temporary workspace). The returned NamedTuple's
`dKdC_cell`, `d2KdCdD_cell` and `ritzStruct.normFacs`/`M_con`/`M_conU` ARE workspace arrays:
they are valid only until the next `ritzSAWPrecompute!` on `ws`. Use them at once (typically:
one precompute per wavelength, then [`ritzSensitivities!`](@ref) on the same `ws` for each
angle), never keep them across wavelengths. `nullConstrMat` and `gaugeBack` are the sweep's
own (shared, not copied). A warmed-up call allocates almost nothing.

`ArgumentError` when the structure's Legendre orders differ from the sweep's, or when `ws` was
built for other Legendre orders or another `nCon`, or without its precompute part.
"""
function ritzSAWPrecompute!(ws::RitzWorkspace, structure::Structure, sweep::RitzWavelengthSweep)
    cps = structure.calcParams
    (length(cps) == length(sweep.legOrders) &&
     all(cps[i].legOrder == sweep.legOrders[i] for i in eachindex(cps))) ||
        throw(ArgumentError("ritzSAWPrecompute:sweepMismatch -- the structure's Legendre " *
                            "orders differ from those the RitzWavelengthSweep was built for."))
    ws.hasPrecompute ||
        throw(ArgumentError("ritzSAWPrecompute!: the workspace has no precompute part " *
                            "(precompute = false)"))
    N = sweep.nullConstrMat
    _check_workspace(ws, cps, N, "ritzSAWPrecompute!")
    for il in eachindex(cps)       # = derivsByDepth(structure, sweep.rawByDepthCell), in place
        _contract_depth_rho!(ws.dKdC[il], ws.d2KdCdD[il], sweep.rawByDepthCell[il],
                             cps[il].relDepth, structure.mater[il].rho)
    end
    layered = length(cps) > 1
    ritzStruct = _structure_from_gram!(ws.normFacs, layered ? ws.M_con : nothing,
                                       layered ? ws.M_conU : nothing, structure, N,
                                       sweep.gaugeBack, sweep.gram)
    return (; nullConstrMat = N, dKdC_cell = ws.dKdC, d2KdCdD_cell = ws.d2KdCdD, ritzStruct)
end

"""
    RitzWorkspace(sweep::RitzWavelengthSweep; nModes, precompute = true, solve = true)

A [`RitzWorkspace`](@ref) for the structures of `sweep` (its Legendre orders and the number of
columns of its `nullConstrMat`): the form to use in a pool factory, e.g.
`forEachElementPooled(Qs, true, () -> RitzWorkspace(sweep; nModes = 8), stepFcn)`.
"""
RitzWorkspace(sweep::RitzWavelengthSweep; kw...) =
    RitzWorkspace(sweep.legOrders, size(sweep.nullConstrMat, 2); kw...)

"""
    completeRitzPrecompute(structure, precomp) -> NamedTuple

`ritzSAWPrecompute` of a possibly partial precompute (free when complete).
"""
completeRitzPrecompute(structure::Structure, precomp) =
    ritzSAWPrecompute(structure, precomp === nothing ? NamedTuple() : precomp)

# ---- Christoffel ------------------------------------------------------------------------

"""
    LimitingAndBulk

`limitingAndBulk` result: `longLimVel`, `longLimPol` [3], `longLimTheta`, `shearLimVels`
(both shear branches concatenated, variable count), `shearLimPols` [3 x n], `shearLimIwt` (2 or
3 per entry), `shearLimTheta` (sagittal root angles), `vS` / `vTf` (slowest / fastest shear
limit whatever the branch: the Rayleigh boundary and the half-space candidate cap),
`bulkVels` [vL, vT1, vT2], `bulkPols` [3 x 3].
"""
struct LimitingAndBulk
    longLimVel::Float64
    longLimPol::Vector{Float64}
    longLimTheta::Float64
    shearLimVels::Vector{Float64}
    shearLimPols::Matrix{Float64}
    shearLimIwt::Vector{Int}
    shearLimTheta::Vector{Float64}
    vS::Float64
    vTf::Float64
    bulkVels::Vector{Float64}
    bulkPols::Matrix{Float64}
end

"""
    limitingAndBulk(mat) -> LimitingAndBulk

Christoffel limiting velocities (LimitingVels.findRoots: longitudinal = lowest root of branch
1, shear = branches 2 and 3 concatenated) and the bulk velocities along [1, 0, 0], for ONE
material of an ALREADY-ROTATED structure (no further rotation). CalcSingleHS passes the
surface layer (the whole medium), CalcSingleLayered the substrate.
"""
function limitingAndBulk(mat::Material)
    rotElTensor = standardElastTensor(mat)
    rho = mat.rho
    rootsByType = findRoots(rotElTensor, rho)
    longRoots = rootsByType[1]
    isempty(longRoots) && throw(ArgumentError("CalcCommon:limitingAndBulk:noLongitudinalRoot -- " *
        "the longitudinal branch scan found zero roots; expected exactly one"))
    longVels = [r.vel for r in longRoots]
    iLong = argmin(longVels)                       # unique(...)(1): the lowest
    shearAll = vcat(rootsByType[2], rootsByType[3])
    shearLimVels = [r.vel for r in shearAll]
    shearLimPols = isempty(shearAll) ? zeros(3, 0) : reduce(hcat, (r.pol for r in shearAll))
    shearLimIwt = vcat(fill(2, length(rootsByType[2])), fill(3, length(rootsByType[3])))
    shearLimTheta = [r.theta for r in shearAll]
    vS = isempty(shearLimVels) ? NaN : minimum(shearLimVels)
    vTf = isempty(shearLimVels) ? NaN : maximum(shearLimVels)
    eigVals, eigVecs, _ = matrix(rotElTensor, [1.0, 0.0, 0.0], rho)
    return LimitingAndBulk(longRoots[iLong].vel, longRoots[iLong].pol, longRoots[iLong].theta,
                           shearLimVels, shearLimPols, shearLimIwt, shearLimTheta, vS, vTf,
                           sqrt.(eigVals), eigVecs)
end

"""
    shearRange(lim::LimitingAndBulk) -> (vS, vTf)
"""
shearRange(lim::LimitingAndBulk) = (lim.vS, lim.vTf)

"""
    surfaceVRef(structure) -> Float64

`sqrt(G_Hill / rho)` of the surface layer (`mater[nLay]`) [km/s]: the score gates' time scale.
"""
function surfaceVRef(structure::Structure)
    top = structure.mater[nlay(structure)]
    _, gHill = HillAverages(top.elConsts.voigt)
    return sqrt(gHill / top.rho)
end

# ---- Green spectra ----------------------------------------------------------------------

function _checkVelAxis(velAxis, fname)
    (!isempty(velAxis) && all(isfinite, velAxis) && all(>(0), velAxis)) ||
        throw(ArgumentError("$fname:velAxis -- velAxis must be a real, finite, positive vector"))
    return nothing
end

"""
    greenScoreSpectrum(structure, velAxis, fname; threaded = true)
        -> (velAxis, G13Eff::Vector{ComplexF64}, absG13Eff::Vector{Float32})

The one Green solve of a 'full' angle: G13_eff = (G13_+ - G31_+)/2 at `SawSpectrum.ETA` on the
ALREADY-ROTATED structure (angle 0). The magnitude is single precision (display); the complex
values stay double for the score.
"""
function greenScoreSpectrum(structure::Structure, velAxis::AbstractVector{<:Real},
                            fname::AbstractString; threaded::Bool = true)
    _checkVelAxis(velAxis, fname)
    va = collect(Float64, velAxis)
    G = green_component(structure, 0.0, va, 1, 3; side = :eff, eta = SawSpectrum.ETA,
                        threaded = threaded)
    return (velAxis = va, G13Eff = G, absG13Eff = Float32.(abs.(G)))
end

"""
    parseGreenComponent(component) -> (I, J, isEff)

`"GIJ"` (one-sided G_IJ(+k)) or `"GIJeff"` (±k symmetrised), I, J in 1..3.
"""
function parseGreenComponent(component)
    m = match(r"^G([123])([123])(eff)?$", String(component))
    m === nothing && throw(ArgumentError("CalcCommon:parseGreenComponent:badComponent -- " *
        "Green component must be 'GIJ' or 'GIJeff' with I, J in 1..3 (e.g. 'G13eff', 'G13', " *
        "'G33'); got '$component'."))
    return parse(Int, m[1]), parse(Int, m[2]), m[3] !== nothing
end

"""
    greenComponentSpectrum(structure, velAxis, component, fname; threaded = true) -> Vector{Float32}

|G| of any one component on demand (never stored), same damping (`SawSpectrum.ETA`) and
solve as `greenScoreSpectrum`.
"""
function greenComponentSpectrum(structure::Structure, velAxis::AbstractVector{<:Real}, component,
                                fname::AbstractString; threaded::Bool = true)
    _checkVelAxis(velAxis, fname)
    I, J, isEff = parseGreenComponent(component)
    G = green_component(structure, 0.0, collect(Float64, velAxis), I, J;
                        side = isEff ? :eff : :plus, eta = SawSpectrum.ETA, threaded = threaded)
    return Float32.(abs.(G))
end

"""
    greenSpectrum(structure, velAxis, fname) -> (velAxis, absGFull, absG33)

One isoWavevec solve at GreenSolver's DEFAULT damping (3e-4) on any axis: the |G33|
reference the isoFreq tests and validations compare against. Not the Calc pipeline.
"""
function greenSpectrum(structure::Structure, velAxis::AbstractVector{<:Real}, fname::AbstractString)
    _checkVelAxis(velAxis, fname)
    va = collect(Float64, velAxis)
    absGFull = abs.(green_gfull(structure, 0.0, va))
    return (velAxis = va, absGFull = absGFull, absG33 = vec(absGFull[3, 3, :]))
end

"""
    buildAdaptiveVelAxis(candidateVels, velRange, denseHalfWidth = 0.03, nDensePerPeak = 150,
                         nSparse = 150) -> Vector{Float64}

Standalone helper (not used by the Calc types, whose axis is uniform): a dense window around
every candidate velocity unioned with a sparse fill, sorted and de-duplicated.
"""
function buildAdaptiveVelAxis(candidateVels, velRange, denseHalfWidth::Real = 0.03,
                              nDensePerPeak::Integer = 150, nSparse::Integer = 150)
    vMin, vMax = Float64(velRange[1]), Float64(velRange[2])
    cv = filter(v -> vMin <= v <= vMax, Float64.(vec(collect(candidateVels))))
    pts = collect(range(vMin, vMax, length = nSparse))
    for v0 in cv
        append!(pts, range(max(vMin, v0 - denseHalfWidth), min(vMax, v0 + denseHalfWidth),
                           length = nDensePerPeak))
    end
    return sort!(unique(pts))
end

# ---- retention, staleness, classification steps -----------------------------------------

"""
    pruneRitzTo(ritzSolver, classification, retentionLevel, errPrefix) -> RitzSolverSAW

SAW-family retention: `:full` keeps `eigVecs_con`, `:saw_only` keeps the columns of the
significant modes (multi-column in general), `:none` clears them.
"""
function pruneRitzTo(ritzSolver::RitzSolverSAW, classification, retentionLevel,
                     errPrefix::AbstractString)
    if retentionLevel === :full
        return ritzSolver
    elseif retentionLevel === :saw_only
        isempty(ritzSolver.eigVecs_con) && return ritzSolver
        return pruneEigVecs(ritzSolver, [m.modeNo for m in classification.significantModes])
    elseif retentionLevel === :none
        return pruneEigVecs(ritzSolver, nothing)
    end
    throw(ArgumentError("$errPrefix:badRetentionLevel -- retentionLevel must be " *
        "$(orList((:full, :saw_only, :none))); got '$retentionLevel'."))
end

"""
    isStale(result) -> Bool

Classified with a `classificationVersion` other than `CLASSIFICATION_VERSION` (a sweep: any
of its angles).
"""
function isStale(cs)
    hasproperty(cs, :cs) && return any(isStale, cs.cs)
    return cs.classification !== nothing && cs.classificationVersion != CLASSIFICATION_VERSION
end

function refuseStale(cs, errPrefix::AbstractString)
    isStale(cs) && throw(ArgumentError("$errPrefix:stale -- This result was classified by an " *
        "older rule set (classification stamp other than $CLASSIFICATION_VERSION) and cannot be " *
        "reclassified: recompute it."))
    return nothing
end

"""
    greenComponentOf(cs, component) -> Vector{Float32}

|G| of any component on a single-angle result's greenVelAxis: `"G13eff"` returns the stored
`greenAbsG13Eff`, any other component is one solve.
"""
function greenComponentOf(cs, component)
    isempty(cs.greenVelAxis) && throw(ArgumentError("$(nameof(typeof(cs))):greenComponent:noGreenAxis -- " *
        "calcGreen() (or calcAll(:full)) must run first -- greenVelAxis is empty."))
    String(component) == "G13eff" && !isempty(cs.greenAbsG13Eff) && return cs.greenAbsG13Eff
    return greenComponentSpectrum(cs.structure, cs.greenVelAxis, component,
                                  "$(nameof(typeof(cs))).greenComponent")
end

"""
    greenComponentMapOf(sweep, component) -> Matrix{Float32}

Any component over a sweep on its greenVelAxis: `"G13eff"` the stored map, any other one
solve per angle (parallel over the angles per `sweep.useParallel`).
"""
function greenComponentMapOf(sweep, component)
    isempty(sweep.greenVelAxis) && throw(ArgumentError("$(nameof(typeof(sweep))):greenComponentMap:notCalculated -- " *
        "calc() must be run with calcMode=:full first."))
    String(component) == "G13eff" && !isempty(sweep.greenAbsG13EffMap) && return sweep.greenAbsG13EffMap
    return greenComponentMapFrom(sweep.cs, component, sweep.useParallel)
end

function greenComponentMapFrom(csCell, component, useParallel)
    useP = resolveParallel(useParallel, length(csCell))
    cols = forEachElement(Any[c for c in csCell], useP, (c, _) -> greenComponentOf(c, component))
    return reduce(hcat, [Vector{Float32}(c) for c in cols])
end

"""
    classifyGreenPath(cs, thr, christoffelResult, greenResult) -> (classification, limitingEntries)

Classify one angle on a Green path: the limiting entries are associated on the first
classification after a calculation (they need the eigenvectors) and carried afterwards.
"""
function classifyGreenPath(cs, thr, christoffelResult, greenResult)
    mc = cs.limitingDone ?
        ModeClassification(threshold = thr, associateLimiting = false, carryEntries = cs.limitingEntries) :
        ModeClassification(threshold = thr)
    classification = calc(mc, cs.ritzSolver, cs.structure, christoffelResult, greenResult)
    limitingEntries = cs.limitingDone ? cs.limitingEntries :
        SignificantMode[m for m in classification.significantModes if m.type === :limiting]
    return classification, limitingEntries
end

"""
    greenStep(cs, greenVelAxis, retentionLevel, greenScoreRef = nothing, greenThreshold = nothing,
              precomp = nothing) -> cs

Phase body of a 'full' sweep for one angle: the Green solve on the common axis ->
existence score -> Ritz candidates -> classify 'green' -> prune. Prune comes AFTER the
classification (it needs the full eigenvector set).
"""
function greenStep(cs, greenVelAxis, retentionLevel, greenScoreRef = nothing,
                   greenThreshold = nothing, precomp = nothing)
    setRunMode(cs, :full)
    calcGreen(cs, greenVelAxis)
    if cs isa CalcSingleLayered
        calcGreenScore(cs)
    else
        calcGreenScore(cs, greenScoreRef)
    end
    calcRitzScore(cs, precomp)
    identifyModes(cs, :green, greenThreshold)
    applyRetention(cs, retentionLevel)
    cs.retentionLevel = retentionLevel
    return cs
end

"""
    ritzOnlyStep(cs, retentionLevel, ritzThreshold, precomp = nothing) -> cs

`calcAll(cs, :ritzOnly, retentionLevel, precomp)` with the run's method ('ritz') and standard.
"""
function ritzOnlyStep(cs, retentionLevel, ritzThreshold, precomp = nothing)
    cs.sawClassMethod = :ritz
    cs.sawThreshold = normaliseThreshold(ritzThreshold)
    return calcAll(cs, :ritzOnly, retentionLevel, precomp)
end

"""
    sawMapsFrom(csCell, sawTypes) -> (sawVelMap, sawModeNoMap, sawTypeMap, significantModesMap)

Per-angle SAW maps from classified single-point objects (ragged per angle). `sawTypes`:
`(:Rayleigh, :pseudoRayleigh)` half-space, `(:SAW,)` layered.
"""
function sawMapsFrom(csCell, sawTypes)
    n = length(csCell)
    sawVelMap = Vector{Vector{Float64}}(undef, n)
    sawModeNoMap = Vector{Vector{Int}}(undef, n)
    sawTypeMap = Vector{Vector{Symbol}}(undef, n)
    significantModesMap = Vector{Vector{SignificantMode}}(undef, n)
    for i in 1:n
        sm = csCell[i].classification.significantModes
        significantModesMap[i] = sm
        s = [m for m in sm if m.type in sawTypes]
        sawVelMap[i] = [m.vel for m in s]
        sawModeNoMap[i] = [m.modeNo for m in s]
        sawTypeMap[i] = [m.type for m in s]
    end
    return sawVelMap, sawModeNoMap, sawTypeMap, significantModesMap
end

"""
    greenG13EffMapFrom(csCell) -> Matrix{Float32}

`[nVel x nAngles]` |G13_eff| map (every angle on one common axis); empty if any angle lacks it.
"""
function greenG13EffMapFrom(csCell)
    any(c -> isempty(c.greenAbsG13Eff), csCell) && return zeros(Float32, 0, 0)
    return reduce(hcat, [c.greenAbsG13Eff for c in csCell])
end
