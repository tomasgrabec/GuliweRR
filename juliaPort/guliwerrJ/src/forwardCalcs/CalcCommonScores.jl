# =========================================================================================
# CalcCommonScores.jl -- Julia-side infrastructure the SAW classification building blocks
# share. The MATLAB CalcCommon functions (thresholdValue, surfaceVRef, shearRange) live in
# CalcCommon.jl -- one home per helper; this file keeps only what has no MATLAB CalcCommon
# counterpart.
#
#   MATLAB                                 Julia (this file, module ForwardCalcs)
#   CalcCommon.thresholdValue              thresholdValue(t, high, low, errId) -> (value, scale)
#   CalcCommon.surfaceVRef                 surfaceVRef(structure)
#   CalcCommon.limitingAndBulk .vS/.vTf    shearRange(lim::LimitingAndBulk) -> (vS, vTf)
#   (solver call sites of the scores)      scoreGreenSolve(sr, v; eta, forceDir, displDir,
#                                          computeMinus, fullOutput) -- ONE adapter for every
#                                          Green solve the scores make (see its docstring)
#   SawSpectrum.memo                       ScoreMemo (a small FIFO memo with a lock)
#
# Also two MATLAB-compat micro-helpers the score maths needs bit-for-bit:
#   mround(x)  -- MATLAB round (ties away from zero; Julia's round is ties-to-even)
#   _with(x; kw...) -- copy of an immutable record with some fields replaced (the MATLAB
#                      "pk(j).field = value" staged fill, on immutable candidate structs)
# =========================================================================================

"""
    mround(x) -> Float64

MATLAB's `round`: half away from zero (`round(2.5) = 3`, `round(-2.5) = -3`). Julia's default
`round` rounds half to even, which differs exactly on the `.5` ties MATLAB code hits when it
rounds a ratio like `FIT_SPAN / dv` or a velocity to 1e-9 km/s.
"""
mround(x::Real) = round(x, RoundNearestTiesAway)

"""
    _with(x::T; kw...) -> T

A copy of the immutable record `x` with the named fields replaced. The scores fill their
candidate records in stages (candidates -> probe slope -> score -> lowScore), as the MATLAB
struct arrays are filled field by field; immutable records keep a stored result from being
changed through an alias.
"""
@inline function _with(x::T; kw...) where {T}
    vals = ntuple(Val(fieldcount(T))) do i
        f = fieldname(T, i)
        haskey(kw, f) ? convert(fieldtype(T, i), kw[f]) : getfield(x, i)
    end
    return T(vals...)
end

"""
    standardElastTensor(mat::Material) -> Array{Float64,4}

MATLAB `mat.GetStandardElastTensor(0)`: the material's tensor in its standard frame, no
further rotation (as `limitingAndBulk` builds it).
"""
standardElastTensor(mat::Material) =
    GetRotElastTensor(SetDirMat(mat, GetStandardDirMat(mat.dirMat)), 0)

"""
    scoreGreenSolve(sr::Structure, v; eta, forceDir = 1, displDir = 3,
                    computeMinus = true, fullOutput = false) -> GreenSolver

The one adapter through which every score module (GreenScore probes, GreenScoreLayered's
polarisation solve, GreenSpotScore / SpotGreen direction solves, the references) calls the
Green solver: `GreenSolver.calc(sr, 0, v.', forceDir, displDir, 'Eta', eta, ...)` in MATLAB.
`sr` is ALREADY rotated, so the angle passed is 0 (the double-rotation rule).

Follows the MATLAB semantics (damping `Eta` as an argument, `G_eff` = `(G13_+ - G31_+)/2`
by reciprocity); `eta` is required and is passed straight to the Julia Green solver.
"""
function scoreGreenSolve(sr::Structure, v::AbstractVector{<:Real}; eta::Real,
                         forceDir::Integer = 1, displDir::Integer = 3,
                         computeMinus::Bool = true, fullOutput::Bool = false)
    return calc(GreenSolver, sr, 0.0, collect(Float64, v); forceDir = forceDir,
                displDir = displDir, computeMinus = computeMinus, fullOutput = fullOutput,
                eta = eta)
end

"""
    ScoreMemo{V}(maxSize)

Port of `SawSpectrum.memo`: a per-session FIFO memo (oldest entry dropped beyond `maxSize`),
keyed by content (`isequal`/`hash` of the key tuple). Guarded by a lock, so threads can read
and write it; MATLAB does not memoise on pool workers, the Julia memo is simply shared.
"""
struct ScoreMemo{V}
    keys::Vector{Any}
    vals::Vector{V}
    maxSize::Int
    lock::ReentrantLock
end
ScoreMemo{V}(maxSize::Integer) where {V} = ScoreMemo{V}(Any[], V[], Int(maxSize), ReentrantLock())

function memoGet(m::ScoreMemo{V}, key) where {V}
    lock(m.lock) do
        for i in eachindex(m.keys)
            isequal(m.keys[i], key) && return m.vals[i]
        end
        return nothing
    end::Union{Nothing,V}
end

function memoPut!(m::ScoreMemo{V}, key, val::V) where {V}
    lock(m.lock) do
        push!(m.keys, key); push!(m.vals, val)
        if length(m.keys) > m.maxSize
            popfirst!(m.keys); popfirst!(m.vals)
        end
    end
    return nothing
end

memoClear!(m::ScoreMemo) = lock(() -> (empty!(m.keys); empty!(m.vals); nothing), m.lock)
memoSize(m::ScoreMemo) = lock(() -> length(m.keys), m.lock)

"""
    structureKey(structure) -> Tuple

Content key of a Structure for the reference memos (MATLAB `SawSpectrum.structureKey`): the
physical content only -- every layer's Voigt matrix and symmetry class, density, direction
matrix, and calcParams' relDepth / legOrder / realThickness; no names, no nMaxModes, so the
same stack under another name hits the memo. Arrays compare by content in a Tuple key.
"""
function structureKey(s::Structure)
    return (Tuple((copy(m.elConsts.voigt), m.elConsts.symmClassNo, m.rho, copy(m.dirMat))
                  for m in s.mater),
            Tuple((c.relDepth, c.legOrder, c.realThickness) for c in s.calcParams))
end
