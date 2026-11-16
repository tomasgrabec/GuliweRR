# =========================================================================================
# SpotGreen.jl -- port of guliwerr/forwardCalcs/SpotGreen.m: the half-space
# Green's function as a finite measurement sees it -- the coherent sum of the complex Green's
# function over the wavevectors a ProbeKernel couples to. A picture tool
# (calcSpotGreenMap / calcSpotGreen); it never changes a classification.
# EXPERIMENTAL (as in MATLAB): tried on Fe3Al(316) only, not validated against measurements;
# the interface may change. See the MATLAB header for the formulae.
#
# The kernel enters as its SAMPLE (MATLAB `kernel.sample(dirStep, relKStep)`: offsets,
# dirIndex, delta, kIndex, w, dirStep, relKStep) plus the component name.
# `compute(structureBase, angles, kernel, opts)` (the MATLAB signature) is a thin wrapper:
# `q = sample(kernel; ...)`, then `computeSampled`; `options(kernel; ...)` is the matching
# method taking the kernel. The accumulation (`accumulate`) is separated from the direction
# solves so it can be checked on MATLAB's own solves. The sweep-level
# entry points are `calcSpotGreen` (CalcSingleHS) and `calcSpotGreenMap` (CalcAngSweepHS).
#
# Efficiency and precision (as MATLAB): the direction solves are cheap;
# the spot sum dominates. The Keys-cubic resampling is a 4-tap stencil per output point (no
# sparse matrix, summed in double); every direction's resampled block T is formed once, and
# the weight columns of every (angle, node) pair reading it -- plus the two quadrature
# checks' columns -- go through ONE BLAS matrix product (multithreaded by BLAS itself), in
# ComplexF32: |X| is a picture stored as Float32 anyway and moves by ~1e-6 of its peak
# against a double sum. Columns are then added in the fixed order, so a single angle equals
# the same angle in any sweep (bit for bit as long as BLAS rounds a product column
# independently of the others, as Accelerate does).
# =========================================================================================

module SpotGreen

using LinearAlgebra: mul!
using ...CoreTypes: Structure, nlay, ProbeKernel, sample
using ...Solvers: matrix
using ..ForwardCalcs: rotateStructure, forEachElement, resolveParallel, scoreGreenSolve,
                      standardElastTensor, mround
using ..SawSpectrum: SawSpectrum, latticeAxis
using ..GreenScore: GreenScore

const VEL_STEP       = 0.004
const END_FACTOR     = 1.05
const CONV_WARN      = 0.05
const ETA_STEPS      = 4
const MAX_DIRECTIONS = 20000

"""
    Options

`SpotGreen.options` result: `dirStep` [deg], `relKStep`, `eta`, `velStep` [km/s].
"""
struct Options
    dirStep::Float64
    relKStep::Float64
    eta::Float64
    velStep::Float64
end

"""
    options(kernelDirStep, kernelRelKStep; DirStep, RelKStep, Eta, VelStep) -> Options

Defaults: DirStep = the kernel's dirStep; RelKStep = min(kernel relKStep, Eta / ETA_STEPS);
Eta = SawSpectrum.ETA; VelStep = VEL_STEP (<= GreenScore.MAX_STEP). An explicit RelKStep is
used as given. (MATLAB passes the kernel; its two steps are what is read.)
"""
function options(kernelDirStep::Real, kernelRelKStep::Real; DirStep = kernelDirStep,
                 RelKStep = nothing, Eta = SawSpectrum.ETA, VelStep = VEL_STEP)
    for (nm, x) in (("DirStep", DirStep), ("Eta", Eta), ("VelStep", VelStep))
        (isfinite(x) && x > 0) || throw(ArgumentError("SpotGreen.options: $nm must be positive and finite"))
    end
    VelStep <= GreenScore.MAX_STEP || throw(ArgumentError("SpotGreen.options: VelStep must be <= $(GreenScore.MAX_STEP)"))
    rk = RelKStep === nothing ? min(kernelRelKStep, Eta / ETA_STEPS) : Float64(RelKStep)
    (isfinite(rk) && rk > 0) || throw(ArgumentError("SpotGreen.options: RelKStep must be positive and finite"))
    return Options(DirStep, rk, Eta, VelStep)
end

"""
    options(kernel::ProbeKernel; DirStep, RelKStep, Eta, VelStep) -> Options

MATLAB `SpotGreen.options(kernel, Name-Value...)`: the kernel's `dirStep` and `relKStep` are
the defaults (see the two-number method).
"""
options(kernel::ProbeKernel; kw...) = options(kernel.dirStep, kernel.relKStep; kw...)

"""
    componentIndices(component) -> (I, J)

`"GIJ"` -> (I, J).
"""
function componentIndices(component::AbstractString)
    (length(component) == 3 && component[1] == 'G' && component[2] in '1':'3' &&
     component[3] in '1':'3') || throw(ArgumentError("SpotGreen: component must be 'GIJ', I, J in 1:3; got '$component'"))
    return (component[2] - '0', component[3] - '0')
end

"""
    velocityEnd(structureBase) -> Float64

END_FACTOR x the largest bulk longitudinal velocity over in-plane directions (0:0.5:179.5
deg): an upper bound of every direction's longLimVel.
"""
function velocityEnd(structureBase::Structure)
    mat = structureBase.mater[nlay(structureBase)]
    C = standardElastTensor(mat)
    vMax = 0.0
    for k in 0:359
        phi = k * 0.5
        ev = first(matrix(C, [cosd(phi), sind(phi), 0.0], mat.rho))
        vMax = max(vMax, sqrt(maximum(real.(ev))))
    end
    return END_FACTOR * vMax
end

"""
    solveAxis(vOutEnd, delta, velStep) -> Vector{Float64}

The structure-wide solve axis VEL_MIN + n velStep covering every v / (1 + delta).
"""
function solveAxis(vOutEnd::Real, delta::AbstractVector{<:Real}, velStep::Real)
    v0 = SawSpectrum.VEL_MIN
    nLo = floor(Int, (v0 / (1 + maximum(delta)) - v0) / velStep) - 2
    nHi = ceil(Int, (vOutEnd / (1 + minimum(delta)) - v0) / velStep) + 3
    return latticeAxis(nLo, nHi, velStep)
end

"""
    Resampler

Keys-cubic (a = -1/2) resampling of a solve-axis spectrum at vOut / (1 + delta_m): for row
(m-1) nOut + i, the solve indices `i0[row] - 1 .. i0[row] + 2` (1-based start `col1`) and the
four weights `c[:, row]` (exact zeros skipped, as the sparse matrix of MATLAB stores none).
"""
struct Resampler
    col1::Vector{Int}
    c::Matrix{Float64}
    nOut::Int
    nD::Int
end

"""
    resampler(vOut, delta, vS1, velStep, nS) -> Resampler

MATLAB `SpotGreen.resampleMatrix` as a 4-tap stencil.
"""
function resampler(vOut::AbstractVector{<:Real}, delta::AbstractVector{<:Real}, vS1::Real,
                   velStep::Real, nS::Integer)
    nOut = length(vOut); nD = length(delta)
    col1 = Vector{Int}(undef, nOut * nD); c = Matrix{Float64}(undef, 4, nOut * nD)
    for m in 1:nD, i in 1:nOut
        row = (m - 1) * nOut + i
        t = (vOut[i] / (1 + delta[m]) - vS1) / velStep
        tr = mround(t)
        abs(t - tr) < 1e-9 && (t = tr)
        i0 = floor(Int, t); f = t - i0
        c[1, row] = -0.5 * f + f^2 - 0.5 * f^3
        c[2, row] = 1 - 2.5 * f^2 + 1.5 * f^3
        c[3, row] = 0.5 * f + 2 * f^2 - 1.5 * f^3
        c[4, row] = -0.5 * f^2 + 0.5 * f^3
        col1[row] = i0                       # 0-based i0 - 1 + 1 -> 1-based first column
        (i0 < 1 || i0 + 3 > nS) && throw(ArgumentError(
            "SpotGreen:resampleMatrix:outOfAxis -- the solve axis does not cover the resampling."))
    end
    return Resampler(col1, c, nOut, nD)
end

"Apply the resampler to one solve-axis spectrum R, writing T ([nOut x nD]) in place."
function resample!(T::AbstractMatrix{<:Complex}, rs::Resampler, R::AbstractVector{<:Number})
    @inbounds for row in eachindex(rs.col1)
        j = rs.col1[row]
        acc = zero(ComplexF64)
        for k in 1:4
            ck = rs.c[k, row]
            ck != 0 && (acc += ck * R[j+k-1])
        end
        T[row] = acc                         # double sum, stored in T's precision
    end
    return T
end

"""
    accumulate(angles, q, component, opts, vOut, vS, Gp, Gm) -> (absX, info)

The spot sum on given direction solves: `Gp`/`Gm` ([nS x nDir], one-sided G_IJ(+k) /
G_IJ(-k) of the sorted unique folded directions on the solve axis `vS`), `q` a kernel
sample (offsets, dirIndex, delta, kIndex, w, dirStep, relKStep). Fixed accumulation order:
direction, lobe (+k first), angle, node -- a single angle equals the same angle in any sweep.
"""
function accumulate(angles::AbstractVector{<:Real}, q, component::AbstractString, opts::Options,
                    vOut::AbstractVector{<:Real}, vS::AbstractVector{<:Real},
                    Gp::AbstractMatrix{<:Number}, Gm::AbstractMatrix{<:Number})
    nA = length(angles); nO = length(q.offsets); nD = length(q.delta); nOut = length(vOut)
    rs = resampler(vOut, q.delta, vS[1], opts.velStep, length(vS))
    key180 = [Int64(mround((angles[a] + q.offsets[j]) * 1e6)) for a in 1:nA, j in 1:nO]
    key = mod.(key180, Int64(180_000_000))
    flip = [isodd(Int(floor(Float64(k) / 180e6))) for k in key180]
    uKey = sort!(unique(vec(key)))
    dirOf = [searchsortedfirst(uKey, k) for k in key]
    I, J = componentIndices(component)
    s = (-1, -1, 1); sIJ = s[I] * s[J]
    w = q.w
    evenO = [iseven(d) for d in q.dirIndex]
    evenD = [iseven(d) for d in q.kIndex]
    wDir = w .* reshape(evenO, 1, :); wDir ./= sum(abs, wDir)
    wK = w .* evenD; wK ./= sum(abs, wK)
    checkDir = any(!=(0), q.dirIndex); checkK = any(!=(0), q.kIndex)
    # single precision from here on (see the header): the weights once, each T as it is formed
    wS = ComplexF32.(w); wDirS = ComplexF32.(wDir); wKS = ComplexF32.(wK)
    X = zeros(ComplexF32, nOut, nA); XD = zeros(ComplexF32, nOut, nA); XK = zeros(ComplexF32, nOut, nA)
    rows = [(dirOf[a, j], flip[a, j], a, j) for a in 1:nA for j in 1:nO]
    sort!(rows)
    T = Matrix{ComplexF32}(undef, nOut, nD)
    R = Vector{ComplexF64}(undef, length(vS))
    maxGrp = 0; g = 1                                  # largest group: W, P hold 3 columns per row
    while g <= length(rows)
        h = g
        while h < length(rows) && rows[h+1][1] == rows[g][1] && rows[h+1][2] == rows[g][2]; h += 1; end
        maxGrp = max(maxGrp, h - g + 1); g = h + 1
    end
    W = Matrix{ComplexF32}(undef, nD, 3 * maxGrp)
    P = Matrix{ComplexF32}(undef, nOut, 3 * maxGrp)
    g = 1
    while g <= length(rows)
        d, fl = rows[g][1], rows[g][2]
        if fl
            @inbounds for i in eachindex(R); R[i] = sIJ * Gm[i, d]; end
        else
            @inbounds for i in eachindex(R); R[i] = Gp[i, d]; end
        end
        resample!(T, rs, R)
        h = g
        while h < length(rows) && rows[h+1][1] == d && rows[h+1][2] == fl; h += 1; end
        # every (angle, node) pair reading T, plus the checks' columns: one matrix product
        c = 0
        for r in g:h; c += 1; W[:, c] .= view(wS, :, rows[r][4]); end
        if checkDir
            for r in g:h; evenO[rows[r][4]] && (c += 1; W[:, c] .= view(wDirS, :, rows[r][4])); end
        end
        if checkK
            for r in g:h; c += 1; W[:, c] .= view(wKS, :, rows[r][4]); end
        end
        Pc = view(P, :, 1:c)
        mul!(Pc, T, view(W, :, 1:c))
        c = 0
        for r in g:h; c += 1; view(X, :, rows[r][3]) .+= view(Pc, :, c); end
        if checkDir
            for r in g:h; evenO[rows[r][4]] && (c += 1; view(XD, :, rows[r][3]) .+= view(Pc, :, c)); end
        end
        if checkK
            for r in g:h; c += 1; view(XK, :, rows[r][3]) .+= view(Pc, :, c); end
        end
        g = h + 1
    end
    absX = abs.(X)
    peak = Float64.(vec(maximum(absX; dims = 1)))
    dirConv = fill(NaN, nA); relKConv = fill(NaN, nA)
    checkDir && (dirConv = Float64.(vec(maximum(abs.(absX .- abs.(XD)); dims = 1))) ./ peak)
    checkK && (relKConv = Float64.(vec(maximum(abs.(absX .- abs.(XK)); dims = 1))) ./ peak)
    worst = maximum(filter(!isnan, [dirConv; relKConv; 0.0]))
    if worst > CONV_WARN
        @warn "SpotGreen:notConverged -- the quadrature check changes |X| by up to " *
              "$(round(100 * worst)) % of its peak (direction step $(q.dirStep) deg, |k| step " *
              "$(q.relKStep)) -- refine DirStep/RelKStep."
    end
    info = (dirStep = q.dirStep, relKStep = q.relKStep, eta = opts.eta, velStep = opts.velStep,
            nDirections = length(uKey), nNodes = (nD, nO), dirConvergence = dirConv,
            relKConvergence = relKConv)
    return absX, info, uKey
end

"""
    solveDirections(structureBase, dirs, vS, eta, component; useParallel = :auto)
        -> (Gp, Gm)

One Green solve per direction (each alone on its rotated copy, angle 0): one-sided G_IJ(+k)
and G_IJ(-k) on the solve axis.
"""
function solveDirections(structureBase::Structure, dirs::AbstractVector{<:Real},
                         vS::AbstractVector{<:Real}, eta::Real, component::AbstractString;
                         useParallel = :auto)
    I, J = componentIndices(component)
    useP = resolveParallel(Symbol(useParallel), length(dirs))
    sol = forEachElement(Any[x for x in dirs], useP, (x, _) -> begin
        sr = rotateStructure(structureBase, x)
        g = scoreGreenSolve(sr, vS; eta = eta, forceDir = I, displDir = J)
        (Vector{ComplexF64}(g.G), Vector{ComplexF64}(g.G_minus))
    end)
    Gp = reduce(hcat, [s[1] for s in sol]); Gm = reduce(hcat, [s[2] for s in sol])
    return Gp, Gm
end

"Throw `SpotGreen:compute:layered` unless the structure is a half-space."
checkHalfSpace(structureBase::Structure) =
    nlay(structureBase) == 1 || throw(ArgumentError("SpotGreen:compute:layered -- SpotGreen " *
        "is for a half-space (nLay = 1): on a layered stack |k| sets kd."))

"""
    computeSampled(structureBase, angles, q, component, opts; useParallel = :auto)
        -> (absX, vOut, info)

|X| for a set of displayed directions from a kernel sample `q` (MATLAB `SpotGreen.compute`
with `q = kernel.sample(opts.dirStep, opts.relKStep)`). Half-space only.
"""
function computeSampled(structureBase::Structure, angles::AbstractVector{<:Real}, q,
                        component::AbstractString, opts::Options; useParallel = :auto)
    checkHalfSpace(structureBase)
    all(isfinite, angles) || throw(ArgumentError("SpotGreen.compute: angles must be finite"))
    vEnd = velocityEnd(structureBase)
    vOut = GreenScore.gridAxis(vEnd, opts.velStep)
    vS = solveAxis(vOut[end], q.delta, opts.velStep)
    key = [mod(Int64(mround((a + o) * 1e6)), Int64(180_000_000)) for a in angles, o in q.offsets]
    uKey = sort!(unique(vec(key)))
    length(uKey) > MAX_DIRECTIONS && throw(ArgumentError("SpotGreen:compute:tooManyDirections " *
        "-- $(length(uKey)) direction solves needed (limit $MAX_DIRECTIONS)."))
    Gp, Gm = solveDirections(structureBase, Float64.(uKey) ./ 1e6, vS, opts.eta, component;
                             useParallel)
    absX, info, _ = accumulate(angles, q, component, opts, vOut, vS, Gp, Gm)
    return absX, vOut, merge(info, (vEnd = vEnd,))
end

"""
    compute(structureBase, angles, kernel::ProbeKernel, opts::Options; useParallel = :auto)
        -> (absX, vOut, info)

|X| for a set of displayed directions (MATLAB `SpotGreen.compute`): the kernel is sampled on
`opts.dirStep` / `opts.relKStep` (`sample(kernel; ...)`), then [`computeSampled`](@ref) solves
the directions and sums. `absX` is `[nVel x nAngles]` `Float32`, `vOut` the output axis [km/s];
`info` a `NamedTuple` of `kernel`, `dirStep`, `relKStep`, `eta`, `velStep`, `vEnd`,
`nDirections`, `nNodes` (`(nD, nO)`), `dirConvergence` / `relKConvergence` (`[nAngles]`, `NaN`
where the kernel has a single node in that coordinate). Half-space only. A warning
(`SpotGreen:notConverged`) reports a quadrature step too coarse for the branch slopes.
"""
function compute(structureBase::Structure, angles::AbstractVector{<:Real}, kernel::ProbeKernel,
                 opts::Options; useParallel = :auto)
    checkHalfSpace(structureBase)
    q = sample(kernel; dirStep = opts.dirStep, relKStep = opts.relKStep)
    absX, vOut, info = computeSampled(structureBase, angles, q, kernel.component, opts;
                                      useParallel)
    return absX, vOut, merge((kernel = kernel,), info)
end

end # module SpotGreen
