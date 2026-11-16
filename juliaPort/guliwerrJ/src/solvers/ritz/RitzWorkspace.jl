# =========================================================================================
# RitzWorkspace.jl -- preallocated arrays for repeated lowest-mode Ritz solves (Julia only)
#
# Why: a Q sweep with Hellmann-Feynman sensitivities (ritzSensitivities on a
# RitzWavelengthSweep precompute) repeats the same few-MB array allocations thousands of
# times: the per-wavelength dK/dc tables, the K and dK/dD layer blocks, the projection, the
# n x n eigenvector matrix that LAPACK.syevr! allocates even for 8 eigenpairs. On a
# four-layer stack (legOrders 40/20/20/8) one sweep allocated ~11 GB and spent 30-40 % of its
# wall time in Julia's stop-the-world garbage collector on 8-12 threads.
#
# A RitzWorkspace holds every one of those arrays, sized once from the Legendre orders, the
# constrained dimension nCon and nModes. The in-place ("!") functions write into it:
#   ritzSAWPrecompute!(ws, structure, sweep)     (ForwardCalcs, CalcCommon.jl) -- tables,
#                                                 M_con and its Cholesky factor;
#   ritzSensitivities!(ws, structure, angle, pre) (RitzSensitivity.jl) -- everything else.
# The allocating ritzSAWPrecompute(structure, sweep) and ritzSensitivities(...) are thin
# wrappers (a temporary workspace + the "!" call), so there is one implementation.
#
# Results are BIT-IDENTICAL to a plain allocating implementation: every BLAS/LAPACK call is
# the same routine on arrays of the same shape and leading dimension (buffers are sized
# exactly, and a shorter result uses a column prefix of a buffer, which has the same memory
# layout as a fresh array), and every elementwise expression is the same; checked with `==`
# against verbatim copies of the allocating functions.
#
# THREADS: a workspace is mutable scratch memory, so it must never be used by two tasks at
# once, and the arrays a "!" call returns live only until the next call on that workspace.
# Use one workspace per task, e.g. through ForwardCalcs.forEachElementPooled (a Channel pool,
# never indexed by threadid(): tasks can migrate between threads).
#
# LOWEST-MODE EIGENSOLVE without hidden allocations: `_syevr_lowest!` calls LAPACK
# dsyevr/zheevr directly (range 'I', abstol = 2*floatmin, exactly as solveConstrainedEig's
# `lowest` path through LinearAlgebra.LAPACK.syevr!), with the work arrays from a workspace
# query done once at construction. LinearAlgebra's wrapper queries and allocates them on
# every call and allocates an n x n eigenvector matrix for any k.
# =========================================================================================

using LinearAlgebra: BlasInt
using LinearAlgebra.LAPACK: chklapackerror, chkuplofinite

export RitzWorkspace

"Work arrays of one LAPACK `?syevr`/`?heevr` (range 'I', `1:k`) on an `n x n` matrix."
struct _SyevrWork{T}
    n::Int
    k::Int
    w::Vector{Float64}           # eigenvalues (LAPACK writes up to n)
    Z::Matrix{T}                 # n x k eigenvectors (ldz = n)
    isuppz::Vector{BlasInt}
    work::Vector{T}
    rwork::Vector{Float64}       # complex only (empty for real)
    iwork::Vector{BlasInt}
end

"Complex-only buffers (non-gauged stacks), allocated on first use."
mutable struct _RitzComplexBuffers
    Kim::Matrix{Float64}         # Im(N^T K N) while projecting
    K_con::Matrix{ComplexF64}
    Kh::Matrix{ComplexF64}
    UmC::Matrix{ComplexF64}      # complex copy of the Cholesky factor (zhegst)
    eig::_SyevrWork{ComplexF64}
    Wr::Matrix{Float64}
    Wi::Matrix{Float64}
    V::Matrix{ComplexF64}
    MVr::Matrix{Float64}
    MVi::Matrix{Float64}
    MV::Matrix{ComplexF64}
end

"""
    RitzWorkspace(legOrders, nCon; nModes, precompute = true, solve = true)
    RitzWorkspace(structure; nModes, precompute = true, solve = true)
    RitzWorkspace(sweep::RitzWavelengthSweep; nModes, precompute = true, solve = true)

Preallocated arrays for repeated Ritz precomputes and lowest-mode solves of one layered
structure family (any wavelength, constants or densities; fixed Legendre orders):
- `precompute` part, filled by [`ritzSAWPrecompute!`](@ref): the per-layer `dKdC`/`d2KdCdD`
  tables, `normFacs`, `M_con` and its upper Cholesky factor `M_conU`;
- `solve` part, used by [`ritzSensitivities!`](@ref): K / dK/dD / real-gauged layer blocks,
  the projection, the symmetrised and reduced matrices, the `?syevr` work arrays and the
  `nCon x nModes` eigenvectors, the back-map and full-basis eigenvectors, the quotient arrays
  (`Q`, `kin`, `qK`, `qD`), and a memo of the angle-dependent constant-rotation matrices
  (`ElastRotReductMat`, keyed by angle). Complex-only buffers (non-gauged stacks) are added on
  first use.

`nCon` is the number of columns of the constraint null space `N` (the sweep form reads it off
`sweep.nullConstrMat`; the structure form builds `N` once, an SVD, so prefer the sweep form in
a pool). `nModes` is the number of eigenpairs every solve computes; using the workspace with
another `nModes`, other Legendre orders or another `nCon` is an `ArgumentError`.

Not thread-safe: one workspace per task (see `forEachElementPooled`), and arrays returned by a
`!` call are valid only until the next call on the same workspace. Julia only (MATLAB has no
equivalent; its `eig` allocates per call).
"""
mutable struct RitzWorkspace
    legOrders::Vector{Int}
    nCon::Int
    nModes::Int
    # ---- precompute part (ritzSAWPrecompute!)
    hasPrecompute::Bool
    dKdC::Vector{Array{ComplexF64,3}}
    d2KdCdD::Vector{Array{ComplexF64,3}}
    normFacs::Vector{Float64}
    M_con::Matrix{Float64}
    M_conU::Matrix{Float64}
    # ---- solve part (ritzSensitivities!)
    hasSolve::Bool
    gauge::Vector{ComplexF64}
    cC::Vector{ComplexF64}                 # 15 rotated constants, as `*` converts them
    Kb::Vector{Matrix{ComplexF64}}
    KDb::Vector{Matrix{ComplexF64}}
    Kr::Vector{Matrix{Float64}}            # real-gauged block (or Re / Im parts)
    tmp::Matrix{Float64}                   # max(bVL) x nCon
    K_con::Matrix{Float64}                 # real-gauge N^T K N (or its real part)
    Kh::Matrix{Float64}
    Ms::Matrix{Float64}
    Uchol::Matrix{Float64}                 # Cholesky factor when none is cached
    eig::_SyevrWork{Float64}
    MV::Matrix{Float64}
    Bidx::Matrix{Float64}                  # |idx| x nModes, sized on first back-map
    PB::Matrix{Float64}
    A::Matrix{ComplexF64}                  # nCon x nModes, original-basis coordinates
    Ar::Matrix{Float64}
    Ai::Matrix{Float64}
    Ur::Matrix{Float64}                    # strBVL x nModes
    Ui::Matrix{Float64}
    U::Matrix{ComplexF64}
    Uj::Vector{Matrix{ComplexF64}}         # bVL_j x nModes per layer
    UjBuf::Vector{Matrix{ComplexF64}}
    Q::Array{Float64,3}                    # 15 x nModes x nLay
    kin::Matrix{Float64}                   # nModes x nLay
    qK::Matrix{Float64}
    qD::Matrix{Float64}
    cplx::Union{Nothing,_RitzComplexBuffers}
    # ---- memo of ElastRotReductMat per angle (see _rot_reduct!)
    rotDirs::Vector{Matrix{Float64}}
    rotCls::Vector{Int}
    rotCache::Dict{Float64,Vector{Matrix{Float64}}}
end

function RitzWorkspace(legOrders::AbstractVector{<:Integer}, nCon::Integer; nModes::Integer,
                       precompute::Bool = true, solve::Bool = true)
    nL = length(legOrders)
    nL >= 1 && all(>=(1), legOrders) ||
        throw(ArgumentError("RitzWorkspace: legOrders must be >= 1, got $legOrders"))
    bVL = 3 .* Int.(legOrders)
    n0  = sum(bVL)
    1 <= nCon <= n0 ||
        throw(ArgumentError("RitzWorkspace: nCon = $nCon outside 1:$n0 (sum of 3*legOrder)"))
    nModes >= 1 || throw(ArgumentError("RitzWorkspace: nModes must be >= 1, got $nModes"))
    m = Int(nCon); k = min(Int(nModes), m)
    cube(b, d) = d ? Array{ComplexF64,3}(undef, b, b, 15) : Array{ComplexF64,3}(undef, 0, 0, 0)
    mat(T, r, c, d) = d ? Matrix{T}(undef, r, c) : Matrix{T}(undef, 0, 0)
    p, s = precompute, solve
    dKdC    = [cube(b, p) for b in bVL]
    d2KdCdD = [cube(b, p) for b in bVL]
    gauge = ComplexF64[]
    if s
        for L in legOrders
            append!(gauge, ones(ComplexF64, 2L)); append!(gauge, fill(ComplexF64(im), L))
        end
    end
    return RitzWorkspace(collect(Int, legOrders), m, Int(nModes),
        p, dKdC, d2KdCdD, zeros(p ? nL : 0), mat(Float64, m, m, p), mat(Float64, m, m, p),
        s, gauge, zeros(ComplexF64, 15),
        [mat(ComplexF64, b, b, s) for b in bVL], [mat(ComplexF64, b, b, s) for b in bVL],
        [mat(Float64, b, b, s) for b in bVL],
        mat(Float64, maximum(bVL), m, s), mat(Float64, m, m, s), mat(Float64, m, m, s),
        mat(Float64, m, m, s), mat(Float64, m, m, s),
        s ? _SyevrWork(Float64, m, k) : _SyevrWork(Float64, 0, 0),
        mat(Float64, m, k, s), Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0),
        mat(ComplexF64, m, k, s), mat(Float64, m, k, s), mat(Float64, m, k, s),
        mat(Float64, n0, k, s), mat(Float64, n0, k, s), mat(ComplexF64, n0, k, s),
        [mat(ComplexF64, b, k, s) for b in bVL], [mat(ComplexF64, b, k, s) for b in bVL],
        s ? zeros(15, k, nL) : zeros(0, 0, 0), zeros(s ? k : 0, nL), zeros(s ? k : 0, nL),
        zeros(s ? k : 0, nL), nothing, Matrix{Float64}[], Int[],
        Dict{Float64,Vector{Matrix{Float64}}}())
end

RitzWorkspace(calcParams::AbstractVector{CalcParams}, nCon::Integer; kw...) =
    RitzWorkspace([cp.legOrder for cp in calcParams], nCon; kw...)

function RitzWorkspace(structure::Structure; kw...)
    N = buildConstraintMatrix(structure.calcParams, SAWBoundary())
    return RitzWorkspace(structure.calcParams, size(N, 2); kw...)
end

function Base.show(io::IO, ws::RitzWorkspace)
    parts = join(filter(!isempty, [ws.hasPrecompute ? "precompute" : "",
                                   ws.hasSolve ? "solve" : ""]), "+")
    print(io, "RitzWorkspace(legOrders = ", ws.legOrders, ", nCon = ", ws.nCon,
          ", nModes = ", ws.nModes, ", ", parts, ")")
end

"`ArgumentError` unless `ws` was sized for these Legendre orders and this `N`."
function _check_workspace(ws::RitzWorkspace, calcParams, N::AbstractMatrix, who::AbstractString)
    (length(calcParams) == length(ws.legOrders) &&
     all(calcParams[i].legOrder == ws.legOrders[i] for i in eachindex(calcParams))) ||
        throw(ArgumentError("$who: the workspace was built for legOrders $(ws.legOrders); " *
                            "the structure has $([cp.legOrder for cp in calcParams])"))
    size(N, 2) == ws.nCon ||
        throw(ArgumentError("$who: the workspace was built for nCon = $(ws.nCon); " *
                            "the constraint basis has $(size(N, 2)) columns"))
    return nothing
end

const _ROT_CACHE_MAX = 256

"""
    _rot_reduct!(ws, structure, angleDeg) -> Vector{Matrix{Float64}}

`ElastRotReductMat(LayerRotMat(structure, j, angleDeg), symmClassNo_j)` for every layer,
memoised in `ws` by angle. These matrices depend only on each layer's `dirMat`, its symmetry
class and the angle -- not on the wavelength, the constants or the densities -- and a Q sweep
revisits the same few angles (e.g. 20 distinct ones among 978 nodes of one sweep),
while one evaluation costs ~20 us and ~180 kB (a 3^4 tensor rotation per independent
constant), 4 per layer and node. The memo is emptied whenever the structure's `dirMat`s or
classes differ from those it was built for, and when it holds `_ROT_CACHE_MAX` angles. The
values are the very matrices `ElastRotReductMat` returns, so results are unchanged bit for
bit (`RedElastCellFor` is `ElastRotReductMat(...) * ec.vector`).
"""
function _rot_reduct!(ws::RitzWorkspace, structure::Structure, ang::Float64)
    mats = structure.mater
    nL   = length(mats)
    valid = length(ws.rotDirs) == nL &&
            all(i -> ws.rotCls[i] == mats[i].elConsts.symmClassNo &&
                     ws.rotDirs[i] == mats[i].dirMat, 1:nL)
    if !valid
        empty!(ws.rotCache)
        ws.rotDirs = [copy(m.dirMat) for m in mats]
        ws.rotCls  = [m.elConsts.symmClassNo for m in mats]
    end
    r = get(ws.rotCache, ang, nothing)
    r === nothing || return r
    length(ws.rotCache) >= _ROT_CACHE_MAX && empty!(ws.rotCache)
    r = [ElastRotReductMat(LayerRotMat(structure, j, ang), mats[j].elConsts.symmClassNo)
         for j in 1:nL]
    ws.rotCache[ang] = r
    return r
end

"Allocate the complex-only buffers on first use."
function _complex_buffers!(ws::RitzWorkspace)
    ws.cplx === nothing || return ws.cplx
    m = ws.nCon; k = size(ws.eig.Z, 2)
    ws.cplx = _RitzComplexBuffers(Matrix{Float64}(undef, m, m), Matrix{ComplexF64}(undef, m, m),
        Matrix{ComplexF64}(undef, m, m), Matrix{ComplexF64}(undef, m, m),
        _SyevrWork(ComplexF64, m, k), Matrix{Float64}(undef, m, k), Matrix{Float64}(undef, m, k),
        Matrix{ComplexF64}(undef, m, k), Matrix{Float64}(undef, m, k),
        Matrix{Float64}(undef, m, k), Matrix{ComplexF64}(undef, m, k))
    return ws.cplx
end

# -----------------------------------------------------------------------------------------
# ?syevr / ?heevr with preallocated work arrays
# -----------------------------------------------------------------------------------------

const _ABSTOL_LOWEST = 2 * floatmin(Float64)    # as solveConstrainedEig(...; lowest = true)

# argument types exactly as LinearAlgebra.LAPACK.syevr! (the same call serves the size query)
function _syevr_call!(A::StridedMatrix{Float64}, n::Int, k::Int, abstol::Float64,
                      m::Base.RefValue{BlasInt}, w::Vector{Float64}, Z::StridedMatrix{Float64},
                      isuppz::Vector{BlasInt}, work::Vector{Float64}, lwork::Integer,
                      ::Vector{Float64}, ::Integer, iwork::Vector{BlasInt}, liwork::Integer,
                      info::Base.RefValue{BlasInt})
    ccall((LinearAlgebra.BLAS.@blasfunc(dsyevr_), LinearAlgebra.libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
           Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{BlasInt},
           Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
           Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
           Ref{BlasInt}, Clong, Clong, Clong),
          'V', 'I', 'U', n, A, max(1, stride(A, 2)), 0.0, 0.0, 1, k, abstol, m,
          w, Z, max(1, n), isuppz, work, lwork, iwork, liwork, info, 1, 1, 1)
    return nothing
end
function _syevr_call!(A::StridedMatrix{ComplexF64}, n::Int, k::Int, abstol::Float64,
                      m::Base.RefValue{BlasInt}, w::Vector{Float64}, Z::StridedMatrix{ComplexF64},
                      isuppz::Vector{BlasInt}, work::Vector{ComplexF64}, lwork::Integer,
                      rwork::Vector{Float64}, lrwork::Integer, iwork::Vector{BlasInt},
                      liwork::Integer, info::Base.RefValue{BlasInt})
    ccall((LinearAlgebra.BLAS.@blasfunc(zheevr_), LinearAlgebra.libblastrampoline), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ref{Float64}, Ref{Float64},
           Ref{BlasInt}, Ref{BlasInt}, Ref{Float64}, Ptr{BlasInt},
           Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
           Ptr{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong, Clong),
          'V', 'I', 'U', n, A, max(1, stride(A, 2)), 0.0, 0.0, 1, k, abstol, m,
          w, Z, max(1, n), isuppz, work, lwork, rwork, lrwork, iwork, liwork, info, 1, 1, 1)
    return nothing
end

"""
    _SyevrWork(T, n, k)

Work arrays for `?syevr`/`?heevr` over index range `1:k` of an `n x n` matrix, sized by the
LAPACK workspace query (`lwork = -1`), once. The optimal `lwork` the query returns is the one
LinearAlgebra's wrapper resizes to on every call, so the blocked reduction (`?sytrd`, block
size from `ilaenv`) and back-transformation run the same code: results are bit-identical.
"""
function _SyevrWork(::Type{T}, n::Int, k::Int) where {T<:Union{Float64,ComplexF64}}
    n == 0 && return _SyevrWork{T}(0, 0, Float64[], Matrix{T}(undef, 0, 0), BlasInt[], T[],
                                   Float64[], BlasInt[])
    A = zeros(T, n, n); w = Vector{Float64}(undef, n); Z = Matrix{T}(undef, n, k)
    isuppz = Vector{BlasInt}(undef, 2n)
    work = Vector{T}(undef, 1); rwork = Vector{Float64}(undef, 1); iwork = Vector{BlasInt}(undef, 1)
    m = Ref{BlasInt}(); info = Ref{BlasInt}()
    _syevr_call!(A, n, k, _ABSTOL_LOWEST, m, w, Z, isuppz, work, -1, rwork, -1, iwork, -1, info)
    chklapackerror(info[])
    lwork = BlasInt(real(work[1]))
    T <: Complex && (rwork = Vector{Float64}(undef, BlasInt(rwork[1])))
    return _SyevrWork{T}(n, k, w, Z, isuppz, Vector{T}(undef, lwork), rwork,
                         Vector{BlasInt}(undef, iwork[1]))
end

"""
    _syevr_lowest!(sw::_SyevrWork, A) -> nFound

Lowest `sw.k` eigenpairs of the Hermitian `A` (upper triangle; destroyed): eigenvalues in
`sw.w[1:nFound]`, eigenvectors in `sw.Z[:, 1:nFound]`. The same LAPACK call with the same
arguments as `LinearAlgebra.LAPACK.syevr!('V', 'I', 'U', A, 0.0, 0.0, 1, k, 2*floatmin)`,
including its finiteness check of the input.
"""
function _syevr_lowest!(sw::_SyevrWork{T}, A::StridedMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(A)
    n == sw.n || throw(DimensionMismatch("_syevr_lowest!: matrix $n x $n, work sized for $(sw.n)"))
    LinearAlgebra.LAPACK.chkstride1(A)
    chkuplofinite(A, 'U')
    m = Ref{BlasInt}(); info = Ref{BlasInt}()
    _syevr_call!(A, n, sw.k, _ABSTOL_LOWEST, m, sw.w, sw.Z, sw.isuppz, sw.work,
                 length(sw.work), sw.rwork, length(sw.rwork), sw.iwork, length(sw.iwork), info)
    chklapackerror(info[])
    return Int(m[])
end

"""
    _lowest_eig!(ws, K_con, M_con, k, cholU) -> (λ, V)

`solveConstrainedEig(K_con, M_con; nev = k, lowest = true, cholM = cholU)` in workspace
memory: the same steps in the same order (symmetrise, Cholesky or the cached factor,
`?sygst`/`?hegst`, `?syevr` over `1:k`, `V = U^-1 W`, M-renormalisation), bit-identical.
`λ` and `V` are views into `ws`. Cases the workspace does not cover -- `M_con` the identity
(standard path), `k >= n` (full solve), an `M_con` that is not positive definite (QZ
fallback) -- go to `solveConstrainedEig` itself, which allocates.
"""
function _lowest_eig!(ws::RitzWorkspace, K_con::Matrix{T}, M_con, k::Int,
                      cholU) where {T<:Union{Float64,ComplexF64}}
    n = size(K_con, 1)
    if M_con === nothing || isempty(M_con) || _near_identity(M_con) || k >= n
        return solveConstrainedEig(K_con, M_con; nev = k, lowest = true, cholM = cholU)
    end
    cb = T <: Complex ? _complex_buffers!(ws) : nothing
    Kh = T <: Complex ? cb.Kh : ws.Kh
    @inbounds for j in 1:n, i in 1:n
        Kh[i, j] = (K_con[i, j] + conj(K_con[j, i])) / 2
    end
    Ms = ws.Ms
    @inbounds for j in 1:n, i in 1:n
        Ms[i, j] = real(M_con[i, j] + M_con[j, i]) / 2
    end
    Um = cholU
    if Um === nothing
        Um = ws.Uchol
        copyto!(Um, Ms)
        try
            cholesky!(Hermitian(Um))
        catch e
            e isa PosDefException || rethrow()
            return solveConstrainedEig(K_con, M_con; nev = k, lowest = true, cholM = nothing)
        end
        @inbounds for j in 1:n, i in (j + 1):n    # as Matrix(F.U): zeros below the diagonal
            Um[i, j] = 0.0
        end
    end
    Ut = UpperTriangular(Um)
    if T <: Real
        _dsygst_upper!(Kh, Um)
        sw = ws.eig
        nF = _syevr_lowest!(sw, Kh)
        V  = nF == sw.k ? sw.Z : view(sw.Z, :, 1:nF)
        ldiv!(Ut, V)
        MV = nF == sw.k ? ws.MV : view(ws.MV, :, 1:nF)
        mul!(MV, Ms, V)
    else
        copyto!(cb.UmC, Um)
        _zhegst_upper!(Kh, cb.UmC)
        sw = cb.eig
        nF = _syevr_lowest!(sw, Kh)
        cols = nF == sw.k ? Colon() : (1:nF)
        W  = view(sw.Z, :, cols)
        Wr = view(cb.Wr, :, cols); Wi = view(cb.Wi, :, cols)
        Wr .= real.(W); Wi .= imag.(W)
        ldiv!(Ut, Wr); ldiv!(Ut, Wi)
        V = view(cb.V, :, cols)
        V .= complex.(Wr, Wi)
        # real.(V) == Wr and imag.(V) == Wi exactly, so these are solveConstrainedEig's products
        MVr = view(cb.MVr, :, cols); MVi = view(cb.MVi, :, cols)
        mul!(MVr, Ms, Wr); mul!(MVi, Ms, Wi)
        MV = view(cb.MV, :, cols)
        MV .= complex.(MVr, MVi)
    end
    # renormalise so a^H M a = 1, exactly as solveConstrainedEig
    @inbounds for j in axes(V, 2)
        s = 0.0
        for i in axes(V, 1)
            s += real(conj(V[i, j]) * MV[i, j])
        end
        sc = 1 / sqrt(s)
        for i in axes(V, 1)
            V[i, j] *= sc
        end
    end
    λ = view(sw.w, 1:nF)
    if !issorted(λ)
        ord = sortperm(λ)
        return λ[ord], V[:, ord]
    end
    return λ, V
end
