# =========================================================================================
# RitzSensitivity.jl -- exact first derivatives of the lowest Ritz modes of a layered stack
#
# MATLAB counterpart: RitzSolverSAW.ritzSensitivities (same outputs; it computes all
# eigenpairs and keeps the lowest, ~1e-8 on stiff thin-film stacks where the lowest-mode path
# here is ~1e-11). Agrees with finite differences of the solved velocities to 1e-5...1e-8.
#
# Hellmann-Feynman for K a = lambda M a with a M-normalised (a^H M a = 1):
#     d lambda / d p = a^H (dK/dp - lambda dM/dp) a,    v = sqrt(lambda)/(2 pi),
# so dln v / dln p = p * d lambda/dp / (2 lambda). Per parameter:
#   * elastic constants of layer j: K_j = normFac_j * sum_i c_red,i dK/dc_i (dK/dc from
#     derivsByDepth, already / rho_j), c_red = ElastRotReductMat(R(theta), class) * vec, so
#     d lambda / d vec_k = sum_i R[i, k] * q_i,  q_i = normFac_j * a_j^H (dK/dc_i) a_j;
#   * layer depth D_j: the block (K_j, M_j) is normFac_j * (K_j(D_j), I) with normFac_j ~ D_j,
#     so d lambda / dln D_j = a_j^H (K_j + D_j dK_j/dD) a_j - lambda * a_j^H M_j a_j (with
#     normFac, i.e. the blocks of gammaMatAll / dGammaDMatAll / normMatAll). The overall
#     normalisation by D_1 rho_1 scales K and M alike and drops out, so this holds for the
#     substrate too (its D_1 is the truncation depth: a bound mode should not feel it);
#   * density rho_j: normFac_j K_j does not depend on rho_j, M_j ~ rho_j, so
#     dln v / dln rho_j = -1/2 * (kinetic share of layer j);
#   * propagation angle: dc_red/dtheta by a +-0.01 deg central difference of the constants
#     alone, as angleSlope (no extra solve);
#   * |k| (wavelength): ApplyWavelength scales the films' D_j with k and keeps the substrate's,
#     so dln v / dln k = sum over films of dln v / dln D_j, and v_g / v_p = 1 + that.
# Exact for non-degenerate modes; at a crossing the pair's derivatives are those of whichever
# basis the solver returned.
#
# Workspaces (Julia only): the computation lives in ritzSensitivities!, which uses
# RitzWorkspace memory only (RitzWorkspace.jl); ritzSensitivities is a thin wrapper around it
# with a temporary workspace.
# =========================================================================================

export ritzSensitivities, ritzSensitivities!

"""
    ritzSensitivities(structure, angleDeg, precomp = NamedTuple(); nModes, vMax = Inf)
        -> NamedTuple

The lowest `nModes` Ritz modes of a layered `structure` at `angleDeg` (only those slower than
`vMax`, e.g. the substrate's slowest shear limit for bound modes) and the exact logarithmic
derivatives of their phase velocities with respect to every elastic constant, density and
depth of every layer, the propagation angle and `|k|`. One Ritz solve, no finite differences
except `dc/dtheta`; see the file header for the formulae. Experimental.

`precomp` -- the RitzSolverSAW precompute ([`ritzSAWPrecompute`](@ref), e.g. from a
`RitzWavelengthSweep` at this wavelength); missing fields are built. Only the lowest `nModes`
eigenpairs are computed (`solveConstrainedEig(...; lowest = true)`), so keep `nModes` small
(tens at most).

Fields of the result (`nB <= nModes` modes, rows = modes, ascending velocity):
- `vels` [nB], `eigVals` [nB] (`omega^2`, normalised units);
- `eigVecs` -- full-basis eigenvectors `N * a` [strBasisVecL x nB], `a^H M a = 1`, original
  (not real-gauged) basis, as `N * eigVecs_con` of a RitzSolverSAW;
- `dlnv_dlnC` -- one `nB x nIndepConsts` matrix per layer, `dln v / dln c_k` for the layer's
  `ElastConsts.vector` entries;
- `dlnv_dlnRho`, `dlnv_dlnD` -- `nB x nLay` (layer 1 = substrate; its depth is the truncation
  depth, so `dlnv_dlnD[:, 1]` measures truncation sensitivity rather than physics);
- `dlnv_dTheta` [nB] -- per degree;
- `vgOverVp` [nB] -- `1 + sum_films dlnv_dlnD` (= `dln f / dln k`);
- `kinShare` -- `nB x nLay`, kinetic-energy share of each layer.

Checks that hold to round-off: `sum(dlnv_dlnC) = 1/2` and `sum(dlnv_dlnRho) = -1/2` per mode
(`v^2` is homogeneous of degree 1 in the constants and -1 in the densities).

For many calls (a Q sweep, an inversion loop) use [`ritzSensitivities!`](@ref) with a
[`RitzWorkspace`](@ref): the same numbers (bit-identical), without the few MB of scratch
arrays this function allocates per call.
"""
function ritzSensitivities(structure::Structure, angleDeg::Real,
                           precomp::NamedTuple = NamedTuple();
                           nModes::Integer, vMax::Real = Inf)
    nModes >= 1 || throw(ArgumentError("ritzSensitivities: nModes must be >= 1"))
    N = precomp_field(precomp, :nullConstrMat,
                      () -> buildConstraintMatrix(structure.calcParams, SAWBoundary()))
    ws = RitzWorkspace(structure.calcParams, size(N, 2); nModes, precompute = false)
    return ritzSensitivities!(ws, structure, angleDeg, merge(precomp, (; nullConstrMat = N));
                              nModes, vMax)
end

"""
    ritzSensitivities!(ws::RitzWorkspace, structure, angleDeg, precomp = NamedTuple();
                       nModes = ws.nModes, vMax = Inf) -> NamedTuple

[`ritzSensitivities`](@ref) in the memory of `ws` (built with the structure's Legendre orders,
its `nCon` and this `nModes`, else `ArgumentError`). Scratch arrays are reused; the result is
the same NamedTuple of **fresh, small copies** (vels, eigVals, eigVecs `strBVL x nB`, the
derivative arrays, kinShare), never views into `ws`, so it stays valid after the next call.
With `precomp` from [`ritzSAWPrecompute!`](@ref) on the same `ws` (the usual pattern: one
precompute per wavelength, then one call per angle) a warmed-up call allocates only those
outputs (tens of kB). Missing `precomp` fields are built as in `ritzSensitivities`
(allocating).

Not thread-safe: one workspace per task (`forEachElementPooled`).
"""
function ritzSensitivities!(ws::RitzWorkspace, structure::Structure, angleDeg::Real,
                            precomp::NamedTuple = NamedTuple();
                            nModes::Integer = ws.nModes, vMax::Real = Inf)
    isfinite(angleDeg) ||
        throw(ArgumentError("ritzSensitivities: angleDeg must be a real finite scalar"))
    nModes >= 1 || throw(ArgumentError("ritzSensitivities: nModes must be >= 1"))
    isnan(vMax) && throw(ArgumentError("ritzSensitivities: vMax must be a real scalar"))
    nModes == ws.nModes ||
        throw(ArgumentError("ritzSensitivities!: the workspace was built for nModes = " *
                            "$(ws.nModes), got $nModes"))
    ws.hasSolve ||
        throw(ArgumentError("ritzSensitivities!: the workspace has no solve part (solve = false)"))
    N = precomp_field(precomp, :nullConstrMat,
                      () -> buildConstraintMatrix(structure.calcParams, SAWBoundary()))
    dKdC = precomp_field(precomp, :dKdC_cell)
    d2KdCdD = precomp_field(precomp, :d2KdCdD_cell)
    if dKdC === nothing || d2KdCdD === nothing
        dKdC, d2KdCdD = derivsByDepth(structure)
    end
    sp = precomp_field(precomp, :ritzStruct)

    ang    = Float64(angleDeg)
    cps    = structure.calcParams
    nL     = length(cps)
    nL > 1 || throw(ArgumentError("ritzSensitivities: a layered structure (nLay > 1) is expected"))
    _check_workspace(ws, cps, N, "ritzSensitivities!")
    ranges = _layer_ranges(cps)

    # K and dK/dD layer blocks only (the dense block-diagonal GammaMatrix is mostly zeros);
    # each block is the same product as in calc(structure, SAWBoundary(), ...). `M * c` with a
    # real c converts c to complex first (LinearAlgebra's `*`), hence the complex copy cC.
    # The rotated constants are RedElastCellFor's, from the memoised rotation matrices.
    h    = 0.01                              # dc/dtheta step, as angleSlope
    Rang = _rot_reduct!(ws, structure, ang)
    RP   = _rot_reduct!(ws, structure, ang + h)
    RM   = _rot_reduct!(ws, structure, ang - h)
    nfs  = Float64[LayerNormFac(structure, j) for j in 1:nL]
    Kb, KDb, cC = ws.Kb, ws.KDb, ws.cC
    for j in 1:nL
        bVL = cps[j].basisVecL
        copyto!(cC, Rang[j] * structure.mater[j].elConsts.vector)
        mul!(vec(Kb[j]),  reshape(dKdC[j],    bVL * bVL, 15), cC)
        mul!(vec(KDb[j]), reshape(d2KdCdD[j], bVL * bVL, 15), cC)
        Kb[j]  .*= nfs[j]
        KDb[j] .*= nfs[j]
    end
    K_con, gauge = _project_blocks!(ws, N, ranges, cps)
    useCache = sp !== nothing && sp.N === N && sp.M_con !== nothing && sp.normFacs == nfs
    M_con = useCache ? sp.M_con :
            _project_mass(N, cat([Matrix{Float64}(I, length(r), length(r)) .* nfs[j]
                                  for (j, r) in pairs(ranges)]...; dims = (1, 2)), ranges)
    cholU = useCache && hasproperty(sp, :M_conU) ? sp.M_conU : nothing
    λall, Vall = _lowest_eig!(ws, K_con, M_con, Int(nModes), cholU)

    nk   = min(Int(nModes), length(λall))
    vAll = sqrt.(abs.(λall[1:nk])) ./ (2π)
    nB   = min(nk, count(<(vMax), vAll))

    Vk = view(Vall, :, 1:nB)
    A = if gauge === nothing
        Vk
    elseif sp !== nothing && sp.N === N && sp.gaugeBack !== nothing && eltype(Vk) === Float64
        _gauge_back_cached!(ws, sp.gaugeBack, Vk)
    else
        _gauge_back(N, gauge, Matrix(Vk))
    end
    U   = _full_eigvecs!(ws, N, A)
    lam = λall[1:nB]

    kin = view(ws.kin, 1:nB, :); qK = view(ws.qK, 1:nB, :); qD = view(ws.qD, 1:nB, :)
    Q   = view(ws.Q, :, 1:nB, :)             # q_i per layer: normFac * a_j^H dK/dc_i a_j
    for (j, r) in pairs(ranges)
        Uj  = view(ws.Uj[j], :, 1:nB)
        Uj .= view(U, r, :)
        nf  = nfs[j]
        buf = view(ws.UjBuf[j], :, 1:nB)
        mul!(buf, Kb[j], Uj);  _col_quotients!(view(qK, :, j), Uj, buf)
        mul!(buf, KDb[j], Uj); _col_quotients!(view(qD, :, j), Uj, buf)
        for m in 1:nB
            kin[m, j] = nf * real(dot(view(Uj, :, m), view(Uj, :, m)))
        end
        for i in 1:15
            mul!(buf, view(dKdC[j], :, :, i), Uj)
            _col_quotients!(view(Q, i, :, j), Uj, buf, nf)
        end
    end
    mtot = vec(sum(kin, dims = 2))
    den  = 2 .* lam .* mtot

    dlnv_dlnD   = (qK .+ [cps[j].relDepth for j in 1:nL]' .* qD .- lam .* kin) ./ den
    dlnv_dlnRho = -0.5 .* kin ./ mtot
    dlnv_dlnC   = Vector{Matrix{Float64}}(undef, nL)
    dlnv_dTheta = zeros(nB)
    for j in 1:nL
        ec  = structure.mater[j].elConsts
        Rm  = Rang[j]
        dl  = transpose(Rm) * view(Q, :, :, j)            # nIndep x nB: d lambda / d vec_k
        dlnv_dlnC[j] = Matrix(transpose(ec.vector .* dl ./ den'))
        dc  = (RP[j] * ec.vector .- RM[j] * ec.vector) ./ (2h)
        dlnv_dTheta .+= vec(transpose(dc) * view(Q, :, :, j)) ./ den
    end
    vgOverVp = 1 .+ vec(sum(view(dlnv_dlnD, :, 2:nL), dims = 2))

    return (; vels = vAll[1:nB], eigVals = lam, eigVecs = Matrix(U), dlnv_dlnC, dlnv_dlnRho,
            dlnv_dlnD, dlnv_dTheta, vgOverVp, kinShare = kin ./ mtot)
end

"""
    _project_blocks!(ws, N, ranges, calcParams) -> (K_con, gauge)

`N^T K N` from the layer blocks `ws.Kb` of a block-diagonal `K`, in the real gauge when every
block passes it (as `_project_for_solve`: `K_con::Matrix{Float64}`, `gauge = d`), else complex
(`gauge = nothing`). `N`'s row blocks are used as strided views, not copied. `K_con` is a
workspace array.
"""
function _project_blocks!(ws::RitzWorkspace, N::AbstractMatrix{Float64}, ranges, calcParams)
    Kb, Kr = ws.Kb, ws.Kr
    im2 = 0.0; all2 = 0.0
    for (j, Kj) in pairs(Kb)
        L  = calcParams[j].legOrder
        nb = size(Kj, 1)
        Rj = Kr[j]
        @inbounds for jj in 1:nb
            jU3 = jj > 2L
            for ii in 1:nb
                iU3 = ii > 2L
                z = Kj[ii, jj]
                g = iU3 == jU3 ? z : (iU3 ? complex(imag(z), -real(z)) : complex(-imag(z), real(z)))
                Rj[ii, jj] = real(g)
                im2  += abs2(imag(g))
                all2 += abs2(g)
            end
        end
    end
    real_ok = sqrt(im2) <= REAL_GAUGE_TOL * sqrt(all2)
    tmp = ws.tmp
    if real_ok
        out = fill!(ws.K_con, 0.0)
        for (j, r) in pairs(ranges)
            Nb = view(N, r, :); t = view(tmp, 1:length(r), :)
            mul!(t, Kr[j], Nb)
            mul!(out, transpose(Nb), t, 1.0, 1.0)
        end
        return out, ws.gauge
    end
    cb = _complex_buffers!(ws)
    re = fill!(ws.K_con, 0.0); ii = fill!(cb.Kim, 0.0)
    for (j, r) in pairs(ranges)
        Nb = view(N, r, :); t = view(tmp, 1:length(r), :)
        Kr[j] .= real.(Kb[j]); mul!(t, Kr[j], Nb); mul!(re, transpose(Nb), t, 1.0, 1.0)
        Kr[j] .= imag.(Kb[j]); mul!(t, Kr[j], Nb); mul!(ii, transpose(Nb), t, 1.0, 1.0)
    end
    cb.K_con .= complex.(re, ii)
    return cb.K_con, nothing
end

"""
    _gauge_back_cached!(ws, gb, B) -> view into ws

[`_gauge_back_cached`](@ref) (`N^T (d .* (N B)) = (I - Pi) B + im Pi B`) into workspace
arrays; `B[idx, :]` and `Pi_b * B[idx, :]` use `|idx| x nModes` buffers sized on first use.
"""
function _gauge_back_cached!(ws::RitzWorkspace, gb, B::AbstractMatrix{Float64})
    idx = gb.idx
    nI = length(idx); k = size(ws.A, 2); nB = size(B, 2)
    if size(ws.Bidx, 1) != nI
        ws.Bidx = Matrix{Float64}(undef, nI, k)
        ws.PB   = Matrix{Float64}(undef, nI, k)
    end
    Bi = view(ws.Bidx, :, 1:nB)
    @inbounds for j in 1:nB, (kk, i) in pairs(idx)
        Bi[kk, j] = B[i, j]
    end
    PB = view(ws.PB, :, 1:nB)
    mul!(PB, gb.Pib, Bi)
    out = view(ws.A, :, 1:nB)
    out .= B
    @inbounds for j in 1:nB, (kk, i) in pairs(idx)
        out[i, j] = complex(B[i, j] - PB[kk, j], PB[kk, j])
    end
    return out
end

"`N * A` ([`_full_eigvecs`](@ref): real BLAS on the real and imaginary parts) into `ws.U`."
function _full_eigvecs!(ws::RitzWorkspace, N::AbstractMatrix{Float64},
                        A::AbstractMatrix{ComplexF64})
    nB = size(A, 2)
    U  = view(ws.U, :, 1:nB)
    if _is_unconstrained(N)
        U .= A
        return U
    end
    Ar = view(ws.Ar, :, 1:nB); Ai = view(ws.Ai, :, 1:nB)
    Ur = view(ws.Ur, :, 1:nB); Ui = view(ws.Ui, :, 1:nB)
    Ar .= real.(A); Ai .= imag.(A)
    mul!(Ur, N, Ar); mul!(Ui, N, Ai)
    U .= complex.(Ur, Ui)
    return U
end

"`out[m] = scale * Re(U[:, m]^H AU[:, m])` for every column `m`."
function _col_quotients!(out, U::AbstractMatrix{ComplexF64}, AU::AbstractMatrix{ComplexF64},
                         scale::Real = 1.0)
    @inbounds for m in axes(U, 2)
        acc = 0.0
        for i in axes(U, 1)
            acc += real(conj(U[i, m]) * AU[i, m])
        end
        out[m] = scale * acc
    end
    return out
end
