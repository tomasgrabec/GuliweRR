# =========================================================================================
# ZgvNewton.jl -- port of guliwerr/solvers/ritz/ZgvNewton.m: Newton refinement of
# zero-group-velocity (ZGV) points of a free plate. MATLAB's static class becomes a module,
# `ZgvNewton.solve` / `ZgvNewton.seedVector`, the same call names (as DispersionBranches).
#
# Solves, for real k, real omega and a unit vector u,
#     (k^2 L2 + k L1 + L0 - omega^2 M) u = 0      (the plate's Ritz eigenproblem)
#     u^H (2k L2 + L1) u                 = 0      (group velocity along k is zero)
#     u^H u                              = 1
# by Newton's method on p = [u; k; omega^2] with the bordered Jacobian and the complex
# correction beta of Lu & Su, which keeps k and omega^2 real. The operator is
# `RitzSolverLamb`'s [`kPolynomial`](@ref).
#
# Adapted from GEWtool's utilities/ZGVNewtonBeta.m (v2.1.1) by Daniel A. Kiefer; the
# algorithm was designed by Bor Plestenjak. Differences from the original: the operator sign
# convention (GULIWERR's L1 already carries the i, omega^2 enters with a minus sign), the
# unknowns are scaled by the seed (k/k0, omega^2/omega0^2) so one tolerance fits every plate,
# an ill-conditioned or non-finite Newton step stops the iteration as not converged, a
# solution with a singular Jacobian is not an isolated ZGV, and the eigenvector seed comes
# from the eigenproblem at k0 rather than from an SVD.
# Reference: D. A. Kiefer, B. Plestenjak, H. Gravenkamp, C. Prada, "Computing
# zero-group-velocity points in anisotropic elastic waveguides: Globally and locally
# convergent methods", J. Acoust. Soc. Am. 153(2), 1386-1398 (2023), doi:10.1121/10.0017252.
#
# GEWtool licence (applies to the adapted parts of this file):
# ISC License
#
# Copyright (c) 2022-2026 Daniel A. Kiefer
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
# REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
# AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
# INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
# LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# =========================================================================================

module ZgvNewton

using LinearAlgebra

"Newton steps (GEWtool default)."
const MAX_STEPS = 10
"Residual tolerance, relative to `norm([L0 L1 L2 M])` of the scaled operator (GEWtool default)."
const TOL = 1e-14
"Largest relative anti-Hermitian part accepted in L2, L1, L0, M."
const HERM_TOL = 1e-10
"A Jacobian with `cond(J, 1)^-1` below this stops the iteration (MATLAB `rcond`)."
const RCOND_MIN = 1e-15
"`k` below `KMIN_REL * k0` is a cut-off (k -> 0), not a ZGV point (GEWtool default)."
const KMIN_REL = 1e-4

_asmat(A) = A isa Number ? fill(complex(float(A)), 1, 1) : Matrix{ComplexF64}(A)
_herm(A) = (A + A') / 2

"""
    luRcond(J) -> (F, rc)

LU factorisation of `J` and its reciprocal 1-norm condition estimate, MATLAB `rcond` (LAPACK
`xGECON` on the LU factors), so that one factorisation serves the estimate and the Newton
step; `rc = 0` for a singular or non-finite matrix.
"""
function luRcond(J::AbstractMatrix)
    all(isfinite, J) || return nothing, 0.0
    F = lu(J; check = false)
    issuccess(F) || return F, 0.0
    return F, LinearAlgebra.LAPACK.gecon!('1', copy(F.factors), opnorm(J, 1))
end

"""
    seedVector(L2, L1, L0, M, k0, omega0) -> (u0, omegaSeed)

Seed vector: the eigenvector of the plate's eigenproblem at `k0` whose eigenvalue is nearest
`omega0^2`, and its angular frequency. Port of `ZgvNewton.seedVector`.
"""
function seedVector(L2, L1, L0, M, k0::Real, omega0::Real)
    A  = _herm(k0^2 .* _asmat(L2) .+ k0 .* _asmat(L1) .+ _asmat(L0))
    Ms = real(_herm(_asmat(M)))
    ev = eigen(Hermitian(A), Hermitian(complex.(Ms)))
    lam = real.(ev.values)
    j = argmin(abs.(lam .- omega0^2))
    u0 = ev.vectors[:, j] ./ norm(ev.vectors[:, j])
    return u0, sqrt(max(lam[j], floatmin(Float64)))
end

"""
    solve(L2, L1, L0, M, k0, omega0, u0 = nothing) -> (k, omega, u, isConverged, nSteps)

Newton iteration for one ZGV point from a seed. Port of `ZgvNewton.solve`.

`L2, L1, L0, M` -- `n x n` Hermitian operator matrices (or scalars for `n = 1`):
`(k^2 L2 + k L1 + L0) u = omega^2 M u`. `k0`, `omega0` -- real positive seed. `u0` -- seed
vector; `nothing`: [`seedVector`](@ref)'s eigenvector, and `omega0` is then replaced by its
own frequency.

Returns `k`, `omega` (`NaN` when not converged), the unit vector `u` (`NaN`s), `isConverged`
(residual below `TOL`, finite steps, the Jacobian regular at the solution -- an isolated ZGV --
and `k >= KMIN_REL * k0`) and the number of Newton steps. Throws `ArgumentError`
`"ZgvNewton:solve:nonHermitian"` for a non-Hermitian operator (a lossy plate would need
GEWtool's ZGVNewtonComplex, not ported).
"""
function solve(L2, L1, L0, M, k0::Real, omega0::Real, u0 = nothing)
    (isfinite(k0) && k0 > 0) || throw(ArgumentError("ZgvNewton.solve: k0 must be real, finite, positive"))
    (isfinite(omega0) && omega0 > 0) ||
        throw(ArgumentError("ZgvNewton.solve: omega0 must be real, finite, positive"))
    mats = (_asmat(L2), _asmat(L1), _asmat(L0), _asmat(M))
    for (im, A) in pairs(mats)
        if norm(A - A') > HERM_TOL * max(norm(A), floatmin(Float64))
            throw(ArgumentError("ZgvNewton:solve:nonHermitian -- the operator matrices must " *
                                "be Hermitian (lossless plate); matrix $im is not."))
        end
    end
    if u0 === nothing
        u0, omega0 = seedVector(L2, L1, L0, M, k0, omega0)
    end

    # scaled unknowns kap = k/k0, m = omega^2/omega0^2: A(kap, m) = kap^2 B2 + kap B1 + B0 - m BM
    mu0 = omega0^2
    B2 = (k0^2 / mu0) .* _herm(mats[1])
    B1 = (k0 / mu0) .* _herm(mats[2])
    B0 = (1 / mu0) .* _herm(mats[3])
    BM = _herm(mats[4])
    n  = size(B0, 1)
    relTol = TOL * norm(hcat(B0, B1, B2, BM))

    kap = 1.0
    m   = 1.0
    u   = ComplexF64.(vec(u0)) ./ norm(u0)
    ek  = zeros(ComplexF64, n + 2); ek[n + 1] = 1
    isConverged = false
    nSteps = 0
    while true
        A  = kap^2 .* B2 .+ kap .* B1 .+ B0 .- m .* BM
        Ak = 2kap .* B2 .+ B1
        F  = vcat(A * u, dot(u, Ak * u), dot(u, u) - 1)
        J  = [A              Ak * u                -BM * u;
              2 .* (u' * Ak) 2 * dot(u, B2 * u)    0;
              2 .* u'        0                     0]
        FJ, rc = luRcond(J)
        isRegular = rc >= RCOND_MIN
        if norm(F) < relTol
            # a solution with a singular Jacobian is not an isolated ZGV point
            isConverged = isRegular
            break
        end
        (nSteps >= MAX_STEPS || !isRegular) && break
        nSteps += 1
        kor = -(FJ \ hcat(F, ek))
        q = kor[:, 1]
        s = kor[:, 2]
        beta = -imag(q[n + 1]) / real(s[n + 1])    # Lu & Su: keeps kap real
        dp = (isfinite(beta) && abs(beta) > 1e-10) ? q .+ im * beta .* s : q
        all(isfinite, dp) || break
        u   = u .+ dp[1:n]
        kap = kap + real(dp[n + 1])                 # imaginary parts are round-off
        m   = m + real(dp[n + 2])
        u   = u ./ norm(u)
    end

    k = kap * k0
    omega = sqrt(max(m, 0.0) * mu0)
    if !isConverged || !(kap >= KMIN_REL) || m <= 0
        return NaN, NaN, fill(complex(NaN), n), false, nSteps
    end
    return k, omega, u, true, nSteps
end

end # module ZgvNewton
