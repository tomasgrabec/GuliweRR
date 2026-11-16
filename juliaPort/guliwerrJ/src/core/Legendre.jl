# src/core/Legendre.jl
#
# L2-normalised Legendre polynomials and the two derived quantities the Ritz-Rayleigh
# solver's `CalcParams` builds from them.
#
# Port of a "half hand-built, half MATLAB-builtin" pair:
#
# - `guliwerr/core/@CalcParams/BasisDecomposition2D.m` (called `legendre(n, x, 'norm')`,
#   MATLAB's own builtin) → `basis_decomposition_2d`, driven by `legendre_norm`
#   below instead — no Julia package reproduces MATLAB's `'norm'` scaling directly,
#   and the Bonnet three-term recurrence is twelve lines and needs no dependency.
# - `guliwerr/core/@CalcParams/LegendreMatrixDiagNondiag.m` (already hand-built, closed form,
#   no builtin) → `legendre_matrix_diag_nondiag`, unchanged in spirit.
#
# Every implementation here was checked against MATLAB dumps; each function's docstring gives
# the measured tolerance.
#
# Plain top-level definitions, not `module Legendre` — see MatlabCompat.jl's header for why
# (one module per top-level folder, not per file).

export legendre_norm, basis_decomposition_2d, legendre_matrix_diag_nondiag

# -----------------------------------------------------------------------------------------
# legendre_norm
# -----------------------------------------------------------------------------------------

"""
    legendre_norm(n::Integer, x::Real) -> Float64

L2-normalised Legendre polynomial of degree `n` on `[-1, 1]`: `sqrt((2n+1)/2) * P_n(x)`.

Equivalent to row 1 (the `m = 0` associated Legendre function) of MATLAB's
`legendre(n, x, 'norm')` — confirmed at `x = 1` (`n=0 → 0.7071 = sqrt(1/2)`,
`n=1 → 1.2247 = sqrt(3/2)`, `n=2 → 1.5811 = sqrt(5/2)`) and by comparison with a MATLAB
dump at `x ∈ {-1, -0.97, -0.5, -0.123, 0, 0.3, 0.7, 0.999, 1}` and
`n ∈ {0, 1, 2, 3, 5, 10, 20, 50, 100, 149}`.

Uses the Bonnet three-term recurrence `(k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}`, forward
and numerically stable for Legendre polynomials at every order tested (worst measured
absolute error 1.3e-13 at `n = 100` — the two implementations' own rounding,
not drift). **No Julia package is used**: `ClassicalOrthogonalPolynomials.jl`,
`FastTransforms.jl` and `SpecialFunctions.jl` were all considered and rejected
— none exposes the `'norm'` scaling directly, each drags in a large dependency tree for
twelve lines of recurrence, and this is not the bottleneck
(`basis_decomposition_2d(150, 0.01)` runs in ~134 µs).
"""
function legendre_norm(n::Integer, x::Real)
    p0 = one(float(x))                      # P_0
    n == 0 && return sqrt(0.5) * p0
    p1 = float(x)                           # P_1
    for k in 1:(n - 1)
        p0, p1 = p1, ((2k + 1) * x * p1 - k * p0) / (k + 1)
    end
    return sqrt((2n + 1) / 2) * p1
end

# -----------------------------------------------------------------------------------------
# basis_decomposition_2d
# -----------------------------------------------------------------------------------------

"""
    basis_decomposition_2d(legO::Integer, meshSpacing::Real = 0.02) -> Matrix{Float64}

Sample the L2-normalised Legendre basis functions of degree `0:legO-1` on a uniform mesh of
`[-1, 1]` with spacing `meshSpacing`. Port of `BasisDecomposition2D.m`; row `k` holds degree
`k - 1`, matching the MATLAB layout exactly.

One recurrence pass per mesh point (computed in-line, not via repeated
[`legendre_norm`](@ref) calls) instead of MATLAB's `legO` separate `legendre(...)` calls.
Verified against a MATLAB dump: `max|pLeg_Julia - pLeg_MATLAB| = 3.3e-16` for a small case,
and for `legO = 150`, `meshSpacing = 0.01` at `atol = 1e-14`.

**Keep both mesh spacings**: this function's default is `0.02`
(MATLAB's own `BasisDecomposition2D` default), but [`CalcParams`](@ref
GuliwerrJ.CoreTypes.CalcParams)'s constructor always calls this with `0.01`, and
`RitzSolver.calcDisplacement` special-cases `dscrtSize == 0.01` to reuse `cp.polyLeg` — do
not "simplify" by unifying the two defaults.
"""
function basis_decomposition_2d(legO::Integer, meshSpacing::Real = 0.02)
    ti   = collect(-1.0:meshSpacing:1.0)
    pLeg = Matrix{Float64}(undef, legO, length(ti))
    for (j, x) in pairs(ti)
        pLeg[1, j] = sqrt(0.5)
        legO >= 2 && (pLeg[2, j] = sqrt(1.5) * x)
        p0, p1 = 1.0, x
        for k in 1:(legO - 2)
            p0, p1 = p1, ((2k + 1) * x * p1 - k * p0) / (k + 1)
            pLeg[k + 2, j] = sqrt((2(k + 1) + 1) / 2) * p1
        end
    end
    return pLeg
end

# -----------------------------------------------------------------------------------------
# legendre_matrix_diag_nondiag
# -----------------------------------------------------------------------------------------

"""
    legendre_matrix_diag_nondiag(legO::Integer) -> (d::Matrix{Float64}, n::Matrix{Float64})

Port of `CalcParams.LegendreMatrixDiagNondiag`.

`n[i, j] = sqrt((2i-1)(2j-1))` for `i > j` and `i + j` odd, zero elsewhere (lower-triangular,
checkerboard-sparse — nonzero only where the two Legendre degrees have opposite parity);
`d = n * transpose(n)` (symmetric, nonzero only where `i + j` is even).

`d[i,j] = ∫ P̃_i' P̃_j' dx`, `n[i,j] = ∫ P̃_i' P̃_j dx` for the L2-normalised polynomials
`P̃_k`. `d = n * n'` follows because the `n[i, k]` are exactly the coefficients of `P̃_i'` in
the orthonormal basis `{P̃_k}`, so orthonormality collapses the `d` integral to
`Σ_k n[i,k] n[j,k]`.

**MATLAB's parity factor cancels its own factor of ½** — worth flagging because it is easy
to get wrong: MATLAB
builds `A = sqrt((2i-1)(2j-1))/2` and multiplies by `1 - (-1)^(i+j)`, which is **2** when
`i+j` is odd and 0 when even, so the product is `sqrt((2i-1)(2j-1))` — the /2 and the ×2
cancel, they do not compound. Keeping the `/2` without the parity factor's implicit ×2
gives `max|n_Julia - n_MATLAB| = 6.982`.

Verified against MATLAB dumps: `n` is exact (`max|n_J - n_M| = 0`, same closed form, same
rounding); `d` carries ~7e-15 from the matrix product's summation order (`legO = 150`).
"""
function legendre_matrix_diag_nondiag(legO::Integer)
    n = zeros(Float64, legO, legO)
    for j in 1:legO, i in 1:legO
        if i > j && isodd(i + j)
            n[i, j] = sqrt((2i - 1) * (2j - 1))
        end
    end
    return n * transpose(n), n
end
