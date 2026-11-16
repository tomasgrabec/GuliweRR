# ProbeKernel.jl -- port of guliwerr/core/ProbeKernel.m
#
# The wavevectors a measurement couples to: its k-space weight W(k).
# EXPERIMENTAL (as in MATLAB): tried on Fe3Al(316) only, not validated against measurements;
# the interface may still change. See the MATLAB class header for the physics; in short
#
#     s(omega) = integral d^2k  W(k) G(k, omega),     W(k) = D^(k) S^(k),
#
# W a function of kappa = k / k0 in the frame of the displayed direction (kappa_x along it),
# over the whole plane (both the +k0 and the -k0 lobe). This file only DESCRIBES W; the
# consumer is SpotGreen (forwardCalcs), which calls `sample`.
#
# Julia API (MATLAB static methods -> plain functions with a `probe_` prefix, so that short
# generic names like `point`/`gaussian` are not exported into the user's namespace):
#
#   MATLAB                                          Julia
#   ProbeKernel(f, ang, relK, dS, rS, comp, name, p) ProbeKernel(f, ang, relK, dS, rS; component, name, params)
#   ProbeKernel()                                   ProbeKernel()            (= probe_point())
#   ProbeKernel.point(comp)                         probe_point(comp = "G13")
#   ProbeKernel.gaussian(sR, sA, 'Component', ...)  probe_gaussian(sR, sA; component, dirStep, relKStep)
#   ProbeKernel.tgs(L, src, prb, ...)               probe_tgs(L, src, prb; component, dirStep, relKStep)
#   ProbeKernel.sras(L, nF, 'DutyCycle', ...)       probe_sras(L, nF; dutyCycle, detectorOffset, ...)
#   ProbeKernel.fromRealSpace(S, D, px, L, ...)     probe_from_real_space(S, D, px, L; lobeRadius, ...)
#   k.sample(dirStep, relKStep)                     sample(k; dirStep, relKStep)
#   ProbeKernel.validateComponent(c)                probe_validate_component(c)
#   ProbeKernel.defaultDirStep / defaultRelKStep    probe_default_dir_step / probe_default_relk_step
#   ProbeKernel.<CONSTANT>                          PROBE_<CONSTANT>
#
# The weight function is a SCALAR function `(kx, ky) -> Number`, broadcast over the node
# arrays by `sample` (MATLAB: an elementwise array function). The built-in kernels use
# small callable structs (not closures) so that `ProbeKernel{F}` is concrete and `sample`
# infers; a user-supplied anonymous function works the same way.

using LinearAlgebra: Tridiagonal
using SpecialFunctions: besselj0, besselj1
using QuadGK: gauss

export ProbeKernel, sample, probe_point, probe_gaussian, probe_tgs, probe_sras,
       probe_from_real_space, probe_validate_component, probe_default_dir_step,
       probe_default_relk_step

# ---------------------------------------------------------------------------------------
# Constants (MATLAB `properties (Constant)`)
# ---------------------------------------------------------------------------------------

"Gaussian constructors: nodes out to this many sigma."
const PROBE_TRUNC_SIGMA = 3
"Default quadrature: at least this many nodes per sigma."
const PROBE_STEPS_PER_SIGMA = 4
"Default direction step ceiling [deg]."
const PROBE_MAX_DIR_STEP = 0.5
"Default direction steps are taken from this list (1/step an integer) [deg]."
const PROBE_NICE_DIR_STEPS = (0.5, 0.25, 0.2, 0.125, 0.1, 0.05, 0.025, 0.02, 0.0125, 0.01, 0.005)
"Zeros of J1: the rings of a disc's transform."
const PROBE_JINC_ZEROS = (3.831705970207512, 7.015586669815619, 10.173468135062722,
                          13.323691936314223, 16.470630050877634, 19.615858510468243)
"sras: mask harmonics |n| <= this are summed."
const PROBE_SRAS_HARMONICS = 7
"sras / fromRealSpace: default quadrature step = shortest oscillation / this."
const PROBE_SAMPLES_PER_OSC = 8

# ---------------------------------------------------------------------------------------
# The struct
# ---------------------------------------------------------------------------------------

"""
    ProbeKernel{F,P}

A measurement's k-space weight W(kappa), kappa = k/k0 in the displayed direction's frame
(immutable; MATLAB value class). Fields (MATLAB property names):

- `name::String`            -- short label
- `weightFcn::F`            -- scalar `(kx, ky) -> Number` (complex W), broadcast by `sample`
- `angHalfWidth::Float64`   -- deg: nodes within +-this of 0 and of 180 deg
- `relKHalfWidth::Float64`  -- nodes at |kappa| - 1 within +-this
- `dirStep::Float64`        -- deg, default direction step of the quadrature
- `relKStep::Float64`       -- default relative |k| step of the quadrature
- `component::String`       -- Green component `"GIJ"`, own propagation frame
- `params::P`               -- NamedTuple of the constructor's physical inputs (MATLAB struct)

Construct with `ProbeKernel(weightFcn, angHalfWidth, relKHalfWidth, dirStep, relKStep;
component = "G13", name = "custom", params = (;))` or one of `probe_point`,
`probe_gaussian`, `probe_tgs`, `probe_sras`, `probe_from_real_space`.
"""
struct ProbeKernel{F,P<:NamedTuple}
    name::String
    weightFcn::F
    angHalfWidth::Float64
    relKHalfWidth::Float64
    dirStep::Float64
    relKStep::Float64
    component::String
    params::P
end

_pk_err(id::AbstractString, msg::AbstractString) = throw(ArgumentError(string(id, " -- ", msg)))

function _pk_check_scalar(x, fname, argname; lo = -Inf, hi = Inf, loIncl = true, hiIncl = true,
                          finite = true)
    x isa Real || _pk_err("$fname:$argname", "$argname must be a real scalar.")
    finite && !isfinite(x) && _pk_err("$fname:$argname", "$argname must be finite, got $x.")
    isnan(x) && _pk_err("$fname:$argname", "$argname must not be NaN.")
    okLo = loIncl ? x >= lo : x > lo
    okHi = hiIncl ? x <= hi : x < hi
    (okLo && okHi) || _pk_err("$fname:$argname",
        "$argname must be in $(loIncl ? '[' : '(')$lo, $hi$(hiIncl ? ']' : ')'), got $x.")
    return Float64(x)
end

"""
    ProbeKernel(weightFcn, angHalfWidth, relKHalfWidth, dirStep, relKStep;
                component = "G13", name = "custom", params = (;))

A kernel from any scalar weight function `(kx, ky) -> Number` (kappa = k/k0, whole plane,
include the -k0 lobe if the measurement has one). `angHalfWidth` in deg, [0, 90);
`relKHalfWidth` in [0, 1); `dirStep` (deg) and `relKStep` > 0.
"""
function ProbeKernel(weightFcn, angHalfWidth::Real, relKHalfWidth::Real, dirStep::Real,
                     relKStep::Real; component::AbstractString = "G13",
                     name::AbstractString = "custom", params::NamedTuple = NamedTuple())
    fn = "ProbeKernel"
    a = _pk_check_scalar(angHalfWidth, fn, "angHalfWidth"; lo = 0, hi = 90, hiIncl = false)
    r = _pk_check_scalar(relKHalfWidth, fn, "relKHalfWidth"; lo = 0, hi = 1, hiIncl = false)
    d = _pk_check_scalar(dirStep, fn, "dirStep"; lo = 0, loIncl = false)
    s = _pk_check_scalar(relKStep, fn, "relKStep"; lo = 0, loIncl = false)
    c = probe_validate_component(component)
    return ProbeKernel{typeof(weightFcn),typeof(params)}(String(name), weightFcn, a, r, d, s,
                                                         c, params)
end

"Default kernel = `probe_point()` (MATLAB `ProbeKernel()` with no arguments)."
ProbeKernel() = probe_point()

function Base.show(io::IO, k::ProbeKernel)
    print(io, "ProbeKernel(\"", k.name, "\", ", k.component, ", ang +-", k.angHalfWidth,
          " deg, relK +-", k.relKHalfWidth, ", steps ", k.dirStep, " deg / ", k.relKStep, ")")
end

# ---------------------------------------------------------------------------------------
# sample
# ---------------------------------------------------------------------------------------

"""
    sample(k::ProbeKernel; dirStep = k.dirStep, relKStep = k.relKStep) -> NamedTuple

Quadrature nodes and normalised complex weights of the kernel (MATLAB `k.sample(dirStep,
relKStep)`; `nothing` = the kernel's default, as MATLAB's `[]`). Fields, as MATLAB's `q`:

- `offsets::Vector{Float64}`  [nO] deg -- node directions relative to the displayed one:
                               (-m:m) x dirStep around 0, then the same around 180
- `dirIndex::Vector{Int}`     [nO] integer -m..m of each offset within its lobe
- `delta::Vector{Float64}`    [nD] |k|/k0 - 1 = (-n:n) x relKStep
- `kIndex::Vector{Int}`       [nD] integer -n..n
- `w::Matrix{T}`              [nD x nO], sum(abs.(w)) = 1; T = the weight function's
                               (floated) return type: Float64 for point / gaussian / tgs,
                               ComplexF64 for sras / fromRealSpace
- `dirStep`, `relKStep`       the steps used

Weight = W(kappa) x dirStep x relKStep at kappa = (1 + delta)(cos o, sin o), normalised
by the total |w|; directions or |k| rows below 1e-9 of the largest |w| are dropped.
Column-major layout as in MATLAB (rows = delta, columns = offsets).
"""
function sample(k::ProbeKernel; dirStep::Union{Real,Nothing} = nothing,
                relKStep::Union{Real,Nothing} = nothing)
    ds = dirStep === nothing ? k.dirStep :
         _pk_check_scalar(dirStep, "ProbeKernel.sample", "dirStep"; lo = 0, loIncl = false)
    rs = relKStep === nothing ? k.relKStep :
         _pk_check_scalar(relKStep, "ProbeKernel.sample", "relKStep"; lo = 0, loIncl = false)
    m = floor(Int, k.angHalfWidth / ds + 1e-9)
    n = floor(Int, k.relKHalfWidth / rs + 1e-9)
    nO1 = 2m + 1
    offsets = Vector{Float64}(undef, 2nO1)
    dirIndex = Vector{Int}(undef, 2nO1)
    @inbounds for (j, i) in enumerate(-m:m)
        offsets[j] = i * ds
        offsets[j+nO1] = 180 + i * ds
        dirIndex[j] = i
        dirIndex[j+nO1] = i
    end
    kIndex = collect(-n:n)
    delta = kIndex .* rs
    dA = 1.0
    m > 0 && (dA *= ds)
    n > 0 && (dA *= rs)
    co = cosd.(offsets)
    so = sind.(offsets)
    r = 1 .+ delta
    W = _pk_eval_weights(k.weightFcn, r, co, so, dA)      # [nD x nO]
    all(isfinite, W) || _pk_err("ProbeKernel:sample:badWeight",
        "weightFcn must return finite values (scalar weight, broadcast over the nodes).")
    # column / row maxima of |W| in one pass
    colMax = zeros(Float64, size(W, 2))
    rowMax = zeros(Float64, size(W, 1))
    @inbounds for j in axes(W, 2), i in axes(W, 1)
        a = abs(W[i, j])
        a > colMax[j] && (colMax[j] = a)
        a > rowMax[i] && (rowMax[i] = a)
    end
    amax = maximum(colMax)
    amax == 0 && _pk_err("ProbeKernel:sample:zeroKernel", "The kernel is zero at every node.")
    keepO = colMax .> 1e-9 * amax
    keepD = rowMax .> 1e-9 * amax
    Wk = (all(keepO) && all(keepD)) ? W : W[keepD, keepO]
    tot = sum(abs, Wk)
    Wk ./= tot
    return (offsets = offsets[keepO], dirIndex = dirIndex[keepO], delta = delta[keepD],
            kIndex = kIndex[keepD], w = Wk, dirStep = ds, relKStep = rs)
end

"Node count from which `sample` evaluates the weight function on several threads."
const PROBE_THREAD_MIN_NODES = 4096

# W[i, j] = f(r_i cos o_j, r_i sin o_j) * dA. Each element is independent, so the threaded
# path (columns split over tasks, disjoint output slots) is bitwise identical
# to the serial one; it pays off for the Bessel-heavy sras kernel (~15 J1 per node). The
# weight function must therefore be thread-safe (pure), as every built-in one is.
function _pk_eval_weights(f::F, r::Vector{Float64}, co::Vector{Float64}, so::Vector{Float64},
                          dA::Float64) where {F}
    T = Base.promote_op((x, y) -> f(x, y) * dA, Float64, Float64)
    nD = length(r); nO = length(co)
    if !isconcretetype(T) || nD * nO < PROBE_THREAD_MIN_NODES || Threads.nthreads() == 1
        return f.(r .* co', r .* so') .* dA                   # one fused pass
    end
    W = Matrix{T}(undef, nD, nO)
    Threads.@threads for j in 1:nO
        cj = co[j]; sj = so[j]
        @inbounds for i in 1:nD
            W[i, j] = f(r[i] * cj, r[i] * sj) * dA
        end
    end
    return W
end

# ---------------------------------------------------------------------------------------
# Static helpers
# ---------------------------------------------------------------------------------------

"""
    probe_validate_component(component) -> String

`"GIJ"` with I, J in 1..3 (own propagation frame, one-sided notation); throws
`ArgumentError` ("ProbeKernel:validateComponent:badComponent") otherwise.
"""
function probe_validate_component(component)
    ok = (component isa AbstractString || component isa Symbol) &&
         occursin(r"^G[123][123]$", String(component))
    ok || _pk_err("ProbeKernel:validateComponent:badComponent",
        "component must be 'GIJ' with I, J in 1..3 (own propagation frame; the kernel's " *
        "two lobes make the +-k combination) -- got '$(component)'.")
    return String(component)
end

"The largest PROBE_NICE_DIR_STEPS entry <= min(target, PROBE_MAX_DIR_STEP); the smallest if none."
function _probe_nice_dir_step(target::Real)
    lim = min(PROBE_MAX_DIR_STEP, target) + 1e-12
    for s in PROBE_NICE_DIR_STEPS
        s <= lim && return s
    end
    return PROBE_NICE_DIR_STEPS[end]
end

"""
    probe_default_dir_step(sigmaAngDeg)

The largest `PROBE_NICE_DIR_STEPS` entry <= sigma / PROBE_STEPS_PER_SIGMA (and <=
PROBE_MAX_DIR_STEP); 1 for sigma = 0.
"""
probe_default_dir_step(sigmaAngDeg::Real) =
    sigmaAngDeg == 0 ? 1.0 : _probe_nice_dir_step(sigmaAngDeg / PROBE_STEPS_PER_SIGMA)

"""
    probe_default_relk_step(sigmaRelK)

sigma / PROBE_STEPS_PER_SIGMA; 1 for sigma = 0 (SpotGreen caps it at Eta / ETA_STEPS).
"""
probe_default_relk_step(sigmaRelK::Real) =
    sigmaRelK == 0 ? 1.0 : Float64(sigmaRelK / PROBE_STEPS_PER_SIGMA)

"Default steps from the shortest oscillation step dq (kappa), MATLAB oscSteps."
_probe_osc_steps(dq, dirStep, relKStep) =
    (dirStep === nothing ? _probe_nice_dir_step(rad2deg(dq)) : Float64(dirStep),
     relKStep === nothing ? Float64(dq) : Float64(relKStep))

function _pk_step_opt(x, fname, argname)
    x === nothing && return nothing
    return _pk_check_scalar(x, fname, argname; lo = 0, loIncl = false)
end

# %.Ng-like formatting for kernel names (avoids a Printf dependency; labels only)
function _pk_fmtg(x::Real, p::Int = 6)
    isfinite(x) || return string(x)
    x == 0 && return "0"
    xr = round(Float64(x); sigdigits = p)
    e = floor(Int, log10(abs(xr)))
    if e < -4 || e >= p
        mant = round(xr / 10.0^e; digits = p - 1)
        ms = _pk_strip(string(mant))
        return string(ms, "e", e < 0 ? "-" : "+", lpad(string(abs(e)), 2, '0'))
    end
    d = max(p - 1 - e, 0)
    s = string(round(xr; digits = d))
    occursin('e', s) && (s = _pk_fixed(xr, d))
    return _pk_strip(s)
end
_pk_strip(s) = occursin('.', s) ? rstrip(rstrip(s, '0'), '.') : s
function _pk_fixed(x, d)
    i = round(BigInt, abs(x) * big(10)^d)
    s = lpad(string(i), d + 1, '0')
    s = d > 0 ? string(s[1:end-d], ".", s[end-d+1:end]) : s
    return x < 0 ? "-" * s : s
end

# ---------------------------------------------------------------------------------------
# Weight functions (callable structs: scalar (kx, ky) -> W)
# ---------------------------------------------------------------------------------------

"1 at kappa = (+-1, 0), 0 elsewhere (MATLAB pointWeight)."
struct ProbePointWeight end
(::ProbePointWeight)(kx::Real, ky::Real) =
    Float64(abs(abs(kx) - 1) < 1e-12 && abs(ky) < 1e-12)

"exp(-x^2/(2 s^2)); a delta (1 at x = 0) for s = 0."
@inline _pk_gauss1(x::Real, s::Float64) =
    s == 0 ? Float64(abs(x) < 1e-9) : exp(-x^2 / (2 * s^2))

"Polar Gaussian around both lobes (dphi from the nearer of 0 / 180 deg)."
struct ProbeGaussPolarWeight
    sR::Float64
    sA::Float64
end
function (g::ProbeGaussPolarWeight)(kx::Real, ky::Real)
    r = hypot(kx, ky)
    phi = atand(ky, kx)
    dphi = min(abs(phi), 180 - abs(phi))
    return _pk_gauss1(r - 1, g.sR) * _pk_gauss1(dphi, g.sA)
end

"Cartesian Gaussian clouds at kappa = (+1, 0) and (-1, 0)."
struct ProbeGaussCartesianWeight
    sx::Float64
    sy::Float64
end
(g::ProbeGaussCartesianWeight)(kx::Real, ky::Real) =
    (_pk_gauss1(kx - 1, g.sx) + _pk_gauss1(kx + 1, g.sx)) * _pk_gauss1(ky, g.sy)

"2 J1(x) / x, = 1 at x = 0."
struct ProbeJinc end
(::ProbeJinc)(x::Real) = x == 0 ? 1.0 : 2 * besselj1(x) / x

"Knife-edge slope i k.n (NaN direction: none), force |k| (kFactors); the lobe cut is the caller's."
@inline function _pk_detection_factors(W::ComplexF64, kx, ky, cKn, sKn, knife::Bool,
                                       kFactors::Bool)
    kr = hypot(kx, ky)
    if knife
        kn = kx * cKn + ky * sKn
        kFactors || (kn = kn / max(kr, eps()))
        W = W * Complex(0.0, kn)
    end
    kFactors && (W = W * kr)
    return W
end

@inline _pk_outside_lobes(kx, ky, rho) = min(hypot(kx - 1, ky), hypot(kx + 1, ky)) > rho

"sras's W (MATLAB srasWeight); Rt = k0 R, xdt = k0 x_d, rdt = k0 r_d, lobe(qR)."
struct ProbeSrasWeight{L}
    Rt::Float64
    n::Vector{Float64}
    c::Vector{Float64}
    lobe::L
    xdt1::Float64
    xdt2::Float64
    damp::Float64            # blurt^2 / 2 + rdt^2 / 8
    knife::Bool
    cKn::Float64
    sKn::Float64
    kFactors::Bool
    rho::Float64
end
function (w::ProbeSrasWeight)(kx::Real, ky::Real)
    _pk_outside_lobes(kx, ky, w.rho) && return zero(ComplexF64)
    S = 0.0
    @inbounds for i in eachindex(w.n)
        S += w.c[i] * w.lobe(w.Rt * hypot(kx - w.n[i], ky))
    end
    k2 = kx^2 + ky^2
    W = (S * exp(-k2 * w.damp)) * cis(kx * w.xdt1 + ky * w.xdt2)
    return _pk_detection_factors(W, kx, ky, w.cKn, w.sKn, w.knife, w.kFactors)
end

# ---------------------------------------------------------------------------------------
# Gaussian-beam disc lobe: tabulated Hankel transform + not-a-knot cubic spline
# ---------------------------------------------------------------------------------------

"Not-a-knot cubic spline on a grid, 'nearest' extrapolation (MATLAB griddedInterpolant 'spline', 'nearest')."
struct ProbeSpline1D
    x::Vector{Float64}
    y::Vector{Float64}
    coefs::Matrix{Float64}       # [n-1 x 4]: ((c1 t + c2) t + c3) t + c4, t = q - x_j
    h::Float64                   # nominal uniform step (interval lookup)
end
function ProbeSpline1D(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n >= 4 || throw(ArgumentError("ProbeSpline1D: need >= 4 points"))
    dx = diff(x)
    divdif = diff(y) ./ dx
    b = zeros(n)
    @inbounds for i in 2:n-1
        b[i] = 3 * (dx[i] * divdif[i-1] + dx[i-1] * divdif[i])
    end
    x31 = x[3] - x[1]
    xn = x[n] - x[n-2]
    b[1] = ((dx[1] + 2 * x31) * dx[2] * divdif[1] + dx[1]^2 * divdif[2]) / x31
    b[n] = (dx[n-1]^2 * divdif[n-2] + (2 * xn + dx[n-1]) * dx[n-2] * divdif[n-1]) / xn
    dl = zeros(n - 1); d = zeros(n); du = zeros(n - 1)
    d[1] = dx[2]; du[1] = x31
    @inbounds for i in 2:n-1
        dl[i-1] = dx[i]
        d[i] = 2 * (dx[i] + dx[i-1])
        du[i] = dx[i-1]
    end
    dl[n-1] = xn; d[n] = dx[n-2]
    s = Tridiagonal(dl, d, du) \ b
    coefs = Matrix{Float64}(undef, n - 1, 4)
    @inbounds for j in 1:n-1
        dzzdx = (divdif[j] - s[j]) / dx[j]
        dzdxdx = (s[j+1] - divdif[j]) / dx[j]
        coefs[j, 1] = (dzdxdx - dzzdx) / dx[j]
        coefs[j, 2] = 2 * dzzdx - dzdxdx
        coefs[j, 3] = s[j]
        coefs[j, 4] = y[j]
    end
    return ProbeSpline1D(x, y, coefs, (x[n] - x[1]) / (n - 1))
end
function (sp::ProbeSpline1D)(q::Real)
    x = sp.x
    n = length(x)
    q <= x[1] && return sp.y[1]
    q >= x[n] && return sp.y[n]
    j = clamp(floor(Int, (q - x[1]) / sp.h) + 1, 1, n - 1)
    @inbounds begin
        while j > 1 && q < x[j]; j -= 1; end
        while j < n - 1 && q >= x[j+1]; j += 1; end
        t = q - x[j]
        return ((sp.coefs[j, 1] * t + sp.coefs[j, 2]) * t + sp.coefs[j, 3]) * t + sp.coefs[j, 4]
    end
end

"""
Transform of the envelope exp(-2 r^2 / w^2) circ(r / R) as a function of x = qR, normalised
to 1 at x = 0; aR = R / w (MATLAB gaussDiscLobe). Tabulated on x = 0:0.05:xMax+1 by the
Hankel integral int_0^1 exp(-2 aR^2 s^2) J0(x s) s ds -- MATLAB uses adaptive `integral`;
here a Gauss-Legendre rule converged to round-off for each x (the integrand is entire), and
the normaliser int_0^1 exp(-2 aR^2 s^2) s ds = -expm1(-2 aR^2) / (4 aR^2) in closed form.
"""
function _probe_gauss_disc_lobe(aR::Float64, xMax::Float64)
    x = collect(0:0.05:xMax+1)
    E = similar(x)
    a2 = 2 * aR^2
    nodesFor(xi) = 48 + 16 * cld(ceil(Int, 0.75 * xi), 16)
    rules = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
    for N in unique(nodesFor.(x))                     # built serially, read-only below
        s0, w0 = gauss(N, 0, 1)
        rules[N] = (s0, w0 .* exp.(-a2 .* s0 .^ 2) .* s0)
    end
    Threads.@threads for i in eachindex(x)
        s, wq = rules[nodesFor(x[i])]
        acc = 0.0
        @inbounds for j in eachindex(s)
            acc += wq[j] * besselj0(x[i] * s[j])
        end
        E[i] = acc
    end
    norm0 = a2 == 0 ? 0.5 : -expm1(-a2) / (2 * a2)
    E ./= norm0
    return ProbeSpline1D(x, E)
end

# ---------------------------------------------------------------------------------------
# fromRealSpace: separable DFT on a Cartesian grid + cubic-convolution interpolation
# ---------------------------------------------------------------------------------------

"""
2-D cubic-convolution interpolant (Keys, a = -1/2; MATLAB griddedInterpolant 'cubic' on a
uniform grid), extrapolation 'none' -> NaN. Rows of `V` go with `ky`, columns with `kx`.
"""
struct ProbeCubicGrid2D
    ky::Vector{Float64}
    kx::Vector{Float64}
    V::Matrix{ComplexF64}
end

@inline function _pk_cell(g::Vector{Float64}, q::Float64)
    n = length(g)
    h = (g[n] - g[1]) / (n - 1)
    j = clamp(floor(Int, (q - g[1]) / h) + 1, 1, n - 1)
    @inbounds begin
        while j > 1 && q < g[j]; j -= 1; end
        while j < n - 1 && q >= g[j+1]; j += 1; end
        return j, (q - g[j]) / (g[j+1] - g[j])
    end
end

# Keys weights of the 4 neighbours j-1, j, j+1, j+2 at fractional position t in [0, 1]
@inline function _pk_keys(t::Float64)
    t2 = t * t; t3 = t2 * t
    return (0.5 * (-t3 + 2t2 - t), 0.5 * (3t3 - 5t2 + 2), 0.5 * (-3t3 + 4t2 + t), 0.5 * (t3 - t2))
end

# fold out-of-range neighbour weights onto the boundary rows (f_0 = 3 f_1 - 3 f_2 + f_3)
@inline function _pk_fold(j::Int, n::Int, w::NTuple{4,Float64})
    w1, w2, w3, w4 = w
    i1, i2, i3, i4 = j - 1, j, j + 1, j + 2
    if i1 < 1          # j == 1: neighbour 0 = 3 f1 - 3 f2 + f3
        w2 += 3w1; w3 -= 3w1; w4 += w1; w1 = 0.0; i1 = j
    end
    if i4 > n          # j == n-1: neighbour n+1 = 3 fn - 3 f(n-1) + f(n-2)
        w3 += 3w4; w2 -= 3w4; w1 += w4; w4 = 0.0; i4 = j + 1
        i1 < 1 && (i1 = 1)
    end
    return (i1, i2, i3, i4), (w1, w2, w3, w4)
end

function (g::ProbeCubicGrid2D)(qy::Float64, qx::Float64)
    ky = g.ky; kx = g.kx
    (qy < ky[1] || qy > ky[end] || qx < kx[1] || qx > kx[end]) &&
        return ComplexF64(NaN, NaN)
    jy, ty = _pk_cell(ky, qy)
    jx, tx = _pk_cell(kx, qx)
    iy, wy = _pk_fold(jy, length(ky), _pk_keys(ty))
    ix, wx = _pk_fold(jx, length(kx), _pk_keys(tx))
    V = g.V
    acc = zero(ComplexF64)
    @inbounds for b in 1:4
        wx[b] == 0 && continue
        col = zero(ComplexF64)
        for a in 1:4
            col += wy[a] * V[iy[a], ix[b]]
        end
        acc += wx[b] * col
    end
    return acc
end

"fromRealSpace's W (MATLAB gridWeight): the interpolated DFT of the nearer lobe, then the factors."
struct ProbeGridWeight
    lobePlus::ProbeCubicGrid2D
    lobeMinus::ProbeCubicGrid2D
    knife::Bool
    cKn::Float64
    sKn::Float64
    kFactors::Bool
    rho::Float64
end
function (w::ProbeGridWeight)(kx::Real, ky::Real)
    _pk_outside_lobes(kx, ky, w.rho) && return zero(ComplexF64)
    kxf = Float64(kx); kyf = Float64(ky)
    W = kxf >= 0 ? w.lobePlus(kyf, kxf) : w.lobeMinus(kyf, kxf)
    isnan(W) && (W = zero(ComplexF64))
    return _pk_detection_factors(W, kxf, kyf, w.cKn, w.sKn, w.knife, w.kFactors)
end

# ---------------------------------------------------------------------------------------
# Constructors (MATLAB static methods)
# ---------------------------------------------------------------------------------------

"""
    probe_point(component = "G13") -> ProbeKernel

No smearing: the nominal +-k0 only (the reference, = the plain Green map). With "G13" it
gives G13_eff.
"""
probe_point(component::AbstractString = "G13") =
    ProbeKernel(ProbePointWeight(), 0.0, 0.0, 1.0, 1.0; component = component, name = "point")

"""
    probe_gaussian(sigmaRelK, sigmaAngDeg; component = "G13", dirStep = nothing, relKStep = nothing)

A +-k symmetric Gaussian kernel in |k| and direction (polar form). sigmaRelK in [0, 0.3),
sigmaAngDeg in [0, 20) deg; a zero sigma is a delta in that coordinate. Nodes out to
PROBE_TRUNC_SIGMA sigma. Default steps: `probe_default_dir_step` / `probe_default_relk_step`.
"""
function probe_gaussian(sigmaRelK::Real, sigmaAngDeg::Real; component::AbstractString = "G13",
                        dirStep = nothing, relKStep = nothing)
    fn = "ProbeKernel.gaussian"
    sR = _pk_check_scalar(sigmaRelK, fn, "sigmaRelK"; lo = 0, hi = 0.3, hiIncl = false)
    sA = _pk_check_scalar(sigmaAngDeg, fn, "sigmaAngDeg"; lo = 0, hi = 20, hiIncl = false)
    dS = _pk_step_opt(dirStep, fn, "DirStep")
    rS = _pk_step_opt(relKStep, fn, "RelKStep")
    dS === nothing && (dS = probe_default_dir_step(sA))
    rS === nothing && (rS = probe_default_relk_step(sR))
    return ProbeKernel(ProbeGaussPolarWeight(sR, sA), PROBE_TRUNC_SIGMA * sA,
                       PROBE_TRUNC_SIGMA * sR, dS, rS; component = component,
                       name = "gaussian($(_pk_fmtg(sR, 3)), $(_pk_fmtg(sA, 3)) deg)",
                       params = (sigmaRelK = sR, sigmaAngDeg = sA))
end

function _probe_radii_pair(r, argName)
    v = r isa Real ? [Float64(r)] : collect(Float64, r)
    (isempty(v) || !all(isfinite, v) || !all(>(0), v)) &&
        _pk_err("ProbeKernel.tgs:$argName", "$argName must be a vector of finite positive reals.")
    length(v) == 1 && return v[1], v[1]
    length(v) == 2 && return v[1], v[2]
    _pk_err("ProbeKernel:tgs:badRadii", "$argName must be a scalar or [along, across].")
end

"""
    probe_tgs(gratingPeriod, sourceRadii, probeRadii = nothing; component = "G13",
              dirStep = nothing, relKStep = nothing)

Transient-grating kernel from the measurement geometry: 1/e^2 intensity radii of the
excitation (`sourceRadii`) and detection (`probeRadii`; `nothing` or empty = much wider
than the source) spots, scalar or [along k, across k], in the unit of `gratingPeriod`.
sigmaRelK = Lambda / (pi sqrt(a_s^2 + a_p^2)), sigmaAng = Lambda / (pi sqrt(b_s^2 + b_p^2))
[rad]; exact Cartesian Gaussians at +-k0. Throws "ProbeKernel:tgs:tooFewFringes" when a
sigma reaches 0.3.
"""
function probe_tgs(gratingPeriod::Real, sourceRadii, probeRadii = nothing;
                   component::AbstractString = "G13", dirStep = nothing, relKStep = nothing)
    fn = "ProbeKernel.tgs"
    L = _pk_check_scalar(gratingPeriod, fn, "gratingPeriod"; lo = 0, loIncl = false)
    aS, bS = _probe_radii_pair(sourceRadii, "sourceRadii")
    if probeRadii === nothing || (!(probeRadii isa Real) && isempty(probeRadii))
        aP = 0.0; bP = 0.0; prb = Float64[]
    else
        aP, bP = _probe_radii_pair(probeRadii, "probeRadii")
        prb = probeRadii isa Real ? [Float64(probeRadii)] : collect(Float64, probeRadii)
    end
    dS = _pk_step_opt(dirStep, fn, "DirStep")
    rS = _pk_step_opt(relKStep, fn, "RelKStep")
    sx = L / (pi * sqrt(aS^2 + aP^2))
    sy = L / (pi * sqrt(bS^2 + bP^2))
    (sx >= 0.3 || sy >= 0.3) && _pk_err("ProbeKernel:tgs:tooFewFringes",
        "The spot is too small for its grating period (sigma = $(_pk_fmtg(sx, 2)), " *
        "$(_pk_fmtg(sy, 2)) of k0).")
    sigAngDeg = rad2deg(sy)
    dS === nothing && (dS = probe_default_dir_step(sigAngDeg))
    rS === nothing && (rS = probe_default_relk_step(sx))
    p = (gratingPeriod = L, sourceRadii = [aS, bS], probeRadii = prb, sigmaRelK = sx,
         sigmaAngDeg = sigAngDeg)
    return ProbeKernel(ProbeGaussCartesianWeight(sx, sy), PROBE_TRUNC_SIGMA * sigAngDeg,
                       PROBE_TRUNC_SIGMA * sx, dS, rS; component = component,
                       name = "tgs($(_pk_fmtg(L, 3)), $(_pk_fmtg(100 * sx, 3))% / " *
                              "$(_pk_fmtg(sigAngDeg, 3)) deg)",
                       params = p)
end

function _pk_detector_args(knifeEdgeDir, kFactors, fn)
    knifeEdgeDir isa Real || _pk_err("$fn:KnifeEdgeDir", "KnifeEdgeDir must be a real scalar.")
    (kFactors isa Bool || kFactors isa Real) ||
        _pk_err("$fn:KFactors", "KFactors must be a logical scalar.")
    knife = !isnan(knifeEdgeDir)
    cKn = knife ? cosd(Float64(knifeEdgeDir)) : 0.0
    sKn = knife ? sind(Float64(knifeEdgeDir)) : 0.0
    return knife, cKn, sKn, kFactors != 0
end

"""
    probe_sras(gratingPeriod, nFringes; dutyCycle = 0.5, detectorOffset = nothing,
               detectorRadius = 0, knifeEdgeDir = 0, kFactors = true, edgeBlur = 0,
               beamRadius = Inf, rings = 2, component = "G13", dirStep = nothing,
               relKStep = nothing) -> ProbeKernel

SRAS kernel: binary fringes in a disc (R = nFringes Lambda / 2), point-like detector
beside the pattern (default `detectorOffset` [1.1 R, 0]), knife edge measuring the slope
along `knifeEdgeDir` (deg from k; NaN = displacement), optional Gaussian beam
(`beamRadius`, then each lobe is the tabulated Hankel transform). EXPERIMENTAL, see the
MATLAB help. Throws "ProbeKernel:sras:tooFewFringes"; warns
"ProbeKernel:sras:evenHarmonics" for a duty cycle other than 0.5.
"""
function probe_sras(gratingPeriod::Real, nFringes::Real; dutyCycle::Real = 0.5,
                    detectorOffset = nothing, detectorRadius::Real = 0.0,
                    knifeEdgeDir::Real = 0.0, kFactors = true, edgeBlur::Real = 0.0,
                    beamRadius::Real = Inf, rings::Integer = 2,
                    component::AbstractString = "G13", dirStep = nothing, relKStep = nothing)
    fn = "ProbeKernel.sras"
    D = _pk_check_scalar(dutyCycle, fn, "DutyCycle"; lo = 0, hi = 1, loIncl = false, hiIncl = false)
    if detectorOffset !== nothing && !isempty(detectorOffset)
        (length(detectorOffset) == 2 && all(isfinite, detectorOffset)) ||
            _pk_err("$fn:DetectorOffset", "DetectorOffset must be [along, across], finite.")
    end
    rd = _pk_check_scalar(detectorRadius, fn, "DetectorRadius"; lo = 0, finite = false)
    bl = _pk_check_scalar(edgeBlur, fn, "EdgeBlur"; lo = 0, finite = false)
    wB = _pk_check_scalar(beamRadius, fn, "BeamRadius"; lo = 0, loIncl = false, finite = false)
    (0 <= rings < length(PROBE_JINC_ZEROS)) ||
        _pk_err("$fn:Rings", "Rings must be an integer in 0..$(length(PROBE_JINC_ZEROS) - 1).")
    knife, cKn, sKn, kF = _pk_detector_args(knifeEdgeDir, kFactors, fn)
    dS = _pk_step_opt(dirStep, fn, "DirStep")
    rS = _pk_step_opt(relKStep, fn, "RelKStep")
    L = _pk_check_scalar(gratingPeriod, fn, "gratingPeriod"; lo = 0, loIncl = false)
    nF = _pk_check_scalar(nFringes, fn, "nFringes"; lo = 0, loIncl = false)
    R = nF * L / 2
    xd = (detectorOffset === nothing || isempty(detectorOffset)) ? [1.1 * R, 0.0] :
         Float64[detectorOffset[1], detectorOffset[2]]
    k0 = 2 * pi / L
    Rt = k0 * R
    rho = max(PROBE_JINC_ZEROS[rings+1] / Rt, 4 * L / (pi * wB))
    rho >= 0.9 && _pk_err("ProbeKernel:sras:tooFewFringes",
        "$(_pk_fmtg(nF)) fringes are too few for $rings rings (or the beam too narrow): " *
        "the lobe would reach $(_pk_fmtg(rho, 2)) k0 from k0 (limit 0.9).")
    abs(D - 0.5) > 1e-12 && @warn("ProbeKernel:sras:evenHarmonics -- Duty cycle " *
        "$(_pk_fmtg(D, 3)): the mask has even harmonics; the 2k0 lobe (apparent velocity " *
        "2v) is not part of the kernel.")
    nH = collect(Float64, -PROBE_SRAS_HARMONICS:PROBE_SRAS_HARMONICS)
    c = [ni == 0 ? D : sin(ni * pi * D) / (ni * pi) for ni in nH]
    xdt1 = k0 * xd[1]; xdt2 = k0 * xd[2]
    rdt = k0 * rd; blurt = k0 * bl
    damp = blurt^2 / 2 + rdt^2 / 8
    p = (gratingPeriod = L, nFringes = nF, radius = R, dutyCycle = D, detectorOffset = xd,
         detectorRadius = rd, knifeEdgeDir = Float64(knifeEdgeDir), kFactors = kF,
         edgeBlur = bl, beamRadius = wB, rings = Int(rings), lobeRadius = rho)
    # shortest oscillation of W in kappa: ring spacing ~ pi / Rt, phase ramp 2 pi / |xdt|
    dq = min(pi / Rt, 2 * pi / max(hypot(xdt1, xdt2), eps())) / PROBE_SAMPLES_PER_OSC
    dS, rS = _probe_osc_steps(dq, dS, rS)
    nm = "sras($(_pk_fmtg(L)), $(_pk_fmtg(nF)) fringes, det $(_pk_fmtg(xd[1] / R, 2)) R)"
    mk(lobe) = ProbeKernel(ProbeSrasWeight(Rt, nH, c, lobe, xdt1, xdt2, damp, knife, cKn, sKn,
                                           kF, rho),
                           asind(rho), rho, dS, rS; component = component, name = nm, params = p)
    if isinf(wB)
        return mk(ProbeJinc())
    else
        return mk(_probe_gauss_disc_lobe(R / wB, Rt * (PROBE_SRAS_HARMONICS + 1 + rho)))
    end
end

"""
    probe_from_real_space(sourceImg, detectorImg, pixelSize, gratingPeriod; lobeRadius = 0.4,
                          knifeEdgeDir = 0, kFactors = true, component = "G13",
                          dirStep = nothing, relKStep = nothing) -> ProbeKernel

Kernel from real-space images of the source pattern S(y, x) and of the detection weight
D(y, x) (rows = y, columns = x, k0 along +x; same size, one pixel grid). W(k) = S^(k) D^(-k)
by a separable direct DFT on a Cartesian grid around each lobe (matrix products, exact at
the grid points), interpolated by cubic convolution; knife-edge and |k| factors analytic,
as `probe_sras`. EXPERIMENTAL (checked against sras only).
"""
function probe_from_real_space(sourceImg::AbstractMatrix{<:Real},
                               detectorImg::AbstractMatrix{<:Real}, pixelSize::Real,
                               gratingPeriod::Real; lobeRadius::Real = 0.4,
                               knifeEdgeDir::Real = 0.0, kFactors = true,
                               component::AbstractString = "G13", dirStep = nothing,
                               relKStep = nothing)
    fn = "ProbeKernel.fromRealSpace"
    rho = _pk_check_scalar(lobeRadius, fn, "LobeRadius"; lo = 0, hi = 0.9, loIncl = false,
                           hiIncl = false)
    knife, cKn, sKn, kF = _pk_detector_args(knifeEdgeDir, kFactors, fn)
    dS = _pk_step_opt(dirStep, fn, "DirStep")
    rS = _pk_step_opt(relKStep, fn, "RelKStep")
    (isempty(sourceImg) || !all(isfinite, sourceImg)) &&
        _pk_err("$fn:sourceImg", "sourceImg must be a nonempty finite real matrix.")
    size(detectorImg) == size(sourceImg) ||
        _pk_err("MATLAB:ProbeKernel:fromRealSpace:incorrectSize",
                "detectorImg must be of size $(size(sourceImg)), got $(size(detectorImg)).")
    all(isfinite, detectorImg) || _pk_err("$fn:detectorImg", "detectorImg must be finite.")
    px = _pk_check_scalar(pixelSize, fn, "pixelSize"; lo = 0, loIncl = false)
    L = _pk_check_scalar(gratingPeriod, fn, "gratingPeriod"; lo = 0, loIncl = false)
    S = Matrix{Float64}(sourceImg)
    D = Matrix{Float64}(detectorImg)
    mS = maximum(abs, S); mD = maximum(abs, D)
    (mS == 0 || mD == 0) &&
        _pk_err("ProbeKernel:fromRealSpace:emptyImage", "Both images must be nonzero somewhere.")
    ny, nx = size(S)
    ixMin = typemax(Int); ixMax = 0; iyMin = typemax(Int); iyMax = 0
    @inbounds for jx in 1:nx, jy in 1:ny
        if abs(S[jy, jx]) > 1e-3 * mS || abs(D[jy, jx]) > 1e-3 * mD
            ixMin = min(ixMin, jx); ixMax = max(ixMax, jx)
            iyMin = min(iyMin, jy); iyMax = max(iyMax, jy)
        end
    end
    Lsup = px * hypot(ixMax - ixMin + 1, iyMax - iyMin + 1)
    k0 = 2 * pi / L
    dq = L / Lsup / PROBE_SAMPLES_PER_OSC
    x = collect((0:nx-1) .* px)
    y = collect((0:ny-1) .* px)
    nq = ceil(Int, rho / dq)
    ky = collect((-nq:nq) .* dq)
    Ey = cis.((-k0 .* ky) .* y')                 # [nky x ny]
    lobes = map((1.0, -1.0)) do sgn
        kx = sgn .+ ky
        Ex = cis.((-k0 .* x) .* kx')             # [nx x nkx]
        W = ((Ey * S) * Ex) .* ((conj.(Ey) * D) * conj.(Ex))   # S^(k) D^(-k), D real
        ProbeCubicGrid2D(ky, kx, W)
    end
    dS, rS = _probe_osc_steps(dq, dS, rS)
    p = (pixelSize = px, gratingPeriod = L, imageSize = [ny, nx], supportDiagonal = Lsup,
         gridStep = dq, knifeEdgeDir = Float64(knifeEdgeDir), kFactors = kF, lobeRadius = rho)
    wf = ProbeGridWeight(lobes[1], lobes[2], knife, cKn, sKn, kF, rho)
    return ProbeKernel(wf, asind(rho), rho, dS, rS; component = component,
                       name = "fromRealSpace($(_pk_fmtg(L)), $(ny)x$(nx) px)", params = p)
end
