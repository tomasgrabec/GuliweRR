# =========================================================================================
# BlasBackend.jl -- which BLAS/LAPACK library the port's dense linear algebra runs on.
# No MATLAB analogue: MATLAB ships its own MKL/Accelerate build and does
# not expose the choice.
#
# Julia's `LinearAlgebra` calls BLAS/LAPACK through libblastrampoline (LBT), a dispatch table
# that can be pointed at any conforming library at run time. Julia's default is OpenBLAS.
# On macOS GuliwerrJ points LBT at Apple's Accelerate framework when it loads (its
# `$NEWLAPACK` symbols, LP64 and ILP64, macOS >= 13.4). The forwarding keeps OpenBLAS
# underneath (`clear = false`), so a routine Accelerate lacks (e.g. `?gemmt`) still falls
# through to OpenBLAS. Measured gain: the Hermitian eigensolvers 1.3-2x. Everything
# that goes through LBT -- `eigen`, `cholesky`, `mul!`, and the hand-written `zhegst`
# `ccall` in RitzSolver.jl -- follows the switch; nothing else in the package changes.
#
# Non-Apple systems: `_init_blas_backend!` does nothing (`@static if Sys.isapple()`), and
# OpenBLAS runs exactly as before.
#
# Opt-out: environment variable `GULIWERRJ_BLAS=openblas`, read when GuliwerrJ loads, skips
# the forwarding -- for debugging and for regression comparison. Results differ between the
# two backends at round-off level only (eigenvalues ~1e-11 relative).
#
# WHY NOT AppleAccelerate.jl (the obvious package, which does the same two `lbt_forward`
# calls in its `__init__`): its AbstractFFTs extension, active whenever FFTW is loaded,
# adds `plan_fft`/`plan_bfft` methods for `Array{ComplexF64}` that are more specific than
# FFTW's and accept power-of-2 lengths only. Loading it made `FFTW.fft` of a length-6 complex
# vector throw "vDSP FFT requires power-of-2 dimensions" for the whole session -- GuliwerrJ's
# own SawSpectrum (length 2987) and any FFT in the user's scripts. Forwarding LBT here
# directly gives the same BLAS/LAPACK behaviour without that side effect.
#
# THREAD CONTROL -- measured, not assumed:
# - `BLAS.set_num_threads` goes through LBT, which knows OpenBLAS's thread setter but NOT
#   Accelerate's: after `BLAS.set_num_threads(1)` Accelerate still reported 12 threads.
# - Accelerate's own switch (`BLASSetThreading`, macOS >= 15) is a single/multi toggle that
#   acts on the CALLING OS THREAD only: set to single on the main thread, a 3000x3000 `dgemm`
#   there took 0.142 s against 0.063 s, while the same `dgemm` on another Julia thread still
#   ran multi-threaded (0.072 s) and reported 12 threads.
# So pinning BLAS around a threaded loop needs two parts: `with_blas_threads` (OpenBLAS + the
# calling thread's Accelerate mode) around the loop, and `with_task_blas_single` inside each
# task's body, which sets and restores the Accelerate mode of the thread the task runs on.
# =========================================================================================

using LinearAlgebra: BLAS
using Base.Libc.Libdl: dlopen_e, dlsym_e

"Environment variable read at load time: `openblas` opts out of Apple Accelerate on macOS."
const BLAS_ENV_VAR = "GULIWERRJ_BLAS"

const _LIBACCELERATE = "/System/Library/Frameworks/Accelerate.framework/Accelerate"
const _BLAS_BACKEND = Ref(:openblas)
# true only when Accelerate is the active backend AND it has the threading API (macOS >= 15)
const _ACCEL_THREAD_API = Ref(false)
const _ACCEL_MULTI_THREADED  = Cint(0)   # vecLib thread_api.h: BLAS_THREADING_MULTI_THREADED
const _ACCEL_SINGLE_THREADED = Cint(1)   #                      BLAS_THREADING_SINGLE_THREADED

"""
    blasBackend() -> Symbol

The BLAS/LAPACK library GuliwerrJ's dense linear algebra runs on in this process:
`:accelerate` (Apple Accelerate, the default on macOS 13.4+) or `:openblas` (every other
platform, older macOS, or the opt-out `ENV["GULIWERRJ_BLAS"] = "openblas"` set before
`using GuliwerrJ`). Results of the two agree to round-off (~1e-11 relative on eigenvalues),
not bitwise.
"""
blasBackend() = _BLAS_BACKEND[]

_accelerate_forwarded() =
    any(l -> occursin("Accelerate", l.libname), BLAS.get_config().loaded_libs)

"macOS product version from the system plist (`nothing` if it cannot be read)."
function _macos_version()
    try
        txt = read("/System/Library/CoreServices/SystemVersion.plist", String)
        m = match(r"<key>ProductVersion</key>\s*<string>([\d.]+)</string>", txt)
        return m === nothing ? nothing : VersionNumber(m.captures[1])
    catch
        return nothing
    end
end

"""
Point LBT at Accelerate (LP64 + ILP64 `\$NEWLAPACK` symbols, OpenBLAS kept as fallback).
Returns a reason string when it cannot, `nothing` on success. macOS 13.4 is the floor
because 13.3's first `\$NEWLAPACK` release had a `dsptrf` bug (as AppleAccelerate.jl).
"""
function _forward_accelerate!()
    v = _macos_version()
    (v === nothing || v < v"13.4") && return "macOS $(something(v, "unknown")) (needs 13.4+)"
    h = dlopen_e(_LIBACCELERATE)
    h == C_NULL && return "the Accelerate framework could not be opened"
    dlsym_e(h, "dgemm\$NEWLAPACK\$ILP64") == C_NULL && return "no ILP64 \$NEWLAPACK symbols"
    BLAS.lbt_forward(_LIBACCELERATE; clear = false, suffix_hint = "\x1a\$NEWLAPACK")
    BLAS.lbt_forward(_LIBACCELERATE; clear = false, suffix_hint = "\x1a\$NEWLAPACK\$ILP64")
    _accelerate_forwarded() || return "libblastrampoline did not accept the forwarding"
    return nothing
end

"Accelerate's per-thread single/multi switch exists (macOS >= 15)."
function _accel_has_threading_api()
    h = dlopen_e(_LIBACCELERATE)
    return h != C_NULL && dlsym_e(h, "BLASSetThreading") != C_NULL &&
           dlsym_e(h, "BLASGetThreading") != C_NULL
end

"Called from `CoreTypes.__init__`: settle the backend once per process."
function _init_blas_backend!()
    @static if Sys.isapple()
        req = lowercase(strip(get(ENV, BLAS_ENV_VAR, "")))
        if !(req in ("", "accelerate", "openblas"))
            @warn "GuliwerrJ: $BLAS_ENV_VAR='$(ENV[BLAS_ENV_VAR])' not recognised " *
                  "(use 'accelerate' or 'openblas'); using Apple Accelerate."
            req = ""
        end
        if req == "openblas"
            if _accelerate_forwarded()   # another package (e.g. AppleAccelerate.jl) did it
                @warn "GuliwerrJ: $BLAS_ENV_VAR=openblas, but Apple Accelerate was already " *
                      "forwarded by another package; it stays active."
                _BLAS_BACKEND[] = :accelerate
                _ACCEL_THREAD_API[] = _accel_has_threading_api()
            else
                _BLAS_BACKEND[] = :openblas
            end
            return _BLAS_BACKEND[]
        end
        why = _accelerate_forwarded() ? nothing : _forward_accelerate!()  # idempotent
        if why === nothing
            _BLAS_BACKEND[] = :accelerate
            _ACCEL_THREAD_API[] = _accel_has_threading_api()
            _ACCEL_THREAD_API[] ||
                @warn "GuliwerrJ: Apple Accelerate is active, but this macOS has no " *
                      "Accelerate threading API (needs 15+), so threaded sweeps cannot pin " *
                      "its BLAS to one thread. Set $BLAS_ENV_VAR=openblas if they run slowly."
        else
            @warn "GuliwerrJ: Apple Accelerate not used ($why); BLAS/LAPACK stay on OpenBLAS."
            _BLAS_BACKEND[] = :openblas
        end
    end
    return _BLAS_BACKEND[]
end

_accel_get_threading() = ccall((:BLASGetThreading, _LIBACCELERATE), Cint, ())
_accel_set_threading(mode::Cint) =
    ccall((:BLASSetThreading, _LIBACCELERATE), Cint, (Cint,), mode)

"Run `f()` with the calling OS thread's Accelerate BLAS single- (`n == 1`) or multi-threaded."
function _with_accel_threads(f, n::Integer)
    (_ACCEL_THREAD_API[] && _BLAS_BACKEND[] === :accelerate) || return f()
    want = n == 1 ? _ACCEL_SINGLE_THREADED : _ACCEL_MULTI_THREADED
    previous = _accel_get_threading()
    previous == want && return f()
    _accel_set_threading(want)
    try
        return f()
    finally
        _accel_set_threading(previous)
    end
end

"""
    with_task_blas_single(f)

Run `f()` with BLAS single-threaded on the OS thread the current task runs on. Call it
INSIDE the body of a threaded loop (each task, each element), inside an enclosing
[`with_blas_threads`](@ref)`(1)`. It matters for Apple Accelerate only, whose threading
switch is per OS thread (this file's header); under OpenBLAS (one process-wide thread count,
already set by `with_blas_threads`) it is just `f()`.
"""
with_task_blas_single(f) = _with_accel_threads(f, 1)

export blasBackend, with_task_blas_single, BLAS_ENV_VAR
