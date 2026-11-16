# =========================================================================================
# ParallelPolicy.jl
#
# Port of guliwerr/core/ParallelPolicy.m -- shared parallel-execution policy for GuliweRR
# sweep classes. MATLAB's version makes two decisions in one place: whether a loop of
# independent iterations runs in parallel ("auto"/"always"/"never"), and acquisition of a
# MATLAB THREAD pool (`ensureThreadPool`). Only the first decision survives unchanged; the
# second has no Julia analogue at all -- see below.
#
# NO per-file module (see CoreTypes.jl's header): this
# file's definitions live directly in `CoreTypes`.
#
# STATELESS -- MATLAB's `ParallelPolicy` is a `classdef` with only `methods (Static)` and one
# `properties (Constant)`, i.e. a namespace, not a value/handle class at all. No struct is needed here --
# `should_use_parallel`, `threads_available` and `with_blas_threads` are plain module-level
# functions, and `AUTO_THRESHOLD` is a plain `const`.
#
# `ensureThreadPool` HAS NO JULIA ANALOGUE: MATLAB's `parpool`/thread-pool machinery acquires a runtime
# resource that Julia has no equivalent of -- Julia's threads are fixed at PROCESS START via
# `-t auto`/`-t N` or the `JULIA_NUM_THREADS` environment variable, never acquired mid-session.
# `CalcCommon.resolveParallel`'s two-part MATLAB contract ("shouldUseParallel and
# ensureThreadPool must travel together") collapses to a single predicate:
# `should_use_parallel(setting, n) && threads_available()`. The `false` branch's MEANING ("run
# serial") is unchanged; only the acquisition step disappears, because there is nothing left
# to acquire -- `threads_available()` only ever REPORTS whether `-t` was passed at startup.
#
# BLAS-THREAD HELPER (kept here, not in the sweep classes): a hazard with NO MATLAB
# equivalent. MATLAB's thread pool and its BLAS calls never compete for the same cores in a
# way ParallelPolicy.m had to manage; Julia's `Threads.@threads` tasks and BLAS's own
# internal threading do, because both draw from the same OS thread budget. Measured:
# 8 Julia threads x 32 complex 300x300 `eigen` solves, BLAS threads = 1 ->
# 0.082s, BLAS threads = 8 -> 0.142s (1.7x SLOWER) -- 8x8 = 64 OS threads oversubscribing 8
# cores. Rule: pin BLAS to 1 thread for the duration of a THREADED BATCH of solves (the sweep
# case); leave BLAS multithreaded for a single large solve (one big `eigen`, not a batch).
# `with_blas_threads` below is the save/set/restore helper, generic over
# both cases via its `nthreads` argument.
# =========================================================================================

using LinearAlgebra: BLAS

"""
    AUTO_THRESHOLD

`"auto"` mode goes parallel once the loop item count is `>= AUTO_THRESHOLD`. A single
benchmarked threshold, not scaled by `legO`/`nVelGreen` -- port of
`guliwerr/core/ParallelPolicy.m`'s `properties (Constant) AUTO_THRESHOLD = 20`
(`ParallelPolicy.m`), turned into a `const` (MATLAB `Constant` property -> Julia `const`).
"""
const AUTO_THRESHOLD = 20

"""
    should_use_parallel(useParallelSetting, nItems::Integer) -> Bool

Decide whether a loop of `nItems` independent iterations should run in parallel.

INPUT: `useParallelSetting` -- `:auto`, `:always` or `:never` (a `Symbol`; an `AbstractString`
such as `"auto"` also works, converted via `Symbol` -- MATLAB accepted `char` and `string`
interchangeably via `string(useParallelSetting)`, Symbol is the direct Julia analogue of a
fixed small set of named modes); `nItems` -- number of independent loop iterations (must be
`>= 0`).

Throws `ArgumentError` if `nItems` is negative, or if `useParallelSetting` is none of
`:auto`/`:always`/`:never` (mirroring `ParallelPolicy.m`'s `validateattributes` guard and its
`otherwise` branch's `error(...)`, respectively -- Julia has no error-identifier system, see
`Material.jl`'s `SetElastConsts` for the same "message only" pattern elsewhere in this port).

Port of `guliwerr/core/ParallelPolicy.m`.
"""
function should_use_parallel(useParallelSetting::Symbol, nItems::Integer)
    nItems >= 0 ||
        throw(ArgumentError("ParallelPolicy:shouldUseParallel:badInput -- nItems must be " *
            "a nonnegative scalar, got $nItems"))
    return useParallelSetting === :never  ? false :
           useParallelSetting === :always ? true  :
           useParallelSetting === :auto   ? nItems >= AUTO_THRESHOLD :
           throw(ArgumentError("ParallelPolicy:shouldUseParallel:badSetting -- " *
               "useParallelSetting must be :auto, :always or :never; got :$useParallelSetting"))
end

"""
    should_use_parallel(useParallelSetting::AbstractString, nItems::Integer) -> Bool

Convenience overload: `"auto"`/`"always"`/`"never"` (or any other string) is converted to a
`Symbol` and dispatched to [`should_use_parallel(::Symbol, ::Integer)`](@ref) — an invalid
string still reaches that method's `ArgumentError` branch. Mirrors MATLAB's `char`/`string`
interchangeability (`ParallelPolicy.m`'s `useParallelSetting = string(useParallelSetting)`).
"""
should_use_parallel(useParallelSetting::AbstractString, nItems::Integer) =
    should_use_parallel(Symbol(useParallelSetting), nItems)

"""
    threads_available() -> Bool

`true` if this Julia process was started with more than one thread (`-t auto`/`-t N` with
`N > 1`, or `JULIA_NUM_THREADS`).

Replaces `ParallelPolicy.ensureThreadPool` (see this file's header comment for why there is
no acquisition step to port, only a report of what was fixed at process start).
`should_use_parallel(setting, n) && threads_available()` is the Julia equivalent of MATLAB's
"decide, then ensure a thread pool exists" two-step.
"""
threads_available() = Threads.nthreads() > 1

"""
    with_blas_threads(f, nthreads::Integer=1)

Run `f()` with `LinearAlgebra.BLAS.set_num_threads(nthreads)` in effect, restoring the prior
BLAS thread count afterwards -- even if `f` throws.

Use `with_blas_threads(1) do ... @threads for ... ; eigen(...) ; end ... end` around a
THREADED BATCH of independent solves (otherwise a measured 1.7x slowdown, from
Julia threads and BLAS's own internal threading oversubscribing the same cores). Do NOT wrap
a single large solve this way -- leave BLAS multithreaded there, since it is not competing
with anything.

No MATLAB analogue: MATLAB's thread pool and its BLAS calls do not compete for cores the way
Julia's `Threads.@threads` and BLAS's internal threading do (see this file's header comment).

Apple Accelerate (the macOS backend, `BlasBackend.jl`): `BLAS.set_num_threads` does not reach
it, and its own single/multi switch acts on the calling OS thread only. This function sets
that switch for the calling thread as well; the body of every task in a threaded loop must
additionally be wrapped in [`with_task_blas_single`](@ref) (as `forEachElement` does).
"""
function with_blas_threads(f, nthreads::Integer = 1)
    previous = BLAS.get_num_threads()
    BLAS.set_num_threads(nthreads)
    try
        return _with_accel_threads(f, nthreads)
    finally
        BLAS.set_num_threads(previous)
    end
end

export AUTO_THRESHOLD, should_use_parallel, threads_available, with_blas_threads
