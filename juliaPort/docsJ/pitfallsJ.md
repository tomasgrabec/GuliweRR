# GuliwerrJ — pitfalls

Short warnings specific to the Julia version of GuliweRR. Read this in full before writing a
new calculation script. Every item here is a mistake that fails, or silently does the wrong
thing, in GuliwerrJ. Differences that are only for reference (another optimiser, another
partial-wave order, round-off between linear-algebra backends, the viewers) and the Julia
status of every item in `docs/pitfalls.md` (the MATLAB version) are in `matlabVsJuliaJ.md`.

The items are grouped as in `docs/pitfalls.md` — method-related first, then standard numerics, then other specifics (no Julia-only item is about the physics) — and ordered by severity within each group: 🔴 **silent wrong result** · 🟠 **misreading** (the result is right, its meaning or comparison easy to get wrong) · ⚪ **loud or cost only**. The IDs (`jMTH1`, `jNUM2`, …) are stable: a new item gets the next free number in its group.

## Method-related differences

### jMTH1 🟠 Raw Ritz eigenvectors are gauge-dependent in both languages — `calcDisplacement` is canonical in both

Never compare `ritzSolver.eigVecs_con` element-wise against a MATLAB dump (it depends on an
arbitrary null-space basis); compare phase-invariant quantities, or `calcDisplacement`, which
agrees between the languages for every non-degenerate mode. Detail: `matlabVsJuliaJ.md` §3.3.

### jMTH2 🟠 `GreenSolver`'s partial-wave ordering is canonically sorted in Julia; MATLAB's is raw LAPACK order

Only matters if you read the layered diagnostic fields `k3OutLayered`/`polarizationULayered`/
`partialAmplitudeLayered`: their order can differ from a MATLAB dump even when every number
matches. Detail: `matlabVsJuliaJ.md` §3.2.

### jMTH3 🟠 `ElasticityInversion.calc!` uses a different optimiser from MATLAB's `lsqnonlin`

A projected Levenberg–Marquardt replaces `lsqnonlin` (same tolerances and bounds): converged
results agree, intermediate iterates do not — compare only converged results. Detail:
`matlabVsJuliaJ.md` §3.1.

## Numerics: threads, workspaces, plotting data, linear algebra

### jNUM1 🔴 Julia does not stop you from racing a shared mutable field — MATLAB's `parfor` did

Every `Calc*` result is a `mutable struct`, so it behaves like a MATLAB `handle` object:
`calc(sweep)` fills the object you pass in and returns it. GuliwerrJ's own parallelised loops
cannot race, because each task writes into its own disjoint slot of a preallocated `Vector`,
never a shared scalar field. A loop of your own is another matter: if you write your *own*
`Threads.@threads` loop around GuliwerrJ calls — e.g. accumulating a summary statistic across
several `calc()` results into one mutable struct — Julia gives no equivalent of MATLAB's
compile-time refusal to accept a `parfor` loop that writes a broadcast variable. Two tasks
incrementing the same mutable field can silently lose updates. If you need per-iteration
scratch state in your own parallel code, allocate it fresh inside the loop body, or write into
disjoint elements of a preallocated array (the pattern GuliwerrJ's own sweep classes use) —
never index a shared buffer by `Threads.threadid()`, since a task can migrate threads across a
yield point.

**Your own threaded loop around GuliwerrJ calls** (e.g. `Threads.@threads` over orientations,
materials or directions) has three further traps, none of which raises an error:

- **BLAS threads.** Every solve's BLAS starts threads of its own on top of yours, and the two
  compete for the same cores (in one test, a batch of `eigen` solves on 8 threads ran 1.7×
  slower with 8 BLAS threads than with 1). Wrap the loop in `with_blas_threads(1) do ... end`
  and each task's body in `with_task_blas_single(() -> ...)`; on macOS the second is needed because
  Apple Accelerate's thread switch acts per OS thread and ignores `BLAS.set_num_threads`.
  `forEachElement`/`forEachElementPooled` do both for you; `largeCalculationsJ.md` §5 shows
  the pattern.
- **`ProbeKernel` with a weight function of your own.** From 4096 quadrature nodes, `sample`
  (and so every spot-averaged Green calculation) calls the weight function from several
  threads at once — inside or outside your own loop. A function that writes to a shared
  variable (a counter, a cache, a buffer) races silently; keep it pure. The built-in kernels
  are safe.
- **`RitzWorkspace`**: one per task, never shared — see jNUM2.

Nesting itself is safe (Julia's tasks compose), but put the threads where the work is: give a
sweep inside your loop `useParallel = :never`, and a direct `calc(GreenSolver, ...)`
`threaded = false` to save its task overhead (results are bitwise the same either way).

Reference detail (MATLAB vs Julia threading, side by side): `matlabVsJuliaJ.md` §4.

### jNUM2 🔴 A `RitzWorkspace` is scratch memory: never share one between tasks, and `!` outputs live only until its next call

`RitzWorkspace` (with `ritzSAWPrecompute!` / `ritzSensitivities!`, Julia only) avoids the few
MB of scratch arrays every Ritz precompute and solve allocates. MATLAB has no counterpart, so
nothing in MATLAB habits warns about its two rules:
- **One workspace per task.** Two tasks writing the same workspace at once corrupt each other's
  matrices silently — or crash inside LAPACK. In a threaded loop use `forEachElementPooled`
  (a `Channel` pool; each element takes a workspace and gives it back). Never pick a workspace
  by `Threads.threadid()`: a task can migrate between threads (jNUM1).
- **`ritzSAWPrecompute!` returns workspace arrays.** Its `dKdC_cell`, `d2KdCdD_cell`,
  `ritzStruct.M_con`/`M_conU` are overwritten by the next `ritzSAWPrecompute!` on that
  workspace. Use them at once (one precompute per wavelength, then the angles), never keep the
  NamedTuple across wavelengths. `ritzSensitivities!` returns fresh copies, which you may keep.
- **In your own threaded loop, mark the task's workspace variable `local`.** If the same name
  is assigned anywhere else in the enclosing function (e.g. a serial branch's
  `ws = RitzWorkspace(...)`), the `Threads.@threads` body's `ws = take!(pool)` assigns that one
  captured variable, shared by every task — a data race that can end in a `PosDefException`
  or a crash inside LAPACK.

### jNUM3 🔴 Cell-centring of maps: Makie centres Cartesian maps itself, polar maps need `pcolorCenteredGrid`, and the matrix orientation differs from Makie's

MATLAB's `pcolor` treats a sample as its cell's lower-left corner (NUM3 of
`docs/pitfalls.md`). Makie's `heatmap!(x, y, C)`, when `length(x) == size(C,1)` (sample
*centres*, the natural way to call it), already performs exactly the edge construction
`pcolorCenteredGrid` exists to force onto MATLAB's `pcolor` — interior edges at consecutive
midpoints, endpoints linearly extrapolated; this holds for increasing, decreasing, and
single-sample axes. So `plotGreenMap` and friends pass angle/velocity/kParRange arrays
straight to `heatmap!`, and calling `pcolorCenteredGrid` first would be actively wrong (it
would double-apply the correction). The helper is **not** dead code, though: `plotPolarHeatmap`
still needs it, because drawing a polar map as a hand-built quad mesh (`poly!`) genuinely
needs explicit cell *edges*, not centres, to build the quads from.

The part to be careful of if you call the lower-level `PlotData`/`VisData` functions directly
and draw your own `heatmap!`: this package's data convention (matching MATLAB's
`pcolor(x,y,C)`/`meshgrid(x,y)`) stores a map as `C[i,j] <-> (x[j], y[i])` — row is the
*y*-index, column is the *x*-index. Makie's `heatmap!(x,y,z)` convention is the opposite,
`z[i,j] <-> (x[i], y[j])`. Every built-in `plot*` function already transposes before calling
`heatmap!`; a script drawing its own heatmap from `ritzSpectrumMapData` or similar and
skipping the transpose will get a map that runs sideways, not an error.

### jNUM4 🟠 On macOS the linear algebra runs on Apple Accelerate — results differ from OpenBLAS at round-off level

Mac and Linux results (or runs with `GULIWERRJ_BLAS=openblas`) agree to ~1e-11, not bitwise:
never compare them by `==`. Do not `using AppleAccelerate` in a session with FFTW — it breaks
complex FFTs of non-power-of-2 length in your own code. Detail and opt-out:
`matlabVsJuliaJ.md` §3.4; thread control: jNUM1.

### jNUM5 ⚪ Julia threads are a process-start flag — a forgotten `-t auto` runs correctly, just silently serial

MATLAB's `parfor` acquires a pool at runtime (`parpool`); Julia's thread pool is sized once,
when the process starts, and cannot be grown afterwards. **Always launch with `-t auto`** (or
a specific count): `julia --project=... -t auto`. A session started without it gives
`Threads.nthreads() == 1`, and `useParallel = :auto`/`:always` in every sweep constructor
then simply runs serially — correct results, no error, just slower, which is a much easier bug
to miss than a crash. Check `Threads.nthreads() > 1` once per session if a sweep feels
unexpectedly slow.

## Other specifics: calling conventions and setup

### jOTH1 ⚪ Constructors are keyword-only — a positional call fails with a confusing error, not a helpful one

`Structure(mat, cp, "name", 500)` (the natural translation of MATLAB's
`Structure(mat, cp, 'name', 500)`) does not raise "use keywords instead" — it silently matches
the struct's own auto-generated **positional** inner constructor (four fields: `name`,
`nMaxModes`, `mater`, `calcParams`, in that field order, not constructor-argument order) and
fails deep inside type conversion:

```
MethodError: Cannot `convert` an object of type Material to an object of type String
```

The same trap exists for `CalcParams`, `Material` and every `CalcSingleXxx`/`CalcAngSweepXxx`
family member — all of them are keyword-only on purpose. Always use the keyword form:
`Structure(materIn = mat, calcParamsIn = cp, name = "name", nMaxModes = 500)`.

`ElastConsts(consts, symmNo)` is the exception: it takes its two arguments positionally, as in
MATLAB. So do the inverse launchers (`measured`, `association`, `sensitivity`, `elasticity`)
and `freqToOmega0`.

### jOTH2 ⚪ `calcMode`/`retentionLevel`/`useParallel`/`waveType` options are `Symbol`s, not `String`s — a quoted MATLAB literal fails or is silently accepted, depending on the argument

Where MATLAB takes `'ritzOnly'`, Julia takes `:ritzOnly`. Passing the quoted string instead
behaves differently from one argument to the next:

- `angSweepHS(structure = str, calcMode = "ritzOnly")` throws `ArgumentError:
  CalcAngSweepHS:constructor:badMode -- mode must be 'full' or 'ritzOnly'; got 'ritzOnly'`.
  The message prints the offending value in quotes, so it can look as if you had typed the
  right thing.
- Where the keyword is typed `Symbol` (e.g. `retentionLevel` of `angSweepPlate`), the call
  throws `TypeError: in keyword argument retentionLevel, expected Symbol, got a value of type
  String`.
- `useParallel = "auto"` on `angSweepHS` is silently accepted (it is converted to a `Symbol`).

Always pass Symbols: `:full`/`:ritzOnly`, `:full`/`:saw_only`/`:none`,
`:auto`/`:always`/`:never` — the MATLAB spelling with a colon in front and no quotes. The one
deliberate exception is `MeasuredVelocities`'s `waveType`, which is a plain `String`
(`"rayleigh"`, `"pseudoRayleigh"`, ...), because it is stored directly in a struct field typed
`String`.

### jOTH3 ⚪ Dispersion `curveType` names use underscores, not hyphens — a hyphen is not legal inside a `Symbol` literal

MATLAB's `curveType` strings are hyphenated (`'wavelength-freq'`, `'k-vphase'`,
`'freq-vgroup'`, ...). Julia's `Symbol` cannot express that directly — `:wavelength-freq`
parses as *subtraction* of the symbol `:wavelength` and the name `freq` — so every dispersion
`curveType` accepted by `plotDispersionCurves`/`plotDispersionCurvesPlate`/`dispersionViewer`
uses an underscore instead: `:wavelength_freq`, `:wavelength_vphase`, `:k_freq`, `:k_vphase`,
`:freq_vphase`, `:freq_vgroup`. Typing the hyphenated form does not "just work" the way copying
a MATLAB string literal might tempt you to assume: it parses, then fails when run with an
`UndefVarError` (or a `MethodError`) instead of a message that names the mistake.

### jOTH4 ⚪ GUI launchers need a loaded Makie backend

`structureControl`, the viewers and the `plot*` functions throw a plain `MethodError` until `using GLMakie` (or `using CairoMakie`) has been run in the same session; load it after `using GuliwerrJ`.
