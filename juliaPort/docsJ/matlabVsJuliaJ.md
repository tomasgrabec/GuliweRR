# GuliwerrJ — differences from MATLAB GuliweRR

One reference page for someone who knows the MATLAB version: everything a user of GuliwerrJ
meets that is spelled, computed, threaded or displayed differently. The physics, the
defaults and the results are the same (round-off apart); this page lists only what differs.
The mistakes that fail or silently do the wrong thing are in `pitfallsJ.md` — read that
first; this page is the background to it. Each function's docstring (`?name`) has the full
detail.

## 1. Calling conventions

| Topic | MATLAB | Julia |
|---|---|---|
| Constructors | Positional, `nargin`-guarded: `Material(ec, 8.96)`, `Structure(mat, cp, 'cu', 500)` | **Keyword-only**: `Material(elConsts = ec, rho = 8.96)`, `Structure(materIn = mat, calcParamsIn = cp, name = "cu", nMaxModes = 500)`. A positional call fails with a confusing `MethodError` (`pitfallsJ.md` jOTH1) |
| Positional exceptions | — | `ElastConsts(consts, symmNo)`, the inverse launchers (`measured`, `association`, `sensitivity`, `elasticity`), `freqToOmega0(freqHz, lambdaRefNm)`, `ProbeKernel(weightFcn, angHalfWidth, relKHalfWidth, dirStep, relKStep; ...)`, and the `calc`/`calcAll`/`classify` step functions (`calc(sweep, mode, retentionLevel)`) |
| Launchers | `Calc.angSweepHS(...)`, `Vis.plotGreenMap(...)`, `Gui.angSweepViewer(...)`, `Inv.association(...)` | Bare exported names, no namespace prefix: `angSweepHS(; ...)`, `plotGreenMap(...)`, `angSweepViewer(...)`, `association(...)`. `angSweepHS()` with no arguments works, as `Calc.angSweepHS()` does |
| Method calls | `sweep.calc()`, `cs.classify('ritz', 'high')` | `calc(sweep)`, `classify(cs, :ritz, :high)` — the object is the first argument |
| Option values | Quoted strings: `'ritzOnly'`, `'saw_only'`, `'high'`, `'log'` | `Symbol`s: `:ritzOnly`, `:saw_only`, `:high`, `:log`. A string fails or is accepted depending on the argument (`pitfallsJ.md` jOTH2). `useParallel` also accepts a string. `MeasuredVelocities`' `waveType` stays a `String` (`"rayleigh"`) |
| `curveType` names | Hyphenated: `'wavelength-freq'`, `'k-vphase'` | Underscored: `:wavelength_freq`, `:k_vphase` (`pitfallsJ.md` jOTH3) |
| "Not given" | `[]` or `''` | `nothing` (e.g. `sawClassMethod = nothing`) |
| Cell argument `{consts, symmNo}` for `Material`'s elastic constants | Cell array | `Tuple`: `(consts, symmNo)` or `(consts,)` |
| Value vs handle objects | `ElastConsts`, `Material`, `CalcParams`, `Structure`, the Ritz and Green solvers are value classes; the `Calc*` results are handles | The same split: the first group are immutable `struct`s, the `Calc*` results (and `ElasticityInversion`) are `mutable struct`s. `calc(sweep)` fills `sweep` in place and returns it; a setter or a solver's `calc` returns a **new** object, so assign it: `mat = SetDensity(mat, 2.7)`, `cp = set_rel_depth(cp, 30)`, `s = calc(s)` |
| `CalcParams` setters | `cp.SetRelDepth(30)` style | snake_case free functions: `set_rel_depth`, `set_real_thickness`, `set_leg_order` |
| Derived quantities | Properties: `str.nLay`, `str.strBasisVecL` | Functions: `nlay(str)`, `str_basis_vec_l(str)` |
| `ElastConsts.Get*` instance wrappers | `ec.GetReducedElastConsts(rotMat)`, `ec.GetYoungsModulus(...)`, `ec.GetHillAverages`, `ec.GetAnisotropyRatio`, ... | Not ported; call the free function on the struct's fields: `ElastRotReductMat(rotMat, ec.symmClassNo) * ec.vector`, `Tensor2Youngs(...)`, `HillAverages(...)`. The anisotropy ratio is the field `ec.shearAnisotropyRatio` (cubic/isotropic closed form) or `calc_anisotropy_ratio(normElTen)` (general angular scan; the one renamed function). `Material`'s `GetRotMat`, `GetRotElastTensor`, `GetRotElastVoigt`, `GetRotDirMat` exist, as functions taking the `Material` first |
| Static methods | `ElastConsts.Voigt2Vector(...)` | Plain functions: `Voigt2Vector(...)` |
| Results with several outputs | `nargout`-dependent: `[vg, solver] = fastCalcGroupVel(...)`; `GreenSolver.calc` returns the displacement field when `nargout > 1` | Fixed return values: `fastCalcGroupVel` always returns `(velGroup, solver)`; `GreenSolver`'s `fullOutput = true` keyword fills `dispField` (and the extra field `depthAxis`) |
| `GammaMatrix` | `calc` (half-space) and `calcLamb` (plate) | One `calc`, chosen by a first argument `SAWBoundary()` or `LambBoundary()` |
| `findRoots` | Three parallel cell arrays `{iwt}` (velocities, polarisations, angles) | `Vector{LimitingRoot}`, each with its own `iwt` (1 longitudinal, 2 fast shear, 3 slow shear) |
| `MeasuredVelocities.fromTable` | A MATLAB `table` | A `NamedTuple` of columns: `fromTable((angleDeg = a, velocity = v))` |
| Data shapes | Row/column vectors, cell arrays, `struct`s | Plain `Vector`s (no row/column orientation), `Vector`s of `Vector`s or `Matrix{Matrix{Float64}}` for cell arrays (displacement fields are indexed `[iDir, ilay]`), `NamedTuple`s for `struct`s (e.g. `latticeParams`). Indices are 1-based and use square brackets |
| `guliwerrVersion()` | The GuliweRR version | The same string; the package's `Project.toml` version tracks the packaging only |

## 2. Errors and warnings

| Topic | MATLAB | Julia |
|---|---|---|
| Error type | `error('Class:method:id', ...)` with an identifier you can catch by ID | A Julia exception (`ArgumentError`, `DomainError`, ...) whose message starts with the same identifier: `ArgumentError("ModeAssociation:fromSweep:angleNotFound -- ...")`. Catch by type and test the message (`occursin("angleNotFound", sprint(showerror, e))`) |
| Wrong argument type | `validateattributes` error | Usually a `MethodError` (no method for that argument type) or a `TypeError` on a typed keyword |
| `phaseVels` on a tensor whose Christoffel matrix has a negative eigenvalue (a non-positive-definite, physically unstable tensor) | Returns complex velocities silently | Throws `DomainError` from `sqrt`. A physically valid tensor never triggers it |
| Warnings | `warning('Class:method:id', ...)`, silenced with `warning('off', id)` | `@warn` with the same identifier at the start of the message. There are no warning IDs to switch off: silence a block with `using Logging; with_logger(NullLogger()) do ... end`. Most warnings fire on every call, as in MATLAB; the ones MATLAB shows once per session (`legOrder` too low for the Ritz-only classification, Green axis step too coarse, near-isotropic Green notice) are also once per session here |

## 3. Numerical differences

Results agree with MATLAB to round-off. The differences below are in the path taken, in the
order of intermediate quantities, or at the round-off level itself; none changes a physical
result.

### 3.1 Inversion optimiser

MATLAB's `lsqnonlin` (`'trust-region-reflective'`) has no open-source Julia equivalent.
`ElasticityInversion.calc!` uses a projected Levenberg–Marquardt on the residual vector and
its Jacobian: the same Gauss–Newton family, the same tolerances and the same box bounds
(applied by projection). The two converge to the same result, but the iterates differ at
every step — compare converged results, never intermediate optimiser state.

The plate ZGV search (`locateZGVs`, `locateZGV`) uses no optimiser: both languages run the
same Newton iteration (`ZgvNewton.solve`) step for step and return the same points to ~1e-12.

### 3.2 `GreenSolver` partial-wave ordering

MATLAB keeps the six roots `k3` of each layer's sextic in the raw order `roots()` returns,
which depends on the LAPACK the MATLAB process links. GuliwerrJ sorts each velocity column's
six roots into one canonical `(real, imag)` order straight after the root-find. `G`, `GFull`
and `dispField` are sums over partial waves and do not change. Only the order inside the
layered diagnostic fields `k3OutLayered`, `polarizationULayered` and
`partialAmplitudeLayered` (with `fullOutput = true`) can differ from a MATLAB dump of the same
calculation. For an isotropic material the sextic has exact double roots, and which member of
such a pair comes first is undetermined in either language.

The near-singular rescaling of the boundary matrix uses a power-of-two row/column
equilibration instead of MATLAB's `equilibrate` (MC64 permutation and scaling, which Julia
lacks). Both give the same inverse (`G` agrees to ~1.5e-16 relative on a relDepth-160 case);
the count of near-singular points in the `GreenSolver:calc:nearSingular` warning can differ.
In layered calculations Julia can flag points where MATLAB flags none (tutorial 07 prints such
warnings in Julia only); the layered Green maps themselves are compared with MATLAB's in the
golden tests.

### 3.3 Eigenvector sign and gauge

| Quantity | Reproducible between MATLAB and Julia? | Compare as |
|---|---|---|
| Christoffel eigenvectors (`matrix`, `SortedEig`) | Sign: no, even away from a degeneracy. In one comparison on Cu(001) (19 angles, eigenvalues well separated), 15 of the 19 angles had at least one column sign-flipped, with eigenvalues agreeing to 1e-12. A degenerate shear pair (isotropic material) has an arbitrary basis | Rays (up to sign); a degenerate pair as a subspace |
| Raw Ritz eigenvectors `eigVecs_con` | No. They live in a constraint null-space basis, and MATLAB's `null()` and Julia's `nullspace()` return different, equally valid orthonormal bases; velocities agree to ~6e-15 regardless | Phase-invariant quantities only: `abs(u' * v)`, outer products, or a comparison that divides out the phase |
| `calcDisplacement`, modal shapes and every plot built from them | Yes, up to round-off, for every non-degenerate mode: both languages rotate each mode so that its largest in-plane Legendre coefficient is real and positive | Element-wise |
| Velocities, scores, classifications | Yes, to round-off | Element-wise, with a tolerance |

### 3.4 Linear-algebra backend (Apple Accelerate)

On macOS 13.4 or later, `using GuliwerrJ` points Julia's BLAS/LAPACK at Apple's Accelerate
framework, for the whole session (your own `eigen`/`*` calls included); everywhere else it
stays on OpenBLAS. `blasBackend()` says which is active.

| Topic | Detail |
|---|---|
| Agreement | Round-off, not bitwise: eigenvalues move by ~1e-11 relative between the backends (and between macOS and Linux). Never compare a Mac result with a Linux one, or with a run under the other backend, by `==` |
| Opt-out | `GULIWERRJ_BLAS=openblas` in the environment before `using GuliwerrJ` (read once, at load). Any value other than `openblas`/`accelerate` warns and keeps Accelerate |
| Thread control | `BLAS.set_num_threads` does not reach Accelerate, whose single/multi-thread switch acts per OS thread. GuliwerrJ's own loops handle it; for your own loops see §4 |
| Old macOS | Below 13.4: OpenBLAS, with a one-time warning. 13.4–14.x: Accelerate without its threading API (threads cannot be pinned), with a one-time warning |
| `AppleAccelerate` package | Do not load it in a session with FFTW: its FFT extension overrides FFTW's complex `plan_fft` and throws for any length that is not a power of 2. GuliwerrJ does not use that package and protects its own FFTs; FFTs in your own code are not protected |

### 3.5 Remaining round-off and path differences

| Where | Difference | What to do |
|---|---|---|
| Limiting velocities (`findRoots`, `refineRoot`) | MATLAB's `fminbnd` is replaced by a Brent minimiser (`fminbnd_brent`) with the same tolerance; it converges to the same minimum by a slightly different iterate path, so refined limiting velocities agree to ~1e-8 relative, not to the Ritz eigenvalues' precision | Anything that thresholds on a limiting velocity: compare the resulting mask or label, which is stable, not the boundary value |
| `ritzSensitivities` | Julia computes only the lowest `nModes` eigenpairs; MATLAB computes all and keeps the lowest. On stiff thin-film stacks the two differ at the ~1e-8 level, the lowest-mode path being the more accurate (~1e-11) | Keep `nModes` small; compare with a ~1e-8 relative tolerance |
| Contours (`extractPlanarCut`, the plane cuts of the 3D surfaces) | Contour.jl replaces `contourc`: the same segments and vertices to machine precision, but a closed loop may start at a different point and either loop or open segment may run the other way | Compare point sets order-insensitively |
| `rcond_matlab` | Reproduces MATLAB's `rcond`: `0.0` on an exactly singular matrix and `NaN` on non-finite input, where Julia's own `lu`/`gecon!` would throw | Nothing — it behaves as MATLAB |
| `matlab_roots` | Reproduces MATLAB's `roots`: unsorted companion-matrix eigenvalues in LAPACK order. Exactly or nearly degenerate roots carry 1e-5 to 1e-8 of conditioning noise, in both languages alike | Do not compare degenerate roots tighter than that |
| Sweep execution order | Same results, different scheduling: a `:full` half-space sweep runs Ritz then Green per angle in one threaded loop (MATLAB runs every angle's Ritz solve before any Green solve, so holds all eigenvectors at once); the 2D layered sweep nests its wavelength and angle loops (no `flattenSweep` option); the plate 2D sweep threads over every (wavevector, angle) point; the first inversion pass is angle-parallel | Nothing — results are independent of the order |

## 4. Threading model

| Topic | MATLAB | Julia |
|---|---|---|
| Thread pool | Acquired at run time (`parpool`, `ParallelPolicy.ensureThreadPool`) | Fixed when the process starts: `julia -t auto` (or `-t N`, or `JULIA_NUM_THREADS`). Without it `Threads.nthreads() == 1` and every sweep runs serially, correctly and silently (`pitfallsJ.md` jNUM5) |
| `useParallel` | `"auto"` (parallel from 20 items), `"always"`, `"never"` | `:auto` (parallel from `AUTO_THRESHOLD = 20` items), `:always`, `:never` |
| Nested parallel loops | Forbidden (`parfor` inside `parfor`), so outer loops stay serial | Allowed and composable: tasks nest without oversubscribing threads |
| Shared variables in a loop | `parfor` refuses a loop that writes a broadcast variable | Nothing stops a data race; write only to disjoint slots (`pitfallsJ.md` jNUM1) |
| Worker memory | Each `parfor` worker gets its own copies | Threads share memory: no copies, and no protection either |
| BLAS threads | Not managed | A threaded batch of solves runs fastest with BLAS on one thread: wrap the loop in `with_blas_threads(1) do ... end` and each task's body in `with_task_blas_single(() -> ...)` (the second matters on Accelerate, §3.4). Leave a single large solve alone. GuliwerrJ's own loops (`forEachElement`, `forEachElementPooled`) do this already |
| `GreenSolver` | No threading inside one call | `calc(GreenSolver, ...)` spreads velocity chunks over the threads by default; the result is bitwise independent of the thread count. Inside your own threaded loop over directions, `threaded = false` saves the (small) task overhead; the default is safe either way |
| `ProbeKernel` weight function | Called from one thread | With 4096 quadrature nodes or more, `sample` (used by every spot-averaged Green calculation) calls the weight function from several threads at once: a weight function of your own must not write to shared state. The built-in kernels are safe |
| `RitzWorkspace`, `ritzSAWPrecompute!`, `ritzSensitivities!` | No counterpart (`eig` allocates per call) | Julia only: one workspace per task (`forEachElementPooled`), and arrays returned by `ritzSAWPrecompute!` live only until the next call on that workspace (`pitfallsJ.md` jNUM2) |
| Batches of independent sweeps | `parfor` over the batch | Worker processes (`Distributed`) scale better than threads for long batches; see `largeCalculationsJ.md` |

## 5. Not ported, and options that differ

| MATLAB | Julia |
|---|---|
| `Calc.saveResult`, `Calc.loadResult`, `.mat` result files | Not ported. A result lives in memory for the session; MATLAB-saved results cannot be read. Store what you need with `Serialization` (same Julia version only) or your own format |
| `Gui.viewSaved`, `Gui.viewResult` | Not ported; open a viewer on the in-memory result (`angSweepViewer(sweep)`, `combinedViewer(sweep)`, ...) |
| `ElastConsts.Get*` wrappers, `SetNewSymmClass` | Not ported (§1) |
| `CalcWavevecAngSweep2Layered.flattenSweep` | Not needed: Julia nests the two loops (§3.5) |
| `ParallelPolicy.ensureThreadPool` | No counterpart: threads are fixed at start-up (§4) |
| `retentionLevel` defaults | Same as MATLAB: `:full` on single angles and angular sweeps, `:none` on the 2D sweeps (set `:full` for branch tracking and the dispersion viewer); plates accept only `:full`/`:none` |
| `greenComponentData(result, 'spot')` and other components | Same as MATLAB: only `G13eff` is stored, any other component is one Green solve per angle; `"spot"` is experimental, as in MATLAB |
| Modal-shape field normalisation | Same as MATLAB: by the positive maximum only, so a negative-dominant mode can exceed `[-1, 1]`; the Julia plot sets its colour range from the data (§6) |
| — | Julia only: `RitzWorkspace`, `ritzSAWPrecompute!`, `ritzSensitivities!`, `forEachElementPooled`, `with_blas_threads`, `with_task_blas_single`, `blasBackend()` |

## 6. Plots and viewers

Plots and viewers are a Makie package extension: they exist only after `using GLMakie` (an
interactive window) or `using CairoMakie` (files, no display) in the same session; before that a
call throws a plain `MethodError` (`pitfallsJ.md` jOTH4).

### 6.1 Plot functions

| Topic | MATLAB | Julia |
|---|---|---|
| Return values | Graphics handles (`[hPoints, hAx, hBg, hLim]`, line handles `hOut`, ...) | Makie objects: an `Axis` for most plots; `(Axis, Colorbar)` for `plotGreenMap` and `plotRitzSpectrumMap`; a `Figure` for the 3D surfaces; `(Axis, Axis)` for `plotSensitivitySpectrum`; a `PolarAxis` for `plotAngularMap(..., :polar)`. Individual plot objects are in `ax.scene.plots` |
| Drawing into existing axes | `hAx` argument | `ax` (or `fig`) keyword |
| Colour bars | `colorbar(hAx)` beside any axes | Added only when the function creates its own figure; with your own `ax`/`fig` no colour bar is drawn (the diagnostic plots), since Makie needs a layout cell to place one |
| Colormap of the 3D surfaces | `parula` (MATLAB default) | `:viridis` |
| Modal-shape colour range | `clim([-1, 1])`, which clips a negative-dominant mode | Symmetric about zero and set from the data (largest absolute value), with a colour bar; one `Axis` per direction, the layers as stacked bands |
| `plotSensitivitySpectrum` heat map | Row 1 at the top (`imagesc`) | Row 1 at the bottom (Makie's y axis increases upwards); cell contents and labels are the same |
| Map orientation for your own heat maps | `pcolor(x, y, C)`: `C[i, j] ↔ (x[j], y[i])` | Data keep MATLAB's orientation; Makie's `heatmap!(x, y, z)` expects `z[i, j] ↔ (x[i], y[j])` — transpose (`pitfallsJ.md` jNUM3) |

### 6.2 Viewers

The Julia viewers are **experimental** (`quickStartJ.md`, "The GUI layer"); the MATLAB ones are the reference for interactive work.

| Topic | MATLAB | Julia |
|---|---|---|
| Names | One launcher per viewer and family: `Gui.angSweepViewer`, `Gui.angSweepViewerPlate`, `Gui.isofreqAngSweepViewer`, `Gui.freqAngSweep2LayeredViewer`, ... | One launcher per viewer, choosing the family from the result type: `structureControl`, `singleAngleViewer`, `angSweepViewer`, `modalShapeViewer`, `combinedViewer`, `doubleSweepViewer`, `dispersionViewer` |
| Embedding | `parentContainer` | `fig` keyword: `nothing` for a new window, or a `GridPosition` |
| `show()` | Builds hidden, `show()` reveals | No `show()`: a Makie `Figure` is never built hidden |
| Modal-shape viewer windows | `windowMode` `'combined'` (one window) or `'split'` | Always two `Figure`s (controls, shapes); no `windowMode` |
| `combinedViewer` windows | `windowMode` `'combined'`, `'split'`, `'fullSplit'` | `:combined` (default, two windows) or `:split` (three windows; MATLAB's `'fullSplit'` folded in) |
| Modal-shape viewer QC panel | SAW/limiting mode list with "Mark as ..." overrides; angle stepping follows the physical wave (`followPhysicalWave`) | Not ported: angle stepping keeps the raw mode index. Use `singleAngleViewer` on the same result to see which modes are classified |
| Dropdowns | Plain dropdowns | `<<`/`>>` stepper buttons on every dropdown |
| Controls while drilled in | Graph-type/criterion dropdowns disabled | Never disabled (Makie has no disabled state); a change takes effect when you return to the map |
| Criterion dropdown (single angle) | Shown only for the Ritz-mode curve | Always visible; it acts only on the Ritz-mode curve |
| Dispersion viewer | Axis-range labels follow the axis names; checkboxes greyed out when not applicable | Static "X/Y min/max" labels; toggles never greyed out (the effective setting is recomputed on every redraw) |
| Crosshair from the shape viewer into a single-angle curve | Snapped onto the displayed curve | Placed at the correct velocity (or k∥) with y fixed at 0 |
| Kept from MATLAB as they are | isoFreq `:polar` map: a click is resolved in (angle, k∥) although the map is drawn in Cartesian coordinates; velocity-range fallback 0–10 km/s when a sweep has no Green axis | Same |
| Structure editor | One tab per layer, radio-button groups, Cij fields enabled per symmetry class, "open viewer" / "assign to workspace" | One shared layer editor with a layer menu, menus instead of radio groups, all 21 Cij fields always editable (a value typed into a cell the symmetry class does not use is ignored, as the disabled cell is in MATLAB); the result is returned as an `Observable` (`win.result[]`), the hand-off point a viewer can watch |
| Axis ranges from a script | `viewer.setAxisRange(axisName, lims)` | `setAxisRange!(viewer, axis, lims)`, same arguments; ranges kept across slice steps and mode switches of `doubleSweepViewer`, as in MATLAB |

## 7. Status of every MATLAB pitfall in Julia

`docs/pitfalls.md` has 30 items (PHY1–PHY6, MTH1–MTH13, NUM1–NUM6, OTH1–OTH5). "Applies" means the MATLAB item holds for GuliwerrJ as
written (read it there); "differs" says what is different; "not applicable" says why.

| ID | MATLAB pitfall | Julia status |
|---|---|---|
| PHY1 | Existence is a threshold on a continuum | **Applies.** Same scores and standards, with `Symbol` arguments (`:high`, `:low`); every candidate, accepted or not, is in `greenScorePeaks` of the single-angle result. |
| PHY2 | Layered stacks: several branches coexist, and a substrate's second shear limit can look like a SAW | **Applies.** Same layered classification as MATLAB: every accepted wave has the single type `:SAW`, with no Rayleigh/Sezawa split. |
| PHY3 | A leaky pseudo-SAW may have no credible Ritz mode | **Applies.** |
| PHY4 | A tall Green peak can be the crest of a broad shoulder | **Applies.** `GreenPeaks.find` is the Julia peak helper. |
| PHY5 | A true Rayleigh peak's height and width depend on η | **Applies.** Same η default as MATLAB (`SawSpectrum.ETA = 3e-3` in the `Calc` classes; `GreenSolver` default 3e-4). |
| PHY6 | Layered classification: `:low` by default, Love waves need `:high` | **Applies.** `classify(sweep, :green, :high)` for Love waves; the layered Green axis has a fixed 4 m/s step. |
| MTH1 | `legOrder >= 2*relDepth`, tiers | **Applies.** Same tiers and defaults (20/50 easy, 30/70, 40/100). The Ritz-only half-space classification needs `legOrder >= 2.3 * relDepth` and warns once per session otherwise (MTH9); `CalcParams(relDepth = 20, legOrder = 50)` satisfies both. |
| MTH2 | Ritz's lowest mode can sit a hair on the wrong side of `minShear` | **Applies.** `RAYLEIGH_LABEL_TOL` is exported; use `greenPeakVel`, not a Ritz velocity. |
| MTH3 | Decide branch membership from one quantity | **Applies.** |
| MTH4 | A peak can be invisible to `findpeaks` if its window clips it | **Applies** to any peak finder that needs a lower neighbour on both sides, including GuliwerrJ's MATLAB-compatible `findpeaks_matlab`. |
| MTH5 | Existence is decided by the score alone | **Applies.** `vel` is the deciding method's velocity; `modeVel` is the bound Ritz mode's eigenvelocity (use it for synthetic inversion data). |
| MTH6 | Numeric thresholds are per method | **Applies.** `classify(sweep, :green, 0.3)` reads the number on the `:green` scale; a numeric `sawThreshold` needs `sawClassMethod` set (`...:numericThresholdNeedsMethod`). Prefer `:high`/`:low`. |
| MTH7 | A score is only comparable within its own method | **Applies.** |
| MTH8 | The `Calc` classes' Green data are at η = 3e-3 on a coarse axis | **Applies.** Same η (`SawSpectrum.ETA = 3e-3`), same `greenVelStep` default (6 m/s on a half-space). |
| MTH9 | The `:ritz` method is a cheaper approximation of `:green` | **Applies.** The half-space rule needs `legOrder >= 2.3 * relDepth` (`RitzScore.checkBasisHS` warns once per session); layered `:ritz` reports trapped waves only. |
| MTH10 | The `:refined` method models the probe spot | **Applies.** Half-space `:full` only; same limits. |
| MTH11 | A per-point SAW pick can jump between branches in a layered sweep | **Applies.** Use `DispersionBranches.track` (`retentionLevel = :full`) or `trackSeededBranch(sweep)` to follow one continuous branch. |
| MTH12 | `classify`/`reclassify` switch without recomputation | **Applies**, with Symbols; `sawClassMethod` and `sawThreshold` default to `nothing` (not `''`/`[]`). `:ritzOnly` results support `:ritz` only; `:saw_only` retention cannot keep eigenvectors it never saw. |
| MTH13 | A stale result cannot be reclassified | **Applies in form only.** `isStale` and the `...:stale` refusal exist, but without save/load a Julia result is always current. |
| NUM1 | `'` conjugates, `.'` doesn't | **Differs in form.** The same distinction exists (`'` is the adjoint, `transpose(...)` the plain transpose), but the single invisible punctuation mark cannot cause it: `transpose` is a spelled-out call and `A.'` is a syntax error in Julia. Use `transpose(A)` for a complex quantity unless conjugation is intended. |
| NUM2 | `ElastConsts` without `symmNo` guesses the symmetry class | **Applies.** `ElastConsts(consts)` infers the class the same way; pass `symmNo` for Voigt or tensor input (`SymmClassLibrary` lists the numbers). |
| NUM3 | `pcolor` cell-centring | **Differs.** Makie centres Cartesian maps itself; polar maps still use `pcolorCenteredGrid`; matrix orientation differs from Makie's (`pitfallsJ.md` jNUM3). |
| NUM4 | Eigenvector retention RAM cost | **Applies.** A `[nCompModes]`-sized eigenvector `Matrix` costs the same. `retentionLevel = :none` for any non-final sweep step is still the right call; the Plate family has only `:full`/`:none`. |
| NUM5 | `useParallel`, thread pool, default `"auto"` | **Differs in mechanism.** Same purpose and threshold (`AUTO_THRESHOLD = 20` angles), but Julia threads are a start-up flag (`pitfallsJ.md` jNUM5, jNUM1 above). `useParallel = :never` for debugging still applies. |
| NUM6 | `parfor` over single-angle objects: pass the Ritz precompute in | **Differs.** Julia threads share memory, so there is no worker copy to rebuild. In a loop of your own over single-angle objects, build the precompute once (`ritzPrecompute(sweep)`) and pass it as `precompStruct`; otherwise each `singleHS` builds its own, which costs more than a Ritz solve. |
| OTH1 | Two coordinate conventions (`dirMat` vs standard solver frame) | **Applies.** `GetStandardDirMat` is still the one place the two meet; `Euler2DirMat`/`DirMat2Euler` reorder internally, `Euler2DCM`/`DCM2Euler` do not. |
| OTH2 | Never `addpath(genpath(...))` | **Not applicable.** Julia's package system has no path-order-dependent class lookup. |
| OTH3 | A saved result is tied to the class definitions | **Not applicable.** `Calc.saveResult`/`Calc.loadResult` exist only in MATLAB; Julia results live in memory for the session. |
| OTH4 | `Inv.association` needs every measured angle to coincide with a sweep angle | **Applies.** `association(sweep, measured, velocityTolerance, angleTolerance)`, default `angleTolerance = 1e-6`; otherwise `ArgumentError` with `ModeAssociation:fromSweep:angleNotFound`. |
| OTH5 | `plotAngularMap(…, 'polar')` needs polar axes; the Green-map plots need plain ones | **Applies.** `plotAngularMap(sweep, :polar)` draws on (and, given `ax`, expects) a Makie `PolarAxis`; `plotGreenMap`/`plotDispersionMap` with `:polar` draw on a plain `Axis`. Do not reuse one axis between the two kinds. |
