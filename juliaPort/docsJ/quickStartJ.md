# GuliwerrJ — quick start

*New to Julia? Start with `gettingStartedJ.md` (install, environment, first plot).*
*Coming from MATLAB GuliweRR? Read `matlabVsJuliaJ.md` — every difference, by topic.*

GuliwerrJ is the Julia version of GuliweRR. Read this, `functionReferenceJ.md`, and
`pitfallsJ.md` before using it. If you want to modify the package's own internals instead,
read `juliaPort/guliwerrJ/README.md` (layout, threading) and the header comments of the
source files — they record the design decisions (eigensolver conventions, threading model,
module layout) this package follows.

**The physics is the same.** GuliwerrJ computes the same Ritz/Christoffel/Green pipeline as
the MATLAB version, with the same conventions — see `docs/background.md` for what the toolbox
is for and the maths behind it; nothing here duplicates that. Only the code is in a different
language.

## Environment setup

Julia has no `run('init.m')` equivalent — there is no global path state to set up, only a
per-project *environment* (a `Project.toml`/`Manifest.toml` pair pinning exact dependency
versions) that you activate with the `--project=` flag. Two ways to get `GuliwerrJ` usable:

**Working inside the package itself** (extending GuliwerrJ). On a fresh
clone, download and precompile the dependencies once:

```
julia --project=juliaPort/guliwerrJ -e 'using Pkg; Pkg.instantiate()'
```

then start Julia in the package environment:

```
julia --project=juliaPort/guliwerrJ -t auto
```

and run `using GuliwerrJ` at the REPL. `-t auto` starts Julia with a thread pool sized to your
machine — **do this every time**; see `pitfallsJ.md` jNUM5 for what happens if you forget it.

**Using GuliwerrJ from somewhere else** (a calculation script that lives outside the package,
the equivalent of a MATLAB script that just needs `guliwerr/` on the path): add it to your own
project as a *development dependency*, which points Julia at the local source tree instead of
a registry:

```julia
using Pkg
Pkg.develop(path = "/absolute/path/to/juliaPort/guliwerrJ")
Pkg.instantiate()
using GuliwerrJ
```

Run once per project environment (it edits that environment's own `Project.toml`); after
that, `using GuliwerrJ` works from any script started with
`julia --project=<that environment> -t auto`.

**BLAS/LAPACK backend.** On macOS (13.4 or later) `using GuliwerrJ` switches Julia's linear
algebra from OpenBLAS to Apple's Accelerate framework — a Ritz step ran 1.4–2.7× faster
on the Apple-silicon machine it was measured on. Every other platform keeps OpenBLAS. `blasBackend()` reports which one is active
(`:accelerate` or `:openblas`). To compare against OpenBLAS on a Mac, set the environment
variable before Julia loads the package:

```
GULIWERRJ_BLAS=openblas julia --project=juliaPort/guliwerrJ -t auto
```

The two backends agree to round-off, not bitwise (eigenvalues ~1e-11 relative), so a result
computed on a Mac and one computed on Linux, or with the opt-out, differ in the last digits.
See `matlabVsJuliaJ.md` §3.4.

## Minimal working example — half-space SAW velocity

The direct Julia counterpart of `docs/quickStart.md`'s own example: build a `Material`, wrap
it in a `Structure`, run a one-angle `angSweepHS` sweep. Every constructor except
`ElastConsts` and the inverse launchers is **keyword-only** — see `pitfallsJ.md` jOTH1 before
reaching for a positional call out of MATLAB habit.

```julia
using GuliwerrJ

# Material -- isotropic Al (isotropic class takes [c11, c44] only)
ec  = ElastConsts([108.0, 26.4], 1)          # [C11, C44] GPa, symmNo=1 (isotropic)
mat = Material(elConsts = ec, rho = 2.70)    # rho, g/cm^3
cp  = CalcParams(relDepth = 20, legOrder = 50)   # "easy" tier, docs/pitfalls.md MTH1

str = Structure(materIn = mat, calcParamsIn = cp, name = "al_isotropic", nMaxModes = 500)

# Single-angle sweep (angles can be a 1-element vector)
sweep = angSweepHS(structure = str, angles = [0.0], calcMode = :ritzOnly, retentionLevel = :saw_only)
sweep = calc(sweep)

vSAW = sweep.sawVelMap[1][1]   # km/s
println("v_SAW = ", round(vSAW, digits = 4), " km/s")
```

This prints `v_SAW = 2.9183 km/s`. As in MATLAB, `angSweepHS` (not the lower-level
`singleHS`) is the right entry point even for one angle: it builds the Ritz precompute once
and shares it between angles. `singleHS` builds its own precompute, which is why a sweep is
cheaper; it is mainly the building block `angSweepHS` uses per angle. Note the Symbol
arguments (`:ritzOnly`, `:saw_only`) where MATLAB used quoted strings (`'ritzOnly'`,
`'saw_only'`) — see `pitfallsJ.md` jOTH2.

## Full forward example — angular sweep + visualisation

`calcMode = :full` also runs the Green's-function solver, exactly as in MATLAB. The material
below switches to anisotropic Cu (rather than the isotropic Al above) for the same reason the
MATLAB quick-start does: an isotropic material's SAW velocity has no angular dependence, so
the angular Green map would be flat.

Plotting is a **package extension**: nothing renders until a concrete Makie backend is loaded
alongside GuliwerrJ. `using CairoMakie` (headless, software-only) or `using GLMakie` (opens a
real window) both work — either one loads `Makie` itself as a transitive dependency, which is
all the extension checks for.

```julia
using GuliwerrJ
using CairoMakie   # or GLMakie for an interactive window

ec  = ElastConsts([170.0, 122.0, 74.0], 2)   # [C11, C12, C44] GPa, symmNo=2 (cubic)
mat = Material(elConsts = ec, rho = 8.96)     # rho, g/cm^3; default dirMat = Cu(001)
cp  = CalcParams(relDepth = 20, legOrder = 50)   # starting point for Cu(001), see docs/pitfalls.md MTH1
str = Structure(materIn = mat, calcParamsIn = cp, name = "cu_001", nMaxModes = 500)

angles = collect(0.0:2.0:90.0)
sweep  = angSweepHS(structure = str, angles = angles, calcMode = :full,
                     retentionLevel = :full)
sweep  = calc(sweep)

ax1, cb1 = plotGreenMap(sweep)                                       # |G13_eff| map
ax2      = plotClassificationDiagnostic(sweep, :type, true, :full)   # classified points over the Green background
```

This runs the 46-angle Ritz and Green sweep in well under a minute (the first run of a
session also compiles the code it needs).

### SAW classification: method and standard

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface
> waves in strongly anisotropic materials. Thresholds, scores and the interface may change.

Every result carries one active classification (`sweep.activeMethod`,
`sweep.activeThreshold`). `calc` runs `:green` (the Green existence score on |G13_eff|) on
`calcMode = :full` and `:ritz` (the Ritz-only score) on `:ritzOnly`. The default standard is
per method: `:high` for `:green` and `:refined` on a half-space, `:low` for `:ritz` and for
every layered stack. Switching reuses the stored spectra — no new Green or Ritz solve
(`:refined`, the spot-averaged half-space score, solves its spot directions once):

```julia
classify(sweep, :ritz, :high)     # Ritz-only score, stricter standard
classify(sweep, :green, :low)     # Green score, looser standard
classify(sweep, :green, 0.3)      # a numeric threshold on the method's own score scale
sm = sweep.cs[10].classification.significantModes   # per-angle entries
sm[1].type, sm[1].method, sm[1].scoreDb, sm[1].grade  # e.g. :Rayleigh, :green, 36.8, :strong

absG, label, v = greenComponentData(sweep, "G33")    # any other Green component: recomputed
plotGreenMap(sweep; component = "G33")
plotRitzSpectrumMap(sweep)        # the Ritz spectrum |S(v)| from the stored numerators
plotRitzModeMap(sweep)            # every Ritz mode coloured by its Ritz score [dB]
plotModalShape(sweep.cs[10].ritzSolver, 1, :none, nothing, :rgb)   # polarisation-RGB view
```

Methods: `:green` | `:refined` (half-space only) | `:ritz`; standards `:high` | `:low` | a
number. Only |G13_eff| is stored (`greenAbsG13Eff` / `greenAbsG13EffMap`, `Float32`); any
other component is recomputed on request. For a working, runnable version of this pattern —
and the layered/isoFreq equivalents — copy from `juliaPort/tutorialsJ/` rather than retyping.

**Same as MATLAB.** The classification is one-to-one with `docs/quickStart.md`, "Choosing the
classification", with `Symbol`s for strings. `sawClassMethod`, `sawThreshold` and
`greenVelStep` are fields of the (mutable) result and can be set before `calc` (unset means
`nothing`, not `''`/`[]`). `reclassify(sweep, :high)` keeps the method, `isStale(sweep)`
exists, and `refinedSawAgreement(sweep)` pairs the `:green` and `:refined` entries after
`classify(sweep, :refined)`. On a layered stack (`singleLayered`, `angSweepLayered`,
`wavevecAngSweep2Layered`) the methods are `:green` and `:ritz`, both with default `:low`:
Rayleigh- and Sezawa-type waves. Love waves need `reclassify(sweep, :high)`; leaky waves need
`:green` (`docs/pitfalls.md` PHY6 and MTH9).

### Spot-averaged Green maps (EXPERIMENTAL)

`calcSpotGreenMap` / `calcSpotGreen` draw |G| as a measurement with a finite pattern sees it:
the complex Green's function summed coherently over the wavevectors a `ProbeKernel` couples
to (a steep branch broadens and drops, a flat one keeps its height). It is a picture only —
no classification reads it — and, as in MATLAB, it is validated only against the published Fe3Al(316) map (UTGS paper, arXiv:2510.10696, Figs. 3–4).
Half-space only, any `calcMode`; each result is stored on its own velocity axis:

```julia
k  = probe_gaussian(0.004, 0.3)               # sigma_|k|/k0, sigma_angle [deg]; also probe_tgs, probe_sras
sw = calcSpotGreenMap(sweep, k)               # DirStep, RelKStep, Eta, VelStep keywords; stores spotGreenMap
plotGreenMap(sw; component = "spot")          # its own spotGreenVelAxis
cs = calcSpotGreen(sweep.cs[10], k)           # one angle; equals column 10 of the map bit for bit
plotGreenCurve(cs; component = "spot")
sw.spotGreenInfo.dirConvergence               # per-angle quadrature check ("comb"); > 5 % warns
```

An SRAS (spatially resolved acoustic spectroscopy) kernel, `probe_sras`, is one-sided:
compute 0–359°, not 0–180°.

## The GUI layer

> **Experimental** (the interactive viewers only — the `plot*` functions above are not affected). The Julia viewers are tested headlessly (wiring, redraws, layout sizes) but have had little use in a real window and are still awkward to use. For interactive work use the MATLAB viewers (`Gui.viewSaved`, `docs/guiGuide.md`); use these when a calculation already lives in Julia and a quick look is enough.

GuliwerrJ ships six interactive viewers plus a structure-editor/job-launcher, all built on
Makie. Each is exported from `GuliwerrJ` itself (as an empty stub, `src/visual/GuiStubs.jl`,
gaining its real implementation from the extension — the same pattern the `plot*` functions
above use), so they are callable directly once Makie is loaded:

```julia
using GuliwerrJ
using GLMakie   # or CairoMakie for a headless smoke-test

win = structureControl()          # build a Structure and launch a calculation, GUI-driven
```

| Launcher | Opens | Needs |
|---|---|---|
| `structureControl()` | Structure editor + job launcher | nothing — builds a `Structure` interactively |
| `singleAngleViewer(result)` | Single-angle curve view (classification / Green / Ritz-mode / Ritz spectrum) | a `calcAll()`'d `CalcSingleHS`/`CalcSingleLayered`/`CalcSinglePlate`, or a `calc()`'d `CalcIsofreqSingleLayered` |
| `angSweepViewer(sweepResult)` | Angular-sweep 2D map | a `calc()`'d `CalcAngSweepHS`/`CalcAngSweepLayered`/`CalcAngSweepPlate`/`CalcIsofreqAngSweepLayered` |
| `modalShapeViewer(result)` | Modal-shape depth-panel view | any of the above single or sweep results |
| `combinedViewer(result)` | Two-pane shell (shape + map/curve, wired together) | any single or sweep result |
| `doubleSweepViewer(sweepResult)` | Outer-axis slice stepper, with a dispersion-mode switch for the Layered/Plate cases | a genuine double sweep (`CalcWavevecAngSweep2Layered`/`...Plate`/`CalcFreqAngSweep2Layered`, `calc()`'d) |
| `dispersionViewer(sweepResult)` | Modal shapes + dispersion-curve map at one fixed angle | a `CalcWavevecAngSweep2Layered`/`...Plate` result with `retentionLevel = :full` |

Every launcher accepts a `fig` keyword (`nothing` by default, opening a standalone window; or
a `GridPosition` to embed into a figure you already own). `GLMakie` is what actually opens an
interactive window on your desktop; `CairoMakie` renders to an in-memory/file target only and
suits machines without a display. Neither is a dependency of GuliwerrJ: install the backend in
an environment on your load path, e.g. your default one (`] activate @v1.13`, `] add GLMakie`),
so that `using GLMakie` works next to `--project=juliaPort/guliwerrJ`.

**Running GuliwerrJ on a remote machine.** The extension triggers on `Makie` alone, not a
specific backend (`Project.toml`'s `[extensions]` table names `Makie`, not `GLMakie`), so
other Makie backends such as `WGLMakie` load the same GUI code. Only `GLMakie` and
`CairoMakie` are tested.

## MATLAB-only features

`Calc.saveResult`, `Calc.loadResult`, `Gui.viewSaved` and `Gui.viewResult` exist only in
MATLAB. A GuliwerrJ result lives in memory for the session; to inspect one interactively, use
the viewers above.

## Where to go next

- **`functionReferenceJ.md`** — one-line-per-function index of every entry point, grouped by
  category, built from the export list and source signatures.
- **`pitfallsJ.md`** — traps specific to the Julia version: mistakes that fail or silently do
  the wrong thing.
- **`matlabVsJuliaJ.md`** — everything that differs from MATLAB GuliweRR (calling
  conventions, errors, numerics, threading, what is not ported, plots and viewers), and which
  items of `docs/pitfalls.md` apply unchanged, differ, or do not apply.
- **`largeCalculationsJ.md`** — batches of many independent sweeps (orientations, materials,
  parameter grids): threads or worker processes, how many workers, restartable runs.
- **`docs/background.md`** — what GuliweRR is and the physics behind the Ritz/Christoffel/Green
  pipeline (identical for both languages).
- **`juliaPort/tutorialsJ/`** — narrated, copy-paste-ready scripts, one per calculation type,
  mirroring `examples/tutorials/`'s numbered scripts one for one (01–09, plus `S1` material
  input and `S2` free-plate Lamb waves), with the Julia-specific concepts explained inline.
- **`juliaPort/tutorialsCodeJ/`** — the same tutorials, same file names, stripped to code:
  one-line section headers and comments, every argument a named keyword/variable so each call
  can be copied directly (the Julia counterpart of `examples/tutorialsCode/`).
- **Recovering elastic constants from measured velocities?** That is the *inverse* problem —
  see `functionReferenceJ.md`'s "Inverse characterisation" section (`measured`/`association`/
  `sensitivity`/`elasticity`).
- Any function's own Julia docstring (`?functionName` at the REPL, or read the source
  directly) — denser than this file, with the reasoning behind any deliberate divergence from
  MATLAB.
