# GuliwerrJ

Julia version of GuliweRR — a guided-wave (Rayleigh/Lamb/pseudo-surface-acoustic-wave)
forward and inverse elasticity solver, also available in MATLAB (`guliwerr/` at the
repository root).

**Read the header comments of the source files before touching them** — they record
empirically verified decisions (eigensolver sorting/phase conventions, the null-space
non-reproducibility, the threading model, six MATLAB-builtin compatibility shims, ...) that are
easy to silently regress if re-derived from scratch. Source comments cite the port's internal
design notes, which are not distributed with this release.

## Status

GuliwerrJ covers the whole of GuliweRR's numerics: the core types (`ElastConsts`, `Material`,
`CalcParams`, `Structure`), the Christoffel, Ritz, and Green solvers, the half-space, layered,
and plate forward calculations (single angles, angular sweeps, and the 2D sweeps), the
surface-wave classification, spot-averaged Green maps, and the inverse characterisation.
Plotting and the GUI are a Makie package extension: one plot function for every MATLAB
`Vis.plot*`, plus a structure editor and six viewers. Results agree with the MATLAB version to
round-off (checked against golden values dumped from MATLAB). The
classification, spot-averaged Green maps, and the inverse module are experimental, as in
MATLAB (`juliaPort/docsJ/functionReferenceJ.md`). `Calc.saveResult`, `Calc.loadResult`,
`Gui.viewSaved`, and `Gui.viewResult` exist only in MATLAB.

`Project.toml`'s `version` (`0.1.0`) tracks the packaging only; `guliwerrVersion()` reports the
GuliweRR version, the same string MATLAB reports.

New to Julia: `juliaPort/docsJ/gettingStartedJ.md`. Coming from MATLAB GuliweRR:
`juliaPort/docsJ/matlabVsJuliaJ.md`. See `juliaPort/docsJ/quickStartJ.md` for how to
actually USE this package (environment setup, a minimal example, the full GUI-launch table) —
this README stays focused on developing the package itself. The comments in the source files
record each module's scope, every deliberate divergence from MATLAB, and the Julia-specific
techniques used.

## Quick start (developing this package)

On a fresh clone, download and precompile the dependencies once (the package has no committed
`Manifest.toml`, so this resolves and installs them):

```sh
cd juliaPort/guliwerrJ
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then:

```sh
julia --project=. -t auto
```

```julia
julia> using GuliwerrJ

julia> Threads.nthreads()   # must be > 1 -- see "Threading" below
8

julia> angSweepHS(; structure = Structure(; materIn = Material(; elConsts = ElastConsts([108.0, 26.4], 1), rho = 2.70),
                                             calcParamsIn = CalcParams(; relDepth = 20, legOrder = 50)),
                    angles = 0.0, calcMode = :ritzOnly, retentionLevel = :saw_only) |> calc |> r -> r.sawVelMap[1][1]
2.9182586815168117   # SAW velocity, km/s -- see juliaPort/docsJ/quickStartJ.md for the full worked example
```

## Threading — read this before running anything sweep-sized

Julia's thread pool is a **process-start flag**, not something acquired at runtime the way
MATLAB's `parpool` is. **Always start Julia with `-t auto` (or `-t N`)** — a session
started without it silently runs everything single-threaded, correctly but slowly, which
is a much harder bug to notice than a crash. `Threads.nthreads() > 1` is worth checking
once per session.

The corresponding hazard: Julia does **not** stop you from writing the same shared, mutable
field from two threads at once (unlike MATLAB's `parfor`, which refuses to compile such a
loop). See `src/core/ParallelPolicy.jl` and `juliaPort/docsJ/pitfallsJ.md` before writing
anything that parallelises over a sweep axis. For batches of many independent sweeps, worker
processes scale better than threads: `juliaPort/docsJ/largeCalculationsJ.md`.

**BLAS backend.** On macOS 13.4+ the package points libblastrampoline at Apple
Accelerate when it loads (`src/core/BlasBackend.jl`); elsewhere it stays on OpenBLAS.
`blasBackend()` reports which; `GULIWERRJ_BLAS=openblas` opts out. `BLAS.set_num_threads`
does not reach Accelerate, whose single-thread switch is per OS thread, so threaded loops
pin it per task with `with_task_blas_single` inside `with_blas_threads(1)`. Results differ
from OpenBLAS at round-off level (`juliaPort/docsJ/matlabVsJuliaJ.md` §3.4).

## Package layout

```
juliaPort/guliwerrJ/
  Project.toml                    [compat] pins version bounds; Manifest.toml is not committed
  src/
    GuliwerrJ.jl              top-level module: includes + re-exports only
    core/         module CoreTypes    -- ElastConsts, Material, CalcParams, Structure, ...
    solvers/      module Solvers      -- Christoffel, Ritz (SAW+Lamb), Green, two levels
                                         below GuliwerrJ (GuliwerrJ.Solvers.Christoffel, ...)
    forwardCalcs/ module ForwardCalcs -- CalcSingleHS/Layered/Plate, angle/2D sweeps, Calc.jl
                                         launcher aliases
    inverseCalcs/ module InverseCalcs -- MeasuredVelocities, ModeAssociation,
                                         SensitivityAnalysis, ElasticityInversion
    visual/       module Visual       -- VisData/PlotData (no Makie), PlotStubs.jl +
                                         GuiStubs.jl (empty stubs the Makie extension below
                                         adds methods to)
  ext/
    GuliwerrJMakieExt.jl      package extension, activates on `using Makie` (any backend)
    Plots*.jl                 plotting (surface/curve, maps, modal-shape, diagnostic)
    Gui*.jl                   GUI (Widgets, StructureControl, SingleAngleViewer,
                               AngSweepViewer, ModalShapeViewer, CombinedViewer,
                               DoubleSweepViewer, DispersionViewer)
```

`Makie`/`GLMakie` are declared under `[weakdeps]` only; `[extensions] GuliwerrJMakieExt =
"Makie"` triggers on `Makie` alone (backend-agnostic) — `using GuliwerrJ, GLMakie` (desktop),
`using GuliwerrJ, CairoMakie` (headless/testing), or in principle `using GuliwerrJ, WGLMakie`
(browser-hosted) activate the same extension code with zero package changes; only `GLMakie`
and `CairoMakie` are tested.

### Extending the package: one module per top-level folder, not one per file

If you add a new file under `src/`, do not wrap it in its own `module X` unless `X` differs
from every type/struct name it defines. A file `module CalcParams` containing
`struct CalcParams` breaks: once the parent module's `include("CalcParams.jl")` runs, the
parent scope has already bound the name `CalcParams` to the *module object*, and a subsequent
`using .CalcParams` cannot rebind it to the struct — bare `CalcParams(...)` then fails with
`MethodError: objects of type Module are not callable`. The convention is therefore one
`module` per top-level folder (`CoreTypes`, `Solvers`, `ForwardCalcs`, `InverseCalcs`,
`Visual`), with individual files `include`d as plain top-level definitions inside it. A new
file belongs inside one of these existing modules' own `include` list, not in a fresh module
of its own.

## Licence

ISC, as the rest of GuliweRR: see `LICENSE` at the repository root.
