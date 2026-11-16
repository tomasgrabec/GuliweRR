# GuliweRR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24280695.svg)](https://doi.org/10.5281/zenodo.24280695)

*Guided Linear Waves by Ritz–Rayleigh* — a MATLAB toolbox, with a one-to-one Julia version, that calculates elastic-wave propagation in anisotropic, layered media: surface-guided waves (surface acoustic waves, SAWs) on a half-space or film(s)-on-substrate stack, and Lamb-wave dispersion on a free (multi-layered) plate.

## What it does

For a stack of materials (elastic constants, density, layer thicknesses, crystal orientation) and a propagation direction, the forward calculations give:

- bulk-wave velocities and limiting velocities of any anisotropic material (Christoffel equation);
- SAW and pseudo-SAW velocities, mode identities, and modal shapes on a half-space or a film-on-substrate stack (Ritz–Rayleigh method);
- the surface elastodynamic Green's function of the same structures, in which surface waves show up as sharp peaks (frequency-domain solver after Every et al.);
- Lamb-wave dispersion of free-standing homogeneous and multi-layered plates, including zero-group-velocity points;
- surface-wave classification (which spectral peaks are surface waves) and spot-averaged Green maps (*experimental*);
- plotting functions and interactive GUI windows for all of the above.

An inverse procedure recovers the elastic constants from measured velocities (*experimental*).

> **Experimental parts:** the inverse procedure, surface-wave classification, and spot-averaged Green maps work and are tested, but their results and interfaces may change. Each is labelled again in its own documentation.

## Installation (MATLAB)

1. Download or clone this repository.
2. Run the initialisation script once per MATLAB session (from the project root):
   ```matlab
   run('guliwerr/init.m')
   ```
3. To load GuliweRR automatically at every MATLAB startup, add that line to your `startup.m` file:
   ```matlab
   % In MATLAB: edit(fullfile(userpath,'startup.m'))
   run('/full/path/to/Guliwerr/guliwerr/init.m')
   ```

> **Note:** Do not use `addpath(genpath(...))` on the `guliwerr/` folder — `genpath` recurses into `@ClassName/` subfolders and breaks MATLAB's class lookup.

## Quick start (MATLAB)

See [`docs/quickStart.md`](docs/quickStart.md): a minimal SAW-velocity calculation, an angular sweep with plots, and the choice of surface-wave classification. Copy-paste starting points are in [`examples/tutorials/`](examples/tutorials/) (narrated) and [`examples/tutorialsCode/`](examples/tutorialsCode/) (the same scripts, code only). For an overview of the whole package, open the Atlas, [`docs/GuliweRR_Atlas_v6.html`](docs/GuliweRR_Atlas_v6.html).

## Julia version

The same calculations are also available as a Julia package, `GuliwerrJ`, in [`juliaPort/`](juliaPort/). It is released together with the MATLAB version, one to one, and gives the same results within numerical tolerance (checked against MATLAB reference values). MATLAB is the user-facing reference version: MATLAB is more widely used, and its interactive viewers are functional, while the Julia viewers are experimental (Julia's plotting functions work fine). Julia is the computational version — about 1.4 to 26 times faster, depending on the calculation (measured on one Apple-silicon machine). Start with [`juliaPort/docsJ/gettingStartedJ.md`](juliaPort/docsJ/gettingStartedJ.md) and [`juliaPort/docsJ/quickStartJ.md`](juliaPort/docsJ/quickStartJ.md); the differences from MATLAB are in the Atlas, section "Two implementations: MATLAB and Julia".

## Documentation

**Start here**

- **[`docs/GuliweRR_Atlas_v6.html`](docs/GuliweRR_Atlas_v6.html)** — navigable, diagrammed overview of the whole package (architecture, pitfalls, changelog, roadmap) in one page
- [`docs/quickStart.md`](docs/quickStart.md) — path setup and minimal/full worked examples (MATLAB)
- [`juliaPort/docsJ/gettingStartedJ.md`](juliaPort/docsJ/gettingStartedJ.md) — the same for Julia, from installing Julia to a first plot; then [`quickStartJ.md`](juliaPort/docsJ/quickStartJ.md)
- `examples/tutorials/` (MATLAB) and `juliaPort/tutorialsJ/` (Julia) — narrated tutorials, numbered in reading order

**Using GuliweRR**

- [`docs/guiGuide.md`](docs/guiGuide.md) — the click-driven alternative: build a structure, run a calculation, and browse the result in the GUI windows
- [`docs/pitfalls.md`](docs/pitfalls.md) — gotchas that have each caused a real bug (coordinate-axis convention, MATLAB complex-transpose, mesh-resolution rule, etc.)
- [`docs/background.md`](docs/background.md) — what GuliweRR is and the physics behind it

**Reference** (`docs/reference/`)

- [`functionReference.md`](docs/reference/functionReference.md) — one-line-per-function index of every entry point
- [`FullFunctionReference.html`](docs/reference/FullFunctionReference.html) — every public class property/method reachable after `run('guliwerr/init.m')`, with its call signature and MATLAB's own `help` text, generated from the source
- [`architecture/`](docs/reference/architecture/) — one page per module (core, Christoffel, Ritz, Green, plates, forward calculations and mode classification, visualisation, GUI): how it works and why
- Julia: [`functionReferenceJ.md`](juliaPort/docsJ/functionReferenceJ.md), [`matlabVsJuliaJ.md`](juliaPort/docsJ/matlabVsJuliaJ.md), [`pitfallsJ.md`](juliaPort/docsJ/pitfallsJ.md), [`largeCalculationsJ.md`](juliaPort/docsJ/largeCalculationsJ.md)

## Requirements

- MATLAB — developed and tested on R2025b only
- Signal Processing Toolbox (`findpeaks`) required — the surface-wave classification of every half-space and layered calculation uses it
- Optimization Toolbox (`lsqnonlin`) — needed only for `inverseCalcs/`
- Parallel Computing Toolbox optional — used automatically for large angular sweeps (`useParallel="auto"`), falls back to serial execution with a warning if unavailable
- Julia version (optional, `juliaPort/`): Julia 1.13 or newer; packages are installed by Julia itself (`juliaPort/docsJ/gettingStartedJ.md`)
- Platform: both versions are developed and tested on macOS on Apple silicon only. There the Julia version runs its linear algebra on Apple's Accelerate library; on Linux and Windows it uses OpenBLAS, which is untested so far, so speed and round-off-level results may differ

## Folder structure

```
docs/               Documentation: the Atlas, quick start, GUI guide, pitfalls and background on top; reference/ holds the function references and the per-module architecture notes
examples/           Narrated tutorials, code-only tutorials, and case studies
juliaPort/          The Julia package GuliwerrJ with its own docs and tutorials
guliwerr/           The MATLAB package (add to the path with guliwerr/init.m)
├── core/           ElastConsts, Material, CalcParams, Structure, ParallelPolicy, ProbeKernel
├── solvers/
│   ├── christoffel/  Bulk wave velocities and limiting velocities (Christoffel equation)
│   ├── green/        Green's function solver (half-space and layered)
│   └── ritz/         Ritz-Rayleigh forward solver (SAW and Lamb)
├── forwardCalcs/   Angle/wavevector/frequency sweep orchestrators (CalcAngSweepHS, CalcAngSweepLayered, ...)
├── inverseCalcs/   Elasticity inversion from measured velocities (MeasuredVelocities, ModeAssociation, SensitivityAnalysis, ElasticityInversion)
├── visual/         Plotting functions
└── gui/            Interactive GUI windows
```

Each class with many static methods has a main classdef file alongside a `@ClassName/` subfolder holding one `.m` file per static method implementation. MATLAB resolves `ClassName.Method()` transparently.

## Validation

The Ritz–Rayleigh results were checked against the published methods and papers listed under Citation, and against the Green's-function solver, which is an independent calculation of the same surface waves. Modes of multi-layer plates agree with GEWtool, open-source guided-wave software by D. A. Kiefer, to about 1e-9 relative, and the zero-group-velocity points of a silicon wafer match published measurements (Kiefer, Mezil & Prada 2023); details in `docs/reference/architecture/plate.md`. The Julia version is checked against reference values computed with MATLAB.

## Roadmap

See [RoadmapPublic.md](docs/RoadmapPublic.md).

## Changelog

See [ChangelogPublic.md](docs/ChangelogPublic.md).

## Versioning

Versions are numbered MAJOR.MINOR.PATCH: MAJOR marks a large restructuring of the code, MINOR new or changed functionality, PATCH fixes only. A minor release can change or remove an existing function or option; every such change is listed under **Breaking** in [ChangelogPublic.md](docs/ChangelogPublic.md), so check that list before updating. Each part is a whole number, so 6.9 is followed by 6.10. `guliwerrVersion()` returns the version you have; a `-dev` suffix means code between releases. Every major and minor release gets its own DOI; patch releases do not (cite the DOI of the minor release they belong to).

## Citation

If you use GuliweRR in your research, please cite the software itself ([doi:10.5281/zenodo.24280695](https://doi.org/10.5281/zenodo.24280695), see [`CITATION.cff`](CITATION.cff)) and the method paper(s) below — Grabec et al. (2024) is the main, most up-to-date reference and should be cited by default; add the others only if specifically relevant, and add the Green's-function paper as well if you used `GreenSolver`:

- **Main reference (default)** — Grabec, T., Soudná, Z., Repček, K., Lünser, K., Fähler, S., Stoklasová, P., Sedlák, P., Seiner, H. (2024). Guided acoustic waves in thin epitaxial films: Experiment and inverse problem solution for NiTi. *Ultrasonics*, 138, 107211. [doi:10.1016/j.ultras.2023.107211](https://doi.org/10.1016/j.ultras.2023.107211)
- **Lamb (plate) waves** — Grabec, T., Sedlák, P., Seiner, H. (2020). Application of the Ritz–Rayleigh method for Lamb waves in extremely anisotropic media. *Wave Motion*, 96, 102567. [doi:10.1016/j.wavemoti.2020.102567](https://doi.org/10.1016/j.wavemoti.2020.102567)
- **Half-space, original paper / inverse-procedure discussion** — Stoklasová, P., Sedlák, P., Seiner, H., Landa, M. (2015). Forward and inverse problems for surface acoustic waves in anisotropic media: A Ritz–Rayleigh method based approach. *Ultrasonics*, 56, 381–389. [doi:10.1016/j.ultras.2014.09.004](https://doi.org/10.1016/j.ultras.2014.09.004)
- **Green's function (`GreenSolver`), additive** — Every, A. G., Kim, K. Y., Maznev, A. A. (1997). The elastodynamic response of a semi-infinite anisotropic solid to sudden surface loading. *Journal of the Acoustical Society of America*, 102(3), 1346–1355. [doi:10.1121/1.420053](https://doi.org/10.1121/1.420053)

## Acknowledgements

The Ritz–Rayleigh method and its extensions implemented here build on theoretical work by Petr Sedlák and Hanuš Seiner, co-authors of the papers above, whose contributions to the underlying theory and earlier (pre-GuliweRR) implementations made this software possible.

`ZgvNewton` (the Newton refinement of zero-group-velocity points) is adapted from `ZGVNewtonBeta` in GEWtool by Daniel A. Kiefer (ISC licence; algorithm by Bor Plestenjak; Kiefer, Plestenjak, Gravenkamp & Prada, *J. Acoust. Soc. Am.* 153, 1386 (2023)); its licence notice is kept in that file.

GuliweRR was developed by Tomáš Grabec with assistance from Claude Code (Anthropic), an AI coding assistant, in implementation, testing and documentation. All physics, method choices and results were directed and checked by the author.

## Licence

GuliweRR is released under the [ISC License](LICENSE) — a permissive licence functionally equivalent to MIT/BSD-2-Clause. It covers the whole repository, the Julia version in `juliaPort/` included.
