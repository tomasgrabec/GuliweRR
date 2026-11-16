# forwardCalcs — architecture and design rationale

`guliwerr/forwardCalcs/`: the half-space (HS), layered, and isoFreq `Calc*` classes below, the
`@Calc` namespace launcher, the shared helper `CalcCommon`, and the classification classes
(`ModeClassification` and the scoring classes). The free-plate `Calc*Plate` classes are in
`plate.md`. This is the orchestration layer: it drives Christoffel, `RitzSolverSAW`, and
`GreenSolver` (see `green.md`) across single points, angle sweeps, and 2D (angle ×
wavelength/frequency) sweeps, and turns their raw output into a classified list of physical
surface acoustic wave (SAW) modes. This layer does not solve the elastodynamic problem itself —
every Green and Ritz solve is delegated to the solver layer; its job is *which combination of
solves to run, in what order, sharing what, and how to interpret the result*.

The classification detail is in the companion file `forwardCalcs_modeClassification.md`. Read this
file first for the class family and the orchestration pattern, that one for how a Green or Ritz
peak becomes a labelled "Rayleigh"/"pseudoRayleigh" (half-space) or "SAW"/"limiting" (layered)
mode.

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface waves in
> strongly anisotropic materials. Thresholds, scores, and the interface may change.

## Class family map

| Class | Geometry | Swept quantity | Drives | Ritz? |
|---|---|---|---|---|
| `CalcSingleHS` | half-space | — (one angle) | Christoffel + Ritz + Green directly | yes |
| `CalcSingleLayered` | layered | — (one angle) | Christoffel (substrate-only) + Ritz + Green directly | yes |
| `CalcAngSweepHS` | half-space | angle | an array of `CalcSingleHS` | via single |
| `CalcAngSweepLayered` | layered | angle | an array of `CalcSingleLayered` | via single |
| `CalcWavevecAngSweep2Layered` | layered | angle × wavelength | an array of `CalcAngSweepLayered`, one per wavelength | via sweep |
| `CalcIsofreqSingleLayered` | layered/HS | k∥ (fixed ω₀) | `GreenSolver.calc` directly, `'ScanMode', 'isoFreq'` | no |
| `CalcIsofreqAngSweepLayered` | layered/HS | angle (fixed ω₀), k∥ inner | an array of `CalcIsofreqSingleLayered` | no |
| `CalcFreqAngSweep2Layered` | layered/HS | angle × ω₀, k∥ inner | an array of `CalcIsofreqAngSweepLayered`, one per ω₀ | no |
| `ModeClassification` | either | — | consumes Ritz/Christoffel/Green results | — |

Every "drives an array of X" class is composition, not inheritance: no shared base class
exists anywhere in this family, the plate family included (see next section). Each level reuses the level below purely by
calling its already-tested public methods (constructor + `calc()`), and adds exactly one new
concern of its own — see each section below for what that concern is.

## Why half-space and layered are separate class families

Not a stylistic choice: layered structures have an arbitrary number of Sezawa modes and no
system-wide Christoffel limiting velocity (a layered stack has no single bulk material to bound
it), so a `CalcSingleHS` generalised with an `nLay` branch would either silently mishandle the
layered case's genuinely different mode taxonomy or grow ad-hoc branches to cope with it. The two
class families reuse code by calling the same solver APIs (`RitzSolverSAW`, `GreenSolver.calc`,
`Christoffel.matrix`), never by inheritance. No common base class exists, deliberately; the free-plate
family (`plate.md`) is built the same way, by composition.

The practical fallout of the split is visible in the properties themselves:

| | `CalcSingleHS` | `CalcSingleLayered` |
|---|---|---|
| Christoffel reference | system-wide (`structure.mater(nLay)`, the whole HS is one material) | substrate-only (`structure.mater(1)`, prefixed `sub*` throughout so it's never mistaken for a system-wide bound) |
| Green axis end | derived: `1.05 * longLimVel` (+ 0.1 km/s end taper) | `max(vMax, 1.05 * bandVLmax)` (+ 0.1 km/s) (`bandVLmax` = the fastest longitudinal limit of any layer, the end of the band the classification searches); the optional `vMax` only extends the displayed axis beyond it |

`Structure`'s own layer-index convention is `mater(1)` = bottom (substrate) … `mater(nLay)` =
top (surface) — `CalcSingleLayered.calcLim()`'s `subIdx = 1` and `CalcSingleHS.calcLim()`'s
`surfIdx = obj.structure.nLay` both follow this, and the comments in each explicitly flag that
using the other class's index convention would be wrong for that class.

## Single-point containers: `CalcSingleHS` / `CalcSingleLayered`

The atomic unit: one `Structure`, one propagation angle, held as a handle object (RAM-heavy
per-object data — Ritz eigenvectors, Green spectrum — and GUI write-back both need reference,
not value, semantics). The constructor rotates the structure once (`CalcCommon.rotateStructure`, shared by every
single-point class) and stores the *rotated* structure — callers never rotate a `Structure` by
hand anywhere in this layer.

`calcAll(mode, retentionLevel)` runs the full pipeline:

```
'full':     calcRitz() -> calcLim() -> calcGreen() -> calcGreenScore() -> calcRitzScore()
            -> identifyModes() ('green', then sawClassMethod if another) -> applyRetention()
'ritzOnly': calcRitz() -> calcLim() -> calcRitzScore() -> identifyModes() ('ritz')
            -> applyRetention()   (no Green)
```

- **`calcGreen(velAxis)`** is the angle's one Green solve (`CalcCommon.greenScoreSpectrum`),
  at `SawSpectrum.ETA` = 3e-3. By default it runs on the angle's own axis (`greenAxis()`, see "The
  Green velocity axis" below), so it needs `calcLim()` first; a sweep passes its common axis
  instead (below). The solve is +k only; by reciprocity `G13_eff = (G13_+ − G31_+)/2` comes from
  it (`green.md`, "±k symmetrisation"). A result stores only `|G13_eff|`: the display magnitude
  `greenAbsG13Eff` in single precision, and the classification input — the complex `G13_eff` on
  the existence score's band — as `greenScoreG13Eff`/`greenScoreVelAxis`, in double. Any other
  component is recomputed on request, never stored: `greenComponent(component)` returns a
  `[nVel × 1]` single spectrum on `greenVelAxis` (one Green solve, same axis and damping).
  `component` is `'GIJ'` for the one-sided `G_IJ(+k)` (`I` = force direction, `J` = displacement
  direction; 1 = propagation, 3 = surface normal), e.g. `'G13'` (= `G13_+`) or `'G33'`, or
  `'GIJeff'` for the ±k-symmetrised one, e.g. `'G13eff'` (returns the stored data, no solve). The
  parsing and the solve live in `CalcCommon.parseGreenComponent` and
  `CalcCommon.greenComponentSpectrum(structure, velAxis, component, fname)`.
- **`calcGreenScore(ref)`** scores this angle's `G13_eff` candidates
  (`forwardCalcs_modeClassification.md`): half-space `GreenScore` (probe slopes, scores, and the
  `'low'` reference `greenScoreRef` — a sweep passes its own; without one, the memoised reference
  of the unrotated structure), layered `GreenScoreLayered` (probes and polarisation; no
  reference). The scored candidates are stored in `greenScorePeaks`. It solves the band itself
  only when `calcGreen` ran on an axis that does not cover it.
- **`calcRitzScore(precomp)`** stores the Ritz candidates: the modal-sum numerators `u3·conj(u1)`
  of the modes up to 1.2 × the (fastest) longitudinal limiting velocity (`ritzScoreNum` — the Ritz
  spectrum), the score axis (`ritzScoreVelAxis`), this direction's own floor (`ritzScoreRef`) and
  the scored candidates (`ritzScoreCands`). It runs on both calculation modes, in ~10–20 ms (machine-dependent), with no other
  direction solved. It needs the eigenvectors, so it runs before pruning.
- **`classify(method, threshold)`** sets the active classification — what `significantModes` (and,
  on a sweep, the `saw*` maps) hold — from the stored candidates, without recomputing Green or
  Ritz. **`reclassify(threshold)`** = `classify(activeMethod, threshold)`. The properties
  `sawClassMethod` / `sawThreshold`, their defaults, `'refined'`'s first-use solve, the limiting
  entries and the stale stamp are in `forwardCalcs_modeClassification.md`, section 8.
- **`calcSpotGreen` / `calcSpotGreenMap(kernel)`** (experimental, half-space): on request, a
  spot-averaged `|G|` picture for a `ProbeKernel` (`core.md`). It never feeds a classification and
  is kept by a new calculation or a reclassification.

Pruning (`applyRetention`) always runs *after* `identifyModes()`, never before — classification
needs to see the full eigenvector set to know which columns are significant, and the Ritz
candidates are computed from it. A standalone single-point object solves on its own axis; a sweep
needs one *common* axis instead (see below), which is exactly why the sweep classes call
`calcRitz`/`calcLim` and then `CalcCommon.greenStep` (`calcGreen(greenVelAxis)` → `calcGreenScore` →
`calcRitzScore` → `identifyModes` → `applyRetention`) as separate phase calls rather than
`calcAll()`. See "The Green velocity axis" below for the axis itself and for `nVelGreen`.

### Rotation convention: pass 0 to solvers after construction

The constructor already rotates `obj.structure` to `obj.angle`, and `RitzSolverSAW` and
`GreenSolver.calc` **also** rotate internally, by whatever angle they are given. Calling
`RitzSolverSAW(obj.structure, obj.angle)` after construction would therefore rotate the structure
twice: net `2×angle` instead of `angle`, for every angle ≠ 0. `calcRitz()` and `calcGreen()` in both
classes accordingly pass `0`, not `obj.angle`, to every solver call — the same rule `calcLim()`
follows (`GetStandardElastTensor(0)`, on the already-rotated structure).

**Rule for any rotation-touching class in this layer:** pass `0`, never `obj.angle`, to any solver
call made after the constructor has rotated `obj.structure`. `CalcIsofreqSingleLayered.calc()`
follows the same pattern (`GreenSolver.calc(obj.structure, 0, ...)`). The sweep and 2D classes
below only drive these single-point classes and never rotate a `Structure` or call a solver
themselves, so they inherit the rule from the building block one level down (see "2D sweep
orchestrators").

**The same rule applies to GUI and analysis code.** "Never calls a solver directly" holds for the
`Calc.*` classes, not for everything built on them. `ModalShapeViewerGUI`'s on-demand
Green's-response redraw calls `GreenSolver.calc` directly. Passing `sweepResult.structure` plus
the real angle is correct only for a genuine `CalcAngSweepHS`/`CalcAngSweepLayered`, whose
top-level `.structure` is UNROTATED; for a bare `CalcSingleHS`/`CalcSingleLayered` (packaged
internally as a one-angle sweep), `.structure` is already rotated and the same call
rotates it twice. **Rule for any code that reads `.structure` off a `Calc*` result and calls a
solver itself:** go through the per-angle element `cs = sweepResult.cs{i}` and pass
`cs.structure, 0`, which is correct in both cases. The same split exists on the isofreq side:
`CalcIsofreqAngSweepLayered.structure` is the UNROTATED base (each element rotates its own internal
copy), whereas `CalcIsofreqSingleLayered.structure` is ALREADY ROTATED to its own angle (its
constructor calls `rotateStructure(structure, angle)`), which is why `IsofreqModalShapeViewerGUI`'s
adapter for a single-point result tracks a `structureAlreadyRotated` flag and passes `angle=0`
(`gui_viewers.md`).

### `retentionLevel` — the eigenvectors, the largest per-angle array

`applyRetention(retentionLevel)` (public so the sweep classes can call it between phases)
controls how much of `ritzSolver.eigVecs_con` — the Ritz eigenvector matrix in `GammaMatrix`'s
constrained solution space — survives after classification:

| Level | Kept |
|---|---|
| `'full'` | everything (needed for modal-shape visualisation) |
| `'saw_only'` | only the columns `identifyModes()` classified as significant (`[classification.significantModes.modeNo]` — can be multiple columns: a pseudoRayleigh (HS) or SAW (layered) cluster, limiting-mode entries) |
| `'none'` | cleared entirely; velocities/scores stay, eigenvectors don't |

The scaling argument (the eigenvector matrix grows as `(3·Σ legO)²` per angle until `nMaxModes`
caps its columns, ~0.35–7 MB per angle, multiplied by every point a sweeping class drives) and the advice on which level to use are in `docs/pitfalls.md`
NUM4. `CalcSingleHS`/`CalcSingleLayered`/`CalcAngSweepHS`/`CalcAngSweepLayered` default to
`'full'` (affordable for a single point or a 1D sweep); `CalcWavevecAngSweep2Layered`/
`CalcWavevecAngSweep2Plate` default to `'none'`, because RAM multiplies across *both* axes there
(`CalcFreqAngSweep2Layered` runs no Ritz solve and has no retention level).
`'full'`/`'saw_only'` remain an explicit opt-in for a genuinely small 2D grid. See
`ritz.md`/`core.md` for the storage mechanism one layer down.

### `useParallel` from this layer's perspective

The angle-sweep classes (`CalcAngSweepHS`, `CalcAngSweepLayered`, `CalcIsofreqAngSweepLayered`,
and, forwarded unchanged, the 2D sweep classes) accept `useParallel` (`"auto"` (default) |
`"always"` | `"never"`), resolved once per `calc()` call by `ParallelPolicy.shouldUseParallel`
(`guliwerr/core/ParallelPolicy.m`) and reused for every phase loop in that call. `"auto"` switches
on once the swept dimension (`nAngles`, at whichever level is looping) reaches
`ParallelPolicy.AUTO_THRESHOLD` (20). It uses a thread pool (why: `core.md`, "ParallelPolicy"). By default the 2D sweep classes' outer loop (over wavelength/ω₀) stays serial
— MATLAB forbids nested `parfor` — and forwards `useParallel` to each inner
`CalcAngSweepLayered`/`CalcIsofreqAngSweepLayered`, which decides at the angle level. Benchmark
numbers and the thread-versus-process rationale are in `core.md`; `docs/pitfalls.md` NUM5 has
the short version.

Both 2D sweep classes also accept `flattenSweep` (default `false`): when `true`, `calc()`
dispatches across the *flattened* `(wavelength/ω₀, angle)` grid in one `parfor` per phase, so
`ParallelPolicy.shouldUseParallel` sees the *product* of both axes — useful when neither axis
alone reaches `AUTO_THRESHOLD` but their product does (e.g. 8 wavelengths × 15 angles). The
nested path (`calcFullPhases()`) and `CalcWavevecAngSweep2Layered.calcFlat()` share the same
`CalcCommon.forEachElement`/`greenStep` helpers, so the two cannot
diverge. `calcFlat()` still respects the per-wavelength shared-Green-axis barrier described below
(Phases 1 and 2 flatten across the whole grid; the axis is built per wavelength from that
wavelength's angles only, then Phase 3 flattens again — see the method's docblock).
`CalcFreqAngSweep2Layered.calcFlat()` has no barrier (isoFreq side, no phase split) and flattens
in a single pass. **Trade-off:** peak RAM is higher than the nested path, since all
`nWav × nAngles` (or `nOmega × nAngles`) elements are in flight through Phases 1–2 together — see
each class's `flattenSweep` docblock, especially when combining it with `retentionLevel='full'`
(required, for example, by `DispersionBranches.track`).

## Angle-sweep orchestrators: `CalcAngSweepHS` / `CalcAngSweepLayered`

These exist to use one identical Green velocity axis at every angle, producing genuine
`[nVel × nAngles]` matrices (`greenAbsG13EffMap`) rather than a ragged per-angle
structure. That is why they are not simply "call `CalcSingleHS.calcAll()` in a loop": each
angle's own axis ends at its own longitudinal limiting velocity, so the common axis needs Phase 2
(Christoffel) to have run across *all* angles first.
`calcFullPhases()` therefore runs explicit phases, each a `parfor`-or-`for` loop over the angle
array, with `obj.cs{iAngle}` extracted to a local variable before each loop and written back
afterwards (`parfor` needs a plain sliced variable, not a chained `obj.cs{i}` property access):

1. **Ritz**, all angles — `cs{iAngle}.calcRitz(precomp)`, with the sweep's shared precompute
   (see below).
2. **Christoffel**, all angles — `cs{iAngle}.calcLim()`.
3. **Common axis preamble**: the axis of the fastest angle, in one step for every angle — HS
   `GreenScore.gridAxis(1.05*max(longLimVel across angles), step)`; layered the lattice axis up to
   the largest bound `max(vMax, 1.05*bandVLmax)` (`CalcAngSweepLayered`; the sweep hands its
   current `vMax` to every angle first). Every angle's own axis is a leading part of it.
4. **Green**, all angles, via `CalcCommon.greenStep`: `calcGreen(greenVelAxis)` →
   `calcGreenScore(greenScoreRef)` → `calcRitzScore(precomp)` → `identifyModes()` ('green') →
   `applyRetention()`, one Green solve per angle, run *immediately* per angle inside the same loop
   iteration — not in a separate post-loop — to keep peak RAM at one angle's eigenvector matrix
   rather than all angles' at once.

Before Phase 3 the half-space sweep computes **one `GreenScore` `'low'` reference for the whole
structure** (`GreenScore.reference`, stored as the sweep's `greenScoreRef`, memoised per structure)
and hands it to every angle, so every angle — and a standalone angle of the same structure —
classifies against the same reference. `GreenScoreLayered` needs no reference. If
the run's method is not `'green'`, the sweep then calls `classify(method, sawThreshold)`.
An unset `sawClassMethod` (`''`, the default) is resolved first, on every run, by
`CalcCommon.resolveSawClassMethod` — `'green'` on `'full'`, `'ritz'` on `'ritzOnly'` — into a
local value (the property stays unset); it refuses another method than `'ritz'` on `'ritzOnly'`
and a numeric `sawThreshold` with the method unset, before any solve. Every angle is given the
run's resolved method and standard, and `greenStep` stamps its angle `calcMode = 'full'`
(`setRunMode`), so the angles match the sweep after a re-run in the other mode. `classify` never
writes `sawClassMethod` / `sawThreshold`.

**`'ritzOnly'` mode** has no shared velocity axis and no phases: each angle runs
`CalcCommon.ritzOnlyStep` (`calcAll('ritzOnly', retentionLevel)` with the sweep's
`sawThreshold`; the method is `'ritz'`). Nothing is shared between angles: each scores against its
own floor, so a standalone angle classifies exactly as inside the sweep.

`classify(method, threshold)` on a sweep sets the active classification of every angle from its
stored candidates and refreshes the `saw*`/`significantModesMap` maps.
`CalcWavevecAngSweep2Layered` resolves `sawClassMethod` the same way, forwards it and
`sawThreshold` to every wavelength's sweep and classifies per wavelength.

The constructor builds the angle-independent Ritz precompute (`nullConstrMat`, `dKdC_cell`,
`d2KdCdD_cell`) once, from the *unrotated* structure, and shares it across
every `cs{iAngle}`. This is correct, not merely convenient: these quantities depend only on
`calcParams.relDepth`/`mater.rho` (see `GammaMatrix.derivsByDepth`), never on `dirMat` or
rotation. Building it costs more than one Ritz solve (~0.15–0.2 s at `CalcParams(30,70)` — `relDepth` 30, the layer depth in wavelengths, and `legOrder` 70, the Legendre order of the depth basis — against
~30 ms per solve, machine-dependent), so sharing it is what makes a Ritz sweep cheap.

**The precompute is passed to every Ritz call explicitly — never left on the per-angle object.**
`precompStructForRitz` is `Transient` on every `Calc*` class, and a thread-pool `parfor`
serialises each object it hands to a worker (as `save` would), so a worker's copy of
`cs{iAngle}` arrives *without* it. The sweep therefore passes its own copy into
`calcRitz(precomp)` / `calcAll(..., precomp)` / `CalcCommon.ritzOnlyStep(..., precomp)`; inside
`CalcCommon.forEachElement`'s step handle it is broadcast once per worker, not copied per angle.
`Transient` is kept for speed as well as for small saved files: as an ordinary property, the
21 MB precompute (`CalcParams(30,70)`) would be serialised into every worker transfer, per angle
and per phase. Each sweep's `ritzPrecompute()` returns it, rebuilding it once after `load()`;
`CalcWavevecAngSweep2Layered.calcFlat()` hands each grid point its own wavelength's copy. A
single-angle object that has to build its precompute *on a worker* warns
(`CalcCommon.warnRitzPrecomputeBuiltOnWorker`, `<Class>:ritzPrecomputeBuiltOnWorker`) — a lost
speed-up, never a wrong result.

### The 0–180° angle grid, and Green's function on it

Both sweep classes default to `(0:2:180)'` (91 angles), not a full circle, for an exact reason:
Christoffel's `Γ_il(n) = C_ijkl n_j n_k` is quadratic in the propagation direction `n`, so
`v(θ) = v(θ+180°)` identically for every velocity Phases 1–2 compute (Ritz SAW, Christoffel
limiting, bulk) — for any elastic tensor, with no symmetry assumption (derived in `green.md`'s
"±k symmetrisation" section). Green's function does *not* share this symmetry, so a one-sided map
(`greenComponentMap('G13')`, `'G33'`, ..., `|G_+|`) must not be mirrored across 180°. The
±k-symmetrised `greenAbsG13EffMap` (the only Green map a sweep stores) needs no extra range or
pass, though: each angle's one solve gives it, on the same axis and at the same damping
(`green.md`, "±k symmetrisation").

For `CalcAngSweepLayered` the 0–180° default is *inherited* rather than re-derived for a layered
assembly: the quadratic-in-`n` argument holds per layer, and nothing in assembling several layers'
stiffness/mass matrices or the SAW boundary condition obviously breaks it, but it is a physics
*generalisation*, checked numerically at θ+180° rather than derived in full.

### The Green velocity axis

Every `'full'` angle runs **one** Green solve at `SawSpectrum.ETA` = 3e-3 on a uniform axis:

- **Half-space:** from 0.3 km/s to 1.05 × `longLimVel` + 0.1 km/s, in `greenVelStep` steps
  (default `GreenScore.DEFAULT_STEP` = 6 m/s), or the step that gives `nVelGreen` points; built as
  0.3 + n × step (`GreenScore.gridAxis`, `CalcSingleHS.greenAxis`). `GreenScore` reads only the
  band [0.8 vS, 1.03 vTf], fits each pole through samples ±16 m/s apart and gates analytically,
  so its labels do not depend on the step up to `GreenScore.MAX_STEP` (16 m/s; larger steps are
  capped, warning once). The step only sets the picture resolution.
- **Layered:** 4 m/s steps on the lattice 0.3 km/s + n × 4 m/s (`GreenScoreLayered.gridAxis`,
  `CalcSingleLayered.greenAxis`), from 0.3 × `bandVSmin` if that is lower than 0.3 km/s, up to
  `max(vMax, 1.05 × bandVLmax)` + 0.1 km/s. The layered sweep constructors have no `nVelGreen`:
  that positional argument is unused and must be empty.

Display and classification share the solve. The displayed `|G13_eff|` is therefore the magnitude of
the very spectrum the classifier scored, at the same damping.

**Sweeps: one common axis, scored per angle.** A sweep solves every angle on the axis of its
fastest angle, in one step for all angles, so the maps are rectangular. Every angle's own axis is a
leading part of the common one, and the existence scores read each angle's band, which is the same
set of grid points (bit for bit) in the common axis as in a standalone angle's own axis. GreenSolver
treats every velocity point independently, so a sweep angle's samples equal a solve on the shorter
axis to rounding. Neither the probe slope, nor the memoised `'low'` reference, nor the layered
direction-own normaliser involves another angle of the sweep: a standalone angle classifies
exactly as inside the sweep.

**Other components, on request.** `CalcAngSweepHS`/`CalcAngSweepLayered.greenComponentMap(component)`
returns a `[nVel × nAngles]` single map of any component (strings as for `greenComponent` above)
on the stored `greenVelAxis`: one Green solve per angle, same axis, same damping
(`SawSpectrum.ETA`). The map is returned, not stored. `CalcCommon.greenComponentMapFrom` assembles
it; `CalcCommon.greenG13EffMapFrom` assembles the stored `|G13_eff|` map (and returns `[]` if any
angle lacks it).

**Cost and RAM.** Timings are in `forwardCalcs_modeClassification.md`, section 10. In a sweep the
slower angles solve a little further, to the fastest angle's bound (a few per cent more points).
The maps are stored in single precision: a 181-angle half-space sweep holds ≈1 MB for
`greenAbsG13EffMap`, the only Green map it stores (other components are recomputed by
`greenComponentMap`, not kept). This is small next to the Ritz eigenvectors (`retentionLevel`,
above). Each `cs{i}` also keeps its own `|G13_eff|` and the complex band `greenScoreG13Eff`
(double, a few tens of KB per angle). The saving multiplies in `CalcWavevecAngSweep2Layered`,
which holds one sweep per wavelength.

**Plotting.** The axis starts at 0.3 km/s, below any surface wave, so `Vis.plotGreenMap` shows an
empty band at the bottom; the viewers default their velocity range to start at 0.9 × the slowest
Ritz mode.

### Convenience maps

`assembleMaps()` (private, called at the end of every `calc()` and `classify()`) is a plain
cellfun/cat extraction from `cs{}` with no classification logic of its own; it reads each angle's
active classification. `sawVelMap`/`sawModeNoMap`/`sawTypeMap` filter `significantModes` to the
SAW-type entries: HS uses `{'Rayleigh', 'pseudoRayleigh'}`; layered uses `{'SAW'}` (both layered
methods emit `'SAW'`). There is no `'Sezawa'` velocity-rank split — see the companion
ModeClassification doc, "No Love type, no Sezawa split". `bulkVelMap`/`greenAbsG13EffMap` exist
only in full mode; `'ritzOnly'` leaves the Green fields `[]`, which would make the extraction
cellfuns error on ragged input, so that block is skipped entirely.

## 2D sweep orchestrators: `CalcWavevecAngSweep2Layered` / `CalcFreqAngSweep2Layered`

Both add exactly one outer loop (over wavelength or ω₀) around the 1D sweep classes above:
neither reimplements the Ritz→Christoffel→Green phase loop, and neither constructs
`RitzSolverSAW`/`GreenSolver`/`CalcSingleLayered` directly. Because they never rotate a
`Structure` themselves, rotation correctness is inherited from the building block one level down
(see "Rotation convention: pass 0 to solvers after construction").

### Precompute sharing differs between the two — and that difference is physical, not incidental

`CalcWavevecAngSweep2Layered` builds `rawByDepthCell` (`GammaMatrix.rawByDepth`'s
Legendre-integral tensor) once, in its constructor, and passes it to every `CalcAngSweepLayered`
it constructs. This works because `rawByDepth` is independent of material and `relDepth` — it
depends only on `calcParams.legOrder`/`basisVecL`/`legDiag`/`legNondiag`, all fixed across the
sweep; `relDepth` enters only `derivsByDepth`'s cheap rescaling step. The expensive tensor is
therefore built `nLay` times for the *whole* 2D sweep, not `nLay × nWav` times. This relies on
`legOrder` staying constant across the sweep, which nothing enforces: varying `legOrder` by
wavelength without rebuilding `rawByDepthCell` would silently invalidate the cache.

`CalcFreqAngSweep2Layered` has **no** equivalent optimisation because none is available:
`GreenSolver.calc` has no analogue of Ritz's expensive/cheap `rawByDepth`/`derivsByDepth` split,
so every `(angle, ω₀)` point is an independent, full-cost call. It is a plain nested loop with no
shared precompute. The asymmetry between the two 2D classes follows from which solver each drives.

### `retentionLevel='none'` default, and no top-level aggregation

`CalcWavevecAngSweep2Layered.retentionLevel` defaults to `'none'` (unlike
`CalcAngSweepLayered`'s `'full'`) — see the RAM discussion above. For the same reason,
the Green maps are **not** aggregated into a top-level 3D/4D array: each `sweeps{iWav}` (a full
`CalcAngSweepLayered`) already holds its own, so read `sweeps{iWav}.greenAbsG13EffMap` (or call
`sweeps{iWav}.greenComponentMap(component)` for another component) for a given wavelength slice. `CalcFreqAngSweep2Layered` follows the
same pattern for `sweeps{iOmega}.absGMap`, for consistency, even though Green spectra alone are
far smaller than eigenvectors.

### Real-unit awareness lives here, not below

`CalcWavevecAngSweep2Layered` is the one real-unit-aware class in this family: every film
layer's `CalcParams.realThickness` must be set (checked in the constructor before any work
starts), and `structure.ApplyWavelength(wavelengths(iWav))` converts it to `relDepth` at each
wavelength. The substrate (layer 1) is **not** swept — its `relDepth`/`legOrder` stay as set in
the reference structure (see `ApplyWavelength`). It is the caller's responsibility to make the
substrate `relDepth` conservatively deep relative to the *smallest* wavelength in the sweep
(small wavelength → large film `relDepth` → the substrate must be deeper still); only a
convergence study shows how deep is enough.

## isoWavevec vs. isoFreq — the orchestration-level view

`green.md` covers the solver-level physics of the two scan modes (fixed `|k∥|=1` sweeping ω, vs.
fixed ω₀ sweeping `|k∥|`); this section covers how the Calc classes expose them.

Everything above (`CalcSingle*`, `CalcAngSweep*HS/Layered`, `CalcWavevecAngSweep2Layered`) is
**isoWavevec-side**: fixed `relDepth` (equivalently, a fixed real wavelength once
`ApplyWavelength` is applied), Green swept over velocity. The isoFreq side is a simpler, parallel
family:

| | isoWavevec side | isoFreq side |
|---|---|---|
| Atomic single-point class | `CalcSingleHS`/`CalcSingleLayered` | `CalcIsofreqSingleLayered` |
| Angle-sweep class | `CalcAngSweepHS`/`CalcAngSweepLayered` | `CalcIsofreqAngSweepLayered` |
| 2D sweep class | `CalcWavevecAngSweep2Layered` | `CalcFreqAngSweep2Layered` |
| Class kind | handle | handle for the sweeps; **value** for the single-point `CalcIsofreqSingleLayered` (the `cs` elements) |
| Ritz involved? | yes | **no** |
| Phased/shared-axis structure? | yes (3-phase, shared Green axis) | no — flat, independent loop |
| Retention policy? | yes | not applicable (no eigenvectors) |

The isoFreq side has no Ritz, for a physical reason: Ritz fixes `|k∥|=1` by construction, so it
has no free parameter that can hold real frequency fixed while `k∥` varies — the isoFreq question
is meaningless to Ritz, not merely unimplemented. `CalcIsofreqSingleLayered` therefore has no
`RitzSolverSAW`, no eigenvector retention, and no `ModeClassification` matching (there is no Ritz
result to confirm a Green peak against): it calls `GreenSolver.calc` once with
`scanMode='isoFreq'`, converting the physically meaningful `kParRange` into the `velRange`
argument `GreenSolver.calc` expects (`v = ω₀/k∥`, elementwise, order-preserving). With no large
per-object state, `CalcIsofreqSingleLayered` is a **value class**, matching
Christoffel/GreenSolver/RitzSolverSAW's default-construct-then-`calc()` pattern. The two sweeps are
handle classes (with `copyElement()`) for consistency with the rest of the family;
`CalcFreqAngSweep2Layered` also holds other handles.

The isoFreq sweep loops are the simplest `parfor` case in the project: value-class elements (no
handle aliasing), no shared precompute, and every angle's/ω₀'s calculation fully independent.

`Calc.freqToOmega0(freqHz, lambdaRefNm)` converts a real experimental frequency into
`GreenSolver.calc`'s normalised `omega0` (`v = f·λ` at the reference wavelength also passed to
`ApplyWavelength`) — full parameter docs in `docs/reference/FullFunctionReference.html`.

For a half-space the two scan modes coincide: with no layer depth scale to introduce dispersion,
`v=ω/|k∥|` is the only combination entering either mode's Christoffel solve, so both give the
same `G_IJ(v)` curve (derived in `green.md`).

## `@Calc` — the namespace launcher

`guliwerr/forwardCalcs/@Calc/` is a thin static-method namespace (`Calc.singleHS(...)`,
`Calc.angSweepLayered(...)`, etc.) over the `Calc*` classes (the plate ones included) — the same
pattern as the `@Christoffel`/`@GammaMatrix`/`@GreenSolver` folders in the solver layer. Each
launcher is a one-line pass-through (`obj = CalcAngSweepHS(varargin{:})`). Other launchers:
`Calc.modeClassification()`, `Calc.freqToOmega0(...)`, and `Calc.saveResult` / `Calc.loadResult`
(over `CalcCommon.saveCalc` / `loadCalc`; `docs/reference/functionReference.md`). Do not instantiate `Calc`
itself — it exists so calling code writes `Calc.angSweepHS(...)` rather than qualifying every class
name.

## References

`docs/pitfalls.md` NUM4 and NUM5 (retention/parallel summaries); `forwardCalcs_modeClassification.md`
(the mode classification this layer drives: methods, standards, `classify`, stored fields); `green.md`
(GreenSolver internals — isoWavevec/isoFreq scan-mode physics, ±k/G_eff symmetrisation,
near-isotropic degeneracy fix); `ritz.md` ("Per-mode diagnostics, the Ritz score and the Ritz
spectrum").
