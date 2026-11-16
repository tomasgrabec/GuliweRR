# core — ElastConsts, Material, CalcParams, Structure, ParallelPolicy, ProbeKernel

`guliwerr/core/`. The building blocks every solver (Ritz, Christoffel, GreenSolver) is handed as input: a material's elastic tensor, its orientation, the numerical mesh for one layer, the layer stack itself, a shared parallel-execution policy, and (experimental) the wavevector weight of a finite measurement. None of these classes do any wave physics themselves — they describe the problem; the solvers in `guliwerr/solvers/` consume it.

## Composition: how the four data classes fit together

```
ElastConsts   — one material's stiffness tensor (three equivalent representations)
     │
     ▼
Material      — ElastConsts + density + orientation (dirMat) for one physical layer
     │
     ▼
CalcParams    — Legendre-basis mesh resolution (legOrder, relDepth) for one layer
     │
     ▼
Structure     — a stack of (Material, CalcParams) pairs, bottom layer first
```

`Material` and `CalcParams` describe orthogonal things about a layer — orientation/physics vs. numerical resolution — so they stay separate objects rather than merging into one. `Structure` pairs them up per layer (`obj.mater(i)`, `obj.calcParams(i)`) and is what actually gets passed to a `Calc.*` sweep class. The four data classes are value classes (see "Value semantics" below); `ParallelPolicy` is stateless and sits outside the pipeline, used by the sweep classes rather than by `Structure` itself; `ProbeKernel` is an optional input to the spot-averaged Green map (below).

## ElastConsts — stiffness tensor storage and rotation

Stores one material's elastic stiffness in three equivalent, always-synchronised representations: `vector` (symmetry-class-specific independent constants), `voigt` (6×6 contracted notation), `tensor` (full 3×3×3×3). All three are rebuilt together on construction or on `SetNewSymmClass` — there is no path that lets them drift out of sync, which matters because downstream code reads whichever representation is convenient (`voigt` for moduli, `tensor` for rotation) without re-deriving it.

The constructor accepts any of the three forms: `FindType` dispatches on shape (`ndims==4` → tensor, a 2-D array with any dimension equal to 1 → vector, `6×6` → Voigt, anything else → `'unknown'`, which the constructor rejects), and `Any2Tensor` builds the canonical tensor from whichever was given.

The symmetry class (`symmClassNo`; `SymmClassLibrary` lists class names and constant counts, from isotropic = 1 with 2 constants to triclinic = 10 with 21) fixes how many independent constants are expected. The constructor errors on a mismatched count rather than letting `Vector2Voigt` silently drop the extra entries.

**`SetNewSymmClass` is safe in one direction only.** It rebuilds from `obj.voigt` under the new class with no compatibility check: `Any2Tensor` maps a Voigt matrix straight into the tensor regardless of `symmNo` (only vector input uses `symmNo`, via `Vector2Voigt`), so the new `vector` is whichever Voigt positions the new class reads (`Tensor2Vector`/`IndepVoigtPositions`). Moving to a **higher class number** (more constants, e.g. hexagonal→triclinic) is always safe — previously tied entries simply become free. Moving to a **lower class number** (e.g. hexagonal→cubic) silently discards the rest, with no check that they were consistent with the tighter symmetry.

`GetYoungsModulus`/`Tensor2Youngs` never forms a compliance *tensor* (which would need the Voigt shear-strain factors 1/2 and 1/4 that stiffness Voigt notation does not): it rotates the *stiffness* tensor to the query direction, converts to Voigt, inverts the 6×6 and reads `compliance(1,1) = 1/E`.

`VoigtAverages`/`ReussAverages`/`HillAverages` and `Voigt2Bulk` give isotropic-equivalent moduli (`[M,G]` = P-wave/S-wave modulus pair, `M=K+4/3·G`): Voigt assumes uniform strain (upper bound), Reuss uniform stress (lower bound, via the compliance matrix — its `K` uses the same `1/sum(compliance 3×3 block)` formula as `Voigt2Bulk`), and Hill is their arithmetic mean (the standard practical polycrystal estimate). Voigt and Reuss coincide exactly for isotropic input. `CubicConversion`/`TetragonalConversion`/`HexagonalConversion` (each with a `*Deconversion` inverse) re-express a class's Voigt constants in a physically transparent basis (e.g. hexagonal `[C11,C33,G1,G2,G3]` = the two uniaxial moduli plus axial/basal/breathing shear) via a fixed linear transform, `vecOut = transfMat*vecIn`; `*Deconversion` applies its matrix inverse.

`RotateTensor` applies a full-tensor rotation, `C'_{ijkl} = R_{ia}R_{jb}R_{kc}R_{ld}C_{abcd}`, via four sequential mode-*n* matrix products (permute the target dimension to the front, left-multiply by the 3×3 rotation matrix, permute back) rather than the mathematically equivalent 8-nested-loop / 3⁸-term sum. `ElastRotReductMat` calls `RotateTensor` once per independent constant (up to 21, for triclinic) and is itself called once per angle inside an angular sweep, so a 360-point sweep of a low-symmetry material is thousands of calls — the vectorised form is what keeps that cheap.

`GetReducedElastConsts`/`ElastRotReductMat` produce the 15-of-21 reduced constant set the Ritz solver's `GammaMatrix` pipeline builds from (`GetIndOf2DReduct_15of21` names the index set: every Voigt entry touching index 2, the in-plane completion axis) — why the reduction is exact, and its coordinate-convention dependency, in `docs/reference/architecture/ritz.md`, "The 21-to-15 elastic-constant reduction".

`CalcAnisotropyRatio` (generalised Zener ratio, `(v_shear_max/v_shear_min)²` scanned over a coarse angular grid) lives here rather than in `christoffel/` because it is conceptually a property of the stiffness tensor, even though it calls `Christoffel.phaseVels`. `init.m` puts every folder on the MATLAB path, so this `core/ → solvers/christoffel/` dependency is harmless. `ElastConsts.GetAnisotropyRatio` (instance method) uses the exact closed-form Zener formula for cubic symmetry (`2*c44/(c11-c12)`) and falls back to the angular scan only for lower symmetries; it feeds the anisotropy-factor tiers for choosing `legO`/`relDepth` (`docs/pitfalls.md` MTH1). The Zener ratio can be below 1; the property `shearAnisotropyRatio` holds the rotation-invariant `(v_shear_max/v_shear_min)²` ≥ 1 instead, set by the constructor where a closed form exists — isotropic 1, cubic `max(A_Z, 1/A_Z)` (the cubic shear extremes are `c44` and `(c11−c12)/2`) — and `NaN` otherwise, because a scan would cost more than the constructor, which finite-difference loops call many times. `GreenSolver.calc`'s near-isotropic screen reads it and scans (`CalcAnisotropyRatio`) only when it is `NaN` (see `green.md`).

## Material — orientation and the coordinate convention

`Material` bundles an `ElastConsts` object, density `rho`, and orientation (`dirMat`) for one physical layer. It is the single place where "which direction is propagation" gets defined.

GuliweRR uses two coordinate conventions, each suited to what it's for:

**Input coordinates — `dirMat`, for specifying an orientation.** A 3×3 matrix of row vectors, each a direction cosine, completed to a right-hand orthonormal frame by `CompleteDirMat` if fewer than 3 rows are given. Rows are ordered by physical importance — what defines the cut first, what's needed to place it in the lab second:

| Row | Meaning |
|---|---|
| `dirMat(1,:)` | surface normal (plane) |
| `dirMat(2,:)` | propagation direction reference |
| `dirMat(3,:)` | in-plane completion — `cross(row1,row2)` |

`Material.CalcRotMat(dirMat, angleDeg)` rotates rows 2 and 3 in-plane by `angleDeg` (pivoting on row 1) and returns a 3×2 matrix `[rotated-propagation | surface-normal]` — this is the format `ElastConsts.RotateTensor` consumes directly, so a `Material` object's own `GetRotElastTensor`/`GetRotElastVoigt` can be called with just an angle.

`Miller2DirMat`/`DirMat2Miller` (Miller/Miller-Bravais index ↔ direction-cosine row vector) and `Euler2DirMat`/`DirMat2Euler` (Bunge-ZXZ Euler angles ↔ a genuine `dirMat`, row 1 = surface normal) live here because they all produce or consume a `dirMat`. `Miller2DirMat` also needs crystal-system geometry (`symmNo`/`latticeParams`, via the private `LatticeParamSpec`/`LatticeParamSpecReport`, which call `ElastConsts.SymmClassLibrary` for the class name).

**`Euler2DirMat`/`DirMat2Euler` are thin `dirMat`-facing wrappers, not the raw conversion.** The Bunge-ZXZ angle↔matrix maths lives in `Euler2DCM`/`DCM2Euler`: a generic conversion between three angles and *any* proper-rotation matrix, in Bunge's own convention (row/column 3 is the distinguished tilt axis), with no `dirMat` meaning. Bunge's frame is not a third GuliweRR frame — it coincides row for row with the computational coordinates' final frame (`[propagation; completion; plane]`, below), whereas `dirMat`'s row 1 (surface normal) sits in Bunge's row/column 3. `Euler2DirMat`/`DirMat2Euler` apply the `[3,1,2]`/`[2,3,1]` reorder (an even, determinant-preserving permutation) around `Euler2DCM`/`DCM2Euler`, so they always take and return a real `dirMat`. For the raw conversion (e.g. checking a `dcm` given in a paper's own Bunge convention, or an arbitrary rotation with no `dirMat` meaning), call `Euler2DCM`/`DCM2Euler` directly.

**Computational coordinates — the internal solver frame** (`x1=propagation, x2=in-plane, x3=surface-normal/depth`, i.e. z=out-of-plane), for computing — the usual convention in the literature this project cites (Every, Maznev). Every solver (Ritz's `GammaMatrix`, Christoffel, `GreenSolver`) computes natively in this frame. `Material.GetStandardDirMat(dirMat, angleDeg)` is the single conversion point: it rotates within the native `dirMat` row layout first (rows 2 & 3, via `RotDirMatFromAngle`), *then* reorders rows `[3,2,1]` to produce `[completion; propagation; plane]`. The order matters — `CalcRotMat`/`RotDirMatFromAngle` always mix rows 2 & 3, so rotating an already-reordered `dirMat` would rotate propagation towards the surface normal instead of within the surface plane. Every solver routes through `GetStandardDirMat` once per layer rather than repeating the permutation inline, so the convention lives in one place; the standard-frame elastic tensor itself comes from `Material.GetStandardElastTensor(angleDeg)` (`GetStandardDirMat`, then `GetRotElastTensor(0)` — the angle applied once), which `GreenSolver`, the limiting velocities (`CalcCommon.limitingAndBulk`), `RitzScore`, `SpotGreen`, `SensitivityAnalysis`, and the Christoffel plots all call.

The `[3,2,1]` reorder is chosen so that feeding the result back through `CalcRotMat`/`RotateTensor` (which read row 1 as pivot/normal and rows 2–3 as the in-plane pair) lands propagation on `R1` (=`x1`), completion on `R2` (=`x2`), and `cross(R1,R2)` on `R3` (=`x3`). Given `CompleteDirMat`'s right-hand convention `completion = cross(normal, propagation)`, `cross(R1,R2)` is the surface normal (via `cross(a, b×a) = b` for orthonormal `a,b`), so the frame is exactly `x1=propagation, x2=in-plane, x3=surface-normal`, with no x2/x3 swap.

**`GetStandardDirMat`'s raw 3×3 output is not itself a proper rotation.** `stdDirMat = rotDirMat([3,2,1],:)` is a single row transposition — an odd permutation — so `det(dirMat)=+1` (right-handed, as `CompleteDirMat` constructs) but `det(GetStandardDirMat(dirMat))=-1`. This is harmless in the chain above: `CalcRotMat` uses row 1 as the pivot and combines rows 2–3 into a 3×2 result, and `RotateTensor` then recomputes the third axis via `cross(R1,R2)`, which is right-handed by construction. It matters only to a consumer that treats the raw output as a 3×3 rotation matrix, as Euler-angle extraction does: `Material.DCM2Euler` (which `DirMat2Euler` wraps) does **not** accept `GetStandardDirMat`'s output — it returns a number, but not a physical rotation. The correct handedness-preserving reorder is the `[2,3,1]` permutation applied inside `DirMat2Euler`; see its docstring.

**Do not confuse the two conventions.** `dirMat`'s row order (`[plane; direction; completion]`) is the user-facing input format and is unrelated, row for row, to the `x1/x2/x3` solver frame produced by `GetStandardDirMat`. Code indexing a raw `dirMat` and code indexing a `GetStandardDirMat` output are reading different orderings of the same three physical directions.

`GetRotMat(angleDeg, relativeToggle)` offers a second axis of choice unrelated to the above: `relativeToggle=0` (default) rotates within the layer's own `dirMat`; `relativeToggle=1` rotates in the absolute `[100]/[010]` crystallographic frame regardless of `dirMat` — useful for sweeping a crystal-relative angle independent of how the layer happens to be mounted.

## CalcParams — Legendre-basis mesh resolution

One `CalcParams` per layer. The Ritz method expands each layer's depth-dependent displacement field in a Legendre-polynomial basis truncated at order `legOrder`, over a domain of relative depth `relDepth` (layer thickness measured in wavelengths). `LegendreMatrixDiagNondiag` precomputes the two integral matrices the Ritz eigenproblem needs — `D_mn = ∫ P'_m P'_n dx` (symmetric) and `N_mn = ∫ P'_m P_n dx` (lower-triangular, nonzero only where `m+n` is odd) — once per `CalcParams` construction, not per solve, since they depend only on `legOrder`.

The two must satisfy `legOrder ≥ 2 × relDepth` — why, in `docs/reference/architecture/ritz.md`, "Legendre order (`legO`) and `relDepth`"; operational tiers in `docs/pitfalls.md` MTH1.

`realThickness` (default `NaN`) is deliberately separate from `relDepth`: `relDepth` is what every solver consumes and must always mean the same dimensionless quantity, whereas `realThickness` is how a real sample is described (e.g. in microns). The same layer is probed at many wavelengths across a sweep, so the conversion happens per wavelength at calculation time, in `Structure.ApplyWavelength` (below).

`SetLegOrder` recomputes every dependent field (`legDiag`, `legNondiag`, `basisVecL`, `polyLeg`), so a `CalcParams` object can never hold a stale mesh for its own `legOrder`.

## Structure — the layer stack

A `Structure` is an ordered stack of `(Material, CalcParams)` pairs — `mater(1)`/`calcParams(1)` is the bottom (substrate) layer, `mater(nLay)`/`calcParams(nLay)` is the top (free-surface) layer. `AddLayer`/`RemoveLayer`/`ModifyLayer` all return the modified object rather than mutating in place (value-class semantics, see below); `strBasisVecL` (the sum of `basisVecL` across all layers, needed to size the assembled Ritz eigenproblem) is kept incrementally consistent by every one of those three methods rather than recomputed from scratch. `isMirrorSymmetric(tol)` (mid-plane mirror: each layer's tensor against the mirror image of its partner's) and `isSagittalMirrorSymmetric(angleDeg, tol)` (the `x2 → −x2` mirror at one angle) are the exact symmetry tests the free-plate mode typing relies on (`plate.md`).

`ApplyWavelength(obj, wavelengthReal)` derives `relDepth = realThickness / wavelengthReal` for every **film** layer (`2:nLay`), turning a fixed real thickness into a wavelength-relative depth that changes at every sweep point. The **substrate (layer 1) is deliberately skipped.** `legOrder` is held fixed across a wavelength sweep (required for `GammaMatrix.rawByDepth` caching), so a substrate `relDepth` derived from a real thickness could grow arbitrarily at long wavelengths and silently break the `legO ≥ 2×relDepth` rule. The substrate's `relDepth`/`legOrder` instead stay at the value chosen for the whole sweep — a truncation-depth parameter, not a real thickness, since a substrate is nominally semi-infinite — and can still be changed via `SetRelDepth`.

## ParallelPolicy — thread pool, not process pool

A stateless static class shared by every angle, wavevector, and frequency sweep class (the `Calc*` classes behind the `Calc.` launchers, including the plate and 2D sweeps; each carries a `useParallel` property), and by `SpotGreen`, `GreenScore`, `GreenSpotScore`, and `SensitivityAnalysis`, with two jobs: (1) `shouldUseParallel(useParallelSetting, nItems)` resolves a sweep's `useParallel` property (`"auto"`/`"always"`/`"never"`) against `AUTO_THRESHOLD=20`; (2) `ensureThreadPool()` starts (or swaps in) a MATLAB **thread** pool, falling back to serial execution with a warning if Parallel Computing Toolbox is unavailable. `isOnWorker()` reports whether the caller already runs on a pool worker; a loop that parallelises on its own (e.g. `GreenScore.reference`) stays serial there, since a worker has no pool to start.

A thread pool shares the client's memory, so broadcast data — elastic tensors, `CalcParams`, and especially the cached `rawByDepthCell` — is not duplicated per worker. The per-angle results are what a sweep's RAM grows with (at `'full'` the eigenvectors are the largest part, ~0.35–7 MB per angle, `docs/pitfalls.md` NUM4), and a process pool would add the broadcast cost once per worker on top. On a production-sized sweep (hundreds of angles, `legOrder` ~70, `relDepth` ~30) a thread pool is several times faster than serial at a modest memory cost, while a process pool of the same width is slower and uses several times the RAM.

`AUTO_THRESHOLD=20` is a single fixed crossover, not scaled by `legO`/`nVelGreen`. Thread-based `parfor` parallelises raw numeric kernels well but classdef method-dispatch overhead (`validateattributes`, string-mode switches) poorly, so the speed-up appears only once the per-angle numeric cost (`GreenSolver.calc` at realistic `nVelGreen`) dominates; small sweeps (a few dozen angles, modest `legO` and `nVelGreen`) show essentially no speed-up. A cost-aware threshold would need a real cost model. Use `useParallel="never"` for debugging, since `parfor` mangles stack traces and breakpoints.

## ProbeKernel — wavevector weight of a finite measurement (experimental)

A measurement couples to a cloud of wavevectors around the nominal ±k₀, weighted by `W(k)` (the product of the Fourier transforms of the source pattern and the detection weight). `ProbeKernel` only describes `W`; it knows nothing about Green's functions. Constructors: `point(component)` (no smearing, the reference), `gaussian`, `tgs`, `sras`, `fromRealSpace`, and the generic `ProbeKernel(weightFcn, ...)`; `sample` returns the quadrature nodes and normalised complex weights. The half-space `calcSpotGreenMap(kernel)` consumes it (`SpotGreen`); the result is a picture only and no classification reads it. Half-space only, and not yet validated beyond one material — see `help ProbeKernel` for the constructors' arguments.

## Value semantics

The four data classes are value classes (none inherit `handle`). Every mutating method (`SetDirMat`, `AddLayer`, `ApplyWavelength`, …) returns the modified object, which the caller must capture (`obj = obj.AddLayer(mat, cp)`; `obj.AddLayer(mat, cp)` alone modifies a temporary copy). In return, a `Structure` handed to a sweep cannot be mutated underneath it by other code holding the same reference — the natural fit for "assemble a `Structure` once, then pass it read-only into a solver".

## References

The `x1=propagation, x2=in-plane, x3=surface-normal` convention matches the layered-medium literature the solvers cite (Every 1997/2002/2015, Maznev 1999 — see `docs/reference/architecture/green.md` for the full reference list). `core/` has no independent literature dependency beyond that shared convention. The Legendre-basis Ritz formulation that `CalcParams` serves is from Stoklasová et al. (2015) and Grabec et al. (2024); full citations in `docs/background.md`, "Key references".
