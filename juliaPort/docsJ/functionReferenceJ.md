# GuliwerrJ — function reference

One-line-per-function index of every public entry point, grouped by category. This is a cheat
sheet, not full documentation — every function below has its own Julia docstring (`?name` at
the REPL) with the complete argument list, defaults, and — where relevant — the reasoning behind
any deliberate divergence from MATLAB. Start here to find *which* function to call.

Built from `juliaPort/guliwerrJ/src/GuliwerrJ.jl`'s own export list and checked against the
constructor and function signatures in the source. All names below are available after
`using GuliwerrJ` (core numerics) or `using GuliwerrJ, CairoMakie`/`GLMakie` (plotting and the
GUI launchers) — see `quickStartJ.md` for environment setup. The forward calculations come
first; the experimental parts (classification, spot-averaged Green maps, the inverse module)
are marked where they appear.

**Every constructor except `ElastConsts` and the inverse launchers is keyword-only.**
`Calc.angSweepHS(structure, angles, 'ritzOnly', ...)` became
`angSweepHS(; structure = ..., angles = ..., calcMode = :ritzOnly, ...)` — no positional
shortcut exists, on purpose (`pitfallsJ.md` jOTH1). Where MATLAB took a quoted string
(`'ritzOnly'`, `'full'`), Julia takes a `Symbol` (`:ritzOnly`, `:full`) — see `pitfallsJ.md`
jOTH2.

**MATLAB-only.** `Calc.saveResult`, `Calc.loadResult`, `Gui.viewSaved` and `Gui.viewResult`
have no Julia counterpart: a GuliwerrJ result lives in memory for the session. Everything else
that differs from MATLAB is listed in `matlabVsJuliaJ.md`.

## Core building blocks (`src/core/`, module `CoreTypes`)

Build order for any calculation: `ElastConsts` → `Material` → `CalcParams` → `Structure` — same
order as MATLAB.

| Name | Purpose | Call | Key inputs |
|---|---|---|---|
| `ElastConsts` | Elastic constants (vector / 6×6 Voigt / 3×3×3×3 tensor forms) for one symmetry class. | `ElastConsts(consts, symmNo=nothing)` | `consts`: vector/Voigt/tensor (positional, not keyword — matches MATLAB's own two-positional-argument shape); `symmNo` (opt, auto-detected if omitted) — see `SymmClassLibrary` for the constant count/order per class (isotropic = `[c11,c44]`, cubic = `[c11,c12,c44]`) |
| `Material` | Elastic constants + density + crystal orientation for one layer. | `Material(; elConsts=nothing, rho=nothing, dirMat=nothing, name=nothing)` | `elConsts`: an `ElastConsts`, or a raw vector/matrix/`Tuple` passed straight to `SetElastConsts`; `rho`: density [g/cm³]; `dirMat` (opt, default a fixed reference orientation, e.g. Cu(001)); `name` (opt) |
| `CalcParams` | Discretisation settings (Legendre order, depth) for one layer. | `CalcParams(; relDepth=20.0, legOrder=50, realThickness=NaN)` | `legOrder` **should satisfy `legOrder >= 2*relDepth`** (`legOrder >= 2.3*relDepth` for the Ritz-only half-space classification), see `docs/pitfalls.md` MTH1 and `matlabVsJuliaJ.md` §7; `realThickness` (opt, real units, for `ApplyWavelength`) |
| `Structure` | A stack of layers (layer 1 = bottom/substrate). | `Structure(; materIn=nothing, calcParamsIn=nothing, name=nothing, nMaxModes=nothing)` | `materIn`/`calcParamsIn` (opt — omitting `materIn` adds no layer, same as MATLAB's `isempty` guard); `nMaxModes` default 500 |

`Structure` layer editing (each returns a **new** `Structure` — immutable, unlike MATLAB's
`handle`): `AddLayer(obj, mater0=nothing, calcParams0=nothing)` (appends on top),
`RemoveLayer(obj, ind)`, `ModifyLayer(obj, indLay, mater0, calcParams0)`,
`ApplyWavelength(obj, wavelengthReal)` (derives film `relDepth` from `realThickness`, substrate
excluded), `isMirrorSymmetric(obj, tol=1e-10)`, `isSagittalMirrorSymmetric(obj, angleDeg, tol=1e-10)`. `nlay(obj)`/`str_basis_vec_l(obj)` are
functions, not properties.

`CalcParams` setters (each returns a new object, rebuilding every derived field): `set_rel_depth`,
`set_real_thickness`, `set_leg_order`.

`ElastConsts` free functions (module-level, no `.` namespace prefix — MATLAB's
`ElastConsts.<Name>(...)` statics become plain functions): `Voigt2Vector`, `Vector2Voigt`,
`Any2Tensor`, `Tensor2Voigt`, `Tensor2Vector`, `SymmClassLibrary`, `FindSymmClass_input`,
`IndepVoigtPositions`, `GetIndOf2DReduct_15of21`, `CubicConversion`/`CubicDeconversion`,
`HexagonalConversion`/`HexagonalDeconversion`, `TetragonalConversion`/`TetragonalDeconversion`,
`RotateTensor`, `ElastRotReductMat`, `RightHandFrame`, `Tensor2Youngs`, `Voigt2Bulk`,
`HillAverages`/`VoigtAverages`/`ReussAverages`, `ZenerRatio`. MATLAB's
`ElastConsts.CalcAnisotropyRatio` is the free function `calc_anisotropy_ratio` (see "Solvers"
below).

`Material` orientation functions: `Euler2DCM`/`DCM2Euler` (raw Bunge-frame rotation, no
`dirMat` meaning), `Euler2DirMat`/`DirMat2Euler` (`dirMat` ↔ Bunge-ZXZ Euler angles, the
reorder applied internally — matches the fixed MATLAB behaviour), `Miller2DirMat`/
`DirMat2Miller`, `CompleteDirMat`, `GetStandardDirMat`, `RotDirMatFromAngle`, `CalcRotMat`.
`SetElastConsts`/`SetDirMat`/`SetDensity`/`SetName` and `GetRotMat`/`GetRotElastTensor`/
`GetRotElastVoigt`/`GetRotDirMat` are the setter/getter functions taking a `Material` as their
first argument, same role as the MATLAB instance methods of the same name.

Materials library: `MaterialFromLibrary(name; dirMat, T, variant)` returns `(Material, entry)`
and `ListMaterialLibrary(; print=false)` the rows -- the ports of `Material.FromLibrary`/
`Material.ListLibrary`, reading the same `guliwerr/data/materialsLibrary.csv` (one file for both
languages; `entry.structure` gives the crystal structure, e.g. `"fcc"`, `"hcp"`). Renamed because Julia's `Material` is a struct and cannot carry static methods.

Utilities: `guliwerrVersion()` (the GuliweRR version string, the same one MATLAB reports; the package's
own `Project.toml` version tracks the packaging only), `blasBackend()`
(`:accelerate` or `:openblas`, see `matlabVsJuliaJ.md` §3.4).

## Half-space forward calculations (`src/forwardCalcs/`, `Calc.jl` launchers)

`singleHS`/`angSweepHS` etc. are thin keyword-forwarding aliases onto the underlying `CalcSingleHS`/
`CalcAngSweepHS` struct constructors (MATLAB's `Calc.<name>` namespace collapses to bare
functions — no `Calc.` prefix needed, these are plain exported names). Use `angSweepHS` even
for a single angle (pass a 1-element `angles` vector) — it builds the shared Ritz precompute
internally; `singleHS` is the lower-level per-angle building block.

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface
> waves in strongly anisotropic materials. Thresholds, scores and the interface may change.

The classification (`classify`, `reclassify`, the `sawClassMethod`/`sawThreshold`/`greenVelStep`
fields, `isStale`, `refinedSawAgreement`) is one-to-one with MATLAB, with `Symbol`s for
strings. The default standard is per method: `:high` for `:green`/`:refined` on a half-space,
`:low` for `:ritz` and for layered stacks.

| Name | Purpose | Call | Key inputs |
|---|---|---|---|
| `singleHS` → `CalcSingleHS` | Single-angle half-space container (Christoffel + Green + Ritz + classification). | `singleHS(; structure=nothing, angle=0, precompStruct=NamedTuple(), calcMode=:full, retentionLevel=:full, nVelGreen=nothing)` | `structure`: **unrotated** `Structure`; `calcMode`: `:full`\|`:ritzOnly`; `retentionLevel`: `:full`\|`:saw_only`\|`:none` |
| ↳ `calcRitz(cs)` / `calcLim(cs)` / `calcGreen(cs)` / `calcGreenScore(cs)` / `calcRitzScore(cs)` / `calcRefined(cs)` / `identifyModes(cs)` / `applyRetention(cs, rl)` | The individual pipeline phases. `CalcSingleHS` is a `mutable struct` (MATLAB handle semantics): each step mutates and returns the object. | | see `calcAll` for the usual one-call path |
| ↳ `classify(cs, method=nothing, threshold=nothing)` | (Re)classify with another method/standard, no new Green/Ritz solve (`:refined` solves its spot directions once). | `classify(cs, :ritz, :high)` | `method`: `:green` (default on `:full`) \| `:refined` (half-space only) \| `:ritz` (default and only choice on `:ritzOnly`); `threshold`: `:high` \| `:low` \| a number on the method's score scale (default `:high` for `:green`/`:refined`, `:low` for `:ritz`). `reclassify(cs, threshold)` keeps the method |
| ↳ `isStale(result)` / `refinedSawAgreement(result)` | `isStale`: the result was classified by another rule set than the current one (`classify` then refuses it). `refinedSawAgreement`: the `:green` and `:refined` entries of each angle, paired (after `classify(result, :refined)`). | | |
| ↳ `greenComponent(cs, "G33")` | \|G\| of any component on `greenVelAxis` — only `"G13eff"` is stored (`greenAbsG13Eff`), any other `"GIJ"`/`"GIJeff"` is one Green solve. | | |
| ↳ `calcAll(cs, mode=nothing, retentionLevel=nothing)` | Runs the full pipeline for this angle and classifies modes. | `calcAll(cs::CalcSingleHS, mode=nothing, retentionLevel=nothing)` | positional, not keyword — mirrors MATLAB's `.calcAll(mode, retentionLevel)` |
| ↳ `calcSpotGreen(cs, kernel; DirStep, RelKStep, Eta, VelStep, useParallel=:auto)` | *EXPERIMENTAL, as in MATLAB (validated only against the published Fe3Al(316) map, UTGS paper arXiv:2510.10696, Figs. 3–4):* the spot-averaged \|G\| curve — the coherent sum of the complex Green's function over the wavevectors a `ProbeKernel` couples to. Fills `spotGreen` (`Float32`), `spotGreenVelAxis`, `spotGreenInfo`; never read by a classification; kept by `calcAll`/`classify`. Equal bit for bit to this angle's column of `calcSpotGreenMap`. Plot with `plotGreenCurve(cs; component="spot")`. | `calcSpotGreen(cs, probe_gaussian(0.004, 0.3); DirStep=0.25)` | `kernel`: `probe_point()` \| `probe_gaussian(sigmaRelK, sigmaAngDeg)` \| `probe_tgs(...)` \| `probe_sras(...)` \| `probe_from_real_space(...)` \| `ProbeKernel(weightFcn, ...)`. Keywords as MATLAB's Name-Value pairs: `DirStep` [deg] (default kernel's), `RelKStep` (default `min(kernel's, Eta/4)`), `Eta` (default `SawSpectrum.ETA`), `VelStep` [km/s] (default 0.004). Half-space only (`ArgumentError` on a layered structure); any `calcMode` |
| `angSweepHS` → `CalcAngSweepHS` | Angular sweep of `singleHS`, sharing one Green velocity axis. | `angSweepHS(; structure=nothing, angles=0:2:180, calcMode=:full, retentionLevel=:full, nVelGreen=nothing, useParallel=:auto)` | `angles` default `0:2:180`; `useParallel`: `:auto`\|`:always`\|`:never` (default `:auto`, parallel ≥20 angles — `pitfallsJ.md` jNUM5, `matlabVsJuliaJ.md` §4) |
| ↳ `calc(sweep, mode=nothing, retentionLevel=nothing)` | Runs the sweep; populates `.cs` (per-angle result), `.sawVelMap`, `.sawTypeMap`, `.significantModesMap`, `.greenVelAxis`, `.greenAbsG13EffMap` (Float32). `classify(sweep, method, threshold)` / `greenComponentMap(sweep, "G33")` as for one angle. | `calc(sweep::CalcAngSweepHS, mode=nothing, retentionLevel=nothing)` | positional, not keyword, same reason as `calcAll` above |
| ↳ `calcSpotGreenMap(sweep, kernel; DirStep, RelKStep, Eta, VelStep, useParallel=sweep.useParallel)` | *EXPERIMENTAL, as in MATLAB:* the spot-averaged \|G\| map (`spotGreenMap` `[nVel x nAngles]` `Float32`, `spotGreenVelAxis`, `spotGreenInfo` with the per-angle quadrature checks `dirConvergence`/`relKConvergence`). One Green solve per distinct direction, shared between angles; `SpotGreen:notConverged` warns when the quadrature is too coarse (the "comb"). Kept by `calc`/`classify`. Plot with `plotGreenMap(sweep; component="spot")`. Lower level: `SpotGreen.options(kernel; ...)`, `SpotGreen.compute(structure, angles, kernel, opts; useParallel)` `-> (absX, vOut, info)`. | `calcSpotGreenMap(sweep, probe_sras(10, 8); DirStep=2.0, Eta=1e-2)` | as `calcSpotGreen`; an SRAS kernel is one-sided — compute 0–359° |

## Layered forward calculations (`src/forwardCalcs/`)

Film(s) on a semi-infinite substrate (`Structure` with `nlay(structure) >= 2`, layer 1 =
substrate). `vMax` is optional: the Green axis ends at 1.05 × the fastest longitudinal limit of
any layer; a larger `vMax` only extends the displayed axis.

| Name | Purpose | Call | Key inputs |
|---|---|---|---|
| `singleLayered` → `CalcSingleLayered` | Single-angle layered container. | `singleLayered(; structure=nothing, angle=0, vMax=nothing, precompStruct=NamedTuple(), calcMode=:full, retentionLevel=:full)` | `vMax`: optional (extends the Green axis); classification `:green` \| `:ritz` (no `:refined` for layered stacks) |
| `angSweepLayered` → `CalcAngSweepLayered` | Angular sweep at a fixed layer stack (isoWavevec family). | `angSweepLayered(; structure=nothing, angles=nothing, vMax=nothing, calcMode=:full, retentionLevel=:full, rawByDepthCell=nothing, useParallel=:auto)` | as for `singleLayered`; `rawByDepthCell` (opt): a pre-built `rawByDepth` cell, as the 2D sweeps pass it |
| `wavevecAngSweep2Layered` → `CalcWavevecAngSweep2Layered` | 2D sweep: angle × real wavelength (film thickness fixed in real units via `ApplyWavelength`, applied internally per wavelength). | `wavevecAngSweep2Layered(; structure=nothing, angles=nothing, wavelengths=nothing, vMax=nothing, calcMode=:full, retentionLevel=:none, useParallel=:auto, sawClassMethod=nothing, sawThreshold=nothing)` | `wavelengths`: required; every film layer's `CalcParams.realThickness` must be set; default `retentionLevel=:none` (RAM multiplies across both axes) |
| `isofreqSingleLayered` → `CalcIsofreqSingleLayered` | Single-angle, single-frequency Green scan (fixed real `omega0`, sweeps k∥). | `isofreqSingleLayered(; structure=nothing, angle=nothing, omega0=nothing, kParRange=nothing)` | `omega0`, `kParRange`: required. No Ritz, no eigenvector retention. |
| `isofreqAngSweepLayered` → `CalcIsofreqAngSweepLayered` | Angular sweep at fixed `omega0` (isoFreq family). | `isofreqAngSweepLayered(; structure=nothing, angles=nothing, omega0=nothing, kParRange=nothing, useParallel=:auto)` | `omega0`, `kParRange`: required |
| `freqAngSweep2Layered` → `CalcFreqAngSweep2Layered` | 2D sweep: angle × real frequency (isoFreq counterpart of `wavevecAngSweep2Layered`). | `freqAngSweep2Layered(; structure=nothing, angles=nothing, omega0List=nothing, kParRange=nothing, useParallel=:auto)` | `omega0List`, `kParRange`: required |

isoWavevec (Ritz+Green) and isoFreq (Green-only) are two different scan families, unchanged
from MATLAB — see `docs/background.md` for which physical question each answers.

### Classifier and unit helper (all families)

| Name | Purpose | Call | Key inputs |
|---|---|---|---|
| `modeClassification` → `ModeClassification` | The classifier the Calc classes run; for custom post-processing. | `modeClassification(; threshold=:high)` | `calc(mc, ritz, structure, christoffelResult, greenResult, ritzScore)`; output `significantModes::Vector{SignificantMode}` (`type, vel, modeNo, modeVel, greenPeakVel, method, score, scoreDb, grade, ...`). `vel` is the deciding method's velocity (the Green / spot peak on `:green` / `:refined`, the mode's eigenvelocity on `:ritz` and limiting entries); `modeVel` is always the bound Ritz mode's eigenvelocity — use it wherever the Ritz model's own velocity is meant (e.g. synthetic inversion data) |
| `freqToOmega0` | Converts a real frequency to `GreenSolver`'s normalised isoFreq `omega0`. | `freqToOmega0(freqHz, lambdaRefNm::Real)` | positional, both required; `lambdaRefNm` must match the reference used in `ApplyWavelength` |

## Free plate (Lamb waves) forward calculations (`src/forwardCalcs/`)

Homogeneous or multi-layered free-standing plate, both surfaces traction-free
(`nlay(structure) >= 1`). No `calcMode` argument anywhere in this family (no Green's-function
analogue for a free plate, same as MATLAB); `retentionLevel` is `:full`/`:none` only (no
`:saw_only` tier).

| Name | Purpose | Call | Key inputs |
|---|---|---|---|
| `singlePlate` → `CalcSinglePlate` | Single-angle, single-wavevector container. | `singlePlate(; structure=nothing, angle=nothing, realThickness=nothing, realWavVec=nothing, precompStruct=NamedTuple(), retentionLevel=:full)` | `realThickness`/`realWavVec` (both opt, default = normalised units only) |
| `angSweepPlate` → `CalcAngSweepPlate` | Angular sweep at fixed thickness/wavevector. | `angSweepPlate(; structure=nothing, angles=nothing, realThickness=nothing, realWavVec=nothing, retentionLevel=:full, useParallel=:auto, precompStruct=NamedTuple())` | dense maps (every angle has the same mode count, not ragged); `precompStruct` (opt): a shared precompute |
| `wavevecAngSweep2Plate` → `CalcWavevecAngSweep2Plate` | 2D sweep: angle × real wavevector (thickness fixed). | `wavevecAngSweep2Plate(; structure=nothing, angles=nothing, wavVecs=nothing, realThickness=nothing, retentionLevel=:none, useParallel=:auto)` | `wavVecs`, `realThickness`: both required |
| `locateZGVs(sweep, angleIdx=1, maxCount=...)` | Points of zero group velocity along k at one fixed angle: candidate scan, then Newton refinement (`ZgvNewton`, no mode number held fixed, adapted from GEWtool's `ZGVNewtonBeta`), kept where an independent solve confirms `|vg| <= 1e-4 vp`; each flagged a true ZGV point (transverse component zero too: symmetry directions, e.g. ⟨100⟩, ⟨110⟩ of a cubic (001) plate) or a TGV wave (tutorial S2, Part 5). | `locateZGVs(sweep::CalcWavevecAngSweep2Plate, angleIdx=1, ...)` | needs `retentionLevel=:full`; returns `Vector{ZGVPoint}` (`branchNo`, `wavVecZGV`, `freqZGV`, `angle`, `velGroupTrans`, `isTrueZGV`) — the same list as MATLAB |
| `classifyBranches(sweep, angleIdx=1, confidenceThresh=0.3)` | Multi-layer-aware mode-type classification with a confidence-driven fallback to the exact parity class or branch identity. | `classifyBranches(sweep::CalcWavevecAngSweep2Plate, angleIdx=1, confidenceThresh=0.3)` | needs `retentionLevel=:full` |
| `reindexByBranch(map3D, angleIdx, ...)` | Re-orders a raw `(mode, angle, wavevec)` map by tracked branch identity. | | internal helper exposed for custom post-processing |
| `isMirrorSymmetric(structure, tol=1e-10)` | Exact test for the mid-plane mirror (each layer's tensor vs the mirror image of its partner's) — when true, every Plate mode's parity is exact. | `isMirrorSymmetric(obj::Structure, tol=1e-10)` | `tol` relative; a single Cu(111) layer is false |
| `isSagittalMirrorSymmetric(structure, angleDeg, tol=1e-10)` | Exact test for the sagittal mirror `x2 → −x2` at one angle — with `isMirrorSymmetric`, all four Plate labels are exact. | `isSagittalMirrorSymmetric(obj::Structure, angleDeg, tol=1e-10)` | angle-dependent |

## Visualisation (`src/visual/` data prep + `ext/GuliwerrJMakieExt.jl` drawing)

Every `plot*` function is declared as an empty stub in the core package (`src/visual/PlotStubs.jl`)
and only gains a real method once a Makie backend is loaded (`using GuliwerrJ, CairoMakie` or
`GLMakie`) — calling one before that throws an ordinary `MethodError`, not a custom message.

| Name | Purpose | Call | Key inputs |
|---|---|---|---|
| `plotGreenMap` | \|G13_eff\| heatmap, angle × velocity (any other component recomputed; `component="spot"`: the stored spot-averaged map of `calcSpotGreenMap`, on its own velocity axis, any `calcMode`). | `plotGreenMap(sweepResult, scaleMode=:log, climRange=(-4.0,-1.0), plotType=:cartesian; ax=nothing, component="G13eff")` | `sweepResult`: `calc()`'d `angSweepHS`/`angSweepLayered` result with `calcMode=:full`; `scaleMode`: `:log`\|`:linear`; `plotType`: `:cartesian`\|`:polar` |
| `plotRitzSpectrumMap` / `plotRitzSpectrumCurve` | The Ritz spectrum \|S(v)\| from the stored modal numerators (no solve); the default background of a `:ritzOnly` result. | `plotRitzSpectrumMap(sweepResult, scaleMode=:log, climRange=(-4.0,-1.0), plotType=:cartesian; ax=nothing, showCands=false)`, `plotRitzSpectrumCurve(singleResult, scaleMode=:log; ax=nothing, showCands=false)` | either `calcMode` |
| `plotDispersionMap` | 2D/1D dispersion map (angle × wavelength/kParRange/omega0); dispatches on the sweep-result type. | `plotDispersionMap(sweepResult, omegaIdx, scaleMode, climRange, plotType; ax=nothing)` | `sweepResult`: any `wavevecAngSweep2Layered`/`freqAngSweep2Layered`/`isofreqAngSweepLayered`/`isofreqSingleLayered` result |
| `plotDispersionCurves` | Fixed-angle k-/wavelength-dispersion curves over tracked branches, optionally over a Green's-\|G\| background. | `plotDispersionCurves(sweepResult, curveType, angleIdx, showRitzBranches, showGreenBackground, scaleMode, climRange; ax=nothing)` | `curveType` one of 6 `Symbol`s, e.g. `:wavelength_vphase` — **underscore, not hyphen**, see `pitfallsJ.md` jOTH3 |
| `plotDispersionCurvesPlate` | The `wavevecAngSweep2Plate` counterpart, coloured by mode type or branch, plus Christoffel bound overlays. | `plotDispersionCurvesPlate(sweepResult, curveType, angleIdx, colorMode, confidenceThresh, showChristoffelBounds; ax=nothing)` | same 6 `curveType` names |
| `plotClassificationDiagnostic` | Classified entries of the active classification over the G13_eff map (Ritz spectrum map on `:ritzOnly`), marker size by score on the method's scale, `'weak'` entries hollow, plus the limiting-velocity reference line. | `plotClassificationDiagnostic(sweepResult, colorBy=:type, showGreenBackground=true, refLevel=:full; ax=nothing, component=nothing)` | `colorBy`: `:type`\|`:score` (scoreDb)\|`:confidence`; `refLevel`: `:sawOnly`\|`:sawLimit`\|`:full`; `component`: `"G13eff"`/`"GIJ"`/`"ritzSpectrum"` |
| `plotClassificationCurve` | The single-angle counterpart of `plotClassificationDiagnostic`. | `plotClassificationCurve(singleResult, colorBy, showGreenBackground, refLevel; ax=nothing, component=nothing)` | |
| `plotAngularMap` | SAW velocity vs. angle. | `plotAngularMap(sweepResult, plotType; ax=nothing)` | `plotType`: `:cartesian`\|`:polar` |
| `plotPolarHeatmap` | Generic angle×radius heatmap drawn polar-style (`plotGreenMap`/`plotDispersionMap`'s `plotType=:polar` delegates here). | `plotPolarHeatmap(angles, radius, C; ax=nothing, colormap=:viridis, colorrange=nothing)` | calls `pcolorCenteredGrid` internally — see `pitfallsJ.md` jNUM3 |
| `plotModalShape` | A Ritz mode's depth×lateral displacement field. | `plotModalShape(ritzSolver, modeNo, substrateCutoff, cutoffParam, dirs, meshSpacing; fig=nothing)` | `ritzSolver`: from `singleResult.ritzSolver`; `dirs`: `1:3` or `:rgb` (one polarisation-RGB panel + key, also on the two Green variants) |
| `plotPolarisationKey` | Colour key of the polarisation-RGB view. | `plotPolarisationKey(; ax=nothing)` | |
| `plotGreenModalShape` | A Green's-response displacement field (half-space only). | `plotGreenModalShape(greenObj, forceDir, cutoffParam, dirs, meshSpacing; fig=nothing)` | |
| `plotGreenModalShapeLayered` | The layered counterpart, one panel per layer. | `plotGreenModalShapeLayered(greenObj, structure, forceDir, substrateCutoff, cutoffParam, dirs, meshSpacing; fig=nothing)` | `nlay(structure) >= 2` |
| `plotGroupVelSurface` | 3D group-velocity ("ray") surface. | `plotGroupVelSurface(mater, nTheta, nPhi; fig=nothing)` | |
| `plotVelSurface` | 3D directional bulk-wave velocity/slowness surface. | `plotVelSurface(mater, surfaceType, nTheta, nPhi; fig=nothing)` | |
| `plotYoungsSurface` | 3D directional Young's-modulus surface. | `plotYoungsSurface(mater, nTheta, nPhi; fig=nothing)` | |
| `plotVelCurve` | 2D bulk-wave phase-velocity vs. angle curve. | `plotVelCurve(mater, angleDeg; ax=nothing)` | |
| `plotYoungsCurve` | The `plotVelCurve` counterpart for Young's modulus. | `plotYoungsCurve(mater, angleDeg; ax=nothing)` | |
| `plotRitzModeMap` / `plotRitzModeCurve` | Scatter/curve of every computed Ritz mode, shaded by a Ritz-mode criterion. | `plotRitzModeMap(sweepResult, criterion=:ritzScore, velCutoffFactor=1.10; ax=nothing)` | `criterion`: `:ritzScore` (per-mode Ritz score, dB) \| `:freqDisp` (with the half-space veto) \| `:depthInsensitivity` |
| `plotRitzModeMapPlate` / `plotRitzModeCurvePlate` | The Plate counterparts, shaded by mode type or confidence. | `plotRitzModeMapPlate(sweepResult, criterion, velCutoffFactor; ax=nothing)` | `criterion`: `:modeType`\|`:confidence` |
| `plotGreenCurve` | \|G13_eff\| (or any component) vs. velocity for one angle (`component="spot"`: the stored `calcSpotGreen` curve). | `plotGreenCurve(singleResult, scaleMode=:log; ax=nothing, component="G13eff")` | |
| `plotAssociationDiagnostic` | `ModeAssociation` points diagnostic. | `plotAssociationDiagnostic(association, overlayMode; ax=nothing)` | |
| `plotInversionFit` | Measured-vs-fitted velocity comparison, before/after inversion. | `plotInversionFit(measured, associationFinal, associationInitial; ax=nothing)` | |
| `plotSensitivitySpectrum` | Sensitivity eigenvalue spectrum + transform-matrix heatmap. | `plotSensitivitySpectrum(sensitivityAnalysis; fig=nothing)` | |

Data-prep helpers feeding the above (not usually called directly, no Makie dependency —
`src/visual/VisData.jl`/`PlotData.jl`): `calcDispField2D`, `calcGreenDispField2D`,
`calcGreenDispField2DLayered`, `calcSubstrateCutoffDepth`, `extractPlanarCut`,
`scaleGreenData(absG, scaleMode, componentLabel)` (label required, e.g. `"G_{13}^{eff}"`),
`greenComponentData(result, component="G13eff") -> (absG, label, velAxis)` (`"spot"`: the stored spot-averaged curve/map, label `"G_{13}^{spot}"`, own axis),
`ritzModeScores(singleResult)` (per-mode Ritz score), `classificationScoreScale(method,
isLayered)` (thresholds and marker-size range of a method), `calcPolarisationRgb(dispField)`,
`calcModalShapeRowHeights`, `pcolorCenteredGrid`, `ColormapRWB`, `ColMapGrayOrange`,
`ColormapRYG`, `plateModeTypeStyle`, plus the `*Data` functions (`velSurfaceData`,
`groupVelSurfaceData`, `youngsSurfaceData`, `velCurveData`, `youngsCurveData`,
`ritzSpectrumMapData`, `ritzSpectrumCurveData`, `dispersionCurvesData`,
`wavevecGreenBackgroundData`) that compute exactly what each `plot*` function above draws.

## GUI (`ext/Gui*.jl`, package extension)

**Experimental** — see `quickStartJ.md`, "The GUI layer"; the MATLAB viewers are the ones to use for interactive work.

Exported from `GuliwerrJ` itself (empty stubs in `src/visual/GuiStubs.jl`, gaining their real
implementation from the extension — the same pattern the `plot*` functions above use), so
every name below is callable directly after `using GuliwerrJ, GLMakie` (or `CairoMakie`) — see
`pitfallsJ.md` jOTH4. Each launcher accepts a `fig` keyword (`nothing` for a standalone window, or a `GridPosition` to embed into).

| Launcher | Opens | Purpose |
|---|---|---|
| `structureControl(; fig=nothing)` | Structure editor + job launcher | Build a `Structure` layer-by-layer and launch a forward calculation on it, GUI-driven. Returns a `NamedTuple` (`fig`, `result::Observable`, `statusLabel`, `layers::Observable`, `widgets`, `run!`) — `result` is the hand-off point a viewer watches. |
| `singleAngleViewer(result; fig=nothing)` | Single-angle curve view | For a `calcAll()`'d `CalcSingleHS`/`CalcSingleLayered`/`CalcSinglePlate`/`CalcIsofreqSingleLayered`. Four graph types in MATLAB's order (classification / Green's-function / Ritz-mode curve / Ritz spectrum). |
| `angSweepViewer(sweepResult; fig=nothing)` | Angular-sweep 2D map | For a `calc()`'d `CalcAngSweepHS`/`CalcAngSweepLayered`/`CalcAngSweepPlate`/`CalcIsofreqAngSweepLayered`. Click a point to select a mode; drills down to a single-angle curve in place. |
| `modalShapeViewer(result; fig=nothing)` | Modal-shape depth-panel view | Accepts any single or sweep result from the families above directly. Always opens its own separate `Figure` for the shape panels (a Green's-response toggle needs the same Figure-clearing machinery `plotGreenModalShape` uses). |
| `combinedViewer(result; windowMode=:combined, fig=nothing)` | Two-pane shell | `modalShapeViewer` (top) wired to `angSweepViewer`/`singleAngleViewer` (bottom, dispatched on `result`'s type), bidirectionally linked. `windowMode`: `:combined` (default) or `:split`. |
| `doubleSweepViewer(sweepResult; fig=nothing)` | Outer-axis slice stepper | For a genuine double sweep (`CalcWavevecAngSweep2Layered`/`...Plate` with `length(angles) > 1` and `length(wavelengths)`/`length(wavVecs)` `> 1`, or `CalcFreqAngSweep2Layered`). Steps an outer wavelength/wavevector/omega0 value, rebuilding an embedded `combinedViewer` at each step; for the Layered/Plate cases a button switches to an embedded `dispersionViewer` ("dispersion mode") — `CalcFreqAngSweep2Layered` has no second mode, matching its own MATLAB source. |
| `dispersionViewer(sweepResult, angleIdx=1; fig=nothing)` | Modal shapes + dispersion-curve map | For a `CalcWavevecAngSweep2Layered`/`...Plate` result with `retentionLevel = :full` (needed by the internal `DispersionBranches.track` call). `sweepResult` may carry more than one angle; `angleIdx` picks which. |
| `setAxisRange!(viewer, axis, lims)` | Set an axis range from a script | `axis` is `:x`, `:y` or `:c` (colour, map viewers only); `lims` is `(lo, hi)`, a `NaN` entry keeps that limit, `[]` keeps both. Fills the min/max fields as typing would. Works on `angSweepViewer`, `singleAngleViewer`, `dispersionViewer`, and forwards from `combinedViewer` (to the map) and `doubleSweepViewer` (to the view of the current mode). As in MATLAB, `doubleSweepViewer` keeps each mode's ranges across a slice step and a mode switch and back. |

## Solvers (`src/solvers/`)

The solver-level exports, for custom post-processing and for the tutorials that call them
directly. Most scripts never need them: the `Calc*` containers above run them.

| Name | Purpose | Call |
|---|---|---|
| `phaseVels` / `phaseVelScan` | Christoffel bulk phase velocities for one direction / for a material over in-plane angle. | `phaseVels(elTensor, dirVec, rho=1.0)`, `phaseVelScan(mater::Material, angleDeg=0.0)` |
| `groupVelVectorsForDir` | Bulk group-velocity vectors for one direction. | `groupVelVectorsForDir(elTensor, dirVec, rho=1.0)` |
| `calc_anisotropy_ratio` | Generalised anisotropy ratio, `(fastest shear / slowest shear)^2` over all directions, of a density-normalised elastic tensor (MATLAB's `ElastConsts.CalcAnisotropyRatio`); prefer `ZenerRatio` for cubic (Zener ratio `2·C44/(C11−C12)`). | `calc_anisotropy_ratio(normElTen, angleStep=4.0)` |
| `LimitingRoot`, `findRoots`, `refineRoot` | The limiting (bulk-wave) velocities of a cut. | `findRoots(rotElTensor, rho=1.0, alphaPrecision=1e-2, relTol=1e-3)` |
| `RitzSolverSAW`, `RitzSolverLamb` | The Ritz solvers for a half-space/layered surface wave and for a free plate. A plate solve also gives `velsGroupTrans`, the group velocity perpendicular to k (positive towards increasing angle). | `RitzSolverSAW(; structure, angleDeg)`; see `?RitzSolverLamb` |
| `locateZGV`, `kPolynomial`, `ZgvNewton.solve` | Plate ZGV point nearest a seed on one raw mode (Newton; `NaN` if not converged); the plate eigenproblem as a polynomial in k; the Newton iteration itself. | `locateZGV(structure, angleDeg, realThickness, modeNo, wavVec0, precomp=NamedTuple())`, `kPolynomial(structure, angleDeg, realThickness)`, `ZgvNewton.solve(L2, L1, L0, M, k0, omega0, u0=nothing)` |
| `calcWithGamma` | `calc` of a `RitzSolverSAW` that also returns the assembled `GammaMatrix` (K, dK/dD, M) it used — MATLAB's `gammaMat` property. | `(s, gm) = calcWithGamma(s, withDiagnostics=true, precomp=NamedTuple())` |
| `ritzSAWPrecompute`, `ritzWavelengthSweep` → `RitzWavelengthSweep` | The angle-independent Ritz precompute of a structure; for a sweep over wavelengths (`ApplyWavelength`) build the wavelength-independent part once (raw Legendre tensor, constraint basis, gauge map, Gram blocks) and complete it per wavelength (about 10× cheaper on a four-layer, legOrder-50 stack; bit-identical). | `ritzSAWPrecompute(structure)`; `ws = ritzWavelengthSweep(structure)`, then `ritzSAWPrecompute(ApplyWavelength(structure, λ), ws)` (same Legendre orders required) |
| `ritzSensitivities` | **Experimental.** The lowest Ritz modes of a layered stack and the exact (Hellmann–Feynman) log-derivatives of their phase velocities with respect to every layer's elastic constants, density and depth, the angle and `|k|` (`vgOverVp`) — one solve, no finite differences. | `ritzSensitivities(structure, angleDeg, precomp=NamedTuple(); nModes, vMax=Inf)`; fields `vels, eigVecs, dlnv_dlnC, dlnv_dlnRho, dlnv_dlnD, dlnv_dTheta, vgOverVp, kinShare` |
| `RitzWorkspace`, `ritzSAWPrecompute!`, `ritzSensitivities!` | **Julia only.** Preallocated arrays for repeated lowest-mode Ritz solves of one stack family (fixed Legendre orders; any wavelength, constants, densities): the in-place precompute and sensitivities give bit-identical results without the few MB of scratch allocation per call. `!` outputs of the precompute live in the workspace until its next call; the sensitivities are fresh copies. One workspace per task — see `forEachElementPooled`, `pitfallsJ.md` jNUM2. | `ws = RitzWorkspace(sweep; nModes)` (or `(structure; nModes)`); `pre = ritzSAWPrecompute!(ws, ApplyWavelength(structure, λ), sweep)`; `ritzSensitivities!(ws, structure, angleDeg, pre; vMax=Inf)` |
| `forEachElementPooled` | **Julia only.** `forEachElement` with a reusable workspace per task: a pool (`Channel`) of `min(threads, n)` workspaces, each element `take!`s one and returns it; BLAS pinned as in `forEachElement`. | `forEachElementPooled(items, useP, () -> RitzWorkspace(sweep; nModes = 8), (x, i, ws) -> ...)` |
| `solveConstrainedEig` | The Ritz eigensolve; `lowest = true` computes only the lowest `nev` pairs (accurate bisection, pays for tens of modes at most). | `solveConstrainedEig(K_con, M_con=nothing; nev=nothing, lowest=false, cholM=nothing)` |
| `GreenSolver` | The Green's function: a spectrum over a velocity axis, optionally with the displacement field. | `calc(GreenSolver, structure, angle, velRange; forceDir=1, displDir=3, fullOutput=false, ...)` |

## Inverse characterisation (`src/inverseCalcs/`, `Inv.jl` launchers)

> **Experimental.** The inverse module is validated on synthetic data only; the interface may
> change. Its pseudo-surface-wave tracking is new.

Recover the elastic constants of one selected layer from measured SAW/pseudo-SAW/bulk-limiting
velocities from one angular sweep, plus optional pulse-echo phase velocities and a soft-fixed
bulk modulus — the same design as in MATLAB, including the object-array joint
multi-orientation inversion (`structure`/`association`/`measured` each accept either one object
or a `Vector` of them, one per crystallographic cut of the same crystal — a single cut is just
a length-1 vector and behaves identically).

| Name | Purpose | Call | Key inputs |
|---|---|---|---|
| `measured` → `MeasuredVelocities` | The measured dataset: angle/velocity/waveType/std per point. | `measured(angleDeg, velocity, waveType, velocityStd)`; `fromTable(cols::NamedTuple)` builds one from columns `angleDeg`, `velocity` and optionally `waveType`, `velocityStd` (MATLAB's `MeasuredVelocities.fromTable`) | `waveType`: member of `WAVE_TYPES` (`"rayleigh"`, `"pseudoRayleigh"`, `"sezawa"`, `"love"`, `"longitudinal"`, `"fastShear"`, `"slowShear"`, `"unknown"` — still plain strings here, not Symbols, matching `MeasuredVelocities`'s own field type) |
| ↳ `addPhaseVelocity(obj, propagationDir, waveType, velocity)` | Adds an independently measured pulse-echo bulk phase velocity. | | `propagationDir` in the inverted layer's **crystal frame** |
| ↳ `setBulkModulus(obj, bulkModulusGPa)` | Sets the soft-fixed bulk-modulus pseudo-datum. | | entered as an artificial `sqrt(B/rho)` pseudo-velocity |
| `association` → `ModeAssociation` (via `fromSweep`) | Builds the initial measured-point ↔ calculated counterpart binding from a forward sweep at the initial guess. | `association(sweep, measured, velocityTolerance=0.3, angleTolerance=1e-6, invertedLayerNo=1)` | `sweep`: `angSweepHS`/`angSweepLayered` result, `retentionLevel=:full`, `calcMode=:full` |
| ↳ `rebind(assoc, angleDeg, eigVecs_con, M_con, vels, ...)` | Re-establishes the `ritzMode`/`pseudoSAW` binding at one angle from a freshly solved eigenproblem. | | MAC-based, phase- and basis-invariant — safe against the null-space non-reproducibility in `pitfallsJ.md` jMTH1 |
| ↳ `rebindLimiting(assoc, angleDeg, roots, ...)` | The `limiting`-class counterpart to `rebind`. | | `roots` from `limitingRootsAt` |
| ↳ `setPinned(assoc, idx, tf)` / `pinnedPoints(assoc)` | Marks a point as human-confirmed. | | |
| ↳ `confidence(assoc)` | Per-point confidence signal. | | |
| `sensitivity` → `SensitivityAnalysis` (via `calc`) | The sensitivity/Jacobian matrix and its eigendecomposition into well-determined linear combinations of elastic constants. | `sensitivity(structure, association, invertedLayerNo, measured; pointWeights=nothing, includePseudoSAW=true, useParallel=:auto, sensitivityCutoff=1/160)` | `structure`/`association`/`measured` each accept a single object or a `Vector` (joint inversion) |
| ↳ `limitingRootsAt(matLay, symmNo, cVec, angleDeg)` | The labelled Christoffel root set `rebindLimiting` re-identifies against. | | |
| `elasticity` → `ElasticityInversion` | The bounded-least-squares optimisation orchestrator, over the well-determined transformed constants. | `elasticity(structure, measured, association, sensitivity, invertedLayerNo; anchorPolicy=:chained, rebindVelocityWindow=1.0, includePseudoSAW=nothing, stepwise=true, pseudoSAWRefinement=:green, lowerBound=Float64[], upperBound=Float64[], ...)` | `includePseudoSAW` defaults to the sensitivity's. The optimiser is a projected Levenberg–Marquardt (Gauss–Newton family, like MATLAB's `lsqnonlin`) with the box bounds applied by projection; the iterate path differs from MATLAB's — `matlabVsJuliaJ.md` §3.1 |
| ↳ `calc!(inv)` | Runs the inversion **in place** (note the `!`; like all `Calc*` results, `ElasticityInversion` is a mutable struct, and it carries best-so-far optimiser state across evaluations). Populates `inv.sensitivity`/updates `inv.structure`/`inv.association`. | `calc!(inv::ElasticityInversion)` | never re-runs the forward sweep — exploits K's linearity in the elastic constants |
| ↳ `recalcForwardAtResult(inv; angles=nothing, vMax=nothing, calcMode=:full, retentionLevel=:full)` | Fresh forward sweep at the inverted elasticity, for display/checking. | | all keywords; deliberately separate and opt-in |

## Where to look next

- New calculation script: start from `juliaPort/tutorialsJ/`, not from scratch.
- Full argument list for any function above: its own Julia docstring (`?name` at the REPL).
- Julia-specific traps: `pitfallsJ.md`; every difference from MATLAB: `matlabVsJuliaJ.md`.
- What the calculation actually means physically: `docs/background.md`.
