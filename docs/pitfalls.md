# GuliweRR — pitfalls

Short warnings for using the published API. Read this in full before writing a new calculation script — every item here has caused a real bug at least once.

The items are grouped by kind, most fundamental first, and ordered by severity within each group:

- **PHY — Physics:** what the waves, peaks and scores mean.
- **MTH — Method-related numerics:** the Ritz and Green solvers and the mode classification.
- **NUM — Standard numerics:** MATLAB, memory and parallel runs.
- **OTH — Other specifics:** coordinates, setup, saved results, inversion and plotting.

Severity: 🔴 **silent wrong result** (the script runs and the numbers look plausible) · 🟠 **misreading** (the result is right, its meaning easy to get wrong) · ⚪ **loud or cost only** (an error, a warning, RAM or run time).

The IDs (`PHY1`, `MTH5`, …) are stable: a new item gets the next free number in its group, and nothing is renumbered.

## Physics

### PHY1 🟠 Existence is a threshold on a continuum — `type` and a branch end are labels, not ground truth

A surface wave is accepted if its existence score reaches the standard of the active classification (`activeMethod`: `'green'`, `'refined'`, or `'ritz'`; standard `activeThreshold`). The physics is continuous — a branch fades gradually, and the Rayleigh→pseudoRayleigh transition is gradual too — so `type` is a label of convenience and the angle where a branch "ends" depends on the standard you choose. In measured maps the end of a fading branch is likewise uncertain by a few degrees. If a branch end matters to you, compare the candidates' scores with the standard across the angles either side (`classification.greenScorePeaks` / `ritzScorePeaks`, and `refinedSaw.peaks` per angle for `'refined'`, list every candidate, accepted or not), or `reclassify` at the other standard; don't chase a sharper end with finer velocity resolution, which does not change it. Each entry's `grade` tells the two standards apart without reclassifying: `'strong'` passes the method's `'low'` standard, `'weak'` passes only the active standard.

What the two standards are for, what the probe-aperture model does to steep branches, and why no threshold should be taken relative to `max(abs(G))` if you write your own classification: [forwardCalcs_modeClassification.md §1 and §2](reference/architecture/forwardCalcs_modeClassification.md#1-overview) (sections "Standards and grades" and "Why each choice").

### PHY2 🔴 Layered stacks: several branches coexist, and a substrate's second shear limit can look like a SAW

The classification is computed independently at *each* point of a sweep — neither the Green-peak matching of `'green'` nor the modal-sum peaks of `'ritz'` have any memory of which physical branch was picked at a neighbouring point. In a layered sweep, higher-order (Sezawa-type) branches can coexist with the fundamental one, so a script that just grabs `sawVelMap{iAngle}(1)` (or any one entry of the classification) at each independent relDepth/wavelength point can silently swap between them, producing a spurious, non-physical peak or jump in the v_SAW curve — or lock onto a substrate bulk-like resonance rather than a surface mode at a single point. Neither `calcMode='full'` nor `'ritzOnly'` removes this: both are evaluated independently per point. How to follow one branch: MTH11.

**Do not expect a layered entry's `type` field to tell you which candidate is "the" fundamental.** Every accepted layered candidate is a single `'SAW'` type — there is no `'Rayleigh'`/`'Sezawa'` split by velocity rank (why: [forwardCalcs_modeClassification.md §7](reference/architecture/forwardCalcs_modeClassification.md#no-love-type-no-sezawa-split)). Use the `vel`, `score`, `matchScore`, `relPeakHeight`, and `aboveSubstrateST` fields on each entry to judge branch order/dominance/leakiness from real numbers yourself, e.g. via `Vis.plotClassificationDiagnostic` (marker size by score).

An anisotropic substrate can have a *second*, distinct bulk-shear limiting velocity beyond the simple isotropic `sqrt(C44/rho)` (e.g. Cu(001) along [100] has two, ≈2.17 and ≈2.87 km/s, computed with GuliweRR via `CalcSingleLayered.calcLim()`'s `subShearLimVels` for the Cu constants of [quickStart.md](quickStart.md)), and a Green's-function cusp right at that threshold can look like a SAW peak. The existence scores' late-time gates suppress such short-lived bulk features, but nothing explicitly excludes candidates coinciding with the substrate's own `'limiting'`-type velocities, and `DispersionBranches.track` follows Ritz modes whatever their classification. If a tracked or classified branch sits suspiciously flat and close to a substrate limiting velocity across a wide relDepth range, treat it with suspicion rather than as a confirmed physical mode.

### PHY3 🟠 A leaky pseudo-SAW may have no individually credible Ritz mode — a weak match is not a missing peak

A leaky pseudo-SAW, or a mode whose weight is smeared over several near-degenerate Ritz modes, has no reason to correspond to one eigenmode of a non-radiating Ritz basis, so a real Green peak's matched Ritz mode can look weak (low `matchScore`/`oopDisplScore`). That is expected, not a sign the entry is wrong: a low-confidence match is still reported, with the plain best-scoring mode ([§6](reference/architecture/forwardCalcs_modeClassification.md#6-identity)). MTH5 covers which field actually decides whether an entry exists, and PHY1 why a genuinely fading branch is not "fixed" by more resolution.

### PHY4 🟠 A tall Green peak can still be the crest of a broad shoulder, not a genuine resonance

Height alone does not show that a Green feature is an isolated, sharp resonance rather than the crest of a wide bulk shoulder or a Van Hove edge (a sharp edge of the bulk-wave continuum), and such a feature can still be matched to a Ritz mode. The existence scores remove most of them by a late-time gate and a relative-width cap ([§5](reference/architecture/forwardCalcs_modeClassification.md#5-sawspectrum--the-shared-late-time-machinery)). If you are peak-searching Green output yourself, compare peak `width`/`prominence` (e.g. `GreenPeaks.find` on that angle's spectrum): genuine peaks are far narrower than shoulders (the widths are given in the same section).

### PHY5 🟠 A true (non-leaky) Rayleigh peak's height and width depend on the damping `η` — its location does not

A genuine, non-leaky Rayleigh peak in `G` is a pole of the undamped problem, regularised only by `GreenSolver`'s limiting absorption (`ω → ω·√(1+i·η)`, default `η = 3e-4`), so its sampled height and width depend on `η` and on how finely the velocity axis samples it; its *location* does not. Only compare absolute heights/widths — or height/width thresholds of your own — between results computed at the same `η` and similar sampling.

The width follows from `η` directly: `ω → ω·√(1+iη) ≈ ω·(1+iη/2)` moves the pole to a complex velocity `v₀(1+iη/2)`, so near it `|G| ∝ 1/|(v−v₀) + iηv₀/2|`, which falls to half its peak value at `|v−v₀| = √3·ηv₀/2`, i.e. the full width at half maximum (FWHM) of `|G|` is √3·η·v₀ (0.52 % of v₀ at η = 3e-3, ≈ 13 m/s at 2.5 km/s); `|G|²` has FWHM η·v₀. Which `η` the `Calc` classes' stored Green data use: MTH8. How the default η is set: `docs/reference/architecture/green.md`.

### PHY6 🟠 Layered stacks: `'low'` by default — Love waves need `'high'`, leaky waves need `'green'`, and the band comes from every layer

The layered classes (`Calc.singleLayered`, `angSweepLayered`, `wavevecAngSweep2Layered`) offer `'green'` (`GreenScoreLayered`, the default on `'full'`) and `'ritz'`; there is no `'refined'` for stacks. Both methods default to `'low'` (`sawThreshold = []`), unlike the half-space `'green'`.

- **`'low'` keeps sagittal-type waves only** — Rayleigh, the Rayleigh cap, the Sezawa families. A surface-polarisation term ranks Love-type waves below the sagittal ones on purpose (`classification.greenScorePeaks.SH`: 0 sagittal … 1 Love). **For Love waves use `classify('green', 'high')` / `reclassify('high')`, or a number below 0.20.** How many Love points and false positives each standard gives on the measured NiTi/MgO maps (Grabec et al. 2024) is discussed in [§11](reference/architecture/forwardCalcs_modeClassification.md#11-validation).
- **Leaky layered waves (above the substrate's shear limit, radiating into it) are never reported by `'ritz'`** — its trapped gate admits trapped waves only, as it cannot tell a leaky wave from the slab continuum; use `'green'` for them (MTH9).
- **The band comes from every layer** (slowest shear and fastest longitudinal limiting velocity of *any* layer, `bandVSmin` / `bandVLmax` after `calcLim`), not the substrate's, so a fast deep layer widens the search. `aboveSubstrateST` still refers to the substrate's minimum shear limit: a label only. The axis: `docs/reference/architecture/forwardCalcs.md`, "The Green velocity axis".
- **Calibrated on one measured stack** (the NiTi film on MgO of Grabec et al. 2024). A strongly leaky sagittal wave between the substrate's shear and longitudinal limits can pass `'low'`; treat an accepted entry with a large pole width above the substrate shear limit with care until a measured case covers it ([§12](reference/architecture/forwardCalcs_modeClassification.md#12-limitations)).
- `Calc.angSweepLayered` argument 6 and `Calc.wavevecAngSweep2Layered` argument 7 are unused and must be empty: the layered Green axis has a fixed 4 m/s step.

Design and thresholds: [Layered stacks](reference/architecture/forwardCalcs_modeClassification.md#layered-stacks).

## Method-related numerics

### MTH1 🔴 `legOrder` should be ≥ 2× `relDepth` — and both scale with anisotropy/plane

`CalcParams(relDepth, legOrder)`'s `legOrder` (Legendre polynomial order — mesh resolution in depth) should satisfy `legOrder >= 2*relDepth`. Under-resolution fails silently: `legOrder=20`/`relDepth=20` (ratio 1) can give velocity errors of up to a few per cent (computed with GuliweRR), depending on material and cut, that look like legitimate physics. The half-space `'ritz'` classification (the default on `calcMode='ritzOnly'`) needs more: `legOrder >= 2.3*relDepth`, and it warns below that (MTH9). Every row of the table below satisfies it.

Starting points, by the Zener ratio `A = 2*C44/(C11-C12)` (for non-cubic materials the generalised ratio, the squared ratio of the fastest to the slowest shear velocity, Grabec et al. 2024, Sect. 4.3) and whether the surface is a principal plane.

| Case | Zener ratio `A` | Plane | `relDepth` | `legOrder` |
|---|---|---|---|---|
| Easy (incl. isotropic) | < 5 | principal | 20 | 50 |
| Complex | < 10 | non-principal | 30 | 70 |
| Extreme | > 10 | non-principal | 40 | 100 |

A combination the table does not list (a Zener ratio below 5 on a non-principal plane, or between 5 and 10 on a principal one) takes the "complex" row, as does anything unclear or borderline. `legOrder` also has its own absolute-magnitude and RAM-cost tiers, independent of this ratio: 50 (standard), 70 (larger systems), 100 (very large). `legOrder` 100 is the largest practical basis; expect memory and run time to grow accordingly.

Plates (Lamb waves) are the exception where more is not safer: at low `kd` a large `legOrder` loses digits to round-off instead of gaining them. In the multi-layer plates checked against GEWtool (`docs/reference/architecture/plate.md`), the lowest flexural mode at kd = 0.1 moved by 1e-6 to 3e-5 (relative) between `legOrder` 20 and 48 per layer, but by up to 4e-3 at `legOrder` 60; the best overall agreement was at 28 per layer. Choose a plate's `legOrder` by convergence of the highest mode you need, not "as high as affordable".

### MTH2 🔴 Ritz's own lowest mode can land a hair on the wrong side of `minShear` — it's still Rayleigh

Ritz eigenvalues are variational approximations while Christoffel's `minShear` is exact, so near the boundary `ritz.vels(1)` can sit ~0.01% on the wrong side of `minShear` — numerical noise, not physics. If you classify results yourself, never decide Rayleigh vs pseudoRayleigh by comparing a Ritz velocity with `minShear`, and do not "fix" such a comparison with a tolerance: use the matched Green peak's own velocity (`greenPeakVel`), which carries no Ritz discretisation error and is what `type` is decided from on the `'green'` method. The Green peak has a resolution limit of its own, though: the damping η = 3e-3 of the existence solve broadens every pole to a width of about η·v (the FWHM of |G|², √3 times narrower than the FWHM of |G|, PHY5), so a Rayleigh wave hugging `minShear` can put its peak a fraction of a m/s *above* it. The classifier therefore labels a Green peak `'Rayleigh'` below `minShear·(1 + ModeClassification.RAYLEIGH_LABEL_TOL)` (= η); do the same if you label Green peaks yourself. The `'ritz'` and `'refined'` labels use their own boundaries, and the `'ritz'` label never uses the jittery modal-sum peak velocity: [§7](reference/architecture/forwardCalcs_modeClassification.md#the-boundary-per-method), "The boundary per method". The general rule behind this: MTH3.

### MTH3 🔴 Two branch-eligibility rules that disagree by a numerical epsilon can lose a peak entirely — decide branch membership from ONE quantity, not two

The general form of MTH2, for any post-processing you write: if one step decides "which side of a physical threshold are we on" from one solver's number (e.g. a Ritz velocity) and another step from a different solver's number for the same physical thing (e.g. the Green peak velocity), the two will disagree by a small, irreducible epsilon right at the threshold. A peak in that sliver can be mislabelled or — when the two disagree in opposite directions, so each step rejects it — lost entirely. Derive the decision once, from the single most exact quantity available, and reuse it everywhere ([§7](reference/architecture/forwardCalcs_modeClassification.md#7-labels), "One search range, one deciding quantity").

### MTH4 🔴 A physically real Green peak can be invisible to `findpeaks` if its search window clips it

`findpeaks` only registers a maximum that has a lower neighbour on both sides *inside the array you pass it*, so peak-searching a masked sub-array (`absG(mask)`) whose edge falls on a peak's falling side misses that peak whatever threshold you use. This bites hardest at physical thresholds such as `minShear` (the slowest shear limiting velocity, `min(shearLimVels)` of the Christoffel result), where a genuine Rayleigh peak can sit just below the boundary as a sharp spike. In your own peak searches, extend the window past any boundary and filter the results afterwards, as the classifier does ([§7](reference/architecture/forwardCalcs_modeClassification.md#7-labels), "One search range, one deciding quantity"); if a smooth progression shows an isolated missing point near a boundary, inspect the raw samples straddling it before assuming a threshold rejected the peak.

### MTH5 🟠 Existence is decided by the score alone — `oopSurfDispl` and `matchScore` describe the match, not whether the entry exists

On the `'green'` and `'refined'` methods every candidate whose score passes the active standard appears in `significantModes`, however poor its Ritz match; there is no `oopSurfDispl` or `matchScore` floor, and `relPeakHeight` (height over the tallest candidate at that angle) never decides existence either. If a real-looking peak is missing from the classification, look at its score against the standard (`classification.greenScorePeaks`, or `reclassify('high')`), not at `oopDisplScore`; `matchScore` tells you how good the match is, not whether the entry exists. An entry's `vel` is its Green peak (`greenPeakVel`, the spot peak on `'refined'`); its Ritz mode's eigenvelocity is `modeVel`, which need not lie within `velMatchTol` of the peak, and two entries can share a `modeNo`. Use `modeVel` wherever the Ritz model's own velocity is meant (e.g. synthetic data for an inversion). Design: [§6](reference/architecture/forwardCalcs_modeClassification.md#6-identity), "Green decides existence, Ritz only assists identity".

### MTH6 🟠 Numeric thresholds are per method — prefer `'high'` / `'low'`

`sawThreshold`, `classify(method, threshold)` and `reclassify(threshold)` take `'high'`, `'low'`, or a number. A number is read on the scale of the method it is applied to, and is a different standard on every other, so a number tuned on one method means nothing on another; on the half-space `'green'` and `'refined'` methods no number reproduces `'low'`, which lives on a second, structure-relative scale. Use the names unless you are deliberately scanning a threshold. A numeric `sawThreshold` therefore needs `sawClassMethod` set explicitly: with the method unset (the default, MTH12) `calc()` / `calcAll()` refuse it (`…:numericThresholdNeedsMethod`), while `'high'` / `'low'` / `[]` take the resolved method's own values. The values of every method's standards are tabulated in [forwardCalcs_modeClassification.md §1](reference/architecture/forwardCalcs_modeClassification.md#standards-and-grades) ("Standards and grades").

### MTH7 🟠 A score is only comparable with scores and standards of its own method

An entry's `method` names the classification that produced it, and its `score` is on that method's own scale (MTH6; the scales are defined in [§§2–4](reference/architecture/forwardCalcs_modeClassification.md#2-green--greenscore-half-space)). Never compare scores across methods or geometries, or a score from one structure with a structure-relative standard of another. No method needs neighbouring angles: a standalone angle classifies exactly as the same angle inside any sweep.

### MTH8 🟠 The `Calc` classes' Green data are at η = 3e-3 on a coarse axis — not at `GreenSolver`'s default

The `Calc` classes' Green data are **not** at the default `η` (PHY5): in `calcMode='full'` each angle runs one Green solve at `SawSpectrum.ETA = 3e-3` on a uniform axis from 0.3 km/s — half-space in `greenVelStep` steps (default 6 m/s; the constructor's `nVelGreen` sets the number of points instead), layered stacks in fixed 4 m/s steps. At this η a pole's FWHM in |G| is ≈ 0.5 % of v (≈ 13 m/s at 2.5 km/s), so at 6 m/s it spans about two samples: plotted peaks look angular and their plotted heights are sampled, not exact. Set `greenVelStep = 0.004` (or smaller) for smoother pictures; the `'green'` classification does not change, since it fits the pole rather than reading sampled heights.

The stored magnitude comes from this solve — `greenAbsG13Eff`/`greenAbsG13EffMap` (`|G13_eff|`, exactly what the existence score read); any other component is recomputed on request (`greenComponent`/`greenComponentMap`) at the same η and axis. A height threshold of your own tuned on a default-`η` `|G33|` does not apply to them: for a default-`η` spectrum call `GreenSolver.calc(structure, angle, velAxis, forceDir, displDir)` yourself (an `'Eta'` name-value pair sets another η). The `'refined'` method solves its own spectra at a larger η (MTH10). Why the classification uses η = 3e-3: [§5](reference/architecture/forwardCalcs_modeClassification.md#5-sawspectrum--the-shared-late-time-machinery).

### MTH9 🟠 The `'ritz'` method is a cheaper approximation of `'green'`, with its own limits

`'ritz'` scores a gated Ritz modal sum (`RitzScore`) the way `'green'` scores the Green's function, from the angle's one Ritz solve — no reference directions, no neighbouring angles. It is the only method of `calcMode='ritzOnly'` (chosen there automatically), and its candidates are stored on every `'full'` result too, where `classify('ritz')` gives the same entries as a `'ritzOnly'` calculation.

- **Half-space:** the gates (non-dispersion veto, dephased score, near-surface window weight, bottom-trapped gate) are calibrated together: do not change one alone, and keep the veto whenever you score Ritz modal sums yourself. `'low'` (the default) is always a subset of `'high'`. The rule needs `legOrder ≥ 2.3 × relDepth` and warns otherwise (`RitzScore.checkBasisHS`).
- **Layered:** `'low'` (default) keeps sagittal-type waves, `'high'` adds Love-type ones, and a trapped gate admits trapped waves only. **Leaky layered waves are never reported by `'ritz'`**; use `'green'` for them (PHY6).

The rules, known limits, and the agreement with measured maps, against `'green'` and `'refined'`: [§4](reference/architecture/forwardCalcs_modeClassification.md#4-ritz--ritzscore-half-space-rule), [§11](reference/architecture/forwardCalcs_modeClassification.md#11-validation) and [§12](reference/architecture/forwardCalcs_modeClassification.md#12-limitations).

### MTH10 🟠 The `'refined'` method models the probe spot — its own scale, its own spectra

`classify('refined')` (half-space `'full'` only) classifies with `GreenSpotScore`, which sums `G13_eff` coherently over the directions a finite probe spot excites and scores the late-time spectra of that sum. Its entries are ordinary `significantModes` entries with `method = 'refined'`; every candidate, accepted or not, stays in each angle's `refinedSaw.peaks`, and `refinedSawAgreement()` pairs the `'green'` and `'refined'` entries of each angle.

- Its standards (or a number on its `'high'` scale) are its own (MTH6). Its `'low'` is not guaranteed to be a subset of its `'high'`, so a point can be `'low'` without being `'high'`.
- It solves its own Green spectra (η = 7e-3, summed over the spot's 11 directions), not the displayed map (η = 3e-3, one direction): a refined peak need not sit on the displayed map's maximum. Its `type` compares the peak with the spot-weighted slowest shear velocity (`refinedSaw.info.minShearSpot`), not `minShear` at the angle itself.
- On a narrow branch that moves faster with angle than its own peak width per 0.5° (observed only for a synthetic material with Zener ratio `2*C44/(C11-C12)` ≈ 20), refined entries are not reliable.

How it is made and why: [§3](reference/architecture/forwardCalcs_modeClassification.md#3-refined--greenspotscore-half-space); the limit: [§12](reference/architecture/forwardCalcs_modeClassification.md#12-limitations).

### MTH11 🔴 To follow one branch across a layered sweep, track it — don't re-pick it at each point

Because each point is classified on its own (PHY2), a per-point pick can jump between coexisting branches. Where you need one continuous physical branch in a layered sweep, use `Calc.wavevecAngSweep2Layered` + `DispersionBranches.track` (`retentionLevel='full'`), which follows one branch by mass-weighted Ritz-eigenvector overlap (the modal assurance criterion, MAC) across sweep points instead of re-guessing from scratch at each. Seed the branch from the classification at the first sweep point: an entry's `modeNo`, e.g. the slowest `'SAW'` entry of `sweep.sweeps{1}.significantModesMap{1}`, which `DispersionBranches.trackSeededBranch(sweep)` does in one call — see `tutorial_06_film_on_substrate.m`.

CAVEAT: even that tracking isn't proof against every failure mode — it follows Ritz modes whatever their classification, including one that sits on a substrate limiting velocity (PHY2, last paragraph).

### MTH12 ⚪ `classify` / `reclassify` switch method and standard without recomputation — but `'saw_only'` retention cannot keep eigenvectors it never saw

`result.classify(method, threshold)` sets the active classification of any single-angle or sweep object (including `Calc.wavevecAngSweep2Layered`, every wavelength) from stored candidates; `reclassify(threshold)` is `classify(activeMethod, threshold)`. Switching costs no recomputation, apart from `'refined'` (half-space `'full'` only), which solves its spot directions on first use; a `'ritzOnly'` result supports `'ritz'` only. With `retentionLevel = 'saw_only'`, eigenvectors were pruned to the modes of the classification active when the calculation ran, so modes newly accepted by another method or a looser standard have no stored eigenvector (modal-shape plots of them need a fresh calculation). With `'none'` there are no eigenvectors anyway; `'full'` keeps everything. To start from another method, set `sawClassMethod` / `sawThreshold` before `calc()` (on `Calc.wavevecAngSweep2Layered` they are forwarded to every wavelength). Left unset (`''`, the default), `sawClassMethod` follows the `calcMode` of each run — `'green'` on `'full'`, `'ritz'` on `'ritzOnly'` — so the same object can be re-run in the other mode; the method in effect is `activeMethod`. A method you set explicitly is kept across runs, and on `'ritzOnly'` anything but `'ritz'` is an error. `classify` / `reclassify` change only the active classification, never `sawClassMethod` / `sawThreshold`. A `Calc.wavevecAngSweep2Layered` object has `classify` and `reclassify` but no `activeMethod`, `activeThreshold`, or `isStale()` of its own; read them from `sweeps{iw}`, the sweep of one wavelength. What is stored and what is carried over: [§8](reference/architecture/forwardCalcs_modeClassification.md#8-classify--reclassify-stored-fields-the-stamp).

### MTH13 ⚪ A stale result cannot be reclassified — recompute it

Every classification carries a stamp (`classificationVersion`), compared with `ModeClassification.CLASSIFICATION_VERSION`, which changes whenever a classification rule or a stored candidate field changes. A saved result with another stamp still loads (OTH3) and shows its maps and modes, but `isStale()` returns `true` (on a `Calc.wavevecAngSweep2Layered`, ask `sweeps{iw}.isStale()`) and `classify` / `reclassify` refuse it (`...:classify:stale`) — recompute it from its structure and options; don't compare its scores or entries with a current result's. Why: [§8](reference/architecture/forwardCalcs_modeClassification.md#the-classification-stamp), "The classification stamp".

## Standard numerics

### NUM1 🔴 MATLAB complex transpose: `'` conjugates, `.'` doesn't

`A'` is the *Hermitian (conjugate)* transpose — for a complex array, `imag(A') = -imag(A)`. `A.'` is a plain shape-only transpose.

For any complex physics quantity (Green's function `G`, eigenvectors `psi`/`U`, wavevector `k3`) always use `.'` unless conjugation is actually intended. The mistake is silent and costly to find: `imag(G_ref(:)')` in a validation script negates the imaginary part, making a correct result look like a 2× discrepancy.

### NUM2 🔴 `ElastConsts` without `symmNo` guesses the symmetry class — pass it explicitly

`ElastConsts(consts)` auto-detects the symmetry class when `symmNo` is omitted. From a vector of independent constants this is reliable; from a Voigt matrix or a tensor it counts the distinct numerical values, which miscounts whenever two independent constants happen to be numerically equal (a cubic crystal with `c12 = c44`, say, looks isotropic). The mistake is silent — the object builds, with the wrong constant set. Pass `symmNo` whenever the constants come as a Voigt matrix or tensor (`ElastConsts.SymmClassLibrary` lists the class numbers).

### NUM3 🔴 `pcolor(x, y, C)` treats a sample as its cell's lower-left corner, not its centre

`pcolor` (used with `shading('flat')` throughout `guliwerr/visual/@Vis/`) draws `C(:,j)` from `x(j)` to `x(j+1)`, not centred on `x(j)` — and separately drops the last row/column of `C` entirely under `'flat'` shading. For a 1°-step angular sweep, this means the cell you'd expect to read as angle 25° actually spans 25.0–26.0°, a full step off from what the axis label implies.

`Vis.pcolorCenteredGrid(x, y, C)` fixes this: it converts sample coordinates to cell-edge coordinates (interior edges at consecutive midpoints, endpoints linearly extrapolated) and pads `C` by one dummy row/column to match, so `pcolor(xEdges, yEdges, CPadded)` draws each sample centred where it belongs. Used by `Vis.plotGreenMap`, `Vis.plotDispersionMap`, `Vis.plotDispersionCurves`, `Vis.plotRitzSpectrumMap`, and `Vis.plotPolarHeatmap`. Route your own `pcolor` plots of discrete per-step data (angle, kParRange, etc.) through it too; dense field meshes such as the modal-shape plots need no centring.

### NUM4 ⚪ Eigenvector retention level — the largest per-angle array a sweep keeps

`Calc.angSweepHS`/`Calc.angSweepLayered`'s `retentionLevel` controls how much of the Ritz eigenvector matrix survives per angle after `calcRitz()` returns:

- `'full'` — keep everything (needed for modal-shape visualisation). `eigVecs_con` is `nCon × min(nMaxModes, nCon)` complex doubles, `nCon ≈ 3·Σ legOrder` (summed over the layers; the basis length minus three constraint rows per layer), so it grows as `(3·Σ legOrder)²` until `nMaxModes` caps the columns. Per angle that is ~0.35 MB for a half-space at `legOrder` 50, 0.9 MB for two layers at `legOrder` 40, ~2.7 MB for two layers at 70, and ~6.6 MB for four layers at 70 (capped by `nMaxModes` = 500). A 91-angle sweep at `'full'` holds ~30–600 MB of eigenvectors; a nested sweep (e.g. `Calc.wavevecAngSweep2Layered`) multiplies that across every outer-loop point (20 wavelengths × 91 angles: ~0.6–12 GB).
- `'saw_only'` — keep only the eigenvectors of the modes classified as significant (enough for modal-shape plots of the classified modes, but not for `Inv.association`, which needs every eigenvector and so `'full'`; see MTH12 for `classify` / `reclassify` afterwards).
- `'none'` — velocities only, minimum RAM; recompute a single angle's eigenvector on demand later if needed.

Use `'none'` for any sweep step that isn't the final, visualisation-bound one — this is why `Calc.wavevecAngSweep2Layered` and `Calc.wavevecAngSweep2Plate` default to `'none'`, not `'full'`. Pruning happens inside the sweep loop, per angle, once that angle's classification candidates and limiting entries are stored, so peak RAM is about one angle's eigenvector matrix rather than all angles' at once, and `classify` never needs the pruned eigenvectors. What each level keeps in detail: `docs/reference/architecture/forwardCalcs.md`, "`retentionLevel`"; the pruning mechanism: `docs/reference/architecture/ritz.md`, "Mode retention".

### NUM5 ⚪ `useParallel` — thread pool, default `"auto"`

The angle-sweep classes (`Calc.angSweepHS`, `Calc.angSweepLayered`, `Calc.isofreqAngSweepLayered`, and forwarded through the 2D sweep classes) accept `useParallel` (`"auto"` | `"always"` | `"never"`), resolved by `guliwerr/core/ParallelPolicy.m`. `"auto"` switches on for ≥20 angles; without the Parallel Computing Toolbox the run falls back to serial execution with a warning. It runs as a MATLAB **thread** pool, not a process pool — deliberately, since a thread pool shares memory rather than duplicating broadcast data (elastic tensors, `CalcParams`) per worker, which matters directly against NUM4. Use `"never"` when debugging (`parfor` mangles stack traces/breakpoints). Don't judge parallel speedup from a small/cheap fixture — the benefit only shows up once the per-angle numeric cost dominates MATLAB's own dispatch overhead (see `docs/reference/architecture/core.md`'s "ParallelPolicy" section for the actual benchmark).

### NUM6 ⚪ Your own `parfor` over single-angle objects: pass the Ritz precompute in

The sweep classes share one angle-independent Ritz precompute across all angles; building it costs more than a Ritz solve (~0.3–0.8 s at `CalcParams(30,70)`, against ~20 ms per solve). It is a `Transient` property, so it does not travel with an object into a `parfor` worker: a `Calc.singleHS`/`singleLayered`/`singlePlate` object arrives on a worker without it and rebuilds it for every angle — correct, but slower than a serial loop. If you write your own parallel loop, take the precompute once on the client (`p = cs{1}.ritzPrecompute()`, or a sweep's `sweep.ritzPrecompute()`) and pass it in: `calcRitz(p)` or `calcAll(mode, retentionLevel, p)` (plate: `calcAll(retentionLevel, p)`). A rebuild on a worker prints `...:ritzPrecomputeBuiltOnWorker`. The sweep classes already do this.

## Other specifics

### OTH1 🔴 Two coordinate conventions, for two different jobs — know which one you're looking at

GuliweRR uses two coordinate frames by design, each suited to what it's for, and they don't share an axis order:

1. **Input coordinates — `dirMat` (`Material` input)**, for *specifying* a crystal cut. Rows are direction-cosine vectors, ordered by physical importance: row 1 = surface normal (the plane you've cut), row 2 = in-plane propagation reference (lab angle = 0), optional row 3 completes the right-handed frame.
2. **Computational coordinates — the standard/literature convention** (internal, solver-facing), for *computing*: x1 = propagation, x2 = in-plane, x3 = surface normal/depth (z = out-of-plane) — the usual convention in the literature this project cites. `Christoffel`, `RitzSolver*`, and `GreenSolver` all natively compute and return components in this convention, not in `dirMat`'s row order.

`Material.GetStandardDirMat` is the one place the two conventions meet: it consumes a `dirMat` and produces the standard frame the solvers actually use internally. If you're indexing displacement/eigenvector components directly (rather than going through a plotting helper that already handles this), check which convention that particular array is in — assuming `dirMat`'s row order applies to a solver's output (or vice versa) silently mislabels which physical direction you're reading.

**`dirMat` means input coordinates everywhere in the public API, with no exceptions — if a function's name contains "DirMat", it takes/returns a genuine `dirMat`.** This includes `Material.Euler2DirMat`/`DirMat2Euler`: the underlying Bunge Euler-angle maths works in Bunge's own convention (row/column 3 = his distinguished tilt axis, not GuliweRR's surface normal), and these two apply the necessary `[3,1,2]`/`[2,3,1]` reorder internally, so no manual reordering is needed around them. The raw, unreordered Bunge-frame conversion is available under names that deliberately do *not* contain "DirMat" — `Material.Euler2DCM`/`DCM2Euler` — for generic proper-rotation maths with no physical `dirMat` meaning (e.g. checking a `dcm` given in a paper's own Bunge row order). See `docs/reference/architecture/core.md`'s `Material` section for the derivation.

### OTH2 ⚪ Path setup — never `addpath(genpath(...))`

Always use `run('guliwerr/init.m')`, which adds each subfolder explicitly. `genpath` recurses into `@ClassName/` subfolders and breaks MATLAB's class lookup.

### OTH3 ⚪ A saved result is tied to the class definitions it was saved with

`Calc.saveResult` stores the result object itself, not a neutral data format. Loading it rebuilds the object against the *current* `Calc*` class definitions: a property added since then gets its default, a removed one is dropped with a warning, and an incompatible change makes MATLAB fall back to a plain struct, which `Calc.loadResult` refuses (`CalcCommon:loadCalc:notAnObject`). The file's `meta.guliwerrVersion` names the GuliweRR release it was saved with — load it with that release, or recompute. `meta.gitCommit` is the repository's HEAD at save time (uncommitted edits are not recorded); it is empty outside a git checkout, such as a downloaded release, where the `commitMismatch` warning on load never fires either. In a checkout the warning is only a reminder; it does not mean anything is wrong. Saved results are working files, not an archive format. A result saved under an older classification rule is stale: MTH13.

### OTH4 ⚪ `Inv.association` needs every measured angle to coincide with a sweep angle

`Inv.association` matches each measured point to the forward sweep at the same angle. It errors (`ModeAssociation:fromSweep:angleNotFound`) if a measured angle has no sweep angle within `angleTolerance`, whose default is 1e-6°. Build the sweep's `angles` array from the measured angles (or round the measured angles onto the sweep grid), or widen `angleTolerance` deliberately.

### OTH5 ⚪ `Vis.plotAngularMap(…, 'polar')` needs polar axes; the Green-map plots need plain ones

`Vis.plotAngularMap` with `plotType='polar'` draws on a `polaraxes` and, given an axes handle `hAx`, expects one. `Vis.plotGreenMap` and `Vis.plotDispersionMap` with `plotType='polar'` draw on plain Cartesian axes, never a `polaraxes` (`Vis.plotPolarHeatmap` builds the polar look itself, as `pcolor` cannot be drawn on polar axes). Do not reuse one axes handle between the two kinds of plot.
