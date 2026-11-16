# ModeClassification — architecture and design rationale

`guliwerr/forwardCalcs/ModeClassification.m` and its scoring classes, all in
`guliwerr/forwardCalcs/`: `GreenScore.m` (half-space `'green'`), `GreenSpotScore.m` (half-space
`'refined'`), `RitzScore.m` (`'ritz'`, both geometries), `GreenScoreLayered.m` (layered `'green'`)
and `SawSpectrum.m` (the late-time machinery they share). Companion to `forwardCalcs.md` (read that
first for how `CalcSingleHS` / `CalcSingleLayered` and the sweep classes drive this). The class
consumes `RitzSolverSAW`'s per-mode diagnostics plus Christoffel and Green results and turns them
into a classified list of surface-wave modes — what a user wants ("where is the Rayleigh mode,
where are the Sezawa modes") rather than a bare list of Ritz eigenvalues.

Value class, default-constructed (optionally with name-value overrides of `threshold`,
`velMatchTol`, `matchScoreThresh`, `associateLimiting`, and `carryEntries`, e.g.
`ModeClassification('threshold','low')`; any other name is an error), then `calc()` is called and
the result reassigned — the same pattern as `Christoffel` / `GreenSolver` / `RitzSolverSAW`. Users
normally never call it directly: the Calc classes do, through `identifyModes()` and `classify()`.

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface waves in
> strongly anisotropic materials. Thresholds, scores, and the interface may change.

> **Terms** used below without further explanation.
>
> - **TGS, SRAS** — transient grating spectroscopy (a standing optical grating excites the surface
>   with an in-plane force; the out-of-plane displacement is detected) and spatially resolved
>   acoustic spectroscopy (a laser-ultrasonic technique that maps the surface-wave velocity).
> - **cut** — a crystal cut: the orientation of the sample surface, written by its Miller indices,
>   e.g. Ni(110). The propagation angle within the surface is what a sweep varies.
> - **limiting velocity** — a bulk-wave velocity (from Christoffel's equation) above which a surface
>   wave can radiate into the bulk. vS / vTf = slowest / fastest shear limiting velocity, vL =
>   longitudinal one; `minShear` = vS, the plain minimum over `shearLimVels` (section 7).
> - **G_Hill, vRef** — G_Hill is the Hill-average shear modulus (the mean of the Voigt and Reuss
>   bounds) of the surface material; vRef = √(G_Hill/ρ) is the reference velocity that sets the gate.
> - **candidate, pole** — a candidate is a local maximum of a spectrum that may be a surface wave;
>   a pole is a damped resonance, which is how each peak is modelled (position, half-width γ,
>   height).
> - **probe spot** — the finite area an experiment excites and detects, which spreads the
>   propagation direction over a fraction of a degree (σθ).
> - **dephasing** — the loss of height of a narrow peak when the probe spot sums directions whose
>   pole velocities differ.
> - **late-time gate** — the spectrum over velocity is Fourier-transformed to a time response, and
>   the first few periods are discarded: surface waves keep ringing, bulk transients die out.
> - **Van Hove edge** — a sharp edge or cusp in the bulk-wave continuum, where the density of bulk
>   states changes abruptly; it is short-lived in time and so removed by the gate.
> - **sagittal, SH** — a sagittal wave moves in the plane of the propagation direction and the
>   surface normal (Rayleigh, Sezawa); an SH (shear-horizontal, Love-type) wave moves perpendicular
>   to it.
> - **relDepth, legOrder** — the depth of a layer of the Ritz domain in wavelengths, and the
>   Legendre order of its depth basis (`docs/pitfalls.md`, MTH1).
> - **MAC** — modal assurance criterion, the normalised overlap of two eigenvectors.

## 1. Overview

### Existence, identity, label

Classification answers three questions, in this order:

1. **Existence** — is there a surface wave here that an experiment would see? Decided by one
   score per spectral peak, from the method's own spectrum. A peak is a surface wave if its score
   reaches the threshold; nothing else decides it.
2. **Identity** — which Ritz mode represents it? Answered from Ritz. The Ritz mode supplies the
   entry's `modeNo`, `modeVel`, and eigenvector. The entry's `vel` is the velocity of the method
   that decided existence: the Green (or spot) peak on `'green'` / `'refined'`, the mode's
   eigenvelocity on `'ritz'`.
3. **Label** — `'Rayleigh'` or `'pseudoRayleigh'` on a half-space, `'SAW'` with an
   `aboveSubstrateST` flag on a layered stack — applied last, from one velocity against the slowest
   shear limit. It has no effect on existence or identity.

### Three methods

| Method | Geometry | Needs | Existence from | Default standard |
|---|---|---|---|---|
| `'green'` (default) | half-space, layered | `calcMode = 'full'` | `GreenScore` / `GreenScoreLayered`: each `G13_eff` peak fitted as one damped pole, corrected for a late-time gate and aperture dephasing | half-space `'high'`, layered `'low'` |
| `'refined'` | half-space only | `calcMode = 'full'`; solved on first use | `GreenSpotScore`: late-time spectra of `G13_eff` summed coherently over the probe spot's directions | `'high'` |
| `'ritz'` | half-space, layered | nothing extra (the only method on `'ritzOnly'`) | `RitzScore`: a late-time Ritz modal sum that mimics `G13` | `'low'` |

The candidates of every method available on a result are scored once, threshold-free, and stored
side by side. `classify(method, threshold)` on a Calc result chooses which of them fills
`significantModes` and the `saw*` maps; switching method or standard only re-thresholds. On a
`'full'` result the Ritz candidates are computed too, so `'ritz'` costs nothing extra there and
equals the `'ritzOnly'` classification of the same Ritz basis.

Layered stacks have no `'refined'`: there is no spot score for stacks.

### Standards and grades

Each method takes a **standard**: `'high'` (sensitive — every mode a good, high-averaging TGS map
shows), `'low'` (distinct modes only — poor TGS data, SRAS), or a number on that method's own
scale. The scales differ between methods: compare a score only with its own method's thresholds,
and prefer the named standards.

| Method | `'high'` | `'low'` | A number is compared with |
|---|---|---|---|
| `'green'`, half-space (`GreenScore`) | `score` ≥ 0.16 | `lowScore` ≥ −0.79 (structure-relative) | `score` |
| `'refined'` (`GreenSpotScore`) | `scoreHigh` ≥ −0.125 | `scoreLow` ≥ −0.55 (structure-relative) | `scoreHigh` (the entry's `score`) |
| `'ritz'`, half-space (`RitzScore`) | `dScore` ≥ 0.16, with its gates | `dScore` ≥ 0.70, the same gates | `dScore`, with the gates |
| `'green'`, layered (`GreenScoreLayered`) | `score` ≥ 0.20 | `score` ≥ 0.70 | `score` |
| `'ritz'`, layered (`RitzScore`) | `dScore + polTerm` ≥ −0.7 | ≥ 0.65 | `dScore + polTerm`, with the trapped gate |

On the half-space `'green'` and `'refined'` methods `'low'` lives on a second, structure-relative
scale, so no number reproduces it: a number always means the `'high'` scale there. A number tuned
on one method means nothing on another, and `classify(method)` without a standard falls back to
the new method's default rather than carrying the old standard over. Scores are log10 amplitude
ratios; an entry's `scoreDb` is the same value in dB (20 × `score`). The constants behind these
values are listed in section 9.

What an experiment observes also depends on the technique: a mode with low out-of-plane surface
displacement is barely visible to an out-of-plane-sensitive technique (e.g. SRAS), while one
sensitive to in-plane motion may still recover it.

Every surface-wave entry carries a **grade**: `'strong'` if its candidate passes the method's
`'low'` standard, `'weak'` if it passes only the active (looser) one. At `'low'` every entry is
`'strong'`.

### Dispatch inside `ModeClassification.calc`

`calc(ritzSolver, structure, christoffelResult, greenResult, ritzScore)` picks its path from which
arguments are present:

| `mode` / `method` | Trigger | Static method | Existence |
|---|---|---|---|
| `'ritzOnly'` / `'ritz'` | `christoffelResult` or `greenResult` omitted or `[]` | `recognizeRitzOnlyScore` | half-space `RitzScore.acceptedHS`; layered `RitzScore.acceptedLayered` |
| `'full'` / `'refined'` | `greenResult.refined` present | `recognizeRefined` | `GreenSpotScore.accepted` |
| `'layered'` / `'green'` | both present, `structure.nLay > 1` | `recognizeLayered` | `GreenScoreLayered.accepted` (`greenResult.greenScoreLayered`) |
| `'full'` / `'green'` | both present, `structure.nLay == 1` | `recognizeFull` | `GreenScore.accepted` (`greenResult.greenScore`) |

The inputs are the stored, already scored candidates (`greenResult.greenScore.peaks`,
`greenScoreLayered.peaks`, `refined.peaks`, `ritzScore.cands`); `ModeClassification` only
thresholds, matches, and labels them. A Green input without its score, or a Ritz-only input without
`ritzScore` or its stored candidates, is an error (`ModeClassification:recognizeFull:noScore`,
`:recognizeLayered:noScore`, `:ritzOnly:noScore`, `:ritzOnly:noCandidates`). `christoffelResult`
is read on the Green paths only (the Rayleigh/pseudoRayleigh label and the limiting entries); the
Ritz-only path reads nothing but `ritzScore.cands`, whose `belowShear` already carries the label.

## 2. `'green'` — `GreenScore` (half-space)

### How the score is made

For one propagation direction (`GreenScore.candidates`, `probeSlope`, `score`, `lowScore`;
`direction` runs them all), on the band [0.80 vS, 1.03 vTf] of the complex `G13_eff` at
η = `SawSpectrum.ETA` = 3e-3 (vS / vTf = slowest / fastest shear limiting velocity,
vRef = √(G_Hill/ρ)):

1. **Candidates.** Local maxima of |G13_eff| over the whole band with half-prominence width
   ≤ 12 % of v, kept at or below vTf(1 + 0.012). The normaliser fBand is the median |G13_eff| over
   [0.80 vS, vTf] of the same direction.
2. **Pole fit** (`GreenScore.fitPeaks`). The complex model (a0 + a1 x)/(x − q) through the peak
   sample and the samples ±`FIT_SPAN` (16 m/s) away gives the pole position, the half-width
   γ = −Im q (floored at the damping width η v/2) and the pole height. A failed fit (no damped pole
   within the stencil) keeps the grid peak height. rawScore = log10((pole height − base)/fBand).
3. **Analytic late-time gate.** A pole of half-width γ (in velocity) rings in time as
   exp(−2πγt/λ). Counting its spectral height only from t0 = 3 periods of vRef multiplies it by
   exp(−6πγ/vRef) (λ cancels): gateTerm = log10 of that. A true SAW (γ set by η only) loses
   ~0.01 dex; a strongly leaky pole (γ/vRef ≈ 0.05) loses ~0.4 dex.
4. **Aperture dephasing** (`GreenScore.dephaseFactor`). A finite probe sums a Gaussian spread of
   directions (σθ = `APERTURE_SIGMA` = 0.4°). On a branch of slope dv/dθ that spreads the pole
   velocity by δv = |dv/dθ| σθ, a Gaussian envelope on the ringing. The gated height ratio to the
   same pole without spread is F = ρ √(π/2) erfcx((ρ + q)/√2) exp(−q²/2), with ρ = γ/δv and
   q = 6π δv/vRef: dephaseTerm = log10 F. Narrow poles on steep branches lose most, broad leaky
   poles almost nothing.
5. **Slope** (`GreenScore.probeSlope`, `combineSlopes`). Two probe solves at θ ± 0.5° on a ±48 m/s
   stencil (16 m/s spacing) around each pole; each reads the *local* maximum of |G| nearest the
   pole, not the stencil's arg-max. Both sides with the same sign and within a factor 2 → their
   mean, otherwise the gentler one; one side only → that one. A pole neither probe finds is
   charged the stencil's bound (0.096 km/s per degree).
6. **Scores.** score = rawScore + gateTerm + dephaseTerm (the `'high'` scale; the entry's `score`).
   logStrength = score + log10 fBand is the gated, dephased prominence in |G13_eff| units;
   lowScore = logStrength − `logStrengthTop` of the structure (`GreenScore.reference`: the median,
   over 12 directions 0:15:165°, of each direction's largest logStrength — the structure's typical
   strongest wave; memoised per structure).

Standards: `'high'` = score ≥ 0.16; `'low'` = lowScore ≥ −0.79 (at least 10^−0.79 ≈ 0.16 of the
structure's typical strongest wave); a number = score ≥ that number. Empirically every `'low'`
point is also a `'high'` point; this is not guaranteed by construction, because the two scales use
different normalisers.

### Why each choice

- **Tuned against measured maps.** The thresholds, σθ and the method itself were tuned against
  measured TGS maps (section 11).
- **Analytic gate instead of an FFT.** For one pole the late-time gate (section 5) is exact in
  closed form, so the score needs neither a fine uniform axis nor padding or tapers; the band
  alone suffices.
- **Dephasing.** It removes the steep, narrow branch ends that the measured maps do not show. A TGS spot averages directions over a fraction of a
  degree, so a branch that moves fast with angle does not ring in the experiment, however sharp it
  is in one direction.
- **Width floor, failed-fit height, two-sided local-maximum probe.** Each removes an estimation
  artefact of a weak pole on a strong neighbour's skirt: a 3-point fit can return a width below
  the damping width, a failed fit's model maximum can exceed the grid peak, and an arg-max probe
  jumps to the neighbour next to a crossing or cusp, reads the stencil bound, and dephases a real,
  strong wave away. A probe that jumps only ever reads a steeper slope, so the gentler side is the
  one to trust. Shared with `GreenScoreLayered`.
- **σθ = 0.4°.** The spot's direction spread is not measured. A larger spread (0.8°) removes
  about half of a bound Rayleigh branch that hugs vS on a steep branch (tens of m/s per degree)
  and is weak in `G13_eff`; at 0.4° it stays, at no measurable cost elsewhere.
- **Probe slope, not neighbouring angles.** Slopes from a sweep's neighbouring angles fail at
  coarse steps and make a single angle differ from a sweep angle. The probe makes every angle
  self-contained: a single angle classifies bit for bit as inside any sweep, at any angular step
  or range.
- **Structure-level `'low'` reference.** A median over the computed sweep would make the `'low'`
  labels depend on which angles were computed; 12 fixed directions do not.
- **Fixed pole-fit span.** Fitting through adjacent samples spans only two steps and misjudges
  broad leaky poles on a fine axis. A fixed ±16 m/s stencil makes the result independent of the
  axis step from 4 to 16 m/s, so the Green axis step (`greenVelStep`, default
  `GreenScore.DEFAULT_STEP` = 6 m/s) is free to serve the pictures.
- **Pole height as numerator, band median as normaliser.** The fitted pole height does not depend
  on how the grid samples the peak; the band median is a smooth property of the direction.
- **No threshold relative to `max(abs(G))`.** The tallest peak at an angle is usually a bound
  Rayleigh pole whose sampled height depends on the grid, so one peak's existence would depend on
  another peak's height and on sampling, and a smoothly growing peak could flip at a single angle.
  Each existence score is normalised within its own direction instead (the band median for the
  Green scores, the late-time floor for `RitzScore`); only the half-space `'low'` standards of
  `'green'` and `'refined'` compare with one reference per structure.

## 3. `'refined'` — `GreenSpotScore` (half-space)

`classify('refined')` on a `'full'` half-space result (`CalcSingleHS`, `CalcAngSweepHS`) solves the
spot directions on first use (printing the cost), stores the scored candidates in `refinedSaw`,
and from then on switches to and from `'refined'` by re-thresholding only. A new calculation clears
them.

### How it is made

For one direction θ (`GreenSpotScore.average`, `direction`, `identify`; vS / vTf / vL = slowest /
fastest shear and longitudinal limiting velocity at θ, vRef = √(G_Hill/ρ)):

1. **Directions.** The probe spot excites a Gaussian spread of directions (`FWHM` 2.5°, cut at
   2.5σ): θ + o, o = −2.5:0.5:2.5° (11 directions). Each is solved once, alone: one `GreenSolver`
   call at η = 7e-3 for `G13_eff` on the axis 0.3 km/s + n × 8 m/s up to 1.05 × *its own* vL
   + 0.1 km/s. Q(θ + 180°) = Q(θ), so directions fold into [0, 180°); a sweep solves each folded
   direction once (`solveDirections`, a direction cache shared by all angles).
2. **Coherent average.** X = Σ w_j G13_eff(θ + o_j), complex, on the shortest of the 11 axes.
3. **Late-time spectra.** `SawSpectrum.lateTime` of X with the gate at 4 periods of vRef
   (detection) and at 2 (position): ×4 zero-padding (2 m/s output), tail taper 0.2. Candidates:
   maxima below vL with half-prominence width ≤ 12 % of v, kept at a reported velocity
   ≤ vTf(1 + 0.012), the same cap as `GreenScore`. Between vTf and vL the spectrum rings at the
   longitudinal branch point (a surface-skimming bulk wave just below vL), which the cap excludes.
4. **Scores.** scoreHigh = log10(prom4 / fU1), fU1 = median |X| (un-gated) over [0.80 vS, vTf].
   scoreLow = log10(prom4 / fT12), fT12 = the median, over the 12 directions 0:15:165°, of each
   direction's largest gate-4 prominence (`GreenSpotScore.reference`, memoised per structure).
   `'high'`: scoreHigh ≥ −0.125; `'low'`: scoreLow ≥ −0.55; a number: scoreHigh ≥ it. The entry's
   `score` is scoreHigh.
5. **Velocity.** The gate-2 peak whose half-prominence half-width contains the gate-4 peak
   (nearest if several), else the gate-4 position (`velDet`, also stored).
6. **Identity and label** (`GreenSpotScore.identify`). The spot peak is mapped to the nearest
   single-direction pole of the angle (a `GreenScore` candidate, accepted or not) within
   max(1.2 % v, 30 m/s) (`POLE_WIN_REL`, `POLE_WIN_ABS`), and that pole is bound by
   `matchRitzToGreenPeak`; with no pole in the window, `matchSpotPeak` at the reported velocity
   (section 6); `'Rayleigh'` below minShearSpot × (1 + 7e-3), else `'pseudoRayleigh'`.
   minShearSpot = Σ w_j minShear(θ + o_j), the same weighted average as the spectrum: a bound
   Rayleigh wave is below minShear in every direction of the spot, so its averaged peak is below
   the averaged minShear. 7e-3 is the method's own η, the resolution of the peak position.

The entries are ordinary `significantModes` entries (`method = 'refined'`, `greenPeakVel` = the
spot peak velocity) and carry the same limiting entries as `'green'`.
`refinedSawAgreement()` pairs the `'green'` and `'refined'` surface-wave entries of each angle
one to one (within max(1.2 % v, 30 m/s)) for comparison.

### Why

- **What it buys.** It uses the spot-averaged velocity and a modified accept decision: the
  averaged spectrum rejects narrow poles on steep branches that the spot dephases, and peaks whose
  single-direction and averaged positions disagree.
- **Coherent sum, then gates.** The spot drives every direction at once, so the response is the
  complex sum; dephasing then shows up in the late-time spectrum without a model (`GreenScore`'s
  analytic dephasing term approximates the same effect from a slope probe).
- **Two gates.** The late gate (4 periods) separates surface poles from short-lived bulk features
  most reliably; the early gate (2) keeps more of the pole's own lineshape and places the peak
  closer to where a measured map puts it (computed with GuliweRR on the published cuts).
- **η = 7e-3, 8 m/s.** Calibrated with the rest. The axis need not resolve the pole
  finely: the late-time spectrum is interpolated by zero-padding.
- **Every direction on its own axis.** A direction's spectrum then does not depend on which other
  directions were solved, so a single angle is bit for bit the same angle in any sweep, at any
  step.
- **A refiner, not the default.** It needs 0.5° directions around every angle: twice the
  Green solves of a 1° `'green'` map, which adds about half the time of the `'full'` sweep itself.
- **`'low'` ⊄ `'high'`.** scoreHigh divides by the direction's own band median, scoreLow by one
  structure-level number, so a point can pass `'low'` and not `'high'`. Not forced by
  intersection: each standard is used as calibrated.

## 4. `'ritz'` — `RitzScore` (half-space rule)

`'ritz'` needs no Green's function. Its candidates come from the angle's one Ritz solve
(`CalcSingleHS.calcRitzScore`, part of `calcAll` on both calculation modes).

### The modal sum (the Ritz spectrum)

With the Ritz phase velocities `v_n` and the free-surface displacements `u_n` (standard frame;
top row of `calcDisplacement`):

```
S(v) = Σ_n u3_n · conj(u1_n) / (v_n² − v²(1 + iη)),   η = RitzScore.ETA = 1e-2
```

summed over the modes with `v_n ≤ 1.2 ×` the longitudinal limiting velocity (`MODE_VEL_FACTOR`;
layered: the fastest of any layer). For M-orthonormal eigenvectors this is the finite Ritz
domain's driven surface response to an in-plane force, read as out-of-plane displacement: the Ritz
analogue of `G13`. Ritz has no radiation boundary condition, so a leaky resonance appears as a
cluster of closely spaced, comparably weighted modes rather than one clean eigenmode; no per-mode
quantity read alone shows it, the coherent sum does. The discrete bulk modes of the finite domain
dephase in time like the bulk continuum of a true half-space, so the same late-time gate
suppresses them, while a surface-wave mode keeps ringing. The `G13`-like numerator `u3·conj(u1)` matters more than any other choice
(an `|u3|²` weighting is clearly worse). η is larger than `SawSpectrum.ETA` because each Ritz mode
is an exact pole: 1e-2 keeps it resolved on the 2 m/s axis and decayed within the FFT window.

The numerators (`RitzScore.modalNumerators`) are taken before `retentionLevel` pruning and stored
(`ritzScoreNum`). `RitzScore.spectrum(vn, num, velAxis)` rebuilds S at any velocity resolution
with no solve; it serves the existence score below and the "Ritz spectrum" view (`Vis.plotRitzSpectrumCurve` / `plotRitzSpectrumMap`).
`RitzScore.decompose(vn, num, v)` gives each mode's share |term_n| / Σ_m |term_m| of S at one
velocity — the contributing modes of a peak.

### The score

S is scored with `SawSpectrum`'s code (section 5): on the uniform axis
`SawSpectrum.velocityAxis(1.05 vL, 2 m/s)`, without padding or tail taper (an analytic spectrum
costs nothing per point), the same gate (3 periods of vRef), candidate search and width cap, and
`score = log10(prominence / floor)`. The floor is this direction's own (`RitzScore.ownReference`:
the median late-time |S| over the band [0.8 vS, vTf]). Above vTf the slab's modes are not converged
in `legOrder` at the usual bases, which would make every score depend on the basis. The own floor
needs no Ritz solve for any other direction, so a single angle classifies exactly as inside a
sweep.

### The half-space rule

Each candidate is scored once, threshold-free, by `RitzScore.direction` and stored
(`ritzScoreCands`), so a reclassification needs no eigenvectors.

1. **Local supporting mode** (`RitzScore.localMode`): the strongest modal-sum term among the modes
   within the peak's half-width (at least 4 m/s), else the nearest mode. The *global* argmax over
   all modes can name a strong Rayleigh pole far below a pseudo-SAW peak, and then the
   veto, slope and velocity would all belong to the wrong mode. The entry's `vel` is this mode's
   eigenvelocity.
2. **Dephased score** `dScore = score + log10 F`, `F = GreenScore.dephaseFactor(γ, δv, vRef)` with
   `γ = ETA·v/2` (the modal sum's own pole half-width) and `δv = |dv/dθ| × APERTURE_SIGMA`
   (0.4°). The slope comes from the same solve, without a second one (`RitzSolver.angleSlope`,
   defined in `ritz.md`: exact Hellmann–Feynman derivative).
3. **Non-dispersion veto**: the supporting mode's `freqDisp ≥ MIN_FREQDISP = 0.3`
   (`freqDisp = exp(−10·((v_p − v_g)/v_p)²)`, defined and motivated in `ritz.md`), so
   `|v_p − v_g|/v_p` up to about 35 % passes. The veto removes slab modes that ring coherently
   near the longitudinal limit (e.g. Ni(110) near 5 km/s), which Green's function does not show as
   surface waves. A higher cut-off than 0.3 drops real branches.
4. **Window mass**, for candidates above `minShear` only (`RitzScore.windowMass`): log10 of
   `Σ r·E1·freqDisp` over the modes within `max(1.2 % v, 30 m/s)` of the peak on its side of
   `minShear` (`r` = surface `|u3|²` over the depth-mean `|u|²`, divided by relDepth; `E1` = energy fraction in
   the top wavelength; modes above 0.98 v_L weigh nothing). A slab artefact rings from one or two
   bulk-like modes; a real leaky wave is carried by several surface-concentrated ones.
   Depth-normalised by + log10(relDepth/30) (`HS_WINDOW_DEPTH_REF`): the summed weight of a real
   leaky wave falls with the slab depth, as more modes share it. Must be
   ≥ `HS_MIN_WINDOW_MASS` = −1.55.
5. **Bottom-trapped gate** (same candidates): rejected if the supporting mode has > 15 % of its
   `|u|²` in the deepest 10 % of the slab and `bottom10/bottom25 > 0.8` (uniform: 0.4). A strongly
   leaky pseudo-SAW that the rigid-bottom slab cannot radiate splits into "U-shaped" modes (energy
   at surface and bottom); the lower member wins the peak tens of m/s below the real branch.
   Values for relDepth ≈ 30.
6. **Basis guard** (`RitzScore.checkBasisHS`): warns once if `legOrder < 2.3 relDepth`; below it
   the slab does not resolve the band up to vL and hundreds of slab modes near vL pass. This is
   stricter than the general rule `legOrder ≥ 2 relDepth` of `docs/pitfalls.md` (MTH1): the
   half-space `'ritz'` classification needs the factor 2.3.
7. **Side of `minShear` and label**: a candidate is on the Rayleigh side if its supporting mode is
   below `minShear + RitzScore.slabResolution` (section 7). The side decides both the label and
   whether the window-mass and bottom gates apply.

| Standard | Gates (`RitzScore.acceptedHS`) |
|---|---|
| `'high'` | veto, `dScore ≥ HS_THRESHOLD_HIGH = 0.16`, Rayleigh side or (window mass ≥ −1.55 and not bottom-trapped) |
| `'low'` (default) | the same, with `dScore ≥ HS_THRESHOLD_LOW = 0.70` — always a subset of `'high'` |
| numeric | a dScore threshold with the `'high'` gates |

The gates and thresholds are calibrated together; change none alone (in particular, do not use a
threshold below 0.9 without the veto).

The Christoffel step (`calcLim`, a few ms) runs in `'ritzOnly'` mode too (search bound, label,
slab resolution). The `'ritz'` classification has no limiting entries.

`classification.ritzScorePeaks` lists every scored candidate at the angle, accepted or not (field
set `RitzScore.emptyPeaks`): `vel` (the modal-sum peak), `modeNo` / `modeVel` (the supporting
mode), `score`, `dScore`, `prominence`, `relWidth`, `relHeight`, `share`, `freqDisp`, `slope`,
`windowMass`, `belowShear`, `vetoed`, `bottom10` / `bottom25`, `SH` / `polTerm`, `accepted`
(`freqDisp` / `windowMass` NaN on the layered rule, `SH` / `polTerm` NaN on the half-space).

## 5. `SawSpectrum` — the shared late-time machinery

`SawSpectrum` holds what makes the three scores comparable: the damping of the default existence
solves (`ETA` = 3e-3), the uniform velocity axis (`velocityAxis`: from `VEL_MIN` = 0.3 km/s in
`VEL_STEP` = 4 m/s steps), the late-time gate (`lateTime`, which also ramps the spectrum with a
0.1 km/s raised cosine at both axis ends) and the peak search (`candidates`). `GreenScore` reads
`ETA` and applies the gate analytically; `GreenSpotScore` and `RitzScore` call `lateTime` /
`candidates`. `latticeAxis` builds every Calc and score axis as 0.3 km/s + n × step, so a point
has the same value on a short axis and on a longer one; the layered Green axis uses the same
lattice (`GreenScoreLayered.gridAxis`).

**`lateTime(G, velAxis, vRef, pad, tailFrac, gatePeriods)`.** At fixed `|k|` the velocity axis is
a frequency axis, so an FFT over it is the impulse response `u(k, t)` of a grating excited at
`t = 0`. The first `GATE_PERIODS` = 3 periods of `vRef = √(G_Hill/ρ)` of the surface layer are
removed (a raised-cosine turn-on of `RAMP_PERIODS` = 0.6 periods, effectively abrupt). The end of
the causal window falls off over its last `TAIL_TAPER` = 20 % (raised cosine), and the gated time
signal is zero-padded `PAD_FACTOR` = 4× before it is transformed back, which interpolates the
late-time spectrum exactly onto a 4× finer axis. `GreenSpotScore` passes its own gates (2 and 4);
`RitzScore` passes `pad = 1`, `tailFrac = 0`.

**`candidates(G, velAxis, vSearchMax, ref, pad, tailFrac)`.** Local maxima of the late-time |G|
below `vSearchMax` with half-prominence width ≤ `MAX_REL_WIDTH` = 12 % of their velocity,
`score = log10(prominence / ref.floor)`. The search runs a little past the bound, so that a peak
at the bound is still interior.

### Why each choice

- **`G13`, not `G33`.** TGS excites with a thermoelastic in-plane force and detects out-of-plane
  displacement (Maznev et al. 1999, p. 2820), so `G13` is the measured component. `G33` over-weights the vertically polarised
  Rayleigh branch relative to pseudo-SAWs and Sezawa modes.
- **The ±k mean `G13_eff`, not the one-sided `G13_+`.** In TGS the standing grating drives the +k
  and −k waves together (Maznev, Akthakul & Nelson 1999, p. 2820). By reciprocity `G(−k) = G(k)ᵀ`,
  so +k and −k share every pole; for a bound SAW `|G13(+k)| = |G13(−k)|` exactly, but at a
  pseudo-SAW the two residues differ (computed with GuliweRR: by a factor 1.1–1.2 on Si(111) and
  2–8 on strongly anisotropic cuts). Maznev et al.'s observation that the two agree near resonances
  therefore holds for bound SAWs and mildly anisotropic cuts only. Scoring `G13_eff` puts the
  classified points on the spectrum that is displayed. It comes from the same single solve by
  reciprocity (`green.md`, "±k symmetrisation", which also says where the two coincide).
- **η = 3e-3, not the solver default 3e-4.** At 3e-4 a bound Rayleigh pole is ~0.5 m/s wide,
  narrower than the axis step, so its sampled height is grid luck and dwarfs a leaky pseudo-SAW of
  the same residue. A measurement-like damping broadens every pole to a resolvable linewidth
  (|G| FWHM ≈ √3·η·v, ≈ 0.5 % of v; derivation in `green.md`, "Unified limiting-absorption
  damping") while conserving its residue, so bound and leaky peaks of equal residue get comparable
  heights. It also makes each pole decay well inside the FFT
  window, which the gate needs.
- **A late-time gate.** Surface-wave poles ring for tens of periods and survive; bulk-continuum
  features — limiting-velocity shoulders, Van Hove edges, cusps, and broad domes — decay within a
  few periods and are suppressed. This is the calculated counterpart of cutting the first
  nanoseconds off a measured TGS signal, which removes the same transient bulk features (UTGS paper,
  arXiv:2510.10696, Fig. 3).
- **Gate counted in periods of `vRef`, so the method is dimensionless.** At fixed `|k|` the
  problem is scale-free; a gate in periods of the material's own shear velocity needs no grating
  wavelength, and the classification is invariant under a uniform scaling of all stiffnesses (a
  fixed gate in nanoseconds is not).
- **Prominence, not height.** Most false candidates are bumps on the flank of a stronger feature.
  A weak mode on such a flank has low prominence, and it is also hard to resolve in an experiment.
- **Zero-padding and a tail taper.** The damping sets how long a pole rings: a pole's half-width
  is `v·η/2`, so its ringing decays with `τ = 1/(π v η)`. The causal window of an FFT over the axis
  lasts `1/(2 dv)`. Zero-padding the gated time signal is exact band-limited interpolation, so
  candidates and prominences are read on a finer axis at no solve cost. On a coarse axis the window
  ends before the ringing has fully died away, and a hard cut there ripples the spectrum into many
  spurious candidates; the raised-cosine tail taper removes that.
- **A relative-width cap of 12 %.** The gate leaves a few broad shoulders on strongly anisotropic
  materials (relative widths of 0.15–0.18); genuine peaks are far narrower. A genuine peak's width (computed with GuliweRR) is a
  few tenths of a percent to a couple of percent of its own velocity, a shoulder's several percent
  to tens of percent — the figure to compare if you peak-search a Green spectrum yourself
  (`GreenPeaks.find` reports each peak's `width` and `prominence`).
- **No angular-persistence criterion.** Classification works from one angle's own spectrum (plus,
  for `'low'` on the Green methods, one structure-level reference), so a single angle is classified
  exactly as inside a sweep.

## Layered stacks

Two methods: `'green'` (`GreenScoreLayered`) and `'ritz'` (the layered Ritz rule). Both default
to `'low'`. Neither singles out a layer as "the substrate": the half-space scores take their band
and gate from the surface material, and on NiTi/MgO the film's shear limit (1.6–2.3 km/s) lies
*below* the Rayleigh wave (2.8 km/s), so a half-space band would exclude it. Which layer acts as
"the substrate" also changes with λ in a multilayer. An all-layer band, an own-period gate, and a
direction-own normaliser need no layer property.

### `'green'` — `GreenScoreLayered`

For one direction (`GreenScoreLayered.candidates`, `normaliser`, `probeSlope`, `polarisation`,
`score`; `direction` runs them all), on the complex `G13_eff` at `SawSpectrum.ETA` and the band
[0.3 vSmin, 1.05 vLmax] — vSmin / vLmax = the slowest shear / fastest longitudinal limiting
velocity over **all** layers (`bandLimits`, stored as `bandVSmin` / `bandVLmax`; a range only,
never a score input):

1. **Candidates.** Local maxima of |G13_eff| inside the band with half-prominence width ≤ 12 % of
   v; `GreenScore.fitPeaks` (±16 m/s, γ floored at η v/2, a failed fit keeps the grid height).
2. **Slope.** `GreenScore.probeSlope`, as on the half-space.
3. **Polarisation.** One full-tensor solve at pole and pole ± 16 m/s gives the residues r11, r22,
   r33 of the diagonal `G_ii` (3-point model a1 + r/(x − q)); SH = |r22| / (|r11| + |r22| + |r33|):
   ≈ 0 for a sagittal (Rayleigh/Sezawa-type) wave, ≈ 1 for a Love wave.
4. **Score.**
   score = log10(pole height − base) − 2πγτ/ln 10 + log10 F(γ, |dv/dθ| σθ) − N + K_POL (1 − 2 SH),
   with the gate τ = `GATE_OWN_PERIODS` = 10 periods of the wave itself (τ = 10/v), `GreenScore`'s
   dephasing factor F (σθ = 0.4°), N = log10 of the median |G13_eff| over the band of **this**
   direction, `K_POL` = 0.75 (SH undetermined → −K_POL). No reference
   over other directions.

Standards (one scale): `'low'` (default) = score ≥ 0.70 — sagittal-type waves (Rayleigh, Rayleigh
cap, both Sezawa families); `'high'` = score ≥ 0.20 — also Love-type waves; a number = score ≥
that number. `'low'` is always a subset of `'high'`.

Why:

- **Own-period gate.** Ten periods of the wave itself charge a leaky pole by its width alone, with
  no layer's shear velocity involved.
- **Polarisation term.** Where the Rayleigh character hands over to the steep Love arms, the
  sagittal-type points that `'low'` must keep are as weak as the strongest Love points, in the
  calculation and in the measured NiTi/MgO map (Grabec et al. 2024). Strength alone cannot separate them; SH can. With K_POL = 0.75, lowering the
  threshold admits Love waves before spurious peaks.
- **Direction-own median normaliser, wide band.** A 12-direction structure reference would cost
  12 extra directions per structure and, in a dispersion sweep, per wavelength; the direction's
  own band median gives the same classification on the NiTi/MgO stack (Grabec et al. 2024).
- **Axis.** The 4 m/s lattice (`forwardCalcs.md`, "The Green velocity axis") also holds a sweep's
  common axis, so a single angle scores bit for bit as inside any sweep.

`recognizeLayered` matches each accepted peak to a Ritz mode with the layered identity score
(section 6), labels every entry `'SAW'` with `aboveSubstrateST` =
`greenPeakVel` above the substrate's slowest shear limit, sorts the entries by velocity (for
readability only) and appends the substrate's limiting entries. `classification.greenScorePeaks`
holds every candidate with its parts (`GreenScoreLayered.emptyPeaks` fields plus `accepted`).

### `'ritz'` — the layered Ritz rule

`RitzScore.directionLayered`, called by `CalcSingleLayered.calcRitzScore`; built like the
half-space rule and consistent with `GreenScoreLayered`.

1. **Band.** Axis (`RitzScore.layeredAxis`, 2 m/s from min(0.3 km/s, 0.3 vSmin)), candidate search
   and floor all span [0.3 vSmin, 1.05 vLmax]; modes are kept up to 1.2 vLmax.
2. **Local supporting mode and dephased score**, exactly as on the half-space (`localMode`,
   `RitzSolver.angleSlope`, `dScore = score + log10 F`). The late-time gate and the dephasing use
   the surface layer's `vRef`: a gate at the slowest shear limit of any layer ranks Love waves far
   worse against spurious peaks.
3. **Polarisation term** `polTerm = LAY_K_POL (1 − 2 SH)`, `SH = |u2|²/|u|²` of the supporting
   mode at the free surface — `GreenScoreLayered`'s term, from the Ritz displacement instead of the
   Green residues. The value thresholded (and the entry's `score`) is `dScore + polTerm`.
4. **Trapped gate.** The supporting mode's energy fraction in the deepest 25 % of the whole domain
   (`bottom25`, density- and thickness-weighted over the layers, `RitzScore.modeProfiles`) must be
   ≤ `LAY_EBOT_MAX` = 1e-2. A guided wave decays into the substrate (1e-22 … 1e-12); a mode of the
   slab continuum reaches the rigid bottom (≥ 0.03). No `freqDisp` veto: guided waves in a stack
   are genuinely dispersive (the film sets a length scale), so `freqDisp` cannot separate them from
   finite-domain artefacts.
5. **Label:** `'SAW'`, with `aboveSubstrateST` = the supporting mode above the substrate's slowest
   shear limit (a label only); `vel` is the supporting mode's eigenvelocity, `matchScore` its share
   of the modal sum.

| Standard | Gates (`RitzScore.acceptedLayered`) |
|---|---|
| `'low'` (default) | `bottom25 ≤ 1e-2` and `dScore + polTerm ≥ LAY_THRESHOLD_LOW = 0.65` — sagittal-type waves |
| `'high'` | the same gate, `≥ LAY_THRESHOLD_HIGH = −0.7` — also Love-type waves |
| numeric | a threshold on `dScore + polTerm` with the same gate |

No threshold keeps all the sagittal-type waves without admitting Love waves: the steep lower Love arm
scores above the Rayleigh-cap ends, because the modal sum's pole half-width `ETA v/2` dephases a
steep branch far less than Green's does. `'low'` gives up those cap ends.

**Leaky layered waves are not reported.** A wave above the substrate's shear limit that radiates
into it is carried by slab modes that reach the bottom, so the trapped gate removes it with the
continuum; no other Ritz measure separates them. Use `'green'` for leaky layered waves. The bottom
energy decides trapping better than a velocity comparison: a sagittal wave above the substrate's
slowest shear limit in a symmetry direction, where it cannot couple to that shear wave, has
`bottom25` ≈ 1e-21 and is kept.

**Not unified with the half-space bottom gate.** On a half-space the leaky candidates are kept and
judged by the window mass, and the bottom gate rejects only U-shaped slab doublets. The layered
gate drops leaky candidates outright.

## 6. Identity

### Green decides existence, Ritz only assists identity

`recognizeFull` and `recognizeLayered` separate *does the physics show a real surface response
here?* — answered by the Green score alone — from *which Ritz mode best represents it?* — answered
by `matchRitzToGreenPeak`. A weak Ritz
representation of a Green-confirmed peak is not evidence that the peak is unreal, only that the
discretised Ritz basis has no single clean eigenmode for it (near-degenerate crowding, a leaky
pseudo-SAW, or an unlucky truncation). The pipeline is therefore:

1. **Existence:** every candidate whose score clears the threshold (`matchGreenScorePeaks` /
   `matchGreenScoreLayeredPeaks`). Nothing after this point can remove a candidate.
2. **Match:** each accepted peak is matched to its best-supporting Ritz mode
   (`matchRitzToGreenPeak`; always succeeds). A low-confidence match is still reported, with the
   plain argmax mode.
3. **Label last** (section 7).
4. **Limiting-velocity entries** from `associateLimitingMode` (below).

No absolute `oopSurfDispl` or `matchScore` floor exists anywhere in this pipeline: each would
reject real modes that the Ritz basis represents poorly. Because a weak match is reported rather
than dropped, an entry's `modeVel` (the matched Ritz mode's velocity) is not guaranteed to lie within
`velMatchTol` of its `greenPeakVel`, and two entries can share a `modeNo`. `matchScore`
describes match quality; it says nothing about whether the entry exists.

The `'refined'` method identifies through the single-direction pole nearest its spot peak
(`matchRitzToGreenPeak`; `matchSpotPeak` when there is none); the `'ritz'`
method by the supporting mode of its own modal-sum peak (`localMode`).

### `matchRitzToGreenPeak(ritz, vGreen, obj, isLayered)` — the `'green'` match

Binds a `'green'` peak (`GreenScore` / `GreenScoreLayered`, this direction's pole-fitted
`G13_eff` peak) to the Ritz mode **nearest** its velocity. Among the modes within
`GREEN_TIE_REL` × vGreen (5e-4, ≈ 1.4 m/s at 2.8 km/s) of the nearest distance, the one with the
largest identity weight wins:

```
identityWeight:  half-space oopSurfDispl .* freqDisp,  layered oopSurfDispl
```

It returns `[modeNo, score, confident]`: `score` is `matchScoreOf` the chosen mode (below),
`confident = score > matchScoreThresh` (default 0.05). Every Ritz mode is a candidate, with no
velocity-based eligibility cut. `isLayered` is passed by the caller, not read from `ritz`, so
plain-struct `ritz` stand-ins work.

**Why nearest.** The `'green'` peak is the pole of the very direction the Ritz basis was solved
for, and it is placed to within a few m/s at most (against a fine, weakly damped `G13_eff` peak). A
trapped wave has a Ritz mode at its peak, so the nearest mode is exact there, including a Rayleigh
wave hugging `minShear` among dense bulk-like slab modes. A leaky (above-shear) wave is carried by
several slab modes 1–2 % of v apart, and a rule that weights proximity loosely against the surface
weight can bind a neighbour tens of m/s away. Within the peak's own placement error the nearer of
two modes is chance, hence the tie window: there the surface weight decides. The mode that
dominates the Ritz modal sum `|u3 conj(u1) / (v_n² − v²(1 + iη))|` at the peak is not used: it
matches leaky peaks well but hands weak (Love-type) layered peaks to a strong neighbour far away in
velocity.

**The entry's velocity is the Green peak, not the mode.** For a leaky entry the bound mode's
eigenvelocity (`modeVel`) can be up to half a slab-mode spacing from the true resonance, where
`greenPeakVel` is within a few m/s, so `'green'` / `'refined'` entries report
`vel = greenPeakVel`. All of them, trapped too: there the eigenvelocity is closer, but the
difference is below Green's own resolution, and one source per method keeps a branch continuous
where it crosses the shear limit. `'ritz'` entries keep their mode's eigenvelocity.

### `'refined'`: through the single-direction pole, `matchSpotPeak` as fallback

A spot peak is not a single-direction pole. `GreenSpotScore` sums `G13_eff` coherently over the
spot's directions θ ± 2.5°, so its peak sits near the spot-weighted average of the branch velocity,
moved by the branch's curvature, or by asymmetry where the branch fades on one side. The Ritz basis
is solved at θ alone, so its modes sit at θ's velocities. The spot peak can therefore lie tens of
m/s from the single-direction pole — as far as the slab-mode spacing — so any rule applied at the
spot velocity can bind a neighbour of the wave's mode on leaky entries.

`GreenSpotScore.identify` therefore maps each spot peak to the nearest single-direction pole of the
angle — a `GreenScore` candidate (`greenScorePeaks`, accepted or not, so no threshold is involved)
within max(`POLE_WIN_REL` v, `POLE_WIN_ABS`) = max(1.2 % v, 30 m/s), the same tolerance that pairs
`'green'` and `'refined'` entries (`refinedSawAgreement`) — and binds that pole with
`matchRitzToGreenPeak`.

With no pole in the window (rare; most often where the spot spectrum forms a comb, section 12)
`matchSpotPeak` decides:

```
matchScore   = velProximity .* identityWeight      (matchScoreOf; half-space)
velProximity = exp(-(relDist / velMatchTol)^2),  relDist = |vels - vSpot| / vSpot
```

`modeNo` is the argmax. Not nearest: beside a Rayleigh wave hugging `minShear` the spot peak falls
among dense bulk-like slab modes, and the surface weight is what finds the wave's mode.

`oopSurfDispl` and `freqDisp` are read from `ritz.diagnostics` (`RitzSolverSAW.modeDiagnostics`,
defined in `ritz.md`). On a half-space `freqDisp` separates a true SAW from bulk-like Ritz modes and
picks the correct mode more often than `depthInsensitivity`. On a layered stack every true guided
mode is dispersive, so `freqDisp` carries no information about which mode represents a peak and is
left out. Either way the identity only ranks Ritz modes for a peak that already exists; it is never
an existence floor.

### `associateLimitingMode(ritz, vLim, polLim, obj)` — limiting-velocity entries

A separate criterion: matches one Christoffel limiting velocity (shear or longitudinal) to the
best-supporting Ritz mode within a ±5-index window around the closest Ritz mode, scored by
`depthInsensitivity × polarisation projection × velProximity`, and reports no match unless the
score clears `matchScoreThresh` — here a genuine existence floor, since no Green peak backs a
limiting association. It is used for every non-degenerate limit: each shear root and the
longitudinal root (the substrate's roots on a layered stack). It associates the grazing mode at
each threshold only, not a general search for bulk modes. It uses the raw `ritz.depthInsensitivity` (log10 form, all modes),
not `ritz.diagnostics.depthInsensitivity`, which is already offset by its minimum over all modes;
re-offsetting it within a small window would double-process it.

**Degenerate limits.** Limits that share one velocity and have orthogonal polarisations (the two
shear roots at an acoustic axis) are grouped first (`degenerateLimitGroups`: velocities equal to
`LimitingVels.DEGENERATE_REL_TOL`, polarisations orthogonal to `LIMIT_ORTH_TOL` = 1e-6). Their
polarisations are only a basis of the degenerate eigenspace, so pairing each with a mode by its own
vector would be decided by that basis. `associateLimitingGroup` therefore associates each group
jointly: the g best modes above `matchScoreThresh`, scored with the projection on the group's span,
are dealt out in mode order — g modes for g roots, never one mode twice.

The limiting entries need the eigenvectors, so a Calc object associates them once per calculation
and every later `'green'` or `'refined'` classification carries them unchanged
(`associateLimiting = false`, `carryEntries`; section 8). The `'ritz'` classification has none.

### Why proximity scoring is Gaussian, not linear

Both `matchSpotPeak` and `associateLimitingMode` score with a Gaussian `velProximity` term
(width `velMatchTol`), not a linear one (`1 − relDist`). A linear term is too weak: a strong but
unrelated match component (`oopSurfDispl × freqDisp`, `oopSurfDispl`, or
`depthInsensitivity × polarisation projection`) can then win the argmax for a candidate several
per cent away in velocity. `associateLimitingMode` additionally requires the combined score to
clear `matchScoreThresh`, not merely be positive: a window centred on the nearest Ritz mode *by
index* can otherwise land on the dominant Rayleigh/SAW mode itself, whose high
`depthInsensitivity` and nonzero polarisation projection would win on magnitude alone. If either
safeguard is weakened, expect a spurious duplicate or mislabelled entry co-located with a genuine
dominant mode, not a crash: both failures are silent under every check except a multi-angle scan
or a visual diagnostic.

## 7. Labels

### One search range, one deciding quantity

- **One contiguous search range, with no split at `minShear`.** MATLAB's `findpeaks` recognises
  only an *interior* maximum — one with a lower neighbour on both sides within the array it is
  given. Near some directions a genuine Rayleigh mode crowds right up against the shear limit,
  giving a sharp, narrow feature just below `minShear`; a search boundary there would cut off its
  falling side and make it structurally invisible. Every method therefore searches one range
  (`GreenScore`: the band up to vTf(1 + 0.012); `GreenSpotScore`: up to vL, capped at
  vTf(1 + 0.012); `RitzScore`: up to vL; layered: [0.3 vSmin, 1.05 vLmax]), and every search runs a
  little past its bound so that a peak at the bound is interior.
- **The label from one quantity: the peak's own velocity.** Ritz eigenvalues are variational
  approximations and Christoffel's `minShear` is exact, so the two never agree to numerical
  precision; near the boundary a Ritz velocity can sit a hair on the wrong side of `minShear`. A
  design that decides the side from Ritz's velocity in one place and the Green peak's in another
  can mislabel a peak, or lose it entirely when the two disagree in opposite directions. On the
  Green methods the Green peak velocity, which carries no Ritz discretisation error, alone decides
  the label, and no velocity-based cut restricts which Ritz modes may match.

The general lesson: when two solvers each supply a number meant to represent the same physical
threshold, derive the "which side are we on" decision from one canonical, most exact quantity and
reuse it everywhere.

### The boundary per method

| Method | `'Rayleigh'` if | Tolerance |
|---|---|---|
| `'green'` | `greenPeakVel < minShear·(1 + RAYLEIGH_LABEL_TOL)`, `RAYLEIGH_LABEL_TOL` = 3e-3 | Green's resolution (below) |
| `'refined'` | spot peak velocity `< minShearSpot·(1 + LABEL_TOL)`, `LABEL_TOL` = 7e-3 | the method's own η; minShearSpot = the spot-weighted minShear |
| `'ritz'` | supporting mode's eigenvelocity `< minShear + slabRes` | the slab's resolution of the slow-shear limit (below) |
| layered (both) | `'SAW'`; `aboveSubstrateST` = velocity above min(substrate shear limits) | none (a flag only) |

**Green's resolution.** The damping η = 3e-3 of the existence solve broadens every pole to a full
width of about η·v, and the branch point at `minShear` pulls the |G| maximum as well, so Green
cannot place a peak on either side of `minShear` more finely than about η·`minShear`. A Rayleigh
wave hugging `minShear` (on a strongly anisotropic test material) puts its peak up to about 1 m/s
above it, and without the tolerance its label flickers with sub-m/s interpolation noise; the
nearest genuine pseudoRayleigh above `minShear` on the cuts tried is tens of m/s away. A leaky
wave within 0.3 % above `minShear` would also be called Rayleigh; the expectation is that a
branch's leakage vanishes as it approaches its limiting velocity, so the difference is a name only
there — not verified in general.

**The slab's resolution** (`RitzScore.slabResolution`). The slab cannot represent the grazing
slow-shear limit exactly: its lowest standing pattern pairs two partial waves whose `k_z` differ by
π/D (rigid bottom, quarter wave), a window of width `1/(2D)` in `q = k_z/k`, whose trace velocity
is at least the larger of the two partial waves'. Ritz eigenvalues are upper bounds, so the error
is one-sided: a Rayleigh mode hugging `minShear` can sit up to that much above it (0.2–0.6 m/s on
the cuts tried; computed from Christoffel alone, < 1 ms). The modal-sum *peak* is not used: on the
2 m/s axis it jitters by 1–3 m/s, which flips labels and would send a vS-hugging Rayleigh wave
through the pseudo-SAW gates.

### `minShear = min(shearLimVels)` — why the plain minimum, regardless of branch, is correct

`shearLimVels` is not "one value per shear branch": `CalcCommon.limitingAndBulk`'s sagittal-plane
scan (`LimitingVels.findRoots`) returns every accepted root per branch — `islocalmin` deliberately
catches tangential/grazing touches as well as crossings — so `numel(shearLimVels)` is often more
than 2 and varies from angle to angle.

This is physical, not a numerical defect. Each limiting velocity is the group velocity of a bulk
state confined to the sagittal plane, and its value follows from a short argument. Let `s = (s1, s3)`
be the slowness of a bulk wave in the sagittal plane; its trace velocity along the surface is `1/s1`.
At a maximum of `s1` along a slowness curve the tangent is parallel to `x3`, so the curve's normal
points along `x1`. The group velocity is normal to the slowness surface (Every et al. 2013, §2.3,
p. 1201, Eqs. 7–8), so it is parallel to the surface, and `v_g · s = 1` (same equations) gives
`v_g = 1/s1,max`: that bulk state carries energy along the surface at exactly the limiting
velocity. The scan varies direction within the plane containing the propagation direction and the
surface normal, not a single fixed k-vector. A shear slowness sheet is generically concave in places
(saddle-shaped regions of negative Gaussian curvature, Every et al. 2013, p. 1201; the origin of
phonon-focusing caustics), so it can have several such stationary points, whereas the longitudinal
sheet is convex and has exactly one — hence `longLimVel` has exactly one root and `shearLimVels`
generally does not. Multiple roots per branch are genuine surface-skimming tangential touches.

**The slowest of them is, by definition, the existence boundary for a non-leaky surface wave:** a
surface wave cannot have a phase velocity above *any* channel through which the bulk can carry
energy away, whichever shear branch that channel belongs to. Branch identity (`shearLimIwt`)
matters for a different question — following one named branch from one calculation to the next.
`minShear = min(shearLimVels)` therefore stays the plain minimum and does not use `shearLimIwt`.

### Why Green owns the label when it is available

Distinguishing a bound Rayleigh mode from a leaky pseudo-Rayleigh one is by definition a comparison
against the exact bulk threshold. Ritz has no radiation boundary condition and so no notion of that
threshold, which is why the `'ritz'` method runs the (cheap) Christoffel step and labels against it
with the slab's own tolerance.

### No Love type, no Sezawa split

No `'Love'` type exists on any path. Love/SH motion has little out-of-plane displacement, so it
couples only weakly into `G13` (and not at all into the Ritz numerator `u3·conj(u1)` for pure SH);
on a layered stack a Love-type wave is a `'SAW'` entry when the standard admits it. Any other Green
component is available on request (`greenComponent(component)`, `greenComponentMap`;
`forwardCalcs.md`). `'other'` exists in the type taxonomy for completeness; nothing emits it.

Every accepted layered candidate is a single `'SAW'` type, with no `'Rayleigh'`/`'Sezawa'` split by
velocity rank: off a principal plane, "lowest velocity" and "tallest, most observable peak"
decouple (the lowest-velocity candidate can be far weaker than the dominant peak). `vel`,
`score`, `relPeakHeight`, and `aboveSubstrateST` are reported per entry so a caller can judge branch
order, dominance and leakiness from the numbers.

## 8. `classify` / `reclassify`, stored fields, the stamp

### On the Calc classes

| Call / property | `CalcSingleHS`, `CalcAngSweepHS` | `CalcSingleLayered`, `CalcAngSweepLayered`, `CalcWavevecAngSweep2Layered` |
|---|---|---|
| `sawClassMethod` (set before `calc`) | `''` (default: `'green'` on `'full'`, `'ritz'` on `'ritzOnly'`, worked out on every run; `classify` leaves it) \| `'green'` \| `'refined'` \| `'ritz'`; on `'ritzOnly'` only `'ritz'` | `''` (default, as half-space) \| `'green'` \| `'ritz'` |
| `sawThreshold` | `[]` = the method's default (`'high'` for `'green'` / `'refined'`, `'low'` for `'ritz'`); a number needs `sawClassMethod` set | `[]` = `'low'` for both; a number needs `sawClassMethod` set |
| `classify(method, threshold)` | sets the active classification; both arguments optional (the active method, the method's default standard) | the same; `CalcWavevecAngSweep2Layered` per wavelength |
| `reclassify(threshold)` | `classify(activeMethod, threshold)` | the same |
| `isStale()` | single and sweep classes | `CalcSingleLayered`, `CalcAngSweepLayered` |

`calcAll('full')` (and a sweep's `calc`) classifies `'green'` first — it associates the limiting
entries while the eigenvectors are in memory — and then the chosen method if another. `classify`
errors with `<Class>:classify:needsFull` for `'green'` / `'refined'` on a `'ritzOnly'` result, and
`CalcCommon.validateSawClassMethod` with `:noRefinedLayered` for `'refined'` on a layered stack.
`'refined'` solves on first use: 11 Green solves per angle, or, on a sweep, every direction the
angles' spot kernels need, once (the reference directions included), in parallel per
`useParallel`.

What every viewer, Vis function, tracker, and drill-down reads — `classification.significantModes`,
the sweep's `significantModesMap` / `sawVelMap` / `sawModeNoMap` / `sawTypeMap` — holds the active
classification. `activeMethod` and `activeThreshold` record it.

Stored per angle (all threshold-free, all surviving eigenvector pruning):

| Field | Content |
|---|---|
| `greenScoreVelAxis`, `greenScoreG13Eff` | the band of the Green axis and the complex `G13_eff` on it |
| `greenScorePeaks`, `greenScoreRef` | the scored `GreenScore` / `GreenScoreLayered` candidates; the half-space `'low'` reference (also on the sweep) |
| `ritzScoreNum`, `ritzScoreVelAxis`, `ritzScoreRef`, `ritzScoreCands` | the modal numerators (the Ritz spectrum), the score axis, the own floor (+ `slabRes` on a half-space), the scored Ritz candidates |
| `refinedSaw` (half-space) | the `GreenSpotScore` candidates with identity and label, their reference and info (incl. `minShearSpot`); `[]` until `classify('refined')`; the sweep keeps `refinedSawRef` |
| `limitingEntries`, `limitingDone` | the limiting entries, associated once per calculation |
| `classification`, `activeMethod`, `activeThreshold`, `classificationVersion` | the active classification and its stamp |

### The `ModeClassification` result

| Property | Content |
|---|---|
| `mode` | `'ritzOnly'` \| `'full'` \| `'layered'` — the path that ran |
| `method` | `'green'` \| `'refined'` \| `'ritz'` |
| `threshold` | the standard it was run at |
| `significantModes` | the entries (below) |
| `greenScorePeaks` | `'green'`: every candidate with its parts and `accepted` (half-space: `vel`, `score`, `lowScore`, `rawScore`, `gateTerm`, `dephaseTerm`, `gamma`, `slope`, `relWidth`, `relHeight`; layered: the `GreenScoreLayered.emptyPeaks` fields). Empty otherwise |
| `ritzScorePeaks` | `'ritz'`: every candidate (section 4). Empty otherwise |

### The `significantModes` entries

One struct shape on every path — `type`, `vel`, `modeNo`, `modeVel`, `matchScore`, `greenPeakVel`,
`oopDisplScore`, `dispersionScore`, `limVelAssoc`, `aboveSubstrateST`, `relPeakHeight`,
`method`, `score`, `scoreDb`, `grade` — with fields that do not apply set to
`NaN` (or `''`) rather than omitted. MATLAB struct-array concatenation (`[sm1, sm2]`) requires
identical field sets and does not backfill a missing field, so a ragged field set would break the
sweep classes' map assembly. `emptyEntries`, `makeSAWEntry`, and `makeLimitingEntry` (local
functions of `ModeClassification.m`) and the `significantModes` property default hold the one
field set.

| Field | `'green'` | `'refined'` | `'ritz'` |
|---|---|---|---|
| `type` | half-space `'Rayleigh'` / `'pseudoRayleigh'`, plus `'limiting'`; layered `'SAW'`, plus `'limiting'` | `'Rayleigh'` / `'pseudoRayleigh'`, plus `'limiting'` | half-space `'Rayleigh'` / `'pseudoRayleigh'`; layered `'SAW'` |
| `vel` | = `greenPeakVel` | = `greenPeakVel` (the spot peak) | the supporting mode's eigenvelocity |
| `modeNo`, `modeVel` | the matched Ritz mode and its eigenvelocity | the same | the supporting Ritz mode (`localMode`) |
| `greenPeakVel` | the `G13_eff` peak velocity | the spot peak velocity | `NaN` (the peak is in `ritzScorePeaks`) |
| `score` | `GreenScore` score / `GreenScoreLayered` score | scoreHigh | half-space dScore; layered dScore + polTerm |
| `matchScore` | `matchScoreOf` the matched mode | the same | half-space the supporting mode's `freqDisp`; layered its share of the modal sum |
| `relPeakHeight` | peak height / tallest candidate at this angle | the same, on the spot spectrum | modal-sum peak height / tallest candidate |

`scoreDb = 20·score` (the scores are log10 amplitude ratios); `grade` is `'strong'` / `'weak'`
(section 1); `method` is `''` and `score` / `scoreDb` / `grade` empty or `NaN` on limiting entries.
`oopDisplScore` / `dispersionScore` are the mode's `oopSurfDispl` / `depthInsensitivity`.
`aboveSubstrateST` is `true` / `false` on layered `'SAW'` entries, `NaN` elsewhere.

`relPeakHeight` is display-only, in [0, 1] (1 = the tallest feature at that angle). Existence does
not establish dominance: a Rayleigh peak can clear the threshold for several degrees after a much
taller pseudoRayleigh peak has become the observable response, and `relPeakHeight` lets a caller tell
the two apart; the classification plots size their markers by `score` instead. It is never a
gate.

### The classification stamp

`ModeClassification.CLASSIFICATION_VERSION` is raised whenever a classification rule
or a stored candidate field changes. Every classification records it on the Calc object
(`classificationVersion`; 0 = unclassified). A result whose stamp differs is **stale**
(`isStale()`, `CalcCommon.isStale`): it still loads and shows its maps, modes, and Green and Ritz data,
but `classify` / `reclassify` refuse it (`<Class>:classify:stale`) — recompute it.

## 9. Parameters

`ModeClassification` properties (name-value overrides in the constructor; the Calc classes pass
`threshold` and, on a reclassification, `associateLimiting` / `carryEntries`):

| Property | Default | Meaning |
|---|---|---|
| `threshold` | `'high'` | `'high'` \| `'low'` \| numeric, read on the scale of the score supplied. The Calc classes pass the method's default |
| `velMatchTol` | 0.02 | relative tolerance for "close"; the Gaussian width of `matchScoreOf` (`matchSpotPeak`'s ranking, every entry's `matchScore`) and of `associateLimitingMode` |
| `matchScoreThresh` | 0.05 | Green and spot matches: sets `confident` only; `associateLimitingMode`: an existence floor |
| `associateLimiting`, `carryEntries` | `true`, `struct([])` | skip the limiting association and carry the given limiting entries over |

`ModeClassification` constants: `RAYLEIGH_LABEL_TOL` = 3e-3, `GREEN_TIE_REL` = 5e-4,
`LIMIT_ORTH_TOL` = 1e-6 (section 6), `CLASSIFICATION_VERSION`.

The scoring classes keep their calibrated values as `Constant` properties — calibrated together
and not meant to be tuned per call:

| Class | Constants |
|---|---|
| `SawSpectrum` | `ETA` 3e-3, `GATE_PERIODS` 3, `RAMP_PERIODS` 0.6, `MAX_REL_WIDTH` 0.12, `VEL_MIN` 0.3 km/s, `VEL_STEP` 0.004 km/s, `PAD_FACTOR` 4, `TAIL_TAPER` 0.2, `EDGE_TAPER` 0.1 km/s, `AXIS_HI_X_VL` 1.05 (half-space axes end at 1.05 × longLimVel) |
| `GreenScore` | `GATE_PERIODS` (= `SawSpectrum.GATE_PERIODS`), `APERTURE_SIGMA` 0.4°, `MAX_REL_WIDTH` (= `SawSpectrum.MAX_REL_WIDTH`), `BAND_LO` 0.80, `BAND_SOLVE_HI` 1.03, `SEARCH_HI_TOL` 0.012, `FIT_SPAN` 0.016 km/s, `PROBE_OFFSET` 0.5°, `PROBE_STEP` 0.016 km/s, `PROBE_HALF` 3, `DEFAULT_STEP` 0.006 km/s, `MAX_STEP` 0.016 km/s, `THRESHOLD_HIGH` 0.16, `THRESHOLD_LOW` −0.79, `N_REF_ANGLES` 12, `REF_CACHE_SIZE` 32 |
| `GreenSpotScore` | `ETA` 7e-3, `VEL_STEP` 0.008 km/s, `FWHM` 2.5°, `FWHM_TO_SIGMA` 2.3548, `TRUNC_SIGMA` 2.5, `DIR_STEP` 0.5°, `GATE_DETECT` 4, `GATE_POS` 2, `BAND_LO` (= `GreenScore.BAND_LO`), `BAND_HI` 1.00, `SEARCH_HI_TOL` (= `GreenScore.SEARCH_HI_TOL`), `THRESHOLD_HIGH` −0.125, `THRESHOLD_LOW` −0.550, `REF_ANGLES` 0:15:165°, `REF_CACHE_SIZE` 32, `LABEL_TOL` 7e-3, `POLE_WIN_REL` 0.012, `POLE_WIN_ABS` 0.030 km/s |
| `RitzScore` | `VEL_STEP` 0.002 km/s, `ETA` 1e-2, `MODE_VEL_FACTOR` 1.2, `MIN_FREQDISP` 0.3; half-space `HS_THRESHOLD_HIGH` 0.16, `HS_THRESHOLD_LOW` 0.70, `HS_MIN_WINDOW_MASS` −1.55, `HS_WINDOW_DEPTH_REF` 30, `HS_BOTTOM10_MAX` 0.15, `HS_BOTTOM_RATIO_MAX` 0.8, `HS_MIN_LEG_PER_DEPTH` 2.3, `HS_FLOOR_LO` 0.8, `WINDOW_REL` 0.012, `WINDOW_ABS` 0.030 km/s, `WINDOW_TOP` 0.98, `LOCAL_MIN_HALF` 0.004 km/s; layered `LAY_THRESHOLD_LOW` 0.65, `LAY_THRESHOLD_HIGH` −0.7, `LAY_K_POL` 0.75, `LAY_EBOT_MAX` 1e-2 |
| `GreenScoreLayered` | `THRESHOLD_LOW` 0.70, `THRESHOLD_HIGH` 0.20, `K_POL` 0.75, `GATE_OWN_PERIODS` 10, `BAND_LO_X_VSMIN` 0.3, `BAND_HI_X_VLMAX` 1.05 |

## 10. Cost

Serial, one thread, per angle unless stated; a strongly anisotropic cut on the half-space, NiTi on
MgO on the layered stack. Times are machine-dependent; read them as orders of magnitude.

| Step | Half-space | Layered |
|---|---|---|
| Ritz solve (baseline basis) | ~30 ms standalone | — |
| The angle's one Green solve + `'green'` scoring (6 m/s axis; probes ~7 points per candidate per side) | ≈ 46 ms | ≈ 0.13 s solve (~2450 points, `vMax` 10 km/s) + ≈ 0.08 s scoring (probes and one full-tensor solve of 3 points per candidate) |
| `GreenScore.reference` (once per structure, memoised) | ≈ 0.56 s | none |
| Ritz candidates (`calcRitzScore`) | ≈ 12 ms on `'ritzOnly'`; +~17 ms per angle on `'full'` | ≈ 9 ms |
| `'refined'` | 11 Green solves (~0.4 s) per standalone angle, plus the memoised reference (up to 132 solves) on a structure's first use; a 0:1:180 sweep needs 360 directions (twice the Green solves of a 1° map), about half the time of the `'full'` sweep itself | — |
| Re-thresholding (`classify` to a computed method, `reclassify`) | no solve | no solve |

A 181-angle `'full'` half-space sweep (Ritz, Green, `'green'`, and the Ritz candidates) takes of
the order of 8 s on 12 threads. The `'refined'` direction cache lives only during the call; a
result stores the scored candidates, not the spectra.

## 11. Validation

Measured TGS maps against which the classification's accepted points can be checked are
published for:

- NiTi/MgO, a layered stack: Grabec et al., "Guided acoustic waves in thin epitaxial films:
  Experiment and inverse problem solution for NiTi", *Ultrasonics* **138**, 107211 (2024),
  doi:10.1016/j.ultras.2023.107211.
- Half-space cuts of Ni(110), Ni(311), and Fe3Al(316): Grabec, Stoklasová, Repček, Kušnír, Mareš, Ševčík,
  Sedlák, and Seiner, "Ultra-transient grating spectroscopy for visualization of surface
  acoustics", *Nat. Commun.* (2026, accepted), arXiv:2510.10696.

The same comparison can be repeated with GuliweRR: compute an angular sweep of a published cut
(`Calc.angSweepHS`) and overlay `classification.significantModes` on the measured map. Each method
is run exactly as a user would run it, at its two named standards.

## 12. Limitations

- **Love-type (shear-horizontal) modes** couple only weakly into `G13`. On a half-space no Love
  type exists. On a layered stack both methods rank Love waves below the sagittal ones on purpose
  (polarisation term): `'low'` excludes them, `'high'` finds most of them at the cost of false
  positives.
- **Branch ends** are a threshold on a continuum: where a branch fades over several degrees,
  the classified end depends on the threshold. Where a branch fades over several degrees its score
  decays smoothly; a finer axis does not reveal a peak further out.
- **`'ritz'` depends on the basis.** Its gates are calibrated on the baseline Ritz basis
  (relDepth ≈ 30, `legOrder ≥ 2.3 relDepth`, stricter than the general `legOrder ≥ 2 relDepth` of
  `docs/pitfalls.md`, MTH1); finer bases change the scores, and the bottom-trapped gate rests on
  one cut. A very leaky pseudo-SAW that shows in TGS can be dropped at `'high'`; a weak Rayleigh
  wave hugging `minShear` on a steep branch can be absent at `'low'`. On a strongly anisotropic
  test material a Rayleigh wave hugging `minShear` keeps a stable label, but at `'low'` only about
  half of that branch is kept (strong dephasing).
- **`'refined'`: the 0.5° kernel is a quadrature.** Where a narrow pole's branch moves more per
  0.5° than the pole's width (|dv/dθ| · 0.5° > η·v), the 11 directions' poles no longer merge into
  one smeared peak but appear as a comb of separate peaks, each at a kernel direction's own
  velocity. Seen only on a strongly anisotropic test material, where some refined Rayleigh points
  are labelled pseudoRayleigh. The steepest narrow branches on the cuts tried sit near that limit.
- **Layered tuning rests on one measured stack** (NiTi on MgO, Grabec et al. 2024: `K_POL`, `LAY_K_POL`, and every
  layered threshold). A strongly leaky sagittal wave between the substrate shear and longitudinal
  limits can pass the `'green'` `'low'`; untested against measurement. An accepted entry with a
  large pole width (`gamma` many times `SawSpectrum.ETA·v/2`) above the substrate shear limit
  deserves care until a measured case covers it: the own-period gate charges it by its width
  only. The layered `'ritz'` never reports leaky waves (trapped gate).
- **`reclassify` to a looser standard with `retentionLevel = 'saw_only'`**: modes newly accepted
  have no stored eigenvector.

## `DispersionBranches.track` — a separate utility, not part of the classification

`DispersionBranches.track(sweepResult, angleIdx)`, a static method of its own class, unrelated to
`ModeClassification`. It tracks each Ritz mode's
*physical identity* across a fixed-angle wavelength/wavevector sweep
(`CalcWavevecAngSweep2Layered` or `CalcWavevecAngSweep2Plate`, which expose the same
`sweeps{iw}.cs{angleIdx}.ritzSolver` shape — see `docs/reference/architecture/plate.md`), for
`numel(angles)==1` or an explicit `angleIdx`, by maximum eigenvector overlap (MAC — modal assurance
criterion, here mass-weighted, see below) with the previous wavelength step, so a plotted
dispersion curve does not jump between unrelated Ritz mode indices. It requires `retentionLevel='full'` on the input sweep: a branch can drift
through classification wobble between steps, so only the full eigenvector matrix (not just the
columns `'saw_only'` keeps) is enough to track it.

One `[nCompModes × nCompModes]` MAC matrix product per step; branches are assigned to candidates
one-to-one by maximising the total MAC over the step (`matchpairs`, a linear assignment), so every
row of `branchModeNoMap` is a permutation of `1:nCompModes`: each raw mode sits on exactly one
branch.

**Mass-weighted MAC.** The overlap is taken on the unconstrained Legendre coefficients `N·α`,
each scaled by the square root of its layer's mass factor (`GammaMatrix.LayerNormFac`, relDepth·ρ
relative to layer 1 — the diagonal of the unconstrained `M`), so it measures the overlap in the
physical inner product ∫ρ u^H u dz, as `M_con` does. Each step uses its own factors, so the
cross-step product carries `sqrt(M_a M_b)`: symmetric in the two steps, and at most 1 without
further normalisation (the eigenvectors are `M`-orthonormal). A plain Euclidean overlap would
weight every layer alike however thin it was: a thin film's coefficients can hold most of a
Rayleigh eigenvector's Euclidean norm but a tiny fraction of its mass norm, so unrelated film modes
reach a MAC near 1 and the one-to-one assignment can hand the Rayleigh branch to a film mode.

**Why one-to-one.** An independent argmax per branch would let several branches claim one
candidate, and a candidate no branch claimed would then fall off every branch for the rest of the
sweep. In a dense, near-degenerate cluster — a layered stack's substrate box modes above its bulk
shear velocity — that happens at almost every step, so a guided mode emerging from the cluster as
k grows would never be plotted. Where an argmax gives distinct matches the
assignment is the same.

**Tracking direction does not matter.** The reference at each step is the previous step's own
eigenvectors, only relabelled, so the pairing between two adjacent steps depends on their MAC
matrix alone (symmetric in the two steps, see above); tracking backwards uses its transpose, whose optimal one-to-one assignment is the
same pairing. Starting from the high-k end gives the same curves with different branch numbers.

**One seeded branch.** `DispersionBranches.trackSeededBranch(sweepResult, angleIdx)` (layered sweeps)
picks the branch through the slowest accepted `'SAW'` entry of the first wavelength step (its
`modeNo`) and returns that branch's velocity at every step — the seed convention the tutorials and
examples use, in one place.

What the tracker cannot do is resolve a pairing the overlaps themselves leave ambiguous: a step too
coarse for a near-degenerate pair, or modes inside the substrate box-mode cluster, which are not
physical branches anyway. Finer steps are the remedy there, not a different direction.

## `GreenPeaks.find(absG, velAxis, mask, heightFrac)` — generic peak finder

Peak-finds a spectrum inside a velocity sub-range (`mask`) via MATLAB's `findpeaks` (Signal
Processing Toolbox), accepting a peak only if its height exceeds `heightFrac × max(absG)` (default
`GreenPeaks.HEIGHT_FRAC` = 0.15, of the whole array, not just the mask), and reports each peak's
`prominence` and half-prominence `width` without filtering on them. It is generic over its
`(absG, xAxis)` inputs — no dependency on Ritz, Christoffel, or the half-space/layered
distinction. An isoFreq spectrum can be searched with `kParRange` in place of `velAxis`. The
classification does not use it.

## References

`forwardCalcs.md` (the orchestration layer that drives this class: `identifyModes()`, `classify`,
`greenStep` / `ritzOnlyStep`, the shared references, the sweep classes' `sawTypeMap` filters, and
"The Green velocity axis"); `green.md` (the `eta` argument, "±k symmetrisation"); `ritz.md`
(`modeDiagnostics`, `angleSlope`, the Ritz score, and the Ritz spectrum); `docs/pitfalls.md`
(classification gotchas); `visual.md` (`plotClassificationDiagnostic` /
`plotClassificationCurve`'s fixed draw order for limiting versus mode points, and
`plotRitzSpectrumCurve` / `plotRitzSpectrumMap`).

- Maznev, A.A., Akthakul, A., and Nelson, K.A., J. Appl. Phys. 86(5), 2818–2824 (1999),
  doi:10.1063/1.371130 (TGS as a horizontal surface dipole; ±q symmetry of G13, Eq. 7, p. 2820).
- Every, A.G., Maznev, A.A., Grill, W., Pluta, M., Comins, J.D., Wright, O.B., Matsuda, O.,
  Sachse, W., and Wolfe, J.P., "Bulk and surface acoustic wave phenomena in crystals: Observation
  and interpretation", Wave Motion 50, 1197–1217 (2013), doi:10.1016/j.wavemoti.2013.02.007 — for
  why shear slowness surfaces are concave in places (caustics) while the longitudinal one is convex.
