# Plate (Lamb waves) — architecture and design rationale

`guliwerr/solvers/ritz/RitzSolverLamb.m` (solver) + `guliwerr/forwardCalcs/CalcSinglePlate.m`, `CalcAngSweepPlate.m`, `CalcWavevecAngSweep2Plate.m` (orchestration, `Calc.singlePlate`/`.angSweepPlate`/`.wavevecAngSweep2Plate`). Computes Lamb-mode (free-plate) velocities, group velocities, mode types, and zero-group-velocity (ZGV) points for a homogeneous or multi-layered anisotropic plate via the same Ritz–Rayleigh variational method `docs/reference/architecture/ritz.md` describes for SAW — this document covers only what is specific to the free-plate boundary condition; read `ritz.md` first for the shared method, `K`/`M` assembly, and Legendre-basis mechanics.

## Why Plate is a separate, leaner class family than HS/Layered

`CalcSinglePlate`/`CalcAngSweepPlate`/`CalcWavevecAngSweep2Plate` mirror `CalcSingleHS`/`CalcAngSweepHS`/`CalcWavevecAngSweep2Layered`'s composition pattern (a sweep is built from per-angle containers, a 2D sweep from per-wavevector angle-sweeps) but have no Christoffel/Green phase at all, for two physical reasons, not omissions:

1. **No Green's-function analogue exists in this codebase.** `GreenSolver` models a half-space or layered-substrate's surface response to a harmonic point force — a different, unbuilt calculation for a *finite* free plate. There is consequently no `'full'` calcMode for Plate — every container in this family only ever runs the Ritz solve.
2. **Every converged Ritz mode is a genuine physical branch.** A free plate's boundary conditions (both surfaces traction-free) are the *true* physical ones, not an artificial truncation of an infinite domain the way SAW's rigid-bottom BC is. There is consequently no "which numerical mode is the real one" question for `ModeClassification`'s `significantModes` machinery to answer — mode *type* (Asymm/Symm/TH-Asymm/TH-Symm; TH marks the shear-horizontal families) is a different question, answered directly by `RitzSolverLamb.modeRecognition` (see below), not by `ModeClassification.calc`. `CalcSinglePlate` never calls it.

**Christoffel limiting velocities are still relevant, despite no Green's function.** Bulk-wave velocities (`LimitingVels`/`Christoffel.matrix`, the same machinery `CalcSingleHS.calcLim`/`CalcSingleLayered.calcLim` use) bound Lamb-mode phase velocities in physically meaningful ways — e.g. at high `k·d` (wavevector times plate thickness; written kd below) a plate's surface modes approach the material's Rayleigh velocity, itself just below the slowest bulk shear velocity, independently of any Green's-function response calculation. `Vis.plotDispersionCurvesPlate`'s `showChristoffelBounds` option draws this: the bottom-layer and (if different) top-layer longitudinal and slowest-shear limiting velocities, as horizontal reference lines — deliberately both surfaces, not one, since a genuinely asymmetric stack's two faces bound its high-kd behaviour independently (see "Why a real (non-mirror-symmetric) multi-layer stack breaks the clean picture" below).

`retentionLevel` is `'full'`/`'none'` only for this family — no `'saw_only'` tier, since (per point 2 above) there is no "significant subset" to prune to; every mode is real. Physical scaling (`realThickness`, `realWavVec`) is `RitzSolverLamb`'s own constructor concept (rescales every layer's `relDepth` by a common factor so that `Σ relDepth = realWavVec·d_total/(2π)`, i.e. `relDepth = d/λ` with `λ = 2π/realWavVec`, or equivalently `k·d_total = 2π·Σ relDepth` with `k = realWavVec`), not `Structure.ApplyWavelength` — the 2D sweep class (`CalcWavevecAngSweep2Plate`) therefore sweeps `realWavVec` directly at one fixed `realThickness`, rather than via Layered's `ApplyWavelength`.

**Precompute sharing differs from the SAW/Layered family in one concrete way:** only `nullConstrMat` and `rawByDepthCell` (the Legendre integrals of `GammaMatrix.rawByDepth`; both depend on `legOrder` alone) are shareable across *both* the angle and wavevector axes. `GammaMatrix.elastByDepth`'s output (the Lamb-path analogue of `derivsByDepth`'s `dK/dc`) bakes in the *rotated* elastic tensor at a specific angle, so it is not angle-shareable — only wavevector-shareable at one fixed angle. `CalcAngSweepPlate` builds `rawByDepthCell` once and hands it to every angle (in MATLAB `CalcWavevecAngSweep2Plate` therefore builds it once per wavevector; the Julia port, whose grid is one flat loop, once per sweep). Each build saved is about 40 % of a single-layer Lamb solve at legOrder 20.

## Mode-type classification: Asymm/Symm/TH-Asymm/TH-Symm

`RitzSolverLamb.modeRecognition` (run by `calc(true)`, i.e. `withModeTypes`) labels every computed mode by **mid-plane parity**: each displacement component is split into its even and odd parts about the plate's mid-plane, and the mode's kinetic energy is divided into four fractions (`RitzSolverLamb.parityFractions`, rows in `modesType` order):

| `modesType` | label | energy fraction |
|---|---|---|
| 1 | Asymm (A0, A1, … — flexural family) | u1 odd + u3 even |
| 2 | Symm (S0, S1, … — extensional family) | u1 even + u3 odd |
| 3 | TH-Asymm (SH1, SH3, …) | u2 odd |
| 4 | TH-Symm (SH0, SH2, …) | u2 even |

The four fractions sum to 1; `modesType` is the largest, **`modeTypeConfidence` = largest − second largest**, in `[0,1]` (1 = a pure type, low = a mixed mode).

**Why parity.** A free plate whose stack is mirror-symmetric about its mid-plane (and whose crystals have a mirror plane parallel to the surface — cubic (001)/(110), isotropic, …) has the mirror `z → D − z` as an exact symmetry, so every mode is either symmetric (u1, u2 even, u3 odd) or antisymmetric (u1, u2 odd, u3 even). The split is a physical definition of the four classical labels, valid for every mode, higher ones included.

**SH-Lamb coupling.** In an anisotropic plate off the symmetry axes the SH and Lamb motions couple. The symmetry classes are then {u1 even, u2 even, u3 odd} and {u1 odd, u2 odd, u3 even}, so a symmetric mode shares its energy between Symm and TH-Symm (never an antisymmetric class while the mirror holds): the label names the dominant polarisation and a low confidence says the mode is genuinely mixed (e.g. S0 and SH0 of a Cu(001) plate propagating 25° off the cube axis are mixed, with a confidence clearly below 1).

**Asymmetric stacks.** Without a mirror plane, parity is not conserved and the fractions measure how far a mode is from either class — e.g. a wave bound to one face of a bilayer is half even, half odd, confidence ≈ 0. This is the meaning `classifyBranches` relies on (below).

**How it is computed** (cheap compared with the Lamb solve itself):
- **One layer.** The normalised Legendre polynomials are orthonormal on the layer and `P_n(−x) = (−1)^n P_n(x)` about its centre, which is the plate's mid-plane, so the even/odd energies are the sums of `|coefficient|²` over even/odd `n`, read off the eigenvector. No depth grid.
- **Several layers.** The displacement is evaluated at depth `z` and at its mirror `D − z` on a Gauss–Legendre grid over the lower half of the plate, whose sub-intervals end at every interface and every mirrored interface. Inside a sub-interval both are single-layer polynomials of degree `< max(legOrder)`, so `max(legOrder)` nodes per sub-interval integrate the even/odd energies exactly. The weight is the symmetrised density `(ρ(z) + ρ(D − z))/2`: even and odd parts are orthogonal under any mid-plane-even weight, so the fractions sum to 1 exactly, and it is the kinetic-energy weight itself whenever the stack is mirror-symmetric. One general path was preferred over an exact layer-mirror map for mirror-symmetric stacks: it costs the same order and covers every stack. A uniform plate split into two unequal layers reproduces the single-layer labels at confidence 1.

**Exact degeneracies.** Where two modes of different classes have the same velocity (e.g. an isotropic plate at `k·d = π`, where the Lamé mode and SH1 both travel at `√2 vT`, with `vT` the transverse bulk velocity), any mixture is an eigenmode; the eigensolver returns an arbitrary one and its confidence is low. That is correct, not a bug.

### Why a real (non-mirror-symmetric) multi-layer stack breaks the clean picture

The parity split above is exact for a mirror-symmetric plate, but a real multi-layer structure generally is not mirror-symmetric (different materials and/or thicknesses on the two faces — e.g. a thin metal transducer coating on one side only, common in laser-ultrasonics for optical generation, or any film-plus-substrate free plate). Without a mirror plane, the true elastodynamic stiffness/mass operators do not block-diagonalise by parity at all: every real mode is an inherent mixture of what would be "symmetric-like" and "antisymmetric-like" behaviour in the homogeneous limit. Labelling such a mode by its largest parity fraction then names only its dominant character, of *locally varying* trustworthiness, not a rigorous classification — and how much it matters is not a fixed property of the stack geometry alone; it also depends on kd, since mixing can be negligible at low kd and grow near a mode crossing at higher kd for the very same stack.

### The classification-fallback design

Rather than a single geometric threshold decided once from layer thicknesses (which cannot track the kd-dependence above), classification uses two tiers:

- **`Structure.isMirrorSymmetric(tol)`** — the exact condition for the mid-plane mirror to be a symmetry of the plate (`core.md`): for every layer pair `i`, `j = nLay+1−i` (the middle layer of an odd stack paired with itself), `relDepth` and density match and layer `j`'s standard-frame tensor is the *mirror image* of layer `i`'s, `C_j = M C_i M` with `M = diag(1,1,−1)` (relative `tol`). Pairing a layer with its mirror image rather than an identical copy covers both requirements in one test — the stack must be symmetric *and* each crystal must look the same from either face: a single Cu(111) layer fails (3-fold axis along the normal, C14-type entries), (001)/(110)/isotropic layers pass, and so does a bilayer of two (111) layers rotated 180° about the normal (mirror twins). For elasticity a mirror normal to the surface is equivalent to a 2-fold axis along the normal, since the tensor is unchanged by inversion. When true, every mode's **parity** is exact. Known conservative case, not worth handling while plates are not the focus: a plate symmetric only after merging identical adjacent layers (one material split 0.1/0.3) returns false, and its labels then follow confidence — safe, since only a false *true* can mislabel.
- **`Structure.isSagittalMirrorSymmetric(angleDeg, tol)`** — every layer's standard-frame tensor at that propagation angle is unchanged by `x2 → −x2` (`M = diag(1,−1,1)`). Then SH (u2) decouples from the sagittal (u1, u3) motion and the **Lamb/TH** half of the label is exact. Angle-dependent: Cu(001) passes along ⟨100⟩ and ⟨110⟩, not at 30°. With both mirrors the tensor is orthotropic in the propagation frame and all four labels are rigorous.
- **`CalcWavevecAngSweep2Plate.classifyBranches(angleIdx, confidenceThresh)`** — trusts the label when both mirrors hold, and otherwise when `modeTypeConfidence` reaches the threshold, per (wavevector step, tracked branch), evaluated *locally* (a branch can be confidently labelled away from a crossing and ambiguous right at one — a single global "is this stack asymmetric enough" decision cannot capture that). Combines `DispersionBranches.track` (which accepts `CalcWavevecAngSweep2Plate`) with the per-mode confidence to produce, for every (step, branch), the physical label (`'Asymm'` etc.); below the threshold, the exact parity class (`'Asymm-class'`/`'Symm-class'`) where the mid-plane mirror holds, or a stable `'branch N'` fallback string where it does not — feed the latter directly to a qualitative (`lines`-style) colormap keyed by branch number, since branch identity itself is exact (chained from step to step by the modal assurance criterion, MAC) even when the physical label is not.

`confidenceThresh` defaults to 0.3: the top parity fraction must exceed the runner-up by 0.3 (e.g. 65 % vs 35 %). **This default is unvalidated** — not tuned against real multi-layer data or literature, because the threshold needs measured or independently classified multi-layer mode types, which the GEWtool comparison of the next section (frequencies, group velocities, ZGV points, parity classes of a mirror-symmetric stack) does not provide. Treat it as a reasonable starting point, not a derived constant.

## Multi-layer stacks: validated against GEWtool

The null-space/constraint-matrix projection `RitzSolverLamb.calc` uses for `nLay>1` is the identical pattern `docs/reference/architecture/ritz.md`'s "Multi-layer stacks" section describes for SAW — inter-layer displacement-compatibility rows on top of the boundary condition (free surfaces here, rigid bottom there), generic over layer count. Nothing in the eigensolve path restricts a Plate structure to a single layer or to a mirror-symmetric one.

**Single-layer validation:**

- **Isotropic:** A0 and S0 of an Al plate at several `relDepth` values agree with the classical Rayleigh–Lamb frequency equations to ~1e-9 relative or better.
- **Anisotropic:** the S1 ZGV points of a [001]-cut Si wafer agree with published values computed by two other codes and with laser-ultrasonic measurements (see "Zero-group-velocity (ZGV) points" below).

**Multi-layer validation.** No analytic reference exists for a multi-layer anisotropic plate, so GuliweRR was compared with an independent numerical code: GEWtool v2.1.1, open-source software by D. A. Kiefer (ISC licence, doi:10.5281/zenodo.10114243; spectral elements on Gauss–Lobatto–Legendre nodes, one element per layer, discretised differently from the Legendre Ritz basis used here). The two codes share only the material data. Four stacks, each 1 mm thick in total, propagation 25° off [100] unless stated, k from kd = 0.1 to 15 (120 wavevectors), lowest 25 modes. Materials: Cu cubic, C11/C12/C44 = 170/122/74 GPa, ρ = 8.96 g/cm³; Al isotropic, C11/C44 = 108/26.4 GPa, ρ = 2.70 g/cm³; the strongly anisotropic layer an illustrative cubic crystal, 129/109/82 GPa, ρ = 8.0 g/cm³, at Bunge angles (17°, 41°, 63°).

| Stack (bottom to top, mm) | Why |
|---|---|
| Cu(001) 0.9 + Al 0.1, 25° | coated plate, coupled SH–Lamb |
| Al 0.8 + Cu(001) 0.2, 25° | the same materials in the other order and thickness ratio |
| Cu(001) 0.3 + strongly anisotropic cubic crystal (Zener ratio 8.2) 0.45 at a general orientation + Al 0.25, 31° | three layers, no symmetry axis, no mirror plane |
| Al 0.2 / Cu(001) 0.6 / Al 0.2, 25° | mirror-symmetric sandwich, exact parity classes |

Discretisation is `legOrder` = 28 per layer in GuliweRR and 28 nodes per layer in GEWtool (both converged for 25 modes at kd = 15; 20 is not: 1e-4 in frequency at the highest mode; GEWtool's rotated tensors reproduce GuliweRR's standard-frame tensors to 4e-10 GPa). Results:

- **Frequencies**, matched by value (both lists sorted at each k): relative difference at most 3e-10 to 9e-10 for kd > 1, over all four stacks and all 25 modes. At kd ≤ 1 it grows as the frequency falls — below 1e-9 above 0.2 MHz, about 1e-8 between 0.05 and 0.2 MHz, up to 1e-5 below 10 kHz. That is round-off in both codes (each code's own change from one `legOrder`/node count to another is of the same size there), not a discrepancy between them.
- **Group velocity** `velsGroup` against GEWtool's axial group velocity: at most 1e-8 km/s (relative 2e-8, against |vg| floored at 50 m/s) for kd > 1; GEWtool's axial energy velocity equals its group velocity to 1e-14 km/s, so the lossless identity holds in both.
- **ZGV points** (`locateZGVs` against `computeZGVScan`, run to 2.5 times the frequency band, f up to about 11 MHz and k up to 15 rad/mm): 11, 4, 14 (15 in GuliweRR, see below), and 5 points in the four stacks, and 10 for a single Cu(001) layer at 0°, all in both lists, k and f identical to 2e-9 relative or better.
- **Parity classes** of the sandwich: the symmetric/antisymmetric class of every mode (modes 1–30 at 120 wavevectors) equals GEWtool's exact symmetric and antisymmetric subspaces; 0 mismatches in 3572 mode labels (near-degenerate pairs excluded).

The check is code against code, not against a measured multi-layer plate, and it does not test the `confidenceThresh` default (above) or mixed-class labels in an asymmetric stack. It covers the free plate only. In the three-layer stack GuliweRR found one ZGV point near a cut-off (k = 1.1452 rad/mm, 11.3179 MHz) that GEWtool's global `computeZGVScan` did not; GEWtool's own Newton step started from it confirms the point. Where `locateZGVs` itself misses points on coarse wavevector grids, see "Known limitations" below.

## Zero-group-velocity (ZGV) points

A zero-group-velocity (ZGV) point is where a dispersion branch's group velocity (the speed at which wave *energy* travels) vanishes at a finite, non-zero wavevector — producing a sharp, long-lived resonance with strong elastic-constant sensitivity, and hence of direct interest for laser-ultrasonic/TGS-style measurement.

### Group velocity along k and across k

The group velocity of a plate mode is an in-plane vector, `c_g = ∇_k ω`. GuliweRR computes both of its components at every solve, without an extra eigen-solve:

- **Along k:** `velsGroup` = `∂ω/∂k` = `vp + Σᵢ dᵢ·∂vp/∂dᵢ` (`vp + D·dvp/dD` for a single layer; `ritz.md`, "Multi-layer stacks").
- **Across k:** `velsGroupTrans` = `c_gθ = (1/k)·∂ω/∂θ = ∂vp/∂θ` (θ the propagation angle in radians, k fixed), positive towards increasing angle. `RitzSolverLamb.calcTransGroupVelocity` evaluates it by the same Hellmann–Feynman (first-order perturbation) step as the thickness derivative: `∂ω²/∂θ = αᴴ (∂K/∂θ) α` for the solved eigenvector α. `K` is linear in the 15 rotated reduced constants, so `∂K/∂θ` is exact once `∂c/∂θ` is known; that is a ±0.01° central difference of the rotated constants alone (no solve), as in `RitzSolver.angleSlope`, which it equals to round-off.

From the two, the energy (group) velocity vector is `c_g = velsGroup·k̂ + velsGroupTrans·θ̂`, and the skew angle between `c_g` and `k` is `α = atan2(velsGroupTrans, velsGroup)`. `velsGroupTrans` vanishes on a direction with a sagittal mirror (e.g. ⟨100⟩ and ⟨110⟩ of a cubic (001) plate) and for an isotropic plate at every angle. Checked against a ±0.001° finite difference of whole solves (~1e-7 km/s, the O(h²) of the constants' difference; one layer and a Cu/Al bilayer) and against GEWtool's energy-velocity vector `energyVelVec` (≤ 1.2e-7 km/s over the eight lowest modes of a Cu(001) plate and of a non-mirror-symmetric Cu/Al bilayer at 25°, same sign). It is computed on every plate solve, at about 20–26 % of the solve's cost (one layer and a bilayer at legOrder 20). Stored as `CalcAngSweepPlate.velsGroupTransMap` and `CalcWavevecAngSweep2Plate.velsGroupTransMap3D`.

### Locating points of zero group velocity (`locateZGVs`)

`CalcWavevecAngSweep2Plate.locateZGVs(angleIdx, maxCount)` works at one fixed angle, in three steps:

1. **Candidate scan** — pure post-processing on a `DispersionBranches.track`-tracked branch's `velsGroupMap3D`, no extra forward solves: (a) a sign change of the group velocity between two wavevector steps; (b) a *dip*: a step where |vg| is below half its value at both neighbours and below 1 % of vp (`ZGV_DIP_REL_TOL`) without a sign change next to it. Two ZGV points of one branch between two grid steps (a minimum and a maximum of ω(k) close together) leave vg with the same sign at both steps; the dip catches that case.
2. **Newton refinement** — from the seed (k₀, f₀), the linear interpolation of the sign change along the branch (for a dip: the dip step and its two neighbours, one candidate each), `ZgvNewton.solve` solves the plate's eigenproblem written as a quadratic polynomial in k, `(k² L2 + k L1 + L0) u = ω² M u` (`RitzSolverLamb.kPolynomial`, the same matrices `calc` solves, to round-off, built once per call), together with the zero-group-velocity condition `uᴴ(2k L2 + L1)u = 0` and `uᴴu = 1`, by Newton's method on `[u; k; ω²]` (bordered Jacobian; a complex correction keeps k and ω² real). It converges in 3–6 steps, each one linear solve of the size of the Ritz basis, and needs **no mode number** — so a branch crossing inside the bracket does not matter. When it does not converge from the interpolated seed, or converges outside the bracket widened by one step (onto a neighbouring point), the two bracket ends are tried, each starting from the sweep's own eigenvector of the branch there: where two modes nearly coincide in frequency the interpolated seed can start on the neighbouring mode (the three-layer stack of the GEWtool comparison above: 3 of 82 points were lost without the retries). The method is GEWtool's `ZGVNewtonBeta` (D. A. Kiefer; algorithm by B. Plestenjak; Kiefer, Plestenjak, Gravenkamp & Prada, *J. Acoust. Soc. Am.* 153, 1386 (2023)), adapted to GuliweRR's operator and with guards: a non-Hermitian operator is refused (`ZgvNewton:solve:nonHermitian`), and an ill-conditioned or non-finite step, a singular Jacobian at the solution (a ZGV that is not isolated, e.g. a band flat in k) or a run into k → 0 (a cut-off) count as not converged.
3. **Acceptance** — a Newton point is kept if it lies inside the swept wavevector range and an independent `RitzSolverLamb` solve there has a mode at that frequency with `|vg| ≤ 1e-4·vp` (`ZGV_VG_REL_TOL`; a genuine point sits at ~1e-12). Two candidates that converge onto the same point are reported once. `branchNo` is the candidate's branch when the point lies inside its bracket, else the tracked branch nearest in frequency there.

Each point carries `velGroupTrans` and **`isTrueZGV`** = `|velGroupTrans| ≤ 1e-4·vp` (the acceptance bound above; on Si(001) |c_gθ|/vp grows by about 1e-3 per degree off [100] and 7e-4 per degree off [110], so a direction within about 0.1° of a symmetry axis still counts as a true ZGV point):

- **Symmetry directions** (a sagittal mirror, e.g. ⟨100⟩ and ⟨110⟩ of a cubic (001) plate): the perpendicular component vanishes, so the point is a true ZGV point, a minimum or a saddle point of the dispersion surface ω(k_x, k_y). `isTrueZGV = true`.
- **Other angles:** the perpendicular component generally does not vanish. The point then carries energy sideways and is not a ZGV resonance: Kiefer, Mezil & Prada (2023) call it a TGV (transverse-group-velocity) wave. It is returned, flagged `isTrueZGV = false`. They showed that a cubic (001) plate (Si) has exactly eight S1 ZGV points: four minima on the ⟨100⟩ directions and four saddle points on ⟨110⟩.

`RitzSolverLamb.locateZGV(structure, angle, thickness, modeNo, wavVec0)` runs the same Newton refinement from one raw mode at one wavevector (NaN when it does not converge).

**Checks of the refinement.** A ZGV point past a branch crossing inside a bracket is found (Cu(001) plate, 1 mm, legOrder 20, wavevectors 1–4 rad/mm: branch 12 at k = 3.9479 rad/mm). GEWtool's global `computeZGVScan` finds the same points (Cu(001) at 0°: 15 of 15 below 12 MHz; Si(001) at 0°, 45°, 18.4°, 22.5°; Cu/Al bilayer at 25°: 6 of 6) to ≤ 1e-10 in frequency and ≤ 8e-10 in k at matched discretisation (GuliweRR legOrder 30, GEWtool 30 nodes per layer). Isotropic Al plate (C11 = 108, C44 = 26.4 GPa, ρ = 2.70 g/cm³): the S1 ZGV frequency equals the Rayleigh–Lamb value, f·d = 2.856744100769 MHz·mm, to 1e-15, at legOrder 16 already.

**Known limitations.**

- *Grid resolution.* Candidates come from the sampled group velocity, so a ZGV point is found only if the wavevector grid shows it: a sign change, or a dip of |vg| below half its neighbours. Two points closer together than about one grid step, on branches that also approach each other (an avoided crossing of SH-like and Lamb-like modes), can still be missed: a Cu(001) 0.9 + Al 0.1 mm bilayer at 25° has two ZGV points at k = 1.612 and 1.746 rad/mm, missed on a grid of k = 1–4 rad/mm with 4 or 7 points and found with 13. Refine the grid where ZGV points are expected.
- *Along k at the swept angles only.* The search finds points where the axial component vanishes and flags whether the transverse one does too. A true ZGV point off the symmetry directions of a low-symmetry plate (a stationary point of ω(k_x, k_y) at an angle that is not scanned) would need a 2-D search in (k, θ); it is not done.

**Validation (anisotropic).** Tutorial S2, Part 5, takes a [001]-cut Si wafer: d = 524.6 µm, C11 = 165.6, C12 = 63.9, C44 = 79.5 GPa, ρ = 2330 kg/m³, the inputs given by Kiefer et al. GuliweRR finds the S1 points below; the published computed values are from GEWtool.

| Direction | GuliweRR | Published, computed | Measured |
|---|---|---|---|
| [100] (a minimum) | 7.7079 MHz at k = 3.4208 rad/mm, ZGV | 7.7079 MHz at k = 3.421 rad/mm | 7.7072 MHz |
| [110] (a saddle point) | 7.7551 MHz at k = 3.1423 rad/mm, ZGV | 7.7551 MHz at k = 3.142 rad/mm | 7.7545 MHz |
| 18.4° from [100] (26.6° from [110]) | 7.7265 MHz, TGV, c_gθ = 0.1778 mm/µs | 7.727 MHz, c_gθ = −0.177 mm/µs | ±0.170 mm/µs |

- GuliweRR gives the same values in MATLAB and Julia, and they are converged at Legendre order 16.
- Kiefer et al. measure the angle from [110] towards [010], so their direction is 71.6° from [100] in GuliweRR's angle; 18.4° is its mirror image across [110], an equivalent wave with the transverse component reversed. GuliweRR gives −0.178 mm/µs at 71.6°.
- Prada et al. (2009) computed the [100]/[110] points with DISPERSE as f·d = 4.044 and 4.069 MHz·mm, using ρ = 2329 kg/m³. GuliweRR gives 4.0445 and 4.0692 MHz·mm with that density.
- The A2 ZGV points near 15 MHz (Prada et al., their Fig. 9) are found as well, but the paper gives them only as a plot.

## GUI

Plate viewer family: see `docs/reference/architecture/gui_viewers.md`, "The Plate viewer family" (and `docs/reference/architecture/gui.md` for `StructureControlGUI`'s Free plate mode).

## References

Grabec, Sedlák, Seiner 2020, *Wave Motion* 96, 102567 (Ritz–Rayleigh method for Lamb waves in extremely anisotropic media) — the primary method reference, `docs/background.md`. Rose, *Ultrasonic Waves in Solid Media*, 1999 §6.2, and Achenbach, *Wave Propagation in Elastic Solids*, 1973 §6.4 — the classical Rayleigh–Lamb frequency equations that the single-layer isotropic validation compares against.

ZGV points in anisotropic plates, with the measurements of the anisotropic validation:

- Prada, Clorennec, Murray & Royer, *J. Acoust. Soc. Am.* 126, 620 (2009), doi:10.1121/1.3167277 — S1 and A2 ZGV resonances measured on a Si(001) wafer.
- Kiefer, Mezil & Prada, *Sci. Adv.* 9, eadk6846 (2023), doi:10.1126/sciadv.adk6846 — the eight true ZGV points, TGV waves, and the measured ZGV frequencies.

ZGV refinement: Kiefer, Plestenjak, Gravenkamp & Prada, *J. Acoust. Soc. Am.* 153(2), 1386–1398 (2023), doi:10.1121/10.0017252 — the Newton method `ZgvNewton` is adapted from, as implemented in GEWtool's `ZGVNewtonBeta.m` (D. A. Kiefer, ISC licence; algorithm design B. Plestenjak).

GEWtool v2.1.1, D. A. Kiefer, doi:10.5281/zenodo.10114243 (open source, ISC licence) — the external reference for the multi-layer validation, the ZGV points and the transverse group velocity.
