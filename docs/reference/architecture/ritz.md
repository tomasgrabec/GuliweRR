# Ritz — architecture and design rationale

`guliwerr/solvers/ritz/` (`RitzSolver` base, `RitzSolverSAW` subclass, `@GammaMatrix/`). The primary forward solver: computes SAW mode velocities and eigenvectors for a half-space or layered anisotropic stack via the Ritz–Rayleigh variational method. `RitzSolverLamb` (plate/Lamb waves) is a sibling subclass sharing the same base and `GammaMatrix.calcLamb` path — see `docs/reference/architecture/plate.md` for everything specific to it (mode-type classification, ZGV).

Ritz reaches the same **plane-wave** assumption Christoffel and Green's function make by definition (`docs/background.md`, "The physics") from a different direction: the field is homogeneous in the in-plane transverse direction (x2) — see "The 21-to-15 elastic-constant reduction" below.

## The Ritz–Rayleigh method

The depth-dependent displacement field is expanded in a **single high-order Legendre polynomial basis** spanning the whole depth domain — one element, not many small low-order elements as in FEM (Stoklasová 2015; Grabec 2024). Each of the three displacement components (u1, u2, u3) gets its own `legOrder`-length Legendre expansion per layer, so one layer's basis vector length is `basisVecL = 3*legOrder`.

Energy minimisation via Hamilton's/Lagrange's principle turns the elastic-wave PDE into a generalised eigenvalue problem:

**K c = ω² M c**

— K the stiffness matrix (elastic constants contracted with Legendre-derivative integrals), M the density-normalised mass matrix, c the Legendre coefficient vector, ω² the eigenvalue (not velocity directly). In the normalised convention used throughout (`|k∥|=1`, relative length λ=1), phase velocity comes out numerically equal to frequency: `vels = freqs = sqrt(eigVals)/(2π)` (`RitzSolver.solveEigenvalueProblem`). λ=1 is hard-coded in the solvers, not a `CalcParams` property.

Rigid-bottom (SAW) or free-surface (Lamb) boundary conditions, plus inter-layer displacement continuity for multi-layer stacks, are encoded as linear constraints on c. `RitzSolver.buildConstraintMatrix` returns their null space N, and the unconstrained problem is projected: `K_con = N.'*K*N`, `M_con = N.'*M*N`, solved for the reduced coefficient vector α (`c = N*α`).

## Class hierarchy

```
RitzSolver          base value class — solveEigenvalueProblem, calcGroupVelocity,
│                    calcDisplacement, angleSlope, pruneEigVecs; static
│                    solveConstrainedEig, realGauge, buildConstraintMatrix
├── RitzSolverSAW    half-space and layered SAW: calc, modeDiagnostics; static
│                    ritzSensitivities (below)
└── RitzSolverLamb   plate — see docs/reference/architecture/plate.md
ZgvNewton           Newton refinement of plate ZGV points (plate.md)
```

`RitzSolverSAW.calc(withDiagnostics, precompStruct)` builds a `GammaMatrix`, projects K/M into the constrained basis, solves the eigenvalue problem, computes group velocity, and (if requested) the per-mode diagnostics (`modeDiagnostics`, below). `precompStruct` lets an orchestrator (e.g. an angular sweep) supply pre-computed `nullConstrMat` and `dKdC_cell`/`d2KdCdD_cell` — optional, computed on demand when absent (the two cells together or not at all).

The `GammaMatrix` is a local of `calc`, not a property of the solver: it is passed to `projectToConstrained` and `calcGroupVelocity` and dropped when the solve returns. Nothing reads it afterwards, and kept per angle it would be most of a sweep's memory — for a layered stack, several times the size of the *full* eigenvector set. Anyone who needs K, dK/dD, or M after the solve rebuilds them with `GammaMatrix().calc(rs.structure, rs.angle)` (`calcLamb` for a plate), which takes about a millisecond.

## The eigen-solve

`RitzSolver.solveConstrainedEig(K_con, M_con)` is the only Ritz eigen-solve: `solveEigenvalueProblem` and `ElasticityInversion`'s fit loop both call it. Two things keep it on MATLAB's fast paths:

- **Exactly Hermitian input.** `eig` takes the Hermitian solver (`eig(K)`) or Cholesky (`eig(K,M)`) only if its input is Hermitian *bit for bit*; `N.'*K*N` misses that by round-off and would fall back to complex Schur / QZ (up to 12× slower on stacks, eigenvectors only approximately M-orthogonal). So the helper symmetrises, `(A+A')/2`, first. The near-identity `M_con` of a single-layer half-space is still solved as `eig(K)`.
- **Real gauge.** If every layer's elastic tensor is invariant under a 2-fold rotation about the surface normal (isotropic, cubic (001)/(110), hexagonal, tetragonal/orthorhombic with a principal axis normal, ...), the imaginary part of K sits only in its in-plane/u3 blocks, and `D = diag(1, 1, i)` per Legendre coefficient makes `DᴴKD` real symmetric with M and N unchanged. `RitzSolver.realGauge` tests this numerically (`‖Im DᴴKD‖ ≤ 1e-13 ‖DᴴKD‖`, no symmetry bookkeeping; the tolerance must stay at round-off — computed with GuliweRR: cubic (111) has a residual of ~3e-5 that would shift eigenvalues by 9–19 % if forced real) and `projectToConstrained` projects the gauged K when it passes. The real symmetric solve is another 1.6–3× faster. `solveEigenvalueProblem` maps the eigenvectors back, `α = N.'(d .* (Nβ))`, so `eigVecs_con` is always in the original basis.

Each eigenvector is defined only up to a phase factor. `calcDisplacement` fixes it per mode (largest in-plane Legendre coefficient real and positive), so displacement profiles and the 2D fields built from them are reproducible; for the C2z structures above u1, u2 then come out real and signed and u3 as i × real. Everything else that reads `eigVecs_con` (MAC, Rayleigh quotients, `RitzScore`'s `u3·conj(u1)`) is phase-independent.

## GammaMatrix — building K and M

A three-step static pipeline factors the K-matrix build into what depends only on the mesh (`CalcParams`), what depends only on the Structure (layer geometry + material), and what changes at every angle:

| Step | Method | Depends on | Cached across |
|---|---|---|---|
| 1 | `rawByDepth` | `legOrder` only | all angles, materials and wavelengths |
| 2 | `derivsByDepth` | Structure (relDepth, ρ per layer) | the angle loop |
| 3 (Lamb path) | `elastByDepth` | Structure + angle | the frequency-sweep loop |

`rawByDepth` builds a `[bVL × bVL × 15 × 3]` tensor per layer: the 15 index runs over the reduced elastic constants (see below), the 3 index over a depth-factor group — `K0` (no depth derivative, factor `(2π)²`, from the in-plane/propagation-direction term `C_{i1,k1}`), `K1` (one depth derivative, Hermitian, factor `±2π`, cross term `C_{i1,k3}`/`C_{i3,k1}`), `K2` (two depth derivatives, factor 1, the pure depth-gradient term `C_{i3,k3}`). `derivsByDepth` contracts these against `relDepth` (`K = K0 + K1/(relDepth/2) + K2/(relDepth/2)²`) and divides by ρ to produce `dK/dc` (feeds `calc`) and `d²K/(dc·dD)` (feeds group velocity / depth sensitivity). The instance method `calc` then contracts `dK/dc` with the layer's rotated reduced elastic constants (`C = Σᵢ (dK/dcᵢ)·cᵢ`) at the requested angle.

The rest of `@GammaMatrix/` is small helpers around these steps:

- `ContractDepth(K0, K1, K2, relDepth)` — the depth factors, `K = K0 + K1/(d/2) + K2/(d/2)²` and `dK/dD = −2·K1/d² − 8·K2/d³` with `d = relDepth`; shared by the SAW path (`derivsByDepth`) and the Lamb path (`calcLamb`).
- `LayerNormFac(structure, ilay)` — the cross-layer normalisation factor described under "Multi-layer stacks".
- `LayerRotMat(structure, ilay, angleDeg)` — the rotation from the layer's crystal frame to the solver frame; it applies the angle once (see "Rotation convention" in `forwardCalcs.md`).
- `RedElastCellFor(structure, angleDeg)` and `RawCellFor(structure)` — the default fill of the rotated reduced constants and of the step-1 tensors per layer; public, so that a sweep can build the step-1 cell once and hand it down.
- `layerBlockDerivs` and `dummyForLayer` — `dK/dc` of one layer as its own diagonal block (or zero-padded to the full stack), used by the inverse module's sensitivity matrices.

**Naming note:** `GammaMatrix` is unrelated to Christoffel's own acoustic tensor Γᵢₖ = Cᵢⱼₖₗnⱼnₗ/ρ (`guliwerr/solvers/christoffel/@Christoffel/matrix.m`) — the name is a coincidence (`christoffel.md`).

## The 21-to-15 elastic-constant reduction

The Ritz field is harmonic in the propagation direction (x1) and explicitly homogeneous in the in-plane transverse direction (x2) — the "line-source" setup described in `docs/background.md`. Since the field never varies along x2, the strain component e22 = ∂u2/∂x2 is identically zero, so every elastic constant that couples to it (any Voigt-index-2 term, C_I2/C_2J, 6 of the 21 independent triclinic constants) contributes nothing to the strain energy and can be dropped before the Ritz computation. `ElastConsts.GetIndOf2DReduct_15of21` returns this 15-element index set; `ElastConsts.ElastRotReductMat` builds the corresponding rotated-constant matrix. Both live in `guliwerr/core/@ElastConsts/`, not in `ritz/`, but `rawByDepth`/`derivsByDepth`/`GammaMatrix.calc` all consume this reduced 15-constant representation directly.

**Coordinate-convention dependency.** This index set is convention-dependent — it must drop whichever Voigt index corresponds to the *current* completion (no-variation) axis. Under the standard convention (x1=propagation, x2=in-plane/completion, x3=surface normal/depth — reached everywhere via `Material.GetStandardDirMat`), the correct exclusion is Voigt index 2. Getting this wrong is silent and severe: it produces a velocity error of tens of per cent for layered structures. The error is self-consistent within Ritz alone — nothing inside Ritz's own maths flags it — and shows up only in an external identity check (same-material two-layer against half-space). Any change to the coordinate convention must update this index set as part of the same change, not as an afterthought.

Separately, `RitzSolverSAW(structure, angle)` rotates the *unrotated* structure internally, so a caller that has already rotated it must pass `angle=0` — see `docs/reference/architecture/forwardCalcs.md`, "Rotation convention: pass 0 to solvers after construction" (the "double-rotation bug" is applying the propagation angle twice).

## Legendre order (`legO`) and `relDepth`

`legO` (`CalcParams.legOrder`) is the truncation order of the Legendre basis per direction per layer; `relDepth` is the domain depth in wavelengths. The K-matrix terms above scale as `1/(relDepth/2)` and `1/(relDepth/2)²` — a deeper domain (larger `relDepth`) compresses more physical depth into the same fixed `[-1,+1]` Legendre argument range, so the true displacement field (decaying roughly like `exp(-|k∥|·z)` in a half-space, or oscillating layer by layer in a film) varies more rapidly across that range and needs a correspondingly higher polynomial order to represent without aliasing. `docs/pitfalls.md` MTH1 gives the resulting operational rule (`legO ≥ 2×relDepth`, with tiers scaling further by anisotropy (Zener ratio `2·C44/(C11−C12)` for cubic materials) and by the cut) — an empirical floor, not a derived convergence bound, and only a starting point for material/orientation combinations far from those validated. Under-resolution fails silently: a ratio of `legO` to `relDepth` of 1 gives velocity errors of up to a few per cent (computed with GuliweRR), depending on material and cut, that look like real physics, e.g. a genuine anisotropy effect. `legO` is therefore not merely a resolution knob: it is what keeps the single-element Legendre expansion from aliasing against the compressed depth scale that `relDepth` implies.

## Multi-layer stacks

`GammaMatrix.calc` loops generically over `structure.nLay`, block-diagonalising each layer's K/M into `gammaMatAll`/`normMatAll` with a cross-layer normalisation factor `normFac = (relDepth_i·ρ_i)/(relDepth_1·ρ_1)` (layer 1's substrate sets the reference scale). `RitzSolver.buildConstraintMatrix` adds inter-layer displacement-compatibility rows (Legendre value at the shared interface must match) on top of the rigid-bottom SAW constraint. A same-material two-layer structure reduces to the single-layer case exactly (block-diagonal cancellation), a machine-precision self-consistency check.

**Layered `depthInsensitivity`** — meant to measure sensitivity to the *artificial/computational* domain boundary — is computed against the substrate's truncation depth only, not against a film layer's real physical thickness (which would wrongly treat real geometry as a numerical-truncation knob): `RitzSolverSAW.calc` calls `calcGroupVelocity(gm, true)`, which takes the substrate block's (layer 1) Rayleigh quotient alone. For `nLay==1` the substrate block is all of `dGammaDMatAll`, so the half-space path is unaffected.

**Group velocity in a stack.** In the normalised convention (|k|=1), the phase velocity depends on the geometry only through `k·d_i`, so `vg = dω/dk = vp + Σᵢ dᵢ·∂vp/∂dᵢ` — each layer's depth derivative weighted by *that layer's own* `relDepth`. `calcGroupVelocity` builds this from per-layer Rayleigh quotients `qᵢ = UᵢᴴKᵢ′Uᵢ` of the block-diagonal `dGammaDMatAll` (which holds `∂Kᵢ/∂dᵢ`) on the full-space eigenvectors `U = Nα`, and weights block *i* by `dᵢ` — the Hellmann–Feynman quotient without ever projecting `dK/dD` into the constrained basis (that product costs O(nCon³) to read O(nCon) numbers). The depth dependence of `normFac` itself needs no explicit term: its contributions cancel exactly in the weighted sum (they sum to `αᴴKα − ω²αᴴMα = 0`), but not layer by layer. For `nLay==1` this reduces to `vp + D·dvp/dD`. Scaling the summed derivative by the total depth instead (`vp + D_total·Σᵢ∂vp/∂dᵢ`) is wrong for `nLay>1`: it over-weights thin films by `D_total/d_film` and gives unphysical (often large negative) layered group velocities. The result agrees with a finite-difference uniform dilation of the stack (SAW and Lamb). `dVel_dDepth` (all layers +1 at once) is unchanged and, for `nLay>1`, is *not* by itself a group-velocity term.

**Angular slope (`RitzSolver.angleSlope`).** The same Hellmann–Feynman projection gives `dv/dθ` of
any mode without a second solve: `K` is linear in each layer's 15 rotated reduced constants, so
`dK/dθ = Σᵢ (dcᵢ/dθ)·dK/dcᵢ` exactly, with only `dcᵢ/dθ` taken as a ±0.01° central difference of the
constants (cheap: no solve), and `dv/dθ = αᴴ(dK_con/dθ)α/(8π²v)`. Exact for a non-degenerate mode;
agrees with a finite difference of two further solves. Used by the
`'ritz'` classification (both geometries) for the aperture dephasing of a branch (see
`forwardCalcs_modeClassification.md`). Needs the full eigenvector set (before retention pruning).

**Sensitivities (`RitzSolverSAW.ritzSensitivities`; new, interface may change).** The same Hellmann–Feynman step, applied to every parameter of a layered stack at once: one solve gives the log-derivatives of the lowest modes' phase velocities with respect to every layer's elastic constants, density and depth, the angle and `|k|`. Two sums hold to round-off and check it: Σ dln v/dln c = ½ and Σ dln v/dln ρ = −½ per mode (v² is homogeneous of degree 1 in the constants and −1 in the densities). Exact for non-degenerate modes. For a sweep over wavelengths, `CalcCommon.ritzWavelengthSweep` builds the wavelength-independent part of the precompute (the step-1 tensors and the constraint null space) once.

## Validation

Ritz surface-wave results are validated in two independent ways. They reproduce the published Ritz–Rayleigh results and the measured velocities those papers compared them with (Stoklasová et al. 2015 for anisotropic half-spaces; Grabec et al. 2024 for NiTi films on MgO). And they agree with the Green's-function solver (`green.md`), which finds the same surface waves by an unrelated route — as peaks of the surface response, not as eigenvalues of a variational problem — so an error in the Ritz assembly, constraints, or rotation of the constants would show up as a disagreement between the two.

## Mode retention

`RitzSolver.pruneEigVecs(keepCols)` is the sanctioned way to discard eigenvector columns after solving — `eigVecs_con` has `SetAccess=protected`, so an external orchestrator cannot trim it directly. This is the *mechanism*; the *policy* (which retention level to use) is not a Ritz-solver concept — `retentionLevel` (`'full'`/`'saw_only'`/`'none'`) lives at the orchestration layer (`Calc.angSweepHS`/`Calc.angSweepLayered`/`CalcSingleHS`, see `docs/pitfalls.md` NUM4), which calls `pruneEigVecs` after `calc()` returns based on what the sweep actually needs downstream. `vels`, `eigVals`, and all mode-recognition scores are untouched by pruning — only the (large, `[strBVL × nModes]`, `strBVL` being the summed basis-vector length of all layers) eigenvector matrix is trimmed.

## Per-mode diagnostics, the Ritz score, and the Ritz spectrum

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface waves in strongly anisotropic materials. Thresholds, scores and the interface may change.

`RitzSolverSAW` computes mode velocities and eigenvectors and, with `withDiagnostics = true`
(default), per-mode diagnostics. It does not classify: which modes are surface waves is decided in
`ModeClassification` (`docs/reference/architecture/forwardCalcs_modeClassification.md`).

### `modeDiagnostics` (`diagnostics`)

`RitzSolverSAW.modeDiagnostics` fills `diagnostics` for every computed mode:

- **`freqDisp`** = `exp(−10·((v − v_g)/v)²)`: 1 for a non-dispersive mode, → 0 for a dispersive
  one. A homogeneous half-space has no length scale, so a true SAW is non-dispersive; a mode of
  the finite-depth Ritz slab is set by the domain depth and is dispersive. Meaningful on a
  half-space (and a Lamb plate) only: in a layered stack the film sets a length scale and true
  guided modes are dispersive too.
- **`oopSurfDispl`** = |out-of-plane displacement| at the top surface.
- **`depthInsensitivity`** = `ritz.depthInsensitivity` minus its minimum over the modes: large for
  modes whose velocity barely changes with the artificial boundary depth (substrate-only `dK/dD`
  for a layered structure, see above), i.e. modes confined near the surface.

`ModeClassification` reads them for the identity of a Green peak (`identityWeight`:
`oopSurfDispl` × `freqDisp` on a half-space, without `freqDisp` on a layered stack — the tie-break
of `matchRitzToGreenPeak`, the ranking of `matchSpotPeak`, and every entry's `matchScore`) and for
the limiting-velocity association (the raw `depthInsensitivity`). The Ritz-mode
viewers offer `freqDisp` and `depthInsensitivity` as colour criteria, beside the Ritz score.

### The Ritz score and the Ritz spectrum (`RitzScore`)

The `'ritz'` classification (the only one on `calcMode = 'ritzOnly'`, and computed on `'full'`
results too) is built from the Ritz modes by `guliwerr/forwardCalcs/RitzScore.m`. Its core, the
**Ritz spectrum** (a damped modal sum with the `G13`-like numerator), its storage, the score, the
supporting-mode logic (`RitzScore.localMode`, `RitzScore.decompose`) and the
half-space and layered rules are all defined in `forwardCalcs_modeClassification.md`, "The modal
sum (the Ritz spectrum)" and sections 4 and 6; the plots are `Vis.plotRitzSpectrumCurve` /
`plotRitzSpectrumMap` (`visual.md`).

## References

Stoklasová 2015 (SAW half-space Ritz baseline); Grabec 2024 (layered media, film-on-substrate generalisation); full citations in `docs/background.md`, "Key references"; `docs/pitfalls.md` MTH1 and NUM4.
