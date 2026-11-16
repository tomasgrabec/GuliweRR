# Christoffel — architecture and design rationale

`guliwerr/solvers/christoffel/` (`Christoffel` classdef + `@Christoffel/` static methods, plus the standalone `LimitingVels` classdef). Solves the Christoffel equation for bulk-wave velocities and polarisations in a homogeneous medium along a given propagation direction — the cheapest of the three solvers, and the source of both an upper bound on SAW velocity and the polarisation data used for mode identification.

The Christoffel equation is, by definition, the equation for a **plane harmonic wave** — the same plane-wave assumption Green's function and Ritz make from their own directions (`docs/background.md`, "The physics").

## What it computes, and why it's used alongside Ritz and Green

For a unit propagation direction **n** and density-normalised elastic tensor C, the Christoffel matrix Γ(n) = C_ijkl n_j n_l has three real eigenpairs (v², **g**): the squared phase velocities and polarisation (displacement) directions of the three bulk-wave branches (quasi-longitudinal L, quasi-shear T1/T2) that can propagate along **n** in an infinite medium. No boundary condition enters — this is the fastest of the three solvers and the only one that doesn't need a depth/layer structure at all.

It feeds the forward pipeline (`CalcSingleHS.calcLim`; `CalcSingleLayered.calcLim` uses the substrate only, as a soft reference) two things Ritz and Green don't directly give:
- **Limiting velocities** (`LimitingVels`) — for the bulk waves with wavevector in the sagittal plane whose energy flux runs parallel to the surface (no surface-normal group-velocity component), the phase velocity along the surface, v/cos θ*; the lowest one bounds the existence domain of true (non-leaky) SAW from above. `CalcSingleHS` uses `1.05 × longLimVel + 0.1 km/s` as the top of its self-built Green/velocity scan range for exactly this reason.
- **Polarisation vectors** at those limiting velocities, used by `ModeClassification` to score candidate Ritz modes by projecting the Ritz surface-displacement eigenvector onto the bulk-wave polarisation (`projection = abs(displNorm * polLim(:))`; for a degenerate group of limiting velocities the projection is onto the span of the group's polarisations, see "Degenerate roots" below) — a physical mode-type indicator independent of the Green spectrum. Applied to both the shear branches and the longitudinal branch — `associateLimitingMode`'s parameters are named generically `vLim`/`polLim`, since the function is genuinely called both ways.

## Calc-level methods: the Christoffel matrix itself

`matrix.m` is the hot path — called once per angle/candidate everywhere else in this folder:

```
Gamma(i,k) = C(i,j,k,l)*n(j)*n(l) / rho
```

implemented as `reshape(reshape(permute(elTensor,[1,3,2,4]),9,9) * kron(n,n), 3,3) / rho`, i.e. the (j,l) index pair contracted as a single 9-element dot product rather than 4 nested loops. `eig()` runs on the symmetrised matrix `(Γ + Γ.')/2` (Γ from a rotated tensor is symmetric only to round-off, and the non-symmetric solver can return an exactly degenerate pair — isotropic shear — as complex), and its output is sorted descending by eigenvalue, so index 1 is always L and 2/3 are T1/T2 (not physically labelled shear/quasi-shear beyond that ordering). `phaseVels` is a thin derivative (`sqrt(eigVals)`).

Unrelated to the Ritz `GammaMatrix` class (`guliwerr/solvers/ritz/@GammaMatrix`), which shares the name by coincidence.

## Scan-level method: Material-based angle sweep

`phaseVelScan` takes a `Material` object and a vector of rotation angles, extracts the propagation direction via `mater.GetRotDirMat(angleDeg)`, and loops `phaseVels` internally. The other methods take tensor + rho directly, so they're testable without constructing a `Material` (the input convention in `Christoffel.m`'s docblock).

## Group velocity and why the out-of-plane component matters

Group velocity: `vG_i = C(i,j,k,l)*g(j)*g(k)*s(l) / rho` (Royer–Dieulesaint Eq. 4.24), where **g** is the polarisation eigenvector and **s** = **n**/v_P is the slowness vector. `groupVelVectorsForDir` returns the full 3-vector for an arbitrary 3D direction, as full-sphere group-velocity ("ray") surfaces need (`Vis.plotGroupVelSurface`), since group velocity **V** = ∇_k ω(**k**) can leave the sagittal plane except on a mirror-symmetry plane of the material (Every et al. 2013, Wave Motion 50, §2.3). Its x3 (surface-normal) component at a sagittal direction is zero where a bulk wave's energy flux has no depth component, i.e. runs parallel to the surface. That's the physical condition `LimitingVels` searches for (through the cheaper `phaseVelAngleDeriv`, below).

Note: this is **not** the same quantity as `oopSurfDispl` in `RitzSolverSAW`/`ModeClassification` — that's a Ritz eigenvector's out-of-plane surface displacement, one of the per-mode diagnostics the mode classification reads (`RitzSolverSAW.modeDiagnostics`). The names are similar ("out-of-plane") but they belong to different solvers and measure different things (bulk-wave energy-flux direction here vs. a Ritz eigenvector's surface-displacement component there); do not conflate them.

## LimitingVels: root-finding algorithm

`LimitingVels.findRoots` (called via `CalcCommon.limitingAndBulk` by `CalcSingleHS`/`CalcSingleLayered.calcLim`) locates directions in the sagittal half-plane θ ∈ (−90°, 90°) where the out-of-plane group-velocity component of a given wave type crosses or grazes zero.

**Coarse scan, then per-candidate refinement:**
1. Evaluate `vg2(θ) = v(θ)·sinθ + v'(θ)·cosθ` (`vg2` is the out-of-plane, x3, component of the group velocity and `vg1` the in-plane, x1, one) on a uniform grid `θ = -π/2+1e-4 : alphaPrecision : π/2-1e-4` (`findRoots` default `alphaPrecision = 1e-2` rad, ≈ 300 samples over the half-plane). `v(θ)` and `v'(θ)` come from `Christoffel.phaseVelAngleDeriv` — a Hellmann–Feynman sandwich `dv/dθ = g' · (dΓ/dθ) · g / (2v)` that reuses the eigenvector already computed for `v(θ)`, cheaper than a full group-velocity contraction (`groupVelVectorsForDir`) per sample, and equal to it (`vg1 = v cosθ − v' sinθ`, `vg2 = v sinθ + v' cosθ`). It takes its angle in **radians**, unlike `phaseVelScan`'s degrees, since the derivative formulas are radian-native and a per-call-site deg/rad conversion is exactly the kind of silent factor-of-(180/π) bug this asymmetry is designed to avoid.
2. `islocalmin(abs(vg2))` flags candidate points — this catches both genuine sign-crossing roots and tangential touches (a bulk wave that grazes zero out-of-plane velocity without crossing sign — a real surface-skimming feature, not scan noise).
3. Each candidate is refined with `fminbnd` inside its neighbouring grid cell, minimising `|vg2(θ)|` to sub-grid precision, with `TolX = 1e-12` rad rather than `fminbnd`'s default 1e-4 rad. The precision actually reached is about `sqrt(eps)·|θ|` ≈ 1e-8 rad (the relative term in `fminbnd`'s stopping rule dominates, and a `TolX` below ~1e-9 changes nothing), which is still four orders better than the default. At a steep sign crossing the residual left by a θ error is slope × error; with the default tolerance a genuine root of a strongly anisotropic cubic (slope of tens of km/s/rad) would fail the accept test, taking the slowest shear limit and with it the Rayleigh search. The θ error would also shift the accepted velocity, since `vg1` is not stationary at the root. `refineRoot` (below) runs the same refinement and accept test (one private helper).
4. A refined candidate is accepted only if its residual is small relative to the in-plane group-velocity magnitude at that point (`residual < relTol * abs(vg1Root)`, default `relTol = 1e-3`) **and** the in-plane component is positive (outward-propagating). Only then is `Christoffel.matrix` called once more, at the refined angle, to extract the polarisation eigenvector.

**Why coarse-scan-then-refine.** Testing the grid samples directly against a fixed *absolute* threshold would make the survival of a genuine root — or of a tangential touch — depend on whether the grid happens to land a sample inside that window. That is a matter of luck, not a resolution guarantee, and it can drop real roots, including the longitudinal one (which is always exactly one), for generic non-symmetric orientations. Refining every candidate and judging it against a threshold relative to the local group velocity decouples correctness from grid alignment.

The longitudinal branch (`iwt=1`) physically yields exactly one root; `CalcCommon.limitingAndBulk` still applies `unique()` there and takes the lowest velocity, to guard against numerical duplicates from the scan.

`LimitingVels.refineRoot(rotElTensor, rho, thetaGuess, iwt, halfWidth, relTol)` re-finds one already-located root inside a bracket of `±halfWidth` (default 0.05 rad) around `thetaGuess`, with the same accept test as `findRoots`; `ok = false` means the root has disappeared and the caller must fall back to a full `findRoots` scan. Because the caller fixes the bracket, the returned velocity is a smooth function of the elastic constants, which a finite-difference derivative needs; `findRoots` brackets from whichever grid cell `islocalmin` flags and can jump by one cell under a tiny perturbation. `findRoots` also returns each root's θ* (third output) to serve as that bracket centre. Used by `inverseCalcs/SensitivityAnalysis`.

**Degenerate roots.** Where roots of different wave types coincide — an acoustic axis lying in the sagittal plane, e.g. cubic (001) at 0°, or any direction of an isotropic medium — the Christoffel eigenvalue is multiple, so each root's own eigenvector is an arbitrary vector of a 2-D eigenspace, decided by round-off. `findRoots` therefore replaces the polarisations of such a group (velocities within `LimitingVels.DEGENERATE_REL_TOL` relative, θ* within `DEGENERATE_THETA_TOL`, and a phase-velocity eigenvalue that is multiple there) by a canonical orthonormal basis of the eigenspace (`canonicaliseDegenerate`, `canonicalBasis`), so MATLAB and Julia return the same vectors. Because that basis only spans the eigenspace, `ModeClassification` associates a degenerate group with Ritz modes jointly (`degenerateLimitGroups`, `associateLimitingGroup`), using the projection of the surface displacement onto the span of the group rather than onto a single vector.

## Coordinate convention

Standard literature convention throughout (x1 = propagation, x2 = in-plane, x3 = surface normal/depth). `Christoffel.matrix`'s own formula is coordinate-frame agnostic — Γ_ik = C_ijkl n_j n_l is self-consistent in whatever frame C happens to be expressed in, since it only involves a propagation-direction vector, not a surface. The caller side enforces the convention: `CalcCommon.limitingAndBulk` (which `CalcSingleHS`/`CalcSingleLayered.calcLim` call) rotate the material's elastic tensor into the standard frame via `Material.GetStandardElastTensor` before calling into Christoffel — the same one-call pattern used by Ritz and `GreenSolver` (`core.md`).

## References

Royer, D. & Dieulesaint, E., *Elastic Waves in Solids I: Free and Guided Propagation*, Springer (2000 English ed.) — Eq. 4.24 (group velocity). Every, A.G. et al., Wave Motion, 50, 1197–1217 (2013), doi:10.1016/j.wavemoti.2013.02.007, §2.3 (group velocity direction off the sagittal plane; `groupVelVectorsForDir`).
