# GuliweRR — background

> The physics GuliweRR calculates and the literature it rests on.

## What GuliweRR calculates

**GuliweRR** — *Guided Linear Waves by Ritz-Rayleigh* — calculates linear elastic-wave propagation in anisotropic and layered media. Given a material stack (elastic constants, density, layer thicknesses, crystal orientation) and a propagation direction, it returns wave velocities, mode identities and polarisations, and the surface response to a surface force.

| Domain | Waves | Depth geometry | Lower boundary | Solvers |
|---|---|---|---|---|
| Bulk, homogeneous | Quasi-longitudinal and two quasi-shear bulk waves | None | — | Christoffel |
| Half-space, bare | Rayleigh and pseudo-surface (leaky) waves | Semi-infinite | None (Green); truncated, rigid (Ritz) | Ritz, Green |
| Half-space, layered (film-on-substrate) | Rayleigh, Sezawa, pseudo-surface waves; dispersive | N layers on a semi-infinite substrate | None (Green); truncated, rigid (Ritz) | Ritz, Green |
| Plate, homogeneous or multi-layered | Lamb waves (symmetric, antisymmetric, shear-horizontal); ZGV points | Finite | Free (both surfaces) | Ritz |

It was written to support **laser-ultrasonics / transient grating spectroscopy (TGS)** measurements, where forward calculations predict what a given sample should show and an inverse procedure recovers the elastic constants from measured velocities.
However, its use is not limited to that.

## The plane-wave solution

All three solvers look for plane-wave solutions of the elastic equation of motion,

ρ üᵢ = ∂ⱼ(Cᵢⱼₖₗ ∂ₗ uₖ),

in the form uₖ = ûₖ(x₃) e^{i(kx₁ − ωt)}: harmonic along the propagation direction x₁, uniform along the in-plane transverse direction x₂, with phase velocity v = ω/k. This mimics a line source, or a TGS grating with one harmonic direction. The solvers differ only in how they find the depth profile û(x₃), which is why their results are directly comparable at the same propagation direction.

Axes follow the standard literature convention (Every, Maznev): x₁ = propagation, x₂ = in-plane transverse, x₃ = surface normal / depth.

Because the field never varies along x₂, the strain ε₂₂ is identically zero, and the six elastic constants that couple only through it drop out: a guided wave along a given direction depends on 15 of the 21 independent constants.

### Christoffel — bulk waves

The simplest case is a bulk wave in an unbounded crystal. With no boundaries, the profile is a plane wave too, and the problem reduces to a 3×3 eigenvalue problem,

det[Cᵢⱼₖₗ nⱼ nₗ − ρv² δᵢₖ] = 0,

whose roots are the three bulk velocities (one quasi-longitudinal, two quasi-shear) along the direction n, with their polarisations.

For guided waves, Christoffel supplies the **limiting velocities**: the bulk velocities whose group velocity lies along the surface, in the sagittal plane. The slowest shear limiting velocity bounds true Rayleigh waves from above. A surface wave faster than it cannot stay confined to the surface: it couples to a bulk wave and leaks energy into the depth — a pseudo-surface (leaky) wave, lying between the slowest shear and the longitudinal limits (Lim & Farnell 1969). Group velocity off the sagittal plane: Every 2013.

### Ritz — the modal spectrum

Guided waves need a bounded domain — a free surface, layer interfaces, or both — and the profile is no longer a single plane wave. The Ritz (Ritz–Rayleigh) method expands it in Legendre polynomials in each layer — one high-order basis spanning the whole layer, not many small low-order elements as in FEM — and makes the energy stationary (Hamilton's principle). This gives one generalised matrix eigenvalue problem,

K α = ω² M α,

for all modes of the domain at once: K from the elastic constants and Legendre-derivative integrals, M the mass matrix, α the Legendre coefficients. The eigenvalue is ω², and the phase velocity follows as v = ω/k. Each eigenvector gives the mode's depth profile and polarisation at negligible extra cost.

- **Half-space (bare or layered).** The semi-infinite substrate is truncated at a finite depth with a rigid bottom. Displacement continuity at interfaces is imposed; traction continuity and the traction-free surface come out of the energy principle only approximately. A finite slab cannot carry energy away into the bulk, so true surface waves come out cleanly, while pseudo-surface waves need a deeper domain and more careful mode characterisation (Stoklasová 2015; Grabec 2024).
- **Plate.** Both surfaces are free — a genuine physical boundary on an already finite domain — and the eigenproblem gives the Lamb-wave dispersion (Grabec 2020).

With the full modal set known, the surface response to a surface force can be rebuilt as a damped sum over the modes — the Ritz spectrum, an approximation to what Green's function computes directly.

The same Ritz framework underlies resonant ultrasound spectroscopy (RUS) of finite samples (Landa 2009; Sedlák 2014).

### Green's function — the surface response

Green's function builds the profile exactly from Christoffel's partial waves in each layer. For a fixed in-plane slowness, the Christoffel determinant is a sextic in the depth wavenumber k₃; in the substrate, the three physical partial waves are kept — those that decay into the depth, or radiate away from the surface (Farnell 1970). All boundary conditions are imposed exactly, and the result is the surface Green's function G_IJ(k∥, ω): the displacement in direction J at the free surface in response to a harmonic surface force in direction I (Every 1997, 1998; layered: Maznev 1999, Every 2002).

Its peaks sit at the zeros of the free-surface boundary determinant, i.e. at the free surface modes: a singularity for a true surface wave, a finite broadened peak for a pseudo-surface wave, and smaller features at the bulk-wave thresholds. The depth is handled analytically, with no truncation. Supersonic surface waves on cubic crystals: Every 2015.

**±k symmetrisation.** A force grating that drives the +k and −k waves together (a standing TGS grating, for example) excites the effective Green's function

G_eff ≡ (G(+k) + (−1)ⁿ G(−k)) / 2, where n is how many of I, J equal 1.

Elastodynamic reciprocity gives G_IJ(−k) = G_JI(+k) for any elastic tensor, cut, and layer stack, so one +k solve gives both sides; for TGS, G₁₃^eff = (G₁₃ − G₃₁)/2. The two directions share every pole, but for a pseudo-surface wave their residues differ, so v(θ) = v(θ + 180°) holds for Christoffel and Ritz but not for the Green's-function response.

### Dispersion in layered stacks

A homogeneous half-space has no length scale, so its surface-wave velocities do not depend on wavelength. A film of thickness h introduces one: velocities depend on h/λ, and higher modes (Sezawa waves) appear as the film grows thicker relative to the wavelength. A dispersive spectrum can be scanned in two ways: at fixed wavelength (fixed k, solving for ω — Ritz and Green), or at fixed frequency (fixed ω, scanning k — Green). The two coincide only for a non-dispersive medium.

### Plates and zero-group-velocity points

In a free plate both surfaces are traction-free and the modes are Lamb waves, all dispersive in frequency × thickness. In an isotropic plate, or along a symmetry direction of an anisotropic one, they split into symmetric, antisymmetric, and shear-horizontal families; otherwise the shear-horizontal modes couple with the Lamb modes, and a homogeneous plate keeps only the symmetric and antisymmetric families (Grabec 2020). On some branches the group velocity vanishes at non-zero wavenumber — the zero-group-velocity (ZGV) points, where energy stays under the source and sharp local resonances appear (Prada 2009). In an anisotropic plate, a true ZGV point needs the full in-plane group velocity to vanish; off the symmetry directions only its component along k vanishes (TGV waves) (Kiefer, Mezil & Prada 2023). ZGV points are refined with a Newton method (Kiefer, Plestenjak, Gravenkamp & Prada 2023).

Green's function for a finite free plate is not implemented. The Christoffel limiting velocities do not bound the plate modes: a free plate has no bulk to leak into, and each higher-order mode starts at a cut-off frequency with k = 0, so its phase velocity ω/k falls from infinity there (Prada 2009). Multi-layer plates have been checked against GEWtool (Kiefer), an independent spectral-element code.

### The inverse problem

Given measured velocities, the elastic constants are found by iterative least-squares minimisation of (measured − calculated)². A sensitivity analysis first transforms the elastic-coefficient space (PCA-like) to find the directions the data constrain well, and only those are optimised (Stoklasová 2015 for the half-space; Grabec 2024 for layered systems).

## Key references

### Ritz method

- **Grabec, T. et al., Ultrasonics, 138 (2024). doi:10.1016/j.ultras.2023.107211** — main, most up-to-date reference: layered media (film on substrate, generalised to N layers), the current complex-exponential-harmonic formulation, and the inverse-procedure discussion for layered systems. Which papers to cite: README, "Citation".
- **Stoklasová, P. et al., Ultrasonics, 56, 381–389 (2015). doi:10.1016/j.ultras.2014.09.004** — original Ritz paper, SAW half-space, real-valued (sin/cos) harmonic formulation — the difference from Grabec 2024's complex-exponential harmonics is marginal. Cited for attribution and its half-space inverse-procedure discussion.
- **Grabec, T., Sedlák, P. & Seiner, H., Wave Motion, 96, 102567 (2020). doi:10.1016/j.wavemoti.2020.102567** — plate/Lamb waves, same real-valued formulation as Stoklasová 2015. The method reference for Lamb waves.
- **Landa, M., Sedlák, P., Seiner, H. et al., Appl. Phys. A, 96, 557–567 (2009). doi:10.1007/s00339-008-5047-4** — modal resonant ultrasound spectroscopy (RUS); the Ritz method described in more depth here than in the SAW-specific papers above.
- **Sedlák, P., Seiner, H., Zídek, J. et al., Experimental Mechanics (2014). doi:10.1007/s11340-014-9862-6** — determination of all 21 triclinic elastic coefficients by RUS.

### Surface waves and Green's function

- **Farnell, G.W., Properties of elastic surface waves, in Physical Acoustics, Vol. VI, eds. W.P. Mason & R.N. Thurston (Academic Press, 1970), pp. 109–166. doi:10.1016/B978-0-12-395666-8.50017-8** — the partial-wave method for surface waves on anisotropic crystals: the sextic for the depth wavenumber, three partial waves, and the free-surface boundary determinant whose zeros are the surface waves. The Green's function inverts the same boundary matrix, so its peaks sit at these free modes.
- **Lim, T.C. & Farnell, G.W., J. Acoust. Soc. Am., 45(4), 845–851 (1969). doi:10.1121/1.1911556** — pseudo-surface (leaky) waves: one partial wave radiates into the bulk, and the wave has a finite attenuation.
- **Every, A.G., Kim, K.Y. & Maznev, A.A., J. Acoust. Soc. Am., 102(3), 1346–1355 (1997). doi:10.1121/1.420053** — half-space surface Green's function, primary method reference.
- **Every, A.G., Kim, K.Y. & Maznev, A.A., Ultrasonics, 36, 349–353 (1998). doi:10.1016/S0041-624X(97)00039-5** — companion half-space Green's-function paper.
- **Every, A.G., Meas. Sci. Technol., 13, R21 (2002). doi:10.1088/0957-0233/13/5/201** — near-surface/thin-film elastic properties; layered Green's function.
- **Maznev, A.A., Akthakul, A. & Nelson, K.A., J. Appl. Phys., 86(5), 2818–2824 (1999). doi:10.1063/1.371130** — surface acoustic modes in thin films on anisotropic substrates; layered (film-on-substrate) Green's function, G₃₃(k∥,ω).
- **Every, A.G., J. Acoust. Soc. Am., 138, 2937 (2015). doi:10.1121/1.4934557** — supersonic SAW on cubic-crystal surfaces.
- **Every, A.G. et al., Wave Motion, 50, 1197–1217 (2013). doi:10.1016/j.wavemoti.2013.02.007** — bulk/surface acoustic wave review; group velocity off the sagittal plane.

### Plates and ZGV points

- **Prada, C., Clorennec, D., Murray, T.W. & Royer, D., J. Acoust. Soc. Am., 126(2), 620–625 (2009). doi:10.1121/1.3167277** — S1 and A2 zero-group-velocity (ZGV) resonances measured on a [001]-cut Si wafer; their frequency depends on direction.
- **Kiefer, D.A., Mezil, S. & Prada, C., Sci. Adv., 9, eadk6846 (2023). doi:10.1126/sciadv.adk6846** — the same Si wafer: eight true ZGV points, four minima on ⟨100⟩ and four saddle points on ⟨110⟩; at other angles the group velocity along k vanishes while the perpendicular component does not (TGV waves). Measured ZGV frequencies used for the anisotropic plate validation (tutorial S2, Part 5).
- **Kiefer, D.A., Plestenjak, B., Gravenkamp, H. & Prada, C., J. Acoust. Soc. Am., 153(2), 1386–1398 (2023). doi:10.1121/10.0017252** — computing ZGV points by a Newton method (GEWtool's `ZGVNewtonBeta`, D. A. Kiefer, ISC licence; algorithm by B. Plestenjak); `locateZGVs` refines its candidates with an adaptation of it.
- **Kiefer, D.A., GEWtool (v2.1.1), doi:10.5281/zenodo.10114243** — open-source (ISC licence) spectral-element code for guided waves; the external reference for the multi-layer plate validation.

Implementation details of each solver: `docs/reference/architecture/` (`christoffel.md`, `ritz.md`, `green.md`, `plate.md`, `forwardCalcs.md`).
