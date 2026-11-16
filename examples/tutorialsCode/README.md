# Code-oriented tutorials

The same tutorials as [`../tutorials/`](../tutorials/), same file names, with most of the explanatory prose removed. Use these to copy calls. Read the narrated versions for *why* each parameter has the value it does.

Conventions used in every file:

- one-line section headers (`%% Section -- what it does`);
- inline comments on one line (two lines at most); valid option values are listed as `'a' | 'b' | 'c'`;
- every function argument is a named variable, so each call's argument order can be read and copied directly (`Calc.angSweepHS(structure, angles, calcMode, retentionLevel, nVelGreen, useParallel)`).

Run `run('guliwerr/init.m')` once per MATLAB session first, from the project root (see `docs/quickStart.md`). Sections that open a GUI window need an interactive MATLAB session, not `-batch` mode.

| File | Topic |
|---|---|
| `tutorial_01_isotropic_halfspace.m` | minimal SAW velocity, isotropic Al |
| `tutorial_02_angular_map.m` | Cu(001) angular sweep, Green map, classification, G13_eff |
| `tutorial_03_modal_shape.m` | Ritz modal shape, Green's response, Green spectrum, combined viewer |
| `tutorial_04_3d_surfaces.m` | bulk velocity / group-velocity / Young's surfaces and in-plane curves |
| `tutorial_05_layered_green.m` | film on substrate, layered `'full'` sweep |
| `tutorial_06_film_on_substrate.m` | SAW velocity vs. film thickness, branch tracking |
| `tutorial_07_dispersion_map_and_viewer.m` | angle sweep, wavelength sweep, double sweep, viewers |
| `tutorial_08_isofreq_scan.m` | isoFreq scans at one and several fixed frequencies |
| `tutorial_09_inversion_synthetic.m` | elastic-constant inversion from synthetic data |
| `tutorial_S1_material_from_literature.m` | Voigt input, conversions, orientation (Miller/Euler/rotation), averaging |
| `tutorial_S2_lamb_waves_plate.m` | free plate: Lamb modes, dispersion, mode types, ZGV, coated plate, Si(001) ZGV check against measurements |

Tutorial 09 uses the inverse module, which is experimental: it is validated on synthetic data only, and the interface may change.

The Julia equivalents are in `juliaPort/tutorialsCodeJ/`.
