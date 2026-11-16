# Changelog

For users of the public GuliweRR release. `Unreleased` collects changes not yet in a tagged
release. Versioning (MAJOR.MINOR.PATCH, DOIs, the **Breaking** heading) is described in the
[README](../README.md#versioning).

Format: [Keep a Changelog](https://keepachangelog.com/) style (Breaking / Added / Changed / Fixed / Removed).

## [6.0.0] – 2026-11-16

First public release of GuliweRR, together with its Julia version, `GuliwerrJ`.

### Added

- **Forward calculations** for anisotropic half-spaces, layered structures (one or more films on
  a substrate), and free plates: the Ritz–Rayleigh solver (surface and Lamb waves), the
  Green's-function solver (half-space and layered), and Christoffel bulk and limiting velocities.
- **Calculation types**: single direction, angular sweep, wavelength–angle and frequency–angle
  maps, and isofrequency scans; modal shapes; dispersion-branch tracking; zero-group-velocity
  points of free plates, with the group velocity perpendicular to k and a flag that tells a true
  ZGV point from a TGV wave (the Newton refinement is adapted from GEWtool's `ZGVNewtonBeta`,
  D. A. Kiefer, ISC licence; Kiefer, Plestenjak, Gravenkamp & Prada, JASA 153, 1386 (2023)).
- **Results**: plotting functions (`Vis`), interactive viewers for every calculation type (`Gui`,
  with axis ranges settable from a script), and saving and reloading of finished results
  (`Calc.saveResult`, `Gui.viewSaved`).
- **Sensitivities**: exact (Hellmann–Feynman) derivatives of the layered Ritz modes' velocities
  with respect to every layer's elastic constants, density and thickness, the angle and |k|
  (`RitzSolverSAW.ritzSensitivities`; new, interface may change).
- **Julia version** (`juliaPort/`): the same calculations, one to one, checked against MATLAB.
- **Experimental** — working and tested, but their results and interfaces may still change:
  - surface-wave classification (`'green'`, `'refined'`, `'ritz'`), which decides which spectral
    peaks count as surface waves. Its thresholds are checked against the measured maps
    published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no
    general agreement on this question for strongly anisotropic materials;
  - spot-averaged Green maps (`ProbeKernel`, `calcSpotGreenMap`);
  - the inverse module (`Inv`): elasticity from measured velocities, for one or several
    crystallographic cuts, with pseudo-surface-wave tracking.
