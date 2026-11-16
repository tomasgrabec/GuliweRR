# GuliweRR — quick start

This page shows the shortest working calculations. Read `docs/reference/functionReference.md` and `docs/pitfalls.md` before writing your own scripts. For design notes on the internals, see `docs/reference/architecture/`.

## Path setup

Once per MATLAB session, from the project root:

```matlab
run('guliwerr/init.m')
```

To load it automatically every session, add the same line with a full path to your `startup.m` (see the README, "Installation"). **Never** use `addpath(genpath(...))` on `guliwerr/` — see `docs/pitfalls.md` OTH2.

Run it at the prompt, not from inside a script started with `run('somescript.m')`: `run` changes into that script's folder first, which breaks the relative path. This is why the example scripts do not call it themselves.

## Minimal working example — half-space SAW velocity

The shortest path from "a material" to "a SAW velocity": build a `Material`, wrap it in a `Structure`, run a one-angle `Calc.angSweepHS` sweep.

```matlab
run('guliwerr/init.m');

%% Material -- isotropic Al (isotropic class takes [c11, c44] only)
alConsts = [108, 26.4];   % [C11, C44] GPa
alSymm   = 1;              % isotropic
alRho    = 2.70;           % g/cm^3

ec  = ElastConsts(alConsts, alSymm);
mat = Material(ec, alRho);
cp  = CalcParams(20, 50);                % relDepth=20, legOrder=50 -- "easy" tier, docs/pitfalls.md MTH1

str = Structure(mat, cp, 'al_isotropic', 500);   % material, calcParams, name, nMaxModes

%% Single-angle sweep (angles can be a 1-element array)
sweep = Calc.angSweepHS(str, 0, 'ritzOnly', 'saw_only', []);
sweep = sweep.calc();

vSAW = sweep.sawVelMap{1}(1);   % km/s
fprintf('v_SAW = %.4f km/s\n', vSAW);
```

For a standard crystal, take the constants from the materials library instead of typing them in (`Material.ListLibrary()` shows what is there, with crystal structure and sources: 48 common metals, semiconductors and oxides at room temperature):

```matlab
matCu = Material.FromLibrary('Cu', [0 0 1; 1 1 0]);   % Cu(001), propagation reference [110]
matMgO = Material.FromLibrary('magnesia', [0 0 1; 1 0 0]);  % name or alias, any case; 'T'/'variant' pick among several rows
```

`Calc.angSweepHS` (not the lower-level `Calc.singleHS`) is the right entry point even for one angle. `Calc.singleHS` builds its own Ritz precompute when none is passed in, which costs more than the solve itself when repeated per angle; the sweep builds it once and shares it.

## Full forward example — angular sweep + visualisation

`calcMode='full'` also runs the Green's-function solver. The material is now anisotropic Cu: an
isotropic material's SAW velocity does not depend on angle, so its angular map would be flat.

```matlab
cuConsts = [170, 122, 74];   % [C11, C12, C44] GPa, cubic
ec  = ElastConsts(cuConsts, 2);   % symmNo=2: cubic
mat = Material(ec, 8.96);         % rho, g/cm^3; default dirMat = Cu(001)
cp  = CalcParams(20, 50);         % "easy" tier (Zener ratio 2*C44/(C11-C12) = 3.1, principal (001) plane), docs/pitfalls.md MTH1
str = Structure(mat, cp, 'cu_001', 500);

angles = (0:2:90)';
sweep  = Calc.angSweepHS(str, angles, 'full', 'full');   % 'full' calcMode runs Christoffel+Green too
sweep  = sweep.calc();

Vis.plotGreenMap(sweep);                 % |G13_eff| map (what the classification read); 6th arg 'G33' for |G33| (recomputed)
Vis.plotClassificationDiagnostic(sweep, 'type', true, 'full');   % classified points over the Green background; refLevel: 'sawOnly' | 'sawLimit' | 'full'
```

### Choosing the classification

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface waves in strongly anisotropic materials. Thresholds, scores and the interface may change.

Which peaks count as surface waves is decided by a score against a standard. Three methods compute
the score:

- `'green'` (default, `calcMode='full'`): `GreenScore` on the Green's function \|G13_eff\|, the map
  plotted above;
- `'refined'` (half-space `'full'` only): `GreenSpotScore`, the Green's function averaged over the
  probe spot's directions — closer to a measured map; for a 1° sweep it costs about
  half the time of the `'full'` calculation itself;
- `'ritz'`: `RitzScore` on the Ritz modal sum — the only method of `calcMode='ritzOnly'` (chosen
  there automatically), and free on a `'full'` result.

Each has two named standards: `'high'` is sensitive (every mode a good map from transient grating spectroscopy, TGS, shows; see `docs/background.md`), `'low'`
keeps only the distinct modes (poor TGS data, or spatially resolved acoustic spectroscopy, SRAS). The default is `'high'` for `'green'` and
`'refined'`, `'low'` for `'ritz'`. Choose before `calc()`, or switch afterwards without
recomputing anything (only `'refined'` solves extra spectra, once, on first use):

```matlab
sweep = Calc.angSweepHS(str, angles, 'full', 'full');
sweep.sawClassMethod = 'green';    % '' (default: 'green' on 'full', 'ritz' on 'ritzOnly') | 'green' | 'refined' | 'ritz'
sweep.sawThreshold   = 'low';      % 'high' | 'low' | [] (= the method's default) | a number (needs sawClassMethod) -- set before calc()
sweep.greenVelStep   = 0.004;      % optional: finer Green map (km/s; default 0.006) -- labels unchanged
sweep = sweep.calc();
sweep = sweep.reclassify('high');            % later: another standard, same method
sweep = sweep.classify('refined');           % another method, at its default standard
sweep = sweep.classify('ritz', 'low');       % method and standard together
T = sweep.refinedSawAgreement();             % 'green' vs 'refined' entries, paired per angle
```

`sweep.activeMethod` / `sweep.activeThreshold` say which classification the maps and
`significantModesMap` hold. Each entry carries its `method`, `score`, and a `grade` (`'strong'`
passes `'low'`; `'weak'` only the looser active standard); the fields are listed in
`docs/reference/functionReference.md`, "Mode classification". Prefer the names `'high'` / `'low'` to a number
as a standard: a number is read on the method's own scale (`docs/pitfalls.md` MTH6). A single
angle classifies exactly as the same angle inside any sweep.

On a layered stack (`Calc.singleLayered`, `Calc.angSweepLayered`, `Calc.wavevecAngSweep2Layered`)
the methods are `'green'` (`GreenScoreLayered`) and `'ritz'`, both with default `'low'`: Rayleigh-
and Sezawa-type waves. Love waves need `sweep.reclassify('high')`; leaky waves need `'green'`
(`docs/pitfalls.md` PHY6 and MTH9).

A result whose classification stamp is older than the current rules (`sweep.isStale()`) still
loads and plots, but cannot be reclassified: recompute it (`docs/pitfalls.md` MTH13). Read
`docs/pitfalls.md` PHY1, MTH5–MTH7, MTH9, and MTH12 before reading much into a single entry.

For runnable versions of this pattern (and the layered and isoFreq equivalents), copy from `examples/tutorials/` rather than retyping.

## Where to go next

- **`docs/guiGuide.md`** — the click-driven alternative to everything above: build a `Structure` and run a calculation through `Gui.structureControl()` instead of a script, then browse the result in a viewer window. Good for exploring one material or checking a single point; script it instead (see the examples above and `examples/tutorials/`) for anything you'll repeat.
- **`examples/tutorials/`** — narrated, copy-paste-ready scripts, one per calculation type (isotropic half-space, angular map, modal shapes, 3D surfaces, layered Green, film on substrate, dispersion/double sweeps, isoFreq, inversion, material input `S1`, free-plate Lamb waves `S2`). Start here for a new script.
- **`examples/tutorialsCode/`** — the same tutorials, same file names, stripped to code: one-line section headers and comments, every argument a named variable so each call's argument order can be copied directly. Use it when you know the physics and just want the calls.
- **Julia version** — `juliaPort/` holds the same calculations as the Julia package `GuliwerrJ` (the same results within numerical tolerance, faster; MATLAB stays the reference version, and its interactive viewers are functional while Julia's are experimental — Julia's plotting functions work fine). Start with `juliaPort/docsJ/gettingStartedJ.md` and `juliaPort/docsJ/quickStartJ.md`; the code-only tutorials are in `juliaPort/tutorialsCodeJ/`, with the same file names as `examples/tutorialsCode/`.
- **Recovering elastic constants from measured velocities?** That's the *inverse* problem, not this file's forward one. It is experimental (the inverse module is validated on synthetic data only, the interface may change, and its pseudo-surface-wave tracking is new). See `examples/tutorials/tutorial_09_inversion_synthetic.m` and `docs/reference/functionReference.md`'s "Inverse characterisation" section (`Inv` namespace).
- **`examples/caseStudies/`** — real research scripts for specific published/measured materials (verbose, with literature provenance) — read these for *why* particular parameters were chosen for a real case, not as templates.
- **`docs/reference/functionReference.md`** — one-line-per-function index of every entry point, grouped by category.
- **`docs/pitfalls.md`** — gotchas that have each caused a real bug once.
- **`docs/background.md`** — what GuliweRR is and the physics behind the Ritz/Christoffel/Green pipeline.
- Any function's own `help <name>` — complete INPUT/OUTPUT documentation, denser than this file, with design-rationale prose included.
