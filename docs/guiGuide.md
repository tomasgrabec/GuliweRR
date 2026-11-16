# GuliweRR — GUI walkthrough

An interactive, click-driven alternative to the scripting workflow in `docs/quickStart.md` — build a `Structure`, run one forward calculation, and browse the result, all without writing a script. Good for exploring a material or checking a single point; for parameter sweeps, batch runs, or anything you'll want to repeat, script it instead (`examples/tutorials/`) — the GUI runs exactly one calculation at a time.

Everything below assumes `run('guliwerr/init.m')` has already been run this session (`docs/quickStart.md` §"Path setup").

## Launching

```matlab
win = Gui.structureControl();       % blank start
% win = Gui.structureControl(str);  % or: every tab pre-filled from an existing Structure
win.show();
```

Builds one window: a tab per layer (Material/CalcParams input) on the left, a control panel that launches a forward calculation on the right. Like every window in this guide it is built hidden; `show()` puts it on screen. An existing `Structure` passed in fills one tab per layer (direction matrix, density, symmetry, Cij, relDepth or thickness, legOrder); without one, the window starts with a blank Substrate tab.

## Step 1 — Build a `Structure`

The **Substrate** tab is always present and can't be removed. Each tab has:

- **Direction Matrix** — surface normal + first (propagation-reference) direction, 2 rows × 3 columns. Non-orthogonal input is auto-corrected (status line notes when this happens).
- **density [g/cm³]**
- **symmetry** — a dropdown (isotropic, cubic, tetragonal I/II, trigonal I/II, hexagonal, orthorhombic, monoclinic, triclinic). Choosing one enables exactly the Cij grid cells that symmetry class needs (`ElastConsts.IndepVoigtPositions`) and greys out the rest — dependent entries are echoed back automatically as you type.
- **Elasticity: Cij (Voigt matrix) [GPa]** — upper-triangle 6×6 grid, only the enabled cells are editable.
- **relDepth** (Substrate) or **thickness [µm]** (film layers) + **legOrder** — see `docs/pitfalls.md` MTH1 for the `legOrder ≥ 2×relDepth` rule; film-layer thickness converts to `relDepth` automatically at Run time via `Structure.ApplyWavelength` once a wavelength is known.

**[+ Add layer]** appends a new film tab on top of the stack (capped at 5 film layers — build more programmatically via `Structure.AddLayer` if you need it). **[- Remove active layer]** removes whichever tab is currently open (never the Substrate). **[Build Structure]** assigns the live `Structure` to a base-workspace variable (default name `guliStruct`, editable in the field beside the button) — useful if you want to keep building via script from this point, but not required before Run: the structure is already live internally.

## Step 2 — Configure the calculation

The right-hand control panel:

| Control | Choices | Meaning |
|---|---|---|
| **Calculation target** | Single layer (HS) \| Full stack (Layered) \| Free plate (Lamb) | HS treats *one selected layer* (via the dropdown that appears) as an infinite half-space on its own — a quick way to check one material in isolation. Full stack requires ≥2 layers. Free plate treats the *whole* current stack as a free-standing plate, both surfaces traction-free, and computes its Lamb modes (Ritz only — no Green's function exists for a free plate here); see Step 7 below. |
| **Family** (Layered only) | isoWavevec \| isoFreq | isoWavevec = fixed wavelength, Ritz+Green (the everyday case). isoFreq = fixed real frequency, Green-only, normalised-`omega0` convention — see `docs/background.md` for when to reach for which. Hidden for HS and Free plate. |
| **Calc type** | Single \| AngSweep \| 2-axis sweep | One angle, an angle sweep, or a 2D sweep (angle × wavelength/frequency for Layered, angle × wavevector for Free plate). 2-axis sweep is disabled in HS mode (no second axis to sweep in a bare half-space). |

The three generic fields below relabel themselves for whichever combination is active (e.g. `angle [deg] (opt, 0)` for HS Single, `wavelength` + `angle range` + `vMax (opt)` for isoWavevec AngSweep, `realWavVec (opt)` + angle + `realThickness (opt)` for Free plate Single/AngSweep, with `wavVecs` and `realThickness` both required for a Free plate 2-axis sweep) — read the live label, not this table, for what to type.

**Advanced** (collapsed by default): `calcMode` (`full` runs Christoffel+Green too; `ritzOnly` is cheaper), `retentionLevel` (eigenvector RAM policy), and the classification the run starts with — **Classification** (`Green` / `Refined` / `Ritz`), **Standard** (`method default` / `high` / `low` / `custom`) and **Custom value** (used with `custom` only). The methods on offer follow the other settings: only `Ritz` on `ritzOnly`, no `Refined` on a layered stack, and no classification at all for isoFreq and Free plate. The method you picked comes back when the settings allow it again (switching `calcMode` from `ritzOnly` back to `full`, or from a layered stack back to a half-space for `Refined`); the **Standard** tooltip lists each method's default. You can switch method and standard later in any viewer without recomputing — see "Choosing and rerunning the classification" below. There is no Green-axis resolution setting: in `full` mode the Green's function is computed on a fixed velocity grid (6 m/s on a half-space, 4 m/s on a layered stack). Free plate ignores `calcMode` and accepts only `retentionLevel` `'full'`/`'none'` (see "Known limitations").

## Step 3 — Run

**[Run]** dispatches to the matching `Calc.*` class (mirrors `docs/reference/functionReference.md`'s half-space/layered tables exactly) and assigns the result to a base-workspace variable (default `guliResult`, editable beside the button).

Next to Run is **"Open viewer after Run"**, checked and enabled whenever the current combination has a viewer to open:

| Combination | Result class | Opens |
|---|---|---|
| HS Single, HS AngSweep | `CalcSingleHS`, `CalcAngSweepHS` | `Gui.combinedViewer` |
| isoWavevec Layered Single, AngSweep | `CalcSingleLayered`, `CalcAngSweepLayered` | `Gui.combinedViewer` |
| isoWavevec Layered 2-axis sweep | `CalcWavevecAngSweep2Layered` | `Gui.dispersionViewer` |
| isoFreq Layered Single, AngSweep | `CalcIsofreqSingleLayered`, `CalcIsofreqAngSweepLayered` | `Gui.isofreqCombinedViewer` |
| isoFreq Layered 2-axis sweep | `CalcFreqAngSweep2Layered` | `Gui.freqAngSweep2LayeredViewer` |
| Free plate Single, AngSweep | `CalcSinglePlate`, `CalcAngSweepPlate` | `Gui.combinedViewerPlate` |
| Free plate 2-axis sweep | `CalcWavevecAngSweep2Plate` | `Gui.dispersionViewerPlate` |

The label and enabled state update live as you change the toggles above, so it always reflects what Run is actually about to produce. If disabled, `lastResult`/`lastViewer` on the `StructureControlGUI` handle are still there to inspect or hand to a viewer manually.

## Step 4 — The combined viewer (`CombinedViewerGUI`)

Opened for any HS or isoWavevec Layered Single/AngSweep result (isoFreq and Free plate results have their own viewers, Steps 6–7). One window, opened maximised, split by a thin divider:

**Top half — modal shapes.** A depth×lateral displacement-field plot per layer, plus:
- **Angle** dropdown + `<< >>` (AngSweep results only — a Single result has one angle, so these are no-ops)
- **Mode** `<< >>` — steps through computed Ritz modes at the current angle, with the selected mode's velocity beside it (`vel: 2.91 km/s`); stepping angle tries to follow the *same physical wave* (SAW/limiting-type identity), not just the same mode index, since Ritz mode ordering can shuffle with angle
- **Source**: `Ritz modal shape` (the eigenvector, unconstrained and re-rendered) or `Green's response` (a `Force dir` row appears below — continuous velocity rather than a discrete mode; half-spaces and layered stacks)
- **Display**, the dropdown beside Source: `Components (u_x, u_y, u_z)` — one signed panel per displacement direction — or `Polarisation (RGB)` — one panel per layer, coloured by the direction of the displacement (red u_x, blue u_y, green u_z; key on the right) with brightness = \|u\|, white where nothing moves. The same dropdown is in every viewer that shows modal shapes (Plate, isoFreq, dispersion viewers)
- **SAW-like modes** / **Limiting-velocity modes** panels — a **[Go To]** button per automatically-classified candidate (each SAW-like entry also shows its method, score in dB and grade), plus **Mark current mode as: [SAW] [lim]** on the panel's bottom line to override the automatic classification by hand when it picked the wrong branch (session-local, not saved back to the result object)

**Bottom half — angular map** (AngSweep results) or **single-angle curve** (Single results). Its top bar has two rows: the graph type (and criterion) on the first, the classification bar on the second (see "Choosing and rerunning the classification" below).

| Graph type | Shows |
|---|---|
| **Classification diagnostic** / **curve** (default) | The classified surface waves on a response background, with the Christoffel limiting-velocity reference. Marker colour = type, marker size = score (bigger = stronger response), filled = `strong`, open = `weak`. The background is `\|G13_eff\|` on a `full` result and the Ritz spectrum on a `ritzOnly` result |
| **Green's-function map** / **curve** | `\|G13_eff\|`, the ±k-symmetrised `G13` — the spectrum the Green classification is decided on. There is no component dropdown; other components come from a script, e.g. `Vis.plotGreenMap(sweep, [], [], [], [], 'G33')`. Not available on a `ritzOnly` result (the status line says so) |
| **Ritz-mode map** / **curve** | Every computed Ritz mode, shaded (map) or placed (curve) by the **Criterion** chosen in the dropdown beside it (below) |
| **Ritz spectrum map** / **curve** | The Ritz spectrum \|S(v)\| — the `G13`-like response summed over the Ritz modes, which the Ritz classification scores. The Ritz-side counterpart of the Green's-function view, on the same log-normalised scale for direct comparison; on a `ritzOnly` result it is the only response view. It is rebuilt from data stored with the result, so drawing it costs no solve. The curve marks the Ritz candidates (▼ at the spectrum's peaks, ○ at the Ritz modes that carry them) |

**Criterion** (Ritz-mode map/curve only):
- **Ritz score** (default) — each mode coloured by the Ritz classification's score (dB) of the candidate it carries; modes that carry no candidate are drawn faint. On the curve the scored modes are joined by a line with small markers, and the others sit as grey dots on a baseline below them. This shows where the Ritz rule looked, accepted or not.
- **freqDisp** — how non-dispersive the mode is (1 = phase and group velocity equal). The Ritz classification vetoes a half-space candidate whose supporting mode has `freqDisp` below 0.3: the map shades that range red and the curve draws the level as a dashed red line, and a click there reads "vetoed". On a layered result the item is greyed and the status line notes that it is not used on layered stacks — the film makes true guided modes dispersive.
- **depthInsensitivity** — how little the mode's velocity changes with the depth of the computational domain (large for modes confined near the surface); a diagnostic kept for comparison.

Navigating:

- Click the map/curve to snap to the nearest computed mode — this drives the modal-shape half above automatically, and stepping the modal-shape half in return marks a crosshair below (bidirectional, both halves stay in sync). For a classified mode the status line also reads its method, score in dB, grade and supporting mode, and on the Ritz-mode view the criterion value
- **Show types** checkboxes filter the classified points by type in the classification view; on the other graph types, ticking a type overlays its classified points as rings (faint rings for `weak` entries) — off by default
- Axis-range fields (angle or velocity on x, velocity or value on y, colour range for the 2D map) apply on edit, and keep their values across a drill-down, a mode switch or a wavelength step; a script can set them too (`setAxisRange`, see [Standalone use](#standalone-use))
- In the single-angle view the selected mode stays marked when you change the graph type or criterion; the **Mode crosshair** checkbox hides the marker
- **[Export figure]**, beside "Show types", opens the map or curve as shown (hidden types stay hidden, your axis ranges kept) in an ordinary MATLAB figure window, to save from its File menu or with `exportgraphics` in any format

**Drill down to one angle, from an AngSweep map:** the bottom half's **[View single angle ->]** button swaps the 2D map for a 1D curve at the angle you last clicked (or the first angle) — the single-angle result is already part of the sweep, so nothing is recalculated. Click-to-select keeps driving the modal-shape half while drilled in, and the **Graph type** and **Criterion** dropdowns switch the curve (Green's-function map → Green's-function curve, and so on); what you pick there is what the map shows when you go back. **[<< Back to angular sweep]** restores the values of the angle/velocity/colour range fields from before the drill-down (a zoom done with the axes toolbar alone is not captured).

`Gui.combinedViewer(sweepResult, windowMode)` also takes an explicit `windowMode`: `'combined'` (default, everything above in one window), `'split'` (modal-shape window + map/curve window, two separate windows), `'fullSplit'` (modal-shape controls, modal-shape plots, and the map/curve each in their own window — three total). Useful on a multi-monitor setup or when the combined window feels cramped.

## Step 5 — The dispersion viewer (`DispersionViewerGUI`)

Opened for a `CalcWavevecAngSweep2Layered` result (isoWavevec 2-axis sweep, one fixed angle × a wavelength range), run with `retentionLevel='full'` (the branch tracking needs the eigenvectors). Modal shapes on one side (stepped by **Wavelength** dropdown + `<< >>` and a numeric **Branch** field, tracking a physical branch continuously across wavelength steps via `DispersionBranches.track` — no "follow the same index" jump), a dispersion-curve map on the other (**Graph type** selects the curve set, click snaps the branch/wavelength selection). A classification bar above the curves (Green | Ritz) reclassifies every wavelength of the sweep at once and highlights the newly classified branches; the branches themselves do not change. Same `windowMode` (`'combined'`/`'split'`) toggle as the combined viewer.

## Choosing and rerunning the classification

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface waves in strongly anisotropic materials. Thresholds, scores and the interface may change.

Which modes count as surface waves is decided by a classification: one of three **methods**, at one **standard**. Every viewer of an HS or isoWavevec Layered result has a classification bar for it (in the combined viewer, in the bottom half): **Classification method**, **Standard** (with a value field for `custom` and the value in dB beside it), **[Reclassify]**, and a status line naming the active classification. Choose a method and standard, press **[Reclassify]**, and the map or curve, the modal-shape panels and any other linked view redraw with the new classification. The classification is re-thresholded from data already stored with the result, so this takes no recomputation. (How the bar is built and refreshes the linked views: `docs/reference/architecture/gui_viewers.md`, "The classification bar and the refresh chain".)

**Methods:**

- **Green** (default) — scores each peak of the Green's function `|G13_eff|`; needs `calcMode` `full`. The reference method.
- **Refined** — the Green score of the response averaged over the directions a finite probe spot excites, as in a measurement; half-space `full` results only. Its first use solves the extra directions (for a 1° sweep about half the time of the calculation itself); switching back and forth after that is free.
- **Ritz** — scores the Ritz spectrum, so it needs no Green's function: the only method on a `ritzOnly` result, and the cheapest. An approximation of Green, with its own limits (`docs/pitfalls.md` MTH9); it reports no limiting-velocity entries, and on a layered stack only trapped waves (leaky waves need Green).

Methods a result cannot use are shown but greyed, and selecting one says why (Green and Refined on a `ritzOnly` result: recompute with `calcMode` `full`). On a layered result Refined is not in the list at all — it does not exist for stacks.

**Standards:**

- **high** — everything a good, high-averaging map from transient grating spectroscopy (TGS) shows. The default for Green and Refined on a half-space.
- **low** — the distinct, stronger waves only; for poor TGS data or spatially resolved acoustic spectroscopy (SRAS). The default for Ritz and for every layered method. On a layered stack, Love-type waves need `high`.
- **custom** — a number on the method's own score scale (a log10 amplitude ratio; the label shows it in dB). The scales differ between methods, so a number tuned for one means nothing for another, and changing the method resets the standard to that method's default. Prefer `high` / `low` unless you are deliberately scanning a threshold (`docs/pitfalls.md` MTH6; the values: `docs/reference/architecture/forwardCalcs_modeClassification.md`, section 1).

**When to use which.** Start with Green at its default. Use `low` when comparing with SRAS or a noisy TGS map, `high` with a good TGS map. Use Refined to compare with a measured map in detail, for instance where a branch moves quickly with angle and its end depends on the probe spot. Use Ritz when a quick `ritzOnly` run is enough, or to see what the cheaper method would report. How each method compares with measured maps is in the validation table of `docs/reference/architecture/forwardCalcs_modeClassification.md`.

**Strong and weak entries.** Each classified surface wave carries a grade: `strong` if it also passes the method's `low` standard, `weak` if it passes only the looser active one. In the classification view weak entries are open markers, strong ones filled (faint and full rings in the overlay); at `low` every entry is strong. Weak entries are the ones a poorer measurement would probably miss. The marker size grows with the score, so the strongest waves stand out within one classification — sizes are not comparable between methods.

**"Stale classification".** A result saved with an older classification rule set still opens, with its maps, modes and Green and Ritz data, but the status line shows this banner and **[Reclassify]** is disabled: its stored scores were made by other rules and cannot be re-thresholded consistently. Recompute it from its structure (re-run it from the structure window, or `calc()` / `calcAll()` in a script) to classify it again.

## Step 6 — The isoFreq viewers (`IsofreqCombinedViewerGUI` / `FreqAngSweep2LayeredViewerGUI`)

The isoFreq family (Green's-function only — no Ritz solve, so no Ritz/classification controls anywhere in this family) has its own parallel set of viewers, structurally mirroring Steps 4–5 above but simplified: one graph (the Green `|G₃₃|` field) instead of a three-way graph-type choice, and the modal-shape half shows only the on-demand **Green's response** displacement field (no Ritz side to toggle to — `Force dir` dropdown is always visible).

**isoFreq Single/AngSweep** → `Gui.isofreqCombinedViewer` (`IsofreqCombinedViewerGUI`): same top/bottom split as `CombinedViewerGUI` — modal shapes (**Angle** stepper, a **kParRange** field with `<<`/`>>` stepping one kParRange sample and the matching velocity `v = ω₀/k∥` shown beside it, in place of the Ritz-mode field, **Force dir**) above, a **kParRange vs. `|G₃₃|`** heatmap (AngSweep) or curve (Single) below. Click the map/curve to drive the shape above; the shape's own steppers mark a crosshair in return. AngSweep's **[View single angle ->]** drill-down works exactly as it does for `CombinedViewerGUI`.

**isoFreq 2-axis sweep** → `Gui.freqAngSweep2LayeredViewer` (`FreqAngSweep2LayeredViewerGUI`): an **omega0** stepper on top rebuilds the same shape+map pair against a different frequency slice on every step (no recomputation — each slice was already `calc()`'d as part of the sweep). Single window only; every step resets to the map view (there is no second "mode" to preserve, unlike `DoubleSweepViewerGUI`'s angular/dispersion switch).

## Step 7 — The Plate viewers (`CombinedViewerGUIPlate` / `DispersionViewerGUIPlate`)

Opened for a **Free plate (Lamb)** result. The Plate family mirrors Steps 4–5 window for window, but is a separate set of classes, not a mode of the Layered ones, because two things are genuinely different for a free plate:

- **No Green's function exists for a free plate**, so every Green-related control — the Source toggle, Force dir, Green's-function graph type, Green background — is absent (not merely greyed out).
- **Every computed Ritz mode is a real physical branch** (both plate surfaces are genuinely free, unlike SAW's artificial rigid bottom), so there is no SAW-recognition or classification step. What each mode gets instead is a **mode type** — `Asymm`, `Symm`, `TH-Asymm` or `TH-Symm` (TH = shear-horizontal-like) — plus a **confidence** in [0, 1] for how clearly it beat the runner-up label.

**Free plate Single/AngSweep** → `Gui.combinedViewerPlate` (`CombinedViewerGUIPlate`): the same top/bottom split and `windowMode` (`'combined'`/`'split'`/`'fullSplit'`) as `CombinedViewerGUI`.
- **Top half — modal shapes:** **Angle** dropdown + `<< >>`, **Mode** field + `<< >>`, and a line showing the current mode's type and confidence (in place of the SAW-like/limiting panels). Stepping angle keeps the same *mode number* rather than following the same physical wave as Step 4 does — mode order by velocity can shuffle with angle, so re-check the type/confidence after a step.
- **Bottom half — one graph only** (Ritz-mode map for AngSweep, Ritz-mode curve for Single), with a **Criterion** dropdown: **Mode type** \| **Confidence**. On the AngSweep map the criterion sets the colour — mode type as a four-colour categorical legend, or confidence on a continuous scale with a colour bar (the colour-range fields are enabled only then). On the Single curve (velocity on x) it sets the y-axis instead, while marker colour always shows the mode type. There is no per-type "Show types" filter. Click-to-select, the crosshair link to the modal-shape half, and the **[View single angle ->]** drill-down (with its axis ranges saved and restored) work exactly as in Step 4.

**Free plate 2-axis sweep** → `Gui.dispersionViewerPlate` (`DispersionViewerGUIPlate`): the Plate counterpart of Step 5's dispersion viewer, at one fixed angle (shown in the top bar). Modal shapes are stepped by **Wavevector** dropdown + `<< >>` and a numeric **Branch** field; the curve side has a **Curve type** dropdown (e.g. `k-vphase`, `freq-vgroup`), **Colour by** `type` \| `branch`, and **Show Christoffel bounds** (bulk-wave limiting velocities as horizontal reference lines, on by default). Needs the sweep run with `retentionLevel='full'` (the Advanced default), since branch tracking follows each branch by eigenvector overlap.

`Colour by: branch` matters for a stack that isn't mirror-symmetric about its mid-plane (different materials or thicknesses on the two faces): mode-type labels there are a heuristic that can become ambiguous near branch crossings, and a low-confidence point falls back to a stable "branch N" label instead. Separately, multi-layer Plate frequencies, group velocities, and ZGV points agree with GEWtool, an independent code, but have not been compared with a measured multi-layer plate (`docs/reference/architecture/plate.md`).

As with the Layered family, a genuine double sweep (several angles *and* several wavevectors) still auto-opens only the dispersion viewer; its angular-map ↔ dispersion-curve switcher is reached directly, via `Gui.doubleSweepViewerPlate(result)`.

## Standalone use

Every viewer above is independently constructible from a script, without going through `StructureControlGUI` at all — useful if you already have a result sitting in the workspace:

```matlab
win = Gui.combinedViewer(sweepResult);   % or Gui.dispersionViewer / Gui.angSweepViewer / Gui.modalShapeViewer / Gui.singleAngleViewer /
                                          %    Gui.isofreqCombinedViewer / Gui.freqAngSweep2LayeredViewer / Gui.isofreqAngSweepViewer /
                                          %    Gui.combinedViewerPlate / Gui.dispersionViewerPlate / Gui.doubleSweepViewerPlate / ...
win.show();
```

Two standalone windows (e.g. `Gui.modalShapeViewer(sweep)` + `Gui.angSweepViewer(sweep)`, opened separately rather than via `Gui.combinedViewer`) aren't linked automatically — wire them by hand:

```matlab
win1 = Gui.modalShapeViewer(sweep);
win2 = Gui.angSweepViewer(sweep);
win1.show();  win2.show();
win2.OnModeSelected     = @(angle, modeNo, ~, rawVel) win1.goToSelection(angle, modeNo, rawVel);
win1.OnSelectionChanged = @(angle, vel) win2.showCrosshairAt(angle, vel);
```

Axis ranges can be set from the script as well, e.g. to open a figure at the same zoom every time. `setAxisRange` fills the side-panel min/max fields exactly as typing would; `NaN` keeps that limit, `'c'` is the colour range of a map:

```matlab
win = Gui.doubleSweepViewer(sw2);
win.setAxisRange('y', [NaN 7]);    % velocity up to 7 km/s, lower limit unchanged
win.setAxisRange('c', [-3 0]);     % colour range
win.show();
```

The combined and double-sweep viewers pass the call on to their map, or to the dispersion curves in dispersion mode, where `'y'` is whatever the selected curve type plots.

## Viewing a result computed elsewhere

A calculation run without a display (a `-batch` run, a remote machine, an overnight script) cannot open any window — and the result disappears when that MATLAB session ends. Save it instead, and open the viewer later from any interactive session:

```matlab
% in the calculation script
sweep = sweep.calc();
path  = Calc.saveResult(sweep, [], 'Cu(001), legO 40');   % [] -> results/CalcAngSweepHS_<name>_<timestamp>.mat

% later, in your own MATLAB
Gui.viewSaved(path)        % or Gui.viewSaved() for a file dialog in results/
```

`Gui.viewSaved` opens the default viewer for the result's class (`Gui.viewResult` does the same for a result already in memory). It is the window the **Open viewer after Run** checkbox would open, except for a double sweep (several angles and several wavelengths), where it opens the double-sweep viewer and Run opens only the dispersion viewer. Saving drops only the angle-independent Ritz precompute, which is rebuilt automatically the first time a recompute needs it (e.g. a modal shape under `retentionLevel='saw_only'`/`'none'`); everything the viewers read is stored. A result saved before the `Calc*` classes changed may not load — see `docs/pitfalls.md` OTH3.

## Known limitations

- Running a *film* layer on its own through **Single layer (HS)** uses a fixed `relDepth=20` (with that tab's own `legOrder`), not an editable field — the film tab's thickness field has no meaning for a half-space. The Substrate tab's own `relDepth` is used as typed.
- In Free plate mode the Advanced `retentionLevel` dropdown still offers `'saw_only'`, which the plate classes do not accept (the status line shows the error) — choose `'full'` or `'none'`.

## Where to go next

- **`docs/quickStart.md`** — the scripting equivalent of everything above, for parameter sweeps and repeatable work.
- **`docs/reference/functionReference.md`** — one-line-per-function index, including the full `Gui.*`/`Calc.*`/`Vis.*` namespaces this guide's windows sit on top of.
- **`docs/pitfalls.md`** — gotchas that have each caused a real bug once (several apply just as much to GUI-entered values as scripted ones, e.g. `legOrder`/`relDepth`).
