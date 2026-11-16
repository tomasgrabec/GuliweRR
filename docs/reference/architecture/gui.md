# GUI — architecture and design rationale

`guliwerr/gui/` (`StructureControlGUI` plus the viewer handle classes, fronted by the `@Gui`
namespace launcher). This file covers the launcher pattern and `StructureControlGUI` — the entry
point: a click-driven alternative to scripting against `Calc.*` directly, for exploring one
material or checking a single point (`docs/quickStart.md`/`docs/guiGuide.md`). The viewer windows
that `StructureControlGUI` hands a finished calculation to are covered in `gui_viewers.md` and
noted here only where the hand-off itself matters.

## `@Gui` — the namespace launcher

`guliwerr/gui/@Gui/Gui.m` is a thin static-method namespace (`Gui.structureControl(...)`,
`Gui.combinedViewer(...)`, etc.) over the GUI handle classes, which live directly in
`guliwerr/gui/` — the same `@ClassName/` pattern as `@Christoffel`, `@GammaMatrix`,
`@GreenSolver` and `@Calc` (see `forwardCalcs.md`, "`@Calc` — the namespace launcher"). Each
launcher file is a one-line pass-through — `obj = StructureControlGUI(varargin{:})` — so calling
code writes `Gui.structureControl(...)` rather than qualifying each class name. Do not
instantiate `Gui` itself. The two exceptions are not pass-throughs: `Gui.viewResult(calcResult)`
picks the default viewer for any `Calc*` result, and `Gui.viewSaved(file)` loads a result saved by
`Calc.saveResult` and hands it to `Gui.viewResult`.

`Gui`'s launchers pass `varargin{:}` through to constructors with fixed, named argument lists
(`StructureControlGUI(initialStructure)`, `AngSweepViewerGUI(sweepResult, parentContainer, classifyTarget)`,
...). As in `@Calc/Calc.m`, the classdef `Gui.m` declares each static method's argument-count
signature separately from its one-line implementation file, and MATLAB validates a call against
the classdef's declared signature, not against what the implementation does with `varargin`.
A new viewer-constructor argument therefore has to be added to `Gui.m`'s declared signature as
well, or the call fails with "too many input arguments" (the same applies on the `@Calc` side,
where a missing declaration can leave an argument silently not forwarded).

## `StructureControlGUI` — combined builder + launcher window

One `uifigure` window with two halves that share a single object's state: a tab per layer
(Substrate, always present, plus film layers added/removed on demand) for building a
`Structure`, and a control panel that dispatches that `Structure` to whichever `Calc.*`
orchestrator matches a set of toggles. There is no separate "hand the structure to the
launcher" step — the same handle object owns both halves, and every field edit on the left
writes straight into the `Structure` the right side will use.

### Why `uifigure`, not classic `figure`/`uicontrol`

The window is built programmatically from `uifigure` components (there is no App Designer `.mlapp` file). `uigridlayout` (declarative row/column layout) replaces manual normalised-position arithmetic, and `uitabgroup`/`uitab` is a genuine tabbed-panel widget replacing manual `uipanel`-visibility toggling. Callback logic can be exercised without a display, but layout correctness — do labels clip, does a title bar overlap a button, is the Run button still inside the window — needs the window rendered.

**Checking a layout from a batch run.** Run MATLAB with `-batch` but without `-nodisplay` (under `-nodisplay`, `exportapp` fails for a `uifigure`). `exportapp(fig, 'shot.png')` then captures the whole window, widgets included, even if the figure was never made visible, so nothing appears on screen. Keep the viewer hidden, give it the size to check, let it settle, and export, e.g. `w = Gui.viewSaved(path, 'combined', false); w.fig.Position = [0 0 1800 1100]; drawnow; pause(3); exportapp(w.fig, 'shot.png'); close(w.fig)`.

### State management: `structureModel` as the single source of truth

```matlab
properties (SetObservable)
    structureModel  Structure = Structure()
end
```

Every widget edit — a Cij cell, density, symmetry, direction-matrix row, depth/thickness,
`legOrder` — calls `syncPanelToStructure(k)`, which rebuilds that tab's `Material`/
`CalcParams` from its current field values and calls `Structure.AddLayer`/`ModifyLayer` on
`structureModel`. There is no separate "commit" step: `structureModel` is always the live,
current structure, and `[Build Structure]` (`onBuild`) does nothing but `assignin('base', ...)`
a copy of whatever is already there — it exists purely to hand a reachable variable to the
command line, not to trigger a build. `structureModel` is declared `SetObservable`, and the
constructor registers a `PostSet` listener (`addlistener(obj, 'structureModel', 'PostSet',
@(~,~) obj.onStructureChanged())`) that repopulates the half-space (HS) branch layer dropdown and toggles
the "Full stack" mode's availability (disabled, forced back to HS, below 2 layers) every time
the structure changes shape — this propagation is driven by the property system, not manually
threaded through every one of the edit callbacks that can change `nLay`.

`recomputeAndRefreshCij`/`onSymmetryChanged`/`onDirMatEdited` all *rebuild the Material from
scratch* on every edit (re-run `ElastConsts(vec, symmNo)`, re-run `Material.CompleteDirMat`)
rather than patching a specific field — cheap enough per keystroke, and it means invalid
intermediate states never propagate: `syncPanelToStructure` returns early without writing to
`structureModel` if `ElastConsts` construction throws (an inconsistent partial Cij entry, for
instance), so the model only ever holds physically valid materials.

`panels` (one entry per open tab: direction-matrix fields, `dirMat`, density field, symmetry
dropdown, the 6×6 Cij field grid, depth field, `legOrder` field) is declared as a plain
`cell`, not a `(1,:) struct` array, because a size/class-constrained struct property starts as a
0-field struct and rejects the first assignment of a struct with fields.

### Building a layer's `Material`

Each tab's symmetry dropdown drives which of the 6×6 Cij grid's upper-triangle cells are
editable, via `ElastConsts.IndepVoigtPositions(symmNo)` (`guliwerr/core/@ElastConsts/`) —
choosing a symmetry class enables exactly the independent constants that class needs and
disables (greys out) the rest; `recomputeAndRefreshCij` then re-derives the dependent cells
from the current independent ones and writes them back, so a triclinic-to-cubic switch
immediately shows the cubic-consistent values rather than stale unrelated numbers.

Direction-matrix completion (`onDirMatEdited`) calls `Material.CompleteDirMat` on the user's
2×3 input (surface normal + first/propagation direction), which orthogonalises the second row
against the first if needed. `CompleteDirMat` itself has no `wasAdjusted` output, and the GUI
deliberately does not add one: it repeats the same `dot(n1, n2) > 1e-3` check (matching
`CompleteDirMat`'s internal tolerance) purely to decide whether to show a "First direction
adjusted for orthogonality" status message — a GUI-only presentation concern kept out of
`Material`'s signature.

### Adding and removing layers

`[+ Add layer]` always appends a new film-layer tab on top of the stack — `Structure.AddLayer`
is itself append-only, with no insert-in-the-middle method, so the widget-side operation
matches what the underlying model can actually do. `maxGuiLayers` (default 5, film layers only,
Substrate is separate) is a soft usability cap, enforced by disabling the button with a status
message, not a construction-time limit — `uitabgroup` builds tabs on demand, so no pre-built
maximum is needed.

`[- Remove active layer]` removes whichever tab is *currently open* (`activeTabIndex`, matched
against `tabGroup.SelectedTab`), not always the top layer — silently removing the top tab while
looking at a different one would be surprising. The Substrate (index 1) is never removable,
checked directly against the active index. Removing a middle layer shifts every layer above it
down by one in `structureModel` (`Structure.RemoveLayer` supports an arbitrary index), and
`renumberLayerTabs` keeps the film-layer tab titles ("Layer N") contiguous.

### Pre-populating from an existing `Structure`

`StructureControlGUI(initialStructure)` fills every tab from a non-empty `Structure`: one extra
film-layer tab per layer beyond the Substrate (`loadFromStructure`), then each tab's
dirMat/density/symmetry/Cij/relDepth-or-realThickness/legOrder fields
(`populatePanelFromLayer`). That method finishes by calling `onSymmetryChanged`/`onDirMatEdited`
— the same two calls `addPanel` ends with — so the recompute and sync into `structureModel`
happen exactly as for a manual edit, rather than by writing to `structureModel` directly.

## The control panel — user-driven dispatch, not automatic on `nLay`

The **Calculation target** toggle (Single layer (HS) / Full stack (Layered) / Free plate (Lamb))
is a deliberate user choice, not something the GUI infers from
`structureModel.nLay`. Picking any one layer — Substrate or a film — through the HS-branch
dropdown runs that layer's `Material` alone through `CalcSingleHS`/`CalcAngSweepHS`, as if it
were an infinite half-space on its own; nothing about having built a 3-layer stack forces
treating it as layered. This mirrors `forwardCalcs.md`'s own point that half-space and layered
are genuinely separate class families with separate mode taxonomies — the GUI toggle is where
that split becomes a literal user decision. Free plate mode uses `structureModel` directly and
in full, like the Layered branch's isoFreq/2-axis paths do (never `buildHSStructure`'s
single-layer extraction). It hides both the HS layer-selector dropdown and the Layered-only
isoWavevec/isoFreq **Family** toggle — Plate has no isoFreq family, since there is no Green's
function to scan (`docs/reference/architecture/plate.md`) — and, unlike HS, keeps the **2-axis sweep**
calc type enabled, since `CalcWavevecAngSweep2Plate` exists. "Full stack" is disabled whenever
`nLay < 2`, matching `CalcSingleLayered`/`CalcWavevecAngSweep2Layered`'s own
constructor-level validation.

`buildHSStructure` packages the selected layer as a synthetic one-layer `Structure`, reusing
the Substrate tab's own `relDepth` directly if the substrate itself was selected, but falling
back to a **hardcoded `relDepth=20`** if a film layer was selected — a film layer's natural
parameter is `realThickness` (a real-unit quantity converted via `ApplyWavelength`), which has
no meaning for a bare half-space run with no wavelength in play. 20 is a plausible generic depth,
not one derived for the material in that tab; there is no editable field for it.

### The three generic calc fields

Rather than a dozen mostly-hidden named fields (one per possible parameter across every
`(layerMode, family, calcType)` combination), the control panel exposes exactly three generic
`uieditfield`s whose labels and visibility are rewritten by `refreshCalcFields` for whichever
combination is currently active — e.g. `angle [deg] (opt, 0)` for HS Single, or
`{wavelength, angle range, vMax}` for an isoWavevec AngSweep. `refreshCalcFields` is wired to
all three toggle groups (`onLayerModeToggle`, `onFamilyToggle`, `onCalcTypeToggle`), so any
change to any of the three re-evaluates all three labels together rather than patching one at
a time. Each label is stacked directly above its field (a `[6×1]` grid) rather than beside it in
a fixed-width column, where longer labels like `angle range, e.g. 0:2:180` would clip illegibly
(a layout defect only visible on screen, not to any headless check).

`onLayerModeToggle` sets `layerModeHSBtn.Value = true` directly in the constructor path
(`buildControlPanel`), which does **not** raise `SelectionChangedFcn` — MATLAB only fires that
callback on a user-driven or programmatic *selection change* through the group's own API, not
a direct property write on a member button. `calcType2AxisBtn` (invalid in HS mode) is disabled
inside the `onLayerModeToggle` callback, so `buildControlPanel` calls `obj.onLayerModeToggle()`
explicitly once at the end, after every toggle group exists, to sync the initial state. The
general rule: any `Value = true` set directly on a `uiradiobutton` needs an explicit follow-up
call to whatever `SelectionChangedFcn` would otherwise have done.

The Layered branch additionally toggles a **Family** group (isoWavevec / isoFreq), mirroring
`forwardCalcs.md`'s isoWavevec/isoFreq split at the orchestration layer: isoWavevec fields are
`{wavelength(s), angle(s), vMax}`, isoFreq fields are `{omega0(s), angle(s), kParRange}` — the
GUI never lets the two families' fields blend, because they map to structurally different
`Calc.*` constructors with no shared parameter meaning. Free plate mode keeps the isoWavevec
field roles (1 = wavelength/wavevector axis, 2 = angle axis, 3 = the extra parameter, `vMax`
there), substituting its own real-unit parameters: `{realWavVec, angle(s), realThickness}` for
Single/AngSweep (all optional), `{wavVecs, angle range, realThickness}` for the 2-axis sweep
(`wavVecs`/`realThickness` required).

### `vMax` is optional

`vMax` (field 3 of an isoWavevec layered calculation, labelled `vMax (opt)`) only extends the Green
velocity axis beyond its default end, 1.05 × the fastest longitudinal limiting velocity of any
layer; an empty field passes `[]`. `kParRange` (field 3 of an isoFreq calculation) is required.

### Advanced options

`calcMode` (`full`/`ritzOnly`), `retentionLevel` (`full`/`saw_only`/`none`) and the
classification controls (below) sit in a `uipanel` collapsed by default behind an
"Advanced ▸/▾" disclosure button — no dedicated MATLAB collapsible-panel component exists, so
toggling the panel's `Visible` property and flipping the button's own text is the whole
mechanism. There is no Green-axis field: the step is the Calc classes' `greenVelStep`
(default `GreenScore.DEFAULT_STEP` = 6 m/s on a half-space; a fixed 4 m/s on a layered stack),
which the classification does not depend on and the GUI does not expose. The disclosure glyphs are built at
runtime via `char(9656)`/`char(9662)`, since single-quoted MATLAB char arrays do not process
`\x`-style escapes.

### Classification method and standard

> **Experimental.** Surface-wave classification is checked against the measured maps published in Grabec et al. 2024 and the UTGS paper (arXiv:2510.10696), and there is no general agreement on which spectral peaks count as surface waves in strongly anisotropic materials. Thresholds, scores and the interface may change.

Three controls in the Advanced panel choose the classification a run starts with
(`forwardCalcs_modeClassification.md`): **Classification** (`Green` / `Refined` / `Ritz`),
**Standard** (`method default` / `high` / `low` / `custom`) and **Custom value** (a number on the
method's own score scale, enabled for `custom` only). `onRunClicked` passes them to the Calc
object as `sawClassMethod` / `sawThreshold` before it runs (`classificationSettings`,
`applyClassification`): `method default` becomes `sawThreshold = []`, the method's own default
standard (`'high'` for Green and Refined, `'low'` for Ritz and for every layered method). The item
is called `method default`, not `default`, because a graphics `Value` of `'default'` means "reset
this property to its default" to MATLAB and cannot be selected as an ordinary item.

`refreshClassificationOptions` keeps the controls consistent with the rest of the panel, and runs
whenever `calcMode`, the standard, or one of the three toggle groups changes:

- **`calcMode = 'ritzOnly'`:** the method shows `Ritz`, disabled — the only method without a
  Green's function. Run passes `'ritz'` whatever the dropdown held.
- **Layered:** `Green` / `Ritz` only (there is no spot score for stacks, so no `Refined`).
- **Free plate and isoFreq:** both controls disabled — neither family has a SAW classification.

The method lists come from `CalcCommon.classificationMethods`, the list the Calc classes validate
against. The user's own pick is remembered (`userMethod`) and shown again as soon as the controls
offer it: back from `ritzOnly` to `full`, or from a layered stack to a half-space (`Refined`), the
dropdown returns to the method chosen before, so a run never silently classifies with a method the
GUI substituted. The **Standard** tooltip lists each method's default standard, read off
`CalcCommon.defaultThreshold`.

This choice only sets the starting point. Any viewer's classification bar switches method and
standard afterwards without recomputing (`gui_viewers.md`, "The classification bar and the
refresh chain"), so the panel offers no more than the classes' own options.

## Running a calculation and handing off to a viewer

`onRunClicked` reads the three toggle groups and the three generic fields, dispatches to the
matching `Calc.*` static launcher (`Calc.singleHS`, `Calc.angSweepLayered`,
`Calc.isofreqAngSweepLayered`, …, covering the full grid `forwardCalcs.md`'s class-family table
already lists), runs it (`calcAll()`/`calc()` as appropriate for that class), and stores the
result on `obj.lastResult`. As with `[Build Structure]`, the result is also pushed to a
base-workspace variable (`runVarNameField`, default `guliResult`) via `assignin('base', ...)`, so
the computed object is reachable from the command line.

If the **"Open viewer after Run"** checkbox is checked, `openResultViewer` dispatches on
`class(obj.lastResult)` to the matching `Gui.*` viewer launcher and calls `.show()` on it,
storing the window on `obj.lastViewer` (a previously opened viewer is left open, not replaced
or closed — just as `StructureControlGUI`'s own window is not closed by `Build`). The viewer
routes:

| Result class | Viewer opened |
|---|---|
| `CalcSingleHS`, `CalcSingleLayered`, `CalcAngSweepHS`, `CalcAngSweepLayered` | `Gui.combinedViewer` |
| `CalcSinglePlate`, `CalcAngSweepPlate` | `Gui.combinedViewerPlate` (`gui_viewers.md`) |
| `CalcWavevecAngSweep2Layered` | `Gui.dispersionViewer` |
| `CalcWavevecAngSweep2Plate` | `Gui.dispersionViewerPlate` (`gui_viewers.md`) |
| `CalcIsofreqSingleLayered`, `CalcIsofreqAngSweepLayered` | `Gui.isofreqCombinedViewer` (`gui_viewers.md`) |
| `CalcFreqAngSweep2Layered` | `Gui.freqAngSweep2LayeredViewer` (`gui_viewers.md`) |

`CalcSingleHS`/`CalcSingleLayered` (a single propagation angle, no angle axis) route to the
*same* `Gui.combinedViewer` as the AngSweep classes: `CombinedViewerGUI`'s bottom half dispatches
internally on the result's class (a 1D curve viewer for a single-point result, a 2D map viewer
for a sweep; see `gui_viewers.md`), so `StructureControlGUI` need not distinguish Single from
AngSweep when choosing a launcher. `CalcSinglePlate`/`CalcAngSweepPlate` route to
`Gui.combinedViewerPlate`, and `CalcIsofreqSingleLayered`/`CalcIsofreqAngSweepLayered` to
`Gui.isofreqCombinedViewer`, in the same way — three families share this "one combined-viewer
entry point, dispatched internally" shape. Every `(layerMode, family, calcType)` combination has
a route. `describeViewerForCurrentSelection` mirrors `onRunClicked`'s branching *before* Run is
clicked, and `refreshOpenViewerOption` uses it to set the checkbox's label; its
disabled/unchecked branch is unreachable through any current toggle combination and is kept only
as a defensive fallback. `StructureControlGUI` never opens the double-sweep viewers
(`Gui.doubleSweepViewer`/`Gui.doubleSweepViewerPlate`), even for a genuine double sweep: call the
launcher directly, or `Gui.viewResult` / `Gui.viewSaved`, which choose the viewer by the result's
class and, for a 2-axis sweep, by its shape (`gui_viewers.md`, "DoubleSweepViewerGUI").

## Known limitations

- Running a film layer standalone through the HS branch uses the hard-coded `relDepth=20`
  described above rather than a dedicated editable field.
- The Advanced panel's `retentionLevel` dropdown offers `'full'`/`'saw_only'`/`'none'` in every
  mode, but Plate accepts only `'full'`/`'none'` (every Plate mode is real, so there is no
  significant subset to prune to — `docs/reference/architecture/plate.md`). Selecting `'saw_only'` in
  Plate mode surfaces `Calc.*Plate`'s own validation error in the status line (caught by
  `onRunClicked`'s try/catch) rather than being disabled outright.

## References

`docs/guiGuide.md` (user-facing walkthrough); `docs/quickStart.md` (where the GUI sits relative
to scripting); `forwardCalcs.md` (the `Calc.*` classes this window dispatches to, the `@Calc`
launcher pattern, the isoWavevec/isoFreq split); `gui_viewers.md` (the viewer windows
`StructureControlGUI` hands a finished result to).
