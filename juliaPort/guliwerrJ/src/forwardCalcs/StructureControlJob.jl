# =========================================================================================
# StructureControlJob.jl -- the headless half of the structure-control GUI
#
# Port of the NON-widget half of `guliwerr/gui/StructureControlGUI.m` (35 methods, 11
# `Calc.*` entry points): "given these form values, build a Structure and call the right
# Calc.* entry point". The MATLAB controller is split here into this headless job-builder and
# a form. This file has NO Makie/GUI-handle dependency at all and is `include`d directly into
# `module ForwardCalcs`, so it is usable (and testable) without the plotting extension.
#
# The GUI FORM (widgets, a "Run" button, a status label) is `ext/GuiStructureControl.jl`
# (Makie-gated), which builds a `StructureControlSpec` from live
# widget state and calls `build_job` from THIS file -- see that file's own header for the
# hand-off contract.
#
# -----------------------------------------------------------------------------------------
# WHAT THE MATLAB SOURCE DOES
# -----------------------------------------------------------------------------------------
# `StructureControlGUI` is two things wearing one class:
#   (a) a live structure EDITOR -- one tab per layer (Substrate + film layers), each holding
#       a symmetry-class dropdown, an independent-Cij-entry grid (gated by
#       `ElastConsts.IndepVoigtPositions`), a density field, a 2-row direction-matrix input
#       (surface normal + first/propagation direction, completed via `Material.CompleteDirMat`
#       exactly like `Material`'s own `SetDirMat` already does internally), a relDepth field
#       (Substrate) or a realThickness field (film layers), and a legOrder field -- kept
#       live-synced into one `Structure` (`syncPanelToStructure`, `StructureControlGUI.m:
#       477-498`) as the user types.
#   (b) a control panel that reads THREE generic, relabelled calc-input fields (angle/
#       wavelength/omega0/wavVecs/vMax/kParRange/realThickness -- the exact role of each of
#       the three slots changes with the (layerMode, family, calcType) toggle combination,
#       `refreshCalcFields`, `StructureControlGUI.m`) and, on "Run", dispatches to ONE
#       of 11 `Calc.*` entry points (`onRunClicked`, `StructureControlGUI.m`).
#
# THE HEADLESS SPLIT, made here: (a) becomes [`LayerSpec`](@ref) (one raw, still-unvalidated
# layer's worth of form fields) plus [`build_structure`](@ref)/[`build_hs_structure`](@ref)
# (the `syncPanelToStructure`/`buildHSStructure` logic, collapsed into a one-shot "build the
# WHOLE stack from scratch" function, since a job-builder builds fresh every call rather than
# live-editing one `Structure` value incrementally the way the GUI's own event-driven
# `AddLayer`/`ModifyLayer` calls do -- there is no meaningful difference in the RESULT, only
# in whether it is built incrementally or in one pass). (b) becomes [`StructureControlSpec`](@ref)'s
# remaining fields plus [`build_job`](@ref) (the `onRunClicked` dispatch table).
#
# NAMED FIELDS INSTEAD OF MATLAB'S THREE RELABELLED SLOTS (a deliberate deviation). MATLAB's
# three generic `uieditfield`s are relabelled text boxes whose parsed meaning (angle?
# wavelength? omega0? vMax? kParRange?) depends on which of 8 (layerMode, family, calcType)
# combinations is toggled (`docs/reference/architecture/gui.md`, "The three generic calc fields").
# `StructureControlSpec` below gives every physical quantity its own typed field
# (`wavelength::Union{Nothing,Float64}`, `angles::Union{Nothing,Vector{Float64}}`,
# `kParRange::Union{Nothing,Vector{Float64}}`, ...): there is no slot-renaming table to keep
# in sync with `build_job`'s branches, and a value of the wrong kind fails at construction
# instead of being mis-parsed when the job runs.
#
# REQUIRED-VS-OPTIONAL FIELDS -- MOSTLY NOT RE-VALIDATED HERE, DELIBERATELY. MATLAB's
# `parseField` (throws if blank) vs `parseFieldOrDefault` (falls back to a default) choice per
# slot is how each field's requiredness is expressed there. In this port, EVERY downstream
# `CalcXxx` constructor and/or its `calc`/`calcAll` ALREADY implements the identical
# required/optional contract, throwing the exact same class of `ArgumentError` MATLAB's own
# `parseField` throws (e.g. `CalcIsofreqSingleLayered`'s constructor already throws
# `"...missingOmega0"`/`"...missingKParRange"`; any requirement that depends on the effective
# calculation mode is likewise enforced downstream, so the mode-conditional check that bit
# the MATLAB GUI is not re-implemented here at all). `build_job` therefore passes `nothing`
# straight through for every such field and lets the ALREADY-TESTED downstream check fire --
# duplicating it here would be the same "hand-copied set of checks that must agree" trap
# `StructureControlGUI.selectionState`'s own docblock warns about, one level down. The ONE
# exception is `wavelength` (isoWavevec `:single`/`:angSweep`): it is consumed by
# [`ApplyWavelength`](@ref) BEFORE any `CalcXxx` constructor ever sees it (mirroring MATLAB's
# own `obj.structureModel.ApplyWavelength(wavelength)` call inside `onRunClicked`, BEFORE
# `Calc.singleLayered`/`Calc.angSweepLayered`), so a missing value has no downstream
# constructor to throw a clear error for it -- `build_job` checks it explicitly, for exactly
# this reason and no other.
#
# NO "OPEN VIEWER"/"ASSIGN TO BASE WORKSPACE" HERE -- BOTH ARE INTERACTIVE-SESSION CONCEPTS
# WITH NO HEADLESS ANALOGUE. `onRunClicked`'s trailing `assignin('base', varName, ...)` and
# `openResultViewer` (launching a MATLAB GUI viewer keyed on `class(obj.lastResult)`) both
# belong to the GUI FORM half, not this one -- `build_job` simply returns whichever `CalcXxx`
# result it produced; what a caller (the GUI form's "Run" button, a script, a test) does with
# that value is entirely up to it. The GUI form's own result-`Observable` (see
# `ext/GuiStructureControl.jl`) is this port's replacement for "assign it somewhere a later
# step can find it".
# =========================================================================================

# -----------------------------------------------------------------------------------------
# LayerSpec -- one layer's worth of raw form fields (StructureControlGUI.m's per-tab panel
# state: dirMatFields, rhoField, symmetryDropdown, cijFields, depthField, legOrderField)
# -----------------------------------------------------------------------------------------

"""
    LayerSpec

One layer's worth of raw `StructureControlGUI` form fields -- the per-tab panel state
(`StructureControlGUI.m`'s `obj.panels{k}`) needed to build one `(Material, CalcParams)`
pair. Immutable, plain data -- no widget/handle dependency.

# Fields
- `dirMat` -- 2x3 (or 3x2/3x3) direction-row input, row 1 = surface normal, row 2 = first/
  propagation direction (`StructureControlGUI.m`'s `dirMatFields`, `dirRowLabels =
  {'Surface normal', 'First direction'}`) -- passed straight to [`Material`](@ref)'s own
  `dirMat` keyword, which completes it via [`CompleteDirMat`](@ref) exactly as
  `onDirMatEdited` does; no need to complete it a second time here.
- `rho` -- density [g/cm^3] (`rhoField`).
- `symmClassNo` -- symmetry class, 1-10 (`symmetryDropdown`; see [`SymmClassLibrary`](@ref) --
  this IS "material selection" in this GUI; there is no separate materials-preset library,
  as in MATLAB).
- `cij` -- independent Voigt constants, in [`IndepVoigtPositions`](@ref)`(symmClassNo)` order
  (`cijFields`' enabled upper-triangle cells, read by `readPanelCijVector`). Default: zeros,
  matching every `cijFields` cell's own MATLAB default value of `0` (a blank/all-zero material
  is exactly as constructible -- and exactly as physically meaningless until filled in -- as
  the freshly-opened MATLAB form).
- `relDepth` -- domain depth relative to wavelength; used ONLY when this layer is the
  substrate (position 1 in a [`StructureControlSpec`](@ref)'s `layers`) or is selected as the
  sole HS material (`depthField` on the Substrate tab).
- `realThickness` -- real physical thickness, same units [`ApplyWavelength`](@ref) is called
  with; used ONLY when this layer is a film layer (position > 1) (`depthField` on a film-layer
  tab).
- `legOrder` -- max Legendre order (`legOrderField`).

Construct via the keyword form; `cij`'s default is sized to `symmClassNo`'s own independent
constant count (so passing only `symmClassNo` still gives a validly-shaped, if physically
blank, `LayerSpec`).
"""
struct LayerSpec
    dirMat::Matrix{Float64}
    rho::Float64
    symmClassNo::Int
    cij::Vector{Float64}
    relDepth::Float64
    realThickness::Float64
    legOrder::Int
end

function LayerSpec(; dirMat::AbstractMatrix{<:Real} = [1.0 0.0 0.0; 0.0 0.0 1.0],
                      rho::Real = 1.0,
                      symmClassNo::Integer = 1,
                      cij::AbstractVector{<:Real} = zeros(Float64, length(IndepVoigtPositions(symmClassNo))),
                      relDepth::Real = 20.0,
                      realThickness::Real = 1.0,
                      legOrder::Integer = 50)
    return LayerSpec(Matrix{Float64}(dirMat), Float64(rho), Int(symmClassNo),
                      Float64.(collect(cij)), Float64(relDepth), Float64(realThickness),
                      Int(legOrder))
end

# -----------------------------------------------------------------------------------------
# build_structure / build_hs_structure -- StructureControlGUI.m's syncPanelToStructure
# (:477-498) and buildHSStructure (:1044-1061), collapsed into one-shot builders
# -----------------------------------------------------------------------------------------

"""
    build_structure(layers::AbstractVector{LayerSpec}) -> Structure

Build the FULL layer stack (`layers[1]` = substrate, `layers[2:end]` = film layers, bottom to
top) from raw form fields, via [`ElastConsts`](@ref) -> [`Material`](@ref) -> [`CalcParams`](@ref)
-> [`AddLayer`](@ref) for each layer in order.

Port of the "build a fresh `Structure`" half of `StructureControlGUI.syncPanelToStructure`
(`StructureControlGUI.m`) -- collapsed into one pass (a job-builder builds once per
`Run`, rather than live-editing one `Structure` value incrementally the way the GUI's own
per-keystroke `AddLayer`/`ModifyLayer` calls do; the RESULT is identical either way).

Layer 1 gets `CalcParams(relDepth = layers[1].relDepth, legOrder = ...)`; every other layer
gets `CalcParams(realThickness = layers[k].realThickness, legOrder = ...)` -- mirroring
`syncPanelToStructure`'s own `if p.isSubstrate ... relDepth ... else ... realThickness`
branch (`StructureControlGUI.m`) exactly (MATLAB's film-layer branch passes `relDepth
= []`, i.e. `CalcParams`'s own default, since a film layer's `relDepth` is meaningless until
[`ApplyWavelength`](@ref) derives it later).

Throws `ArgumentError` if `layers` is empty (no substrate to build from), or propagates
whatever [`ElastConsts`](@ref)/[`Material`](@ref)'s own constructors throw for an invalid
`cij`/`dirMat` on any one layer (no extra "layer k:" context is added -- see this file's
header for why: MATLAB's own `catch err; GuiCommon.setStatusError(obj.statusLine,
err.message)` shows the raw underlying message unchanged, so this stays a faithful port
rather than an embellished one).
"""
function build_structure(layers::AbstractVector{LayerSpec})
    isempty(layers) &&
        throw(ArgumentError("StructureControlSpec: at least one layer (the substrate) is required"))
    structure = Structure()
    for (k, layer) in enumerate(layers)
        elC = ElastConsts(layer.cij, layer.symmClassNo)
        mat = Material(; elConsts = elC, rho = layer.rho, dirMat = layer.dirMat)
        cp = k == 1 ? CalcParams(; relDepth = layer.relDepth, legOrder = layer.legOrder) :
                      CalcParams(; realThickness = layer.realThickness, legOrder = layer.legOrder)
        structure = AddLayer(structure, mat, cp)
    end
    return structure
end

"""
    build_hs_structure(layers::AbstractVector{LayerSpec}, hsLayerIndex::Integer) -> Structure

Extract ONE layer (1 = substrate, matching a [`StructureControlSpec`](@ref)'s `layers`
indexing) as a fresh 1-layer half-space `Structure`, named `"HS"`, `nMaxModes = 500`.

Port of `StructureControlGUI.buildHSStructure`.

**KNOWN LIMITATION, carried forward verbatim from the MATLAB source's own class docblock**
(`StructureControlGUI.m`): running a FILM layer (`hsLayerIndex > 1`) standalone through
the HS branch uses a hardcoded `relDepth = 20` truncation depth, not an editable field --
`layer.relDepth` is used only when `hsLayerIndex == 1` (the substrate, which DOES have its own
editable `relDepth`); this is a first-pass simplification inherited from MATLAB, not
introduced by this port.

Throws `ArgumentError` if `hsLayerIndex` is out of range for `layers`.
"""
function build_hs_structure(layers::AbstractVector{LayerSpec}, hsLayerIndex::Integer)
    (1 <= hsLayerIndex <= length(layers)) ||
        throw(ArgumentError("StructureControlSpec: hsLayerIndex = $hsLayerIndex is out of " *
            "range (structure has $(length(layers)) layer(s))"))
    layer = layers[hsLayerIndex]
    elC = ElastConsts(layer.cij, layer.symmClassNo)
    mat = Material(; elConsts = elC, rho = layer.rho, dirMat = layer.dirMat)
    relDepthHS = hsLayerIndex == 1 ? layer.relDepth : 20.0   # see docstring, "KNOWN LIMITATION"
    return Structure(; materIn = mat,
                        calcParamsIn = CalcParams(; relDepth = relDepthHS, legOrder = layer.legOrder),
                        name = "HS", nMaxModes = 500)
end

# -----------------------------------------------------------------------------------------
# StructureControlSpec -- the control-panel half (StructureControlGUI.m's toggle groups +
# three generic calc fields + Advanced panel)
# -----------------------------------------------------------------------------------------

"""
    StructureControlSpec

Every `StructureControlGUI` form input needed to build a `Structure` and dispatch to one of
11 `Calc.*` entry points, as plain, immutable, always-correctly-typed data -- see this file's
header for why each physical quantity gets its OWN named field rather than reusing MATLAB's
three generic relabelled slots.

# Fields
- `layers::Vector{LayerSpec}` -- `layers[1]` = substrate, `layers[2:end]` = film layers.
- `layerMode::Symbol` -- `:hs` (Single layer, half-space) | `:layered` (Full stack) | `:plate`
  (Free plate / Lamb) -- `layerModeGroup`.
- `family::Symbol` -- `:isoWavevec` | `:isoFreq` -- `familyGroup`; consulted ONLY when
  `layerMode === :layered` (mirrors `StructureControlGUI.selectionState`'s own comment: the
  toggle exists but is meaningless for HS/Plate).
- `calcType::Symbol` -- `:single` | `:angSweep` | `:twoAxis` -- `calcTypeGroup`. `:hs` +
  `:twoAxis` has no target class (there is no `CalcWavevecAngSweep2HS`) and
  [`build_job`](@ref) throws `ArgumentError` for that combination, mirroring
  `onLayerModeToggle`'s own `calcType2AxisBtn.Enable = 'off'` for HS.
- `hsLayerIndex::Int` -- which entry of `layers` acts as the HS material (`layerMode === :hs`
  only; 1 = substrate) -- `layerSelectDropdown`.
- `angle`/`angles` -- propagation angle(s) [deg]. `angle` (`calcType === :single`) and
  `angles` (`:angSweep`/`:twoAxis`) are BOTH optional (`nothing` -> each target `CalcXxx`
  constructor's own default: `0.0` / `0:2:180`, applied identically across HS/Layered/Plate --
  see this file's header, "REQUIRED-VS-OPTIONAL FIELDS").
- `wavelength`/`wavelengths` -- layered, `family === :isoWavevec` only. `wavelength`
  (`:single`/`:angSweep`) is REQUIRED -- checked explicitly by [`build_job`](@ref), since it
  feeds [`ApplyWavelength`](@ref) directly, before any `CalcXxx` constructor exists to
  validate it (see this file's header). `wavelengths` (`:twoAxis`) is REQUIRED, validated by
  `CalcWavevecAngSweep2Layered`'s own constructor.
- `vMax` -- layered, `family === :isoWavevec` only; OPTIONAL (`nothing` =
  the band end, 1.05 x the fastest longitudinal limit; it only extends the Green axis).
- `omega0`/`omega0List` -- layered, `family === :isoFreq` only. `omega0`
  (`:single`/`:angSweep`) and `omega0List` (`:twoAxis`) are REQUIRED, validated by the target
  `CalcIsofreqXxx`/`CalcFreqAngSweep2Layered` constructor.
- `kParRange` -- layered, `family === :isoFreq` only (ALL calc types); REQUIRED, validated by
  the target constructor.
- `realWavVec` -- plate, `:single`/`:angSweep` only; OPTIONAL (`nothing` -> `0.0`, the target
  constructor's own default).
- `wavVecs` -- plate, `:twoAxis` only; REQUIRED, validated by `CalcWavevecAngSweep2Plate`'s
  own constructor.
- `realThickness` -- plate, ALL calc types; OPTIONAL for `:single`/`:angSweep` (`nothing` ->
  `Float64[]`), REQUIRED for `:twoAxis` (both contracts enforced by the respective target
  constructor -- `StructureControlGUI.m`'s own comment at `onRunClicked:851-860` notes
  `retentionLevel` is "passed through UNCHECKED" for Plate for the same reason: the target
  class's own validation is the single source of truth).
- `calcMode::Symbol` -- `:full` | `:ritzOnly`; read only for `layerMode in (:hs, :layered)`
  with `family === :isoWavevec` (HS has no family toggle, but DOES have `calcMode` --
  `advCalcMode` is always visible; isoFreq and Plate have no Green's-function scan and simply
  never read this field, exactly as `onRunClicked` never reads `advCalcMode`/`advNVelGreen`
  for Plate, or `calcMode` for isoFreq).
- `retentionLevel::Symbol` -- `:full` | `:saw_only` | `:none`; passed through to every target
  regardless of `layerMode` -- Plate/`:twoAxis` targets reject `:saw_only` via their OWN
  `validateRetentionLevel` call (not pre-checked here, see above).
- `sawClassMethod` -- `nothing` (default: the calc's own per-mode default, 'green' on
  `:full`) | `:green` | `:refined` (half-space only) | `:ritz`; forced to `:ritz` on
  `calcMode = :ritzOnly` (as `StructureControlGUI.classificationSettings`). HS and layered
  isoWavevec targets only (plate / isoFreq results carry no classification).
- `sawThreshold` -- `nothing` (the method's default standard) | `:high` | `:low` | a number
  on the method's score scale (`advStandard` 'method default' / 'high' / 'low' / 'custom').
  There is no `advNVelGreen` field (as in MATLAB): the Green axis is each result's own.

Construct via the keyword form; every field except `layers` defaults to
`StructureControlGUI`'s OWN freshly-opened-window defaults (`layerMode = :hs`, `family =
:isoWavevec`, `calcType = :single`, `hsLayerIndex = 1`, `calcMode = :full`, `retentionLevel =
:full`, `sawClassMethod = sawThreshold = nothing`, every angle/wavelength/etc. field
`nothing`).
"""
struct StructureControlSpec
    layers::Vector{LayerSpec}

    layerMode::Symbol
    family::Symbol
    calcType::Symbol

    hsLayerIndex::Int

    angle::Union{Nothing,Float64}
    angles::Union{Nothing,Vector{Float64}}

    wavelength::Union{Nothing,Float64}
    wavelengths::Union{Nothing,Vector{Float64}}
    vMax::Union{Nothing,Float64}

    omega0::Union{Nothing,Float64}
    omega0List::Union{Nothing,Vector{Float64}}
    kParRange::Union{Nothing,Vector{Float64}}

    realWavVec::Union{Nothing,Float64}
    wavVecs::Union{Nothing,Vector{Float64}}
    realThickness::Union{Nothing,Float64}

    calcMode::Symbol
    retentionLevel::Symbol
    sawClassMethod::Union{Nothing,Symbol}
    sawThreshold::Union{Nothing,Symbol,Float64}
end

_scj_tovec(x) = x === nothing ? nothing : Float64.(vec(collect(x)))
_scj_toscalar(x) = x === nothing ? nothing : Float64(x)

function StructureControlSpec(; layers::AbstractVector{LayerSpec} = LayerSpec[LayerSpec()],
                                 layerMode::Symbol = :hs,
                                 family::Symbol = :isoWavevec,
                                 calcType::Symbol = :single,
                                 hsLayerIndex::Integer = 1,
                                 angle = nothing,
                                 angles = nothing,
                                 wavelength = nothing,
                                 wavelengths = nothing,
                                 vMax = nothing,
                                 omega0 = nothing,
                                 omega0List = nothing,
                                 kParRange = nothing,
                                 realWavVec = nothing,
                                 wavVecs = nothing,
                                 realThickness = nothing,
                                 calcMode::Symbol = :full,
                                 retentionLevel::Symbol = :full,
                                 sawClassMethod = nothing,
                                 sawThreshold = nothing)
    return StructureControlSpec(collect(layers), layerMode, family, calcType, Int(hsLayerIndex),
        _scj_toscalar(angle), _scj_tovec(angles),
        _scj_toscalar(wavelength), _scj_tovec(wavelengths), _scj_toscalar(vMax),
        _scj_toscalar(omega0), _scj_tovec(omega0List), _scj_tovec(kParRange),
        _scj_toscalar(realWavVec), _scj_tovec(wavVecs), _scj_toscalar(realThickness),
        calcMode, retentionLevel,
        sawClassMethod === nothing ? nothing : Symbol(sawClassMethod),
        normaliseThreshold(sawThreshold))
end

# StructureControlGUI.classificationSettings + applyClassification: the run's method
# ('ritz' on ritzOnly) and standard onto a freshly built HS / layered isoWavevec object.
function _scj_setClass(cs, spec::StructureControlSpec)
    cs.sawClassMethod = spec.calcMode === :ritzOnly ? :ritz : spec.sawClassMethod
    cs.sawThreshold = spec.sawThreshold
    return cs
end

# keyword arguments with the `nothing` entries dropped (the target's own default applies)
_scj_kw(; kw...) = (; (k => v for (k, v) in kw if v !== nothing)...)

# -----------------------------------------------------------------------------------------
# build_job -- StructureControlGUI.onRunClicked (StructureControlGUI.m)
# -----------------------------------------------------------------------------------------

function _scj_validateChoice(value::Symbol, allowed::Tuple, name::AbstractString)
    value in allowed ||
        throw(ArgumentError("StructureControlSpec: $name must be $(orList(allowed)); " *
            "got '$value'."))
    return value
end

"""
    build_job(spec::StructureControlSpec)

Validate `spec` and dispatch to the matching `Calc.*` entry point, running it to completion
and returning the result -- the headless equivalent of `StructureControlGUI.onRunClicked`
(`StructureControlGUI.m`), minus the two purely-interactive tails (assigning the
result to a base-workspace variable, opening a viewer window) -- see this file's header for
why neither has a place here.

Returns whichever `CalcXxx` result the (`layerMode`, `family`, `calcType`) combination
selects -- one of `CalcSingleHS`, `CalcAngSweepHS`, `CalcSingleLayered`, `CalcAngSweepLayered`,
`CalcSinglePlate`, `CalcAngSweepPlate`, `CalcIsofreqSingleLayered`,
`CalcIsofreqAngSweepLayered`, `CalcWavevecAngSweep2Layered`, `CalcWavevecAngSweep2Plate`, or
`CalcFreqAngSweep2Layered` (the same 11 types `onRunClicked` itself can produce). The return type is therefore genuinely a big `Union` /
runtime-decided -- exactly as untyped as MATLAB's own `obj.lastResult`; this is a one-shot GUI
job launcher, not a hot loop, so the resulting type instability is immaterial.

Throws `ArgumentError` -- with a clear, specific message -- for: an unrecognised
`layerMode`/`calcType`/`family`; `layerMode === :hs` with `calcType === :twoAxis` (no
`CalcWavevecAngSweep2HS` exists); an out-of-range `hsLayerIndex`; a missing `wavelength` for
`layered`/`isoWavevec` `:single`/`:angSweep`; an empty `layers`; an invalid `dirMat`/`cij` on
any layer (propagated from [`Material`](@ref)/[`ElastConsts`](@ref)); or any of the
required-field/mode-mismatch errors the target `CalcXxx` constructor or `calc`/`calcAll` itself
raises (`missingOmega0`, `missingKParRange`, `missingWavVecs`,
`missingRealThickness`, a `retentionLevel` a target rejects, ...) -- see this file's header,
"REQUIRED-VS-OPTIONAL FIELDS", for why most of these are NOT re-checked here.
"""
function build_job(spec::StructureControlSpec)
    _scj_validateChoice(spec.layerMode, (:hs, :layered, :plate), "layerMode")
    _scj_validateChoice(spec.calcType, (:single, :angSweep, :twoAxis), "calcType")
    spec.layerMode === :layered &&
        _scj_validateChoice(spec.family, (:isoWavevec, :isoFreq), "family")

    if spec.layerMode === :hs
        structureHS = build_hs_structure(spec.layers, spec.hsLayerIndex)
        if spec.calcType === :single
            cs = singleHS(; _scj_kw(structure = structureHS, angle = spec.angle, calcMode = spec.calcMode,
                            retentionLevel = spec.retentionLevel)...)
            return calcAll(_scj_setClass(cs, spec))
        elseif spec.calcType === :angSweep
            cs = angSweepHS(; _scj_kw(structure = structureHS, angles = spec.angles, calcMode = spec.calcMode,
                              retentionLevel = spec.retentionLevel)...)
            return calc(_scj_setClass(cs, spec))
        else   # :twoAxis
            throw(ArgumentError("StructureControlSpec: layerMode = :hs has no 2-axis sweep " *
                "(there is no CalcWavevecAngSweep2HS) -- mirrors StructureControlGUI's own " *
                "disabled '2-axis sweep' radio button when 'Single layer (HS)' is selected."))
        end

    elseif spec.layerMode === :plate
        structure = build_structure(spec.layers)
        if spec.calcType === :single
            cs = singlePlate(; structure = structure, angle = spec.angle,
                                realThickness = spec.realThickness, realWavVec = spec.realWavVec,
                                retentionLevel = spec.retentionLevel)
            return calcAll(cs)
        elseif spec.calcType === :angSweep
            cs = angSweepPlate(; structure = structure, angles = spec.angles,
                                  realThickness = spec.realThickness, realWavVec = spec.realWavVec,
                                  retentionLevel = spec.retentionLevel)
            return calc(cs)
        else   # :twoAxis
            cs = wavevecAngSweep2Plate(; structure = structure, angles = spec.angles,
                                          wavVecs = spec.wavVecs, realThickness = spec.realThickness,
                                          retentionLevel = spec.retentionLevel)
            return calc(cs)
        end

    else   # :layered
        if spec.family === :isoWavevec
            structure = build_structure(spec.layers)
            if spec.calcType === :single
                spec.wavelength === nothing && throw(ArgumentError(
                    "StructureControlSpec: wavelength is required for layerMode = :layered, " *
                    "family = :isoWavevec, calcType = :single"))
                appliedStruct = ApplyWavelength(structure, spec.wavelength)
                cs = singleLayered(; _scj_kw(structure = appliedStruct, angle = spec.angle, vMax = spec.vMax,
                                      calcMode = spec.calcMode, retentionLevel = spec.retentionLevel)...)
                return calcAll(_scj_setClass(cs, spec))
            elseif spec.calcType === :angSweep
                spec.wavelength === nothing && throw(ArgumentError(
                    "StructureControlSpec: wavelength is required for layerMode = :layered, " *
                    "family = :isoWavevec, calcType = :angSweep"))
                appliedStruct = ApplyWavelength(structure, spec.wavelength)
                cs = angSweepLayered(; _scj_kw(structure = appliedStruct, angles = spec.angles, vMax = spec.vMax,
                                        calcMode = spec.calcMode, retentionLevel = spec.retentionLevel)...)
                return calc(_scj_setClass(cs, spec))
            else   # :twoAxis -- wavelength varies per point, ApplyWavelength runs INSIDE calc()
                cs = wavevecAngSweep2Layered(; _scj_kw(structure = structure, angles = spec.angles,
                                                wavelengths = spec.wavelengths, vMax = spec.vMax,
                                                calcMode = spec.calcMode, retentionLevel = spec.retentionLevel)...)
                return calc(_scj_setClass(cs, spec))
            end
        else   # :isoFreq
            structure = build_structure(spec.layers)
            if spec.calcType === :single
                cs = isofreqSingleLayered(; structure = structure, angle = spec.angle,
                                             omega0 = spec.omega0, kParRange = spec.kParRange)
                return calc(cs)
            elseif spec.calcType === :angSweep
                cs = isofreqAngSweepLayered(; structure = structure, angles = spec.angles,
                                               omega0 = spec.omega0, kParRange = spec.kParRange)
                return calc(cs)
            else   # :twoAxis
                cs = freqAngSweep2Layered(; structure = structure, angles = spec.angles,
                                             omega0List = spec.omega0List, kParRange = spec.kParRange)
                return calc(cs)
            end
        end
    end
end
