# =========================================================================================
# GuiStructureControl.jl -- structure editor / job launcher form
#
# Port of the WIDGET half of `guliwerr/gui/StructureControlGUI.m` -- see
# `src/forwardCalcs/StructureControlJob.jl`'s own header for the split rationale (the MATLAB
# class is both a structure editor and a job launcher, so the job-building half is headless)
# and for the headless `StructureControlSpec`/`build_job` this
# file builds on. `include`d directly into the EXISTING `module GuliwerrJMakieExt`, after
# `GuiWidgets.jl`.
#
# `using GuliwerrJ.ForwardCalcs: StructureControlSpec, LayerSpec, build_job` is the one
# explicit qualified import of the headless job-builder names this file needs; everything
# else it touches
# (`IndepVoigtPositions`, the 11 `CalcXxx` result types for [`StructureControlResult`](@ref),
# ...) reaches this file via `GuliwerrJMakieExt`'s existing `using GuliwerrJ`.
# The GUI form itself never calls
# `singleHS`/`angSweepHS`/.../`wavevecAngSweep2Layered` directly -- `build_job` is the ONE
# entry point it dispatches through, exactly mirroring MATLAB's `onRunClicked` calling into
# ONE of 11 `Calc.*` launchers, except here that branching lives entirely on the headless side.
#
# -----------------------------------------------------------------------------------------
# WIDGET-LAYER SIMPLIFICATIONS FROM THE MATLAB SOURCE -- DOCUMENTED, NOT ACCIDENTAL
# -----------------------------------------------------------------------------------------
# 1. ONE REUSABLE LAYER-EDITOR PANEL, NOT ONE TAB PER LAYER. MATLAB's `uitabgroup` gives each
#    layer its OWN persistent set of widgets (`addPanel`/`obj.panels{k}`); Makie has no
#    built-in tab-group widget. Here, ONE set of layer-field widgets (direction matrix,
#    density, symmetry, Cij grid, depth, legOrder) is shared across every layer: a "Layer"
#    `Menu` (mirroring `layerSelectDropdown`) plus Add/Remove buttons manage an internal
#    `Observable{Vector{LayerSpec}}`, and selecting a different layer REPOPULATES the shared
#    widgets from that layer's own stored `LayerSpec` (mirroring `populatePanelFromLayer`),
#    while editing a widget writes back into the CURRENTLY SELECTED layer's slot (mirroring
#    `syncPanelToStructure`). The `layers` Observable is exactly `StructureControlGUI`'s own
#    `structureModel`, one level less baked (a `Vector{LayerSpec}`, not a Structure) -- see
#    `build_job`'s own header for why raw form fields, not a `Structure`, is the right level of
#    "still validatable" to keep in the Observable.
# 2. NAMED CALC-FIELD TEXTBOXES, NOT THREE RELABELLED SLOTS. `StructureControlSpec` already
#    gives every physical quantity (`angle`, `angles`, `wavelength`, `wavelengths`, `vMax`,
#    `omega0`, `omega0List`, `kParRange`, `realWavVec`, `wavVecs`, `realThickness`) its own
#    named field (`StructureControlJob.jl`'s own header explains why); this form mirrors that
#    directly with one always-existing, always-correctly-mapped, ALWAYS VISIBLE `Textbox` per
#    quantity, each with its own label -- there is no per-combination RELABELLING function to
#    keep in sync with `build_job`'s own branches, unlike MATLAB's `refreshCalcFields`. Point 6
#    below records why "shown/hidden per combination" (`refreshCalcFields`'s OTHER role) is not
#    reproduced either.
# 3. ALL 21 UPPER-TRIANGLE Cij CELLS ALWAYS EXIST; NONE ARE VISUALLY DISABLED PER SYMMETRY
#    CLASS. `onSymmetryChanged`'s `Enable = 'on'/'off'` toggle in
#    `StructureControlGUI.m` is a pure display affordance -- [`build_structure`](@ref)/[`build_job`](@ref)
#    only ever READ the cells [`IndepVoigtPositions`](@ref) says are independent for the
#    layer's current `symmClassNo`, exactly as `readPanelCijVector` does; a value typed into a
#    cell that symmetry does not use is simply never read. The MATLAB
#    enable/disable affordance is not reproduced; it is a display nicety, not a
#    correctness gap. The same holds for the symmetry LOCK while a library entry or the
#    Hill average is selected (MATLAB greys the dropdown out, `updateMaterialLock`): a
#    `Menu` cannot be disabled, so a change made there is reverted instead, and the
#    "Hill avg" `Toggle` ignores a switch-on while the symmetry is already isotropic.
# 4. `layerMode`/`family`/`calcType` (MATLAB `uibuttongroup` radio-button stacks) become
#    `Menu`s -- Makie has no built-in radio-button-GROUP widget either, and `Menu` is the
#    substitute for "a small, fixed, mutually-exclusive choice
#    set" (`GuiWidgets.jl`'s own `guiDropdownStepperRow`).
# 5. NO "OPEN VIEWER" / "ASSIGN TO BASE WORKSPACE" -- see `StructureControlJob.jl`'s header.
#    The result [`Observable`](@ref) below (returned as `result` alongside `fig`) IS this
#    port's hand-off point: a caller or viewer watches it instead of `Gui.combinedViewer`
#    being invoked by name from inside `onRunClicked`.
# 6. NO PER-COMBINATION SHOW/HIDE OF CALC FIELDS OR MODE MENUS -- UNAVAILABLE in Makie,
#    not merely skipped for convenience. `refreshCalcFields`'s OTHER role (beyond
#    relabelling, see point 2) is show/hide: MATLAB toggles `.Visible` on whichever of the
#    three generic fields do not apply. Every Makie PLOT OBJECT exposes `.visible`, but
#    Blocks mostly do not (`:visible in fieldnames(Makie.Textbox)`, `Makie.Menu`,
#    `Makie.Button` all `false`; only `Makie.Label` has one) in this Makie version. Toggling only the field LABEL's
#    visibility while leaving its (unhideable) `Textbox`/`Menu` sitting there would be a worse
#    affordance than leaving both always visible with a clear label, so every calc field and
#    every mode-dependent `Menu` (`hsLayerMenu`, `familyMenu`, `calcModeMenu`) stays always
#    visible -- a scoped-down layout simplification, not a
#    correctness gap (`build_job` still only ever reads the fields relevant to the current
#    combination).
# =========================================================================================

using GuliwerrJ.ForwardCalcs: StructureControlSpec, LayerSpec, build_job

# -----------------------------------------------------------------------------------------
# Result hand-off type -- the 11 possible build_job results (one per `Calc.*` entry point),
# spelled out as an explicit Union rather than `Any` so a
# viewer watching `result` can dispatch on `result[]`'s concrete type without
# this file needing to define a shared abstract supertype none of these 11 structs otherwise
# share.
# -----------------------------------------------------------------------------------------

const StructureControlResult = Union{
    CalcSingleHS, CalcAngSweepHS, CalcSingleLayered, CalcAngSweepLayered,
    CalcSinglePlate, CalcAngSweepPlate, CalcIsofreqSingleLayered, CalcIsofreqAngSweepLayered,
    CalcWavevecAngSweep2Layered, CalcWavevecAngSweep2Plate, CalcFreqAngSweep2Layered,
}

# -----------------------------------------------------------------------------------------
# Small pure helpers -- parsing / layer bookkeeping. Kept free of any Makie object so they can
# be exercised directly.
# -----------------------------------------------------------------------------------------

"""
    gscParseNumericField(s) -> Union{Nothing, Float64, Vector{Float64}}

Parse one calc-field `Textbox`'s text as MATLAB's `str2num` would for the scalar/colon-range/
list syntax `StructureControlGUI.parseField` supports:
blank -> `nothing`; `"start:stop"` or `"start:step:stop"` -> the expanded range as a
`Vector{Float64}`; a comma/whitespace-separated list -> a `Vector{Float64}`; a single number ->
a bare `Float64`. Returns `nothing` (never throws) if the text is unparseable by any of these
rules -- the caller decides what a bad/blank field means (required vs. optional), exactly the
division of responsibility [`guiParseRange`](@ref) uses for the axis-range
fields.
"""
function gscParseNumericField(s::AbstractString)
    str = strip(s)
    isempty(str) && return nothing
    parts = split(str, ':')
    if length(parts) in (2, 3)
        nums = tryparse.(Float64, strip.(parts))
        if all(!isnothing, nums)
            return length(nums) == 2 ? collect(nums[1]:nums[2]) : collect(nums[1]:nums[2]:nums[3])
        end
    end
    toks = split(str, r"[,\s]+"; keepempty = false)
    isempty(toks) && return nothing
    nums = tryparse.(Float64, toks)
    all(!isnothing, nums) || return nothing
    return length(nums) == 1 ? nums[1] : nums
end

"""
    gscParseScalar(s) -> Union{Nothing, Float64}

[`gscParseNumericField`](@ref), collapsed to a scalar (first element, if a vector resulted) --
for the SINGLE-valued calc fields (`angle`, `wavelength`, `vMax`, `omega0`, `realWavVec`,
`realThickness`).
"""
function gscParseScalar(s::AbstractString)
    v = gscParseNumericField(s)
    v === nothing && return nothing
    return Float64(v isa AbstractVector ? first(v) : v)
end

"""
    gscParseVector(s) -> Union{Nothing, Vector{Float64}}

[`gscParseNumericField`](@ref), lifted to a `Vector{Float64}` (a bare scalar becomes a
1-element vector) -- for the VECTOR-valued calc fields (`angles`, `wavelengths`, `omega0List`,
`kParRange`, `wavVecs`).
"""
function gscParseVector(s::AbstractString)
    v = gscParseNumericField(s)
    v === nothing && return nothing
    return v isa AbstractVector ? Float64.(v) : [Float64(v)]
end

# 7. THE LIBRARY ENTRY'S NAME IS NOT CARRIED INTO THE LAYER. MATLAB's library dropdown also
#    names the layer's Material after the chosen entry (`panelMaterialName`); `LayerSpec`
#    has no name field, so here the dropdown only fills density, symmetry and Cij. The
#    name is cosmetic (no calculation reads it).

_gscLayerMenuItems(layers::AbstractVector{LayerSpec}) =
    isempty(layers) ? String["Substrate"] :
    vcat("Substrate", ["Layer $i" for i in 1:(length(layers) - 1)])

# -----------------------------------------------------------------------------------------
# Layer editor panel
# -----------------------------------------------------------------------------------------

function _gscBuildLayerPanel(parent, layers::Observable{Vector{LayerSpec}}, populating::Ref{Bool})
    gl = guiGridLayout(parent)

    # ---- layer selector row: Menu + Add/Remove --------------------------------------------
    selRow = guiGridLayout(gl[1, 1])
    Label(selRow[1, 1], "Layer:"; halign = :left, tellwidth = false)
    layerMenu = Menu(selRow[1, 2], options = _gscLayerMenuItems(layers[]))
    addBtn = Button(selRow[1, 3], label = "+ Add layer")
    removeBtn = Button(selRow[1, 4], label = "- Remove")

    # ---- direction matrix (2x3): row 1 = surface normal, row 2 = first direction ----------
    Label(gl[2, 1], "Direction matrix (surface normal / first direction)"; halign = :left, tellwidth = false)
    dirGl = guiGridLayout(gl[3, 1])
    defaultDir = [1.0 0.0 0.0; 0.0 0.0 1.0]
    dirMatBoxes = Matrix{Textbox}(undef, 2, 3)
    for r in 1:2, c in 1:3
        dirMatBoxes[r, c] = Textbox(dirGl[r, c], stored_string = string(defaultDir[r, c]), validator = Float64)
    end

    # ---- materials library: picking an entry fills density, symmetry and Cij (orientation is
    # left as it is); a hand edit of any of those sets it back to "custom" (item 1). While an
    # entry is selected the symmetry is locked; the "Hill avg" toggle (isotropic Hill
    # averages c11 = M, c44 = G; off restores the anisotropic constants) is the only
    # symmetry switch then -- it also works on custom anisotropic Cij ----------------------
    libRows = ListMaterialLibrary()
    libItems = vcat("custom", [e.default ? e.name : "$(e.name) ($(e.variant), T=$(e.T))"
                               for e in libRows])
    libGl = guiGridLayout(gl[4, 1])
    Label(libGl[1, 1], "library material"; halign = :left, tellwidth = false)
    libraryMenu = Menu(libGl[1, 2], options = libItems)
    hillToggle = Toggle(libGl[1, 3], active = false)
    Label(libGl[1, 4], "Hill avg"; halign = :left, tellwidth = false)
    hillSource = Ref{Any}(nothing)   # (symmNo, voigt) the Hill average replaced, for "off"

    # ---- density / symmetry ---------------------------------------------------------------
    rowGl = guiGridLayout(gl[5, 1])
    Label(rowGl[1, 1], "density [g/cm3]"; halign = :left, tellwidth = false)
    rhoBox = Textbox(rowGl[1, 2], stored_string = "1.0", validator = Float64)
    Label(rowGl[1, 3], "symmetry"; halign = :left, tellwidth = false)
    symmItems = ["isotropic", "cubic", "tetragonal I", "tetragonal II", "trigonal I",
                 "trigonal II", "hexagonal", "orthorhombic", "monoclinic", "triclinic"]
    symmetryMenu = Menu(rowGl[1, 4], options = symmItems)

    # ---- Cij grid: full upper triangle of the 6x6 Voigt matrix, always present -------------
    Label(gl[6, 1], "Elasticity: Cij (Voigt matrix) [GPa]"; halign = :left, tellwidth = false)
    cijGl = guiGridLayout(gl[7, 1])
    cijBoxes = Matrix{Union{Nothing,Textbox}}(nothing, 6, 6)
    for r in 1:6, c in r:6
        cijBoxes[r, c] = Textbox(cijGl[r, c], stored_string = "0.0", validator = Float64)
    end

    # ---- depth (relDepth for the substrate, realThickness for a film layer) ---------------
    depthRow = guiGridLayout(gl[8, 1])
    depthLabel = Label(depthRow[1, 1], "relDepth"; halign = :left, tellwidth = false)
    depthBox = Textbox(depthRow[1, 2], stored_string = "20.0", validator = Float64)

    # ---- legOrder ---------------------------------------------------------------------------
    legRow = guiGridLayout(gl[9, 1])
    Label(legRow[1, 1], "legOrder"; halign = :left, tellwidth = false)
    legOrderBox = Textbox(legRow[1, 2], stored_string = "50", validator = Int)

    w = (; layerMenu, addBtn, removeBtn, dirMatBoxes, libraryMenu, hillToggle, rhoBox,
           symmetryMenu, cijBoxes, depthLabel, depthBox, legOrderBox)

    # ---- read the CURRENTLY DISPLAYED widgets into a LayerSpec -----------------------------
    readCurrent = function ()
        dirMat = [something(tryparse(Float64, dirMatBoxes[r, c].stored_string[]), 0.0) for r in 1:2, c in 1:3]
        rho = something(tryparse(Float64, rhoBox.stored_string[]), 1.0)
        symmClassNo = symmetryMenu.i_selected[]
        indepPos = IndepVoigtPositions(symmClassNo)
        cij = [something(tryparse(Float64, cijBoxes[p[1], p[2]].stored_string[]), 0.0) for p in indepPos]
        depthVal = something(tryparse(Float64, depthBox.stored_string[]), 20.0)
        legOrder = something(tryparse(Int, legOrderBox.stored_string[]), 50)
        return LayerSpec(; dirMat = dirMat, rho = rho, symmClassNo = symmClassNo, cij = cij,
                            relDepth = depthVal, realThickness = depthVal, legOrder = legOrder)
    end

    # ---- push the currently-displayed widgets into layers[][selected index] ---------------
    syncToLayers! = function ()
        populating[] && return
        idx = layerMenu.i_selected[]
        (idx === nothing || idx < 1 || idx > length(layers[])) && return
        newLayers = copy(layers[])
        newLayers[idx] = readCurrent()
        layers[] = newLayers
    end

    # ---- populate the widgets FROM layers[][idx], guarded against re-triggering the sync
    # listeners above (mirrors `populatePanelFromLayer`) -------------------------------------
    populateFrom! = function (idx::Integer)
        (idx < 1 || idx > length(layers[])) && return
        layer = layers[][idx]
        populating[] = true
        libraryMenu.i_selected[] = 1   # a stored layer may differ from any library entry
        hillToggle.active[] = false
        hillSource[] = nothing
        for r in 1:2, c in 1:3
            dirMatBoxes[r, c].stored_string[] = string(layer.dirMat[r, c])
        end
        rhoBox.stored_string[] = string(layer.rho)
        symmetryMenu.i_selected[] = layer.symmClassNo
        for r in 1:6, c in r:6
            cijBoxes[r, c].stored_string[] = "0.0"
        end
        for (i, p) in enumerate(IndepVoigtPositions(layer.symmClassNo))
            cijBoxes[p[1], p[2]].stored_string[] = string(layer.cij[i])
        end
        isSubstrate = idx == 1
        depthLabel.text[] = isSubstrate ? "relDepth" : "realThickness [um]"
        depthBox.stored_string[] = string(isSubstrate ? layer.relDepth : layer.realThickness)
        legOrderBox.stored_string[] = string(layer.legOrder)
        populating[] = false
    end

    # ---- write a full ElastConsts into the symmetry menu + Cij grid (callers hold populating)
    showElastConsts! = function (ec)
        symmetryMenu.i_selected[] = ec.symmClassNo
        for r in 1:6, c in r:6
            cijBoxes[r, c].stored_string[] = string(ec.voigt[r, c])
        end
    end

    # ---- Hill average of the displayed Cij, remembering them in hillSource (mirrors
    # applyHillAverage; callers hold populating) --------------------------------------------
    applyHillAverage! = function ()
        spec = readCurrent()
        ec = ElastConsts(spec.cij, spec.symmClassNo)
        hillSource[] = (spec.symmClassNo, ec.voigt)
        mHill, gHill = HillAverages(ec.voigt)
        showElastConsts!(ElastConsts([mHill, gHill], 1))
    end

    # ---- fill density/symmetry/Cij from the chosen library row (mirrors onLibrarySelected) --
    fillFromLibrary! = function (i)
        populating[] && return
        if i === nothing || i <= 1
            # "custom": keep whatever is in the fields, as plain hand-typed values
            populating[] = true
            hillToggle.active[] = false
            populating[] = false
            hillSource[] = nothing
            return
        end
        e = libRows[i - 1]
        populating[] = true
        rhoBox.stored_string[] = string(e.rho)
        showElastConsts!(ElastConsts(e.consts, e.symmNo))
        hillToggle.active[] && applyHillAverage!()
        populating[] = false
        syncToLayers!()
    end

    # ---- Hill avg toggle (mirrors onHillToggled) -------------------------------------------
    hillToggled! = function (on)
        populating[] && return
        populating[] = true
        if on && symmetryMenu.i_selected[] == 1
            hillToggle.active[] = false   # already isotropic: nothing to average
        elseif on
            applyHillAverage!()
        elseif hillSource[] !== nothing
            symmNo, voigt = hillSource[]
            showElastConsts!(ElastConsts(voigt, symmNo))
            hillSource[] = nothing
        end
        populating[] = false
        syncToLayers!()
    end

    # ---- a hand edit of density/symmetry/Cij: back to "custom" and Hill off, the edited
    # values stay (mirrors onMaterialFieldEdited) ----------------------------------------------
    materialEdited! = function ()
        populating[] && return
        populating[] = true
        libraryMenu.i_selected[] = 1
        hillToggle.active[] = false
        populating[] = false
        hillSource[] = nothing
        syncToLayers!()
    end

    # ---- symmetry is locked while a library entry or the Hill average is selected: revert
    # the change (a Menu cannot be disabled; mirrors updateMaterialLock) ----------------------
    shownSymm = Ref(symmetryMenu.i_selected[])
    symmetryEdited! = function (i)
        if !populating[] && (libraryMenu.i_selected[] != 1 || hillToggle.active[])
            populating[] = true
            symmetryMenu.i_selected[] = shownSymm[]
            populating[] = false
            return
        end
        shownSymm[] = i
        materialEdited!()
    end

    on(i -> fillFromLibrary!(i), libraryMenu.i_selected)
    on(a -> hillToggled!(a), hillToggle.active)
    on(_ -> materialEdited!(), rhoBox.stored_string)
    on(i -> symmetryEdited!(i), symmetryMenu.i_selected)
    on(_ -> syncToLayers!(), depthBox.stored_string)
    on(_ -> syncToLayers!(), legOrderBox.stored_string)
    for r in 1:2, c in 1:3
        on(_ -> syncToLayers!(), dirMatBoxes[r, c].stored_string)
    end
    for r in 1:6, c in r:6
        on(_ -> materialEdited!(), cijBoxes[r, c].stored_string)
    end

    on(i -> populateFrom!(i), layerMenu.i_selected)

    on(addBtn.clicks) do _
        newLayers = push!(copy(layers[]), LayerSpec())
        layers[] = newLayers
        layerMenu.options[] = _gscLayerMenuItems(newLayers)
        layerMenu.i_selected[] = length(newLayers)
    end

    on(removeBtn.clicks) do _
        idx = layerMenu.i_selected[]
        if idx === nothing || idx == 1 || length(layers[]) <= 1
            return   # Substrate is never removable -- mirrors onRemoveLayer's own early return
        end
        newLayers = copy(layers[])
        deleteat!(newLayers, idx)
        layers[] = newLayers
        layerMenu.options[] = _gscLayerMenuItems(newLayers)
        layerMenu.i_selected[] = min(idx, length(newLayers))
    end

    return w
end

# -----------------------------------------------------------------------------------------
# Calc-mode toggles + named calc fields + Advanced panel
# -----------------------------------------------------------------------------------------

function _gscBuildControlPanel(parent, layers::Observable{Vector{LayerSpec}})
    gl = guiGridLayout(parent)

    modeRow = guiGridLayout(gl[1, 1])
    Label(modeRow[1, 1], "Target:"; halign = :left, tellwidth = false)
    layerModeMenu = Menu(modeRow[1, 2],
        options = ["Single layer (HS)", "Full stack (Layered)", "Free plate (Lamb)"])
    Label(modeRow[1, 3], "HS layer:"; halign = :left, tellwidth = false)
    hsLayerMenu = Menu(modeRow[1, 4], options = _gscLayerMenuItems(layers[]))
    on(ls -> (hsLayerMenu.options[] = _gscLayerMenuItems(ls)), layers)

    familyRow = guiGridLayout(gl[2, 1])
    Label(familyRow[1, 1], "Family:"; halign = :left, tellwidth = false)
    familyMenu = Menu(familyRow[1, 2], options = ["isoWavevec", "isoFreq"])

    calcTypeRow = guiGridLayout(gl[3, 1])
    Label(calcTypeRow[1, 1], "Calc type:"; halign = :left, tellwidth = false)
    calcTypeMenu = Menu(calcTypeRow[1, 2], options = ["Single", "AngSweep", "2-axis sweep"])

    # ---- named calc fields, one row each, ALWAYS visible (see this file's header, point 6:
    # Makie's `Textbox`/`Menu` Blocks expose no `visible` attribute in this Makie version --
    # unlike `Label`/plot objects, so per-combination
    # show/hide is not reproduced; every field's own label says what it is) ------------------
    fieldsGl = guiGridLayout(gl[4, 1])
    fieldNames = (:angle, :angles, :wavelength, :wavelengths, :vMax,
                  :omega0, :omega0List, :kParRange, :realWavVec, :wavVecs, :realThickness)
    fieldLabels = Dict(
        :angle => "angle [deg]", :angles => "angles [deg]", :wavelength => "wavelength [um]",
        :wavelengths => "wavelengths [um]", :vMax => "vMax", :omega0 => "omega0",
        :omega0List => "omega0 list", :kParRange => "kParRange", :realWavVec => "realWavVec",
        :wavVecs => "wavVecs", :realThickness => "realThickness")
    fieldBoxes = Dict{Symbol,Textbox}()
    fieldLabelObjs = Dict{Symbol,Label}()
    for (i, name) in enumerate(fieldNames)
        rowGl = guiGridLayout(fieldsGl[i, 1])
        fieldLabelObjs[name] = Label(rowGl[1, 1], fieldLabels[name]; halign = :left, tellwidth = false)
        fieldBoxes[name] = Textbox(rowGl[1, 2], stored_string = "")
    end

    advRow1 = guiGridLayout(gl[5, 1])
    Label(advRow1[1, 1], "calcMode:"; halign = :left, tellwidth = false)
    calcModeMenu = Menu(advRow1[1, 2], options = ["full", "ritzOnly"])
    Label(advRow1[1, 3], "retentionLevel:"; halign = :left, tellwidth = false)
    retentionLevelMenu = Menu(advRow1[1, 4], options = ["full", "saw_only", "none"])

    # Classification method / standard (as MATLAB StructureControlGUI): "auto" follows
    # calcMode ('green' on full, 'ritz' on ritzOnly). No nVelGreen: the 2D sweeps use a
    # fixed-step Green axis.
    advRow2 = guiGridLayout(gl[6, 1])
    Label(advRow2[1, 1], "Classification method:"; halign = :left, tellwidth = false)
    sawClassMethodMenu = Menu(advRow2[1, 2], options = ["auto", "green", "refined", "ritz"])
    Label(advRow2[1, 3], "standard:"; halign = :left, tellwidth = false)
    sawThresholdMenu = Menu(advRow2[1, 4], options = ["high", "low"])

    w = (; layerModeMenu, hsLayerMenu, familyMenu, calcTypeMenu, fieldBoxes, fieldLabelObjs,
           calcModeMenu, retentionLevelMenu, sawClassMethodMenu, sawThresholdMenu)

    return w
end

# -----------------------------------------------------------------------------------------
# structureControl -- the top-level entry point (mirrors Gui.structureControl())
# -----------------------------------------------------------------------------------------

"""
    structureControl(; fig = nothing) -> (fig, result, widgets)

Build the structure-editor / job-launcher form -- the Makie/Observable GUI form
half of `guliwerr/gui/StructureControlGUI.m` (`Gui.structureControl()`,
`guliwerr/gui/@Gui/structureControl.m`), built on [`StructureControlSpec`](@ref)/
[`build_job`](@ref) (`src/forwardCalcs/StructureControlJob.jl`) -- see this file's header for
the widget-layer simplifications made porting the tab-per-layer MATLAB editor onto Makie.

`fig`: `nothing` (default) builds a standalone `Figure`; otherwise a `GridPosition` to build
into (via [`guiHostLayout`](@ref)), for embedding in a caller-supplied window.

# Returns a `NamedTuple`
- `fig::Figure` -- the (possibly ancestor) `Figure`.
- `result::Observable{Union{Nothing,StructureControlResult}}` -- `nothing` until the first
  successful `Run`; thereafter the latest [`build_job`](@ref) result. THIS is the hand-off
  point a viewer watches (`on(result) do r; r === nothing || ...; end`) --
  the Julia replacement for MATLAB's `onRunClicked` calling `Gui.combinedViewer(obj.lastResult)`
  by name. Left UNCHANGED by a failed `Run` (a bad spec never clobbers the last good result --
  see `run!`'s own docstring below).
- `statusLabel::Label` -- `StructureControlGUI`'s own status line
  ([`guiSetStatusOk!`](@ref)/[`guiSetStatusError!`](@ref)).
- `layers::Observable{Vector{LayerSpec}}` -- the live layer-stack model
  (`StructureControlGUI.structureModel`, one level less baked -- see this file's header).
- `widgets::NamedTuple` -- every interactive `Block` this form built (layer-panel widgets,
  mode/calc-type `Menu`s, the 11 named calc-field `Textbox`es, the Advanced controls, the Run
  `Button`), for a caller that needs to drive the form programmatically -- set
  `.stored_string[]`/`.i_selected[]` directly rather than simulating a click/keystroke event.
- `run!::Function` -- the EXACT zero-argument function the Run button's own click listener
  calls; exposed so a caller can trigger a run WITHOUT going through a simulated
  Makie click event.
"""
function GuliwerrJ.structureControl(; fig = nothing)
    figOut, layout = guiHostLayout(fig)
    contentGL, controlsGL, _ = guiContentControlsSplit(layout[1, 1])

    layers = Observable(LayerSpec[LayerSpec()])
    populating = Ref(false)

    layerW = _gscBuildLayerPanel(contentGL[1, 1], layers, populating)

    controlW = _gscBuildControlPanel(controlsGL[1, 1], layers)

    runRow = guiGridLayout(controlsGL[2, 1])
    runBtn = Button(runRow[1, 1], label = "Run")
    statusLabel = guiStatusLabel(controlsGL[3, 1])

    result = Observable{Union{Nothing,StructureControlResult}}(nothing)

    run! = function ()
        try
            layerMode = (:hs, :layered, :plate)[controlW.layerModeMenu.i_selected[]]
            family = (:isoWavevec, :isoFreq)[controlW.familyMenu.i_selected[]]
            calcType = (:single, :angSweep, :twoAxis)[controlW.calcTypeMenu.i_selected[]]
            calcMode = (:full, :ritzOnly)[controlW.calcModeMenu.i_selected[]]
            retentionLevel = (:full, :saw_only, :none)[controlW.retentionLevelMenu.i_selected[]]
            sawClassMethod = (nothing, :green, :refined, :ritz)[controlW.sawClassMethodMenu.i_selected[]]
            sawThreshold = sawClassMethod === nothing ? nothing :
                           (:high, :low)[controlW.sawThresholdMenu.i_selected[]]

            fb = controlW.fieldBoxes
            spec = StructureControlSpec(;
                layers = layers[],
                layerMode = layerMode, family = family, calcType = calcType,
                hsLayerIndex = something(controlW.hsLayerMenu.i_selected[], 1),
                angle = gscParseScalar(fb[:angle].stored_string[]),
                angles = gscParseVector(fb[:angles].stored_string[]),
                wavelength = gscParseScalar(fb[:wavelength].stored_string[]),
                wavelengths = gscParseVector(fb[:wavelengths].stored_string[]),
                vMax = gscParseScalar(fb[:vMax].stored_string[]),
                omega0 = gscParseScalar(fb[:omega0].stored_string[]),
                omega0List = gscParseVector(fb[:omega0List].stored_string[]),
                kParRange = gscParseVector(fb[:kParRange].stored_string[]),
                realWavVec = gscParseScalar(fb[:realWavVec].stored_string[]),
                wavVecs = gscParseVector(fb[:wavVecs].stored_string[]),
                realThickness = gscParseScalar(fb[:realThickness].stored_string[]),
                calcMode = calcMode, retentionLevel = retentionLevel,
                sawClassMethod = sawClassMethod, sawThreshold = sawThreshold)

            jobResult = build_job(spec)
            result[] = jobResult
            guiSetStatusOk!(statusLabel, "Done: $(typeof(jobResult))")
        catch err
            msg = err isa ArgumentError ? err.msg : sprint(showerror, err)
            guiSetStatusError!(statusLabel, msg)
        end
        return nothing
    end

    on(_ -> run!(), runBtn.clicks)

    widgets = merge(layerW, controlW, (; runBtn = runBtn))

    return (fig = figOut, result = result, statusLabel = statusLabel, layers = layers,
            widgets = widgets, run! = run!)
end
