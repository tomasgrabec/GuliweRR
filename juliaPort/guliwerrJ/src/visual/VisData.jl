# =========================================================================================
# VisData.jl -- numeric data preparation behind the `@Vis` plotting functions
#
# Port of the NUMERIC half of `guliwerr/visual/@Vis/`: the 8 public data-prep methods
# (`calcDispField2D`, `calcGreenDispField2D`, `calcGreenDispField2DLayered`,
# `extractPlanarCut`, `scaleGreenData`, `calcSubstrateCutoffDepth`,
# `calcModalShapeRowHeights`, `pcolorCenteredGrid`), the 13 numeric `private/` helpers
# (`visDispField2D`, `visDecayCutoffIdx`, `visDepthMagnitude`, `visSphereDirGrid`,
# `visFlattenSignificantModes`, `visClassificationColorVals`, `visPlateRitzModeData`,
# `visPropOrEmpty`, `visValidateChoice`, `visDispersionCurveTypes`,
# `visClassificationTypesAllowed`, `visStandardFrameTensor`) and the two solver-touching
# helpers (`visPlateLimVels`, `visPlateLongLimVel`). Colormap functions live in
# `Colormaps.jl`; the compute halves of the plotting functions (`plotVelSurface` etc.) live in
# `PlotData.jl`. `Visual.jl` `include`s this file directly into `module Visual`; none of
# these MATLAB files is a `classdef`, so no Julia type of the same name exists.
#
# No new result types: every function here returns a plain Julia container (`Vector`,
# `Matrix`, `Tuple`, or a `NamedTuple` for `visFlattenSignificantModes`'s MATLAB
# `struct(...)`-of-parallel-arrays return). `visFlattenSignificantModes` returns a
# `NamedTuple` rather than a one-off struct: it is built once, returned, and never mutated by
# any caller, as with MATLAB's own plain `struct(...)` literal.
#
# `dispField`'s MATLAB SHAPE, PORTED LITERALLY: `Vis.calcDispField2D`/`calcGreenDispField2D`/
# `calcGreenDispField2DLayered` all return a MATLAB `{3 x nLay}` CELL array (`dispField{iDir,
# ilay}` = one `[nMeshDepth x nMeshLateral]` real field). Ported as `Matrix{Matrix{Float64}}`
# sized `(3, nLay)` -- NOT collapsed to a single `Array{Float64,4}` (`(nMeshDepth,
# nMeshLateral, 3, nLay)`), even though today's depth mesh happens to be uniform across
# layers (same `meshSpacing` everywhere) and a 4D array would therefore "work": the MATLAB
# cell's ragged-by-construction shape is what every one of its own unit tests asserts against
# (`tc.verifySize(dispField, [3, nLay])`, `dispField{iDir,ilay}`), and other consumers
# (a GUI's `ModalShapeViewerGUI` panel grid) index it the same `dispField[iDir, ilay]` way. `Matrix{Matrix{Float64}}` is the direct Julia analogue of a
# MATLAB `{M x N}` cell of same-family arrays used throughout this port already (e.g.
# `GreenSolver.k3OutLayered::Vector{Matrix{ComplexF64}}`).
#
# `scaleMode`/`surfaceType`/`substrateCutoff`/`colorBy`/`refLevel`/`curveType` -- MATLAB
# `char` option arguments become `Symbol`s throughout, per this port's established
# convention (`calcMode`, `retentionLevel`, `velAxisMode`, ... -- see `CalcCommon.jl`).
# `curveType`'s six hyphenated MATLAB strings (`'wavelength-freq'`, ...) become
# UNDERSCORE-separated Symbols (`:wavelength_freq`, ...) since a hyphen is not valid inside
# a bare Julia `Symbol` literal (`:wavelength-freq` parses as a subtraction, not one token) --
# `visDispersionCurveTypes` is the single place this translation table is spelled out.
#
# `visValidateChoice` reuses `ForwardCalcs.orList` (imported below) rather than re-deriving
# its own or-list formatter -- one formatter, matching the "message-only" `ArgumentError`
# convention used throughout the port (`validateCalcMode`, `validateRetentionLevel`,
# `validateVelAxisMode`) instead of MATLAB's auto-generated `'MATLAB:<fn>:<id>'`
# identifiers, which have no Julia equivalent and are not reproduced.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# visDispField2D -- private/visDispField2D.m
# -----------------------------------------------------------------------------------------

"""
    visDispField2D(depthProfiles, meshSpacing) -> Matrix{Matrix{Float64}}

Shared depth-profile -> signed 2D (depth x lateral) field assembly.

Port of `guliwerr/visual/@Vis/private/visDispField2D.m`. Shared behind both
[`calcDispField2D`](@ref) (Ritz eigenmode, `nLay` profiles) and [`calcGreenDispField2D`](@ref)
(Green's response, a single half-space profile): they differ only in where the complex depth
profile comes from, not in how it becomes a field.

`depthProfiles[ilay]` is `[nMeshDepth x 3]` COMPLEX, row 1 = deepest (`ti = -1`), row end =
shallowest (`ti = +1`). Returns a `(3, nLay)` matrix of REAL `[nMeshDepth x nMeshLateral]`
fields, normalised by ONE global max across every direction and layer TOGETHER.

**CONTRACT (matches MATLAB, do not "fix"): that max is over POSITIVE values ONLY, not
`max(abs(...))`**, so a negative-dominant mode can fall outside `[-1, 1]`. `lateralFactor` stays genuinely COMPLEX through
the multiplication: the depth profile is complex (Hermitian eigenproblem), so simplifying to
a real cosine would silently drop its imaginary part.
"""
function visDispField2D(depthProfiles::AbstractVector{<:AbstractMatrix{<:Complex}},
                         meshSpacing::Real)
    nMeshLateral = round(Int, 2 / meshSpacing, RoundNearestTiesAway) + 1
    ti = collect(range(-1.0, 1.0, length = nMeshLateral))
    lateralFactor = exp.(1im * pi .* (ti .+ 1))        # [nMeshLateral], complex

    nLay = length(depthProfiles)
    dispField = Matrix{Matrix{Float64}}(undef, 3, nLay)
    maxPos = 0.0
    for ilay in 1:nLay
        for iDir in 1:3
            field2D = real.(depthProfiles[ilay][:, iDir] * transpose(lateralFactor))
            dispField[iDir, ilay] = field2D
            maxPos = max(maxPos, maximum(field2D))
        end
    end

    maxPos = max(maxPos, eps())
    for k in eachindex(dispField)
        dispField[k] = dispField[k] ./ maxPos
    end
    return dispField
end

# -----------------------------------------------------------------------------------------
# visDecayCutoffIdx -- private/visDecayCutoffIdx.m
# -----------------------------------------------------------------------------------------

"""
    visDecayCutoffIdx(dispMag, thresholdFrac) -> Int

First depth index (scanning inward from the surface) at which a displacement profile has
decayed below a fraction of its surface value.

Port of `guliwerr/visual/@Vis/private/visDecayCutoffIdx.m`. `dispMag[end]` is the reference
surface value (floored at `eps()`); scans CONTIGUOUSLY from `dispMag[end]` down to `dispMag[1]`
and returns the first index below `thresholdFrac * reference`, or `1` if the profile never
drops below threshold (keep everything).
"""
function visDecayCutoffIdx(dispMag::AbstractVector{<:Real}, thresholdFrac::Real)
    reference = max(dispMag[end], eps())
    threshold = thresholdFrac * reference

    cutoffIdx = 1
    for k in length(dispMag):-1:1
        if dispMag[k] < threshold
            cutoffIdx = k
            break
        end
    end
    return cutoffIdx
end

# -----------------------------------------------------------------------------------------
# visDepthMagnitude -- private/visDepthMagnitude.m
# -----------------------------------------------------------------------------------------

"""
    visDepthMagnitude(fieldsPerDir) -> Vector{Float64}

Per-depth-row displacement magnitude of a signed 2D field triplet.

Port of `guliwerr/visual/@Vis/private/visDepthMagnitude.m`. `fieldsPerDir` is a 3-element
collection, each `[nMeshDepth x nMeshLateral]` (one [`calcDispField2D`](@ref)/
[`calcGreenDispField2D`](@ref) direction slice, i.e. `dispField[:, ilay]`). Each direction
contributes its own lateral MAXIMUM (not RMS/mean): the lateral factor is a harmonic sweeping
the full `2*pi`, so a row's peak value equals the true depth-profile amplitude at that row
regardless of phase -- see this function's own MATLAB docstring for the derivation.
"""
function visDepthMagnitude(fieldsPerDir::AbstractVector{<:AbstractMatrix{<:Real}})
    nRows = size(fieldsPerDir[1], 1)
    rowMagSq = zeros(Float64, nRows)
    for iDir in 1:3
        rowMagSq .+= vec(maximum(fieldsPerDir[iDir], dims = 2)) .^ 2
    end
    return sqrt.(rowMagSq)
end

# -----------------------------------------------------------------------------------------
# visSphereDirGrid -- private/visSphereDirGrid.m
# -----------------------------------------------------------------------------------------

"""
    visSphereDirGrid(nTheta, nPhi) -> (theta, phi, Nx, Ny, Nz)

Full-sphere `(theta,phi)` direction grid and its unit direction cosines.

Port of `guliwerr/visual/@Vis/private/visSphereDirGrid.m`. Shared by the
[`PlotData`](@ref)-family sphere-surface computations (`velSurfaceData`, `groupVelSurfaceData`,
`youngsSurfaceData`) so all three sweep exactly the same directions -- a precondition for
[`extractPlanarCut`](@ref), which takes the SAME `theta`/`phi` vectors back and would silently
mis-locate crossings against a grid built even slightly differently. `Nx`/`Ny`/`Nz` are
`[nPhi x nTheta]`, matching MATLAB's `meshgrid(theta, phi)` orientation exactly (`theta` varies
along dim 2, `phi` along dim 1).
"""
function visSphereDirGrid(nTheta::Integer, nPhi::Integer)
    theta = collect(range(0.0, pi, length = nTheta))
    phi = collect(range(0.0, 2pi, length = nPhi))

    Nx = [sin(th) * cos(ph) for ph in phi, th in theta]
    Ny = [sin(th) * sin(ph) for ph in phi, th in theta]
    Nz = repeat(transpose(cos.(theta)), length(phi), 1)   # cos(theta), broadcast over phi rows
    return theta, phi, Nx, Ny, Nz
end

# -----------------------------------------------------------------------------------------
# visFlattenSignificantModes -- private/visFlattenSignificantModes.m
# -----------------------------------------------------------------------------------------

"""
    visFlattenSignificantModes(smCell, angles, typeAllowed) -> NamedTuple

Flatten per-angle `SignificantMode` vectors into flat point arrays.

Port of `guliwerr/visual/@Vis/private/visFlattenSignificantModes.m`.
`smCell[ia]` is a `Vector{SignificantMode}` (`ModeClassification`'s per-angle output);
`angles[ia]` the angle of that entry (`[NaN]` for a single-angle caller with no angle axis);
`typeAllowed` a collection of `Symbol` type names to keep
([`visClassificationTypesAllowed`](@ref)).

Returns a `NamedTuple` of flat vectors, one entry per kept point: `angle`, `vel`, `modeNo`,
`modeVel`, `type`, `method`, `grade` (`Vector{Symbol}`), `greenPeakVel`, `score`, `scoreDb` -- the
fields the classification plots read; nothing is recomputed here.
"""
function visFlattenSignificantModes(smCell::AbstractVector{<:AbstractVector{SignificantMode}},
                                     angles::AbstractVector{<:Real},
                                     typeAllowed)
    angle = Float64[]; vel = Float64[]; modeNo = Int[]; modeVel = Float64[]
    type = Symbol[]; method = Symbol[]; grade = Symbol[]
    greenPeakVel = Float64[]; score = Float64[]; scoreDb = Float64[]

    for ia in eachindex(smCell)
        for m in smCell[ia]
            m.type in typeAllowed || continue
            push!(angle, Float64(angles[ia])); push!(vel, m.vel); push!(modeNo, m.modeNo)
            push!(modeVel, m.modeVel)
            push!(type, m.type); push!(method, m.method); push!(grade, m.grade)
            push!(greenPeakVel, m.greenPeakVel); push!(score, m.score); push!(scoreDb, m.scoreDb)
        end
    end
    return (angle = angle, vel = vel, modeNo = modeNo, modeVel = modeVel, type = type,
            method = method, grade = grade, greenPeakVel = greenPeakVel, score = score, scoreDb = scoreDb)
end

# -----------------------------------------------------------------------------------------
# visClassificationColorVals -- private/visClassificationColorVals.m
# -----------------------------------------------------------------------------------------

"""
    visClassificationColorVals(colorBy, vels, scoreDb, greenPeakVel) -> (cVals, cLabel)

Continuous colour values for a classification plot's `:score`/`:confidence` colouring.

Port of `guliwerr/visual/@Vis/private/visClassificationColorVals.m`.
`:score` -- the entry's own `scoreDb` (the score its classification method thresholded, in
dB; `NaN` on limiting entries; comparable only within one method). `:confidence` --
`abs(modeVel - greenPeakVel)`, the bound Ritz mode's distance from the Green peak it was
matched to: 0 is a perfect match, `NaN` means "never matched". The callers pass the entries'
`modeVel` as `vels`, not `vel`: a Green entry's `vel` IS its `greenPeakVel` (MATLAB
likewise).
"""
function visClassificationColorVals(colorBy::Symbol, vels::AbstractVector{<:Real},
                                     scoreDb::AbstractVector{<:Real},
                                     greenPeakVel::AbstractVector{<:Real})
    if colorBy === :score
        return Float64.(scoreDb), "score [dB]"
    else
        return abs.(vels .- greenPeakVel), "|v_{mode} - v_{Green peak}| [km/s]"
    end
end

# -----------------------------------------------------------------------------------------
# visPlateRitzModeData -- private/visPlateRitzModeData.m
# -----------------------------------------------------------------------------------------

"""
    visPlateRitzModeData(ritz::RitzSolverLamb) -> (vels, types, conf)

Per-mode velocity, mode type and type confidence from a Lamb Ritz solver.

Port of `guliwerr/visual/@Vis/private/visPlateRitzModeData.m`. Works at any
`retentionLevel`, since `pruneEigVecs` only ever touches `eigVecs_con`, never these fields.

MATLAB's own docstring flags its `reshape(..., [], 1)` calls as load-bearing: `modesType`/
`modeTypeConfidence` are MATLAB ROW vectors while `vels` is a COLUMN, and naive indexing
silently preserved row orientation, corrupting a caller's vertical concatenation. **That
distinction does not exist in Julia** -- every field here is already a plain `Vector` (no
row/column orientation), so the three returned vectors need no reshape; this is a genuine
simplification the port gets for free, not a missing feature.
"""
function visPlateRitzModeData(ritz::RitzSolverLamb)
    nM = ritz.nCompModes
    return ritz.vels[1:nM], ritz.modesType[1:nM], ritz.modeTypeConfidence[1:nM]
end

# -----------------------------------------------------------------------------------------
# visPropOrEmpty -- private/visPropOrEmpty.m
# -----------------------------------------------------------------------------------------

"""
    visPropOrEmpty(obj, names) -> value or `nothing`

Read the first field in a preference list that `obj` actually has.

Port of `guliwerr/visual/@Vis/private/visPropOrEmpty.m`. Several result properties carry the
same physical meaning under different names on the HS and Layered result types (`longVelMap`
vs `subLongVelMap`, `longLimVel` vs `subLongLimVel`), and one may simply be absent. `names` is
a collection of `Symbol` field names, most-preferred first. MATLAB's `isprop`/dynamic `.()`
duck-typing becomes `hasfield`/`getfield` on the concrete struct type -- the direct Julia
analogue, since every candidate here is a plain immutable `struct`, not a MATLAB handle
object with runtime-variable properties.

Returns `nothing` (MATLAB's `[]`) if no candidate name is present. The two families' names
mean genuinely different physical quantities (HS's system-wide longitudinal limiting velocity
vs. Layered's SUBSTRATE-only one) -- interchangeable only for the display purposes that use
this helper (a velocity-axis cutoff), not in general; a caller that needs to KNOW which family
it has should branch on `hasfield` itself rather than call this.
"""
function visPropOrEmpty(obj, names)
    T = typeof(obj)
    for name in names
        hasfield(T, name) && return getfield(obj, name)
    end
    return nothing
end

# -----------------------------------------------------------------------------------------
# visValidateChoice -- private/visValidateChoice.m
# -----------------------------------------------------------------------------------------

"""
    visValidateChoice(value, allowed, errPrefix, paramName)

Reject a `Symbol`-option argument that isn't one of the allowed values.

Port of `guliwerr/visual/@Vis/private/visValidateChoice.m`. Nearly every data-prep/plotting
function in `visual/` takes one or more `Symbol` options (`:log`/`:linear`,
`:velocity`/`:slowness`, `:score`/`:confidence`, ...); this is the one shared check, reusing
[`ForwardCalcs.orList`](@ref) for the "or"-joined message (the same formatter
`validateCalcMode`/`validateRetentionLevel`/`validateVelAxisMode` already use) rather than
re-deriving one. MATLAB's own error identifier convention
(`'Vis:<fnName>:bad<ParamName>'`) is folded into the `ArgumentError` message text, per this
port's message-only convention -- there is no Julia analogue of a MATLAB error identifier to
preserve.
"""
function visValidateChoice(value, allowed, errPrefix::AbstractString, paramName::AbstractString)
    value in allowed && return nothing
    badName = isempty(paramName) ? paramName : uppercase(paramName[1]) * paramName[2:end]
    throw(ArgumentError("Vis:$errPrefix:bad$badName -- $paramName must be $(orList(allowed)); " *
        "got '$value'."))
end

# -----------------------------------------------------------------------------------------
# visDispersionCurveTypes -- private/visDispersionCurveTypes.m
# -----------------------------------------------------------------------------------------

"""
    visDispersionCurveTypes() -> NTuple{6,Symbol}

The six `curveType` names the dispersion-curve data-prep functions accept.

Port of `guliwerr/visual/@Vis/private/visDispersionCurveTypes.m`. MATLAB's hyphenated
strings (`'wavelength-freq'`, ...) become underscore-separated `Symbol`s -- a hyphen cannot
appear in a bare Julia `Symbol` literal (`:wavelength-freq` parses as subtraction), so this is
the one place the exact translation table is spelled out:

| MATLAB | Julia |
|---|---|
| `'wavelength-freq'` | `:wavelength_freq` |
| `'wavelength-vphase'` | `:wavelength_vphase` |
| `'k-freq'` | `:k_freq` |
| `'k-vphase'` | `:k_vphase` |
| `'freq-vphase'` | `:freq_vphase` |
| `'freq-vgroup'` | `:freq_vgroup` |
"""
visDispersionCurveTypes() = (:wavelength_freq, :wavelength_vphase, :k_freq, :k_vphase,
                              :freq_vphase, :freq_vgroup)

# -----------------------------------------------------------------------------------------
# visClassificationTypesAllowed -- private/visClassificationTypesAllowed.m
# -----------------------------------------------------------------------------------------

"""
    visClassificationTypesAllowed(refLevel::Symbol) -> NTuple{N,Symbol}

Which classified mode types a `refLevel` tier displays.

Port of `guliwerr/visual/@Vis/private/visClassificationTypesAllowed.m`. `refLevel` is
`:sawOnly` | `:sawLimit` | `:full`. `:other`-type points (unclassified/rejected candidates)
are never shown at any tier; `:limiting` points (Ritz modes sitting AT a bulk-wave limiting
velocity) appear at `:full` only.
"""
function visClassificationTypesAllowed(refLevel::Symbol)
    if refLevel === :full
        return (:Rayleigh, :pseudoRayleigh, :SAW, :limiting)
    else
        return (:Rayleigh, :pseudoRayleigh, :SAW)
    end
end

# -----------------------------------------------------------------------------------------
# visStandardFrameTensor -- private/visStandardFrameTensor.m
# -----------------------------------------------------------------------------------------

"""
    visStandardFrameTensor(mater::Material) -> (elTensor, rho)

A `Material`'s elastic tensor in the project's standard frame, plus density.

Port of `guliwerr/visual/@Vis/private/visStandardFrameTensor.m`. `mater.dirMat` already
carries whatever orientation the caller wants evaluated (a propagation-angle rotation, if
any, is the caller's business). Standard frame, traced directly
(`GetStandardElastTensor(mater, 0)`: `GetStandardDirMat -> CalcRotMat -> RotateTensor`): axis 1 = `dirMat`'s in-plane
propagation reference (`dirMat` row 2), axis 2 = the remaining in-plane direction (right-hand
completion, row 3), axis 3 = the SURFACE NORMAL (row 1) -- the same frame Ritz/Christoffel/
GreenSolver use throughout the project.
"""
function visStandardFrameTensor(mater::Material)
    return GetStandardElastTensor(mater, 0.0), mater.rho
end

# -----------------------------------------------------------------------------------------
# visPlateLimVels / visPlateLongLimVel -- private/visPlateLimVels.m, visPlateLongLimVel.m
# -----------------------------------------------------------------------------------------

"""
    visPlateLimVels(structure::Structure, layerIdx::Integer, angleDeg) -> (longLimVel, shearLimVel)

Christoffel longitudinal + slowest-shear limiting velocities for one layer (surface) of a
Plate structure.

Port of `guliwerr/visual/@Vis/private/visPlateLimVels.m`. `angleDeg` is `Union{Real,Nothing}`:
a real value rotates `structure.mater[layerIdx]` to that propagation angle first; `nothing`
(MATLAB's `[]`) assumes `structure`'s `dirMat` ALREADY carries the propagation-angle rotation.
This distinction is LOAD-BEARING, not a convenience -- `CalcWavevecAngSweep2Plate`'s and
`CalcAngSweepPlate`'s `.structure` are kept UNROTATED (so those Plate-family callers
must rotate here), whereas `CalcSinglePlate.structure` is already rotated for its single
propagation angle and rotating again would double-apply it.

**NOT deduplicated against `CalcCommon.limitingAndBulk`'s own longitudinal-root selection**:
`limitingAndBulk` takes `unique(...)` then the LOWEST-velocity root (defensive against a
spurious numerically-distinct duplicate); this function instead takes `vLim_sA[1][1]` --
literally the FIRST root `findRoots` returns for the longitudinal branch (in ascending-theta
SCAN order, not sorted by velocity). Ported exactly as MATLAB has it, not "fixed" to match
`limitingAndBulk`'s own convention -- flagged as a genuine, pre-existing
MATLAB fragility (in the ordinary case there is exactly one longitudinal root over the coarse
scan, so the two conventions agree; they would disagree only if the coarse scan spuriously
found more than one).
"""
function visPlateLimVels(structure::Structure, layerIdx::Integer, angleDeg::Union{Real,Nothing})
    mater = structure.mater[layerIdx]
    if angleDeg !== nothing
        mater = SetDirMat(mater, GetRotDirMat(mater, angleDeg))
    end

    elTensor, rho = visStandardFrameTensor(mater)
    vLim_sA = findRoots(elTensor, rho)

    longLimVel = vLim_sA[1][1].vel
    shearLimVel = minimum(r.vel for r in Iterators.flatten((vLim_sA[2], vLim_sA[3])))
    return longLimVel, shearLimVel
end

"""
    visPlateLongLimVel(structure::Structure, angleDeg) -> Float64

Longitudinal limiting velocity of a Plate structure, taken over BOTH free surfaces.

Port of `guliwerr/visual/@Vis/private/visPlateLongLimVel.m`. `longLimVel` -- max over the
bottom (layer 1) and, when the stack has more than one layer, the top (layer `nlay`) surface,
so a genuinely asymmetric stack's cutoff is never set by the slower surface alone (which
would silently crop real modes belonging to the faster one). `angleDeg` -- see
[`visPlateLimVels`](@ref) for the `nothing`-vs-real distinction.
"""
function visPlateLongLimVel(structure::Structure, angleDeg::Union{Real,Nothing})
    longLimVel = visPlateLimVels(structure, 1, angleDeg)[1]
    if nlay(structure) > 1
        longLimVel = max(longLimVel, visPlateLimVels(structure, nlay(structure), angleDeg)[1])
    end
    return longLimVel
end

# -----------------------------------------------------------------------------------------
# calcDispField2D -- calcDispField2D.m
# -----------------------------------------------------------------------------------------

"""
    calcDispField2D(ritzSolver::RitzSolver, modeNo::Integer, meshSpacing::Real = 0.01)
        -> Matrix{Matrix{Float64}}

Signed 2D (depth x lateral) displacement field for one Ritz mode.

Port of `guliwerr/visual/@Vis/calcDispField2D.m`. `ritzSolver` needs `eigVecs_con`/
`nullConstrMat` populated for the requested mode (`retentionLevel = :full`; there is no
recompute-on-demand for a pruned solver). `modeNo` selects which mode (column of
`eigVecs_con`) to evaluate.

Returns a `(3, nLay)` matrix, `dispField[iDir, ilay]` = one real signed 2D field
`[nMeshDepth x nMeshLateral]`, normalised by one global max across every direction and layer
together -- see [`visDispField2D`](@ref) for the exact (asymmetric, positive-values-only)
normalisation contract. Built on [`calcDisplacement`](@ref); this function only slices
that output into per-layer depth profiles.

**Phase**: `calcDisplacement` is phase-fixed exactly as MATLAB's (`phaseFix`, largest
in-plane coefficient real; C2z structures give real u1/u2 and imaginary u3), so this output
matches MATLAB's `Vis.calcDispField2D` element-wise on every non-degenerate mode.
"""
function calcDispField2D(ritzSolver::RitzSolver, modeNo::Integer, meshSpacing::Real = 0.01)
    displDir = calcDisplacement(ritzSolver, modeNo, meshSpacing)   # [depthPts*nLay x modeNo x 3]
    depthProfile = displDir[:, end, :]                              # mode modeNo only, [depthPts*nLay x 3]

    nLay = nlay(ritzSolver.structure)
    displSize = round(Int, 2 / meshSpacing, RoundNearestTiesAway) + 1

    depthProfiles = Vector{Matrix{ComplexF64}}(undef, nLay)
    for ilay in 1:nLay
        depthProfiles[ilay] = depthProfile[((ilay - 1) * displSize + 1):(ilay * displSize), :]
    end

    return visDispField2D(depthProfiles, meshSpacing)
end

# -----------------------------------------------------------------------------------------
# calcGreenDispField2D -- calcGreenDispField2D.m
# -----------------------------------------------------------------------------------------

"""
    calcGreenDispField2D(greenObj::GreenSolver, forceDir::Integer, meshSpacing::Real = 0.01)
        -> Matrix{Matrix{Float64}}

Signed 2D (depth x lateral) displacement field for one Green's-function response (one
velocity, one force direction) -- the Green's-response analogue of [`calcDispField2D`](@ref).

Port of `guliwerr/visual/@Vis/calcGreenDispField2D.m`. `greenObj` needs `.dispField`
populated: call `calc(GreenSolver, structure, angle, velSelected; fullOutput = true, ...)`
first (MATLAB's `nargout > 1`). Full output is half-space only. `velSelected` should be a
SCALAR velocity (a 1-element `velRange`) -- `dispField` is built over every point of whatever
`velRange` is given, so a single value keeps the expensive full-partial-wave solve to the one
velocity being displayed.

Returns a `(3, 1)` matrix (half-space only, so `nLay = 1` always) -- the same shape/
normalisation convention as `calcDispField2D`'s `(3, nLay)` output. Row 1 = deepest, row end =
the free surface, RE-ORIENTED from `GreenSolver`'s own storage (index 1 = surface) to match
`calcDispField2D`'s row convention.

No physical depth-axis reconstruction is attempted here: this mirrors MATLAB's own choice
(only a normalised row-index fraction is ever needed downstream) even though `GreenSolver`
in this port additionally stores `.depthAxis` (a deliberate divergence from MATLAB) -- that field is
simply not read by this function, exactly as MATLAB's `xx3` axis is never read by its own
`calcGreenDispField2D.m` either.
"""
function calcGreenDispField2D(greenObj::GreenSolver, forceDir::Integer, meshSpacing::Real = 0.01)
    isempty(greenObj.dispField) &&
        throw(ArgumentError("calcGreenDispField2D:noFullOutput -- greenObj.dispField is " *
            "empty -- call calc(GreenSolver, ...; fullOutput = true) (full output) first; " *
            "full output is available for half-space structures only."))
    (1 <= forceDir <= 3) ||
        throw(ArgumentError("calcGreenDispField2D:notLessEqual -- forceDir must be an " *
            "integer in 1:3, got $forceDir"))

    depthProfile = greenObj.dispField[:, :, forceDir, 1]   # [nxx3 x 3], complex
    depthProfile = reverse(depthProfile, dims = 1)          # row 1 = deepest now

    # one "layer" only (half-space), so the shared helper's global normalisation reduces to
    # "across the 3 directions"
    return visDispField2D([depthProfile], meshSpacing)
end

# -----------------------------------------------------------------------------------------
# calcGreenDispField2DLayered -- calcGreenDispField2DLayered.m
# -----------------------------------------------------------------------------------------

"""
    calcGreenDispField2DLayered(greenObj::GreenSolver, structure::Structure,
                                 forceDir::Integer, meshSpacing::Real = 0.01)
        -> Matrix{Matrix{Float64}}

Signed 2D (depth x lateral) displacement field for a LAYERED Green's-function response (one
velocity, one force direction) -- the layered counterpart of [`calcGreenDispField2D`](@ref).

Port of `guliwerr/visual/@Vis/calcGreenDispField2DLayered.m`. Built from `GreenSolver.calc`'s
`k3OutLayered`/`polarizationULayered`/`partialAmplitudeLayered`/`hScaleOut` (`fullOutput =
true`, `nlay(structure) > 1`) instead of the half-space-only `.dispField`. `structure` must be
the SAME `Structure` passed to `calc(GreenSolver, ...)` -- `GreenSolver` stores no `Structure`
copy (deliberately, to avoid growing the value type), so the caller must supply it again here;
only `.calcParams[:].relDepth` is read, Structure order (index 1 = substrate, index `nlay` =
outermost film).

Returns a `(3, nlay(structure))` matrix, Structure order, normalised by ONE global max across
every direction and layer together -- the SAME convention as `calcDispField2D`'s own output,
including that the substrate's FULL `relDepth[1]` extent is returned uncropped (no
`substrateCutoff` applied here, matching `calcDispField2D`).

**Reconstruction.** Every layer's solved amplitude vector is referenced to the SAME global
`x3` coordinate as the substrate (`x3 = 0` at the substrate interface, `x3 = -H_total` at the
free surface, Every 2002 convention):
`u(x3,ci,cj) = sum_n A[n,cj] * U[ci,n] * exp(i*k3[n]*x3)`, `ci` = displacement component,
`cj` = force direction. Mapped to `RitzSolver`'s own `ti in [-1,1]` convention (`ti = -1`
deep, `ti = +1` shallow) linearly per layer.

**Substrate display extent**: `GreenSolver`'s substrate is analytically semi-infinite --
`structure.calcParams[1].relDepth` has no meaning to `GreenSolver`'s own physics, but is
reused here PURELY as a finite display extent, exactly as if the substrate were another film
of that thickness.

**Evanescent-overflow limitation** (inherited from MATLAB, not fixed here): `exp(i*k3*x3)`
for a single partial wave can overflow to `Inf` for a thick film/deep substrate extent even
though the true summed field is finite -- see `docs/reference/architecture/green.md`, "Conditioning".
"""
function calcGreenDispField2DLayered(greenObj::GreenSolver, structure::Structure,
                                      forceDir::Integer, meshSpacing::Real = 0.01)
    isempty(greenObj.k3OutLayered) &&
        throw(ArgumentError("calcGreenDispField2DLayered:noFullOutput -- " *
            "greenObj.k3OutLayered is empty -- call calc(GreenSolver, ...; fullOutput = true) " *
            "on a LAYERED structure (nlay >= 2) first."))
    nlay(structure) == length(greenObj.k3OutLayered) ||
        throw(ArgumentError("calcGreenDispField2DLayered:structureMismatch -- " *
            "structure has $(nlay(structure)) layer(s) but greenObj.k3OutLayered has " *
            "$(length(greenObj.k3OutLayered)) -- pass the SAME Structure used to build greenObj."))
    (1 <= forceDir <= 3) ||
        throw(ArgumentError("calcGreenDispField2DLayered:notLessEqual -- forceDir must be " *
            "an integer in 1:3, got $forceDir"))

    nLay = nlay(structure)
    hScale = greenObj.hScaleOut
    nMesh = round(Int, 2 / meshSpacing, RoundNearestTiesAway) + 1
    ti = collect(range(-1.0, 1.0, length = nMesh))

    # film thickness per layer, Structure order (index 1 = substrate, left at 0/unused;
    # index 2..nLay = films), matching calc.m's own h = 2*pi*relDepth convention
    filmThickness = zeros(Float64, nLay)
    for vIdx in 2:nLay
        filmThickness[vIdx] = 2 * pi * structure.calcParams[vIdx].relDepth
    end

    depthProfiles = Vector{Matrix{ComplexF64}}(undef, nLay)
    for vIdx in 1:nLay
        k3 = greenObj.k3OutLayered[vIdx][:, 1]              # [nWaves], first (only) velocity
        U = greenObj.polarizationULayered[vIdx][:, :, 1]    # [3 x nWaves]
        A = greenObj.partialAmplitudeLayered[vIdx][:, :, 1] # [nWaves x 3]
        nWaves = length(k3)

        x3 = if vIdx == 1
            # substrate: ti=-1 (deep) -> x3 = full extent, ti=+1 (shallow) -> x3 = 0 (the
            # substrate/first-film interface)
            substrateExtent = 2 * pi * structure.calcParams[1].relDepth * hScale
            substrateExtent .* (1 .- ti) ./ 2
        else
            # film: depth below the substrate interface (x3=0) of this film's lower (deep,
            # ti=-1) and upper (shallow, ti=+1) boundaries, summed from the films strictly
            # BELOW it (Structure order: index 2 = the bottom film, touching the substrate)
            dBelow = sum(@view filmThickness[2:(vIdx - 1)]) * hScale    # 0 when vIdx==2
            dThis = filmThickness[vIdx] * hScale
            @. -dBelow - (ti + 1) / 2 * dThis
        end

        profile = zeros(ComplexF64, nMesh, 3)
        for ci in 1:3, n in 1:nWaves
            coef = A[n, forceDir] * U[ci, n]
            @. profile[:, ci] += coef * exp(im * k3[n] * x3)
        end
        depthProfiles[vIdx] = profile
    end

    return visDispField2D(depthProfiles, meshSpacing)
end

# -----------------------------------------------------------------------------------------
# extractPlanarCut -- extractPlanarCut.m
# -----------------------------------------------------------------------------------------

"""
    _interp2_linear(x, y, V, xq, yq) -> Vector{Float64}

Bilinear interpolation of `V[i,j]` (defined on the grid `x[j]`/`y[i]`, MATLAB's
`interp2(x,y,V,...)` orientation) at the query points `(xq[k], yq[k])`.

Internal helper for [`extractPlanarCut`](@ref); not a port of a named MATLAB function (MATLAB
calls the builtin `interp2` directly). Every query point here is guaranteed to lie strictly
inside `[minimum(x),maximum(x)] x [minimum(y),maximum(y)]` -- a `Contour.jl` vertex sits on a
grid cell's edge by construction -- so no extrapolation behaviour is needed or provided.
"""
function _interp2_linear(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                          V::AbstractMatrix{<:Real},
                          xq::AbstractVector{<:Real}, yq::AbstractVector{<:Real})
    itp = interpolate((collect(Float64, y), collect(Float64, x)), Float64.(V), Gridded(Linear()))
    return [itp(yq[k], xq[k]) for k in eachindex(xq)]
end

"""
    extractPlanarCut(X, Y, Z, theta, phi, planeNormal, inPlaneRef = nothing)
        -> (uCoord, wCoord, ptsXYZ)

Extract the planar cut of a `(theta,phi)`-parametrised 3D surface.

Port of `guliwerr/visual/@Vis/extractPlanarCut.m`. `X`,`Y`,`Z` are `[nPhi x nTheta]` surface
coordinates for ONE sheet/wave-type (call once per sheet for a multi-sheet surface, e.g.
`PlotData.velSurfaceData`'s/`PlotData.groupVelSurfaceData`'s 3 wave types); `theta`,`phi` are
the parameter grids used to build them (the SAME `nTheta`/`nPhi` -- see
[`visSphereDirGrid`](@ref)); `planeNormal` is a 3-vector, normal of the cutting plane through
the origin; `inPlaneRef` (optional) is one of the two in-plane basis directions to project
onto (default: an arbitrary vector perpendicular to `planeNormal`, chosen as whichever
standard-basis axis is LEAST aligned with it).

Returns `uCoord`,`wCoord` (2D in-plane coordinates of the cut curve, possibly several
disconnected segments, `NaN`-separated -- matching MATLAB's own plot/contour convention for
broken curves) and `ptsXYZ` (`N x 3`, the same points in original 3D coordinates, same
`NaN`-row separation).

**`contourc` -> `Contour.jl`.** `Contour.jl`'s
`levels`/`lines`/`coordinates` API replaces MATLAB's hand-parsed packed contour matrix. The
port was checked against MATLAB on a fixture with one closed contour, one contour
crossing a phi-wrap boundary, and two disjoint segments, with this result: **segment grouping and vertex LOCATIONS match MATLAB's
`contourc` to machine precision, ORDER-INSENSITIVELY** (matched by bounding box, each segment's
point sequence coincides with MATLAB's own up to a direction reversal, and for a CLOSED loop
additionally up to a cyclic rotation of its starting point) -- but **NOT order-sensitively**:
`Contour.jl` may trace a closed loop starting at a different point and in the opposite
direction, and may trace a boundary-touching open segment in the opposite direction, from
MATLAB's `contourc`. Nothing in the port (or MATLAB's own test suite) reads segment order
or per-segment traversal direction -- both are drawn as a plain broken polyline downstream,
which renders identically regardless -- so this is not a functional gap, only a fact to keep in mind when
comparing against MATLAB (match the point sets order-insensitively, not index by index).

`Contour.jl` expects `z[i,j] <-> x[i],y[j]`; MATLAB's `F` here is `[nPhi x nTheta]` (`phi`
rows, `theta` columns, matching `X`/`Y`/`Z`'s own layout) exactly as `contourc(theta,phi,F,
[0 0])` expects it (`theta` as x, `phi` as y) -- so `Contour.contour` is called on
`permutedims(F)` (`[nTheta x nPhi]`), not `F` itself.
"""
function extractPlanarCut(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real},
                           Z::AbstractMatrix{<:Real},
                           theta::AbstractVector{<:Real}, phi::AbstractVector{<:Real},
                           planeNormal::AbstractVector{<:Real},
                           inPlaneRef::Union{AbstractVector{<:Real},Nothing} = nothing)
    n = planeNormal ./ norm(planeNormal)
    if inPlaneRef === nothing
        # pick whichever standard-basis axis is LEAST aligned with n, to avoid a
        # near-degenerate (n, candidate) pair
        candidates = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        proj = [abs(dot(candidates[:, j], n)) for j in 1:3]
        idx = argmin(proj)
        v = candidates[:, idx]
    else
        v = collect(Float64, inPlaneRef)
    end
    u = v .- dot(v, n) .* n
    u = u ./ norm(u)
    w = cross(n, u)

    F = X .* n[1] .+ Y .* n[2] .+ Z .* n[3]   # signed distance from the plane

    if all(>(0), F) || all(<(0), F)
        # surface never crosses the plane -- Contour.jl would find nothing; return empty
        # rather than let a degenerate level through
        return Float64[], Float64[], zeros(Float64, 0, 3)
    end

    z = permutedims(F)   # [nTheta x nPhi], matching Contour.jl's z[i,j] <-> x[i],y[j]
    cl = Contour.contour(theta, phi, z, 0.0)

    uCoord = Float64[]; wCoord = Float64[]
    ptsXYZ = zeros(Float64, 0, 3)
    for c in Contour.lines(cl)
        segTheta, segPhi = Contour.coordinates(c)

        segX = _interp2_linear(theta, phi, X, segTheta, segPhi)
        segY = _interp2_linear(theta, phi, Y, segTheta, segPhi)
        segZ = _interp2_linear(theta, phi, Z, segTheta, segPhi)

        segU = segX .* u[1] .+ segY .* u[2] .+ segZ .* u[3]
        segW = segX .* w[1] .+ segY .* w[2] .+ segZ .* w[3]

        append!(uCoord, segU); push!(uCoord, NaN)
        append!(wCoord, segW); push!(wCoord, NaN)
        ptsXYZ = vcat(ptsXYZ, hcat(segX, segY, segZ), fill(NaN, 1, 3))
    end
    return uCoord, wCoord, ptsXYZ
end

# -----------------------------------------------------------------------------------------
# scaleGreenData -- scaleGreenData.m
# -----------------------------------------------------------------------------------------

"""
    _nanmax(A) -> Float64

`maximum`, excluding `NaN` entries (MATLAB's `max(A(:))` default behaviour) -- `Inf`/`-Inf`
pass through unchanged; `NaN` if every entry is `NaN`. Internal helper for
[`scaleGreenData`](@ref); not itself a port of a named MATLAB function.

Julia's own `maximum`/`max` PROPAGATE `NaN` (`max(NaN, x) === NaN`, unlike MATLAB's two-array
`max`, which treats `NaN` as absent whenever the other operand is not also `NaN`) -- this is
guide-flagged MATLAB-vs-Julia divergence #1 relevant to this file; see `scaleGreenData`'s own
docstring for why it does not otherwise matter here.
"""
function _nanmax(A)
    m = -Inf
    seen = false
    for x in A
        if !isnan(x)
            seen = true
            m = max(m, x)
        end
    end
    return seen ? m : NaN
end

"""
    scaleGreenData(absG, scaleMode::Symbol, componentLabel::AbstractString) -> (plotData, dataLabel)

Apply the shared Green's-magnitude display scaling, with its axis label.

Port of `guliwerr/visual/@Vis/scaleGreenData.m`. `absG` -- `|G|` data of
any shape (curve, map, or a partially-`NaN` interpolated map), any real element type (the
Calc results store `Float32`; the scaling is done in `Float64`, as MATLAB's `double(absG)`);
`scaleMode` -- `:log` | `:linear` (anything but `:log` is linear, as MATLAB);
`componentLabel` -- REQUIRED TeX label of the quantity, used verbatim: `"G_{13}^{eff}"`,
`"G_{33}"`, `"S"` (the Ritz spectrum), ... (see [`greenComponentData`](@ref)).

`:log` plots `log10(|G| / max|G|)`, zeros floored to `eps()`; `NaN` in -> `NaN` out, and
`NaN`s are excluded from the maximum (`_nanmax`).
"""
function scaleGreenData(absG::AbstractArray{<:Real}, scaleMode::Symbol,
                        componentLabel::AbstractString)
    A = Array{Float64}(absG)
    if scaleMode === :log
        normVal = A ./ _nanmax(A)
        plotData = log10.(max.(normVal, eps()))
        plotData[isnan.(normVal)] .= NaN
        return plotData, "log_{10}(|$componentLabel| / max|$componentLabel|)"
    else
        return A, "|$componentLabel|"
    end
end

# -----------------------------------------------------------------------------------------
# calcSubstrateCutoffDepth -- calcSubstrateCutoffDepth.m
# -----------------------------------------------------------------------------------------

"""
    calcSubstrateCutoffDepth(substrateCutoff::Symbol, cutoffParam, relDepths,
                              substrateDispMag = nothing) -> Float64

Effective substrate `relDepth` to DISPLAY for `plotModalShape`'s axis layout.

Port of `guliwerr/visual/@Vis/calcSubstrateCutoffDepth.m`. `substrateCutoff` -- `:none` |
`:decay` | `:fixedMultiple`; `cutoffParam` -- `nothing` (MATLAB's `[]`) or a threshold
fraction (`:decay`, default `0.05`) / multiple (`:fixedMultiple`, default `4`), ignored for
`:none`; `relDepths` -- every layer's `relDepth`, layer 1 = substrate, layer `nlay` =
outermost film (nearest the free surface); `substrateDispMag` -- required only for `:decay`:
displacement magnitude across the SUBSTRATE's own depth mesh, index 1 = deepest (`ti = -1`),
index end = the substrate/first-film interface (`ti = +1`, or the free surface itself when
`nlay(structure) == 1`).

Returns the effective substrate `relDepth` to use for axis-height allocation
(`<= relDepths[1]`).
"""
function calcSubstrateCutoffDepth(substrateCutoff::Symbol, cutoffParam::Union{Real,Nothing},
                                   relDepths::Union{Real,AbstractVector{<:Real}},
                                   substrateDispMag::Union{AbstractVector{<:Real},Nothing} = nothing)
    rd = relDepths isa Real ? [Float64(relDepths)] : Float64.(relDepths)
    all(>(0), rd) ||
        throw(ArgumentError("calcSubstrateCutoffDepth: relDepths must be positive"))
    nLay = length(rd)

    if substrateCutoff === :none
        return rd[1]

    elseif substrateCutoff === :fixedMultiple
        nLay >= 2 ||
            throw(ArgumentError("calcSubstrateCutoffDepth:noFilmLayers -- 'fixedMultiple' " *
                "needs at least one film layer to reference -- use :none for a plain half-space."))
        cp = cutoffParam === nothing ? 4.0 : Float64(cutoffParam)
        cp > 0 ||
            throw(ArgumentError("calcSubstrateCutoffDepth: cutoffParam must be positive"))
        filmStackDepth = sum(@view rd[2:end])
        return min(rd[1], cp * filmStackDepth)

    elseif substrateCutoff === :decay
        nLay >= 2 ||
            throw(ArgumentError("calcSubstrateCutoffDepth:noFilmLayers -- 'decay' needs at " *
                "least one film layer to reference -- use :none for a plain half-space."))
        substrateDispMag === nothing &&
            throw(ArgumentError("calcSubstrateCutoffDepth:missingDispMag -- 'decay' requires " *
                "substrateDispMag (the computed mode's displacement magnitude across the " *
                "substrate depth mesh)."))
        cp = cutoffParam === nothing ? 0.05 : Float64(cutoffParam)
        (0 < cp < 1) ||
            throw(ArgumentError("calcSubstrateCutoffDepth: cutoffParam must be in (0,1)"))
        all(>=(0), substrateDispMag) ||
            throw(ArgumentError("calcSubstrateCutoffDepth: substrateDispMag must be nonnegative"))

        # scan from the interface (index end, ti=+1) down toward the deepest point (index 1,
        # ti=-1)
        nMesh = length(substrateDispMag)
        cutoffIdx = visDecayCutoffIdx(substrateDispMag, cp)

        ti = collect(range(-1.0, 1.0, length = nMesh))
        return rd[1] * (1 - ti[cutoffIdx]) / 2

    else
        throw(ArgumentError("calcSubstrateCutoffDepth:badMode -- substrateCutoff must be " *
            ":none, :decay, or :fixedMultiple (got '$substrateCutoff')."))
    end
end

# -----------------------------------------------------------------------------------------
# calcModalShapeRowHeights -- calcModalShapeRowHeights.m
# -----------------------------------------------------------------------------------------

"""
    calcModalShapeRowHeights(relDepths, minFrac::Real = 0.3) -> Float64 or Vector{Float64}

Non-proportional per-layer height fractions for the modal-shape depth stack.

Port of `guliwerr/visual/@Vis/calcModalShapeRowHeights.m`. `relDepths` -- every layer's
`relDepth`, uncropped, layer 1 = substrate, layer `nlay` = outermost film; `minFrac` -- floor
fraction of the total displayed height every layer is guaranteed, before the remainder is
split proportionally to true `relDepth` (default `0.3`, auto-clamped down if `nLay * minFrac`
would reach/exceed 1). Returns a value that sums to 1 (a bare `Float64` when `relDepths` is a
scalar/1-layer input, matching MATLAB's own `1x1` return for that case; a `Vector{Float64}`
otherwise).

Every layer is guaranteed a floor share up front so a thin film stays visible next to a much
deeper substrate.

**Deviation from MATLAB (not a functional gap)**: MATLAB returns the OUTPUT in the same
orientation (`1xN` vs `Nx1`) as the INPUT -- a row/column distinction that plain Julia
`Vector`s do not have, so it is not mirrored here.
"""
function calcModalShapeRowHeights(relDepths::Union{Real,AbstractVector{<:Real}}, minFrac::Real = 0.3)
    rd = relDepths isa Real ? [Float64(relDepths)] : Float64.(relDepths)
    all(>(0), rd) ||
        throw(ArgumentError("calcModalShapeRowHeights: relDepths must be positive"))
    (0 < minFrac < 1) ||
        throw(ArgumentError("calcModalShapeRowHeights: minFrac must be in (0,1)"))
    nLay = length(rd)

    mf = Float64(minFrac)
    if nLay * mf >= 1
        mf = 0.9 / nLay
    end

    remainPool = 1 - nLay * mf
    frac = mf .+ remainPool .* (rd ./ sum(rd))
    return relDepths isa Real ? frac[1] : frac
end

# -----------------------------------------------------------------------------------------
# pcolorCenteredGrid -- pcolorCenteredGrid.m
# -----------------------------------------------------------------------------------------

"""
    _local_edges(v::AbstractVector{<:Real}) -> Vector{Float64}

Cell-edge coordinates for one axis: interior edges at midpoints between consecutive samples;
the two open ends linearly extrapolated by their nearest interior half-gap. A single-sample
(degenerate) axis gets an arbitrary unit-width edge pair. Internal helper for
[`pcolorCenteredGrid`](@ref) (MATLAB's local function `localEdges`).
"""
function _local_edges(v::AbstractVector{<:Real})
    if length(v) == 1
        return [v[1] - 0.5, v[1] + 0.5]
    end
    mid = (v[1:(end - 1)] .+ v[2:end]) ./ 2
    return vcat(2 * v[1] - mid[1], mid, 2 * v[end] - mid[end])
end

"""
    _is_monotonic(v::AbstractVector{<:Real}) -> Bool

Non-strict monotonicity: tolerates repeated-equal values, rejects only an actual direction
reversal. Internal helper for [`pcolorCenteredGrid`](@ref) (MATLAB's local function
`isMonotonic`).
"""
function _is_monotonic(v::AbstractVector{<:Real})
    length(v) <= 1 && return true
    d = diff(v)
    return all(>=(0), d) || all(<=(0), d)
end

"""
    pcolorCenteredGrid(x, y, C) -> (xEdges, yEdges, CPadded)

Cell-edge coordinates + padded `C` so a `pcolor`-style draw centres each `C[i,j]` cell on
`(x[j], y[i])` rather than treating it as the cell's lower-left corner (`docs/pitfalls.md`
NUM3).

Port of `guliwerr/visual/@Vis/pcolorCenteredGrid.m`. `x`,`y` -- sample coordinates (a bare
scalar accepted for the degenerate single-sample-axis case, matching MATLAB's own `x(:)'`
handling of a `1x1` input); `C` -- `[length(y) x length(x)]` colour data, `C[i,j]` belongs to
`(x[j], y[i])`. Returns `xEdges` (`length(x)+1`), `yEdges` (`length(y)+1`), and `CPadded`
(`(length(y)+1) x (length(x)+1)`, last row/column duplicated from the last real row/column --
exists only to satisfy `pcolor`'s size requirement, never drawn by `'flat'` shading).

**CONTRACT: `x`/`y` need only be MONOTONIC, not increasing** --
`CalcIsofreqAngSweepLayered`'s `kParRange = omega0 ./ vAxis` is a naturally DECREASING axis,
and this function's own non-strict monotonicity check (`_is_monotonic`, tolerating repeated
values, rejecting only a direction reversal) must keep working on it; ported and tested for
both the increasing and the decreasing case.

**MATLAB had two overlapping implementations of this workaround** (`Vis.pcolorCenteredGrid`
itself, and `private/visPcolorCentered.m`, which bundled the edge construction with an actual
`pcolor(...)`+`shading('flat')` call). Only this one is ported here -- `visPcolorCentered` is
a thin PLOTTING wrapper (it needs a graphics handle) with no numeric content of its own
beyond calling this function, so it is left to the plotting code rather than duplicating
the `CPadded` construction.

Julia has no row/column-vector distinction, so MATLAB's own
`testOutputOrientationMatchesPcolorConvention` (asserting `xEdges` stays `1x(N+1)`, `yEdges`
`(M+1)x1`) does not translate -- both are returned as plain `Vector{Float64}` here.
"""
function pcolorCenteredGrid(x::Union{Real,AbstractVector{<:Real}},
                             y::Union{Real,AbstractVector{<:Real}},
                             C::AbstractMatrix{<:Real})
    xv = x isa Real ? [Float64(x)] : Float64.(vec(x))
    yv = y isa Real ? [Float64(y)] : Float64.(vec(y))

    all(isfinite, xv) ||
        throw(ArgumentError("pcolorCenteredGrid: x must be finite"))
    all(isfinite, yv) ||
        throw(ArgumentError("pcolorCenteredGrid: y must be finite"))
    size(C) == (length(yv), length(xv)) ||
        throw(ArgumentError("pcolorCenteredGrid:incorrectSize -- C must be " *
            "$(length(yv))x$(length(xv)), got $(size(C))"))
    _is_monotonic(xv) ||
        throw(ArgumentError("pcolorCenteredGrid:xNotMonotonic -- x must be monotonic " *
            "(increasing or decreasing; repeated values are fine)."))
    _is_monotonic(yv) ||
        throw(ArgumentError("pcolorCenteredGrid:yNotMonotonic -- y must be monotonic " *
            "(increasing or decreasing; repeated values are fine)."))

    xEdges = _local_edges(xv)
    yEdges = _local_edges(yv)
    CPadded = [C C[:, end]; permutedims(C[end, :]) C[end, end]]
    return xEdges, yEdges, CPadded
end
