# =========================================================================================
# PlotData.jl -- numeric halves of the plotting functions
#
# Several `@Vis` plotting functions in MATLAB run real solves before the actual
# `surf`/`plot`/`pcolor`/axis-label call. This file holds the NUMERIC half of each of them --
# everything up to (but not including) the drawing call -- so the Makie extension receives
# finished arrays and draws them, with no Makie dependency here. Two such call sites
# (`calcDispField2D.m`, `private/visPlateLimVels.m`) live in `VisData.jl` because they are
# public, directly reusable helpers; this file calls back into those rather than duplicating
# them. Covered here: `plotDispersionCurves`, `plotRitzSpectrumMap`/`Curve`,
# `plotVelSurface`, `plotGroupVelSurface`, `plotYoungsSurface` and
# `plotVelCurve`/`plotYoungsCurve`.
#
# `visDispersionCurveAxes` is not ported here: it is pure arithmetic (a `curveType` ->
# `(xAll,yAll,xLabelStr,yLabelStr)` lookup) with no graphics call, but it only decides which
# array goes on which axis, so it belongs to the drawing side. `dispersionCurvesData` below
# stops at the same point `plotDispersionCurves.m` does right before calling it: it returns
# the assembled `wavelength`/`k`/`vel`/`velGroup`/`freq` arrays and the branch-highlight
# mask, and the "which of these goes on which axis for this curveType" reshuffle is done by
# the plotting function.
#
# No new result types: every function below returns a `Tuple` or `NamedTuple` of plain
# arrays/scalars, exactly like `VisData.jl`.
#
# Imports: none beyond `Visual.jl`'s own `using`s -- every solver/forward-calc entry point
# this file calls (`phaseVels`, `phaseVelScan`, `groupVelVectorsForDir`, `Tensor2Youngs`,
# `GetRotMat`, `RitzScore.spectrum`, `DispersionBranches.track`) is already imported there.
# =========================================================================================

# -----------------------------------------------------------------------------------------
# velSurfaceData -- compute half of plotVelSurface.m
# -----------------------------------------------------------------------------------------

"""
    velSurfaceData(mater::Material, surfaceType::Symbol = :velocity,
                    nTheta::Integer = 61, nPhi::Integer = 121)
        -> (theta, phi, Xs, Ys, Zs, V)

Compute half of `guliwerr/visual/@Vis/plotVelSurface.m`: 3D directional bulk-wave velocity
(or slowness) surfaces.

`surfaceType` -- `:velocity` (default) | `:slowness` (`= 1/velocity`, computed slot-wise, per
wave-type branch, not re-sorted after inversion). `nTheta`/`nPhi` -- polar/azimuthal mesh
resolution.

Returns the SAME `theta`/`phi` [`visSphereDirGrid`](@ref) built the sphere from (hand these,
unchanged, to [`extractPlanarCut`](@ref) for a planar cut), `Xs`/`Ys`/`Zs`
(`[nPhi x nTheta x 3]`, one sheet per wave type `[qL, qT1, qT2]` -- each sheet is the unit
wave normal scaled by that wave type's own scalar radius: unlike
[`groupVelSurfaceData`](@ref), where the plotted point is a genuine vector that need not lie
along the wave normal), and `V` (`[nPhi x nTheta x 3]`, the raw velocity/slowness value per
sheet -- what a Makie `surf(...)`'s colour argument would use).

`mater`'s tensor is taken in the project's standard frame (`visStandardFrameTensor`): axis 1
= `dirMat`'s in-plane propagation reference (row 2), axis 2 = the remaining in-plane
direction (row 3), axis 3 = the SURFACE NORMAL (row 1).
"""
function velSurfaceData(mater::Material, surfaceType::Symbol = :velocity,
                         nTheta::Integer = 61, nPhi::Integer = 121)
    visValidateChoice(surfaceType, (:velocity, :slowness), "plotVelSurface", "surfaceType")

    elTensor, rho = visStandardFrameTensor(mater)
    theta, phi, X, Y, Z = visSphereDirGrid(nTheta, nPhi)

    V = zeros(Float64, nPhi, nTheta, 3)
    for it in 1:nTheta, ip in 1:nPhi
        dirVec = [X[ip, it], Y[ip, it], Z[ip, it]]
        V[ip, it, :] = phaseVels(elTensor, dirVec, rho)
    end

    if surfaceType === :slowness
        V = 1 ./ V
    end

    Xs = X .* V
    Ys = Y .* V
    Zs = Z .* V
    return theta, phi, Xs, Ys, Zs, V
end

# -----------------------------------------------------------------------------------------
# groupVelSurfaceData -- compute half of plotGroupVelSurface.m
# -----------------------------------------------------------------------------------------

"""
    groupVelSurfaceData(mater::Material, nTheta::Integer = 61, nPhi::Integer = 121)
        -> (theta, phi, X, Y, Z, mag)

Compute half of `guliwerr/visual/@Vis/plotGroupVelSurface.m`: 3D group-velocity ("ray")
surface -- the direction of energy/signal propagation, which in general deviates from the
wave-normal direction.

Same full-sphere sweep and frame convention as [`velSurfaceData`](@ref). Unlike
`velSurfaceData`, the surface point at `(theta,phi)` is NOT `dirVec` scaled by a scalar: `X`,
`Y`, `Z` (`[nPhi x nTheta x 3]`) are the actual group-velocity VECTOR components
`groupVelVectorsForDir` returns, which can point in a genuinely different direction from the
wave normal that produced them -- a self-intersecting, "multiply-folded" surface is real
physics (Every et al. 2013, Wave Motion 50, Sec. 2.3), not a computation bug. `mag` is the
group-velocity MAGNITUDE per sheet (what a Makie `surf(...)`'s colour argument would use --
the surface radius is NOT itself that magnitude, since the plotted point need not lie along
the wave normal).

For a 2D planar cut (e.g. the sagittal plane), pass `theta`/`phi` and one sheet's `X[:,:,k]`/
`Y[:,:,k]`/`Z[:,:,k]` to [`extractPlanarCut`](@ref). Computing group velocity over the WHOLE
sphere and then filtering is necessary here, unlike `velSurfaceData`, where restricting the
sweep to one plane from the start already gives the right answer.
"""
function groupVelSurfaceData(mater::Material, nTheta::Integer = 61, nPhi::Integer = 121)
    elTensor, rho = visStandardFrameTensor(mater)
    theta, phi, Nx, Ny, Nz = visSphereDirGrid(nTheta, nPhi)

    X = zeros(Float64, nPhi, nTheta, 3)
    Y = zeros(Float64, nPhi, nTheta, 3)
    Z = zeros(Float64, nPhi, nTheta, 3)
    for it in 1:nTheta, ip in 1:nPhi
        dirVec = [Nx[ip, it], Ny[ip, it], Nz[ip, it]]
        vG, _, _ = groupVelVectorsForDir(elTensor, dirVec, rho)   # [component, wavetype]
        X[ip, it, :] = vG[1, :]
        Y[ip, it, :] = vG[2, :]
        Z[ip, it, :] = vG[3, :]
    end

    mag = sqrt.(X .^ 2 .+ Y .^ 2 .+ Z .^ 2)
    return theta, phi, X, Y, Z, mag
end

# -----------------------------------------------------------------------------------------
# youngsSurfaceData -- compute half of plotYoungsSurface.m
# -----------------------------------------------------------------------------------------

"""
    youngsSurfaceData(mater::Material, nTheta::Integer = 61, nPhi::Integer = 121)
        -> (theta, phi, Xs, Ys, Zs, E)

Compute half of `guliwerr/visual/@Vis/plotYoungsSurface.m`: 3D directional Young's-modulus
surface.

Single-panel analogue of [`velSurfaceData`](@ref): Young's modulus is one scalar field (not
three wave-type sheets), so `Xs`/`Ys`/`Zs`/`E` are `[nPhi x nTheta]` (no third, wave-type
dimension) -- the plotted point is the unit direction scaled by `E(direction)`, same
"radius = scalar modulus" convention as `velSurfaceData`'s velocity sheets.

Per-point orthonormal completion `(n1,n2,n3)` is built directly here (a cheap 2-line
`cross`/`normalize`), rather than via `ElastConsts.RightHandFrame`'s MATLAB `disp()` side
effect that would otherwise fire on nearly every one of the `nTheta*nPhi` points -- porting
that side effect would be pure noise in a headless numeric function; MATLAB's own comment at
`plotYoungsSurface.m` already flags this as deliberate, not an oversight this port
introduces.
"""
function youngsSurfaceData(mater::Material, nTheta::Integer = 61, nPhi::Integer = 121)
    elTensor, _ = visStandardFrameTensor(mater)
    theta, phi, X, Y, Z = visSphereDirGrid(nTheta, nPhi)

    E = zeros(Float64, nPhi, nTheta)
    for it in 1:nTheta, ip in 1:nPhi
        n1 = [X[ip, it], Y[ip, it], Z[ip, it]]

        ref = [0.0, 0.0, 1.0]
        if abs(dot(n1, ref)) > 0.9   # n1 too close to ref -- switch axis
            ref = [1.0, 0.0, 0.0]
        end
        n2 = cross(ref, n1); n2 = n2 ./ norm(n2)
        n3 = cross(n1, n2)

        E[ip, it] = Tensor2Youngs(elTensor, vcat(n1', n2', n3'))
    end

    Xs = X .* E
    Ys = Y .* E
    Zs = Z .* E
    return theta, phi, Xs, Ys, Zs, E
end

# -----------------------------------------------------------------------------------------
# velCurveData -- compute half of plotVelCurve.m
# -----------------------------------------------------------------------------------------

"""
    velCurveData(mater::Material, angleDeg = range(0.0, 360.0, length = 361))
        -> (angleDeg, vels)

Compute half of `guliwerr/visual/@Vis/plotVelCurve.m`: bulk-wave phase velocity vs. in-plane
angle, at the fixed surface normal in `mater.dirMat`.

The 2D-curve analogue of [`velSurfaceData`](@ref): that sweeps the full sphere, this sweeps
one great circle -- the in-plane rotation of `mater.dirMat` about its fixed surface normal
(row 1), the same angle convention every other angular sweep in this project uses
(`phaseVelScan`, `CalcAngSweepHS`/`Layered`). `vels` is `[3 x nAngles]`, `[vL;vT1;vT2]` --
directly `phaseVelScan(mater, angleDeg)` (Solvers/Christoffel), already exactly this
compute; this function exists so `PlotData` has one call per MATLAB plotting function to
mirror, and so the Makie extension imports only from `Visual`, never reaching back into
`Solvers` itself for a plot it did not otherwise need.
"""
function velCurveData(mater::Material, angleDeg = range(0.0, 360.0, length = 361))
    return angleDeg, phaseVelScan(mater, angleDeg)
end

# -----------------------------------------------------------------------------------------
# youngsCurveData -- compute half of plotYoungsCurve.m
# -----------------------------------------------------------------------------------------

"""
    youngsCurveData(mater::Material, angleDeg = range(0.0, 360.0, length = 361))
        -> (angleDeg, E)

Compute half of `guliwerr/visual/@Vis/plotYoungsCurve.m`: directional Young's modulus vs.
in-plane angle, at the fixed surface normal in `mater.dirMat`.

The 2D-curve analogue of [`youngsSurfaceData`](@ref) (which sweeps the full sphere): this
sweeps one great circle, the same convention [`velCurveData`](@ref)/`phaseVelScan` use for
bulk waves. `GetRotMat(mater, angleDeg[ia])` returns the `3x2` `[propagation|normal]` format
`Tensor2Youngs` expects directly -- `ElastConsts.GetYoungsModulus`'s MATLAB instance-method
wrapper has no Julia port (it adds nothing over calling the underlying free function
directly on the struct's own fields), so this calls
`Tensor2Youngs(mater.elConsts.tensor, rotMat)` in its place.
"""
function youngsCurveData(mater::Material, angleDeg = range(0.0, 360.0, length = 361))
    E = zeros(Float64, length(angleDeg))
    for (ia, a) in pairs(angleDeg)
        rotMat = GetRotMat(mater, a)
        E[ia] = Tensor2Youngs(mater.elConsts.tensor, rotMat)
    end
    return angleDeg, E
end

# -----------------------------------------------------------------------------------------
# ritzSpectrumMapData / ritzSpectrumCurveData -- compute halves of plotRitzSpectrumMap.m /
# plotRitzSpectrumCurve.m (the Ritz spectrum replaced RitzSurfaceResponse)
# -----------------------------------------------------------------------------------------

"""
    ritzSpectrumMapData(sweepResult, scaleMode::Symbol = :log)
        -> (angles, velAxis, plotData, cbLabel, absMap)

Compute half of `@Vis/plotRitzSpectrumMap.m`: the Ritz spectrum `|S(v)|` as an angle x
velocity map, rebuilt from every angle's stored modal numerators (no solve; works for
`:ritzOnly` results too, where it is the default background). One common axis: the union of
the angles' `ritzScoreVelAxis`, `vLo:RitzScore.VEL_STEP:vHi` (MATLAB colon semantics);
scaled by [`scaleGreenData`](@ref) with label `"S"`. Throws
`Vis:plotRitzSpectrumMap:noRitzSpectrum` when the sweep stores no numerators.
"""
function ritzSpectrumMapData(sweepResult::Union{CalcAngSweepHS,CalcAngSweepLayered},
                             scaleMode::Symbol = :log)
    visValidateChoice(scaleMode, (:log, :linear), "plotRitzSpectrumMap", "scaleMode")
    visHasRitzSpectrum(sweepResult) || throw(ArgumentError("Vis:plotRitzSpectrumMap:" *
        "noRitzSpectrum -- This sweep stores no Ritz modal numerators (calc() stores them)."))
    cs = sweepResult.cs
    vLo = minimum(c.ritzScoreVelAxis[1] for c in cs)
    vHi = maximum(c.ritzScoreVelAxis[end] for c in cs)
    velAxis = SawSpectrum._colon(vLo, RitzScore.VEL_STEP, vHi)
    absMap = Matrix{Float64}(undef, length(velAxis), length(cs))
    for ia in eachindex(cs)
        absMap[:, ia] = first(visRitzSpectrumData(cs[ia], velAxis))
    end
    plotData, cbLabel = scaleGreenData(absMap, scaleMode, "S")
    return sweepResult.angles, velAxis, plotData, cbLabel, absMap
end

"""
    ritzSpectrumCurveData(singleResult, scaleMode::Symbol = :log)
        -> (velAxis, plotData, yLabel, absS, cands)

Compute half of `@Vis/plotRitzSpectrumCurve.m`: one angle's Ritz spectrum `|S(v)|` on its own
`ritzScoreVelAxis`, scaled by [`scaleGreenData`](@ref) (label `"S"`), plus the stored Ritz
candidates (`ritzScoreCands`, for the optional overlay). Throws
`Vis:plotRitzSpectrumCurve:noRitzSpectrum` when no numerators are stored.
"""
function ritzSpectrumCurveData(singleResult::Union{CalcSingleHS,CalcSingleLayered},
                               scaleMode::Symbol = :log)
    visValidateChoice(scaleMode, (:log, :linear), "plotRitzSpectrumCurve", "scaleMode")
    visHasRitzSpectrum(singleResult) || throw(ArgumentError("Vis:plotRitzSpectrumCurve:" *
        "noRitzSpectrum -- This result stores no Ritz modal numerators (calcAll() stores them)."))
    absS, velAxis = visRitzSpectrumData(singleResult)
    plotData, yLabel = scaleGreenData(absS, scaleMode, "S")
    return velAxis, plotData, yLabel, absS, singleResult.ritzScoreCands
end

# -----------------------------------------------------------------------------------------
# dispersionCurvesData -- compute half of plotDispersionCurves.m
# -----------------------------------------------------------------------------------------

"""
    dispersionCurvesData(sweepResult::CalcWavevecAngSweep2Layered, angleIdx::Integer = 1)
        -> (branchModeNoMap, wavelengths, k, vel, velGroup, freq, highlightModeNos1)

Compute half of `guliwerr/visual/@Vis/plotDispersionCurves.m` (the Ritz-branch part; the
Green-background part is [`wavevecGreenBackgroundData`](@ref)): fixed-angle dispersion curves
for a `CalcWavevecAngSweep2Layered` wavelength sweep.

`sweepResult` must have `retentionLevel = :full` -- `DispersionBranches.track`'s own
requirement (this function calls it directly and lets its error propagate).
`branchModeNoMap` is `[nWav x nCompModes]`, the tracked-branch mode-number map;
`wavelengths` is `[nWav]`; `k = 2*pi ./ wavelengths` (`[nWav]`, standard angular wavenumber);
`vel`/`velGroup`/`freq` are `[nWav x nCompModes]` (`freq = vel / wavelength`, broadcast per
column); `highlightModeNos1` lists which of the FIRST wavelength step's branch mode numbers
were classified `:SAW`/`:limiting` there (the layered type list is `{SAW}` on every method
in MATLAB too) -- the branches a caller draws thicker/coloured rather than thin grey.

**Stops exactly where `plotDispersionCurves.m` calls `visDispersionCurveAxes`** -- the
"which of `wavelength`/`k`/`vel`/`velGroup`/`freq` goes on which axis for this `curveType`"
reshuffle is left to the plotting function (see this file's header).
"""
function dispersionCurvesData(sweepResult::CalcWavevecAngSweep2Layered, angleIdx::Integer = 1)
    branchModeNoMap = DispersionBranches.track(sweepResult, angleIdx)
    nWav, nCompModes = size(branchModeNoMap)
    wavelengths = sweepResult.wavelengths

    sm1 = sweepResult.sweeps[1].significantModesMap[angleIdx]
    isSig1Types = (:SAW, :limiting)
    highlightModeNos1 = [m.modeNo for m in sm1 if m.type in isSig1Types]

    vel = zeros(Float64, nWav, nCompModes)
    velGroup = zeros(Float64, nWav, nCompModes)
    for iw in 1:nWav
        ritz = sweepResult.sweeps[iw].cs[angleIdx].ritzSolver
        modeNoRow = @view branchModeNoMap[iw, :]
        vel[iw, :] = ritz.vels[modeNoRow]
        velGroup[iw, :] = ritz.velsGroup[modeNoRow]
    end
    k = 2 * pi ./ wavelengths
    freq = vel ./ wavelengths   # broadcast over columns

    return branchModeNoMap, wavelengths, k, vel, velGroup, freq, highlightModeNos1
end

# -----------------------------------------------------------------------------------------
# wavevecGreenBackgroundData -- compute half of plotDispersionCurves.m's local function
# plotWavevecGreenBackground
# -----------------------------------------------------------------------------------------

"""
    wavevecGreenBackgroundData(sweepResult::CalcWavevecAngSweep2Layered, angleIdx::Integer,
                                curveType::Symbol, scaleMode::Symbol = :log)
        -> (xAxis, velCommon, plotData, cbLabel)

Compute half of `guliwerr/visual/@Vis/plotDispersionCurves.m`'s local function
`plotWavevecGreenBackground`: Green's `|G13_eff|` heatmap (`greenAbsG13EffMap`) over `(wavelength-or-k, velocity)`,
fixed angle -- the wavevec-sweep counterpart of a `plotGreenMap` background.

`curveType` must be `:wavelength_vphase` or `:k_vphase` (velocity y-axis) -- the 4
frequency-based curve types are excluded because `freq = vel/wavelength` rescales the y-axis
per wavelength column, making the background a non-rectangular grid; any other `curveType`
throws `ArgumentError` (an explicit opt-in cannot be silently ignored, matching MATLAB's own
`plotDispersionCurves:backgroundNeedsVphase`). Throws if `sweepResult` was built without Green
data (`calcMode = :ritzOnly`).

Every wavelength step's own differently-ranged `greenVelAxis`/`greenAbsG13EffMap` column is
interpolated onto ONE shared velocity axis (`velCommon`) first -- linear, `NaN` outside a
wavelength's own local range (there is no Green data below that wavelength's cutoff) -- then
[`scaleGreenData`](@ref) applies the shared log/linear scaling, whose `NaN` contract
is exactly what keeps those out-of-range cells blank rather than painted
at the scale floor.
"""
function wavevecGreenBackgroundData(sweepResult::CalcWavevecAngSweep2Layered, angleIdx::Integer,
                                     curveType::Symbol, scaleMode::Symbol = :log)
    curveType in (:wavelength_vphase, :k_vphase) ||
        throw(ArgumentError("Vis:plotDispersionCurves:backgroundNeedsVphase -- " *
            "showGreenBackground needs curveType :wavelength_vphase or :k_vphase (velocity " *
            "y-axis) -- got :$curveType. The 4 frequency-based curve types would need a " *
            "non-rectangular background grid (freq = vel/wavelength rescales per column); " *
            "not supported."))

    nWav = length(sweepResult.sweeps)
    isempty(sweepResult.sweeps[1].greenAbsG13EffMap) &&
        throw(ArgumentError("Vis:plotDispersionCurves:noGreenData -- greenAbsG13EffMap is empty -- " *
            "sweepResult must be built with calcMode = :full to populate Green's function " *
            "data (showGreenBackground needs it)."))

    nVelGreen = length(sweepResult.sweeps[1].greenVelAxis)
    vMin = Inf
    vMax = -Inf
    for iw in 1:nWav
        axLocal = sweepResult.sweeps[iw].greenVelAxis
        vMin = min(vMin, minimum(axLocal))
        vMax = max(vMax, maximum(axLocal))
    end
    velCommon = collect(range(vMin, vMax, length = nVelGreen))

    C = fill(NaN, nVelGreen, nWav)
    for iw in 1:nWav
        axLocal = sweepResult.sweeps[iw].greenVelAxis
        gLocal = Float64.(sweepResult.sweeps[iw].greenAbsG13EffMap[:, angleIdx])
        itp = interpolate((collect(Float64, axLocal),), gLocal, Gridded(Linear()))
        etp = extrapolate(itp, NaN)
        C[:, iw] = etp.(velCommon)
    end

    plotData, cbLabel = scaleGreenData(C, scaleMode, "G_{13}^{eff}")

    xAxis = curveType === :wavelength_vphase ? sweepResult.wavelengths :
            2 * pi ./ sweepResult.wavelengths   # decreasing for increasing wavelengths;
                                                  # pcolorCenteredGrid supports either direction
    return xAxis, velCommon, plotData, cbLabel
end
