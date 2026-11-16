# =========================================================================================
# GreenSolverTypes.jl
#
# Port of the TYPE SHAPE of `guliwerr/solvers/green/@GreenSolver/GreenSolver.m` (177 lines):
# the value class' 19 `SetAccess = protected` properties, its `GreenSolver(forceDir,
# displDir)` constructor, and a handful of predicates that replace MATLAB's ad-hoc
# `isempty(obj.dispField)` guards.
#
# `calc` itself lives in GreenSolver.jl.
#
# Plain top-level definitions, not `module GreenSolverTypes` -- one module per top-level
# folder (`green/`), not one per file. See `src/core/CoreTypes.jl`'s header comment: a
# per-file `module X` containing `struct X` breaks bare `X(...)` construction one scope up,
# and `struct GreenSolver` is exactly that case.
#
# MUTABLE vs IMMUTABLE -- checked, not assumed. `classdef GreenSolver` (no `< handle`) is a
# MATLAB *value* class, and every property is `SetAccess = protected`, so only `calc.m`
# itself may write them. Nothing in the MATLAB code writes any of the 19 properties from
# outside `@GreenSolver/` -- every external site (`CalcIsofreqSingleLayered.m`,
# `CalcCommon.m`, `Vis.calcGreenDispField2D`/`...Layered`, `ModalShapeViewerGUI`,
# `IsofreqModalShapeViewerGUI`) only ever READS. So this
# is an ordinary **immutable `struct` with mutable array fields**, the natural translation
# of a built-once-then-read result class: cheap to pass and return (a struct holding a
# `Matrix` stores a pointer, never a copy), value semantics preserved, and data races on
# the scalar fields impossible by construction.
# =========================================================================================

"""
    GreenSolver

Surface dynamic Green's function `G_IJ(k_par, omega)` for anisotropic solids -- half-space
and layered geometry (Every et al. 1997, Every 2002, Maznev et al. 1999). Result object of
`GreenSolver.calc`.

Port of `guliwerr/solvers/green/@GreenSolver/GreenSolver.m`'s property block. Field names
are kept MATLAB-verbatim (`camelCase`, `G_eff`, `GFull_minus`, ...) so a port can be checked
against the original line by line; only the *functions* are renamed to Julia style.

**"Empty" convention.** MATLAB writes `[]` into the properties a given call does not
populate (`G_minus`/`G_eff` when `computeMinus = false`; the half-space `dispField`/
`depthAxis` for a layered structure; the `*Layered` group for a half-space). The half-space
`k3Out`/`polarizationU`/`partialAmplitude` fields of earlier MATLAB versions are not ported.
Those that remain stay *typed but zero-sized* containers rather than becoming
`Union{Nothing,...}` -- `isempty(obj.G_eff)` is then the exact analogue of MATLAB's
`isempty(obj.G_eff)`, and the struct stays fully concretely typed (no `Union` field to
de-optimise every access). The predicates [`has_half_space_full_output`](@ref) and
[`has_layered_full_output`](@ref) below name the two composite checks `calc.m`'s callers
actually make.

# Fields
- `forceDir`, `displDir` -- I (force) and J (displacement) direction indices; `(1, 3)` is
  G_13, the TGS component. Standard GULIWERR convention: 1 = propagation, 3 = surface normal.
- `velRange` -- phase-velocity axis used in the calculation [km/s].
- `G` -- `G_IJ(+k_par)`, one entry per velocity.
- `G_minus` -- `G_IJ(-k_par) = G_JI(+k_par)` (elastodynamic reciprocity; only `+k_par` is
  ever solved); empty when `computeMinus = false`.
- `G_eff` -- parity-signed mean `(G_+ + (-1)^n G_-)/2`, `n` = how many of `(forceDir,
  displDir)` equal 1 (Maznev 1999): `G13_eff = (G13_+ - G31_+)/2`, `G33_eff = G33_+`. Empty
  when `computeMinus = false`.
- `computeMinus` -- whether the minus/eff fields were filled (free either way).
- `GFull`, `GFull_minus`, `GFull_eff` -- the same three quantities for all nine `(I,J)`
  pairs at once, `3 x 3 x nVel`.
- `dispField` -- displacement-field profile `u(depth, component, forceDir, velocity)`;
  half-space + `nargout > 1` only.
- `depthAxis` -- depth axis for `dispField`, normalised units (`2*pi` = one wavelength).
- `k3OutLayered`, `polarizationULayered`, `partialAmplitudeLayered` -- their layered
  counterparts, one ragged entry per layer in **Structure order** (index 1 = substrate,
  index `nLay` = outermost film): `6 x nVel` / `3 x 6 x nVel` / `6 x 3 x nVel` for a film,
  `3`-wave shapes for the substrate.
- `hScaleOut` -- the `hScale` actually used to build the `*Layered` group (1 for
  `isoWavevec`, `omega0/v` for `isoFreq`); `NaN` when unset.
"""
Base.@kwdef struct GreenSolver
    forceDir::Int = 1
    displDir::Int = 3

    velRange::Vector{Float64} = Float64[]
    G::Vector{ComplexF64} = ComplexF64[]
    G_minus::Vector{ComplexF64} = ComplexF64[]
    G_eff::Vector{ComplexF64} = ComplexF64[]

    computeMinus::Bool = true

    GFull::Array{ComplexF64,3} = Array{ComplexF64,3}(undef, 0, 0, 0)
    GFull_minus::Array{ComplexF64,3} = Array{ComplexF64,3}(undef, 0, 0, 0)
    GFull_eff::Array{ComplexF64,3} = Array{ComplexF64,3}(undef, 0, 0, 0)

    # --- optional full output, half-space (nLay == 1) only ---
    dispField::Array{ComplexF64,4} = Array{ComplexF64,4}(undef, 0, 0, 0, 0)
    depthAxis::Vector{Float64} = Float64[]

    # --- optional full output, layered (nLay > 1) only; Structure order, ragged ---
    k3OutLayered::Vector{Matrix{ComplexF64}} = Matrix{ComplexF64}[]
    polarizationULayered::Vector{Array{ComplexF64,3}} = Array{ComplexF64,3}[]
    partialAmplitudeLayered::Vector{Array{ComplexF64,3}} = Array{ComplexF64,3}[]

    hScaleOut::Float64 = NaN
end

"""
    GreenSolver(forceDir, displDir = 3) -> GreenSolver

Positional convenience constructor mirroring `GreenSolver.m`'s own
`GreenSolver(forceDir, displDir)` -- everything else takes the class-default "not computed"
empty value. MATLAB's `nargin >= 1 && ~isempty(forceDir)` guard is expressed by the
argument simply being optional here, and the keyword form `GreenSolver(; forceDir = ...,
displDir = ...)` (generated by `Base.@kwdef`) covers the `[]`-means-default call shape.

Two positional arguments, not 19 -- deliberately well short of the struct's own field count,
so this cannot collide with `@kwdef`'s all-fields positional constructor during package
precompilation (the same trap as for `Material`).
"""
GreenSolver(forceDir::Integer, displDir::Integer = 3) =
    GreenSolver(; forceDir = Int(forceDir), displDir = Int(displDir))

# -----------------------------------------------------------------------------------------
# Predicates -- the composite `isempty(...)` guards `calc.m`'s own callers make
# -----------------------------------------------------------------------------------------

"""
    nvel(obj::GreenSolver) -> Int

Number of velocity points the object was solved on (`length(obj.velRange)`). MATLAB writes
`numel(obj.velRange)` inline at every site; named here so downstream code reads the same way
whichever array it is about to index.
"""
nvel(obj::GreenSolver) = length(obj.velRange)

"""
    has_minus(obj::GreenSolver) -> Bool

Whether the `k1 = -k_par` branch was solved *and* its results are present. Equivalent to
MATLAB's `obj.computeMinus && ~isempty(obj.G_minus)`; the two can only disagree on a
half-built object.
"""
has_minus(obj::GreenSolver) = obj.computeMinus && !isempty(obj.G_minus)

"""
    has_half_space_full_output(obj::GreenSolver) -> Bool

Whether the half-space `fullOutput` profile is populated -- the guard
`Vis.calcGreenDispField2D` makes before using `dispField`/`depthAxis` (both stay empty for a
layered structure, which populates the `*Layered` group instead).
"""
has_half_space_full_output(obj::GreenSolver) =
    !isempty(obj.dispField) && !isempty(obj.depthAxis)

"""
    has_layered_full_output(obj::GreenSolver) -> Bool

Whether the layered `nargout > 1` per-layer diagnostics are populated -- the guard
`Vis.calcGreenDispField2DLayered:75` makes before reading `k3OutLayered` and friends.
"""
has_layered_full_output(obj::GreenSolver) = !isempty(obj.k3OutLayered)

export GreenSolver, nvel, has_minus, has_half_space_full_output, has_layered_full_output
