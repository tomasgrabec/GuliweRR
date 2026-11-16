## Tutorial: building a Material from real literature/lab data (Julia port)
#
# Same physics as examples/tutorials/tutorial_S1_material_from_literature.m: every other
# tutorial in this series starts from elastic constants already in this codebase's own compact
# vector shorthand ([C11,C12,C44] for cubic, etc.), at the default crystal orientation. Data
# from a paper rarely look like that: a full Voigt matrix, or another parameterisation (e.g.
# bulk modulus and shear constants), and a sample cut that is not the default. This script
# builds a Material from such data: full-Voigt input, symmetry detection, converting
# parameterisations, and the crystal orientation from direction rows, Euler angles, Miller
# indices or an in-plane rotation; then isotropic averaging for a polycrystalline film and the
# structure GUI. The SiC constants below are example input values; they make the round trips checkable.
#
# Copy the parts you need as a starting point for building a new Material from your own data.
#
# See tutorial_01_isotropic_halfspace.jl's header for how to run a script in this series.
#
# Julia concepts introduced in this tutorial:
#   - none new -- this script is pure numerics (no plotting), reusing only `using GuliwerrJ`,
#     positional `ElastConsts(...)`, keyword `Material(...)`, and `@printf` from earlier
#     tutorials. One genuine API difference from MATLAB is called out inline where it occurs
#     (the shear anisotropy ratio is a stored field here, `ec.shearAnisotropyRatio`).

using GuliwerrJ
using Printf

## 1. Full Voigt matrix input, symmetry given EXPLICITLY (the safe,
# recommended path). A paper very often gives you the whole 6x6, not this
# codebase's own compact per-symmetry vector -- ElastConsts's constructor
# accepts a full Voigt matrix directly, no conversion needed on your part.
# SiC (3C, cubic): C11=390, C12=134, C44=253 GPa.
voigtSiC = [390.0  134.0  134.0    0.0    0.0    0.0
            134.0  390.0  134.0    0.0    0.0    0.0
            134.0  134.0  390.0    0.0    0.0    0.0
              0.0    0.0    0.0  253.0    0.0    0.0
              0.0    0.0    0.0    0.0  253.0    0.0
              0.0    0.0    0.0    0.0    0.0  253.0]   # GPa, full symmetric 6x6
# A multi-row matrix literal spanning several SOURCE LINES needs no line-continuation marker at
# all in Julia -- unlike MATLAB's `...`, a newline INSIDE an open `[...]` is simply another row
# separator (exactly like a `;` would be), so this reads almost identically to the MATLAB
# original with the `...`s just dropped. As in tutorial 9: spaces (not commas) between elements
# on the same row.

symmCubic = 2   # see SymmClassLibrary() for every class's number/constant order
ecSiC = ElastConsts(voigtSiC, symmCubic)

@printf("Parsed from full Voigt matrix: [C11,C12,C44] = [%.1f, %.1f, %.1f] GPa\n", ecSiC.vector...)
# `ecSiC.vector...` -- the trailing `...` SPLATS the 3-element vector into three separate
# arguments, one per remaining `%.1f` in the format string. Julia's `@printf` (unlike MATLAB's
# `fprintf`) does not itself loop over a vector argument to fill repeated `%`-specifiers, so the
# vector must be splatted into individual scalar arguments explicitly -- `...` is Julia's own
# "unpack this container as separate function arguments" syntax, used here for exactly that.

## 2. Symmetry auto-detection EXISTS (FindSymmClass_input, used internally
# when symmNo is omitted) but is fragile: it uses `unique` on the Voigt
# matrix, which silently misidentifies the class if two independent
# constants happen to be numerically equal by coincidence.
# ALWAYS supply symmNo explicitly, as in section 1 above, once you know it.
symmDetected = FindSymmClass_input(voigtSiC)
name, _ = SymmClassLibrary(symmDetected)
@printf("Auto-detected symmetry class: %d (%s)\n", symmDetected, name)

## 3. Converting between literature parameterisations. Some papers report
# cubic elasticity as (B, c', c44) -- bulk modulus, tetragonal shear
# constant c'=(C11-C12)/2, c44 -- instead of (C11,C12,C44) directly.
# CubicConversion goes standard -> (B,c',c44); CubicDeconversion goes back.
# Demonstrated as a round trip so the numbers are checkable against
# section 1's known-correct SiC values.
BcPrimeC44, _ = CubicConversion(ecSiC.vector)
@printf("SiC as (B, c', c44)         = [%.2f, %.2f, %.2f] GPa\n", BcPrimeC44...)

## Suppose a DIFFERENT paper reports SiC this way instead -- recover the
# standard (C11,C12,C44) form ElastConsts needs from it:
recovered, _ = CubicDeconversion(BcPrimeC44)
@printf("Recovered [C11,C12,C44]     = [%.1f, %.1f, %.1f] GPa (matches section 1)\n", recovered...)
# Tetragonal and hexagonal have the same pair: TetragonalConversion /
# TetragonalDeconversion and HexagonalConversion / HexagonalDeconversion.

## 4. Real crystallographic orientation -- dirMat, NOT RotateTensor.
# RotateTensor exists (used internally by e.g. plotVelSurface/
# plotGroupVelSurface, which sweep EVERY orientation), but for a single fixed
# crystal cut -- the normal case -- Material's dirMat is the right tool: give
# it the cut's Miller-index-style directions directly IN THE CRYSTAL'S OWN
# AXES, no manual tensor rotation needed.
#
# dirMat row convention:
#   row 1 = surface normal
#   row 2 = in-plane propagation reference (lab angle = 0)
#   (row 3, in-plane, completes the right-handed frame -- can be omitted,
#    a 2x3 dirMat is accepted and completed automatically)
#
# SiC (001) surface, <110>-type in-plane reference direction:
dirMatSiC = [0 0 1    # surface normal  = [001]
             1 1 0]    # propagation ref = [110]

rhoSiC     = 3.21   # g/cm^3
matNameSiC = "SiC (001)/<110>"
matSiC = Material(; elConsts = ecSiC, rho = rhoSiC, dirMat = dirMatSiC, name = matNameSiC)

## 5. Zener ratio -- selects the CalcParams (relDepth/legOrder) tier
# (docs/pitfalls.md MTH1). Cubic: the exact Zener ratio; other symmetries:
# (v_shear_max/v_shear_min)^2 from an angular scan.
A = matSiC.elConsts.shearAnisotropyRatio
@printf("SiC Zener ratio 2*C44/(C11-C12) = %.2f\n", A)
println("-> ratio < 5 on (001): start from the \"easy\" tier")
println("   (relDepth=20, legOrder=50); check convergence (docs/pitfalls.md MTH1).")
# The factor is a stored FIELD of ElastConsts (`ec.shearAnisotropyRatio`), where MATLAB has the
# method `ec.GetAnisotropyRatio()`. It is filled only where a closed form exists (isotropic,
# cubic -- as here); for every other symmetry class the field is NaN and
# `calc_anisotropy_ratio(ec.tensor)` runs the angular scan instead.

## 6. Orientation from EBSD-style Bunge-ZXZ Euler angles (phi1, Phi, phi2) -- the form
# texture/EBSD data arrives in. Euler2DirMat returns a ready-to-use dirMat (row 1 = surface
# normal) -- no manual row reorder. (Convention matches Dryburgh et al. 2020/2022; not yet
# cross-checked against a real EBSD worked example -- docs/reference/functionReference.md.)
phi1Deg = 30.0   # deg
PhiDeg  = 45.0   # deg
phi2Deg = 10.0   # deg
dirMatEuler = Euler2DirMat(phi1Deg, PhiDeg, phi2Deg)
matSiCEuler = Material(; elConsts = ecSiC, rho = rhoSiC, dirMat = dirMatEuler,
                         name = "SiC (Euler 30/45/10)")
phi1Back, PhiBack, phi2Back = DirMat2Euler(dirMatEuler)   # round trip; returns a 3-tuple
@printf("Euler round trip: (%.1f, %.1f, %.1f) deg\n", phi1Back, PhiBack, phi2Back)

## 7. Orientation from Miller indices -- Miller2DirMat turns an index into a direction-cosine
# row. A surface normal is a PLANE normal (hkl), a propagation reference is a real-space
# DIRECTION [uvw] -- the two differ for every non-cubic lattice. Call Miller2DirMat(symmNo)
# alone to list which latticeParams fields that symmetry needs (cubic needs none).
latticeParamsCubic = NamedTuple()   # cubic: no lattice-shape parameters (MATLAB: struct())
normalDirCos = Miller2DirMat([1, 1, 1], symmCubic, latticeParamsCubic, "planeNormal")
propDirCos   = Miller2DirMat([1, -1, 0], symmCubic, latticeParamsCubic, "direction")
dirMatMiller = vcat(permutedims(normalDirCos), permutedims(propDirCos))   # stack as 2 rows
# `permutedims` turns each returned Vector into a 1-row matrix so `vcat` stacks them as ROWS
# (MATLAB's `[a; b]` of two row vectors); a bare `vcat(a, b)` would give one long Vector.
matSiC111 = Material(; elConsts = ecSiC, rho = rhoSiC, dirMat = dirMatMiller,
                       name = "SiC (111)/[1-10]")

## 8. Same cut, in-plane reference rotated -- RotDirMatFromAngle rotates rows 2-3 about the
# surface normal. It needs the COMPLETED 3x3 right-handed frame, not the 2-row input --
# CompleteDirMat builds it.
inPlaneRotationDeg = 30.0   # deg
dirMatFull    = CompleteDirMat(dirMatSiC)
dirMatRotated = RotDirMatFromAngle(dirMatFull, inPlaneRotationDeg)
matSiCRotated = Material(; elConsts = ecSiC, rho = rhoSiC, dirMat = dirMatRotated,
                           name = "SiC (001), 30 deg")
# (A sweep's own `angles` argument already does this rotation per angle -- only rotate dirMat
# by hand to move where "angle = 0" points.)

## 9. Polycrystalline (e.g. sputtered) film: isotropic constants AVERAGED from single-crystal
# ones. Three bounds, each returning the isotropic pair (M, G) = (C11, C44): Voigt (upper,
# uniform strain), Reuss (lower, uniform stress), Hill (their arithmetic mean -- the usual
# practical choice). In this port they are plain functions of the 6x6 Voigt matrix
# (`ec.voigt`), not methods of the ElastConsts object as in MATLAB (`ec.GetHillAverages()`).
matAlCrystal, _ = MaterialFromLibrary("Al")   # Al single crystal from the library
ecAlCubic       = matAlCrystal.elConsts          # [C11, C12, C44] = [107.3, 60.8, 28.3] GPa
alC11Voigt, alC44Voigt = VoigtAverages(ecAlCubic.voigt)
alC11Reuss, alC44Reuss = ReussAverages(ecAlCubic.voigt)
alC11Hill,  alC44Hill  = HillAverages(ecAlCubic.voigt)
@printf("Al isotropic [C11, C44]: Voigt [%.2f, %.2f], Reuss [%.2f, %.2f], Hill [%.2f, %.2f] GPa\n",
        alC11Voigt, alC44Voigt, alC11Reuss, alC44Reuss, alC11Hill, alC44Hill)
symmIsotropic = 1
ecAlHill  = ElastConsts([alC11Hill, alC44Hill], symmIsotropic)
rhoAl     = matAlCrystal.rho   # g/cm^3, the library's single-crystal density
matAlFilm = Material(; elConsts = ecAlHill, rho = rhoAl, name = "Al film (Hill-averaged)")
# dirMat simply omitted: an isotropic material has no orientation to set

## 10. Two-layer Structure from the Materials above -- and the structure-building GUI.
filmThicknessNm = 200.0   # real film thickness [nm]
cpSubstrate = CalcParams(; relDepth = 20.0, legOrder = 50)   # "easy" tier, see section 5
cpFilm      = CalcParams(; relDepth = 1.0, legOrder = 50,    # relDepth placeholder, set per wavelength
                       realThickness = filmThicknessNm)
strDemo = Structure(; materIn = matSiC, calcParamsIn = cpSubstrate, name = "SiC_sub_Al_film",
                      nMaxModes = 500)
strDemo = AddLayer(strDemo, matAlFilm, cpFilm)
# PORT GAP: MATLAB's Gui.structureControl(initialStructure) opens pre-filled from a Structure;
# this port's structureControl() always opens EMPTY (layers are entered in the window). It
# needs GLMakie, which this otherwise plotting-free script does not load -- uncomment to try:
#   using GLMakie
#   gui = structureControl()

## 11. Standard crystals from the materials library -- no typing of constants.
# ListMaterialLibrary(print = true) prints every stored data set (constants, density,
# temperature, source, status); MaterialFromLibrary builds the Material, with the orientation
# given exactly as in sections 4-8, and also returns the row (source, reference, notes). Every
# row names its primary measurement; a material with several rows (another temperature or
# source) is chosen by `T` or `variant`, otherwise its default row. Same data file as MATLAB's
# Material.FromLibrary.
ListMaterialLibrary(; print = true)
matCu001, _ = MaterialFromLibrary("Cu"; dirMat = [0.0 0 1; 1 1 0])   # Cu(001), reference [110]
matMgO, _   = MaterialFromLibrary("MgO"; dirMat = [0.0 0 1; 1 0 0])  # MgO(001) substrate, reference [100]
matTi, entryTi = MaterialFromLibrary("titanium")                      # aliases work too
println("$(matTi.name): rho = $(matTi.rho) g/cm^3, C11 = $(matTi.elConsts.voigt[1, 1]) GPa ",
        "(source: $(entryTi.source))")

## matSiC is now exactly the kind of Material every other tutorial's own
# "Material" section builds directly -- continue from here into
# tutorial_01_isotropic_halfspace.jl/tutorial_02_angular_map.jl's own next step
# (CalcParams -> Structure -> the Calc-family sweep functions).
