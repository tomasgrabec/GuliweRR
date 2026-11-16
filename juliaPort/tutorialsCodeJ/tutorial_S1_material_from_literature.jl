## Tutorial S1 (code) -- building a Material from literature data: Voigt input, conversions, orientation, averaging
# Narrated version: juliaPort/tutorialsJ/tutorial_S1_material_from_literature.jl.
# Run: at the prompt of `julia --project=. -t auto` started in your own environment with GuliwerrJ and GLMakie
#   (juliaPort/docsJ/gettingStartedJ.md, section 2): include("/path/to/GuliweRR/juliaPort/tutorialsCodeJ/tutorial_S1_material_from_literature.jl")

using GuliwerrJ
using Printf

## 1. Full 6x6 Voigt matrix input, symmetry given explicitly (SiC, cubic)
voigtSiC = [390.0 134.0 134.0   0.0   0.0   0.0
            134.0 390.0 134.0   0.0   0.0   0.0
            134.0 134.0 390.0   0.0   0.0   0.0
              0.0   0.0   0.0 253.0   0.0   0.0
              0.0   0.0   0.0   0.0 253.0   0.0
              0.0   0.0   0.0   0.0   0.0 253.0]   # GPa
symmCubic = 2                      # SymmClassLibrary() lists all classes
ecSiC = ElastConsts(voigtSiC, symmCubic)
@printf("[C11, C12, C44] = [%.1f, %.1f, %.1f] GPa\n", ecSiC.vector...)

## 2. Symmetry auto-detection -- exists but fragile, prefer an explicit symmNo
symmDetected = FindSymmClass_input(voigtSiC)
className, _ = SymmClassLibrary(symmDetected)
@printf("Auto-detected class %d (%s)\n", symmDetected, className)

## 3. Cubic (C11, C12, C44) <-> (B, c', c44) round trip
bulkShearC44, _   = CubicConversion(ecSiC.vector)     # [B, (C11-C12)/2, C44]
standardVector, _ = CubicDeconversion(bulkShearC44)
@printf("(B, c', c44) = [%.2f, %.2f, %.2f] GPa\n", bulkShearC44...)
@printf("back: [%.1f, %.1f, %.1f] GPa\n", standardVector...)

## 4. Orientation from direction rows: row 1 = surface normal, row 2 = propagation reference
dirMatSiC  = [0 0 1; 1 1 0]        # (001) surface, [110] reference
densitySiC = 3.21                  # g/cm^3
nameSiC    = "SiC (001)/<110>"
matSiC     = Material(; elConsts = ecSiC, rho = densitySiC, dirMat = dirMatSiC, name = nameSiC)

## 5. Anisotropy factor -- selects the CalcParams tier (cubic: exact Zener formula)
anisotropyFactor = ecSiC.shearAnisotropyRatio   # stored for isotropic/cubic; NaN otherwise -> calc_anisotropy_ratio(ec.tensor)
@printf("Anisotropy factor %.2f\n", anisotropyFactor)

## 6. Orientation from Bunge-ZXZ Euler angles (EBSD form)
phi1Deg = 30.0; PhiDeg = 45.0; phi2Deg = 10.0
dirMatEuler = Euler2DirMat(phi1Deg, PhiDeg, phi2Deg)
matSiCEuler = Material(; elConsts = ecSiC, rho = densitySiC, dirMat = dirMatEuler,
                         name = "SiC (Euler 30/45/10)")
phi1Back, PhiBack, phi2Back = DirMat2Euler(dirMatEuler)
@printf("Euler round trip (%.1f, %.1f, %.1f) deg\n", phi1Back, PhiBack, phi2Back)

## 7. Orientation from Miller indices -- normal = plane normal, reference = direction
latticeParams = NamedTuple()       # cubic needs none; Miller2DirMat(symmNo) lists fields
normalHkl     = [1, 1, 1]
propUvw       = [1, -1, 0]
normalDirCos  = Miller2DirMat(normalHkl, symmCubic, latticeParams, "planeNormal")
propDirCos    = Miller2DirMat(propUvw, symmCubic, latticeParams, "direction")
dirMatMiller  = vcat(permutedims(normalDirCos), permutedims(propDirCos))   # stack as 2 rows
matSiC111     = Material(; elConsts = ecSiC, rho = densitySiC, dirMat = dirMatMiller,
                           name = "SiC (111)/[1-10]")

## 8. Rotate the in-plane reference about the surface normal
inPlaneRotationDeg = 30.0
dirMatFull    = CompleteDirMat(dirMatSiC)   # RotDirMatFromAngle needs the 3x3 frame
dirMatRotated = RotDirMatFromAngle(dirMatFull, inPlaneRotationDeg)
matSiCRotated = Material(; elConsts = ecSiC, rho = densitySiC, dirMat = dirMatRotated,
                           name = "SiC (001), 30 deg")

## 9. Polycrystalline film -- Voigt / Reuss / Hill isotropic averages (C11, C44)
matAlCrystal, _ = MaterialFromLibrary("Al")   # Al single crystal from the materials library
ecAlCubic = matAlCrystal.elConsts              # [C11, C12, C44] GPa
alC11Voigt, alC44Voigt = VoigtAverages(ecAlCubic.voigt)   # upper bound
alC11Reuss, alC44Reuss = ReussAverages(ecAlCubic.voigt)   # lower bound
alC11Hill,  alC44Hill  = HillAverages(ecAlCubic.voigt)    # mean, usual choice
@printf("Al [C11, C44]: Voigt [%.2f %.2f], Reuss [%.2f %.2f], Hill [%.2f %.2f] GPa\n",
        alC11Voigt, alC44Voigt, alC11Reuss, alC44Reuss, alC11Hill, alC44Hill)
symmIsotropic = 1
ecAlHill  = ElastConsts([alC11Hill, alC44Hill], symmIsotropic)
densityAl = matAlCrystal.rho         # g/cm^3
matAlFilm = Material(; elConsts = ecAlHill, rho = densityAl, name = "Al film (Hill)")   # isotropic: no dirMat

## 10. Two-layer Structure from the Materials above
filmThicknessNm     = 200.0
relDepthSub         = 20.0         # "easy" tier ((001) principal plane)
legOrder            = 50
relDepthPlaceholder = 1.0          # film relDepth, set per wavelength
nMaxModes           = 500
cpSubstrate = CalcParams(; relDepth = relDepthSub, legOrder = legOrder)
cpFilm      = CalcParams(; relDepth = relDepthPlaceholder, legOrder = legOrder,
                           realThickness = filmThicknessNm)
structure   = Structure(; materIn = matSiC, calcParamsIn = cpSubstrate, name = "SiC_sub_Al_film",
                          nMaxModes = nMaxModes)
structure   = AddLayer(structure, matAlFilm, cpFilm)
# GUI (needs GLMakie; opens empty in this port, unlike MATLAB's pre-filled one):
#   using GLMakie; gui = structureControl()

## 11. Standard crystals from the materials library (default row, or by T / variant)
ListMaterialLibrary(; print = true)
matCu001, _ = MaterialFromLibrary("Cu"; dirMat = [0.0 0 1; 1 1 0])   # Cu(001), reference [110]
matMgO, _   = MaterialFromLibrary("MgO"; dirMat = [0.0 0 1; 1 0 0])  # MgO(001), reference [100]
matTi, entryTi = MaterialFromLibrary("titanium")                      # name or alias, any case
println("$(matTi.name): rho = $(matTi.rho) g/cm^3, C11 = $(matTi.elConsts.voigt[1, 1]) GPa ",
        "(source: $(entryTi.source))")
