## Tutorial 01 (code) -- isotropic half-space: one SAW velocity, checked against v_SAW ~ 0.92*v_T
# Narrated version: juliaPort/tutorialsJ/tutorial_01_isotropic_halfspace.jl.
# Run: julia --project=juliaPort/guliwerrJ -t auto juliaPort/tutorialsCodeJ/tutorial_01_isotropic_halfspace.jl (after the one-off Pkg.instantiate(), see juliaPort/docsJ/gettingStartedJ.md)

using GuliwerrJ
using Printf

## Material -- isotropic Al (isotropic class takes [C11, C44] only), Hill average of the library's single crystal
matAlCrystal, _ = MaterialFromLibrary("Al")
alC11, alC44 = HillAverages(matAlCrystal.elConsts.voigt)   # GPa
alConsts    = [alC11, alC44]
symmNo      = 1                    # 1 = isotropic (SymmClassLibrary() lists all)
density     = matAlCrystal.rho     # g/cm^3
elastConsts = ElastConsts(alConsts, symmNo)
material    = Material(; elConsts = elastConsts, rho = density)

## Discretisation -- "easy" tier, legOrder >= 2*relDepth (docs/pitfalls.md MTH1)
relDepth   = 20.0                  # domain depth, in wavelengths
legOrder   = 50                    # Legendre basis order across that depth
calcParams = CalcParams(; relDepth = relDepth, legOrder = legOrder)

## Structure -- one half-space layer
structureName = "al_isotropic_halfspace"
nMaxModes     = 500                # max Ritz modes computed/retained
structure     = Structure(; materIn = material, calcParamsIn = calcParams,
                            name = structureName, nMaxModes = nMaxModes)

## Sweep -- one angle, Ritz only
angleDeg       = 0.0               # isotropic: any angle gives the same result
calcMode       = :ritzOnly         # :ritzOnly | :full (adds Christoffel + Green)
retentionLevel = :saw_only         # :full | :saw_only | :none
sweep = angSweepHS(; structure = structure, angles = [angleDeg], calcMode = calcMode,
                     retentionLevel = retentionLevel)
sweep = calc(sweep)

## Result -- SAW velocity vs. transverse bulk velocity
vSAW = sweep.sawVelMap[1][1]       # km/s, first angle, first SAW
vT   = sqrt(alC44 / density)       # km/s
@printf("v_SAW = %.4f km/s, v_T = %.4f km/s, ratio = %.3f\n", vSAW, vT, vSAW / vT)
