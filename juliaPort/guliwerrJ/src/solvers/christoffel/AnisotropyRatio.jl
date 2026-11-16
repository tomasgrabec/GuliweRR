# =========================================================================================
# AnisotropyRatio.jl -- `include`d directly into `module Christoffel` (see Christoffel.jl's
# header for the "one module per folder" rationale).
#
# Port of guliwerr/core/@ElastConsts/CalcAnisotropyRatio.m, relocated here.
# `ElastConsts.CalcAnisotropyRatio` calls `Christoffel.phaseVels`, while `Christoffel`
# itself depends on `ElastConsts` (`SymmClassLibrary`, `Any2Tensor`, `Tensor2Vector`, the
# constructor). MATLAB does not care about this mutual reference (both are static-method
# classdefs resolved at call time); Julia's `include` order does -- keeping this function in
# `ElastConsts.jl` (core/, loaded BEFORE Christoffel) would create a forward reference to a
# not-yet-defined function. Defining the function here, in the module that already depends on
# `ElastConsts` in the other direction, breaks the cycle cleanly: `ElastConsts.jl` has ZERO
# solver dependencies, and `Christoffel`'s dependency on `CoreTypes` stays one-directional.
#
# The function's only other MATLAB caller, `ElastConsts.GetAnisotropyRatio`'s
# cubic/non-cubic dispatch (`ElastConsts.m`), is likewise NOT ported into `ElastConsts.jl`
# for the same reason -- `GreenSolver.calc.m`'s `NEAR_ISO_ANGLE_STEP` near-isotropy test
# is this function's only other real call site, and it already sits above `Christoffel` in
# the dependency graph, so it calls `calc_anisotropy_ratio` directly from here rather than
# through an `ElastConsts` wrapper.
#
# NAMING: `calc_anisotropy_ratio` is the one function in `solvers/christoffel/` that is NOT
# kept MATLAB-verbatim (`CalcAnisotropyRatio`); every other function here keeps its MATLAB
# casing (`matrix`, `phaseVels`, ...). The exact Zener formula for cubic symmetry
# (`ElastConsts.GetAnisotropyRatio`'s cubic branch, which needs no Christoffel scan at all)
# lives with `ElastConsts.jl`; this file only ever holds the general angular-scan method.
# =========================================================================================

"""
    calc_anisotropy_ratio(normElTen, angleStep=4.0) -> Float64

Generalised anisotropy ratio `A`, as `(fastest shear / slowest shear)^2` over all
propagation directions.

INPUT: `normElTen` -- elastic tensor (3x3x3x3), density-normalised; `angleStep` -- angular
step in degrees (default `4`).
OUTPUT: generalised anisotropy ratio (`>= 1`; `= 1` for an isotropic material).

**Prefer the exact Zener formula (`GetAnisotropyRatio`'s cubic branch, ported with
`ElastConsts.jl`) for cubic symmetry** -- this angular scan is the general-symmetry
fallback and does not know about symmetry-specific closed forms.

Port of `guliwerr/core/@ElastConsts/CalcAnisotropyRatio.m`, relocated here -- see this
file's header for why, and for the deliberate `calc_anisotropy_ratio` renaming (the sole
exception to the otherwise MATLAB-verbatim naming in `solvers/christoffel/`).
"""
function calc_anisotropy_ratio(normElTen::AbstractArray{<:Real,4}, angleStep::Real = 4.0)
    size(normElTen) == (3, 3, 3, 3) ||
        throw(ArgumentError("AnisotropyRatio:calc_anisotropy_ratio:badInput -- " *
                             "normElTen must be a 3x3x3x3 tensor"))
    angleStep > 0 ||
        throw(ArgumentError("AnisotropyRatio:calc_anisotropy_ratio:badInput -- " *
                             "angleStep must be a positive scalar (degrees)"))

    alpha = collect(-45.0:angleStep:45.0)
    tanAlpha = tand.(alpha)
    nAlpha = length(alpha)

    # scan all directions -- three orthogonal projections of the +/-45 deg tangent grid
    phaseVelsArr = zeros(Float64, 3, nAlpha, nAlpha, 3)
    dirV = zeros(Float64, 3, 3)

    for ii in 1:nAlpha, jj in 1:nAlpha
        dirV[1, :] = [tanAlpha[ii], tanAlpha[jj], 1.0]
        dirV[2, :] = [tanAlpha[ii], 1.0, tanAlpha[jj]]
        dirV[3, :] = [1.0, tanAlpha[ii], tanAlpha[jj]]
        for idir in 1:3
            d = @view dirV[idir, :]
            d ./= norm(d)
            phaseVelsArr[idir, ii, jj, :] = phaseVels(normElTen, d)
        end
    end

    # shear modes only (mode 1 is longitudinal)
    phaseVels2 = zeros(Float64, 2, nAlpha, 3 * nAlpha)
    for iwt in 2:3
        phaseVels2[iwt - 1, :, :] = hcat(phaseVelsArr[1, :, :, iwt],
                                          phaseVelsArr[2, :, :, iwt],
                                          phaseVelsArr[3, :, :, iwt])
    end

    # velocity^2 ratio -> elastic-stiffness units
    return (maximum(phaseVels2) / minimum(phaseVels2))^2
end

export calc_anisotropy_ratio
