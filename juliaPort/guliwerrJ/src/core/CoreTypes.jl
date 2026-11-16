"""
    CoreTypes

Namespace for GULIWERR's `core/` layer (Legendre basis, `CalcParams`, `ElastConsts`,
`Material`, `Structure`, `ParallelPolicy`) plus the MATLAB-builtin compatibility shims.

One module per top-level folder, not per file: each file below is `include`d as a set of
plain top-level definitions, not wrapped in its own `module`. The reason: a per-file
`module X` containing `struct X` (the common case here, since every ported class is named
after its file) breaks bare construction one scope up. The convention applies to every
module in the tree (`Ritz`, `Green`, `Christoffel`, `ForwardCalcs`, `InverseCalcs`,
`Visual`), not just `CoreTypes`.

Include order is a real dependency constraint here (unlike MATLAB's path-based lookup):
`Legendre.jl` before `CalcParams.jl` (CalcParams's constructor calls Legendre's functions).
"""
module CoreTypes

include("MatlabCompat.jl")
include("Legendre.jl")
include("CalcParams.jl")
include("ElastConsts.jl")
include("ElastConstsSymmetry.jl")
include("ElastConstsDerived.jl")
include("Material.jl")
include("Orientation.jl")
include("MaterialLibrary.jl")
include("Structure.jl")
include("BlasBackend.jl")
include("ParallelPolicy.jl")
include("Version.jl")
include("ProbeKernel.jl")

# BLAS/LAPACK backend (BlasBackend.jl): runs after AppleAccelerate's own __init__ (a
# dependency loads first), so it sees the forwarding Accelerate set up, or undoes it on opt-out.
__init__() = (_init_blas_backend!(); nothing)

end # module CoreTypes
