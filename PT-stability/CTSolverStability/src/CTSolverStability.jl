module CTSolverStability

using ForwardDiff
using LinearAlgebra
using StructArrays

include("nlsolve.jl")
include("types.jl")
include("solvecubic.jl")
include("constants.jl")
include("initials.jl")
include("stability.jl")

include("vanderwaals.jl")



end # module CTSolverStability
