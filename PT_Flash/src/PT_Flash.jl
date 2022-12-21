module PT_Flash

using StructArrays
using LinearAlgebra
using ForwardDiff

include("types.jl")
include("constants.jl")
include("nlsolve.jl")
include("solvecubic.jl")
include("vanderwaals.jl")
include("flashsolver.jl")
include("initials.jl")

end # module PT_Flash
