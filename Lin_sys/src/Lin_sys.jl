module Lin_sys

using LinearAlgebra

include("subs.jl")
include("matrix_checking.jl")
include("tridiag.jl")
include("sys_solving.jl")

greet() = print("Hello world!")

end # module Lin_sys
