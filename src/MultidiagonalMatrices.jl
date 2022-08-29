module MultidiagonalMatrices
using LinearAlgebra
using SparseArrays
using StaticArrays


include("multidiagonalmatrix.jl")
export MultidiagonalMatrix
export istridiagonal,fdsizes,mdrand,mdzeros

include("base.jl")
include("linearalgebra.jl")
include("sparsearrays.jl")

include("solvers.jl")
export tdma,tdma!


end # module MultidiagonalMatrices
