module MultiDiagonalMatrices
using LinearAlgebra
using SparseArrays
using StaticArrays


include("multidiagonalmatrix.jl")
export MultiDiagonalMatrix
export istridiagonal,fdsizes,mdrand,mdzeros

include("base.jl")
include("linearalgebra.jl")
include("sparsearrays.jl")

include("solvers.jl")
export tdma,tdma!


end # module MultiDiagonalMatrices
