## The multidiagonal matrix type
```@docs
MultidiagonalMatrix
```

## Constructors
```@docs
MultidiagonalMatrix(diags::Pair{Int64,Vector{T}}...) where T
MultidiagonalMatrix(A::AbstractMatrix{T},diags=[-1,0,1]; blocksize=1) where T<:Number
```

```@example
using MultidiagonalMatrices # hide
a = rand(5)
b = rand(6)
c = rand(4)
MultidiagonalMatrix( -1 => a, 0 => b, 2 => c)
```


```@example
using MultidiagonalMatrices # hide
A=rand(12,12)
MultidiagonalMatrix(A,[-1,0,2];blocksize=2)
```

```@docs
mdrand
```
```@example
using MultidiagonalMatrices # hide
mdrand(6,6,[-1,0,2])
```

```@docs
mdzeros
```

## Solvers
```@docs
tdma!
tdma
```


## Utilities
```@docs
blocksize
istridiagonal
fdsizes
```

## Base methods
```@docs
Base.eltype
Base.size
Base.getindex
Base.setindex!
Base.Matrix
Base.:*
Base.:\
```


## LinearAlgebra methods
```@docs
LinearAlgebra.Tridiagonal
LinearAlgebra.mul!
```

## SparseArray methods
```@docs
SparseArrays.sparse
SparseArrays.nnz
```

