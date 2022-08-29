MultiDiagonalMatrices.jl
========================


Aims:
- Provide efficient structure for  matrices from finite volume/finite difference discretizations on structured grids
- Block version for systems
- Fallback to sparse direct solvers
- Conversions from/to SparseMatrixCSC
- Fully differentiable
- Compatibility to LinearSolve
- Solvers/preconditioners taking advantage of the structure:
  - TDMA for scalar and block tridiagonal matrices
  - (planned) ILU0, multigrid incl. block versions
