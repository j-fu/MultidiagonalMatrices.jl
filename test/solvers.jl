using Test
using MultiDiagonalMatrices
using SparseArrays
using LinearAlgebra


@testset "Scalar tridiagonal solvers" begin
    n=20
    A=mdrand(n,[-1,0,1])
    f=rand(n)
    x=Tridiagonal(A)\f
    @test Matrix(A)\f≈x
    @test sparse(A)\f≈x
    @test A\f≈x
    @test tdma(A,f)≈x
    @test tdma(Tridiagonal(A),f)≈x

end
