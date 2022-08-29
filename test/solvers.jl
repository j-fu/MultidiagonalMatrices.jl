using Test
using MultidiagonalMatrices
using SparseArrays
using LinearAlgebra
using StaticArrays

@testset "TDMA scalar" begin
    n=20
    A=mdrand(n,n,[-1,0,1])
    x=rand(n)
    f=A*x
    
    @test Matrix(A)\f≈x
    @test sparse(A)\f≈x
    @test tdma(A,f)≈x
    @test tdma(Tridiagonal(A),f)≈x

end


@testset "TDMA block" begin
    n=20
    blocksize=5
    A=mdrand(n*blocksize,n*blocksize,[-1,0,1];blocksize)

    x=rand(blocksize,n)
    f=A*x
    @test tdma(A,f)≈x

    x=rand(blocksize*n)
    f=A*x
    @test tdma(A,f)≈x

    x=[@SVector rand(blocksize) for i=1:n]
    f=A*x
    @test tdma(A,f)≈x
end

@testset "backslash scalar" begin
    n=30
    A=mdrand(n,n,[-5, -1,0,1,7])
    x=rand(n)
    f=A*x
    
    @test Matrix(A)\f≈x
    @test sparse(A)\f≈x
    @test A\f≈x
    
end


@testset "backslash  block" begin
    n=30
    blocksize=5
    A=mdrand(n*blocksize,n*blocksize,[-5, -1,0,1,7];blocksize)


    x=rand(blocksize*n)
    f=A*x
    
    @test Matrix(A)\f≈x
    @test sparse(A)\f≈x
    @test A\f≈x


    x=rand(blocksize,n)
    f=A*x
    @test A\f≈x

    x=[@SVector rand(blocksize) for i=1:n]
    f=A*x
    @test A\f≈x
end
