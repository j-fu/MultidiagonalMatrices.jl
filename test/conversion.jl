using Test
using MultiDiagonalMatrices
using SparseArrays
using LinearAlgebra

function matrixtest(n,diags,blocksize)
    f=rand(blocksize*n)
    A=mdrand(n,diags;blocksize)
    @test A*f ≈ Matrix(A)*f
end

function matrixtest2(n,diags,blocksize)
    f=rand(blocksize,n)
    A=mdrand(n,diags;blocksize)
    @test vec(A*f) ≈ Matrix(A)*vec(f)
end

function matrixtest3(n,diags,blocksize)
    f=[@SVector rand(blocksize) for i=1:n]
    A=mdrand(n,diags;blocksize)
    @test reinterpret(Float64,A*f) ≈ Matrix(A)*reinterpret(Float64,f)
end


@testset "Matrix" begin
    matrixtest(10,[0],1)
    matrixtest(10,[1],1)
    matrixtest(10,[-1,0,1],1)
    matrixtest(10,[-5,-1,0,1],1)
    matrixtest(10,[-5,-1,0,1,5],1)

    matrixtest(10,[0],2)
    matrixtest(10,[1],2)
    matrixtest(10,[-1,0,1],2)
    matrixtest(10,[-5,-1,0,1],2)
    matrixtest(10,[-5,-1,0,1,5],2)

    matrixtest(50,[0],5)
    matrixtest(50,[1],5)
    matrixtest(50,[-1,0,1],5)
    matrixtest(50,[-5,-1,0,1],5)
    matrixtest(50,[-5,-1,0,1,5],5)

    matrixtest2(10,[0],2)
    matrixtest2(10,[1],2)
    matrixtest2(10,[-1,0,1],2)
    matrixtest2(10,[-5,-1,0,1],2)
    matrixtest2(10,[-5,-1,0,1,5],2)

    matrixtest2(50,[0],5)
    matrixtest2(50,[1],5)
    matrixtest2(50,[-1,0,1],5)
    matrixtest2(50,[-5,-1,0,1],5)
    matrixtest2(50,[-5,-1,0,1,5],5)

    matrixtest3(10,[0],2)
    matrixtest3(10,[1],2)
    matrixtest3(10,[-1,0,1],2)
    matrixtest3(10,[-5,-1,0,1],2)
    matrixtest3(10,[-5,-1,0,1,5],2)

    matrixtest3(50,[0],5)
    matrixtest3(50,[1],5)
    matrixtest3(50,[-1,0,1],5)
    matrixtest3(50,[-5,-1,0,1],5)
    matrixtest3(50,[-5,-1,0,1,5],5)


end


function sparsetest(n,diags)
    f=rand(n)
    A=mdrand(n,diags)
    @test A*f ≈ sparse(A)*f
end

@testset "Sparse" begin
    sparsetest(10,[0])
    sparsetest(10,[1])
    sparsetest(10,[-1,0,1])
    sparsetest(10,[-5,-1,0,1])
    sparsetest(10,[-5,-1,0,1,5])
end


function tridiagtest(n,diags)
    f=rand(n)
    A=mdrand(n,diags)
    @test A*f ≈ Tridiagonal(A)*f
end

@testset "Tridiag" begin
    tridiagtest(10,[0])
    tridiagtest(10,[1])
    tridiagtest(10,[-1,0,1])
end
