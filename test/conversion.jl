using Test
using MultiDiagonalMatrices
using SparseArrays
using LinearAlgebra
using StaticArrays

function matrixtest(n,diags,blocksize)
    f=rand(blocksize*n)
    A=mdrand(blocksize*n,blocksize*n,diags;blocksize)
    @test A*f ≈ Matrix(A)*f
end

function matrixtest2(n,diags,blocksize)
    f=rand(blocksize,n)
    A=mdrand(blocksize*n,blocksize*n,diags;blocksize)
    @test vec(A*f) ≈ Matrix(A)*vec(f)
end

function matrixtest3(n,diags,blocksize)
    f=[@SVector rand(blocksize) for i=1:n]
    A=mdrand(blocksize*n,blocksize*n,diags;blocksize)
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


function sparsetest(n,diags,blocksize)
    f=rand(blocksize*n)
    A=mdrand(blocksize*n,blocksize*n,diags;blocksize)
    @test A*f ≈ sparse(A)*f
end

function sparsetest2(n,diags,blocksize)
    f=rand(blocksize,n)
    A=mdrand(blocksize*n,blocksize*n,diags;blocksize)
    @test vec(A*f) ≈ sparse(A)*vec(f)
end

function sparsetest3(n,diags,blocksize)
    f=[@SVector rand(blocksize) for i=1:n]
    A=mdrand(blocksize*n,blocksize*n,diags;blocksize)
    @test reinterpret(Float64,A*f) ≈ sparse(A)*reinterpret(Float64,f)
end



@testset "Sparse" begin
    sparsetest(10,[0],1)
    sparsetest(10,[1],1)
    sparsetest(10,[-1,0,1],1)
    sparsetest(10,[-5,-1,0,1],1)
    sparsetest(10,[-5,-1,0,1,5],1)

    sparsetest(10,[0],2)
    sparsetest(10,[1],2)
    sparsetest(10,[-1,0,1],2)
    sparsetest(10,[-5,-1,0,1],2)
    sparsetest(10,[-5,-1,0,1,5],2)

    sparsetest(50,[0],5)
    sparsetest(50,[1],5)
    sparsetest(50,[-1,0,1],5)
    sparsetest(50,[-5,-1,0,1],5)
    sparsetest(50,[-5,-1,0,1,5],5)
   
    sparsetest(10,[0],2)
    sparsetest(10,[1],2)
    sparsetest(10,[-1,0,1],2)
    sparsetest(10,[-5,-1,0,1],2)
    sparsetest(10,[-5,-1,0,1,5],2)

    sparsetest(50,[0],5)
    sparsetest(50,[1],5)
    sparsetest(50,[-1,0,1],5)
    sparsetest(50,[-5,-1,0,1],5)
    sparsetest(50,[-5,-1,0,1,5],5)

    sparsetest2(10,[0],2)
    sparsetest2(10,[1],2)
    sparsetest2(10,[-1,0,1],2)
    sparsetest2(10,[-5,-1,0,1],2)
    sparsetest2(10,[-5,-1,0,1,5],2)

    sparsetest2(50,[0],5)
    sparsetest2(50,[1],5)
    sparsetest2(50,[-1,0,1],5)
    sparsetest2(50,[-5,-1,0,1],5)
    sparsetest2(50,[-5,-1,0,1,5],5)

    sparsetest3(10,[0],2)
    sparsetest3(10,[1],2)
    sparsetest3(10,[-1,0,1],2)
    sparsetest3(10,[-5,-1,0,1],2)
    sparsetest3(10,[-5,-1,0,1,5],2)

    sparsetest3(50,[0],5)
    sparsetest3(50,[1],5)
    sparsetest3(50,[-1,0,1],5)
    sparsetest3(50,[-5,-1,0,1],5)
    sparsetest3(50,[-5,-1,0,1,5],5)

end


function tridiagtest(n,diags)
    f=rand(n)
    A=mdrand(n,n,diags)
    @test A*f ≈ Tridiagonal(A)*f
end

@testset "Tridiag" begin
    tridiagtest(10,[0])
    tridiagtest(10,[1])
    tridiagtest(10,[-1,0,1])
end


@testset "to MD" begin
    
    A=mdrand(20,20,[-2,0,2])
    B=sparse(A)
    C=MultiDiagonalMatrix(B, [-2,0,2])
    @test A==C

    A=mdrand(20*3,20*3,[-2,0,2]; blocksize=3)
    B=sparse(A)
    C=MultiDiagonalMatrix(B, [-2,0,2]; blocksize=3)
    @test A==C

end
