using Test
using MultiDiagonalMatrices
using SparseArrays
using LinearAlgebra
using StaticArrays


function setindextest(A)
    B=sparse(A)
    colptr=SparseArrays.getcolptr(B)
    rowval=SparseArrays.getrowval(B)
    for i=1:size(B,1)
        for k=colptr[i]:colptr[i+1]-1
            j=rowval[k]
            A[i,j]=1.0
        end
    end
    Bx=sparse(A)
    nonzeros(Bx)==ones(nnz(Bx))
end

@testset "setindex" begin
    @test setindextest(mdrand(10,10))
    @test setindextest(mdrand(60,60,blocksize=3))
    @test setindextest(mdrand(250,250,[-7,-1,0,1,7];blocksize=5))
end


