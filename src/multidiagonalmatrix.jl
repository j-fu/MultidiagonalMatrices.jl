"""
    $(TYPEDEF)
Multidiagonal matrix structure. This is a generalization of a tridiagonal matrix structure
with an arbitrary number of diagonals. No entries outside these diagonals are stored, so this
is not a banded matrix.


### Fields
$(TYPEDFIELDS)
"""
struct MultidiagonalMatrix{T,S} <: AbstractMatrix{T}
    """
    Array of pairs, each  containing the index of a diagonal and its content
    """
    diags::Vector{Pair{Int64,Vector{T}}}

    """
    Array of pairs the index of a diagonal and a mutable view of its content.
    In case that `T` is an immutable `SMatrix`, the mutable view consists
    in a reinterpration of the diagonal content as a rank 3 tensor.
    """
    shadow::Vector{Pair{Int64,S}}
end



"""
    $(TYPEDSIGNATURES)

Create a multidiagonal matrix by specifying its diagonals.
"""
function MultidiagonalMatrix(diags::Pair{Int64,Vector{T}}...) where T
    d=diags[1]
    n=abs(d.first)+length(d.second)
    if !all(d->abs(d.first)+length(d.second)==n,diags)
	error("MultidiagonalMatrix(): diagonal length mismatch")
    end
    diags=sort([diags...],lt=(a,b)->isless(a.first,b.first))
    shadow=[d.first => _shadow(d.second) for d ∈ diags]
    MultidiagonalMatrix(diags,shadow)
end



"""
    $(TYPEDSIGNATURES)

Create a multidiagonal matrix by extracting diagonals from an AbstractMatrix.

"""
function MultidiagonalMatrix(A::AbstractMatrix{T},diags=[-1,0,1]; blocksize=1) where T<:Number
    MA=mdzeros(T,size(A)...,diags;blocksize)
    _setentries!(MA,A)
end

"""
    $(TYPEDSIGNATURES)
Return vector of number as its own shadow.
"""
_shadow(a::Vector{T}) where T<:Number = a

"""
    $(TYPEDSIGNATURES)
Return rank 3 tensor as a mutable view of a vector of `SMatrix`
"""
function _shadow(a::Vector{T}) where T<:SMatrix
    b=size(a[1],1)
    Te=eltype(a[1])
    n=length(a)
    reshape(reinterpret(Te,a),(b,b,n))
end


"""
    $(TYPEDSIGNATURES)

Set entries of a MultidiagonalMatrix from the corresponding ones of an abstract matrix
"""                                 
function _setentries!(MA::MultidiagonalMatrix{T},A::AbstractMatrix) where T<:Number
    for d in MA.diags
	o=d.first
	if o>=0
	    for i ∈ eachindex(d.second)
		d.second[i]=A[i,i+o]
   	    end
	else
	    for i ∈ eachindex(d.second)
		d.second[i]=A[i-o,i]
   	    end
	end
    end
    MA
end

"""
    $(TYPEDSIGNATURES)

Set entries of a MultidiagonalMatrix from the corresponding ones of an abstract matrix
"""                                 
function _setentries!(MA::MultidiagonalMatrix{T},A::AbstractMatrix) where T<:SMatrix
    b=blocksize(MA)
    for d in MA.shadow
	o=d.first
	if o>=0
	    for i = 1:size(d.second,3)
                ib=(i-1)*b
                jb=(i+o-1)*b
	        @views d.second[:,:,i].=A[ib+1:ib+b,jb+1:jb+b] 
   	    end
	else
	    for i = 1:size(d.second,3)
                ib=(i-o-1)*b
                jb=(i-1)*b
	        @views d.second[:,:,i].=A[ib+1:ib+b,jb+1:jb+b] 
   	    end
	end
    end
    MA
end


"""
    $(SIGNATURES)
 
Return the block size.
"""
blocksize(md::MultidiagonalMatrix{T,S}) where {T<: SMatrix,S} =size(md.diags[1].second[1],1)
blocksize(md::MultidiagonalMatrix{T,S}) where {T<: Number,S} = 1


"""
    $(SIGNATURES)

Detect if matrix is tridiagonal
"""
istridiagonal(m::MultidiagonalMatrix)= first.(m.diags)==[-1,0,1]

"""
    $(SIGNATURES)

Detect if matrix could come from a one- two or three dimensional grid
"""
function fdsizes(m::MultidiagonalMatrix)
    diags=first.(md.diags)
    ndiags=length(diags)
    if ndiags==3 && diags==[-1.0,1]
	n=length(diags[2])
	return (n,)
    elseif ndiags==5 && 
	diags[2:4]==[-1,0,1] && 
	diags[1]==-diags[5]
	n=length(diags[3])
        nx=diags[5]
	ny=n/nx
        return(nx,ny)
    elseif ndiags==7 && 
	ndiags[3:5]==[-1,0,1] && 
	diags[2]==-diags[6] &&
	diags[1]==-diags[7]
	n=length(diags[4])
        nxy=diags[7]
	nx=diags[5]
	ny=nxy/nx
	nz=n/nxy
	return (nx,ny,nz)
    end
    nothing
end

"""
    mdrand(n,m,diags=[-1,0,1]; blocksize=1)

Create a strictly diagonally dominant random (block) multidiagonal matrix. 
"""
function mdrand(n,m,diags=[-1,0,1]; blocksize=1)
    n==m || error("Can handle square matrices only")
    T=Float64
    ndiags=length(diags)*blocksize
    pairs=[]
    if blocksize==1
        for d in diags
            if d==0
                v=rand(ndiags+1:0.01:ndiags+2,n)
            else
                v=rand(-1:0.01:0.0,n-abs(d))
            end
            push!(pairs,d=>v)
        end
    else
        mod(n,blocksize)==0 || error("n=$n is no multiple of blocksize $blocksize")
        nb=n÷blocksize
        for d in diags
            if d==0
                v=[SMatrix{blocksize,blocksize}(rand(-1:0.01:0,blocksize,blocksize) + Diagonal(rand(ndiags+1:0.01:ndiags+2,blocksize))) for i=1:nb]
            else
                v=[SMatrix{blocksize,blocksize}(rand(-1:0.01:0.0,blocksize,blocksize)) for i=1:nb-abs(d)]
            end
            push!(pairs,d=>v)
        end

    end
    MultidiagonalMatrix(pairs...)
end

"""
    mdzeros(n,m,diags=[-1,0,1]; blocksize=1)

Create a zero block diagonal matrix with element type `Float64`
"""
function mdzeros(n,m,diags=[-1,0,1], ::Type{T}=Float64; blocksize=1) where T
    n==m || error("Can handle square matrices only")
    ndiags=length(diags)*blocksize
    pairs=[]
    if blocksize==1
        for d in diags
            push!(pairs,d=>zeros(T,n-abs(d)))
        end
    else
        mod(n,blocksize)==0 || error("n=$n is no multiple of blocksize $blocksize")
        nb=n÷blocksize
        for d in diags
            push!(pairs,d=>[@SMatrix zeros(T,blocksize,blocksize) for i=1:nb-abs(d)])
        end
    end
    MultidiagonalMatrix(pairs...)
end

"""
    mdzeros(T, n,m,diags=[-1,0,1]; blocksize=1)

Create a zero block diagonal matrix with element type `T`
"""
mdzeros(::Type{T},n,m,diags=[-1,0,1];blocksize=1)  where T =mdzeros(n,m,diags,T;blocksize)
