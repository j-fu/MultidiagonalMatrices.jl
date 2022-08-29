"""
    $(TYPEDSIGNATURES)

Matrix element type.
"""
Base.eltype(A::MultidiagonalMatrix)=eltype(eltype(A.diags[1].second[1]))


"""
    $(TYPEDSIGNATURES)

Matrix size.
"""
function Base.size(m::MultidiagonalMatrix)
    d=m.diags[1]
    bs=blocksize(m)
    n=abs(d.first)+length(d.second)
    (n*bs,n*bs)
end

"""
    $(TYPEDSIGNATURES)

Matrix size in one dimension
"""
Base.size(m::MultidiagonalMatrix,i)=size(m)[i]

"""
    $(TYPEDSIGNATURES)

Matrix value at index i,j for scalar matrix
"""
function Base.getindex(m::MultidiagonalMatrix{T}, i, j) where T<: Number
    o=j-i
    idiag=findfirst(d->d.first==o, m.diags)
    if idiag==nothing
	0.0
    elseif o>=0
	m.diags[idiag].second[j-o]
    else
	m.diags[idiag].second[i+o]
    end
end

"""
    $(TYPEDSIGNATURES)

Return block numbers and indices within  block
"""
function _blockindices(m::MultidiagonalMatrix{T,S},i,j) where {T<: SMatrix,S}
    bs=blocksize(m)
    ib=(i-1)÷bs+1
    jb=(j-1)÷bs+1
    ii=mod(i-1,bs)+1
    jj=mod(j-1,bs)+1
    ib,jb,ii,jj
end


"""
    $(TYPEDSIGNATURES)


Matrix value at index i,j for block matrix
"""
function Base.getindex(m::MultidiagonalMatrix{T,S}, i, j) where {T<: SMatrix,S}
    ib,jb,ii,jj=_blockindices(m,i,j)
    o=jb-ib
    idiag=findfirst(d->d.first==o, m.diags)
    if idiag==nothing
	0.0
    elseif o>=0
	m.diags[idiag].second[jb-o][ii,jj]
    else
	m.diags[idiag].second[ib+o][ii,jj]
    end
end


"""
    $(TYPEDSIGNATURES)

Set matrix value at index i,j for scalar matrix
"""
function Base.setindex!(m::MultidiagonalMatrix{T}, v, i, j) where T<: Number
    o=j-i
    idiag=findfirst(d->d.first==o, m.diags)
    if idiag==nothing
	throw(ArgumentError("cannot set entry ($i, $j) outside of a defined diagonal to a nonzero value ($v)"))
    elseif o>=0
	m.diags[idiag].second[j-o]=v
    else
	m.diags[idiag].second[i+o]=v
    end
end


"""
    $(TYPEDSIGNATURES)

Set matrix value at index i,j for block matrix. Throws error when
attempting to set value outside of defined diagonals.
"""
function Base.setindex!(m::MultidiagonalMatrix{T}, v, i, j) where T<: SMatrix
    ib,jb,ii,jj=_blockindices(m,i,j)
    o=jb-ib
    idiag=findfirst(d->d.first==o, m.diags)
    if idiag==nothing
	throw(ArgumentError("cannot set entry ($i, $j) outside of a defined block diagonal to a nonzero value ($v)"))
    elseif o>=0
	m.shadow[idiag].second[ii,jj,jb-o]=v
    else
	m.shadow[idiag].second[ii,jj,ib+o]=v
    end
end

"""
    $(TYPEDSIGNATURES)

Set matrix value at index i,j for block matrix. Quietly ingnores
attempt to set value outside of defined diagonals.
"""
function _setindex!(m::MultidiagonalMatrix{T}, v, i, j) where T<: SMatrix
    ib,jb,ii,jj=_blockindices(m,i,j)
    o=jb-ib
    idiag=findfirst(d->d.first==o, m.diags)
    if idiag==nothing
	return 0.0
    elseif o>=0
	m.shadow[idiag].second[ii,jj,jb-o]=v
    else
	m.shadow[idiag].second[ii,jj,ib+o]=v
    end
end


"""
    $(TYPEDSIGNATURES)

Create a dense matrix from scalar multidiagonal matrix.
"""
function Base.Matrix(A::MultidiagonalMatrix{T}) where T<:Number
    m=zeros(T,size(A)...)
    for d in A.diags
	o=d.first
	if o>=0
	    for i ∈ eachindex(d.second)
		m[i,i+o]=d.second[i]
   	    end
	else
	    for i ∈ eachindex(d.second)
		m[i-o,i]=d.second[i]
   	    end
	end
    end
    m
end

"""
    $(TYPEDSIGNATURES)

Create a dense matrix from multidiagonal block matrix
"""
function Base.Matrix(A::MultidiagonalMatrix{T}) where T<:SMatrix
    b=blocksize(A)
    m=zeros(eltype(T),size(A)...)
    for d in A.diags
	o=d.first
	if o>=0
	    for i ∈ eachindex(d.second)
                ib=(i-1)*b
                jb=(i+o-1)*b
		m[ib+1:ib+b,jb+1:jb+b] = d.second[i]
   	    end
	else
	    for i ∈ eachindex(d.second)
                ib=(i-o-1)*b
                jb=(i-1)*b
	        m[ib+1:ib+b,jb+1:jb+b] = d.second[i]
   	    end
	end
    end
    m
end


"""
    $(TYPEDSIGNATURES)

Multiply vector of number by scalar multidiagonal matrix.
"""
function Base.:*(A::MultidiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:Number,Tu<:Number}
    Tv=promote_type(T,Tu)
    v=Vector{Tv}(undef,length(u))
    LinearAlgebra.mul!(v,A,u)
end

"""
    $(TYPEDSIGNATURES)

Multiply vector of SVectors by multidiagonal block matrix.
"""
function Base.:*(A::MultidiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:SMatrix,Tu<:SVector}
    b=blocksize(A)
    if b!=size(u[1],1)
        error("incompatible blocksizes")
    end
    Tv=promote_type(eltype(T),eltype(Tu))
    v=Vector{SVector{b,Tv}}(undef,length(u))
    LinearAlgebra.mul!(v,A,u)
end


"""
    $(TYPEDSIGNATURES)

Multiply vector of numbers by multidiagonal block matrix.
"""
function Base.:*(A::MultidiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:SMatrix,Tu<:Number}
    Tv=promote_type(eltype(T),Tu)
    v=Vector{Tv}(undef,length(u))
    b=blocksize(A)
    vv=reinterpret(SVector{b,Tv},v)
    uu=reinterpret(SVector{b,Tu},u)
    LinearAlgebra.mul!(vv,A,uu)
    v
end


"""
    $(TYPEDSIGNATURES)

Multiply matrix of numbers by multidiagonal block matrix.
"""
function Base.:*(A::MultidiagonalMatrix{T},u::AbstractMatrix{Tu}) where {T<:SMatrix,Tu<:Number}
    b=blocksize(A)
    if b!=size(u,1)
        error("incompatible blocksizes")
    end
    if size(A,1)!=length(u)
        error("incompatible matrix and vector sizes")
    end
    
    Tv=promote_type(eltype(T),Tu)
    v=Matrix{Tv}(undef,size(u)...)
    vv=vec(reinterpret(SVector{b,Tv},v))
    uu=vec(reinterpret(SVector{b,Tu},u))
    LinearAlgebra.mul!(vv,A,uu)
    v
end



"""
    $(TYPEDSIGNATURES)

Solve scalar matrix problem for vector of numbers
"""
function Base.:\(A::MultidiagonalMatrix{T},f::AbstractVecOrMat) where T<:Number
    if istridiagonal(A)
        tdma(A,f)
    else
        sparse(A)\f
    end
end

"""
    $(TYPEDSIGNATURES)

Solve block matrix problem for vector of SVectors
"""
function Base.:\(A::MultidiagonalMatrix{T},f::AbstractVector{Tf}) where {T<:SMatrix,Tf<:SVector}
    b=blocksize(A)
    if b!=size(f[1],1)
        error("incompatible blocksizes")
    end
    if istridiagonal(A)
        tdma(A,f)
    else
        ff=reinterpret(eltype(Tf),f)
        uu=sparse(A)\ff
        reinterpret(SVector{b,eltype(uu)},uu)
    end
end


"""
    $(TYPEDSIGNATURES)

Solve block matrix problem for vector of numbers
"""
function Base.:\(A::MultidiagonalMatrix{T},f::AbstractVector{Tf}) where {T<:SMatrix,Tf<:Number}
    if size(A,1)!=length(f)
        @show size(A,1) ,length(f)
        error("incompatible sizes")
    end
    if istridiagonal(A)
        tdma(A,f)
    else
        sparse(A)\f
    end
end


"""
    $(TYPEDSIGNATURES)

Solve block matrix problem for matrix
"""
function Base.:\(A::MultidiagonalMatrix{T},f::AbstractMatrix{Tf}) where {T<:SMatrix,Tf<:Number}
    b=blocksize(A)
    if b!=size(f,1)
        error("incompatible blocksizes")
    end
    if istridiagonal(A)
        tdma(A,f)
    else
        uu=sparse(A)\vec(f)
        reshape(uu,size(f))
    end
end

