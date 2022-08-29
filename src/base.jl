function Base.size(m::MultiDiagonalMatrix)
    d=m.diags[1]
    bs=blocksize(m)
    n=abs(d.first)+length(d.second)
    (n*bs,n*bs)
end

Base.size(m::MultiDiagonalMatrix,i)=size(m)[i]

function Base.getindex(m::MultiDiagonalMatrix{T}, i, j) where T<: Number
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


function _blockindices(m::MultiDiagonalMatrix{T,S},i,j) where {T<: SMatrix,S}
    bs=blocksize(m)
    ib=(i-1)÷bs+1
    jb=(j-1)÷bs+1
    ii=mod(i-1,bs)+1
    jj=mod(j-1,bs)+1
    ib,jb,ii,jj
end

function Base.getindex(m::MultiDiagonalMatrix{T,S}, i, j) where {T<: SMatrix,S}
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



function Base.setindex!(m::MultiDiagonalMatrix{T}, v, i, j) where T<: Number
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


function Base.setindex!(m::MultiDiagonalMatrix{T}, v, i, j) where T<: SMatrix
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

function _setindex!(m::MultiDiagonalMatrix{T}, v, i, j) where T<: SMatrix
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




function Base.Matrix(A::MultiDiagonalMatrix{T}) where T<:Number
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

function Base.Matrix(A::MultiDiagonalMatrix{T}) where T<:SMatrix
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


function Base.:*(A::MultiDiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:Number,Tu<:Number}
    Tv=promote_type(T,Tu)
    v=Vector{Tv}(undef,length(u))
    LinearAlgebra.mul!(v,A,u)
end

function Base.:*(A::MultiDiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:SMatrix,Tu<:SVector}
    b=blocksize(A)
    if b!=size(u[1],1)
        error("incompatible blocksizes")
    end
    Tv=promote_type(eltype(T),eltype(Tu))
    v=Vector{SVector{b,Tv}}(undef,length(u))
    LinearAlgebra.mul!(v,A,u)
end


function Base.:*(A::MultiDiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:SMatrix,Tu<:Number}
    Tv=promote_type(eltype(T),Tu)
    v=Vector{Tv}(undef,length(u))
    b=blocksize(A)
    vv=reinterpret(SVector{b,Tv},v)
    uu=reinterpret(SVector{b,Tu},u)
    LinearAlgebra.mul!(vv,A,uu)
    v
end


function Base.:*(A::MultiDiagonalMatrix{T},u::AbstractMatrix{Tu}) where {T<:SMatrix,Tu<:Number}
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



function Base.:\(A::MultiDiagonalMatrix{T},f::AbstractVecOrMat) where T<:Number
    if istridiagonal(A)
        tdma(A,f)
    else
        sparse(A)\f
    end
end

function Base.:\(A::MultiDiagonalMatrix{T},f::AbstractVector{Tf}) where {T<:SMatrix,Tf<:SVector}
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


function Base.:\(A::MultiDiagonalMatrix{T},f::AbstractVector{Tf}) where {T<:SMatrix,Tf<:Number}
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


function Base.:\(A::MultiDiagonalMatrix{T},f::AbstractMatrix{Tf}) where {T<:SMatrix,Tf<:Number}
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

