struct MultidiagonalMatrix{T,S} <: AbstractMatrix{T}
    diags::Vector{Pair{Int64,Vector{T}}}
    shadow::Vector{Pair{Int64,S}}
end

_shadow(a::Vector{T}) where T<:Number = a

function _shadow(a::Vector{T}) where T<:SMatrix
    b=size(a[1],1)
    Te=eltype(a[1])
    n=length(a)
    reshape(reinterpret(Te,a),(b,b,n))
end


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


function MultidiagonalMatrix(A::AbstractMatrix{T},diags=[-1,0,1]; blocksize=1) where T<:Number
    MA=mdzeros(T,size(A)...,diags;blocksize)
    _setentries!(MA,A)
end

                                 
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



blocksize(md::MultidiagonalMatrix{T,S}) where {T<: SMatrix,S} =size(md.diags[1].second[1],1)
blocksize(md::MultidiagonalMatrix{T,S}) where {T<: Number,S} = 1

istridiagonal(m::MultidiagonalMatrix)= first.(m.diags)==[-1,0,1]

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

mdzeros(::Type{T},n,m,diags=[-1,0,1];blocksize=1)  where T =mdzeros(n,m,diags,T;blocksize)
