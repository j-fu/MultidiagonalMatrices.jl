struct MultiDiagonalMatrix{T,S} <: AbstractMatrix{T}
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


function MultiDiagonalMatrix(diags::Pair{Int64,Vector{T}}...) where T
    d=diags[1]
    n=abs(d.first)+length(d.second)
    if !all(d->abs(d.first)+length(d.second)==n,diags)
	error("MultiDiagonalMatrix(): diagonal length mismatch")
    end
    diags=sort([diags...],lt=(a,b)->isless(a.first,b.first))
    shadow=[d.first => _shadow(d.second) for d ∈ diags]
    MultiDiagonalMatrix(diags,shadow)
end

function MultiDiagonalMatrix(m::AbstractMatrix,diags::Integer...)
    MultiDiagonalMatrix([ i => diag(m,i) for i ∈ diags]...)
end

blocksize(md::MultiDiagonalMatrix{T,S}) where {T<: SMatrix,S} =size(md.diags[1].second[1],1)
blocksize(md::MultiDiagonalMatrix{T,S}) where {T<: Number,S} = 1

istridiagonal(m::MultiDiagonalMatrix)= first.(m.diags)==[-1,0,1]

function fdsizes(m::MultiDiagonalMatrix)
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


function mdrand(n,diags=[-1,0,1]; blocksize=1)
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
        for d in diags
            if d==0
                v=[SMatrix{blocksize,blocksize}(rand(-1:0.01:0,blocksize,blocksize) + Diagonal(rand(ndiags+1:0.01:ndiags+2,blocksize))) for i=1:n]
            else
                v=[SMatrix{blocksize,blocksize}(rand(-1:0.01:0.0,blocksize,blocksize)) for i=1:n-abs(d)]
            end
            push!(pairs,d=>v)
        end

    end
    MultiDiagonalMatrix(pairs...)
end

