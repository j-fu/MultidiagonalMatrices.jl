_lu(a::Number)=a
_lu(a::AbstractMatrix)=lu(a)

"""
    $(TYPEDSIGNATURES)

Tridiagonal matrix algorithm (метод прогонки) for both scalar and block cases.
In the block case, this uses the non-allocating lu factorization for `SMatrix`.

This algorithm does no pivoting, which in the context of nonliner solver iterations
and in the case of diagonally dominant matrices may of lesser concern. 
"""
function tdma!(u,a,b,c,f)
    N=length(f)
    alpha=similar(a)
    beta=similar(u)
    
    d=_lu(b[1])
    alpha[1] = - (d\c[1])
    beta[1]  = 	d\f[1]  										
    
    for i = 2:1:N-1 		             
	d=_lu(a[i-1]*alpha[i-1] + b[i]) 
	alpha[i] = - (d\c[i])
	beta[i]  = d\(f[i]-a[i-1]*beta[i-1])
    end

    d=_lu(a[N-1]*alpha[N-1]+b[N])
    u[N] = d\(f[N]-a[N-1]*beta[N-1])

    for i =N-1:-1:1 
	u[i] = alpha[i]*u[i+1]+beta[i]
    end		
    u
end


"""
    $(TYPEDSIGNATURES)

Tridiagonal matrix algorithm (метод прогонки) for tridiagonal matrix.
See [`tdma!(u,a,b,c,f)`](@ref) for more explanations.
"""
function tdma!(u,A::Tridiagonal{T},f) where T<:Number
    dl=diag(A,-1)
    d=diag(A,0)
    du=diag(A,1)
    tdma!(u,dl,d,du,f)
end


"""
    $(TYPEDSIGNATURES)

Tridiagonal matrix algorithm (метод прогонки) for scalar tridiagonal multdidiagonal matrix.
See [`tdma!(u,a,b,c,f)`](@ref) for more explanations.
"""
function tdma!(u,A::MultidiagonalMatrix{T},f) where T
    istridiagonal(A)|| error("tdma works only for tridiagonal matrices")
    dl=A.diags[1].second
    d=A.diags[2].second
    du=A.diags[3].second
    tdma!(u,dl,d,du,f)
end

"""
    $(TYPEDSIGNATURES)

Tridiagonal matrix solver for scalar matrices and vector of numbers.
See [`tdma!(u,a,b,c,f)`](@ref) for more explanations.
"""
function tdma(A::Union{Tridiagonal{T},MultidiagonalMatrix{T}},f) where T<:Number
    istridiagonal(A)|| error("tdma works only for tridiagonal matrices")
    Tu=promote_type(T,eltype(f))
    u=similar(f,Tu)
    tdma!(u,A,f)
end


"""
    $(TYPEDSIGNATURES)

Tridiagonal matrix solver for block matrices and vector of SVectors.
See [`tdma!(u,a,b,c,f)`](@ref) for more explanations.
"""
function tdma(md::MultidiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:SMatrix,Tu<:SVector}
    b=blocksize(md)
    if b!=size(u[1],1)
        error("incompatible blocksizes")
    end
    Tv=promote_type(eltype(T),eltype(Tu))
    v=Vector{SVector{b,Tv}}(undef,length(u))
    tdma!(v,md,u)
end


"""
    $(TYPEDSIGNATURES)

Tridiagonal matrix solver for block matrices and vector of numbers.
See [`tdma!(u,a,b,c,f)`](@ref) for more explanations.
"""
function tdma(md::MultidiagonalMatrix{T},u::AbstractVector{Tu}) where {T<:SMatrix,Tu<:Number}
    Tv=promote_type(eltype(T),Tu)
    v=Vector{Tv}(undef,length(u))
    b=blocksize(md)
    vv=reinterpret(SVector{b,Tv},v)
    uu=reinterpret(SVector{b,Tu},u)
    tdma!(vv,md,uu)
    v
end


"""
    $(TYPEDSIGNATURES)

Tridiagonal matrix solver for block matrices and matrix of numbers.
See [`tdma!(u,a,b,c,f)`](@ref) for more explanations.
"""
function tdma(md::MultidiagonalMatrix{T},u::AbstractMatrix{Tu}) where {T<:SMatrix,Tu<:Number}
    b=blocksize(md)
    if b!=size(u,1)
        error("incompatible blocksizes")
    end
    if size(md,1)!=length(u)
        error("incompatible matrix and vector sizes")
    end
    
    Tv=promote_type(eltype(T),Tu)
    v=Matrix{Tv}(undef,size(u)...)
    vv=vec(reinterpret(SVector{b,Tv},v))
    uu=vec(reinterpret(SVector{b,Tu},u))
    tdma!(vv,md,uu)
    v
end

