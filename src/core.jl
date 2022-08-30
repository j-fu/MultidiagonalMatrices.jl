#
# This file contains the "hot loops" of the package.
#

"""
    $(TYPEDSIGNATURES)

Matrix-Vector multiplication.
"""
function LinearAlgebra.mul!(v::AbstractVector{Tv},
                            A::MultidiagonalMatrix{T},
                            u::AbstractVector{Tu}) where {
                                T<: Union{SMatrix,Number},
                                Tu<: Union{SVector,Number},
                                Tv<: Union{SVector,Number}}
    z=zero(Tv)
    for i∈eachindex(v)
        v[i]=z
    end
    for d ∈ A.diags
	o=d.first
	m=d.second 
	if o>=0
	    for i ∈ eachindex(m)
		v[i]+=m[i]*u[i+o] 
   	    end
	else
	    for i ∈ eachindex(m)
		v[i-o]+=m[i]*u[i] 
   	    end
	end
    end
    v
end


"""
    $(TYPEDSIGNATURES)

"LU  factorization" for number - returning number
This function allows to write generic algorithms for SMatrix and number.
"""
_lu(a::Number)=a

"""
    $(TYPEDSIGNATURES)

"LU  factorization" for matrix
"""
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

