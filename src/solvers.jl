_lu(a::Number)=a
_lu(a::AbstractMatrix)=lu(a)

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


function tdma!(u,A::Union{Tridiagonal{T},MultiDiagonalMatrix{T}},f) where T<:Number
    istridiagonal(A)|| error("tdma works only for tridiagonal matrices")
    dl=diag(A,-1)
    d=diag(A,0)
    du=diag(A,1)
    tdma!(u,dl,d,du,f)
end

function tdma(A::Union{Tridiagonal{T},MultiDiagonalMatrix{T}},f) where T<:Number
    istridiagonal(A)|| error("tdma works only for tridiagonal matrices")
    Tu=promote_type(T,eltype(f))
    u=similar(f,Tu)
    tdma!(u,A,f)
end





# function tdma!(u,a::Vector{Tm},b::Vector{Tm},c::Vector{Tm},f) where Tm<:AbstractMatrix
#     N = length(b)
#     bsize=size(b[1],1)
#     alpha=similar(a)
    
    
#     beta  = similar(u)
#     d=lu(b[1])
#     alpha[1]  = - (d\c[1])	
#     @views beta[:,1]=d\f[:,1] 
    
#     for i = 2:1:N-1 
# 	d=lu(a[i-1]*alpha[i-1]+b[i])
# 	alpha[i] = - (d\c[i]) 
# 	@views beta[:,i]=d\(f[:,i]-a[i-1]*beta[:,i-1])
#     end
    
#     d=lu(a[N-1]*alpha[N-1]+b[N])
#     @views u[:,N] = d\(f[:,N]-a[N-1]*beta[:,N-1])
    
#     for i =N-1:-1:1 								
# 	@views u[:,i] = alpha[i]*u[:,i+1]+beta[:,i]				
#     end															
#     u 											
# end
