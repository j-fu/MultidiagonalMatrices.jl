"""
    $(TYPEDSIGNATURES)

Create tridiagonal matrix from multidiagonal matrix.
"""
function LinearAlgebra.Tridiagonal(md::MultidiagonalMatrix{T}) where T
	n=size(md,1)
	idl=findfirst(d->d.first==-1, md.diags)
	id=findfirst(d->d.first==0, md.diags)
	idu=findfirst(d->d.first==1, md.diags)
	dl= idl == nothing ? zeros(T,n-1) : md.diags[idl].second
	d=  id  == nothing ? zeros(T,n)   : md.diags[id].second
	du= idu == nothing ? zeros(T,n-1) : md.diags[idu].second
	Tridiagonal(dl,d,du)
end

istridiagonal(::LinearAlgebra.Tridiagonal)=true


"""
    $(TYPEDSIGNATURES)

Matrix-Vector multiplication.
"""
function LinearAlgebra.mul!(v::AbstractVector{Tv},md::MultidiagonalMatrix{T},u::AbstractVector{Tu}) where {T<: Union{SMatrix,Number}, Tu<: Union{SVector,Number},Tv<: Union{SVector,Number}}
    z=zero(Tv)
    for i∈eachindex(v)
        v[i]=z
    end
    for d ∈ md.diags
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


