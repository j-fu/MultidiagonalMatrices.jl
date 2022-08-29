function SparseArrays.sparse(md::MultiDiagonalMatrix{T}) where T
    I=Int64[]
    J=Int64[]
    V=T[]
    for d in md.diags
	o=d.first
	if o>=0
	    for i ∈ eachindex(d.second)
		push!(I,i)
		push!(J,i+o)
		push!(V,d.second[i])
   	    end
	else
	    for i ∈ eachindex(d.second)
		push!(I,i-o)
		push!(J,i)
		push!(V,d.second[i])
   	    end
	end
    end
    sparse(I,J,V,size(md)...)
end 
 
