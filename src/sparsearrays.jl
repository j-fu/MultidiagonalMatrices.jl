function SparseArrays.sparse(md::MultiDiagonalMatrix{T}) where T<: Number
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


function SparseArrays.sparse(md::MultiDiagonalMatrix{T}) where T<: SMatrix
    b=blocksize(md)
    I=Int64[]
    J=Int64[]
    Tv=eltype(T)
    V=Tv[]
    for d in md.diags
	o=d.first
	if o>=0
	    for i ∈ eachindex(d.second)
                ib=(i-1)*b
                jb=(i+o-1)*b
                for ii=1:b
                    for jj=1:b
		        push!(I,ib+ii)
		        push!(J,jb+jj)
		        push!(V,d.second[i][ii,jj])
                    end
                end
   	    end
	else
	    for i ∈ eachindex(d.second)
                ib=(i-o-1)*b
                jb=(i-1)*b
                for ii=1:b
                    for jj=1:b
		        push!(I,ib+ii)
		        push!(J,jb+jj)
		        push!(V,d.second[i][ii,jj])
                    end
                end
   	    end
	end
    end
    sparse(I,J,V,size(md)...)
end 
 
