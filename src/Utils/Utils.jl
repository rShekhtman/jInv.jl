module Utils
	include("sparseUtils.jl")
	include("testing.jl")
	include("expandPolygon.jl")
	include("sortpermFast.jl")
	include("uniqueidx.jl")
	
	
	
	
	export clear,clear!
	
	function clear!(T)
	end

	function Base.sub2ind(n::Array{Int64,1},ii::Array{Int64,1},jj::Array{Int64,1},kk::Array{Int64,1})	
		return Base.sub2ind((n[1],n[2],n[3]),ii,jj,kk)
	end

	function Base.sub2ind(n::Array{Int64,1},ii::Int64,jj::Int64,kk::Int64)	
		return Base.sub2ind((n[1],n[2],n[3]),ii,jj,kk)
	end
	  		
	function clear!(R::RemoteRef{Channel{Any}})
		p = take!(R)
		clear!(p)
		put(p,1.0)
	end

	function clear!(PF::Array{RemoteRef{Channel{Any}}})
		@sync begin
			for p=workers()
				@async begin
					for i=1:length(PF)
						if p==PF[i].where
							remotecall(p, clear!, PF[i])
						end
					end
				end
			end
		end
	end
	
	
	function clear{T,N}(x::Array{T,N})
		return Array(T,ntuple((i)->0, N))
	end
	
	function clear{T}(x::Vector{T})
		return Array(T,0)
	end
	
	function clear{T}(A::SparseMatrixCSC{T})
		return spzeros(T,size(A,1),size(A,2));
	end	
end
