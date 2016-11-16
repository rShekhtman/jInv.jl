module Utils
	include("sparseUtils.jl")
	include("testing.jl")
	include("expandPolygon.jl")
	include("sortpermFast.jl")
	include("uniqueidx.jl")
	include("initRemoteChannel.jl")



  import Base.clear!
	export clear!

	function Base.sub2ind(n::Array{Int64,1},ii::Array{Int64,1},jj::Array{Int64,1},kk::Array{Int64,1})
		return Base.sub2ind((n[1],n[2],n[3]),ii,jj,kk)
	end

	function Base.sub2ind(n::Array{Int64,1},ii::Int64,jj::Int64,kk::Int64)
		return Base.sub2ind((n[1],n[2],n[3]),ii,jj,kk)
	end

	function clear!(R::RemoteChannel)
		p = take!(R)
		p = clear!(p)
	end

	function clear!(F::Future)
		p = fetch(F)
		p = clear!(p)
  end

	function clear!(PF::Union{Array{RemoteChannel},Array{Future}})
		@sync begin
			for p=workers()
				@async begin
					for i=1:length(PF)
						if p==PF[i].where
							remotecall(clear!, p, PF[i])
						end
					end
				end
			end
		end
	end


	function clear!{T,N}(x::Array{T,N})
		return Array(T,ntuple((i)->0, N))
	end

	function clear!{T}(x::Vector{T})
		return Array(T,0)
	end

	function clear!{T}(A::SparseMatrixCSC{T})
		return spzeros(0,0);
	end

	export getWorkerIds
	function getWorkerIds(A::Array{RemoteChannel})
		Ids = []
		for k=1:length(A)
			push!(Ids,A[k].where)
		end
		return unique(Ids)
	end
end
