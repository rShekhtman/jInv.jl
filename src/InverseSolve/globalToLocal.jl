export GlobalToLocal, getGlobalToLocal
export interpLocalToGlobal, interpGlobalToLocal
export prepareGlobalToLocal

abstract AbstractModel

type GlobalToLocal <: AbstractModel
	PForInv::SparseMatrixCSC # interpolation matrix from fwd mesh to inv mesh
	sigmaBackground::Vector{Float64} #  (# of cells fwd mesh)
end # type GlobalToLocal

# Constructors
getGlobalToLocal(P::SparseMatrixCSC) = GlobalToLocal(P,1e-8*ones(P.n))
getGlobalToLocal(P::SparseMatrixCSC,sigBack::Vector{Float64}) = GlobalToLocal(P,sigBack)
getGlobalToLocal(P::SparseMatrixCSC,sigBack::Vector{Float64},fname) = GlobalToLocal(P,sigBack)



function prepareGlobalToLocal(Mesh2Mesh::SparseMatrixCSC,Iact,sigmaBackground,fname)
	return getGlobalToLocal(Iact'*Mesh2Mesh,interpGlobalToLocal(sigmaBackground,Mesh2Mesh),fname)
end

function prepareGlobalToLocal(Mesh2Mesh::RemoteRef{Channel{Any}},Iact,sigmaBackground,fname)

	M2M = take!(Mesh2Mesh)
	
	G2L = prepareGlobalToLocal(M2M,Iact,sigmaBackground,fname)
   
	put!(Mesh2Mesh,false)
	return G2L
	
end

function prepareGlobalToLocal(Mesh2Mesh::RemoteRef{Channel{Any}},IactRef::RemoteRef{Channel{Any}},
	                           sigmaBackgroundRef::RemoteRef{Channel{Any}},fname)
	
	Iact = fetch(IactRef)
	sigmaBackground = fetch(sigmaBackgroundRef)
	return prepareGlobalToLocal(Mesh2Mesh,Iact,sigmaBackground,fname)
	
end

function prepareGlobalToLocal(Mesh2Mesh::Array{RemoteRef{Channel{Any}},1},Iact,sigmaBackground,fname="")
	pMod               = Array(RemoteRef{Channel{Any}},length(Mesh2Mesh))
	IactRef            = Array(RemoteRef{Channel{Any}},maximum(workers()))
	sigmaBackgroundRef = Array(RemoteRef{Channel{Any}},maximum(workers()))
	@sync begin
		for p=workers()
			@async begin
				# communicate Iact and sigmaBackground
				IactRef[p]            = remotecall_wait(p,identity,Iact)            # send active cells projection matrix
				sigmaBackgroundRef[p] = remotecall_wait(p,identity,sigmaBackground) # send background conductivity
				# call prepareGlobalToLocal with remote refs
				for idx=1:length(Mesh2Mesh)
					if p==Mesh2Mesh[idx].where
						pMod[idx] = remotecall_wait(p,prepareGlobalToLocal,Mesh2Mesh[idx],IactRef[p],sigmaBackgroundRef[p],fname)
					end
				end
			end
		end
	end
	return pMod	
end