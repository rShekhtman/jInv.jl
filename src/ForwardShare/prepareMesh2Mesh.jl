export prepareMesh2Mesh

function prepareMesh2Mesh(pF::ForwardProbType, Minv::AbstractMesh, compact::Bool=true)

	P    = getInterpolationMatrix(Minv,pF.Mesh)'
	if compact
		Ps = SparseMatrixCSC(P.m,P.n,round(UInt32,P.colptr),round(UInt32,P.rowval),round(Int8,log2(P.nzval)/3))
	else
		Ps = SparseMatrixCSC(P.m,P.n,round(UInt32,P.colptr),round(UInt32,P.rowval),P.nzval)
	end
	return Ps
	
end

function prepareMesh2Mesh(pFor::RemoteRef{Channel{Any}}, MinvRef::RemoteRef{Channel{Any}}, compact::Bool=true)
	pF   = fetch(pFor)
	Minv = fetch(MinvRef)
	return prepareMesh2Mesh(pF, Minv, compact)
end

function prepareMesh2Mesh(pFor::Array{RemoteRef{Channel{Any}}},Minv::AbstractMesh,compact::Bool=true)

	Mesh2Mesh = Array(RemoteRef{Channel{Any}},length(pFor))

	# find out which workers are involved
	workerList = []
	for k=1:length(pFor)
		push!(workerList,pFor[k].where)
	end
	workerList = unique(workerList)
	# send sigma to all workers
	MinvRef = Array(RemoteRef{Channel{Any}},maximum(workers()))
	
	tic()
	@sync begin
		for p=workerList
			@async begin
				MinvRef[p] = remotecall_wait(p,identity,Minv)   # send model to workers
			end
		end
	end
	sendTime=toq()
	println("Time for sending out Minv $sendTime")
		
	@sync begin
		for p=workerList
			@async begin
				for idx=1:length(pFor)
					if p==pFor[idx].where
						Mesh2Mesh[idx] = remotecall_wait(p,prepareMesh2Mesh,pFor[idx],MinvRef[p],compact)
					end
				end
			end
		end
	end
	return Mesh2Mesh
end