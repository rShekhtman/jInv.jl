export getData



function getData(sigma::RemoteRef{Channel{Any}},pFor::RemoteRef{Channel{Any}},Mesh2Mesh,doClear::Bool=false)
	#= 
		load a forward problem from RemoteRef
	=#
	# get mesh 2 mesh interpolation matrix
	sig = interpGlobalToLocal(fetch(sigma),fetch(Mesh2Mesh))
	pF  = take!(pFor)
	dobs,pF   = getData(sig,pF,doClear)
	put!(pFor,pF)
	Dobs  = remotecall(myid(),identity,dobs)
	return Dobs,pFor
end


function getData(sigma::Vector,pFor::Array{ForwardProbType},Mesh2Mesh::Array=ones(length(pFor)),doClear::Bool=false,workerList=workers())
	#= 
		load and solve forward problems in parallel
	=#
	i=1; nextidx() = (idx = i; i+=1; idx)
	## Compute Dobs
	Dobs = Array(Any,length(pFor))
	workerList = intersect(workers(),workerList)
	if isempty(workerList)
		error("getData: workers do not exist!")
	end
	@sync begin
		for p=workerList
			@async begin
				while true
					idx = nextidx()
					if idx > length(pFor)
						break
					end
					P = Mesh2Mesh[idx]
					if isa(P,SparseMatrixCSC) && eltype(P.nzval)==Int16
						P = SparseMatrixCSC(P.m,P.n,P.colptr,P.rowval,2.^float(3*P.nzval))
					end
					Dobs[idx],pFor[idx]    = remotecall_fetch(p, getData,P'*sigma,pFor[idx],doClear)
				end
			end
		end
	end
	return Dobs,pFor
end

function getData(sigma::Vector,pFor::Array{RemoteRef{Channel{Any}}},Mesh2Mesh::Array=ones(length(pFor)),doClear::Bool=false)
	#= 
		load and solve forward problems in parallel
	=#
	## Compute Dobs
	Dobs = Array(RemoteRef{Channel{Any}},length(pFor))
	
	# find out which workers are involved
	workerList = []
	for k=1:length(pFor)
		push!(workerList,pFor[k].where)
	end
	workerList = unique(workerList)
	# prepare remoterefs for sigma
	sigmaRef = Array(RemoteRef{Channel{Any}},maximum(workers()))
	@sync begin
		for p=workerList
			@async begin
				# send model to worker
				sigmaRef[p] = remotecall(p,identity,sigma)
				# solve forward problems
				for idx=1:length(pFor)
					if p==pFor[idx].where
						Dobs[idx],pFor[idx]    = remotecall_fetch(p, getData,sigmaRef[p],pFor[idx],Mesh2Mesh[idx],doClear)
					end
				end
			end
		end
	end
	return Dobs,pFor
end