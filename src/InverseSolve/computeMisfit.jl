export computeMisfit


"""
	Dc,F,dF,d2F = computeMisfit(...)
	
	Computes misfit for PDE parameter estimation problem.
	
	computeMisfit has several options.
"""
function computeMisfit(sigma,     # conductivity on inv mesh (active cells only)
                       gloc::GlobalToLocal,
                       misfit::Function,
                       pFor::ForwardProbType,
                       Dobs, Wd,  # observed data, 1/standard deviation
                       doDerivative::Bool=true, doClear::Bool=false)
	#
	#	computeMisfit for a single forward problem
	#
	#	Notes: Everything is stored in memory on the node executing this function.
	#          The gradient and sensitivities returned are with respect to the conductivity on the inv mesh.
	#
	# try
	times = zeros(4)
	
	tic(); 
		sigmaloc = interpGlobalToLocal(sigma,gloc.PForInv,gloc.sigmaBackground); 
	times[1]=toq()
	tic()
		Dc,pFor  = getData(sigmaloc,pFor)      # fwd model to get predicted data
	times[2]=toq()
	tic()
		F,dF,d2F = misfit(Dc,Dobs,Wd)
	times[3]=toq()
		
	if doDerivative
		tic()
		dF = interpLocalToGlobal(getSensTMatVec(dF,sigmaloc,pFor),gloc.PForInv)
		times[4]=toq()
	end

	# free memory
	if doClear
		clear!(pFor.Ainv)
	end

	return Dc,F,dF,d2F,pFor,times
end # function computeMisfit


function computeMisfit(sigmaRef::RemoteRef{Channel{Any}},
                        glocRef::RemoteRef{Channel{Any}},
                         misfit::Function,
                        pForRef::RemoteRef{Channel{Any}},
                        DobsRef::RemoteRef{Channel{Any}},
				      WdRef::RemoteRef{Channel{Any}},
				      dFRef::RemoteRef{Channel{Any}},
                  doDerivative,doClear::Bool=false)
	#
	#	computeMisfit for single forward problem
	#
	#	Note: model (including interpolation matrix) and forward problems are RemoteRefs
	#
	rrlocs = [glocRef.where pForRef.where DobsRef.where WdRef.where dFRef.where]
	if !all(rrlocs .== myid())
		warn("computeMisfit: Problem on worker $(myid()) not all remote refs are stored here, but rrlocs=$rrlocs")
	end
	
	sigma = fetch(sigmaRef)
	gloc  = fetch(glocRef)
	pFor  = take!(pForRef)
	Dobs  = fetch(DobsRef)
	Wd    = fetch(WdRef)
	
	Dc,F,dFi,d2F,pFor,times = computeMisfit(sigma,gloc,misfit,pFor,Dobs,Wd,doDerivative,doClear)
	
	put!(pForRef,pFor)
	# add to gradient
	if doDerivative
		dF = take!(dFRef)
		put!(dFRef,dF += dFi)
	end
	# put predicted data and d2F into remote refs (no need to communicate them)
	Dc  = remotecall(myid(),identity,Dc)
	d2F = remotecall(myid(),identity,d2F)

	return Dc,F,d2F,times
end


function computeMisfit(sigma,
	glocRef::Array{RemoteRef{Channel{Any}},1},
	misfit::Function,
	pForRef::Array{RemoteRef{Channel{Any}},1},
	DobsRef::Array{RemoteRef{Channel{Any}},1},
	WdRef::Array{RemoteRef{Channel{Any}},1},
	doDerivative::Bool=true,
	indCredit=collect(1:length(pForRef)))
	#
	#	computeMisfit for multiple forward problems
	#
	#	This method runs in parallel (iff nworkers()> 1 )
	#
	#	Note: ForwardProblems and Mesh-2-Mesh Interpolation are RemoteRefs
	#								(i.e. they are stored in memory of a particular worker).
	#
	
	F   = 0.0
	dF  = (doDerivative) ? zeros(length(sigma)) : []
	d2F = cell(length(pForRef));
	Dc  = Array(RemoteRef{Channel{Any}}, size(DobsRef))

	indDebit = []
	updateRes(Fi,idx) = (F+=Fi;push!(indDebit,idx))
	updateDF(x) = (dF+=x)
	
	# find out which workers are involved
	workerList = []
	for k=indCredit
		push!(workerList,pForRef[k].where)
	end
	workerList = unique(workerList)
	
	# send sigma to all workers
	sigRef = Array(RemoteRef{Channel{Any}},maximum(workers()))
	dFiRef = Array(RemoteRef{Channel{Any}},maximum(workers()))

	times = zeros(4);
	updateTimes(tt) = (times+=tt)
	
	@sync begin
		for p=workerList
			@async begin
				# communivate model and allocate RemoteRef for gradient
				sigRef[p] = remotecall(p,identity,sigma)   # send conductivity to workers
				dFiRef[p] = remotecall(p,identity,0.0) # get remote Ref to part of gradient
				# solve forward problems
				for idx=1:length(pForRef)
					if pForRef[idx].where==p
						Dc[idx],Fi,d2F[idx],tt =	remotecall_fetch(p,computeMisfit,sigRef[p],glocRef[idx],misfit,pForRef[idx],DobsRef[idx],WdRef[idx],dFiRef[p],doDerivative)
						updateRes(Fi,idx)
						updateTimes(tt)
					end
				end
				
				# sum up gradients 
				if doDerivative
					updateDF(fetch(dFiRef[p]))
				end
			end
		end
	end
	return Dc,F,dF,d2F,pForRef,times,indDebit
end


function computeMisfit(sigma,GLoc::Array,misfit::Function,PF::Array,Dobs::Array,Wd::Array,doDerivative::Bool=true,indCredit=collect(1:length(PF)))
	#
	#	computeMisfit for multiple forward problems
	#
	#	This method runs in parallel (iff nworkers()> 1 )
	#
	#	Note: ForwardProblems and Mesh-2-Mesh Interpolation are stored on the main processor
	#		  and then sent to a particular worker, which returns an updated pFor.
	#
	numFor   = length(PF)
 	F        = 0.0
    dF       = (doDerivative) ? zeros(length(sigma)) : []
 	d2F      = cell(numFor)
 	Dc       = cell(numFor)
	indDebit = []

	# draw next problem to be solved
	nextidx() = (idx = (isempty(indCredit)) ? -1 : pop!(indCredit))

 	updateRes(Fi,dFi,idx) = (F+=Fi; dF= (doDerivative)? dF+dFi : []; push!(indDebit,idx))

	times = zeros(4);
	updateTimes(tt) = (times+=tt)
	
 	@sync begin
 		for p = workers()
 				@async begin
 					while true
 						idx = nextidx()
 						if idx == -1
 							break
 						end
 							Dc[idx],Fi,dFi,d2F[idx],PF[idx],tt = 
 remotecall_fetch(p,computeMisfit,sigma,GLoc[idx],misfit,PF[idx],Dobs[idx],Wd[idx],doDerivative)
 							updateRes(Fi,dFi,idx)
							updateTimes(tt)
 					end
 				end
 		end
 	end
 	
 	return Dc,F,dF,d2F,PF,times,indDebit
 end