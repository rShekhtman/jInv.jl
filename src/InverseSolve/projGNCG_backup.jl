export projGNCG

type projGNCGhis
	Jc::Array
	F::Array
	Dc::Array
	Rc::Array
	alphas::Array
	Active::Array
	stepnorm::Array
	timeMisfit::Array
	timeReg::Array
	timePCG::Array
	hisPCG::Array
	timeGradMisfit::Array
end

function getProjGNCGhis(maxIter,maxIterCG)
	Jc = zeros(maxIter+1)
	F = zeros(maxIter+1)
	Dc = []
	Rc = zeros(maxIter+1)
	alphas = zeros(maxIter+1)
	Active = zeros(maxIter+1)
	stepnorm = zeros(maxIter+1)
	timeMisfit = zeros(maxIter+1,4)
	timeReg = zeros(maxIter+1)
	timePCG = zeros(maxIter+1,1)
	hisPCG = []
	timeGradMisfit = zeros(maxIter+1,2)

	return projGNCGhis(Jc,F,Dc,Rc,alphas,Active,stepnorm,timeMisfit,timeReg,timePCG,hisPCG,timeGradMisfit)
end

function updateHis!(iter::Int64,His::projGNCGhis,Jc::Real,Fc,Dc,Rc::Real,alpha::Real,nActive::Int64,stepnorm::Real,timeMisfit::Vector,timeReg::Real)
	His.Jc[iter+1] = Jc
	His.F[iter+1] = Fc
	push!(His.Dc,Dc)
	His.Rc[iter+1] = Rc
	His.alphas[iter+1] = alpha
	His.Active[iter+1] = nActive
	His.stepnorm[iter+1] = stepnorm
	His.timeMisfit[iter+1,:]+=timeMisfit'
	His.timeReg[iter+1] += timeReg
end


"""
	mc,Dc,outerFlag = projGNCG(mc,pInv::InverseParam,PF,indFor = [])
	
	(Projected) Gauss-NewtonCG
	
	Input:
	
		mc::Vector          - intial guess for model
		pInv::InverseParam  - parameter for inversion
		PF                  - forward problem(s)
		indCredit           - indices of forward problems to work on
	
	Output:
		mc                  - final model
		Dc                  - data
		outerFlag           - flag for convergence
	
"""
function  projGNCG(mc,pInv::InverseParam,PF,indCredit=[],filenamePrefix=[])
#[mc,dpre] = projNewtonCG(mc,param)
#  
	maxIter     = pInv.maxIter      #  Max. no. iterations.
	pcgMaxIter  = pInv.pcgMaxIter   #  Max cg iters.
	pcgTol      = pInv.pcgTol       #  CG stopping tolerance.
	stepTol     = pInv.minUpdate    #  Step norm stopping tol.
	maxStep     = pInv.maxStep
	low         = pInv.boundsLow
	high        = pInv.boundsHigh
	alpha       = pInv.alpha
	Dobs        = pInv.dobs
	Wd          = pInv.Wd

	His = getProjGNCGhis(maxIter,pcgMaxIter)
	#---------------------------------------------------------------------------
	#  Initialization.
	#---------------------------------------------------------------------------
	
	Active = (mc .<=low) | (mc.>=high)  # Compute active set
	
	
	## evaluate function and derivatives
	sig,dsig = pInv.modelfun(mc)
	if isempty(indCredit)
		Dc,F,dF,d2F,PF,tMis = computeMisfit(sig,pInv.model,pInv.misfit,PF,Dobs,Wd,true)
	else
		Dc,F,dF,d2F,PF,tMis,indDebit = computeMisfit(sig,pInv.model,pInv.misfit,PF,Dobs,Wd,true,indCredit)
	end
	dF = dsig'*dF
	
	
	# compute regularizer
	tic()
	R,dR,d2R = pInv.regularizer(mc,pInv.mref,pInv.MInv,Iact=pInv.Iact,C=pInv.regparam)
	tReg = toq()    
	
	# objective function
	Jc  = F  + alpha*R
	gc  = dF + alpha*dR
	
	F0 = F; J0 = Jc
	############################################################################
	##  Outer iteration.                                                        #
	############################################################################
	iter = 0
	outerFlag = 0

	
	outStr = @sprintf("%4s\t%08s\t%08s\t%08s\t%08s\t%08s\n", 
					  "i.LS", "F", "R","alpha","Jc/J0","#Active")
	updateHis!(0,His,Jc,F,Dc,R,alpha,countnz(Active),0.0,tMis,tReg)
	
	print(outStr)
	f = open("jInv.out", "w")
	write(f, outStr)
	close(f)

	while outerFlag == 0
		
		iter += 1
		outStr = @sprintf("%3d.0\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3d\n", 
		         iter, F, R,alpha,Jc/J0,countnz(Active))
		print(outStr)
		f = open("jInv.out", "a")
		write(f, outStr)
		close(f)

		
		#  Set up Hessian and preconditioner.
		if isempty(indCredit)
			Hs(x) = dsig'*HessMatVec(dsig*x,PF,sig,pInv.model,d2F) + alpha*(d2R*x); 
		else
			Hs(x) = dsig'*HessMatVec(dsig*x,PF,sig,pInv.model,d2F,indDebit) + alpha*(d2R*x);
		end
		dd = 1./diag(d2R)
		xt  = zeros(length(mc))
		SSOR(r) = (xt[:]=0.0; xt=ssorPrecTrans!(d2R,xt,r,dd); return xt)
		
		##  inner CG iterations.
		tic()
		delm,hisPCG = projPCG(Hs,gc,Active,SSOR,pcgTol,pcgMaxIter)
		His.timePCG[iter+1]+= toq()
		push!(His.hisPCG,hisPCG)

		# scale step 
		if maximum(abs(delm)) > maxStep; delm = delm./maximum(abs(delm))*maxStep; end
		
		# take gradient direction in the active cells
		ga = -gc[mc .== low]
		if !isempty(ga)
			if maximum(abs(ga)) > maximum(abs(delm)); ga = ga./maximum(abs(ga))*maximum(abs(delm)); end
			delm[mc .== low] = ga
		end
		ga = -gc[mc .== high]
		if !isempty(ga)
			if maximum(abs(ga)) > maximum(abs(delm)); ga = ga./maximum(abs(ga))*maximum(abs(delm)); end
			delm[mc .== high] = ga
		end
		
		## Begin projected Armijo line search
		muLS = 1; lsIter = 1; mt = zeros(size(mc)); Jt = Jc
		while true
			mt = mc + muLS*delm
			mt[mt.<low]  = low[mt.<low]
			mt[mt.>high] = high[mt.>high]
			## evaluate function 
			sigt, = pInv.modelfun(mt)
			if isempty(indCredit)
				Dc,F,dF,d2F,PF,tMis = computeMisfit(sigt,pInv.model,pInv.misfit,PF,Dobs,Wd,false)
			else
				Dc,F,dF,d2F,PF,tMis,indDebit = computeMisfit(sigt,pInv.model,pInv.misfit,PF,Dobs,Wd,false,indCredit)
			end
			His.timeMisfit[iter+1,:]+=tMis'
			
			# compute regularizer(m0,mref,pInv.model.MInv,pInv.model.Iact);
			tic()
			R,dR,d2R = pInv.regularizer(mt,pInv.mref,pInv.MInv,Iact=pInv.Iact,C=pInv.regparam)
			His.timeReg[iter+1] += toq()
			# objective function
			Jt  = F  + alpha*R
			println(@sprintf( "   .%d\t%3.2e\t%3.2e\t\t\t%3.2e",
			           lsIter, F,       R,       Jt/J0))
			
			if Jt < Jc
			    break
			end
			muLS /=2; lsIter += 1
			if lsIter > 6
			    warn("Line search failed!")
			    return mc,Dc,outerFlag
			end
		end
		## End Line search
		
		## Check for termination 
		stepnorm = norm(mt-mc,Inf)
		mc = mt
		Jc = Jt
		
		sig, dsig = pInv.modelfun(mc)
		
		Active = (mc .<=low) | (mc.>=high)  # Compute active set
		  
		#  Check stopping criteria for outer iteration. 
		updateHis!(iter,His,Jc,F,Dc,R,alpha,countnz(Active),stepnorm,tMis,tReg)
	
		if iter >= maxIter
			outerFlag = 1
			println(" \n\n\n  *** Max iterations exceeded ***\n\n\n\n")
			break
		elseif stepnorm < stepTol
			outerFlag = 2
			println("  \n\n\n    *** Step norm tolerance met ***\n")
			break
		end
		# Evaluate gradient
		tic()
		if isempty(indCredit)
			dF = computeGradMisfit(sig,pInv.model,Dc,Dobs,Wd,pInv.misfit,PF)
		else
			dF = computeGradMisfit(sig,pInv.model,Dc,Dobs,Wd,pInv.misfit,PF,indDebit)
		end
		His.timeGradMisfit[iter+1]+=toq()
		
		dF = dsig'*dF
		gc = dF + alpha*dR
		
		if !isempty(filenamePrefix)
			filename = string(filenamePrefix,"_GN",iter,".dat");
			writedlm(filename,reshape(pInv.modelfun(mc)[1],tuple((pInv.MInv.n+1)...)));
		end
		
		
		
	end # while outer_flag == 0
	
	return mc,Dc,outerFlag,His
end  # Optimization code


