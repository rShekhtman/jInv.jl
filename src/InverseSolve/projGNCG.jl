export projGNCG

type projGNCGhis
	Jc::Array
	F::Array
	Dc::Array
	Rc::Array
	alphas::Array
	Active::Array
	stepNorm::Array
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
	stepNorm = zeros(maxIter+1)
	timeMisfit = zeros(maxIter+1,4)
	timeReg = zeros(maxIter+1)
	timePCG = zeros(maxIter+1,1)
	hisPCG = []
	timeGradMisfit = zeros(maxIter+1,2)

	return projGNCGhis(Jc,F,Dc,Rc,alphas,Active,stepNorm,timeMisfit,timeReg,timePCG,hisPCG,timeGradMisfit)
end

function updateHis!(iter::Int64,His::projGNCGhis,Jc::Real,Fc,Dc,Rc::Real,alpha::Real,nActive::Int64,stepNorm::Real,timeMisfit::Vector,timeReg::Real)
	His.Jc[iter+1]            = Jc
	His.F[iter+1]             = Fc
	push!(His.Dc,Dc)
	His.Rc[iter+1]            = Rc
	His.alphas[iter+1]        = alpha
	His.Active[iter+1]        = nActive
	His.stepNorm[iter+1]      = stepNorm
	His.timeMisfit[iter+1,:] += timeMisfit'
	His.timeReg[iter+1]      += timeReg
end

function dummy(mc,Dc,iter,pInv,pMis) 
# this function does nothing and is used as a default for dumpResults(). 
end;



"""
	mc,Dc,outerFlag = projGNCG(mc,pInv::InverseParam,pMis, indFor = [], dumpResults::Function = dummy)
	
	(Projected) Gauss-NewtonCG
	
	Input:
	
		mc::Vector          - intial guess for model
		pInv::InverseParam  - parameter for inversion
		pMis                - misfit terms
		indCredit           - indices of forward problems to work on
		dumpResults			- A function pointer for saving the results throughout the iterations.
							- We assume that dumpResults is dumpResults(mc,Dc,iter,pInv,pMis), 
							- where mc is the recovered model, Dc is the predicted data. 
							- If dumpResults is not given, nothing is done (dummy() is called).
		out::Int            - flag for output (-1: no output, 1: final status, 2: residual norm at each iteration)
							
	Output:
		mc                  - final model
		Dc                  - data
		outerFlag           - flag for convergence
		His                 - iteration history
	
"""
function  projGNCG(mc,pInv::InverseParam,pMis;indCredit=[],dumpResults::Function = dummy,out::Int=2)

	maxIter     = pInv.maxIter      #  Max. no. iterations.
	pcgMaxIter  = pInv.pcgMaxIter   #  Max cg iters.
	pcgTol      = pInv.pcgTol       #  CG stopping tolerance.
	stepTol     = pInv.minUpdate    #  Step norm stopping tol.
	maxStep     = pInv.maxStep
	low         = pInv.boundsLow
	high        = pInv.boundsHigh
	alpha       = pInv.alpha
	
	His = getProjGNCGhis(maxIter,pcgMaxIter)
	#---------------------------------------------------------------------------
	#  Initialization.
	#---------------------------------------------------------------------------
	
	Active = (mc .<=low) | (mc.>=high)  # Compute active set
	
	
	## evaluate function and derivatives
	sig,dsig = pInv.modelfun(mc)
	if isempty(indCredit)
		Dc,F,dF,d2F,pMis,tMis = computeMisfit(sig,pMis,true)
	else
		Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sig,pMis,true,indCredit)
	end
	dF = dsig'*dF
	
	
	# compute regularizer
	tic()
	R,dR,d2R = computeRegularizer(pInv.regularizer,mc,pInv.mref,pInv.MInv,pInv.Iact,alpha,pInv.regparam) 
	tReg = toq()    
	
	# objective function
	Jc  = F  + R
	gc  = dF + dR
	
	F0 = F; J0 = Jc
	############################################################################
	##  Outer iteration.                                                        #
	############################################################################
	iter = 0
	outerFlag = -1; stepNorm=0.0

	outStr = @sprintf("%4s\t%08s\t%08s\t%08s\t%08s\t%08s\n", 
					  	"i.LS", "F", "R","alpha[1]","Jc/J0","#Active")
	updateHis!(0,His,Jc,F,Dc,R,alpha[1],countnz(Active),0.0,tMis,tReg)
	
	if out>=2; print(outStr); end
	f = open("jInv.out", "w")
	write(f, outStr)
	close(f)
	
	while outerFlag == -1
		
		iter += 1
		outStr = @sprintf("%3d.0\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3d\n", 
		         iter, F, R,alpha[1],Jc/J0,countnz(Active))
		if out>=2; print(outStr); end
		f = open("jInv.out", "a")
		write(f, outStr)
		close(f)

		
		#  Set up Hessian and preconditioner.
		if isempty(indCredit)
			Hs(x) = dsig'*HessMatVec(dsig*x,pMis,sig,d2F) + d2R*x; 
		else
			Hs(x) = dsig'*HessMatVec(dsig*x,pMis,sig,d2F,indDebit) + d2R*x;
		end
		
		
		pInv.HesPrec.setupPrec(Hs, d2R,pInv.HesPrec.param);
		
		PC(x) = pInv.HesPrec.applyPrec(Hs,d2R,x,pInv.HesPrec.param);
		
		##  inner CG iterations.
		tic()
		delm,hisPCG = projPCG(Hs,gc,Active,PC,pcgTol,pcgMaxIter)
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
				Dc,F,dF,d2F,pMis,tMis = computeMisfit(sigt,pMis,false)
			else
				Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sigt,false,indCredit)
			end
			His.timeMisfit[iter+1,:]+=tMis'
			
			# compute regularizer(m0,mref,pInv.model.MInv,pInv.model.Iact);
			tic()
			R,dR,d2R = computeRegularizer(pInv.regularizer,mt,pInv.mref,pInv.MInv,pInv.Iact,alpha,pInv.regparam) 
			His.timeReg[iter+1] += toq()
			# objective function
			Jt  = F  + R
			if out>=2;
				println(@sprintf( "   .%d\t%3.2e\t%3.2e\t\t\t%3.2e",
			           lsIter, F,       R,       Jt/J0))
			end
			
			if Jt < Jc
			    break
			end
			muLS /=2; lsIter += 1
			if lsIter > 6
			    outerFlag = -2
				break
			end
		end
		## End Line search
		
		## Check for termination 
		stepNorm = norm(mt-mc,Inf)
		mc = mt
		Jc = Jt
		
		sig, dsig = pInv.modelfun(mc)
		
		Active = (mc .<=low) | (mc.>=high)  # Compute active set
		  
		#  Check stopping criteria for outer iteration. 
		updateHis!(iter,His,Jc,F,Dc,R,alpha[1],countnz(Active),stepNorm,tMis,tReg)
	
		if iter >= maxIter
			break
		elseif stepNorm < stepTol
			outerFlag = 1
			break
		end
		# Evaluate gradient
		tic()
		if isempty(indCredit)
			dF = computeGradMisfit(sig,Dc,pMis)
		else
			dF = computeGradMisfit(sig,Dcp,pMis,indDebit)
		end
		His.timeGradMisfit[iter+1]+=toq()
		
		dF = dsig'*dF
		gc = dF + dR
		
		dumpResults(mc,Dc,iter,pInv,pMis);
		
	end # while outer_flag == 0
	
	if out>=1
		if outerFlag==-1
			println("projGNCG iterated maxIter=$maxIter times but reached only stepNorm of $(stepNorm) instead $(stepTol)." )
		elseif outerFlag==-2
			println("projGNCG stopped at iteration $iter with a line search fail.")
		elseif outerFlag==1
			println("projGNCG reached desired accuracy at iteration $iter.")
		end
	end
	
	
	return mc,Dc,outerFlag,His
end  # Optimization code


