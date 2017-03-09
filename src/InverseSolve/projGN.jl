export projGN
export projGNhis
export normalEqGN


type projGNhis
	Jc::Array
	dJ::Array
	F::Array
	Dc::Array
	Rc::Array
	alphas::Array
	Active::Array
	stepNorm::Array
	lsIter::Array
	timeMisfit::Array
	timeReg::Array
	timePCG::Array
	hisLinSol::Array
	timeLinSol::Array
	timeGradMisfit::Array
end

function getProjGNhis(maxIter,maxIterCG)
	Jc = zeros(maxIter+1)
	dJ = zeros(maxIter+1)
	F  = zeros(maxIter+1)
	Dc = []
	Rc = zeros(maxIter+1)
	alphas = zeros(maxIter+1)
	Active = zeros(maxIter+1)
	stepNorm = zeros(maxIter+1)
	lsIter = zeros(Int,maxIter+1)
	timeMisfit = zeros(maxIter+1,4)
	timeReg = zeros(maxIter+1)
	timePCG = zeros(maxIter+1,1)
	hisLinSol = []
	timeLinSol = zeros(maxIter+1,1)
	timeGradMisfit = zeros(maxIter+1,2)

	return projGNhis(Jc,dJ,F,Dc,Rc,alphas,Active,stepNorm,lsIter,timeMisfit,timeReg,timePCG,hisLinSol,timeLinSol,timeGradMisfit)
end

function updateHis!(iter::Int64,His::projGNhis,Jc::Real,dJ::Real,Fc,Dc,Rc::Real,alpha::Real,nActive::Int64,stepNorm::Real,lsIter::Int,timeMisfit::Vector,timeReg::Real)
	His.Jc[iter+1]            = Jc
	His.dJ[iter+1]            = dJ
	His.F[iter+1]             = Fc
	push!(His.Dc,Dc)
	His.Rc[iter+1]            = Rc
	His.alphas[iter+1]        = alpha
	His.Active[iter+1]        = nActive
	His.stepNorm[iter+1]      = stepNorm
	His.lsIter[iter+1]        = lsIter
	His.timeMisfit[iter+1,:] += timeMisfit
	His.timeReg[iter+1]      += timeReg[]
end



import jInv.InverseSolve.projPCG
function projPCG(gc,pMis,pInv,sig,dsig,d2F,d2R,Active)
		
		#  Set up Hessian and preconditioner.
		Hs(x) = dsig'*HessMatVec(dsig*x,pMis,sig,d2F) + d2R*x; 
		#  build preconditioner
		pInv.HesPrec.setupPrec(Hs, d2R,pInv.HesPrec.param);
		PC(x) = pInv.HesPrec.applyPrec(Hs,d2R,x,pInv.HesPrec.param);		
		# call projPCG and return
		return projPCG(Hs,gc,Active,PC,pInv.pcgTol,pInv.pcgMaxIter)
end

function normalEqGN(gc,pMis,pInv,sig,dsig,d2F,d2R,Active)
	if all(Active)
		return 0*gc,[0.0;0.0]
	end
	# get Hessian of misfit
	tic();
	Hm = getHessian(sig,pMis,d2F)
	timeBuild = toq();
	
	# build overall Hessian
	H  = dsig'*Hm*dsig + d2R
	
	# remove Active constraints
	tic()
	Hr = H[!Active,!Active]
	gr = gc[!Active]
	dr = -(Hr\gr)
	timeSolve=toq();
	dm = 0*gc
	dm[!Active] = dr
	
	# solve and return
	return dm,[timeBuild;timeSolve]
end
	
"""
	mc,Dc,outerFlag = projGN(mc,pInv::InverseParam,pMis, indFor = [], dumpResults::Function = dummy)
	
	(Projected) Gauss-Newton method. 
	
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
function  projGN(mc,pInv::InverseParam,pMis;indCredit=[],
	dumpResults::Function = dummy,out::Int=2,solveGN=projPCG)

	maxIter     = pInv.maxIter      #  Max. no. iterations.
	pcgMaxIter  = pInv.pcgMaxIter   #  Max cg iters.
	pcgTol      = pInv.pcgTol       #  CG stopping tolerance.
	stepTol     = pInv.minUpdate    #  Step norm stopping tol.
	maxStep     = pInv.maxStep
	low         = pInv.boundsLow
	high        = pInv.boundsHigh
	alpha       = pInv.alpha
	
	His = getProjGNhis(maxIter,pcgMaxIter)
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
	R,dR,d2R = computeRegularizer(pInv.regularizer,mc,pInv.mref,pInv.MInv,alpha) 
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

	outStr = @sprintf("\n %4s\t%08s\t%08s\t%08s\t%08s\t%08s\n", 
					  	"i.LS", "F", "R","alpha[1]","Jc/J0","#Active")
	updateHis!(0,His,Jc,norm(projGrad(gc,mc,low,high)),F,Dc,R,alpha[1],countnz(Active),0.0,-1,tMis,tReg)
	
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

		# solve linear system to find search direction
		tic()
		delm,hisLinSol = solveGN(gc,pMis,pInv,sig,dsig,d2F,d2R,Active)
		His.timeLinSol[iter+1]+= toq()
		push!(His.hisLinSol,hisLinSol)
		
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
			His.timeMisfit[iter+1,:]+=tMis
			
			tic()
			R,dR,d2R = computeRegularizer(pInv.regularizer,mt,pInv.mref,pInv.MInv,alpha) 
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
		updateHis!(iter,His,Jc,-1.0,F,Dc,R,alpha[1],countnz(Active),stepNorm,lsIter,tMis,tReg)
	
		dumpResults(mc,Dc,iter,pInv,pMis);
		if stepNorm < stepTol
			outerFlag = 1
			break
		elseif iter >= maxIter
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
		
		
		His.dJ[iter+1] = norm(projGrad(gc,mc,low,high))
		
	end # while outer_flag == 0
	
	if out>=1
		if outerFlag==-1
			println("projGN iterated maxIter=$maxIter times but reached only stepNorm of $(stepNorm) instead $(stepTol)." )
		elseif outerFlag==-2
			println("projGN stopped at iteration $iter with a line search fail.")
		elseif outerFlag==1
			println("projGN reached desired accuracy at iteration $iter.")
		end
	end
	
	
	return mc,Dc,outerFlag,His
end  # Optimization code


