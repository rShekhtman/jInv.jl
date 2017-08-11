export projGN, projGNCG, computeMisfitOutput


"""
	function projGrad

	Projects gradient, i.e.,

					 | gc[i],  				xl[i] < x[i] < xh[i]
	gc[i]  = | max(gc[i],0)	 	xc[i] == xh[i]
					 | min(gc[i],0)   xc[i] == xl[i]


	Input:

		gc 				  - gradient vector
		mc 				  - model
		boundsLow   - lower bounds
		boundsHigh  - upper bounds
"""
function projGrad(gc,mc,boundsLow,boundsHigh)
	pgc = gc[:]
	pgc[mc.==boundsLow] = min(gc[mc.==boundsLow],0)
	pgc[mc.==boundsHigh] = max(gc[mc.==boundsHigh],0)
	return pgc
end


function dummy(mc,Dc,iter,pInv,pMis)
# this function does nothing and is used as a default for dumpResults().
end;


type computeMisfitOutput
   # Used to store the output of computeMisfit() for later.
   Dc
   F
   d2F
   pMis
   tMis
end  # type computeMisfitOutput

	
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
		solveGN::Function   - solver for GN system, default: projPCG, other options: projGNexplicit
		                      The interface for this method is: 
							            dm,his = solveGN(gc,pMis,pInv,sig,dsig,d2F,d2R,Active)
							
	Output:
		mc                  - final model
		Dc                  - data
		outerFlag           - flag for convergence.
		                        = -1 maximum number of iterations.
		                        = -2 could not reduce misfit during line search.
		                        =  1 last step was less than stepTol.
		His                 - iteration history
	
"""
function projGN(mc, pInv::InverseParam, pMis;
                indCredit=[],
                dumpResults::Function = dummy,
                out::Int=2,
                solveGN::Function=projPCG,
                outputcm::Bool=false,  # are we using cm?
                cm::computeMisfitOutput=computeMisfitOutput([],[],[],[],[])  )

   logfile = "jInv.out"
	maxIter     = pInv.maxIter      #  Max. no. iterations.
	pcgMaxIter  = pInv.pcgMaxIter   #  Max cg iters.
	pcgTol      = pInv.pcgTol       #  CG stopping tolerance.
	stepTol     = pInv.minUpdate    #  Step norm stopping tol.
	maxStep     = pInv.maxStep
	low         = pInv.boundsLow
	high        = pInv.boundsHigh
	alpha       = pInv.alpha        # tradeoff parameter
	
	His = getGNhis(maxIter,pcgMaxIter)
	#---------------------------------------------------------------------------
	#  Initialization.
	#---------------------------------------------------------------------------
	
	Active = (mc .<=low) | (mc.>=high)  # Compute active set
	
	
	## evaluate function and derivatives
   sig,dsig = pInv.modelfun(mc)
   if isempty(cm.Dc) || !outputcm
      # Initial misfit
      if isempty(indCredit)
         Dc,F,dF,d2F,pMis,tMis = computeMisfit(sig,pMis, true)
      else
         Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sig,pMis, true,indCredit)
      end
      
   elseif outputcm
      # Use results from previous iteration.
      Dc,F,d2F,pMis,tMis = cm.Dc, cm.F, cm.d2F, cm.pMis, cm.tMis
      
		if isempty(indCredit)
			dF = computeGradMisfit(sig,Dc,pMis)
		else
			dF = computeGradMisfit(sig,Dc,pMis,indCredit)
		end
      
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
	iter = 0   # iteration counter
	MaxlsIter = 6  # maximum number of line search iterations
	outerFlag = -1; stepNorm=0.0

	outStr = @sprintf("\n %4s\t%08s\t%08s\t%08s\t%08s\t%08s\n", 
					  	"i.LS", "F", "R","alpha[1]","Jc/J0","#Active")
	updateHis!(0,His,Jc,norm(projGrad(gc,mc,low,high)),F,Dc,R,alpha[1],countnz(Active),0.0,-1,tMis,tReg)
	
	if out>=2; print(outStr); end
	f = open(logfile, "a")
	write(f, outStr)
	
	while outerFlag == -1
		
		iter += 1
		outStr = @sprintf("%3d.0\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3d\n", 
		         iter, F, R,alpha[1],Jc/J0,countnz(Active))
		if out>=2; print(outStr); end
		write(f, outStr)

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
			mt[mt.<low]  = low[ mt.<low]
			mt[mt.>high] = high[mt.>high]
			## evaluate function 
			sigt, = pInv.modelfun(mt)
			if isempty(indCredit)
				Dc,F,dF,d2F,pMis,tMis = computeMisfit(sigt,pMis, false)
			else
				Dc,F,dF,d2F,pMis,tMis,indDebit = computeMisfit(sigt, false,indCredit)
			end
			His.timeMisfit[iter+1,:]+=tMis
			
			tic()
			R,dR,d2R = computeRegularizer(pInv.regularizer,mt,pInv.mref,pInv.MInv,alpha) 
			His.timeReg[iter+1] += toq()
			# objective function
			Jt  = F  + R
			
			outStr = @sprintf( "   .%d\t%3.2e\t%3.2e\t\t\t%3.2e\n",
			                   lsIter, F,       R,       Jt/J0)
			if out>=2;
				print(outStr)
			end
   		write(f, outStr)
   		flush(f)
			
			if Jt < Jc
			    break
			end
			muLS /=2; lsIter += 1
			if lsIter > MaxlsIter
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
		
	end # while outer_flag == -1


	if out>=1
		if outerFlag==-1
         outStr = @sprintf("projGN iterated maxIter=%i times but reached only stepNorm of %.2e instead %.2e.",
                            maxIter, stepNorm, stepTol)
		elseif outerFlag==-2
			outStr = "projGN stopped at iteration $iter with a line search fail."
		elseif outerFlag==1
			outStr = "projGN reached desired accuracy at iteration $iter."
		end
		println(outStr)
		println(f, outStr)
	end  # out>=1

	close(f)

   if outputcm
      cm = computeMisfitOutput( Dc,F,d2F,pMis,tMis )  # save for next tradeoff parameter
         
      return mc,Dc,outerFlag,His, cm

   else
      return mc,Dc,outerFlag,His
   end
end  # function projGN  Optimization code


projGNCG = projGN
