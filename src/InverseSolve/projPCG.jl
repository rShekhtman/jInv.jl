export projPCG

"""
	dm = projPCG(H,g,Active,Precond,cgTol,maxIter)
	
	Projected Preconditioned Conjugate Gradient method for solving
	
		H*dm = g    subject to    dm[!Active] == 0 
	
	Input:
	
		H::Function       -  computes action of Hessian
		g::Vector         -  right hand side
		Active            -  describes active cells
		Precond::Function - preconditioner
		cgTol             - tolerance
		maxIter             - maximum number of iterations
"""
function projPCG(H::Function,gc::Vector,Active::BitArray,Precond::Function,cgTol::Real,maxIter::Int)
#  PCG over active cells
    his        = zeros(maxIter,3)

	delm       = zeros(size(gc))
	cgiter     = 0
	resid      = - !Active.*gc
	normResid0 = norm(resid)
	while true
		cgiter += 1
		tic()
		dc = !Active.*Precond(resid)
		his[cgiter,2]=toq()
		rd = dot(resid,dc)
	
		#  Compute conjugate direction pc.
		if cgiter == 1
			pc = dc
		else
			betak = rd / rdlast
			pc = dc + betak * pc
		end
	
		#  Form product Hessian*pc.
		tic()
		Hp = H(pc)
		his[cgiter,3] = toq()
		
		Hp = !Active.*Hp
		#  Update delm and residual.
		alphak = rd / dot(pc,Hp)
		delm   = delm + alphak*pc
		resid  = resid - alphak*Hp
		rdlast = rd
	    
	    his[cgiter,1] = norm(resid)/normResid0;
		if (his[cgiter,1] <= cgTol) || (cgiter == maxIter )
	  	  break
		end
	end

	return delm,his[1:cgiter,:]
end