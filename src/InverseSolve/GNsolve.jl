export normalEqGN

"""
function normalEqGN(gc,pMis,pInv,sig,dsig,d2F,d2R,Active)

explicitly builds normal equation, projects it and solves it. 

Inputs:

	gc        - gradient
	pMis      - misfit params
	pInv      - inverse param
	sig,dsig  - current model and derivative
	d2F,d2R   - Hessians os misfit and regularizer
	Active    - indicator of active set

Outputs:

	dm        - search direction
	times     - array containing time to build and solve Hessian
"""
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