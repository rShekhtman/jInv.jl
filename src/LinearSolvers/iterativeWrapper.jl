export IterativeSolver,getIterativeSolver,solveLinearSystem

""" 
type jInv.LinearSolvers.IterativeSolver <: AbstractSolver
	
Fields:

	PC      - symbol, (:ssor, :jac,...)
	maxIter - maximum number of iterations
	tol     - tolerance
	Ainv    - preconditioner
	out     - flag for output
	doClear - flag for deleting preconditioner

Example:

	getIterativeSolver(cg)

"""
type IterativeSolver<: AbstractSolver
	IterMethod::Function
	PC::Symbol
	maxIter::Int
	tol::Real
	Ainv
	out::Int
	doClear::Bool
	nthreads::Int
	nIter::Int
	nBuildPC::Int
	timePC::Real
	timeSolve::Real
	timeMV::Real
end

"""
function jInv.LinearSolvers.getIterativeSolver
	
constructs IterativeSolver

Required Input:

	IterMethod::Function   - function handle for linear solvers 
		Inputs are: (A,b,M), A is matrix, b is right hand side, M is preconditioner
		Outputs are: (x,flag,err,iter), x is approximate solution

Optional Inputs:

	PC::Symbol     - specifies preconditioner, default:ssor
	maxIter        - maximum number of iterations, default:500
	tol            - tolerance on relative residual, default=1e-5
	Ainv           - preconditioner, default=identity
	out            - flag for output, default=-1 (no output)
	doClear        - flag for clearing the preconditioner, default=true
	nthreads       - number of threads to use for matvecs (requires ParSpMatVec.jl), default=4
           
"""
function getIterativeSolver(IterMethod::Function;PC=:ssor,maxIter=500,tol=1e-5,Ainv=identity,out=-1,doClear::Bool=true,nthreads::Int=4)
 	return IterativeSolver(IterMethod,PC,maxIter,tol,Ainv,out,doClear,nthreads,0,0,.0,.0,.0)
 end




function solveLinearSystem!(A,B,X,param::IterativeSolver,doTranspose=0)
	if param.doClear
		# clear preconditioner
		clear!(param)
		param.doClear=false
	end
	
	# build preconditioner
	if param.Ainv == []
		if param.PC==:ssor
			OmInvD = 1./diag(A);
			x      = zeros(eltype(A),size(B,1))
			M(r)   = (x[:]=0.0; tic(); x=ssorPrecTrans!(A,x,r,OmInvD); param.timePC+=toq(); return x);
			param.Ainv= M
		elseif param.PC==:jac
			OmInvD = 1./diag(A)
			M(r)   = (tic(); x=r.*OmInvD; param.timePC+=toq(); return x); 
			param.Ainv= M
		else 
			error("Iterativesolver: preconditioner $(param.PC) not implemented.")
		end
		param.nBuildPC+=1
	end
	
	# solve systems
	y     = zeros(eltype(A),size(X,1))
	if hasParSpMatVec
		Af(x) = (y[:]=0.0; tic(); ParSpMatVec.Ac_mul_B!(one(eltype(A)),A,x,zero(eltype(A)),y,param.nthreads); param.timeMV+=toq(); return y)
	else
		Af(x) = (y[:]=0.0; tic(); Ac_mul_B!(one(eltype(A)),A,x,zero(eltype(A)),y); param.timeMV+=toq(); return y)
	end
	
	tic()
	for i=1:size(X,2)
		bi      = vec(full(B[:,i]))
		# if param.IterMethod ==:cg
			X[:,i],flag,err,iter = param.IterMethod(Af,bi,param.Ainv,tol=param.tol,maxIter=param.maxIter,out=param.out)
		# elseif param.IterMethod ==:bicgstb
			# X[:,i],flag,err,iter = bicgstb(Af,bi, tol=param.tol,maxIter=param.maxIter,M1=param.Ainv,out=param.out)
		# end
		param.nIter+=iter
	end	
	param.timeSolve+=toq();
	return X, param
end # function solveLinearSystem PCGsolver
