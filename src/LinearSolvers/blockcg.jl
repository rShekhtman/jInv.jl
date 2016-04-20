export BlockPCGsolver,getBlockPCGsolver,solveLinearSystem

""" 
type BlockPCGsolver
	
Fields:
PC      - symbol, (:ssor, :jac,...)
maxIter - maximum number of iterations
tol     - tolerance
Ainv    - preconditioner
out     - flag for output
doClear - flag for deleting preconditioner

Example
getBlockPCGsolver()
"""
type BlockPCGsolver<: AbstractSolver
	PC::Symbol
	maxIter::Int
	tol::Real
	Ainv
	out::Int
	doClear::Bool
	ortho::Bool
	nthreads::Int
	nIter::Int
	nBuildPC::Int
	timePC::Real
	timeCG::Real
	timeMV::Real
end

function getBlockPCGsolver(PC=:ssor,maxIter=10,tol=1e-4,Ainv=identity,out=-1,doClear::Bool=true,ortho::Bool=false,nthreads::Int=4)
	return BlockPCGsolver(PC,maxIter,tol,Ainv,out,doClear,ortho,nthreads,0,0,.0,.0,.0)
end


function solveLinearSystem!(A,B,X,param::BlockPCGsolver,doTranspose=0)
	if param.doClear
		# clear preconditioner
		clear!(param)
		param.doClear=false
	end
	
	n = size(B,1)
	nrhs = size(B,2)
	# build preconditioner
	if param.Ainv == []
		if param.PC==:ssor
			OmInvD = 1./diag(A);
			Xt      = zeros(n,nrhs)
			M(R)   = (Xt[:]=0.0; tic(); Xt=ssorPrecTrans!(A,Xt,R,OmInvD); param.timePC+=toq(); return Xt);
			param.Ainv= M
		elseif param.PC==:jac
			OmInvD = 1./diag(A)
			M(r)   = (tic(); for k=1:size(r,2); x[:,k]=r[:,k].*OmInvD; end; param.timePC+=toq(); return x); 
			param.Ainv= M
		else 
			error("PCGsolver: preconditioner $(param.PC) not implemented.")
		end
		param.nBuildPC+=1
	end
	
	# solve systems
	Y    = zeros(n,nrhs)
	if hasParSpMatVec
		Af(X) = (Y[:]=0.0; tic(); ParSpMatVec.Ac_mul_B!(1.0,A,X,0.0,Y,param.nthreads); param.timeMV+=toq(); return Y)
	else
		Af(X) = (Y[:]=0.0; tic(); Ac_mul_B!(1.0,A,X,0.0,Y); param.timeMV+=toq(); return Y)
	end
		
	tic()
	X[:]=0.0
	X,flag,err,iter = KrylovMethods.blockcg(Af,full(B),X=X, tol=param.tol,maxIter=param.maxIter,M=param.Ainv,out=param.out,ortho=param.ortho)
	param.nIter+=iter*nrhs
	param.timeCG+=toq();
	return X, param
end # function solveLinearSystem PCGsolver



