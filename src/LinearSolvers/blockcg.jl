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
	X,flag,err,iter = blockcg!(Af,full(B),X=X, tol=param.tol,maxIter=param.maxIter,M=param.Ainv,out=param.out,ortho=param.ortho)
	param.nIter+=iter*nrhs
	param.timeCG+=toq();
	return X, param
end # function solveLinearSystem PCGsolver



"""
blockcg(A,B,X)

Preconditioned Conjugate Gradient Method for solving

A * X = B

Input:

Output:


Reference: 

"""
function blockcg!(A::Function,B::Array;X=zeros(eltype(B),size(B)),M::Function=identity,maxIter=20,tol=1e-2,ortho::Bool=false,pinvTol =eps(Float64)*size(B,1),out::Int=0,storeInterm::Bool=false)
# [X] = blkPCG(A,B,M,iter,tol)
# 
R = copy(B)
Z = M(R)
if ortho
    P, = gschmidt(Z)
else
    P = Z
end
 
n, nrhs = size(X)
nB      = computeNorm(B)
resmat  = zeros(maxIter,nrhs) 

if storeInterm
    Xout = zeros(n,nrhs,maxIter)	# allocate space for intermediates
end

if out==2
	println("=== blockcg ===")
	println(@sprintf("%4s\t%7s","iter","max(relres)"))
end

# pre-allocate
PTQ = zeros(nrhs,nrhs)
PTR = zeros(nrhs,nrhs)
QTZ = zeros(nrhs,nrhs)

flag = -1;
iter = 0
for iter=1:maxIter
	Q     = A(P)
    # PTQ   = P'*Q
	BLAS.gemm!('T','N',1.0,P,Q,0.0,PTQ)

    # Alpha = (PTQ)\(P'*R);
    PTR = P'*R
    BLAS.gemm!('T','N',1.0,P,R,0.0,PTR)
	pinvPTQ = getPinv!(PTQ,pinvTol)
    Alpha = pinvPTQ*PTR
    
    # X     += P*Alpha
    BLAS.gemm!( 'N','N',1.0, P, Alpha, 1.0, X)
	if storeInterm; Xout[:,:,iter] = x; end
    # R     -= Q*Alpha
    BLAS.gemm!('N','N',-1.0,Q,Alpha,1.0,R)

    resmat[iter,:] = computeNorm(R)./nB
    if out==2
		println(@sprintf("%3d\t%1.2e",iter,maximum(resmat[iter,:])))
	end
    if maximum(resmat[iter,:]) < tol
	    flag = 0;
        break;
    end
    
    Z     = M(R)
    #Beta  = -(PTQ)\(Q'*Z);
	BLAS.gemm!('T','N',1.0,Q,Z,0.0,QTZ)
    Beta  = -pinvPTQ*QTZ
    
    if ortho     
        P,     =  gschmidt(Z + P*Beta)
    else
	    P     = Z  + P*Beta;
	end
end
if out>=0
	if flag==-1
		println(@sprintf("blockg iterated maxIter (=%d) times but reached only residual norm %1.2e instead of tol=%1.2e.",
																							maxIter,maximum(resmat[iter,:]),tol))
	elseif flag==0 && out>=1
		println(@sprintf("blockcg achieved desired tolerance at iteration %d. Residual norm is %1.2e.",iter,maximum(resmat[iter,:])))
	end
end
if storeInterm
    return Xout[:,:,1:iter],flag,resmat[iter,:],iter,resmat[1:iter,:]
else
    return X,flag,resmat[iter,:],iter,resmat[1:iter,:]
end
end

function computeNorm(R)
	n,nrhs = size(R)
	res    = zeros(nrhs)
	for k=1:nrhs
		for i=1:n
			res[k]+=R[i,k]*R[i,k]
		end
	end
	return sqrt(res)
end

function getPinv!(A,pinvTol)
	SVD = svdfact!(A)
	Sinv        = zeros(length(SVD.S))
    index       = SVD.S .> pinvTol*maximum(SVD.S)
    Sinv[index] = 1.0./ SVD.S[index]
    Sinv[find(!isfinite(Sinv))] = 0.0
    return SVD.Vt'scale(Sinv, SVD.U')
end
 
function gschmidt(V)
# Input: V is an m by n matrix of full rank m<=n
# Output: an m-by-n upper triangular matrix R
# and an m-by-m unitary matrix Q so that A = Q*R.
 
m,n      = size(V)
R      = zeros(eltype(V),n,n)
Q      = zeros(eltype(V),m,n)
R[1,1] = norm(V[:,1])
Q[:,1] = V[:,1]./R[1,1]
for k=2:n
    R[1:k-1,k] = Q[:,1:k-1]'*V[:,k]
    Q[:,k]     = V[:,k]-Q[:,1:k-1]*R[1:k-1,k]
    R[k,k]     = norm(Q[:,k])
    Q[:,k]     = Q[:,k]./R[k,k]
end
return Q,R
end