using jInv.Mesh
using jInv.LinearSolvers
using Base.Test
using KrylovMethods

println("===  Example 2D DivSigGrad ====");

domain = [0.0, 1.0, 0.0, 1.0];
n      = [50,50];
Mr     = getRegularMesh(domain,n)
G      = getNodalGradientMatrix(Mr);
m      = spdiagm(exp(randn(size(G,1))));
Ar     = G'*m*G;
Ar     = Ar + 1e-1*norm(Ar,1)*speye(size(Ar,2));
N      = size(Ar,2); 
b      = Ar*rand(N);

sPCG   = getIterativeSolver(KrylovMethods.cg,out=1,sym=1);
x,  = solveLinearSystem(Ar,b,sPCG);
@test norm(Ar*x-b)/norm(b) < sPCG.tol

IterMethod = (A,b;M=M,tol=1e-1,maxIter=10,out=-1) ->
                  bicgstb(A,b,M1=M,tol=tol,maxIter=maxIter,out=out)
sBiCG = getIterativeSolver(IterMethod,out=1);
x, = solveLinearSystem(Ar,b,sBiCG);
@test norm(Ar*x-b)/norm(b) < sBiCG.tol

println("===  Test for nonsymmetric matrices ====");
n = 100
m = 4
sBiCG.PC = :jac
sBiCG.doClear=true
sBiCG.out=0
A = sprandn(n,n,5/n) + 10*speye(n)
B = randn(n,m)
X, = solveLinearSystem(A,B,sBiCG)
for k=1:m
	@test norm(A*X[:,k] - B[:,k])/norm(B[:,k]) < sBiCG.tol
end

sBiCG.doClear=true
sBiCG.isTranspose=true
Xt, = solveLinearSystem(A',B,sBiCG)
for k=1:m
	#@test norm(A*X[:,k] - B[:,k])/norm(B[:,k]) < sBiCG.tol
	@test norm(Xt[:,k] - X[:,k])/norm(X[:,k]) < 1e-13
end
