using jInv.Mesh
using jInv.LinearSolvers
using Base.Test

println("===  Test Julia Wrapper: Symmetric ====");

domain = [0.0, 1.0, 0.0, 1.0];
n      = [50,50];
Mr     = getRegularMesh(domain,n)
G      = getNodalGradientMatrix(Mr);
m      = spdiagm(exp(randn(size(G,1))));
Ar     = G'*m*G;
Ar     = Ar + 1e-1*norm(Ar,1)*speye(size(Ar,2));
N      = size(Ar,2); 
b      = Ar*rand(N);

sSym   = getJuliaSolver(sym = 1);
x,  = solveLinearSystem(Ar,b,sSym);
@test norm(Ar*x-b)/norm(b) < 1e-10


println("===  Test Julia Wrapper: nonsymmetric matrices ====");
n = 100
sNonSym   = getJuliaSolver(sym = 0);
A = sprandn(n,n,5/n) + 10*speye(n)
B = randn(n)

X, = solveLinearSystem(A,B,sNonSym,0);
@test norm(A*X - B)/norm(B) < 1e-10

X, = solveLinearSystem(A,B,sNonSym,1);
@test norm(A'*X - B)/norm(B) < 1e-10

X, = solveLinearSystem(A,B,sNonSym,1);
@test norm(A'*X - B)/norm(B) < 1e-10

X, = solveLinearSystem(A,B,sNonSym,0);
@test norm(A*X - B)/norm(B) < 1e-10

println("===  End Test Julia Wrapper ====");