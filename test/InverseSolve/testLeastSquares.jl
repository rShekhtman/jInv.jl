@everywhere begin
	using jInv.InverseSolve
	using jInv.Mesh
	using jInv.LinearSolvers
	using jInv.ForwardShare
	using Base.Test
end

type LSparam <: ForwardProbType
	A::SparseMatrixCSC
	Ainv
end

function jInv.ForwardShare.getData(m::Vector,pFor::LSparam)
	return pFor.A*m,pFor
end

function jInv.ForwardShare.getSensMatVec(v::Vector,m::Vector,pFor::LSparam)
	return pFor.A*v
end

function jInv.ForwardShare.getSensTMatVec(v::Vector,m::Vector,pFor::LSparam)
	return pFor.A'*v
end

# build domain and true image
domain = [0.0 1.0 0.0 1.0]
n      = [16,16]
Minv   = getRegularMesh(domain,n)
xc     = getCellCenteredGrid(Minv)
xtrue = sin(2*pi*xc[:,1]).*sin(pi*xc[:,2])

# get noisy data 
A     = speye(Minv.nc)
ids   = sort(rand(1:Minv.nc,round(Int64,Minv.nc*.8)))
A     = A[ids,:]
btrue = A*xtrue
bdata = btrue + .1*randn(size(btrue))/norm(btrue,Inf)
Wd    = ones(length(btrue))

# generate misfit param
i1     = (1:round(Int64,size(A,1)/2))
i2     = (round(Int64,size(A,1)/2)+1:size(A,1))
pFor1  = LSparam(A[i1,:],[])
pFor2  = LSparam(A[i2,:],[])
sigmaBack = zeros(length(xtrue))
Iact         = speye(Minv.nc)
gl1          = getGlobalToLocal(Iact,sigmaBack)
gl2          = getGlobalToLocal(Iact,sigmaBack)
bd1          = bdata[i1]
bd2          = bdata[i2]
Wd1          = Wd[i1]
Wd2          = Wd[i2]

pMisRefs    = Array{RemoteRef{Channel{Any}}}(2)
pMisRefs[1] = @spawn getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1) 
pMisRefs[2] = @spawn getMisfitParam(pFor2,Wd2,bd2,SSDFun,fMod,gl2) 

pMis    = Array{MisfitParam}(2)
pMis[1] = getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1) 
pMis[2] = getMisfitParam(pFor2,Wd2,bd2,SSDFun,fMod,gl2) 

# setup pInv
alpha        = 1.
x0           = mean(xtrue)*ones(Minv.nc)
boundsLow    = minimum(xtrue)*ones(Minv.nc)
boundsHigh   = maximum(xtrue)*ones(Minv.nc)
sigmaBack    = zeros(Minv.nc)

#  solve with automatic distribution
pInv         = getInverseParam(Minv,Iact,fMod,diffusionReg,alpha,x0,boundsLow,boundsHigh)
pInv.maxIter = 5
x1, = projGNCG(x0,pInv,pMis)
pInv.maxIter = 5
x2, = projGNCG(x0,pInv,pMisRefs)
@test norm(x1-x2)/norm(x1) < 1e-12