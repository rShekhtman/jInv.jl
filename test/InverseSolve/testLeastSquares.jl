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

# generate forward param
i1     = (1:round(Int64,size(A,1)/2))
i2     = (round(Int64,size(A,1)/2)+1:size(A,1))
pFor1  = @spawn LSparam(A[i1,:],[])
pFor2  = @spawn LSparam(A[i2,:],[])
pFor   = [fetch(pFor1);fetch(pFor2)]
pForp  = [pFor1;pFor2]

# setup pInv
alpha        = 1.
x0           = mean(xtrue)*ones(Minv.nc)
Iact         = speye(Minv.nc)
boundsLow    = minimum(xtrue)*ones(Minv.nc)
boundsHigh   = maximum(xtrue)*ones(Minv.nc)
sigmaBack    = zeros(Minv.nc)
gl1          = @spawnat pForp[1].where getGlobalToLocal(Iact,sigmaBack)
gl2          = @spawnat pForp[2].where getGlobalToLocal(Iact,sigmaBack)
bd1          = @spawnat pForp[1].where identity(bdata[i1])
bd2          = @spawnat pForp[2].where identity(bdata[i2])
Wd1          = @spawnat pForp[1].where identity(Wd[i1])
Wd2          = @spawnat pForp[2].where identity(Wd[i2])

gloc         = [fetch(gl1);fetch(gl2)]
glocp        = [gl1; gl2]
dobs         = Array{Any}(2); 
dobs[1] = fetch(bd1); 
dobs[2] = fetch(bd2)
dobp         = [bd1;bd2]
Wds          = Array{Any}(2); 
Wds[1] = fetch(Wd1); 
Wds[2] = fetch(Wd2)
Wdp          = [Wd1;Wd2]

#  solve with automatic distribution
pInv         = getInverseParam(gloc,Minv,Iact,fMod,diffusionReg,alpha,x0,
                                    [],SSDFun,dobs,Wds,[],boundsLow,boundsHigh)
pInv.maxIter = 5
x1, = projGNCG(x0,pInv,pFor)

# solve with pre-distribution
pInv         = getInverseParam(glocp,Minv,Iact,fMod,diffusionReg,alpha,x0,
                                    [],SSDFun,dobp,Wdp,[],boundsLow,boundsHigh)
pInv.maxIter = 5
x2, = projGNCG(x0,pInv,pForp)

@test norm(x1-x2)/norm(x1) < 1e-12