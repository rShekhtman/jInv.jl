include("../setupTests.jl")

# build domain and true image
domain = [0.0 1.0 0.0 1.0]
n      = [32,32]
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
i1     = (1:round(Int64,size(A,1)/4));
i2     = (round(Int64,size(A,1)/4)+1:round(Int64,2*size(A,1)/4));
i3 	   = (round(Int64,2*size(A,1)/4)+1 : round(Int64,3*size(A,1)/4));
i4 	   = (round(Int64,3*size(A,1)/4)+1 : size(A,1));
pFor1  = LSparam(A[i1,:],[])
pFor2  = LSparam(A[i2,:],[])



sigmaBack = zeros(length(xtrue))
Iact         = speye(Minv.nc);
gl1          = getGlobalToLocal(1.0,sigmaBack)
gl2          = getGlobalToLocal(Iact,sigmaBack)

bd1          = bdata[i1]
bd2          = bdata[i2]
bd3          = bdata[i3]
bd4          = bdata[i4]
Wd1          = Wd[i1]
Wd2          = Wd[i2]
Wd3          = Wd[i3]
Wd4          = Wd[i4]

pMisRefs    = Array{RemoteRef{Channel{Any}}}(6)
pMisRefs[1] = @spawn getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1) 
pMisRefs[2] = @spawn getMisfitParam(pFor2,Wd2,bd2,SSDFun,fMod,gl2)

pForRefs = Array(RemoteRef{Channel{Any}},4)
pForRefs[1] = @spawn LSparam(A[i3,:],[])
pForRefs[2] = @spawn LSparam(A[i4,:],[])
pForRefs[3] = @spawn LSparam(A[i3,:],[])
pForRefs[4] = @spawn LSparam(A[i4,:],[])

Mesh2Mesh = Array(RemoteRef{Channel{Any}},2)	
for i=1:length(Mesh2Mesh)
	k = pForRefs[i].where
	if i==1
		Mesh2Mesh[i] = remotecall(k,speye,Minv.nc);
		wait(Mesh2Mesh[i]);
	elseif i==2
		Mesh2Mesh[i] = remotecall(k,identity,1.0);
		wait(Mesh2Mesh[i]);
	end
end

Wd = Array(Array{Float64},2);
dobs = Array(Array{Float64},2);
dobs[1] = bd3;
dobs[2] = bd4;
Wd[1] = Wd3;
Wd[2] = Wd4;
pMisRefs[3:4] = getMisfitParam(pForRefs[1:2],Wd,dobs,SSDFun,Iact,sigmaBack,Mesh2Mesh);
pMisRefs[5:6] = getMisfitParam(pForRefs[3:4],Wd,dobs,SSDFun,Iact,sigmaBack); # single mesh

M2M3		 = speye(Minv.nc);
M2M4		 = speye(Minv.nc);
pMis    = Array{MisfitParam}(6)
pMis[1] = getMisfitParam(pFor1,Wd1,bd1,SSDFun,fMod,gl1);
pMis[2] = getMisfitParam(pFor2,Wd2,bd2,SSDFun,fMod,gl2);
pMis[3] = fetch(pMisRefs[3]);
pMis[4] = fetch(pMisRefs[4]);
pMis[5] = fetch(pMisRefs[5]);
pMis[6] = fetch(pMisRefs[6]);

# setup pInv
alpha        = 1.
x0           = mean(xtrue)*ones(Minv.nc)
boundsLow    = minimum(xtrue)*ones(Minv.nc)
boundsHigh   = maximum(xtrue)*ones(Minv.nc)
sigmaBack    = zeros(Minv.nc)

#  solve with automatic distribution
pInv         = getInverseParam(Minv,fMod,diffusionReg,alpha,x0,boundsLow,boundsHigh)
pInv.maxIter = 5
x1, = projGNCG(x0,pInv,pMis)

pInv.maxIter = 5
x2, = projGNCG(x0,pInv,pMisRefs)
@test norm(x1-x2)/norm(x1) < 1e-12

clear!(pMisRefs);
for k=1:length(pMis)
	clear!(pMis[k],clearPFor=true, clearData=true,clearMesh2Mesh=true);
end


