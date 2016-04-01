using jInv.InverseSolve
using Base.Test
using jInv.Mesh
using jInv.Utils

# build regular mesh and Iact
domain = [0;5;0;5; 0;3]
n     = [12;13;10]
M     = getRegularMesh(domain,n)
idx   = ones(Int,(n[1],n[2],n[3]))
idx[:,:,1:5] = 0
Iact  = speye(Bool,M.nc)
Iact  = Iact[:,vec(idx).==1] 
mc    = randn(size(Iact,2))

regFuns = [ 
			(m,mref,M)->diffusionReg(m,mref,M,Iact=Iact), 
			(m,mref,M)->wdiffusionReg(m,mref,M,Iact=Iact), 
			smallnessReg, 
			(m,mref,M)->wTVReg(m,mref,M,Iact=Iact)]
for k=1:length(regFuns)
	println("checkDerivative of $(regFuns[k])")
	
	function testFun(x,v=[])
		Sc,dS,d2S = regFuns[k](x,0.*x,M)
		if isempty(v)
			return Sc
		else
			return Sc,dot(dS,v)
		end
	end
	chkDer, = checkDerivative(testFun,mc,out=false)
	@test chkDer
end

regFuns = [wdiffusionRegNodal]
idx   = ones(Int,(n[1]+1,n[2]+1,n[3]+1))
idx[:,:,1:5] = 0
Iact  = speye(Bool,prod(M.n+1))
Iact  = Iact[:,vec(idx).==1] 
mc    = randn(size(Iact,2))

for k=1:length(regFuns)
	println("checkDerivative of $(regFuns[k])")
	
	function testFun(x,v=[])
		Sc,dS,d2S = regFuns[k](x,0.*x,M,Iact=Iact)
		if isempty(v)
			return Sc
		else
			return Sc,dot(dS,v)
		end
	end
	chkDer, = checkDerivative(testFun,mc,out=false)
	@test chkDer
end
