using jInv.InverseSolve
using Base.Test
using jInv.Mesh
using jInv.Utils

# setup forward mesh
domain = [0;5;0;5; 0;3]
n     = [12;13;10]
Minv  = getRegularMesh(domain,n)

modFun = (expMod,boundMod,idMod)
for k=1:length(modFun)
	println("\t\tcheckDerivative of $(modFun[k])")
	function testModFun(m,v=[])
		mc,dm = modFun[k](m)
		if !isempty(v)
			dm = dm*v
			return mc,dm
		end
		return mc
	end
	derivativeDistModel, = checkDerivative(testModFun,rand(prod(n)),out=false)
	@test derivativeDistModel
end


