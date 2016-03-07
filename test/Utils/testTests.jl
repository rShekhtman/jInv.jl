print("   test checkDerivative ...")
using jInv.Utils
using Base.Test

function quadRight(x,v=[])
	f = x.^2
	if !(isempty(v))
		return f, diagm(2*x)*v
	else
		return f
	end
end
x0 = randn(10).+2
p,e,o = checkDerivative(quadRight,x0,out=true)
@test p==true

function quadWrong(x,v=[])
    f = x.^2
    if !(isempty(v))
		return f, -diagm(2*x)*v
	else
		return f
	end
end
x0 = randn(10)
p,e,o = checkDerivative(quadWrong,x0,out=false)
@test p==false

print("passed\n")
