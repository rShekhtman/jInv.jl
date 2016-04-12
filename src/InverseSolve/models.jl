export expMod, boundMod, fMod,dummyMod


"""
	sigma,dsigma = expMod(model)
	
	maps model parameter to conductivity via
	
	sigma(m) = exp(m)
	
"""
function expMod(m)
	return exp(m), sdiag(exp(m))
end

"""
	sigma,dsigma = fMod(model;f::Function=identity,df::Function=m->speye(length(m)))
	
	maps model parameter to conductivity via

	sigma(m) = f(m) and dsigma(m) = sdiag(df(m))
	
"""
function fMod(m;f::Function=identity,df::Function=m->ones(length(m)))
	return f(m),sdiag(df(m))
end

function identityMod(m)
	return m, UniformScaling(1.0)
end


export boundMod
"""
	sigma,dsigma = boundMod(m;boundLow=0.0,boundHigh=1.0)
	
	maps model parameter to conductivity via
	
		sigma = 0.5*(boundHigh-boundLow) * (tanh(m)+1.0 ) + boundLow

"""
function boundMod(m; boundLow = 0.0, boundHigh = 1.0)
	u        = tanh(m)
	d        = 0.5 * (boundHigh - boundLow)
	sigma    = d * (u + 1.0) + boundLow
	dsigmadm = sdiag(d * (1.0 - u .* u))
	return sigma, dsigmadm
end

