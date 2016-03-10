export expMod, boundMod, idMod


"""
	sigma,dsigma = expMod(model)
	
	maps model parameter to conductivity via
	
	sigma(m) = exp(m)
	
"""
function expMod(m)
    sigma    = exp(m)
	dsigmadm = sdiag(sigma)
    return sigma, dsigmadm
end

"""
	sigma,dsigma = idMod(model)
	
	maps model parameter to conductivity via

	sigma(m) = m
	
"""
function idMod(m)
	return m,speye(length(m))
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

