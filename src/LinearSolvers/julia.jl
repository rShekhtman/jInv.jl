type JuliaSolver<: AbstractSolver
	Ainv
end

function getJuliaSolver(Ainv=[])
	return JuliaSolver(Ainv)
end

function solveLinearSystem(A,b,param::JuliaSolver,doTranspose)
	if doTranspose == 0
		return A\b
	else
		return (A')\b
	end
end
