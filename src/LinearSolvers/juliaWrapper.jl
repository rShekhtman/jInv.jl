type JuliaSolver<: AbstractSolver
	Ainv
	facTime::Real
	nSolve::Int
	solveTime::Real
end

function getJuliaSolver()
	return JuliaSolver([],0.0,0,0.0)
end

function solveLinearSystem!(A,B,X,param::JuliaSolver,doTranspose=0)
	if param.doClear == 1
		clear!(param)
	end
	if doTranspose
		clear!(param);
		A = A';
	end
	if param.Ainv==[]
		tic()
		param.Ainv = lufact(A)
		param.facTime+=toq()
		param.nFac+=1
	end
	
	tic()
	U = param.Ainv\B;
	param.solveTime+=toq()
	param.nSolve+=1
	
	return U, param		
end # function solveLinearSystem
			
import jInv.Utils.clear!
function clear!(M::JuliaSolver)
	M.Ainv = [];
end
	
function copySolver(Ainv::JuliaSolver)
	return getJuliaSolver();
end
