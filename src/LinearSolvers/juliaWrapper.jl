export JuliaSolver,getJuliaSolver,copySolver

type JuliaSolver<: AbstractSolver
	Ainv
	sym::Int # 0 = unsymmetric, 1 = symmetric s.t A = A';
	isTransposed::Int
	doClear::Int
	facTime::Real
	nSolve::Int
	solveTime::Real
	nFac::Int
end

function getJuliaSolver(;Ainv = [],sym = 0,isTransposed = 0, doClear = 0)
	return JuliaSolver(Ainv,sym,isTransposed,doClear,0.0,0,0.0,0);
end

solveLinearSystem(A,B,param::JuliaSolver,doTranspose::Int=0) = solveLinearSystem!(A,B,[],param,doTranspose)

function solveLinearSystem!(A,B,X,param::JuliaSolver,doTranspose=0)
	if param.doClear == 1
		clear!(param)
	end
	if param.sym==0
		if doTranspose==1 && param.isTransposed==0
			clear!(param);
			A = A';
			param.isTransposed = 1;
		end
		if doTranspose==0 && param.isTransposed==1
			clear!(param);
		end
	end
	if param.Ainv == []
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
function clear!(param::JuliaSolver)
	param.Ainv = [];
	param.isTransposed = 0;
	param.doClear = 0;
end
	
function copySolver(Ainv::JuliaSolver)
	return getJuliaSolver([],Ainv.sym,0,Ainv.doClear);
end