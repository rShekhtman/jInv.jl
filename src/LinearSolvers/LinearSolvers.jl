module LinearSolvers
	
	abstract AbstractSolver
	export   AbstractSolver

	using KrylovMethods
	
	# check if MUMPS can be used
	const minMUMPSversion = VersionNumber(0,0,1)
	hasMUMPS=false
	try 
		using MUMPS
		hasMUMPS = true
		if myid()==1
			vMUMPS = Pkg.installed("MUMPS")
			if vMUMPS < minMUMPSversion; 
				warn("MUMPS is outdated! Please update to version $(minMUMPSversion) or higher.")
			end
		end
	catch 
	end
	
	# check if Pardiso is installed
	hasPardiso = false
	try
		using Pardiso
		hasPardiso = true
	catch
	end
	
	# check if ParSPMatVec is available
	hasParSpMatVec = false
	const minVerParSpMatVec = VersionNumber(0,0,1)
	try 	
		using ParSpMatVec
		hasParSpMatVec = true
		if myid()==1
			vParSpMatVec = Pkg.installed("ParSpMatVec")
			if vParSpMatVec < minVerParSpMatVec; 
				warn("ParSpMatVec is outdated! Please update to version $(minVerParSpMatVec) or higher.")
			end
		end
	catch 
	end
	
	export solveLinearSystem!,solveLinearSystem
	
	solveLinearSystem(A,B,param::AbstractSolver,doTranspose::Int=0) = solveLinearSystem!(A,B,zeros(eltype(B),size(B)),param,doTranspose)

	include("mumpsWrapper.jl")
	include("iterativeWrapper.jl")
	include("blockIterativeWrapper.jl")
	include("PardisoWrapper.jl")
	
	import jInv.Utils.clear!

	function clear!(M::AbstractSolver)
		M.Ainv = []
	end
	
	export clear!

end # module LinearSolvers
