module LinearSolvers
	
	abstract AbstractSolver
	export   AbstractSolver

	using KrylovMethods
	
	# check if MUMPS can be used
	hasMUMPS=false
	try 
		using MUMPS
		hasMUMPS = true
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
		vParSpMatVec = VersionNumber(0,0,0)
	try 
		vParSpMatVec = Pkg.installed("ParSpMatVec")
		hasParSpMatVec = vParSpMatVec>=minVerParSpMatVec
	catch 
	end
	if hasParSpMatVec
		using ParSpMatVec
	end
	
	export solveLinearSystem!,solveLinearSystem
	
	solveLinearSystem(A,B,param::AbstractSolver,doTranspose::Int=0) = solveLinearSystem!(A,B,zeros(eltype(B),size(B)),param,doTranspose)

	include("mumpsWrapper.jl")
	# include("pcgWrapper.jl")
	include("iterativeWrapper.jl")
	include("blockcg.jl")
	include("PardisoWrapper.jl")
	# include("julia.jl")
	
	import jInv.Utils.clear!

	function clear!(M::AbstractSolver)
		M.Ainv = []
	end
	
	export clear!

end # module LinearSolvers
