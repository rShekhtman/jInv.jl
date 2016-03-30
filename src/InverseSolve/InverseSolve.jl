module InverseSolve
	
	using KrylovMethods
	
	using jInv.Utils
	using jInv.Mesh
	using jInv.ForwardShare
	using jInv.LinearSolvers
	
	# import jInv.ForwardShare.ForwardProbType
	
	export getName,AbstractModel, AbstractMisfit

	abstract AbstractModel
	abstract AbstractMisfit
	
	include("HessianPreconditioners.jl")
	
	export MisfitParam
	"""
	type MisfitParam
		
	Type storing information about one term in the misfit
	
	F(m) = sum_i^n phi_i(pFor(model(m)),dobs,Wd)
	
	Fields:
		pFor::ForwardProbType  - forward problem
		Wd                     - inverse standard deviation
		dobs                   - observed data
		misfit::Function       - misfit function
		modelfun::Function        - model function (evaluated locally)
		gloc	               - mapping from inverse to forward mesh
	   
	Constructors:
		getMisfitParam(pFor,Wd,dobs,misfit,model,gloc=identity)
	"""
	type MisfitParam
		pFor::ForwardProbType
		Wd
		dobs
		misfit::Function
		modelfun::Function
		gloc	
	end
	
	export getMisfitParam
	function getMisfitParam(pFor::ForwardProbType, Wd, dobs, misfit::Function, modelfun::Function=fMod, gloc=identity)
		return MisfitParam(pFor,Wd,dobs,misfit,modelfun,gloc)
	end

	"""
	InverseParam
	
	Type storing parameters for Inversion. 
	
	Fields:
		Minv::AbstractMesh
		modelfun::Function    - model function (evaluated by main worker), see models.jl
		regularizer::Function - regularizer, see regularizer.jl
		alpha::Real           - regularization parameter
		regparam::Vector      - additional parameters for regularizer
		boundsLow::Vector     - lower bounds for model
		boundsHigh::Vector    - upper bounds for model
		maxStep::Real         - maximum step in optimization
    	pcgMaxIter::Int       - maximum number of PCG iterations
		pcgTol::Real          - tolerance for PCG
		minUpdate::Real       - stopping criteria
    	maxIter::Int          - maximum number of iterations
	
	Constructor:
		getInverseParam
	
	Example: 
		pInv = getInverseParam(model,Minv,Iact,modelfun,regularizer,alpha,mref,regparam,misfit,dobs,Wd,misparam)
	"""
	type InverseParam
		MInv::AbstractMesh
		modelfun::Function  # function that goes from model to conductivity
		regularizer::Union{Function,Array{Function}}  # function to calculate WTW
		alpha::Union{Float64,Array{Float64}}  # tradeoff parameter
		mref::Array  # reference model
		boundsLow::Vector
		boundsHigh::Vector
		maxStep::Real  # maximum step size for delta m
		pcgMaxIter::Int
		pcgTol::Real
		minUpdate::Real  #  Step norm stopping tol.
		maxIter::Int
		HesPrec
	end  # type InverseParam
	
	"""
	getInverseParam(...)
	
	Constructs an InverseParam
	
	Required Input:
	
		Minv::AbstractMesh    - mesh of model
		regularizer::Function - regularizer, see regularizer.jl
		alpha::Real           - regularization parameter
		boundsLow::Vector     - lower bounds for model
		boundsHigh::Vector    - upper bounds for model
		
	Optional Inputs:
	
		maxStep::Real=1.0     - maximum step in optimization
	    pcgMaxIter::Int=10    - maximum number of PCG iterations
	    pcgTol::Real          - tolerance for PCG
		minUpdate::Real=1e-4  - stopping criteria
	    maxIter::Int=10       - maximum number of iterations
	"""
	function getInverseParam(MInv,modFun,
							 regularizer,alpha,mref,
							 boundsLow::Vector,boundsHigh::Vector; maxStep=1.0,pcgMaxIter=10,pcgTol=1e-1,
							 minUpdate=1e-4,maxIter=10,HesPrec=getSSORRegularizationPreconditioner(1.0,1e-15,10))
							 
		return 	InverseParam(MInv,modFun,
							 regularizer,alpha,mref,
		                     boundsLow,boundsHigh,
							 maxStep,pcgMaxIter,pcgTol,minUpdate,maxIter,HesPrec)
	end



	export InverseParam, getInverseParam

	include("globalToLocal.jl")
	include("models.jl")
	include("misfit.jl")
	include("regularizers.jl")
	include("projGNCG.jl")
	include("projPCG.jl")
	include("computeMisfit.jl")
	include("computeGradMisfit.jl")
	include("HessMatVec.jl")
	
end
