module InverseSolve
	
	using KrylovMethods
	
	using jInv.Utils
	using jInv.Mesh
	using jInv.ForwardShare
	using jInv.LinearSolvers
	
	# import jInv.ForwardShare.ForwardProbType
	
	export getName,AbstractModel, AbstractMisfit, AbstractRegularizer

	abstract AbstractModel
	abstract AbstractMisfit
	abstract AbstractRegularizer
	
	include("HessianPreconditioners.jl")

	"""
	InverseParam
	
	Type storing parameters for Inversion. 
	
	Fields:
	
		model
		Minv::AbstractMesh
		Iact                  - active cells
		modelfun::Function    - model for conductivities, see models.jl
		regularizer::Function - regularizer, see regularizer.jl
		alpha::Real           - regularization parameter
		regparam::Vector      - additional parameters for regularizer
		misfit::Function      - misfit function
		dobs                  - observed data
		Wd                    - 1./standardDev of noise
		misparam              - additional parameters for misfit
		boundsLow::Vector     - lower bounds for model
		boundsHigh::Vector    - upper bounds for model
		maxStep::Real         - maximum step in optimization
    	pcgMaxIter::Int          - maximum number of PCG iterations
		pcgTol::Real          - tolerance for PCG
		minUpdate::Real       - stopping criteria
    	maxIter::Int          - maximum number of iterations
	
	Constructor:
		getInverseParam
	
	Example: 
		pInv = getInverseParam(model,Minv,Iact,modelfun,regularizer,alpha,mref,regparam,misfit,dobs,Wd,misparam)
	"""
	type InverseParam
		model
		MInv::AbstractMesh
		Iact  # active cells
		modelfun::Function  # function that goes from model to conductivity
		regularizer::Union{Function,Array{Function}}  # function to calculate WTW
		alpha::Union{Float64,Array{Float64}}  # tradeoff parameter
		mref::Array  # reference model
		regparam::Array  # alpha values and/or weights for regularization
		misfit::Function  # calculate data misfit and derivarives
		dobs  # observed data
		Wd    # 1/standard deviation
		misparam
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
	
		model
		Minv::AbstractMesh
		Iact                  - active cells
		modelfun::Function    - model for conductivities, see models.jl
		regularizer::Function - regularizer, see regularizer.jl
		alpha::Real           - regularization parameter
		regparam::Vector      - additional parameters for regularizer
		misfit::Function      - misfit function
		dobs                  - observed data
		Wd                    - 1./standardDev of noise
		misparam              - additional parameters for misfit
		boundsLow::Vector     - lower bounds for model
		boundsHigh::Vector    - upper bounds for model
		
	Optional Inputs:
	
		maxStep::Real=1.0     - maximum step in optimization
	    pcgMaxIter::Int=10    - maximum number of PCG iterations
	    pcgTol::Real          - tolerance for PCG
		minUpdate::Real=1e-4  - stopping criteria
	    maxIter::Int=10       - maximum number of iterations
	"""
	function getInverseParam(model,MInv,Iact,modelfun,
							 regularizer,alpha,mref,regparam,
		                     misfit,dobs,Wd,misparam,
							 boundsLow::Vector,boundsHigh::Vector;
							 maxStep=1.0,pcgMaxIter=10,pcgTol=1e-1,minUpdate=1e-4,maxIter=10,HesPrec=getSSORRegularizationPreconditioner(1.0,1e-15,10))
							 
		return 	InverseParam(model,MInv,Iact,modelfun,
							 regularizer,alpha,mref,regparam,
		                     misfit,dobs,Wd,misparam,
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
