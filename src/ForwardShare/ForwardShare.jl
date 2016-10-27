"""
jInv's module for forward operators

To define your own module, you must at least

1) Create your own type (documented param that describes a forward problem)

   		type MyParam{M}
      		M::T
      		Sources::Array{Complex128}
      		Obs::SparseMatrixCSC
      		getData::Function
      		getSensMatVec::Function
      		getSensTMatVec::Function
      		Ainv::AbstractSolver
	    		fname::AbstractString
	  	end

2) Code a method called 'computeFwd' for your type that solves the forward problem and has
	    the interface

	    function getData(sigma,param::YourType)
	    	return fields
	    end

3) Code methods for matrix vector products involving the Sensitivities. Your methods should look like

      function getSensMatVec(x,sigma,param::YourType)
   			return Sens*x
      end

      function getSensTMatVec(x,sigma,param::YourType)
      	return transpose(Sens)*x
      end

	4) clear method for your type
"""
module ForwardShare


	export ForwardProbType
	abstract ForwardProbType

	using jInv.Mesh
	using jInv.Utils

	export getSensMatVec
	"""
	Jv  = getSensMatVec(v::Vector,m::Vector,param::ForwardProbType)

	Computes matrix-vector product with the Jacobian.

	"""
	getSensMatVec(v::Vector,m::Vector,param::ForwardProbType) = error("nyi")
	export getSensTMatVec
	"""
	JTv  = getSensMatVec(v::Vector,m::Vector,param::ForwardProbType)

	Computes matrix-vector product with the transpose of Jacobian. Implementation
	depends on forward problem.

	"""
	getSensTMatVec(v::Vector,m::Vector,param::ForwardProbType) = error("nyi")

	"""
	(m,n) = getSensMatSize(pFor)
	
	Returns size of sensitivity matrix where m is the number of data points
	and n the number of parameters in the model. 
	 
	Input
	
		pFor - forward problem:: Union{ForwardProbType, Array, RemoteRef}
	
	This is problem dependent and should be implemented in the respective
	packages. 
	"""
	getSensMatSize(pFor::ForwardProbType) = error("nyi")
	
	function getSensMatSize(pFor::RemoteRef)
		if pFor.where != myid()
			return remotecall_fetch(pFor.where,getSensMatSize,pFor)
		else
			return getSensMatSize(fetch(pFor))
		end
	end
	
	function getSensMatSize(pFor::Array)
		n  = length(pFor)
		sz = [0;0]
		for k=1:n
			(s1,s2) = getSensMatSize(pFor[k])
			sz[2]=s2
			sz[1]+=s1
		end
		return tuple(sz...)
	end
	
	export getSensTMatVec,getSensMatVec, getSensMatSize

	# # ===== Methods for parallelization =====
	include("getDataParallel.jl")
	include("prepareMesh2Mesh.jl")
	include("interpLocalToGlobal.jl")
	include("getSensMat.jl")

	import jInv.Utils.clear!
	function clear!(P::ForwardProbType;clearAinv::Bool=true,clearFields::Bool=true, clearMesh::Bool=false, clearSources::Bool=false, clearObs::Bool=false,clearAll::Bool=false)
		if clearAll || clearMesh
			Utils.clear(P.M)
		end
		if clearAll || clearSources
			P.Sources = clear(P.Sources)
		end
		if clearAll || clearObs
			P.Obs     = clear(P.Obs)
		end
		if clearAll || clearFields
			P.Fields = clear(P.Fields)
		end
		if clearAll || clearAinv
			clear!(P.Ainv)
		end
	end


end
