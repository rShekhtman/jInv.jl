
println(" test module InverseSolve")

include("testModels.jl")
include("testMisfit.jl")
include("testRegularizers.jl")
include("testMisfitParams.jl")
include("testLeastSquares.jl")
include("testGetHessian.jl")
include("testHessMatVec.jl")
# include("FWI/runtests.jl")
# include("DivSigGrad/runtests.jl")
println(" InverseSolve: All tests passed!")
