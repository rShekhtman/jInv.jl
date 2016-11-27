using Base.Test

println("\t\tcheckDerivative for ForwardProbTypes")
include("../setupTests.jl")
A    = sprandn(100,10,.1)
pFor = LSparam(A,[])
m0   = randn(10)

passed, = checkDerivative(m0,pFor)

println("\t\tadjoint test")

passed, = adjointTest(m0,pFor)
@test passed

@test passed