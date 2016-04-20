using Base.Test

allPassed = true
try 
	include("Mesh/runtests.jl")
catch
	allPassed = false
	warn("Mesh had test errors")
end
try 
	include("ForwardShare/runtests.jl")
catch
	allPassed = false
	warn("ForwardShare had test errors")
end

try 
	include("Utils/runtests.jl")
catch
	allPassed = false
	warn("Utils had test errors")
end

try 
	include("InverseSolve/runtests.jl")
catch
	allPassed = false
	warn("InverseSolve had test errors")
end

@test allPassed == true

