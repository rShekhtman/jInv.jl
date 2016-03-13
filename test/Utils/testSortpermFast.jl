using jInv.Utils
using Base.Test


v = unique(rand(1:200,54))
w = randn(100)

data = (v,w)
println("testing sortpermFast")

for k=1:length(data)
	
	res1 = sortpermFast(data[k])
	res2 = sortperm(data[k])
	b2   = data[k][res2]
	
	@test all(res1[1] .== res2)
	@test all(res1[2] .== b2)
end