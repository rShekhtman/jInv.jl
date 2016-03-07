using jInv.Mesh
using Base.Test

# setup 3D problem
nc = [5, 7, 8]
x0 = rand(3)
domain = [x0[1], 4, x0[2], 2, x0[3], 6]
h   = (domain[2:2:end]-domain[1:2:end])./nc
h1  = h[1]*ones(nc[1])
h2  = h[2]*ones(nc[2])
h3  = h[3]*ones(nc[3])

Mt = getTensorMesh3D(h1,h2,h3,x0)
Mr = getRegularMesh(domain,nc)
Mt2 = getTensorMesh3D(h1+rand(nc[1]),h2+rand(nc[2]),h3+rand(nc[3]),x0)

print("\ttest DIV*CURL==0...")
Dr = getDivergenceMatrix(Mr)
Cr = getCurlMatrix(Mr)
@test norm(Dr*Cr,Inf)<1e-12

Dt = getDivergenceMatrix(Mt)
Ct = getCurlMatrix(Mt)
@test norm(Dt*Ct,Inf)<1e-12

Dt2 = getDivergenceMatrix(Mt2)
Ct2 = getCurlMatrix(Mt2)
@test norm(Dt2*Ct2,Inf)<1e-12
print("passed!\n")


