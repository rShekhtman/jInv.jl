
function Base.display(Mesh::TensorMesh3D)
	
	n  = Mesh.n
	nc = Mesh.nc
	nf = Mesh.nf
	ne = Mesh.ne
	nn = prod(Mesh.n+1)
	x0 = Mesh.x0
	h1 = Mesh.h1
	h2 = Mesh.h2
	h3 = Mesh.h3
	
	println("Tensor mesh of size $(n[1]) x $(n[2]) x $(n[3])")
	println("Number of cells:   $nc")
	println("Number of faces:   $(nf[1]) + $(nf[2]) + $(nf[3]) = $(sum(nf))")
	println("Number of edges:   $(ne[1]) + $(ne[2]) + $(ne[3]) = $(sum(ne))")
	println("Number of nodes:   $nn")
	println("Coordinate origin: ($(x0[1])m, $(x0[2])m, $(x0[3])m)")
	println("Domain size:       $(sum(h1))m x $(sum(h2))m x $(sum(h3))m")
	println("Minimum cell size: $(minimum(h1))m x $(minimum(h2))m x $(minimum(h3))m")
	println("Maximum cell size: $(maximum(h1))m x $(maximum(h2))m x $(maximum(h3))m")
	
end
