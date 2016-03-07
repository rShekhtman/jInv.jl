export sdiag

function sdiag(a::Vector)
# S = sdiag(s) builds sparse diagonal matrix
	n = length(a)
	i = collect(1:n+1) # colptr
	j = collect(1:n)   # rowval
	return SparseMatrixCSC(n,n,i,j,a)
end