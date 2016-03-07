export checkDerivative

function checkDerivative(f::Function,x0;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3)
# checkDerivative(f::Function,x0;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3)
	if out
		println(@sprintf("%9s\t%9s\t%9s\t%9s\t%9s\t%5s","h","E0","E1","O1","O2","OK?"))
	end
	v = randn(size(x0))
	if eltype(x0) <: Complex
		v += 1im*randn(size(v))
	end
	f0,dvf = f(x0,v)
	nf0 = norm(f0)
	#dvf   = real(dvf)
	Error = zeros(10,2)
	Order = zeros(10,2)
	Success = zeros(10)
	for j=1:10
		ft = f(x0+10.0^(-j)*v)            # function value
		Error[j,1] = norm(f0-ft)/nf0          # Error TaylorPoly 0
		Error[j,2] = norm(f0 .+10.0^(-j)*dvf .- ft)/nf0 # Error TaylorPoly 1
		if j>1
			Order[j,:] = log10(Error[j-1,:]./Error[j,:]);
		end
		if (Order[j,2]>tol) || (Error[j,1]/Error[j,2] > 100); Success[j]=1; end
		if out 
			println(@sprintf("%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%5d",
							10.0^(-j), Error[j,1],Error[j,2], Order[j,1],Order[j,2],Success[j]))
		end
	end
	pass = sum(Success) > nSuccess
	return  pass,Error,Order
end
