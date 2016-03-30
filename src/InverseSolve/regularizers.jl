export diffusionReg, wdiffusionReg, wTVReg, wdiffusionRegNodal,computeRegularizer,smallnessReg

function computeRegularizer(regFun::Function,mc::Vector,mref::Vector,MInv::AbstractMesh,alpha)
	R,dR,d2R = regFun(mc,mref,MInv)
	return alpha*R, alpha*dR, alpha*d2R
end	
	
function computeRegularizer(regFun::Array{Function},mc::Vector,mref::Array,MInv::AbstractMesh,alpha::Array{Float64})
	numReg = length(regFun)
	if size(mref,2)!=numReg; 
		error("computeRegularizer: number of regularizer (=$numReg) does not match number of reference models (=$(size(mref,2))).")
	end
	if length(alpha)!=numReg; 
		error("computeRegularizer: number of regularizer (=$numReg) does not match number of alphas (=$(length(alpha))).")
	end
	
	R = 0.0; dR = zeros(length(mc)); d2R = spzeros(length(mc),length(mc))
	for k=1:numReg
		Rt,dRt,d2Rt = regFun[k](mc,mref[:,k],MInv)
		R   += alpha[k]*Rt
		dR  += alpha[k]*dRt 
		d2R += alpha[k]*d2Rt
	end
	return R, dR, d2R
end


"""
	Rc,dR,d2R = diffusionReg(m,mref,M,Iact=1.0)
	
	Compute diffusion regularizer
		0.5*||GRAD*(m-mref)||_V^2
		
	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - projector on active cells
	
	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function diffusionReg(m::Vector,mref,M::AbstractMesh;Iact=1.0)
	# Rc = .5* || Grad*m ||^2
	dm   = m .- mref
	Div  = getDivergenceMatrix(M)
	Div  = Iact'*Div   # project to the active cells
	V    = getVolume(M)
	Af   = getFaceAverageMatrix(M)
	d2R  = Div * sdiag(Af'*vec(diag(V))) * Div'
	dR   = d2R*dm
	Rc   = 0.5*dot(dm,dR)
	return Rc,dR,d2R
end

"""
	Rc,dR,d2R = smallnessReg(m,mref,M,Iact=1.0)
	
	Compute smallness regularizer (L2 difference to reference model)
		
		R(m) = 0.5*||Iact*(m-mref)||_V^2
		
	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - regularization operator
	
	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function smallnessReg(m::Vector,mref,M::AbstractMesh;Iact=1.0)
	# Rc = .5* || Grad*m ||^2
	dm   = m .- mref
	d2R  = Iact'*Iact
	dR   = d2R*dm
	Rc   = 0.5*dot(dm,dR)
	return Rc,dR,d2R
end

function wdiffusionReg(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])
   # Rc = a1*\\Dx(m-mref)||^2 + .. + a3*\\Dz(m-mref)||^2 + a4*|| m -mref ||^2

   if length(C) == 0
      alpha1 = 1; alpha2 = 1; alpha3 = 1; alpha4 = 1e-4;
      Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])])
   elseif length(C) == 4
      # C = alphax, alphay, alphaz, alphas
      alpha1 = C[1]; alpha2 = C[2]; alpha3 = C[3]; alpha4 = C[4]
      Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])])
   elseif length(C) == 5
      #Af  = getFaceAverageMatrix(M)
      #wt  = Af'*C[5]
      wt   = C[5]
      Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])]) *
             sdiag(wt)
      #Wt = sdiag(C[5])
   elseif length(C) == sum(M.nf) + 1
      # C contains alphas followed by interface weights.
      alpha4 = C[1]
      Wt = sdiag( C[2:end] )
      #println("inteface weights used.")
   else
      error("bad regparams in wdiffusionReg")
   end
   
   dm   = m .- mref
   Div  = getDivergenceMatrix(M)

   V    = getVolume(M)
   #Div  = Iact'*(Div*Wt)   # project to the active cells
   Div  = Iact'*(V*Div*Wt)   # project to the active cells

   Af   = getFaceAverageMatrix(M)

   #d2R  = Div * sdiag(Af'*vec(diag(V))) * Div' + alpha4*Iact'*V*Iact
   d2R  = Div * sdiag(1./(Af'*vec(diag(V)))) * Div' + alpha4*Iact'*V*Iact

   dR   = d2R*dm
   Rc   = 0.5*dot(dm,dR)

   return Rc,dR,d2R
end # function wdiffusionReg

"""
	Rc,dR,d2R = wTVReg(m,mref,M,Iact,C=[])
	
	Compute weighted total variation regularizer
		
	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - projector on active cells
		C     - anisotropy parameters (default: [1 1 1])
		eps   - conditioning parameter for TV norm (default: 1e-3)
	
	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""
function wTVReg(m::Vector,mref,M::AbstractMesh; Iact=1.0, C=[],eps=1e-3)
	# Rc = sqrt(a1*\\Dx(m-mref)||^2 + .. + a3*\\Dz(m-mref)||^2) + a4*|| m -mref ||^2
	if length(C) == 0
		alpha1 = 1; alpha2 = 1; alpha3 = 1;
	elseif length(C) == 1
		alpha1 = 1; alpha2 = 1; alpha3 = 1;
	else
		alpha1 = C[1]; alpha2 = C[2]; alpha3 = C[3];
	end
	dm   = m .- mref
	Div  = getDivergenceMatrix(M)
	Wt   = sdiag([alpha1*ones(M.nf[1]);alpha2*ones(M.nf[2]);alpha3*ones(M.nf[3])])
	Div  = Iact'*(Div*Wt)   # project to the active cells
	V    = getVolume(M); v = diag(V)
	Af   = getFaceAverageMatrix(M)

	Rc   = dot(v,sqrt(Af*(Div'*dm).^2 .+eps))
	d2R  = Div*sdiag(Af'*(v./sqrt(Af*(Div'*dm).^2 .+eps)))*Div'
	dR   = d2R*dm;
	return Rc,dR,d2R
end



"""
	Rc,dR,d2R = wdiffusionRegNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])
	
	Computes weighted diffusion regularizer for nodal model
	
	Input:
		m     - model
		mref  - reference model
		M     - Mesh
		Iact  - projector on active cells
		C     - optional parameters
	
	Output
		Rc    - value of regularizer
		dR    - gradient w.r.t. m
		d2R   - Hessian
"""

function wdiffusionRegNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])	
	dm = m.-mref;
	if isempty(C)
		C = [1,1,1,1,1e-5];
	end
	if M.dim==3
		Wt   = sdiag([C[1]*ones((M.n[3]+1)*(M.n[2]+1)*(M.n[1]));C[2]*ones((M.n[3]+1)*(M.n[2])*(M.n[1]+1));C[3]*ones((M.n[3])*(M.n[2]+1)*(M.n[1]+1))]);
	else
		Wt   = sdiag([C[1]*ones((M.n[2]+1)*(M.n[1]));C[3]*ones((M.n[2])*(M.n[1]+1))]);
	end
	Grad   = Wt*getNodalGradientMatrix(M)*Iact
	d2R = Grad'*Grad;
	d2R += C[4]*speye(length(m));
	dR  = d2R*dm;
	Rc  = 0.5*dot(dm,dR);
	if isnan(Rc)
		dump(m)
	end
   return Rc,dR,d2R
end	 


# function wdiffusionRegNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])
	
	# G   = getNodalGradientMatrix(M)*Iact
	# d2R = G'*G
	# dR  = d2R*m
	# Rc  = 0.5*dot(m,dR)
   # return Rc,dR,d2R
# end	 
   