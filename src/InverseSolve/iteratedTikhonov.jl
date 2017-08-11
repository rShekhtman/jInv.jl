export iteratedTikhonov

"""
    mc,DC = iteratedTikhonov(mc,pInv::InverseParam,pMis,nAlpha,alphaFac,
                            targetMisfit;indCredit=[],dumpResults::Function=dummy)
    
    Perform (Projected) Gauss-NewtonCG using iterated Tikhonov procedure to decrease
    regularization parameter and update reference model after fixed number of GN iterations
    set in pInv.
    
    Input:
           mc::Vector         - Initial guess for model
           pInv::InverseParam - parameter for inversion
           pMis               - misfit terms
           nAlpha             - maximum number of allowed regularization parameter (alpha) values
           alphaFac           - alpha decrease factor. alpha_(i+1) = alpha_i/alphaFac
           targetMisfit       - Termination criterion. iteratedTikhonov will exit 
                                 -when data misfit < targetMisfit
           indCredit          - indices of forward problems to work on
           dumpResults        - A function pointer for saving the results throughout the iterations.
                                 - We assume that dumpResults is dumpResults(mc,Dc,iter,pInv,pMis), 
                                 - where mc is the recovered model, Dc is the predicted data. 
                                 - If dumpResults is not given, nothing is done (dummy() is called).
           
    Output:
           mc                 - final model
           Dc                 - final computed data
           tikhonovFlag       - Data misfit convergence flag
           hist               - Iteration history. This is a vector. Each entry is
                                 - a structure containing the projGNCG history for each
                                 - alpha value.
"""

function iteratedTikhonov(mc,pInv::InverseParam,pMis,nAlpha,alphaFac,
                          targetMisfit;indCredit=[],dumpResults::Function=dummy,solveGN=projPCG)

  if alphaFac <= 1 ; error("alphaFac <= 1") ; end
  changeMref = true  # Should the reference model change after each iteration?
  iter = 0
  tikhonovFlag = -1
  hist = GNhis[]
  Dc = []
  cm = computeMisfitOutput([],[],[],[],[])  # remember results from previous iteration.
  
  while iter < nAlpha 
    iter += 1
    println("Starting projGNCG minimization with alpha $iter of $nAlpha")
    println("alpha = $(pInv.alpha)")
    mc,Dc,GNflag,GNhist, cm = projGN(mc,pInv,pMis,dumpResults=dumpResults, outputcm=true, cm=cm)

    dumpResults(mc,Dc,iter+1000,pInv,cm.pMis)

    push!(hist,GNhist)
    if changeMref
       pInv.mref  = copy(mc)
    end
    pInv.alpha = pInv.alpha / alphaFac
    GNiterIdx  = find(hist[end].F .> 0.0)
    println(hist[end].F)
    if minimum(hist[end].F[GNiterIdx]) <= targetMisfit
      tikhonovFlag = 1
      break
    end
  end  # while iter < nAlpha 
  
  if tikhonovFlag < 0 
    println("iteratedTikhonov iterated maximum number of times, using $nAlpha alpha values but reached only a misfit of $(hist[end].F[end]).")
  else
    println("iteratedTikhonov exiting after reaching desired misfit")
  end
  return mc,Dc,tikhonovFlag,hist
end