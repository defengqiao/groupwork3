
#Hessian function
gethess = function(hess, grad, ..., n, theta, eps, gx) {
  # if hessian is not given, it will be estimated
  # by finite differencing of the gradient vector
  if (is.null(hess)) {
    H = matrix(0, n, n)
    for (i in 1:n) {
      di = theta
      di[i] = di[i] + eps
      H[i,] = (grad(di,...) - gx) / eps
    }
    H = (t(H) + H) / 2   # make estimated Hessian symmetric
  } else{
    # if hessian is given, it will be calculated
    H = hess(theta, ...)
  }
  return(H)
}

# Newton's Method
newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=20,eps=1e-6) {
  
  # calculate function, gradient for a given theta
  fx = func(theta, ...)
  gx = grad(theta, ...)
  # check if fx,gx is finite number
  if (!is.finite(fx)){
   stop("objective is not finite at initial parameters")
    stopifnot()
  }else if (!all(is.finite(gx))){
    stop("derivatives are not finite at initial parameters")
  }
  n = length(theta)     # number of optimization parameters
  ni = 0                # number of Newton iterations to try
  
  # iteration begin, determine if convergence 
  # and if number of iterations is too large
  while (any(abs(gx)>=tol*(abs(fx)+fscale))&(ni<maxit)) {
    #get Hessian matrix
    H=gethess(hess,grad,...,n=n,theta=theta,eps=eps,gx=gx)
    # if Hessian is not positive definite, perturb it
    # by adding M*I to it, where M is min(eigenvalue)+1
    # H2 is Cholesky factorization of Hessian
    while (mode((try(chol(H),silent=TRUE)->H2)) == "character") {
      lam = eigen(H)$values
      H = H + diag(1 - min(lam), n)
    }
    
    # step is a descent direction
    step = -chol2inv(H2) %*% gx  # chol2inv(H2) is H^-1
    
    # if step make value of fun bigger, halving step
    # if halving too many times, stop
    nh = 0              # halving times
    options(warn = -1)
    while ((!is.finite(func(theta + step))) | func(theta + step)>=fx) {
      if (nh <= max.half) {
        step = step / 2
        nh = nh + 1
      } else {
        options(warn = 1)
        return(warning("fail to reduce the objective despite trying ",
             max.half," step halvings"))
        
      }
    }
    options(warn = 1)
    # after finding a suitable step, new parameters is theta+step
    theta = theta + step
    # calculate calculate function, gradient for new theta
    fx = func(theta, ...)
    gx = grad(theta, ...)
    # finish iteration once
    ni = ni + 1
  }
  # iterations are finished
  
  # prepare for return
  result=list(f=fx,theta=theta,iter=ni,g=gx,Hi=NULL)
  # check if convergence
  if (all(abs(gx) < tol * (abs(fx) + fscale))){
    # check if Hessian is positive definite
    H=gethess(hess,grad,...,n=n,theta=theta,eps=eps,gx=gx)
    if (mode((try(chol(H),silent=TRUE)->H2)) == "character"){
      # not positive definite
      warning("Hessian is not positive definite at convergence")
    }else{
      # positive definite
      # calculate inverse of the Hessian
      result$Hi=chol2inv(H2)
    }
  }else{
    # fail to converge <= maxit times
    warning("fail to converge no more than ",maxit," times")
  }
  
  return(result)
  
}

