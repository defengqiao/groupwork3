#Defeng Qiao, s2419769; Tianai Ren, s2329207; YiZhou Chen, s2450877
#
#Contributions:Defeng Qiao: 
#              Tianai Ren:  
#              YiZhou Chen: 

# Overview: 
# Code to perform newton's method for minimization for functions.
# The challenge is that in some situations the Hessian matrix is not given,
# thus, the we need to write a function to create a Hessian function in 
# this situation. 
# In addition to that, the function would not work unless fx and gx provide
# finite value, and warnings would be called.
# Moreover, if the hessian matrix is not positive definite, an obvious approach 
# is to add a multiple of the identity matrix to it, large enough to force 
# positive definiteness, when everything is ready, the newton's method would 
# be applied in the objective function and start the optimazation process.
# However, it would issue warning if halving too many times of iteration which 
# means it larger than max.half and number of Newton iterations is larger than 
# maximum number of Newton iterations.
# Last but not least, the convergence should be judged by seeing whether all 
# elements of the gradient vector have absolute value less than tol times the 
# absolute value of the objective function plus fscale, and another warning 
# would be called when hessian is not positive definite at convergence.
# Finally, the code here would investigates the minimize value of the objective 
# function, the theta value in the optimal place, time of iteration, the final 
# gradient function and the inverse ofthe hessian matrix.


# Hessian function
gethess = function(hess, grad, ..., n, theta, eps, gx) {
  # if hessian is not given,estimate it by finite differencing of the gradient vector.
  # hess is the Hessian matrix function, grad is the gradient function
  # ... is used to pass arguments to a function, n is the number of optimization parameters
  # theta is a vector of initial values, eps is the finite difference intervals
  if (is.null(hess)) { # if hessian is not given
    H = matrix(0, n, n) # initializate hessian
    for (i in 1:n) {
      di = theta 
      di[i] = di[i] + eps 
      H[i,] = (grad(di,...) - gx) / eps #calculate the column vector of ith hessian matrix
    }
    H = (t(H) + H) / 2   # estimated Hessian is symmetric
  } else{
    H = hess(theta, ...) # if hessian is given, output will be calculated directly
  }
  return(H) # return Hessian function
} #gethess

# Newton's Method
newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                maxit=100,max.half=20,eps=1e-6) {
  # Implementing Newton's method for minimization of functions.
  # theta is a vector of initial values, func is the objective function to minimize
  # grad is the gradient function, hess is the Hessian matrix function, ... is used to pass arguments to a function,
  # tol the convergence tolerance, fscale a rough estimate of the magnitude of func near the optimum,
  # maxit the maximum number of Newton iterations, max.half is the maximum number of times a step should be halved,
  # eps is the finite difference intervals.
  
  fx = func(theta, ...) # calculate the objective function at the initial theta
  gx = grad(theta, ...) # calculate the gradient at the initial theta
  
  if ((!is.finite(fx))|(!all(is.finite(gx)))){ # check if fx,gx are finite 
    stop("function is not finite at initial parameters") # issue warnings if fx or gx is not finite
  }
  
  n = length(theta)     # number of optimization parameters
  ni = 0                # number of Newton iterations
  
  # iteration begin, determine if convergence whether all elements of the gradient vector 
  # have the absolute value less than tol times the absolute value of the objective function plus fscale.
  # and if number of iterations is less than the maximum number of Newton iterations
  while (any(abs(gx)>=tol*(abs(fx)+fscale))&(ni<maxit)) { 
    H=gethess(hess,grad,...,n=n,theta=theta,eps=eps,gx=gx) 
    # We also have to ensure that the approximating quadratic has a proper minimum
    # if Hessian is not positive definite, perturb it 
    # by adding M*I to it, where M is min(eigenvalue)+1
    # positive definiteness is tested by seeing if a Cholesky decomposition is possible
    while (mode((try(chol(H),silent=TRUE)->H2)) == "character") { # H2 is Cholesky factorization of Hessian
      #denote the eigrnvalue of Hessian as lam
      lam = eigen(H)$values
      #when hessian is not positive definite, then at least one of the pivot is negative, then eliminatie it by subtraction
      H = H + diag(1 - min(lam), n)
    }
    
    # step is a descent direction which is a sufficiently small step will decrease func
    step = -chol2inv(H2) %*% gx  # chol2inv(H2) is H^-1, step is -H^-1*grad
    
    # if step make value of fun bigger func(theta+step)>func(theta), halving step
    # until a suitable step make func(theta+step)<func(theta)
    # issue warning if halving too many times which means larger than max.half
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
   
    theta = theta + step  # after finding a suitable step, new parameters is theta+step
    
    fx = func(theta, ...) # calculate the objective function for new theta
    gx = grad(theta, ...) # calculate the gradient for new theta
    
    ni = ni + 1 # count the number of time of iteration
  }# iterations are finished
  
  # newt function should return a list containing:
  # the value of the objective function(f), the parameters(theta) at the minimum f
  # the number of iterations(iter), the gradient vector(g) reach the minimum
  # the inverse of the Hessian matrix at the minimum(Hi)
  result=list(f=fx,theta=theta,iter=ni,g=gx,Hi=NULL)
  if (all(abs(gx) < tol * (abs(fx) + fscale))){ # check if convergence
    # check if Hessian is positive definite
    H=gethess(hess,grad,...,n=n,theta=theta,eps=eps,gx=gx)
    if (mode((try(chol(H),silent=TRUE)->H2)) == "character"){
      # issue warning if Hessian is not positive definite at convergence
      warning("Hessian is not positive definite at convergence")
    }else{
      # if Hessian is positive definite, calculate the inverse of it
      result$Hi=chol2inv(H2) 
    }
  }else{
    # issue warning if iterations fail to converge
    warning("fail to converge no more than ",maxit," times")
  }
  
  return(result) # return a list
  
} # newt

