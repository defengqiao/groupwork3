rm(list = ls())

# 几种求step的效率比较
A=matrix(1,10,10)
A=A+diag(1,10,10)
B=c(1,2,3,4,5,6,7,8,9,10)

a=system.time({for (i in 1:200000) {
  Ai=chol2inv(chol(A))
  sol1=Ai%*%B
}
})

b=system.time({for (i in 1:200000) {
  Ac=chol(A)
  sol2=solve(Ac,solve(t(Ac),B))
}
})

c=system.time({for (i in 1:200000) {
  Aci=solve(chol(A))
  sol3=Aci%*%t(Aci)%*%B
}
})
print(a)
print(b)
print(c)

# 非正定
A1=matrix(1,3,3)
A2=A1+diag(-3,3)
A1=A1+diag(3,3)
b1=try(chol(A2),silent = TRUE)
b2=try(chol(A1),silent = TRUE)
print(mode((try(chol(A2),silent = TRUE)->H2))=="character")
print(mode(b2)=="numeric")

# 测试判断语句中是否能进行赋值
i=1
while ((i=i+1)<=10){
  print(i)
}
# 这样proj4中可以写成
#A1=matrix(1,3,3)
#A2=A1+diag(-3,3)
#H=A1+diag(3,3)   #假设是Hessian given theta
while (mode((try(chol(H),silent = TRUE)->H2))=="character"){
  lam=eigen(H)$values
  H=H+diag(1-min(lam),nrow(A))
}

# 步长的计算
# chol2inv(H2) is H^-1, g is the vector of gradient given theta
step=-chol2inv(H2)%*%g

# estimate Hessian by finite differencing of the gradient vector
H=matrix(0,n,n)
for (i in 1:n){
  di=theta
  di[i]=di[i]+eps
  H[i,]=(g-grad(di))/eps
}
H=(t(H)+H)/2      # symmetric

# stop and warning
iter <- 12
try(if(iter > 10) stop("too many iterations"))

f=function(x){
  y=1/x
  if (y==Inf) {
    stop("function cannot be evaluated at initial ",x," parameters")
  }
  print(y)
}

try(f(0))


#Hessian function
gethess = function(hess, grad, ..., n, theta, eps, gx) {
  # if hessian is not given, it will be estimated
  # by finite differencing of the gradient vector
  if (is.null(hess)) {
    H = matrix(0, n, n)
    for (i in 1:n) {
      di = theta
      di[i] = di[i] + eps
      H[i,] = (gx - grad(di, ...)) / eps
    }
    H = (t(H) + H) / 2   # make estimated Hessian symmetric
  } else{
    # if hessian is given, it will be calculated
    H = hess(theta, ...)
  }
  return(H)
}





#
foo=function(x){
  log(x)
}
optim(par=-1,fn=foo)



