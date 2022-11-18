#rm(list = ls())

# test obj fun
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
  #log(x)
}
# test gradient fun
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
  #1/x
}
# test Hessian
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
  #-x^(-2)
}


newt(c(-25,100),rb,gb,hb)
hb(c(-25,100),k=2)

# test if fx, gx is not finite
newt(c(-25,-Inf),rb,gb,hb)

# test when H is not positive definiteness, because not positive, Hi=null, but still calculate
# test if iterations fail to converge
h1 <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- h[2,2] <- th[1]
  h[1,2] <- h[2,1] <- th[2]
  h
  #-x^(-2)
}
newt(c(-25,100),rb,gb,h1)
