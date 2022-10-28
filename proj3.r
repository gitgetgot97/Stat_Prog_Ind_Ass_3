# Scott Jenkins s2451950 
# Practical 3: Smoothing with basis expansions and penalties

## This program ....

### Deadline 12:00 Friday 4th November 2022

## x,y vectors to smooth
## k is the number of basis functions to use
## logsp are the ends of the interval over which to search for smoothing parameters (log lambda scale).
## If only one provided, then no searching, spline returned for the given log lambda value
## bord is the order of b (B-spline order): 3 corresponds to cubic
## pord is the order of difference to use in the penalty. 2 is used in the example sheet
## ngrid is the number of smoothing parameter values to try. Use even spacing on a log scale.

## Sources
# https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-019-0666-3

library(MASS)
x <- mcycle$times ## set x = mcycle times column, y = mcycle accel column
y <- mcycle$accel
# x <- x[1:5]; y <- y[1:5] # Show the first 5 elements of x and y

# Plot the x and y data - we can see it's non-linear!
par(mfrow=c(1,2))
plot(x,y)
model <- lm(y~x)
abline(coef(model), col=2)
plot(predict(model), residuals(model))


# Try to understand knots and X
# k <- 20; bord <- 3; pord <- 2
# dk <- diff(range(x)/(k-bord)) 
# length(dk)
# dk is a real number. range(x) returns the smallest and largest x values
# q: why divide by k-bord instead of just k. k is the number of basis functions.
# 
# knots <- seq(min(x)-dk*bord,by=dk,length = k+bord+1)
# length(knots) # sequence starts with min(x)-dk*bord, increments by dk, legnth is k+bord+1
# X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
# dim(X) # has length(x) rows, and k columns 
# D <- diff(diag(k),differences=pord)
# dim(D)

# logsp <- c(-5,5); ngrid <- 100
# 
# lambda_ss_range <- max(logsp)-min(logsp)
# lambda_search <- seq(min(logsp),by=lambda_ss_range/(ngrid-1),length = ngrid)
# lambda_search[1]

###############################################################################

## Q1. pspline

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
  ## Uses the splineDesign function to fit Pslines to x,y data.
  ## Chooses the smoothing parameter lambda from the logsp which minimises GCV
 
  ## Find n, the length of the data
  n <- length(x)
  
  ## Produce the knots, X, and D
  dk <- diff(range(x))/(k-bord) ## knot spacing
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  print("knots : "); print(knots)
  
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)  # k x k 
  print("X : "); print(X)
  
  D <- diff(diag(k),differences=pord)    # k-2 x k
  print("D : "); print(D)
  
  ## QR Decomposition of X to find R and Q^T*y
  qrx <- qr(X)
  R <- qr.R(qrx) # R # k x k upper diagonal
  R_inv <- solve(R) # R_inverse # k * k  upper diagonal
  print("R, R_inv : "); print(R); print(R_inv)
  QT_y <- qr.qty(qrx,y)  # Q^T*y  # 1 x n
  print("QT*y : "); print(QT_y)
  
  # ## Use eigen decomposition to find U and (capital) Lambda.
  r <- eigen(t(R_inv) %*% crossprod(D) %*% R_inv) 
  U <- r$vectors; # k x k 
  print('U : '); print(U)
  Lambda <- r$values   # 1 x k
  print('Lambda : '); print(Lambda)

  ## Set up the search space for lambda.
  if (length(logsp) == 2) {
    lambda_ss_range <- max(logsp)-min(logsp)
    lambda_search <- seq(min(logsp),by=lambda_ss_range/(ngrid-1),length = ngrid)
  } else {
    lambda_search <- logsp
  }
  # print(lambda_search)
  
  ## For each lambda in our search space, compute kappa, beta_hat, mu_hat, sigma^2_hat and GCV
  lambda <- logsp # Just try 1 value to begin
  
  diag_inv <- solve(diag(lambda*Lambda + 1, k, k)) ## inverse of diagonal matrix: I + lambda*Lambda
  print('Diagonal Inverse : '); print(diag_inv)

  kappa <- sum(diag(diag_inv))  # trace of the inverse of diagonal, a real number
  print("kappa : "); print(kappa)    
  
  # Naively, can compute beta_hat as...
  beta_hat <- solve(crossprod(X) + lambda*crossprod(D)) %*% t(X) %*% y    
  print("Slow Beta Hat : "); print(beta_hat) # But this takes O(k^3) time. We can do it in O(k) as follows...
  # beta_hat <- R_inv %*% U %*% diag_inv %*% t(U) %*% QT_y
  # print("Fast Beta Hat : "); print(beta_hat) # check they match!!
  
  mu_hat <- X %*% beta_hat
  print("mu_hat: "); print(mu_hat)
  sigma_squared_hat <- (norm(y - mu_hat)^2)/(n-kappa)   ## norm(y-mu_hat)^2/(n-kappa)
  print("signma_squared_hat: "); print(sigma_squared_hat)
  GCV <- sigma_squared_hat/(n-kappa)                ## GCV = sigma^2 / (n-kappa),
  print("GCV : "); print(GCV)
  
  output <- list()  # Should contain elements: coef (??^), fitted (µ^) and sig2 (??^2), knots vector, anything else?
  class(output) <- "psline"
  return(output)
  } ## pspline

pspline(x,y,k=40,logsp=5,bord=3,pord=2,ngrid=100) # Test my function

## Error QT_y dimensions not matching up nicely... Think this is messing up my fast beta_hat
## Slow beta hat works. NB: can't set k higher than the number of data points I think


## Q2

print.pspline(m)


## Q3 

predict.pspline(m,x,se=TRUE)

## Q4
plot.pspline(m)




