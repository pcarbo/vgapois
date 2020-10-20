# Illustration of variational Gaussian approximation for Poisson-normal
# model with two unknowns.
library(mvtnorm)
source("../code/vgapois.R")
set.seed(1)

# Simulate data.
n <- 14
b <- c(1.3,1.5)
A <- matrix(rnorm(2*n,mean = -2),n,2)
X <- matrix(rnorm(2*n),n,2)
Y <- matrix(0,n,2)
R <- A + scalecols(X,b)
Y <- matrix(rpois(2*n,exp(R)),n,2)

# Compute Monte Carlo estimate of marginal likelihood.
S0   <- rbind(c(3,2.8),
              c(2.8,3))
ns   <- 1e5
B    <- rmvnorm(ns,sigma = S0)
logw <- rep(0,ns)
for (i in 1:ns)
  logw[i] <- compute_loglik_pois(X,Y,A,B[i,])
d    <- max(logw)
logZ <- log(mean(exp(logw - d))) + d

# Compute importance sampling estimates of mean and covariance.
w     <- exp(logw - d)
w     <- w/sum(w)
mu.mc <- drop(w %*% B)
S.mc  <- crossprod(sqrt(w)*B) - tcrossprod(mu.mc)

# Fit variational approximation.
fit <- vgapois(X,Y,A,S0)
mu  <- fit$mu
S   <- fit$S

cat(fit$message,"\n")
cat(sprintf("Monte Carlo estimate:    %0.12f\n",logZ))
cat(sprintf("Variational lower bound: %0.12f\n",-fit$value))
cat("Monte Carlo mean estimates:\n")
print(mu.mc)
print(S.mc)
cat("Variational covariance estimates:\n")
print(mu)
print(S)

