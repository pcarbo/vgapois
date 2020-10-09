# TO DO: Explain here what this script is for, and how to use it.
library(pracma)
source("../code/vgapois.R")
set.seed(1)

# SIMULATE DATA
# -------------
# TO DO: Explain here what these lines of code do.
n <- 10
b0 <- -1
b <- 1.5
x <- rnorm(n)
r <- exp(b0 + x*b)
y <- rpois(n,r)

# PLOT ELBO SURFACE
# -----------------
# TO DO: Explain here what these lines of code do.
mu   <- seq(-1,2.5,length.out = 50)
s    <- 10^seq(-3,0,length.out = 50)
out  <- meshgrid(s,mu)
S    <- out$X
MU   <- out$Y
n    <- length(S)
elbo <- matrix(0,nrow(S),ncol(S))
for (i in 1:n)
  elbo[i] <- compute_elbo_vgapois1(x,y,b0,s0,MU[i],S[i])
contour(mu,s,elbo,col = "dodgerblue",nlevels = 24,
        xlab = "mu",ylab = "s")

# TO DO: Compute importance sampling estimate of marginal likelihood.

# FIT VARIATIONAL APPROXIMATION
# -----------------------------
# TO DO: Explain here what these lines of code do.
fit <- vgapois1(x,y,b0,s0)
mu  <- fit$par["mu"]
s   <- fit$par["s"]
print(compute_elbo_vgapois1(x,y,b0,s0,mu,s),digits = 6)
points(mu,s,pch = 4,col = "magenta")

# PLOT (EXACT) POSTERIOR SURFACE
# ------------------------------
# TO DO: Explain here what these lines of code do.
s0   <- 10
b    <- seq(-1,3,length.out = 1000)
n    <- length(b)
logp <- rep(0,n)
qb   <- dnorm(b,mu,sqrt(s))
for (i in 1:n)
  logp[i] <- compute_logp_pois(x,y,b0,b[i],s0)
maxlp <- max(logp)
plot(b,exp(logp - maxlp),type = "l",lwd = 2,col = "dodgerblue",
     xlab = "b",ylab = "relative posterior")
lines(b,qb/max(qb),col = "magenta",lwd = 2)
