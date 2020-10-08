# TO DO: Explain here what this script is for, and how to use it.
library(pracma)
source("../code/vgapois.R")
set.seed(1)

# SIMULATE DATA
# -------------
# TO DO: Explain here what these lines of code do.
n <- 20
b <- 1.5
x <- rnorm(n)
r <- exp(x*b)
y <- rpois(n,r)

# PLOT POSTERIOR SURFACE
# ----------------------
# TO DO: Explain here what these lines of code do.
s0 <- 10
b <- seq(-1,3,length.out = 1000)
n <- length(b)
loglik <- rep(0,n)
for (i in 1:n)
  loglik[i] <- compute_logp_pois(x,y,b[i],s0)
plot(b,exp(loglik - max(loglik)),type = "l",lwd = 2,col = "dodgerblue",
     xlab = "b",ylab = "relative likelihood")
