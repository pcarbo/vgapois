# TO DO: Explain here what this script is for, and how to use it.
set.seed(1)

# SIMULATE DATA
# -------------
# TO DO: Explain here what these lines of code do.
n <- 20
b <- 1.5
x <- rnorm(n)
r <- exp(x*b)
y <- rpois(n,r)

# PLOT LIKELIHOOD SURFACE
# -----------------------
# TO DO: Explain here what these lines of code do.
b <- seq(-2,4,length.out = 1000)
n <- length(b)
loglik <- rep(0,n)
for (i in 1:n) {
  r <- exp(x*b[i])
  loglik[i] <- sum(dpois(y,r,log = TRUE))
}
plot(b,exp(loglik - max(loglik)),type = "l",lwd = 2,col = "black",
     xlab = "b",ylab = "relative likelihood")
