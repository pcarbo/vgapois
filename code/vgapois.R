# Fit a variational Gaussian approximation to the univariate
# Poisson-normal model by maximizing the variational lower bound
# ("ELBO"). Under this model, the count data y are Poisson with
# log-rates a + x*b. The unknown b is assigned a normal prior with
# zero mean and variance s0. The intractable posterior for b is
# approximated by a normal with mean mu and variance s.
# 
# The optim L-BFGS-B solver is used to fit the mean and variance of
# the normal approximation. Note that the optimal solution should be
# unique in most cases; see Arridge et al (2018) for details.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
vgapois1 <- function (x, y, a, s0, mu = 0, s = 1, factr = 1e5,
                      maxit = 100, ...) {
  f <- function (par)
    -compute_elbo_vgapois1(x,y,a,s0,par[1],par[2])
  g <- function (par)
    -compute_elbo_grad_vgapois1(x,y,a,s0,par[1],par[2])
  out <- optim(c(mu,s),f,g,method = "L-BFGS-B",lower = c(-Inf,1e-15),
               control = list(factr = factr,maxit = maxit,...))
  names(out$par) <- c("mu","s")
  return(out)
}

# TO DO: Explain here what this function does, and how to use it.
vgapois <- function (X, Y, A, S0, mu = rep(0,ncol(X)), S = diag(ncol(X)),
                     factr = 1e5, maxit = 100, ...) {
  f <- function (mu)
    -compute_elbo_vgapois(X,Y,A,S0,mu,S)
  g <- function (mu)
    -compute_elbo_grad_vgapois(X,Y,A,S0,mu,S)
  out <- optim(mu,f,g,method = "L-BFGS-B",
               control = list(factr = factr,maxit = maxit,...))
  out$mu <- out$par
  return(out)
}

# TO DO: Explain here what this function does, and how to use it,
get_vgapois_params <- function (par, n) {

}

# Compute the log-likelihood of the counts, y, under the univariate
# Poisson-normal model.  See function vgapois1 for details.
compute_loglik_pois1 <- function (x, y, a, b) {
  r <- a + x*b # Poisson log-rates
  return(sum(dpois(y,exp(r),log = TRUE)))
}

# Compute the log-likelihood of the counts, Y, under the multivariate
# Poisson-normal model. See function vgapois for details. Note that
# the special case of one dimension (i.e., the univariate model) is
# also permitted, and should give the same result as compute_loglik_pois1.
compute_loglik_pois <- function (X, Y, A, b) {
  R <- A + scalecols(X,b) # Poisson log-rates
  return(sum(dpois(Y,exp(R),log = TRUE)))
}

# Compute the log-posterior of the counts, y, under the univariate
# Poisson-normal model. See function vgapois1 for details.
compute_logp_pois1 <- function (x, y, a, b, s0)
  compute_loglik_pois1(x,y,a,b) + dnorm(b,sd = sqrt(s0),log = TRUE)

# Compute the variational lower bound ("ELBO") for the variational
# approximation to the univariate Poisson-normal model. See function
# vgapois1 for a description of the input arguments.
compute_elbo_vgapois1 <- function (x, y, a, s0, mu, s) {
  r <- a + x*mu    # mean log-rates
  u <- r + s*x^2/2 # "overdispersed" log-rates
  return(sum(dpois(y,exp(u),log = TRUE) - y*s*x^2/2) - kl_norm1(0,mu,s0,s))
}

# Compute the variational lower bound ("ELBO") for the variational
# approximation to the multivariate Poisson-normal model. See function
# vgapois for details on the input arguments.
compute_elbo_vgapois <- function (X, Y, A, S0, mu, S) {
  n <- nrow(X)
  m <- ncol(X)
  f <- (m - drop(t(mu) %*% solve(S0,mu))
        - sum(diag(solve(S0) %*% S))
        + determinant(S,logarithm = TRUE)$modulus
        - determinant(S0,logarithm = TRUE)$modulus)/2  
  for (i in 1:n) {
    x <- X[i,]
    y <- Y[i,]  
    a <- A[i,]
    r <- a + x*mu  # mean log-rates
    u <- r + diag(diag(x) %*% S %*% diag(x))/2 # "overdispersed" log-rates
    f <- f + sum(y*r) - sum(exp(u)) - sum(lfactorial(y))
  }
  return(f)
}

# Compute the gradient of the univariate Poisson-normal ELBO with
# respect to the mean (mu) and variance (s) of the normal
# approximation. See function vgapois1 for details about the input
# arguments.
compute_elbo_grad_vgapois1 <- function (x, y, a, s0, mu, s) {
  r <- a + x*mu    # mean log-rates
  u <- r + s*x^2/2 # "overdispersed" log-rates
  return(c(sum(x*(y - exp(u))) - mu/s0,
           (1/s - 1/s0 - sum(x^2*exp(u)))/2))
}

# TO DO: Explain here what this function does, and how to use it.
compute_elbo_grad_vgapois <- function (X, Y, A, S0, mu, S) {
  n <- nrow(X)
  m <- ncol(X)
  gmu <- -solve(S0,mu)
  for (i in 1:n) {
    x <- X[i,]
    y <- Y[i,]  
    a <- A[i,]
    r <- a + x*mu  # mean log-rates
    u <- r + diag(diag(x) %*% S %*% diag(x))/2 # "overdispersed" log-rates
    gmu <- g + x*(y - exp(u))
  }
  return(gmu)
}

# Return the Kullback-Leibler divergence KL(p1 || p2) between two
# (univariate) normal distributions p1 and p2, where p1 is normal with
# mean mu1 and variance s1, and p2 is normal with mean mu2 and
# variance s2.
kl_norm1 <- function (mu1, mu2, s1, s2)
  (log(s1/s2) + s2/s1 - 1  + (mu1 - mu2)^2/s1)/2

# TO DO: Explain here what this function does, and how to use it.
kl_norm <- function (mu1, mu2, S1, S2) {
  # TO DO.
}

# scalecols(A,b) scales each column A[,i] by b[i].
scalecols <- function (A, b)
  t(t(A) * b)
