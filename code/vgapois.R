# Fit a variational Gaussian approximation to the univariate
# Poisson-normal model by maximizing the variational lower bound
# ("ELBO"). Under this model, the count data y are Poisson with
# log-rates a + x*b. The unknown b is assigned a normal prior with
# zero mean and variance s0. The intractable posterior for b is
# approximated by a normal with mean mu and variance s.
# 
# The optim L-BFGS-B solver is used to fit the mean and variance of
# the normal approximation. Note that the optimal solution should be
# unique; see Arridge et al (2018) for details.
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

# Fit a variational Gaussian approximation to the multivariate
# Poisson-normal model by maximizing the variational lower bound
# ("ELBO"). Under this model, the count data Y[i,j] are Poisson with
# log-rates A[i,j] + X[i,j]*b[j]. The unknown vector b is assigned a
# multivariate normal prior with zero mean and covariance S0. The
# intractable posterior for b is approximated by a multivariate normal
# with mean mu and covariance S.
#
# Note that the limiting univariate case (b is a scalar) is also
# handled.
#
# The optim L-BFGS-B solver is used to fit the mean and covariance of
# the normal approximation. Note that the optimal solution should be
# unique; see Arridge et al (2018) for details.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
vgapois <- function (X, Y, A, S0, mu = rep(0,ncol(as.matrix(X))),
                     S = diag(ncol(as.matrix(X))),
                     factr = 1e5, maxit = 100, ...) {
  f <- function (par) {
    params <- par2vgapois(par)
    return(-compute_elbo_vgapois(X,Y,A,S0,params$mu,params$R))
  }
  g <- function (par) {
    params <- par2vgapois(par)
    ans <- compute_elbo_grad_vgapois(X,Y,A,S0,params$mu,params$R)
    return(-vgapois2par(ans$mu,ans$R))
  }
  out <- optim(vgapois2par(mu,chol(S)),f,g,method = "L-BFGS-B",
               control = list(factr = factr,maxit = maxit,...))
  params <- par2vgapois(out$par)
  out$mu <- params$mu
  out$S  <- drop(crossprod(params$R))
  return(out)
}

# Given a mean (mu) and R = chol(S), where S is the covariance matrix,
# return a parameter vector passed to optim. It is a vector containing
# mu and the upper triangular portion of R = chol(S). Note that the
# limiting case of one dimension is also handled, in which case the
# return value is a 2-element vector.
vgapois2par <- function (mu, R)
  c(mu,R[upper.tri(R,diag = TRUE)])

# Given an optim "par" vector, output the mean (mu) and Cholesky
# factor (R) of the covariance matrix, such that S = crossprod(R).
# Note that the limiting case of one dimension (n = 1) is also
# handled.
par2vgapois <- function (par) {
  n   <- (sqrt(1 + 8*length(par)) - 1)/2
  mu  <- par[1:n]
  par <- par[-(1:n)]
  R   <- matrix(0,n,n)
  R[upper.tri(R,diag = TRUE)] <- par
  return(list(mu = mu,R = drop(R)))
}

# Compute the log-likelihood of the counts, y, under the univariate
# Poisson-normal model. See function vgapois1 for details.
compute_loglik_pois1 <- function (x, y, a, b) {
  r <- a + x*b # Poisson log-rates
  return(sum(dpois(y,exp(r),log = TRUE)))
}

# Compute the log-likelihood of the counts, Y, under the multivariate
# Poisson-normal model. See function vgapois for details. Note that
# the special case of one dimension is also handled, and should give
# the same result as compute_loglik_pois1.
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
  return(sum(dpois(y,exp(u),log = TRUE) - y*(u - r)) - klnorm1(0,s0,mu,s))
}

# Compute the variational lower bound ("ELBO") for the variational
# approximation to the multivariate Poisson-normal model, in which the
# approximating distribution is parameterized by its mean mu and the
# Cholesky factor of the covariance matrix, R = chol(S). See function
# vgapois for additional details on the input arguments. Note that the
# special case of one dimension is also handled, and should give the
# same result as compute_elbo_vgapois1.
compute_elbo_vgapois <- function (X, Y, A, S0, mu, R) {
  mu0 <- rep(0,length(mu))
  S <- crossprod(R)
  R <- A + scalecols(X,mu)          # mean log-rates
  U <- R + scalecols(X^2,diag(S))/2 # "overdispersed" log-rates
  return(sum(dpois(Y,exp(U),log = TRUE) - Y*(U - R)) - klnorm(mu0,S0,mu,S))
}

# Compute the gradient of the univariate Poisson-normal ELBO with
# respect to the mean (mu) and variance (s) of the normal
# approximation. See function vgapois1 for details about the input
# arguments.
compute_elbo_grad_vgapois1 <- function (x, y, a, s0, mu, s) {
  r   <- a + x*mu    # mean log-rates
  u   <- r + s*x^2/2 # "overdispersed" log-rates
  gmu <- sum(x*(y - exp(u))) - mu/s0
  gs  <- (1/s - 1/s0 - sum(x^2*exp(u)))/2
  return(c(gmu,gs))
}

# Compute the gradient of the multivariate Poisson-normal ELBO with
# respect to the mean (mu) and Cholesky factor of the covariance
# matrix, R = chol(S). See function vgapois for additional details
# about the input arguments.
compute_elbo_grad_vgapois <- function (X, Y, A, S0, mu, R) {
  n <- nrow(X)
  m <- ncol(X)
  S <- crossprod(R)
  gmu <- -solve(S0,mu)
  gR <- solve(t(R)) - R %*% solve(S0)
  for (i in 1:n) {
    x <- X[i,]
    y <- Y[i,]  
    a <- A[i,]
    D <- diag(x)
    r <- a + x*mu                  # mean log-rates
    u <- r + diag(D %*% S %*% D)/2 # "overdispersed" log-rates
    gmu <- gmu + x*(y - exp(u))
    gR <- gR - R %*% diag(x^2*exp(u))
  }
  return(list(mu = gmu,R = gR))
}

# Return the Kullback-Leibler divergence KL(p2 || p1) between two
# (univariate) normal distributions, p1 and p2, where p1 is normal
# with mean mu1 and variance s1, and p2 is normal with mean mu2 and
# variance s2.
klnorm1 <- function (mu1, s1, mu2, s2)
  (log(s1/s2) + s2/s1 + (mu1 - mu2)^2/s1 - 1)/2

# Return the Kullback-Leibler divergence KL(p2 || p1) between two
# multivariate normal distributions, p1 and p2, where p1 is
# multivariate normal with mean mu1 and covariance S1, and p2 is
# multivariate normal with mean mu2 and covariance S2.
klnorm <- function (mu1, S1, mu2, S2)
  drop((logdet(S1) - logdet(S2) + sum(diag(S2 %*% solve(S1)))
        + (mu1 - mu2) %*% solve(S1,mu1 - mu2) - length(mu1))/2)

# Compute the logarithm of the matrix determinant.
logdet <- function (x)
  as.numeric(determinant(as.matrix(x),logarithm = TRUE)$modulus)

# scalecols(A,b) scales each column A[,i] by b[i].
scalecols <- function (A, b)
  t(t(A) * b)
