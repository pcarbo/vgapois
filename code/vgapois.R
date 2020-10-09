# Fit a variational Gaussian approximation to the Poisson-normal model
# by maximizing the variational lower bound ("ELBO"). Under this
# model, the count data y are Poisson with log-rates b0 + x*b. The
# unknown b is assigned a normal prior with zero mean and variance
# s0. The intractable posterior for b is approximated by a normal with
# mean mu and variance s.
# 
# The optim L-BFGS-B solver is used to fit the mean and variance of
# the normal approximation. Note that the optimal solution should be
# unique in most cases; see Arridge et al (2018) for details.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
vgapois1 <- function (x, y, b0, s0, mu = 0, s = 1, factr = 1e5,
                      maxit = 100, ...) {
  f <- function (par)
    -compute_elbo_vgapois1(x,y,b0,s0,par[1],par[2])
  g <- function (par)
    -compute_elbo_grad_vgapois1(x,y,b0,s0,par[1],par[2])
  out <- optim(c(mu,s),f,g,method = "L-BFGS-B",lower = c(-Inf,1e-15),
               control = list(factr = factr,maxit = maxit,...))
  names(out$par) <- c("mu","s")
  return(out)
}

# Compute the log-likelihood of the counts, y, under the Poisson model.
# See function vgapois1 for details.
compute_loglik_pois <- function (x, y, b0, b) {
  r <- b0 + x*b
  return(sum(dpois(y,exp(r),log = TRUE)))
}

# Compute the log-posterior of the counts, y, under the Poisson model.
# See function vgapois1 for details.
compute_logp_pois <- function (x, y, b0, b, s0)
  compute_loglik_pois(x,y,b0,b) + dnorm(b,sd = sqrt(s0),log = TRUE)

# Compute the variational lower bound ("ELBO") for the variational
# approximation to the Poisson model. See function vgapois1 for
# a description of the input arguments.
compute_elbo_vgapois1 <- function (x, y, b0, s0, mu, s) {
  r <- b0 + x*mu
  return(sum(y*r) - sum(exp(r + x^2*s/2)) - sum(lfactorial(y)) +
         (1 - mu^2/s0 - s/s0 + log(s) - log(s0))/2)
}

# Compute the gradient of the ELBO with respect to the mean (mu) and
# variance (s) of the normal approximation. See function vgapois1 for
# details about the input arguments.
compute_elbo_grad_vgapois1 <- function (x, y, b0, s0, mu, s) {
  r <- b0 + x*mu
  u <- exp(r + s*x^2/2)
  return(c(sum(y*x) - sum(x*u) - mu/s0,
           (1/s - 1/s0 - sum(x^2*u))/2))
}
