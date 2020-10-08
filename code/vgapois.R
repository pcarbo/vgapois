# TO DO: Explain here what this function does, and how to use it.
#
# This implementation is meant to be "instructive"---that is, I've
# tried to make the code as simple as possible, with an emphasis on
# clarity. Very little effort has been devoted to making the
# implementation efficient, or the code concise.
vgapois1d <- function (x, y) {

}

# TO DO: Explain here what this function does, and how to use it.
compute_logp_pois <- function (x, y, b, s0)
  sum(dpois(y,exp(x*b),log = TRUE)) + dnorm(b,sd = s0^2,log = TRUE)

# TO DO: Explain here what this function does, and how to use it.
compute_elbo_vgapois1d <- function (x, y, s0, mu, s)
  mu*sum(x*y) - sum(exp(x*mu + x^2*s/2)) - sum(lfactorial(y)) +
    (1 - mu^2/s0 - s/s0 + log(s) - log(s0))/2
