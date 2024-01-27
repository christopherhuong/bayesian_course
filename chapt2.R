

# -------------------------------------------------------------------------


sample <- c("W","L","W","W","W","L","W","L","W")
W <- sum(sample == "W") #observed number of W
L <- sum(sample == "L") #observed number of L
p <- c(0, 0.25, 0.5, 0.75, 1) #proportions of W. Prior probability


#how many ways the proportion(conjecture) can produce the observed sample
ways <- sapply(p, function(q) (q*4)^W * ((1-q)*4)^L)

prob <- ways / sum(ways) #probability of proportion as the true data generating model (posterior distribution/probability)
cbind(p, ways, prob)


# -------------------------------------------------------------------------
# Ways for p to produce W,L = (4p)^W * (4-4p)^L
# function to compute posterior distribution


compute_posterior <- function(sample, poss) {
  W <- 3 #number of W observed
  L <- 3  #number of L observed
  n <- length(poss)-1
  ways <- sapply(poss, function(q) (q*n)^W * ((1-q)*n)^L )
  post <- ways/sum(ways)
  print(data.frame(poss, ways, post=round(post,3)))
  plot(post)
}

p = seq(from=0, to=1, length.out=21)
compute_posterior(sample, p)


# -------------------------------------------------------------------------

dbinom(55, size=120, prob=.5)

# -------------------------------------------------------------------------

# Grid Approximation
#1) define the grid. 
#decide how many points to use in estimating the posterior,
#then make a list of the parameter values on the grid.
n=21
p_grid <- seq(from=0, to=1, length.out=n)
p_grid
# 2) definite prior at each parameter value on the grid
prior <- rep(1, n)
prior
# 3) compute likelihood at each parameter value
likelihood <- dbinom(6, size=9, prob=p_grid)
likelihood
# 4) compute the unstandardized posterior at each parameter value
# by multiplying the prior by the likelihood
unstd.posterior <- likelihood * prior
# 5) standardize the posterior by dividing each value by the sum of all values
posterior <- unstd.posterior / sum(unstd.posterior)
posterior
plot( p_grid , posterior , type="b" )

# try different priors
prior <- ifelse(p_grid >.7, 0, 1)
prior <- exp( -5*abs( p_grid - 0.5 ) )


# quadratic approximation --------------------------------------------------

# A Gaussian approximation is called “quadratic approximation” 
# because the logarithm of a Gaussian distribution forms a parabola. 

d <- rnorm(1000)
hist(d, breaks=50)
hist(log(d), breaks=50)

library(rethinking)

globe.qa <- quap(
  alist(
    W ~ dbinom(W + L, p), # binomial likelihood
    p ~ dunif(0,1) # uniform prior
  
  ),
  data = list(W = 6, L = 3)
)

precis(globe.qa)

#     mean   sd  5.5% 94.5%
#   p 0.67 0.16  0.42  0.92

# Assuming the posterior is Gaussian, it
# is maximized at 0.67, and its standard deviation is 0.16.


# markov chain monte carlo ------------------------------------------------

n_samples <- 1000
p <- rep(NA, n_samples)
p[1] <- 0.5
W <- 6
L <- 3

for (i in 2:n_samples) {
  p_new <- rnorm(1, p[i-1], 0.1)
  if(p_new < 0) p_new <- abs(p_new)
  if(p_new > 1) p_new <- 2-p_new
  q0 <- dbinom(W, W+L, p[i-1])
  q1 <- dbinom(W, W+L, p_new)
  p[i] <- ifelse( runif(1) < q1/q0, p_new, p[i-1])
}
# The values in p are samples from the posterior distribution.
# To compare to the analytical posterior:






























