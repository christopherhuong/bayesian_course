library(ggplot2)


pos <- runif(100000, -1,1)
pos <- replicate(1000, sum(runif(16, -1,1)))


hist(pos, breaks=500)
plot(density(pos, adjust=.1))


#normal distribution emerging from summing values of 
#an underlying uniform distribution


prod(runif(12, 1, 1.1))

growth <- replicate(10000, prod(runif(12, 1, 1.1)))

hist(growth, breaks=500)
plot(density(growth, adjust=0.1))

#small multipicative effects of fluctuations 
#converge to gaussian

big <- replicate(10000, prod(runif(12, 1, 1.5)))
small <- replicate(10000, prod(runif(12, 1.0, 1.1)))

hist(big, breaks=500)
hist(small, breaks=500)


# large multipicative effects produce log-normal distributions
log.big <- replicate(1e4, log(prod(runif(12, 1, 2))))
hist(log.big, breaks=500)


x <- seq(from=-10, to=10, length.out=21)

y <- x^2

plot(y~x)


# -------------------------------------------------------------------------

# Approach for describing statistical models

# 1) Variables: 
#     observable = data
#     unobservable = parameters

# 2) Define each variable in terms of others variables, or in terms of a probability distribution

# 3) Combination of variables and their probability distributions = joint generative model
      # can simulate hypothetical observations and analyze real ones



# -------------------------------------------------------------------------
library(dplyr)

data("Howell1")
d <- Howell1

precis(d)

d2 <- d[d$age >= 18, ]

density(d2$height, adjust=1) %>% plot()

curve(dnorm(x, 178, 20), from=100, to=250)

curve(dunif(x, 0, 50), from=0, to=50)


# prior predictive distribution

sample_mu <- rnorm(10000, 178, 20) 
# sample of means from norm distribution centered at m=178, sd=20
sample_sigma <- runif(10000, 0, 50)
# sample of standard deviations from uniform distribution
prior_h <- rnorm(10000, sample_mu, sample_sigma)
density(prior_h, adjust=1) %>% plot()


















