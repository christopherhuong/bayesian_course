
library(rethinking)

data(Howell1)
d <- Howell1

d2 <- d[d$age >=18, ]

#prior of average height: N(178,20)
curve(dnorm(x, 178, 20), from=100,to=250)


#prior of sigma: Unif(0,50)
curve(dunif(x, 0, 50), from=-10, to=60)

#prior predictive distribution

sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)


# Grid approximation of the posterior distribution -------------------------------------------------------------------------

mu.list <- seq(from=150, to=160, length.out=100)
sigma.list <- seq(from=7, to=9, length.out=100)
post <- expand.grid(mu=mu.list, sigma=sigma.list)

post$LL <- sapply(1:nrow(post), function(i) sum(
  dnorm(d2$height, post$mu[i], post$sigma[i], log=T)))

post$prod <- post$LL + dnorm(post$mu, 178, 20, log=T) +
  dunif(post$sigma, 0, 50, log=T)

post$prob <- exp(post$prod - max(post$prod))

# inspect posterior distribution
contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)


# sample from posterior (2 parameters) in proportion to their prob ------------------------------------

sample.rows <- sample(1:nrow(post), size=1e4, replace=T,
                      prob=post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]

plot(sample.mu, sample.sigma, cex=0.5, pch=16,col=col.alpha(rangi2,0.1))
dens(sample.mu)
dens(sample.sigma)
#confidence / credible intervals from samples
PI(sample.mu)
PI(sample.sigma)


# quap to approximate posterior distribution --------------------------------------------------------------------

flist <- alist(
  height ~ dnorm(mu, sigma),#likelihood function used in bayes theorem
  mu ~ dnorm(178, 20), #prior mu 
  sigma ~ dunif(0,50) #prior mu
)

m4.1 <- quap(flist, data=d2)

#mean and quantile of posterior distribution for mu and sigma
precis(m4.1)



# sampling from the quap posterior distribution ------------------------------------------------------
# posterior w/ 2 parameter dimensions = multidimensional gaussian distribution
# means & variance-covariance matrix to describe multidimensional gaussian distribution

vcov(m4.1)
cov2cor(vcov(m4.1))
#mu and sigma are uncorrelated, which they should be because they're independent

#extract samples from the quap posterior
post <- extract.samples(m4.1, n=1e4)
post <- MASS::mvrnorm(n=1e4, mu=coef(m4.1), Sigma=vcov(m4.1))
head(post)

precis(post)
plot(post)





















