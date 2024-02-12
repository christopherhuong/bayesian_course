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
library(rethinking)

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
# sample of standard deviations from uniform distribution 0 to 50
prior_h <- rnorm(10000, sample_mu, sample_sigma)
# sample heights from the prior distribution

density(prior_h, adjust=1) %>% plot()
# check to see if looks reasonable


# grid approximation of the posterior distribution ------------------------

mu.list <- seq(from=150, to=160, length.out=100)
#list of plausible mean heights
sigma.list <- seq(from=7, to=9, length.out=100)
#list of plausible SDs
post <- expand.grid(mu=mu.list, sigma=sigma.list)
#all possible combinations of height mean + SD
post$LL <- sapply(1:nrow(post), function(i) sum(
  dnorm(d2$height, post$mu[i], post$sigma[i], log=T)))
#sums log-likelihood of each mean+SD combination (row),
# given each observed height

post$prod <- post$LL + dnorm(post$mu, 178, 20, T) +
  dunif(post$sigma, 0, 50, T)
# multiply the prior by the likelihood to get the product
# that is proportional to the posterior density
# (adding logs = multipling the raw densities)

post$prob <- exp(post$prod - max(post$prod))
# subtract maximum posterior probability gets you relative posterior probabilities

#inspect
contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)

plot(post$prob, type="p")

# sampling from the posterior ---------------------------------------------

# now that we have 2 parameters, we sampling combinations of values
# sample from rows of post, in proportion to values in prob

sample.rows <- sample(1:nrow(post), size=10000, replace=T,
                      prob=post$prob)

plot(density(sample.rows, adjust=.1))

sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]

#now we have 10000 samples from the posterior height data
plot(sample.mu, sample.sigma, cex=0.5, pch=16,
     col=col.alpha(rangi2,0.1) )


plot(density(sample.mu, adjust=1))
plot(density(sample.sigma, adjust=1))

PI(sample.mu)
PI(sample.sigma)


# sample size and normality of sigma's posterior --------------------------

d3 <- sample(d2$height, size=20)


mu.list <- seq( from=150, to=170 , length.out=200 )
sigma.list <- seq( from=4 , to=20 , length.out=200 )
post2 <- expand.grid( mu=mu.list , sigma=sigma.list )
post2$LL <- sapply( 1:nrow(post2) , function(i)
  sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] ,
              log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
  dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
                        prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
      col=col.alpha(rangi2,0.1) ,
      xlab="mu" , ylab="sigma" , pch=16 )

dens( sample2.sigma , norm.comp=TRUE )



# approximating the posterior with quadratic approximation ----------------
rm(list=ls())
data(Howell1)

d <- Howell1
d2 <- d[d$age >= 18, ]

flist <- alist(
  height ~ dnorm(mu, sigma),  #height is normally distributed with mean=mu, and sd=sigma
  mu ~ dnorm(178, 20), #mu is normally distributed, with mean=178, sd=20
  sigma ~ dunif(0, 50) #sigma is uniformally distributed from 0 to 50
)

m4.1 <- quap(flist, data=d2)
#look at posterior distribution
precis(m4.1)


start <- list(
  mu=mean(d2$height),
  sigma=sd(d2$height)
)

m4.1 <- quap(flist, data=d2, start=start)
precis(m4.1)

#smaller prior for SD
m4.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0,50)),
    data=d2
  )

precis(m4.2)


# sampling from quap ------------------------------------------------------

vcov(m4.1)

post <- extract.samples(m4.1, n=1e4)
head(post)

precis(post)

library(MASS)
post <- mvrnorm(n=1e4, mu=coef(m4.1), Sigma=vcov(m4.1))
head(post)


# regression --------------------------------------------------------------
rm(list=ls())

data(Howell1); d <- Howell1; d2 <- d[d$age >= 18, ]

plot(d2$height ~ d2$weight)


# height ~ N(ui, sigma)
# ui = alpha + beta * (xi - xbar)
# alpha ~ N(178, 20)
# beta ~ N(0, 10)
# sigma ~ N(0, 50)


# simulate heights from the priors ----------------------------------------
library(rethinking)
set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rnorm(N, 0, 10)



plot(NULL, xlim=range(d2$weight), ylim=c(-100,400),
     xlab="weight", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1,lwd=0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for (i in 1:N) curve( a[i] + b[i] * (x - xbar),
                      from=min(d2$weight), to=max(d2$weight), add=T,
                      col=col.alpha("black",0.2))


# new prior, beta ~ log-normal(0,1)
b <- rlnorm(1e4, 0, 1)
dens(b, xlim=c(0,5), adj=0.1)


set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rlnorm(N, 0, 1)


plot(NULL, xlim=range(d2$weight), ylim=c(-100,400),
     xlab="weight", ylab="height")
abline(h=0, lty=2)
abline(h=272, lty=1,lwd=0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)
for (i in 1:N) curve( a[i] + b[i] * (x - xbar),
                      from=min(d2$weight), to=max(d2$weight), add=T,
                      col=col.alpha("black",0.2))




# finding posterior distribution ------------------------------------------
rm(list=ls())

# height ~ dnorm(mu, sigma)
# mu = a + b * (weight - xbar)
# a ~ dnorm(178, 20)
# b ~ dlnorm(0, 1)
# sigma ~ udunif(0, 50)

data(Howell1); d <- Howell1; d2 <- d[d$age >= 18, ]

xbar <- mean(d2$weight)

m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * (weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0,50)
  ), data=d2
)

precis(m4.3)
round(vcov(m4.3), 3)
pairs(m4.3)



# plotting the posterior --------------------------------------------------

plot(height ~ weight, data =d2, col=rangi2)
post <- extract.samples(m4.3, n=1e4)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map *(x - xbar), add=T)


# plotting uncertainty from the posterior ---------------------------------
N <- 30
dN <- d2[1:N, ]
mN <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * (weight - mean(weight)),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data=dN
)

# extract 20 samples from the posterior
post <- extract.samples( mN, n=20)
#display raw data and sample size
plot(dN$weight, dN$height,
     xlim=range(d2$weight), ylim=range(d2$height),
     col=rangi2, xlab="weight", ylab="height")
mtext(concat("N = ", N))

#plot lines
for (i in 1:20)
  curve(post$a[i] + post$b[i] * (x-mean(dN$weight)),
        col=col.alpha("black", 0.3), add=T)


# plotting regression intervals and contours ------------------------------

post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b * (50 - xbar)

PI(mu_at_50, prob=0.89)


mu <- link(m4.3)
str(mu)

# definte sequence of weights to compute predictions for
# these values will be on horizontal axis
weight.seq <- seq(from=25, to=70, by=1)

#use link to compute mu for each sample from posterior and for each weight in weight.seq
mu <- link(m4.3, data=data.frame(weight=weight.seq))
str(mu)

plot(height ~ weight, d2, type="n")

for (i in 1:100) 
  points(weight.seq, mu[i], pch=16, col=col.alpha(rangi2, 0.1))



# curved lines ------------------------------------------------------------
rm(list=ls())
library(rethinking)
data(Howell1)
d <- Howell1

plot(height ~ weight, d)























