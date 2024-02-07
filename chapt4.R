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




