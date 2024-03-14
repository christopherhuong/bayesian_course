

library(rethinking)


#information = the measured decrease in uncertainty upon observing the outcome
#information entropy = the uncertainty contained in a probability distribution 
#is the average log-probability of an event

#information entropy
p <- c(.2, .8)
-sum(p*log(p))



#divergence = additional uncertainty induced by using one probability distribution
#to describe another probability distribution
p <- c(.70, .10, .20)
q <- c(.10, .10, .80)
# KL divergence
sum(p*(log(p)-log(q)))
sum(p*log(p/q))



sum(q*(log(q)-log(p)))


# brain size --------------------------------------------------------------


sppnames <- c( "afarensis","africanus","habilis","boisei",
              "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame( species=sppnames , brain=brainvolcc , mass=masskg )

d$mass_std <- (d$mass - mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain / max(d$brain)


m7.1 <- quap(
  alist(
    brain_std~ dnorm(mu, log_sigma),
    mu <- a + b*mass_std,
    a ~ dnorm(0.5, 1),
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0, 1)
  ),
  data=d)

precis(m7.1)

set.seed(1)
#log-pointwise predictive-density (bayesian version of the log-probability score)
lppd(m7.1, n=1e4)
#log probability score for each specific observation
#sum these for a total log-probabiltiy score for the model and data
#larger values = better average accuracy
set.seed(1)
logprob <- sim(m7.1, ll=T, n=1e4) #log-probability of each observation
n <- ncol(logprob)
ns <- nrow(logprob)

lse <- function(x) { #compute the log of the sum of exponentiated values
  xmax <- max(x)
  xsum <- sum(exp(x-xmax)) #take all log-probs for given observation
  xmax + log(xsum)  #exponentiates teach, sums them, then takes the log
}

f <- function(i) log_sum_exp(logprob[i,]) - log(ns)

lppd <-sapply(1:n, f) |> print()


#AIC = estimate of average out-of-sample deviance
# more complex models tend to overfit in direct proportion
# to number of parameters

#WAIC = approximates the out-of-sample deviance
#that converges to the cross-validation approximation in a large sample
#it tries to guess the out-of-sample KL divergence


data(cars)
m <- quap(
  alist(
    dist ~ dnorm(mu, sigma),
    mu <- a + b*speed,
    a ~ dnorm(0, 100),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ),data=cars
)
set.seed(94)
post <- extract.samples(m, n=1e3)

#get log-likelihood of each observation i at each sample s from posterior
n_samples <- 1000
logprob <- sapply(1:n_samples,
                  function(s){
                    mu <- post$a[s] + post$b[s]*cars$speed
                    dnorm(cars$dist, mu, post$sigma[s], log=T)
                  })
#now compute lppd, the bayesian deviance
#average sample in each row, take the log, sum the logs
n_cases <- nrow(cars)
lppd <- sapply(1:n_cases,
               function(i) log_sum_exp(logprob[i,])-log(n_samples))

#now the penalty term: variance across samples for each observation
pWAIC <- sapply(1:n_cases, 
                function(i) var(logprob[i,]))

#WAIC
-2*(sum(lppd)-sum(pWAIC))

#standard error of WAIC 
waic_vec <- -2*(lppd - pWAIC)
sqrt(n_cases * var(waic_vec))


WAIC(m7.1)


#outliers and WAIC
data("WaffleDivorce"); d<- WaffleDivorce
d$A <- standardize(d$MedianAgeMarriage)
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)

# D ~ A
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA * A,
    a ~ dnorm(0, 0.2),
    bA ~dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=d)

# D ~ M
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),data=d)


# D ~ M + A
m5.3 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )


compare(m5.1, m5.2, m5.3, func=PSIS)

PSIS_m5.3 <- PSIS(m5.3, pointwise=T)
WAIC_m5.3 <- WAIC(m5.3, pointwise=T)
plot(PSIS_m5.3$k, WAIC_m5.3$penalty, col=rangi2,lwd=2)


# students t (thicker tails, robust to outliers), v=2
m5.3t <- quap(
  alist(
    D ~ dstudent(2, mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    c(bM, bA) ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),data=d)

PSIS(m5.3t)















