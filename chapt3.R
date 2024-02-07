

P_Positive_Vampire <- 0.95
P_Positive_Mortal <- 0.01
P_Vampire <- 0.001

P_Positive <- P_Positive_Vampire*P_Vampire + P_Positive_Mortal*(1-P_Vampire)



# -------------------------------------------------------------------------

#simulate posterior predictive distribution

posterior_samples <- rbeta(1000, 6+1, 3+1)
#each value is a random sample drawn from the 
#beta distribution based on our sample W=6, L=3
rethinking::dens(posterior_samples, lwd=4, col=2, xlab="prop. water", adj=0.1)
curve(dbeta(x, 6+1, 3+1), add=T, lty=2, lwd=3)
#the sample distribution approximates the actual posterior beta distribution
hist(posterior_samples, breaks=50)
plot(posterior_samples)

# -------------------------------------------------------------------------
#function to simulate a N sized sample of W's and L's based on probability p
library(ggplot2)
sim_globe <- function(p, N){
  sample(c("W","L"), size=N, prob=c(p,1-p), replace=T)
}

#simulate drawing N=1000 samples of p 
#from a posterior beta distribution  based on our sample of W=6, L=3
post_samples <- rbeta(1000, 6+1, 3+1)

#posterior distribution
ggplot(data.frame(values=post_samples), aes(x=values)) +
  geom_density() +
  xlim(0,1) + labs(x="probability of W")

#draw samples of 10 based on the 1000 sample of p's
pred_post <- sapply(post_samples, function(p) sum(sim_globe(p,10)=="W"))

#predictive distribution
ggplot(data.frame(values=pred_post), aes(x=values)) +
  geom_histogram(bins=100) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x="number of W samples")


# for (i in 0:10) lines(c(i,i),c(0, tab_post[i+1]), lwd=4, col=4)

#for each possible explanation of the data
#count all the ways the data can happen
#explanations with more ways to produce the data are more plausible


# -------------------------------------------------------------------------

# sampling from the grid-approximate posterior
p_grid <- seq(from=0, to=1, length.out=1000) #possible probabilities of W
prob_p <- rep(1, 1000) #uniform probability of each probability
prob_data <- dbinom(6, size=9, prob=p_grid) #how probable W=6, L=3 are given each possible probability
plot(prob_data)

posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)
plot(posterior)


#draw samples from possible probabilities (0.0-1.0) based on our posterior distribution
samples <- sample(p_grid, prob=posterior, size=10000, replace=T)
plot(samples)
ggplot(data.frame(values=samples), aes(x=values)) +
  geom_histogram(bins=1000) + xlim(0,1) + labs(x="samples from probability of W")

library(rethinking)
dens(samples)

# All you’ve done so far is crudely replicate the posterior density you had already
# computed. That isn’t of much value. But next it is time to use these samples to describe and
# understand the posterior. That is of great value

# -------------------------------------------------------------------------

#find posterior probability that the proportion of water is less than 0.5
#intervals defined by boundaries

#using grid-approximate posterior, add up all probabilities where parameter value < 0.5

sum( posterior[ p_grid < 0.5 ] ) # = 0.1718746
sum(posterior[1:500])

# So about 17% of the posterior probability is below 0.5.

# perform the same calculation, using samples from the posterior
sum( samples < 0.5 ) / 10000

#ask how much posterior probability lies between 0.5 and 0.75
sum( samples > 0.5 & samples < 0.75 ) / 10000


#intervals defined by probability mass

#find boundaries of lower 80% of posterior probability

quantile(samples, 0.8)
#the lower 80% of the posterior probability distribution is between 
# p = 0 and p = 0.7597598 

quantile(samples, 0.1, 0.9)
# the probabilities that lie at the boundaries of the middle 80%
# Intervals of this sort, which assign equal probability mass to each tail, are very common
# in the scientific literature. We’ll call them percentile intervals (PI).


# -------------------------------------------------------------------------
library(rethinking)
#flat prior, posterior computed from 3 W's
p_grid <- seq(from=0, to=1, length.out=1000)
prior <- rep(1, 1000)
likelihood <- dbinom(3, size=3, prob=p_grid)
plot(p_grid, likelihood)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
plot(p_grid, posterior)

samples <- sample(p_grid, size=10000, replace=T, prob=posterior)
plot(samples)
#a sample of probability parameters of W, from the posterior 

library(ggplot2)
ggplot(data.frame(values=samples), aes(x=values)) +
  geom_density(adjust=.1) +
  xlim(0,1) + labs(x="probability of W")

PI(samples, prob=0.5)
#      25%          75% 
#     0.7047047    0.9299299 
#central 50% probability boundaries. Does not describe the most probabable
#probability values (near)

#highest probability density interval 
# = narrowest interval containing the specified probability mass (50%)
HPDI(samples , prob=0.5 )

#       |0.5      0.5| 
#   0.8388388 1.0000000 


# -------------------------------------------------------------------------

#point estimates from posterior distribution

# the parametervalue with highest posterior probability
# maximum a posteriori estimate (MAP)
p_grid[ which.max(posterior) ] # 1
#finds the index for the max value of posterior vector, 
#finds corresponding p_grid value

chainmode( samples , adj=0.01 )
# Returns estimated mode of a density computed from samples




# simulating data implied from model ---------------------------------------------------------

dbinom(0:2, size=2, prob=0.7)
# 0.09   0.42   0.49
# This means that there’s a 9% chance of observing w = 0,
# a 42% chance of w = 1, and a 49% chance of w = 2 with 2 tosses

#simulate observations using these probabilities
# by sampling from the above binomial distribution

dummy_w <- rbinom(1e5, size=9, prob=0.7)
table(dummy_w)/1e5
hist(dummy_w)



# -------------------------------------------------------------------------
p_grid <- seq(from=0, to=1, length.out=1000)
prior <- rep(1, 1000)
likelihood <- dbinom(6, size=9, prob=p_grid)
plot(p_grid, likelihood)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
plot(p_grid, posterior)

samples <- sample(p_grid, size=10000, replace=T, prob=posterior)
plot(samples)
ggplot(data.frame(values=samples), aes(x=values)) +
  geom_density(adjust=.1) +
  xlim(0,1) + labs(x="probabilities of W sampled from posterior")


###########

w <- rbinom(1000, size=9, prob=0.6)
simplehist(w)


w <- rbinom(1000, size=9, prob=samples)
simplehist(w)

# For each sampled value, a random binomial observation 
# is generated. Since the sampled values appear in proportion 
# to their posterior probabilities, the resulting simulated 
# observations are averaged over the posterior











# class demo --------------------------------------------------------------

p_grid <- seq(from=0, to=1, length.out=1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(6, size=9, prob = p_grid)

posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)

plot(p_grid, posterior, type="l")
plot(posterior)


sample <- sample(p_grid, size=10000, prob=posterior, replace=T)

plot(sample, type="p", pch=".")


require(rethinking)

dens(sample)


ggplot(data.frame(values=sample), aes(x=values)) +
  geom_density(adjust=1) +
  xlim(0,1) + labs(x="probability of W")


#how do you get the function for which to integrate
#why not use function defined by your posterior distribution

sum(posterior * p_grid)
mean(sample)
sum(posterior[p_grid < 0.5])
sum(sample<0.5)/10000
mean(sample < 0.5)

mean(prob_data <0.5)
# prior model -> new data -> bayesian updating -> 
# posterior model (distribution) -> 
# simulate implied data -> (because true function is complex??)
# inferences about the world

p_grid[posterior == max(posterior)]

chainmode(sample)


binom <- rbinom(length(sample),size=9, prob=sample)
#size = how many tosses of globe

simplehist(binom)
















































