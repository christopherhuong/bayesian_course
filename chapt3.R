

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

#using grid-approximate posterior, add up all probabilities where parameter value < 0.5

sum( posterior[ p_grid < 0.5 ] ) # = 0.1718746
sum(posterior[1:500])

# So about 17% of the posterior probability is below 0.5.

# perform the same calculation, using samples from the posterior
sum( samples < 0.5 ) / 10000

#ask how much posterior probability lies between 0.5 and 0.75
sum( samples > 0.5 & samples < 0.75 ) / 10000


















