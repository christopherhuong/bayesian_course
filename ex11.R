library(rethinking)

globe.qa <- quap(
  alist(
    W ~ dbinom(W+L, p), #binomial likelihood
    p ~ dunif(0,1)   #uniform prior
    
  ),
  data = list(W=6,L=3)
)

precis(globe.qa)


prior <- rep(1,3)

a <- dbinom(3, 3,)














