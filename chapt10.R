
library(rethinking)



# 10 pebbles, 5 buckets ---------------------------------------------------

# 5 difference distributions of pebbles
p <- list(
  A = c(0, 0, 10, 0, 0),
  B = c(0, 1, 8,  1, 0),
  C = c(0, 2, 6,  2, 0),
  D = c(1, 2, 4,  2, 1),
  E = c(2, 2, 2,  2, 2)
)

# normalize such that it's a probability distribution
p_norm <- lapply(p, function(q) q/sum(q))
# check
lapply(p_norm, function(x) sum(x))
# compute information entropy of each (average log-probability)
(H <- sapply(p_norm, function(q) -sum(ifelse(q==0, 0, q*log(q)))))
# -sum(q*log(q))
H
# the pebble distribution with most number of ways = highest entropy
ways <- c(1, 90, 1260, 37800, 113400)
log(ways) / 10



# binomial distribution ---------------------------------------------------
rm(list=ls())

# p(B) = 0.5
# 4 different candidate probability distributions for 2 flips
p <- list()
          # WW   BW   WB   BB
p[[1]] <- c(1/4, 1/4, 1/4, 1/4)
p[[2]] <- c(2/6, 1/6, 1/6, 2/6)
p[[3]] <- c(1/6, 2/6, 2/6, 1/6)
p[[4]] <- c(1/8, 4/8, 2/8, 1/8)

# compute expected value of each (i dont udnerstand this)
sapply(p, function(p) sum(p*c(0, 1, 1, 2))) 
                              #how many B's
# WW   BW    WB   BB
c(1/4, 1/4, 1/4, 1/4) * 
 c(0,   1,   1,   2)



# compute entropy
sapply(p, function(p) -sum(p*log(p)))
# the binomial distribution A (flattest) has the largest entropy

# what if p=0.7?
p <- 0.7

A <- c(
  (1-p)^2,  p*(1-p), (1-p)*p, p^2
) |> print()

# not a flat distribution, yet still has max entropy
-sum(A*log(A))


sim.p <- function(E = 1.4){ #expected value of BLUE for 2 draws
  x123 <- runif(3, min=0, max=1)
  x4 <- ((E) * sum(x123) - x123[2] - x123[3]) / (2-E)
  z <- sum(c(x123, x4))
  p <- c(x123, x4) / z
  list(H = -sum(p*log(p)), p=p)
}

H <- replicate(1e4, sim.p())
dens(as.numeric(H[1,]), adj=0.1)

entropies <- as.numeric(H[1, ])
distributions <- H[2, ]
#largest observed entropy
max(entropies)
# same as what we calculated from the binomial distribution
# p^y * (1-p)^(n-y)

distributions[which.max(entropies)]

#the most even distribution that satisfies the p=0.7 constraint
#nominated by the binomial distribution, has the greatest entropy


















