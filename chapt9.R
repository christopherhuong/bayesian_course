
library(rethinking)



# markov island example ---------------------------------------------------
num_weeks <- 1e5
positions <- rep(0, num_weeks)

current <- 10
for(i in 1:num_weeks) {
  positions[i] <- current #record current position
  proposal <- current + sample(c(-1, 1), size=1) #flip coin
  if(proposal < 1) proposal <- 10
  if(proposal > 10) proposal <- 1
  prob_move <- proposal/current
  current <- ifelse(runif(n=1,min=0,max=1) < prob_move, proposal, current)
}


hist(positions)
plot(table(positions))



# concentration of measure in high dimensionality -------------------------
rm(list=ls())

D <- 10
T <- 1e3
# multivariate gaussian sample of 10 variables
# means = 0, sd= 10x10 variance-covariance matrix
# variances = 1, covariances = 0
# thus, 10 independent variables
Y <- rmvnorm(T, rep(0, D), diag(D))
cor(Y) |> round(2) # independent vars

rad_dist <- function(Y) sqrt(sum(Y^2)) #square each entry and sum
Rd <- sapply(1:T, function(i) rad_dist(Y[i, ]))

dens(Rd)
# further from zero with more dimensions


# Hamiltonian Monte Carlo -------------------------------------------------------------------------

# 1) function U to return negative log-probability of the data at current position (parameter values)
U <- function(q, a=0, b=1, k=0, d=1){
  mu_y <- q[1]
  mu_x <- q[2]
  U <- sum( dnorm(y, mu_y, 1, log=T)) + sum(x, mu_x, 1, log=T) +
    dnorm(mu_y, a, b, log=T) + dnorm(mu_x, k, d, log=T)
  return(-U)
}


# 2) gradient function 
# need vector of partial derivatives of U with respect to vector q
U_gradient <- function(q, a=0, b=1, k=0, d=1){
  mu_y <- q[1]
  mu_x <- q[2]
  G1 <- sum(y-mu_y) + (a-mu_y)/b^2 # dU/mu_y
  G2 <- sum(y-mu_x) + (k-mu_x)/d^2 # dU/mu_x
  return(c(-G1, -G2)) # negative because energy is neg-log-prob
}

# HMC2 function- 

HMC2_fx <- function(U, grad_U, epsilon, L, current_q){
  q = current_q
  p = rnorm(length(q), 0.1) #random flick - p is momentum
  current_p = p
  # make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # initialize bookkeeping - saves trajectory
  qtraj <- matrix(NA, nrow=L+1, ncol=length(q))
  ptraj <- qtraj
  qtraj[1, ] <- current_q
  ptraj[1, ] <- p
  #alternate full steps for position and momentum
  for(i in 1:L) {
    q = q + epsilon * p #full step for the position
    #make a full step for the momentum, except at end of trajectory
    if(i != L) {
      p = p - epsilon * grad_U(q)
      ptraj[i + 1, ] <- p
    }
    qtraj[i+1, ] <- q
  }
  #make a half step for momentum at the end
  p = p - epsilon * grad_U(q) / 2
  ptraj[L + 1, ] <- p
  # negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # accept or reject the state at end of trajectory, returning
  # either the position at end of trajectory or the initial position
  accept <- 0
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K)){
    new_q <- q #accept
    accept <- 1
  } else new_q <- current_q #reject
  return(list(q = new_q, traj = qtraj, ptraj = ptraj, accept = accept))
}



# test data
set.seed(7)
y <- rnorm(50)
x <- rnorm(50)
y <- as.numeric(scale(y))
x <- as.numeric(scale(x))

library(shape)
Q <- list()
Q$q <- c(-0.1, 0.2)
pr <- 0.3

plot(NULL, ylab="mu_y", xlab="mu_x", xlim=c(-pr,pr), ylim=c(-pr,pr))
step <- 0.03
L <- 11
n_samples <- 4
path_col <- col.alpha("black", 0.5)
points(Q$q[1], Q$q[2], pch=4, col="black")
for(i in 1:n_samples){
  Q <- HMC2_fx(U, U_gradient, step, L, Q$q)
  if(n_samples<10) {
    for(j in 1:L){
      K0 <- sum(Q$ptraj[j,]^2)/2 #kinetic energy
      lines(Q$traj[j:(j+1), 1], Q$traj[j:(j:j+1), 2], col=path_col, lwd=1+2*K0)
    }
    points(Q$traj[1:L+1, ], pch=16, col="white", cex=0.35)
    Arrows(Q$traj[L, 1], Q$traj[L,2], Q$traj[L+1, 1], Q$traj[L+1, 2],
           arr.length=0.35, arr.adj=0.7)
    text(Q$traj[L+1, 1], Q$traj[L+1, 2], i, cex=0.8, pos=4, offset=0.4)
    }
    points(Q$traj[L+1, 1], Q$traj[L+1, 2], pch=ifelse(Q$accept==1, 16, 1),
           col=ifelse(abs(Q$dH) > 0.1, "red", "black"))
}


# -------------------------------------------------------------------------
library(rethinking)

data(rugged); d <- rugged; rm(rugged)

d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000), ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp) #mean = 1
dd$rugged_std <- dd$rugged / max(dd$rugged) #max = 1
dd$cid <- ifelse(dd$cont_africa==1, 1, 2)


m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std - mean(rugged_std)),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ),data=dd
)
precis(m8.3,depth=2)

dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)
str(dat_slim)


m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid]*(rugged_std - 0.215),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ),data=dat_slim, chains=4, cores=4
)

show(m9.1)
precis(m9.1, depth=2)
pairs(m9.1)

stancode(m9.1)

traceplot(m9.1, chains=4)
trankplot(m9.1)




# taming a wild chain -----------------------------------------------------

y <- c(-1, 1)
set.seed(11)
m9.2 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a,
    a ~ dnorm(0, 1000),
    sigma ~ dexp(0.0001)
  ),data=list(y=y), chains=3
)

precis(m9.2)

pairs(m9.2)
pairs(m9.2@stanfit)
traceplot(m9.2)
trankplot(m9.2)

# add weakly informative priors
set.seed(11)
m9.3 <- ulam(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- alpha,
    alpha ~ dnorm(1, 10),
    sigma ~ dexp(1)
  ),data=list(y=y), chains=3
)

precis(m9.3)
pairs(m9.3)
traceplot(m9.3)
trankplot(m9.3)



# non-identifiable parameters ---------------------------------------------
rm(list=ls())
gc()
set.seed(41)
y <- rnorm(100)

# fit whacky model
set.seed(384)
m9.4 <- ulam(alist(
  y ~ dnorm(mu, sigma),
  mu <- a1 + a2,
  a1 ~ dnorm(0, 1000),
  a2 ~ dnorm(0, 1000),
  sigma ~ dexp(1)
),data=list(y=y), chains=3)


precis(m9.4)

traceplot(m9.4)
trankplot(m9.4)


m9.5 <- ulam(alist(
  y ~ dnorm(mu, sigma),
  mu <- a1 + a2,
  a1 ~ dnorm(0, 10),
  a2 ~ dnorm(0, 10),
  sigma ~ dexp(1)
), data=list(y=y), chains=3)

precis(m9.5)
traceplot(m9.5)
trankplot(m9.5)



# wine example ------------------------------------------------------------

rm(list=ls())
data("Wines2012"); d <- Wines2012; rm(Wines2012)

d <- list(
  S = scale(d$score),
  J = as.numeric(d$judge),
  W = as.numeric(d$wine),
  WO = ifelse(d$wine.amer==1, 1, 2),
  JO = ifelse(d$judge.amer==1, 1, 2)
)


mQ <- ulam(
  alist(
    S ~ dnorm(mu, sigma),
    mu <- Q[W], #estimating quality of wine, using wine index variable
    Q[W] ~ dnorm(0, 1),
    sigma~ dexp(1)
  ),data=d, chains=4, cores=4
)

precis(mQ, 2)
traceplot(mQ)
trankplot(mQ)


# Score ~ wine quality + wine origin
mQO <- ulam(
  alist(
    S ~ dnorm(mu, sigma),
    mu <- bQ[W] + bWO[WO],
    bQ[W] ~ dnorm(0, 1),
    bWO[WO] ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),data=d, chains=4, cores=4
)
plot(precis(mQO,2))

# item response model (introduce judge effects)
MQOJ <- ulam(
  alist(
    S ~ dnorm(mu, sigma),
    #subtract judge harshness and multiply by judge discrimination
    mu <- (bQ[W] + bWO[WO] - H[J]) * D[J],
    bQ[W] ~ dnorm(0, 1),
    bWO[WO] ~ dnorm(0, 1),
    H[J] ~ dnorm(0, 1),
    D[J] ~ dexp(1),
    sigma ~ dexp(1)
  ),data=d, chains=4, cores=4
)

plot(precis(MQOJ, 2))


ma <- quap(
  alist(
    S ~ dnorm(mu, sigma),
    mu <- -H[J],
    H[J] ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),data=d
)

plot(precis(ma, 2))














