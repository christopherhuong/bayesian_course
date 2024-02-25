
library(rethinking)


data("WaffleDivorce"); d <- WaffleDivorce

#standardize vars; divorce rate, marriage rate, median age of marriage
d$D <- (d$Divorce - mean(d$Divorce)) / sd(d$Divorce)
d$M <- (d$Marriage - mean(d$Marriage)) / sd(d$Marriage)
d$A <- (d$MedianAgeMarriage - mean(d$MedianAgeMarriage)) / sd(d$MedianAgeMarriage)


#quadratic approximation of the posterior distribution
#of divorce as a linear function of age
m5.1 <- quap(alist(
  D ~ dnorm(mu, sigma),
  mu <- a + bA * A,
  a ~ dnorm(0, 0.2),
  bA ~ dnorm(0, 0.5),
  sigma ~ dexp(1)
), data = d
)

#of divorce as a linear function of marriage rate
m5.2 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

#simulate from prior distribution
prior <- extract.prior(m5.1) #10000 samples from the prior distribution
mu <- link(m5.1, post=prior, data=list(A=c(-2,2)))
plot(NULL, xlim=c(-3,3),ylim=c(-3,3))
for(i in 1:50) lines(c(-2,2), mu[i, ], col=col.alpha("black",0.4) )
#plausible regression lines implied by the prior distribution


#posterior predictions
A_seq <- seq(from=-3, to=3.2, length.out=30) #values of A 
#vector of possible standardized values of median age of marriage
mu <- link(m5.1, data=list(A=A_seq))
str(mu) #1000 samples of each unique value of A
mu.mean <- apply(mu, 2, mean)
#mean of 1000 samples for each unique value of A
mu.PI <- apply(mu, 2, PI, prob=0.89)

precis(m5.1)

#plot
plot(D ~ A, data=d, col=rangi2)
lines(x=A_seq, y=mu.mean, lwd=2)
shade(mu.PI, A_seq)

# approximate the posterior distribution for multiple regression
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2), #because vars are standardized, a should be close to 0
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data = d
)

precis(m5.3)
plot( coeftab(m5.1,m5.2,m5.3), par=c("bA","bM") )
# Once we know median age at marriage for a State, there is little
# or no additional predictive power in also knowing the
# rate of marriage in that State


# plotting multivariate posteriors ----------------------------------------




# predictor residual plots ------------------------------------------------

#to compute predictor residuals, model one predictor with the others
m5.4 <- quap(
  alist(M ~ dnorm(mu, sigma),
        mu <- a + bAM * A,
        a ~ dnorm(0, 0.2),
        bAM ~ dnorm(0, 0.5),
        sigma ~ dexp(1)
  ), data=d
)

# compute the residuals by subtracting the observed
# marriage rate in each State from the predicted rate
mu <- link(m5.4)
mu_mean <- apply(mu,2,mean)
mu_resid <- d$M - mu_mean
#positive residual = observed was higher than model prediction
#residuals are variation in outcome left over, after accounting for
#linear relationship of the predictor

#so this is the variation of marriage rates, accounting for linear
#effect of median age of marraige
plot(d$A, d$M)
abline(a=coef(m5.4)["a"], b=coef(m5.4)["b"])
plot(mu_resid, d$D)




# posterior predictive plot -----------------------------------------------

mu <- link(m5.3)

#summarize samples across cases
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

#simulate observations
D_sim <- sim(m5.3, n=1e4)
D_PI <- apply(D_sim, 2, PI)


plot(mu_mean ~ d$D,
     col=rangi2, ylim=range(mu_PI), xlab="observed divorce rate", ylab="predicted divorce rate")
abline(a=0, b=1, lty=2)
for(i in 1:nrow(d)) lines(rep(d$D[i],2), mu_PI[,i], col=rangi2)
identify( x=d$D , y=mu_mean , labels=d$Loc )




# simulating counterfactual plots ----------------------------------------------------

data("WaffleDivorce")
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)


m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1),
    ## A -> M
    M ~ dnorm(mu_M, sigma_M),
    mu_M <- aM + bAM*A,
    aM ~ dnorm(0, 0.2),
    bAM ~ dnorm(0, 0.5),
    sigma_M ~ dexp(1)
    
  ), data=d
)

precis(m5.3_A)

#M and A are strongly negatively associated. 
#If we interpret this causally, it indicates 
#that manipulating A reduces M.
#The goal is to simulate what would happen, if we manipulate A.
#So next we define a range of values for A.
A_seq <- seq(from=-2, to=2, length.out=30) 
#30 imaginary interventions on A

#use sim() to simulate observations from m5.3_A, 
#but for both M and D in that order.
# Because we have to simulate the influence of A on M
# before we simulate the joint influence of A and M on D.
sim_dat <- data.frame(A=A_seq)
s <- sim(m5.3_A, data=sim_dat, vars=c("M", "D"))

plot(sim_dat$A, colMeans(s$D),
     ylim=c(-2,2), type="l",
     xlab="manipulated A", ylab="counterfactual D")
shade(apply(s$D,2,PI), sim_dat$A)
mtext("total coutnerfactual effect of A on D")


plot(sim_dat$A, colMeans(s$M),
     ylim=c(-2,2), type="l",
     xlab="manipulated A", ylab="counterfactual M")
shade(apply(s$M, 2, PI), sim_dat$A)



sim_dat <- data.frame(M=seq(from=-2, to=2, length.out=30), A=0)
s <- sim(m5.3_A, data=sim_dat, vars="D")

plot(sim_dat$M, colMeans(s), ylim=c(-2,2), type="l",
     xlab="manipulated M", ylab="counterfactual D")
shade(apply(s,2,PI), sim_dat$M)
mtext("total counterfactual effect of M on D")



# simulating counterfactuals without sim() --------------------------------

A_seq <- seq(from=-2, to=2, length.out=30)
# simulate A on M
#extract posterior samples to simulate from
post <- extract.samples(m5.3_A)
M_sim <- with(post, sapply(1:30,
            function(i) rnorm(1e3, aM + bAM*A_seq[i], sigma_M)))

str(M_sim)
# samples of counterfactual M's at different A values(-2 to 2)
plot(A_seq, colMeans(M_sim), type="l",
     ylim=c(-2, 2), xlab="manipulated A", ylab="counterfactual A on M")

#simulate A on D too
D_sim <- with(post, sapply(1:30,
              function(i) rnorm(1e3, a + bA*A_seq[i] + bM*M_sim[,i], sigma)))
plot(A_seq, colMeans(D_sim), type="l",
  ylim=c(-2,2), xlab="manipulated A", ylab="total counterfactual effect of A on D")




# masked relationship -----------------------------------------------------
rm(list=ls())
data("milk"); d<- milk; rm(milk)

d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(d$mass)

# quap error due to missing values in N
dcc <- d[complete.cases(d$K, d$N, d$M),]
#bivariate regression, kcal ~ neocortex
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a + bN * N, # N= standardized neocortex
    a ~ dnorm(0, 1),
    bN ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data=dcc
)

#plot 50 prior regression lines
prior <- extract.prior(m5.5_draft) #extract samples from prior distributions
xseq <- c(-2,2)
mu <- link(m5.5_draft, post=prior, data=list(N=xseq))
# get values of mu for values of N at -2 and 2 SDs

plot(NULL, xlim=xseq, ylim=xseq,
     xlab='neocoretex perc (std)', ylab="kcal milk (std)")
for (i in 1:50) lines(xseq, mu[i, ])

for (i in 1:50) abline(a=prior$a[i], b=prior$bN[i])

# tightening the α prior so that it sticks closer to zero.

m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bN * N,
    a ~ dnorm(0, 0.2),
    bN ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=dcc
)

prior <- extract.prior(m5.5)
mu <- link(m5.5, post=prior, data=list(N=xseq))

plot(NULL, xlim=xseq, ylim=xseq,
     xlab='neocoretex perc (std)', ylab="kcal milk (std)")
for (i in 1:50) lines(xseq, mu[i, ])


# now posterior
precis(m5.5)

post <- extract.samples(m5.5)
plot(K ~ N, data=dcc, xlim=xseq, ylim=xseq,
     xlab='neocoretex perc (std)', ylab="kcal milk (std)")
abline(a=mean(post$a), b=mean(post$bN))


# alternative plot
xseq <- seq(from=min(dcc$N)-0.15, to=max(dcc$N)+0.15, length.out=30)
mu <- link(m5.5, data=list(N=xseq))
mu_mean <- apply(mu,2,mean) #colMeans
mu_PI <- apply(mu,2,PI) #col PI

plot(K ~ N, data=dcc)
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)


# now consider mass -------------------------------------------------------

# by using the logarithm of body mass here, we’re saying
# that we suspect that the magnitude of a mother’s body mass
# is related to milk energy, in a linear fashion


m5.6 <- quap(
  alist(
    K~ dnorm(mu, sigma),
    mu <- a + bM * M,
    a ~ dnorm(0, 0.2), 
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=dcc
)

precis(m5.6)
Mseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)-0.15, length.out=30)
mu <- link(m5.6, data=list(M=Mseq)) 
mu_mean <- colMeans(mu)
mu_pi <- apply(mu, 2, PI)

plot(K ~ M, data=dcc)
lines(Mseq, mu_mean, lwd=2)
shade(mu_pi, Mseq)


# log mass ----------------------------------------------------------------
dcc$log_M <- log(dcc$mass)
dcc$log_M <- standardize(dcc$log_M)

m5.6_log <- quap(
  alist(
    K~ dnorm(mu, sigma),
    mu <- a + bM * log_M,
    a ~ dnorm(0, 0.2), 
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=dcc
)

precis(m5.6_log)
post <- extract.samples(m5.6_log)

plot(K ~ log_M, data=dcc)
abline(a=mean(post$a), b=mean(post$bM))
#cant get shade


#errored way for some reason
Mseq <- seq(from=min(dcc$log_M)-0.15, to=max(dcc$log_M)-0.15, length.out=30)
mu <- link(m5.6_log, data=list(M=Mseq)) 
mu_mean <- colMeans(mu)
mu_pi <- apply(mu, 2, PI)

plot(K ~ log_M, data=dcc)
lines(Mseq, mu_mean, lwd=2)
shade(mu_pi, Mseq)














