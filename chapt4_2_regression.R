
library(rethinking)

data("Howell1")
d <- Howell1
d2 <- d[d$age >= 18, ]


# height by weight --------------------------------------------------------
plot(d2$height ~ d2$weight)

xbar <- mean(d2$weight)


m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),   #likelihood function
    mu <- a + b * (weight-xbar), #linear model
    a ~ dnorm(178, 20),          #a prior    
    b ~ dlnorm(0, 1),            #b prior
    sigma ~ dunif(0, 50)         # sigma prior
  ),
  data = d2)






# tables of marginal distributions ----------------------------------------
precis(m4.3)

# lm(d2$height~d2$weight) |> summary()

vcov(m4.3) |> round(3)

#very little covariation among parameters, results from centering
pairs(m4.3)


# plotting posterior inference against the data ---------------------------

#plotting line implied by mean posterior values of a & b
plot(height ~ weight, data = d2, col=rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map*(x-xbar), add=T)

#adding uncertainty around the mean
post[1:5,]

#see how sample size affects the posterior distribution
N <- 100
dN <- d2[1:N, ]
mN <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a+b*(weight-xbar),
    a ~ dnorm(178,20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data = dN)

#extract 20 smaples to see uncertainty
post <- extract.samples(mN, n=20)

#raw data and sample size
plot(dN$weight, dN$height, 
     xlim=range(d2$weight), ylim=range(d2$height),
     col=rangi2, xlab="weight", ylab="height")
mtext(concat("N= ", N))

#plot the lines
for(i in 1:20)
  curve(post$a[i] + post$b[i]*(x-mean(dN$weight)),
        col=col.alpha("black", 0.3), add=T)


# posterior prediction of mu for individual who weighs 50kg ---------------------

post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b * (50-xbar)
# vector of predicted means, one for each random sample from the posterior
dens(mu_at_50, col=rangi2, lwd=2, xlab="mu|weight=50")
#quadratic approx posterior distribution of
#the mean height mu when weight=50kg
#relative plausibility of difference values of mu
PI(mu_at_50, prob=0.89)
#     5%      94% 
#   158.57 159.66
#the central 89% ways for the model to produce the data
#placed the average height between 159-160cm
#(conditional on the model and data) assuming weight=50kg

#link will sample from quap posterior, then compute mu for 
#each case in the data&sample from the posterior
mu <- link(m4.3)
str(mu)
#matrix mu values.
#each row = sample from posterior. each column = case in the data
#link provides a posterior distribution for each case. 
#so we have a distribution (1000 samples) of mu for each individual (352)

#we want a distribution for mu for each unique weight tho
weight.seq <- seq(from=25, to=70, by=1)
mu <- link(m4.3, data=data.frame(weight=weight.seq))
#compute mu for each sample from posterior, and for each weight in weight.seq
str(mu)
#46 columns = 46 different weight values
plot(height ~ weight, d2, type="n")
for(i in 1:100)
  points(weight.seq, mu[i,], pch=16, col=col.alpha(rangi2,0.1))
#pile of computed mu values are shown for each weight value in weight.seq

#summarize distribution of mu
mu.mean <- apply(mu, 2, mean)
#compute mean of each column of mu matrix (mean of computed means for each weight)
# average µ at each weight value
mu.PI <- apply(mu, 2, PI, prob=0.89)
plot(height~weight, data=d2, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)


#decomposing link()
post <- extract.samples(m4.3)
mu.link <- function(weight) post$a + post$b*(weight-xbar)
weight.seq <- seq(from=25, to=70, by=1)
mu <- sapply(weight.seq, mu.link)



# prediction intervals for actual heights ----------------------------------------------------
# hi ~ N(mu, sigma)

sim.height <- sim(m4.3, data=list(weight=weight.seq), n=10000)
str(sim.height)
#matrix of simulated heights based on computed mu's for each weight

height.PI <- apply(sim.height, 2, PI, prob=0.89)
height.PI

#plot
plot(height ~ weight, d2, col=col.alpha(rangi2, 0.5))
lines(weight.seq, mu.mean) #average line
shade(mu.PI, weight.seq) #PI for average line
shade(height.PI, weight.seq) #PI region for simulated heights

#decompose sim()
post <- extract.samples(m4.3)
# post <- MASS::mvrnorm(n=1e4, mu=coef(m4.3), Sigma=vcov(m4.3))
weight.seq<- (25:70)
sim.height <- sapply(weight.seq, function(weight)
  rnorm(n = nrow(post),
        mean=post$a +post$b*(weight-xbar),
        sd=post$sigma
        ))
str(sim.height)



#curves from lines
####################
# polynomial regression ---------------------------------------------------
rm(list=ls())
library(rethinking)
data("Howell1")
d <- Howell1

plot(height~weight, data=d)
#curved relationship

d$weight_s <- (d$weight - mean(d$weight))/sd(d$weight)
d$weight_s2 <- d$weight_s^2

m4.5 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178,20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,1),
    sigma ~ dunif(0,50)
  ), data=d)

precis(m4.5)
# since the relationship between the outcome height and 
# the predictor weight depends
# upon two slopes, b1 and b2, it isn’t so easy 
# to read the relationship off a table of coefficients

#calculate mean relationship and 89% intervals of mean and predictions
weight.seq <- seq(from=-2.2, to=2.2, length.out=30)
pred_dat <- list(weight_s=weight.seq, weight_s2=weight.seq^2)
mu <- link(m4.5, data=pred_dat)
mu.mean <- apply(mu,2,mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)

sim.height <- sim(m4.5, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

#plot
plot(height ~ weight_s, d, col=col.alpha(rangi2,0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)



# cubic -------------------------------------------------------------------





# b-splines ---------------------------------------------------------------
rm(list=ls())
library(rethinking)
data("cherry_blossoms")
d <- cherry_blossoms

precis(d)
plot(doy ~ year, data=d)

#knots
d2 <- d[complete.cases(d$doy),]
num_knots <- 15
knot_list <- quantile(d2$year, probs=seq(0,1, length.out=num_knots))
#evenly spaced 15 dates

#polynomial degree. 3 = cubic splines
library(splines)

#get basis functions
B <- bs(d2$year,
        knots=knot_list[-c(1,num_knots)], #drop first and last knot
        degree=3, intercept=T)

#row = year corresponding to rows in d2
#columns = basis functions
plot(NULL, xlim=range(d2$year), ylim=c(0,1), xlab="year",ylab="basis")
for(i in 1:ncol(B)) lines(d2$year, B[, i])

#get parameter weights for each basis function
m4.7 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w, 
    a ~ dnorm(100,10),
    w ~ dnorm(0,10),
    sigma ~ dexp(1)),
  data=list(D=d2$doy, B=B),
  start=list(w=rep(0, ncol(B))))


precis(m4.7, depth = 2)

#plot posterior predictions
post <- extract.samples(m4.7)
w <- apply(post$w, 2, mean)
plot(NULL, xlim=range(d2$year), ylim=c(-6,6), 
     xlab="year", ylab="basis * weight")
for (i in 1:ncol(B)) lines(d2$year, w[i]*B[,i])

#97% posterior interval for mu at each year
mu <- link(m4.7)
mu_PI <- apply(mu,2,PI,0.97)
plot(d2$year, d2$doy, col=col.alpha(rangi2,0.3),pch=16)
shade(mu_PI, d2$year, col=col.alpha("black",0.5))


























