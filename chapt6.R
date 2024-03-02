
library(rethinking)


# multicollinearity -------------------------------------------------------

#predicting height from leg length
N <- 100      
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5)  #leg as proportion of height
left_leg <- leg_prop*height + rnorm(N, 0, 0.02) #sim left leg as proportion + error
right_leg <- leg_prop*height + rnorm(N, 0, 0.02)
d <- data.frame(height, left_leg, right_leg)


lm(height ~ left_leg + right_leg, data=d) |> summary()

m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + bl*left_leg+ br*right_leg,
    a ~ dnorm(10,100),
    bl ~ dnorm(2,10),
    br ~ dnorm(2,10),
    sigma ~ dexp(1)
  ),data=d
)
precis(m6.1)
plot(precis(m6.1))

#look at joint posterior distribution for bl and br
post <- extract.samples(m6.1, size=1e4)
plot(bl ~ br, post, col=col.alpha(rangi2,0.1), pch=16)

sum_blbr <- post$bl + post$br
dens( sum_blbr , col=rangi2 , lwd=2 , xlab="sum of bl and br" )


# multicollinear milk -----------------------------------------------------

rm(list=ls())
data(milk);d<-milk;rm(milk)

d$Kcal <- standardize(d$kcal.per.g)
d$Fat <- standardize(d$perc.fat)
d$Lac <- standardize(d$perc.lactose)

#bivariate regressions with both predictors

m6.3 <- quap(
  alist(
    Kcal ~ dnorm(mu, sigma),
    mu <- a + bF * Fat,
    a ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
    
  ), data=d
)


m6.4 <- quap(
  alist(
    Kcal ~ dnorm(mu, sigma),
    mu <- a + bL * Lac,
    a ~ dnorm(0, 0.2),
    bL ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  
  ), data=d
)

precis(m6.3)
precis(m6.4)

#multiregression with 2 correlated predictors

m6.5 <- quap(
  alist(
    Kcal ~ dnorm(mu, sigma),
    mu <- a + bF*Fat + bL*Lac,
    a ~ dnorm(0, 0.2),
    bF ~ dnorm(0, 0.5),
    bL ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  )
  ,data=d
)

precis(m6.5)

pairs( ~ kcal.per.g + perc.fat + perc.lactose , data=d , col=rangi2 )


# simulating collinearity -------------------------------------------------
sim.coll <- function(r=0.9){
  d$x <- rnorm(nrow(d), mean=r * d$perc.fat,
               sd=sqrt((1-r^2)*var(d$perc.fat)))
  m <- lm(kcal.per.g ~ perc.fat + x, data=d)
  sqrt(diag(vcov(m)))[2] #sd of parameter
}

rep.sim.coll <- function(r=0.9, n=100){
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}

r.seq <- seq(from=0, to=0.99, by=0.01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r=z, n=100))

plot(stddev ~ r.seq, type="l", col=rangi2, lwd=2, xlab="correlation")



# post-treatment bias -----------------------------------------------------
rm(list=ls())
#plant exmaple (initial height, treatment, fungus, final height)
N <- 100
h0 <- rnorm(N, 10, 2) #simulate initial heights

#assign treatments and simulate fungus growth
treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size=1, prob=0.5-treatment*0.4) 
#less prob of fungus for treatment=1
# 100 tosses
h1 <- h0 + rnorm(N, 5-3*fungus)

#treatment (1) = taller
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d)

sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))
# 'data.frame': 10000 obs. of 1 variables:
#   mean   sd 5.5% 94.5%     histogram
# sim_p 1.03 0.26 0.67  1.48 ▁▁▃▇▇▃▁▁▁▁▁
#prior expects anything from 40% shrinkage, to 50% growth

m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p ~ dlnorm(0, 0.25),
    sigma ~ dexp(1)
  ),
  data=d
)
precis(m6.6)

#average of 40% growth

m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p ~ a + bT * treatment + bF * fungus,
    a ~ dlnorm(0, 0.25),
    bT ~ dnorm(0, 0.5),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=d
)

priors <- extract.prior(m6.7)
precis(priors)

precis(m6.7)

#now omit the post-treatment variable
m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p ~ a + bT * treatment,
    a ~ dlnorm(0, 0.25),
    bT ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=d
)
precis(m6.8)


# fungus and d-separation -------------------------------------------------
library(dagitty)

plant_dag <- dagitty("dag{
                     H_0 -> H_1
                     F -> H_1
                     T -> F
                     }")
coordinates(plant_dag) <- list(x=c(H_0=0, T=2, F=1.5, H_1=1),
                               y=c(H_0=0,T=0,F=0,H_1=0))
drawdag(plant_dag)

N <- 1000
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each=N/2)
M <- rbinom(N, size=1, prob=0.5) #bernolli
fungus <- rbinom(N, size=1, prob=0.5-treatment*0.4 + M*0.4)
h1 <- h0 + rnorm(N, 5 + 3*M)
d2 <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus,M=M)

cor(d2) |> round(2)
# M influences h1 and fungus (confounder)

# h1 ~ treatment + fungus
m6.7_2 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p ~ a + bT * treatment + bF * fungus,
    a ~ dlnorm(0, 0.25),
    bT ~ dnorm(0, 0.5),
    bF ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=d2
)
precis(m6.7_2)

# h1 ~ treatment
#now omit the post-treatment variable
m6.8_2 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),
    mu <- h0 * p,
    p ~ a + bT * treatment,
    a ~ dlnorm(0, 0.25),
    bT ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=d2
)
precis(m6.8_2)



# collider bias -----------------------------------------------------------

#agent-based simulation of age, happiness, and marriage
# H -> M <- A

# 1) each year, 20 people are born with uniformly distributed happiness
# 2) each year, each person ages one year. happiness does not change
# 3) at age 18, individuals can become married. the odds of marriage
# each year are proportional to an individuals happiness
# 4) once married, an individual remains married
# 5) after age 65, individuals leave the sample

d <- sim_happiness()
precis(d)
cor(d) |> round(2)

lm(happiness ~ age + married, data=d) |> summary()

d2 <- d[d$age>17, ]
d2$A <- (d2$age-18) / (65-18)
#rescales 0=18 to 1=65
d2$mid <- d2$married+1 #index variable for married

m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm(0, 1), #95% of mass from -2 to 2 happiness
    bA ~ dnorm(0, 2),  #95% of mass from -4 to 4 slope of age on happiness
    sigma ~ dexp(1)
  )
  ,data=d2
)

precis(m6.9 ,depth=2)

# mean happiness for unmarried 18 y/o = -0.23
# mean happiness for married 18 y/o = 1.26
# effect of age on happiness = -0.75 (older = unhappier)

# now leave out marriage
m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 1),
    bA ~ dnorm(0, 2),
    sigma ~ dexp(1)
  
  ), data=d2
)
precis(m6.10)
# no effect of age on happiness


# haunted DAGs (unobserved common cause) ----------------------------------
rm(list=ls())
N <- 200
b_GP <- 1 #direct effect of G on P
b_GC <- 0 #direct effect of G on C
b_PC <- 1 #direct effect of P on C
b_U <- 2 #direct effect of U on P and C

set.seed(1)
U <- 2*rbern(N, 0.5) -1 # -1 or 1
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U*U) #P as function of G and U
C <- rnorm(N, b_GC*G + b_PC*P + b_U*U)

d <- data.frame(C=C, P=P, G=G, U=U)

# C ~ G + P  #look at influence of grandparents
m6.11 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G,
    a ~ dnorm(0, 1),
    b_PC ~ dnorm(0, 1),
    b_GC ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data=d
)

precis(m6.11)


# add U:   C ~ G + P + U
m6.12 <- quap(
  alist(
    C ~ dnorm(mu, sigma),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm(0, 1),
    c(b_PC, b_GC, b_U) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data=d
)
precis(m6.12)











