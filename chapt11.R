
library(rethinking)

data("chimpanzees"); d <- chimpanzees ; rm(chimpanzees)


# build index variables for all 4 binary predictor combinations

d$treatment <- 1 + d$prosoc_left + 2*d$condition
xtabs(~ treatment + prosoc_left + condition, d)



# logistic regression
m11.1 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a,
    a ~ dnorm(0, 10) #unreasonably wide prior
  ),data=d
)

set.seed(1999)
prior <- extract.prior(m11.1, n=1e4) # dnorm(0, 10)
# convert parameter to outcome scale
# p <- inv_logit()
p = inv_logit(prior$a)
plot(density(p, adj=0.1))
dens(p, adj=0.1)

# reasonable prior
m11.1.1 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a,
    a ~ dnorm(0, 1)
  ),data=d
)

prior <- extract.prior(m11.1.1, n=1e4)
p = inv_logit(prior$a)
dens(p, adj=0.1)

# prior for treatment effect
m11.2 <- quap(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a + b[treatment],
    a ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 1)
  ),data=d
)
set.seed(1999)
prior <- extract.prior(m11.2, n=1e4)
p <- sapply(1:4, function(k) inv_logit(prior$a + prior$b[, k]))
# prior probability of pulling left for each treatment condition
# a + b[treatment]
# prior differences between treatments
dens(abs(p[,1]-p[,2]), adj=0.1)
mean(abs(p[,1]-p[,2]))

# trimmed data list

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment)
)

# do Markov chain monte carlo approximation of the posterior
m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a[actor] + b[treatment],
    a[actor] ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ),data=dat_list, chains=4, log_lik=T #to compute necessary values for PSIS and WAIC
)
precis(m11.4, 2)

# examine chimp-specific intercepts (lever tendency)
post <- extract.samples(m11.4)
p_left <- inv_logit(post$a)
plot(precis(as.data.frame(p_left)), xlim=c(0,1))

# exmaine treatment effects
labs <- c("R/N","L/N","R/P","L/P")
plot(precis(m11.4, depth=2, pars="b"), labels=labs)
# calculate differences between no-partner/partner
diffs <- list(
  db12 = post$b[,1] - post$b[,3],
  db24 = post$b[,2] - post$b[, 4]
)
# plot contrasts, scale = log-odds of pulling left lever
plot(precis(diffs))

# posterior prediction check
# summarize proportions of left pulls for each chimp in each treatment
# and then plot against posterior predictions
pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
pl[,] # row=chimps, columns=treatment

# plot
plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("actor ",j) , xpd=TRUE )
for ( j in (1:7)[-2] ) {
  lines( (j-1)*4+c(1,3) , pl[j,c(1,3)] , lwd=2 , col=rangi2 )
  lines( (j-1)*4+c(2,4) , pl[j,c(2,4)] , lwd=2 , col=rangi2 )
}
points( 1:28 , t(pl) , pch=16 , col="white" , cex=1.7 )
points( 1:28 , t(pl) , pch=c(1,1,16,16) , col=rangi2 , lwd=2 )
yoff <- 0.01
text( 1 , pl[1,1]-yoff , "R/N" , pos=1 , cex=0.8 )
text( 2 , pl[1,2]+yoff , "L/N" , pos=3 , cex=0.8 )
text( 3 , pl[1,3]-yoff , "R/P" , pos=1 , cex=0.8 )
text( 4 , pl[1,4]+yoff , "L/P" , pos=3 , cex=0.8 )
mtext( "observed proportions\n" )


# posterior predictions
dat <- list( actor=rep(1:7,each=4) , treatment=rep(1:4,times=7) )
p_post <- link( m11.4 , data=dat )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )


# relative effects
# calculate proportional odds by exponentiating the parameter

# calculate proportional odds of switching from treatment 4 to 2
post <- extract.samples(m11.4)
mean(exp(post$b[,4] - post$b[,2]))
# on average, switching from treatment 4 to 2, multiples the odds of left_pull
# by 0.9259718, an 8% reduction in odds



# aggregated binomial -----------------------------------------------------
rm(list=ls())
data("chimpanzees"); d <- chimpanzees; rm(chimpanzees)

d$treatment <- 1 + d$prosoc_left + 2*d$condition
d$side <- d$prosoc_left + 1 #right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2

d_aggregated <- aggregate(d$pulled_left,
                          list(treatment=d$treatment, actor=d$actor,
                               side=d$side, cond=d$cond),
                          sum)

colnames(d_aggregated)[5] <- "left_pulls"


dat <- with(d_aggregated, list(
  left_pulls=left_pulls,
  treatment=treatment,
  actor=actor,
  side=side,
  cond=cond))


m11.6 <- ulam(
  alist(
    left_pulls ~ dbinom(18, p),
    logit(p) <- a[actor] + b[treatment],
    a[actor] ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ),data=dat,chains=4,log_lik=T
)

precis(m11.6, 2)
compare(m11.6, m11.4, func=PSIS)


# -------------------------------------------------------------------------

rm(list=ls())
data(UCBadmit); d <- UCBadmit; rm(UCBadmit)

# model admissions decisions, with gender as predictor

dat_list <- with(d, list(
  admit=admit,
  applications=applications,
  gender = ifelse(applicant.gender == "male", 1, 2)
))

# admit ~ gender
m11.7 <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a[gender],
    a[gender] ~ dnorm(0, 1.5)
  ),data=dat_list, chains=4
)


precis(m11.7,2)
# posterior mean for males > females
# compute contrast on logit scale and absolute (outcome) scale
post <- extract.samples(m11.7)
diff_a <- post$a[, 1]  - post$a[, 2]
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
# males minus females.  values > 0 = male advantage
precis(list(diff_a=diff_a, diff_p=diff_p))
# diff_a = 0.61 = 
# diff_p = 0.14 = women 14% lower adds of admissions

# posterior check
postcheck( m11.7 )
# draw lines connecting points from same dept
for ( i in 1:6 ) {
  x <- 1 + 2*(i-1)
  y1 <- d$admit[x]/d$applications[x]
  y2 <- d$admit[x+1]/d$applications[x+1]
  lines( c(x,x+1) , c(y1,y2) , col=rangi2 , lwd=2 )
  text( x+0.5 , (y1+y2)/2 + 0.05 , d$dept[x] , cex=0.8 , col=rangi2 )
}

# different admit rates for each department, misleading posterior estimates
# add department variable
dat_list$department <- as.integer(d$dept)

m11.8 <- ulam(
  alist(
    admit ~ dbinom(applications, p),
    logit(p) <- a[gender] + b[department],
    a[gender] ~ dnorm(0, 1.5),
    b[department] ~ dnorm(0, 1.5)
  ),data=dat_list,chains=4
)

precis(m11.8, 2)
# now males are lower average odds
# calculate contrasts
post <- extract.samples(m11.8)
diff_a <- post$a[, 1] - post$a[, 2]
diff_p <- inv_logit(post$a[, 1]) - inv_logit(post$a[, 2])
precis(list(diff_a=diff_a, diff_p=diff_p))

# men and women applied to different departments
pg <- with(dat_list, sapply(1:6, function(k)
  applications[department==k]/sum(applications[department==k])))

rownames(pg) <- c("male", "female")
colnames(pg) <- unique(d$dept)
round(pg, 2)



# poisson dist ------------------------------------------------------------

# employed 1000 monks (size), on any day about 1 finished a manuscript (p)
# over 10,000 days (n)
y <- rbinom(n=1e5, size=1000, p=1/1000)
c(mean(y), var(y))
plot(y)
hist(y, breaks=100)


hist(rbinom(n=1e5, size=10, p=0.5))


# poisson example - oceanic tool kits -------------------------------------

data(Kline); d <- Kline; rm(Kline)

d$log_pop <- scale(log(d$population))
d$contact_id <- ifelse(d$contact=="high", 2, 1)

# model form
# m <- quap(alist(
#   total_tools ~ dpois(lambda),
#   log(lambda) <- a[contact_id] + b[contact_id]*log_pop,
#   a[contact_id] ~ ,
#   b[contact_id] ~
# ))

# wide prior for log scale?
# T[i] ~ Poisson(lambda[i])
# log(lambda[i]) = a
# a ~ dnorm(0, 10)
curve(dlnorm(x, mean=0, sd=10), from = 0, to=100, n=200)
# zero tools on average. doesnt make sense
# try narrower prior
curve(dlnorm(x, mean=3, sd=0.5), from=0, to=100, n=200)

# now prior for effect of log_population
N <- 100
a <- rnorm(N, mean=3, sd=0.5)
b <- rnorm(N, mean=0, sd=10)
plot(NULL, xlim=c(-2,2), ylim=c(0,100),
     ylab="total_tools",xlab="z-scores for log_pop")
for(i in 1:N) curve(exp( a[i] + b[i]*x), add=T, col=grau())

# unstandardize log_pop to preserve natural 0
x_seq <- seq(from=log(100), to=log(2e5), length.out=100)
lambda <- sapply(x_seq, function(x) exp(a + b*x))
plot( NULL , xlim=range(x_seq) , ylim=c(0,500) , xlab="log population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( x_seq , lambda[i,] , col=grau() , lwd=1.5 )

# natural population scale
plot( NULL , xlim=range(exp(x_seq)) , ylim=c(0,500) , xlab="population" ,
      ylab="total tools" )
for ( i in 1:N ) lines( exp(x_seq) , lambda[i,] , col=grau() , lwd=1.5 )

# posteriors
dat <- list(
  tools = d$total_tools,
  log_pop = d$log_pop,
  cid = d$contact_id
)

# intercept only
m11.9 <- ulam(
  alist(
    tools ~ dpois(lambda),
    log(lambda) <- a,
    a ~ dnorm(3, 0.5)
  ),data=dat, chains=4, log_lik=T
)

# interaction model
m11.10 <- ulam(
  alist(
    tools ~ dpois(lambda),
    log(lambda) <- a[cid] + b[cid]*log_pop,
    a[cid] ~ dnorm(3, 0.5),
    b[cid] ~ dnorm(0, 0.2)
  ), data=dat, chains=4, log_lik=T
)

compare(m11.9, m11.10, func=PSIS)

# plot posterior predictions to find influential point k
k <- PSIS(m11.10, pointwise=T)$k
plot(dat$log_pop, dat$tools, xlab="log pop (std)", ylab="total tools",
     col=rangi2, pch=ifelse(dat$cid==1, 1, 16), lwd=2,
     ylim=c(0,75), cex=1+normalize(k))

# compute predictions
ns <- 100
p_seq <- seq(from= -1.4, to=3, length.out=ns)

# predictions for cid=1 (low contact)
lambda <- link(m11.10, data=data.frame(log_pop=p_seq, cid=1))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(p_seq, lmu, lty=2, lwd=1.5)
shade(lci, p_seq, xpd=T)

# predictions for cid=2 (high contact)
lambda <- link(m11.10, data=data.frame(log_pop=p_seq, cid=2))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(p_seq, lmu, lty=1, lwd=1.5)
shade(lci, p_seq, xpd=T)



# poisson -----------------------------------------------------------------

num_days <- 30
y <- rpois(num_days, 1.5)

num_weeks <- 4
y_new <- rpois(num_weeks, 0.5*7)

y_all <- c(y, y_new)
exposure <- c(rep(1,30), rep(7, 4))
monastary <- c(rep(0, 30), rep(1, 4))
d <- data.frame(y=y_all, days=exposure, monastary=monastary)

d$log_days <- log(d$days)

m11.12 <- quap(
  alist(
    y ~ dpois(lambda),
    #              log(tau) is the rate of exposure. here, it is days (1, 7)
    log(lambda) <- log_days + a + b*monastary,
    a ~ dnorm(0, 1),
    b ~ dnorm(0, 1)
  ),data=d
)

precis(m11.12)
# a=0.31; average log rate of production
# b= -0.88; monastary 1 is slower than 0 by a log rate of 0.88

# to get parameters back on outcome scale, exp() the samples
post <- extract.samples(m11.12) # samples from the posterior probability dist
lambda_old <- exp(post$a)
lambda_new <- exp(post$a + post$b)
precis(data.frame(lambda_old, lambda_new))
# new monastary produces half as many manuscripts per DAY



# multinomial -------------------------------------------------------------


# simulate career choices for 500 ppl
N <- 500
income <- c(1, 2, 5)
score <- 0.5*income
# convert to probabilities
p <- softmax(score[1], score[2], score[3])

career <- rep(NA, N)
for(i in 1:N) career[i] <- sample(1:3, size=1, prob=p)

# 
# # simulate career choices among 500 individuals
# N <- 500 # number of individuals
# income <- c(1,2,5) # expected income of each career
# score <- 0.5*income # scores for each career, based on income
# # next line converts scores to probabilities
# p <- softmax(score[1],score[2],score[3])
# # now simulate choice
# # outcome career holds event type values, not counts
# career <- rep(NA,N) # empty vector of choices for each individual
# # sample chosen career for each individual
# set.seed(34302)
# for ( i in 1:N ) career[i] <- sample( 1:3 , size=1 , prob=p )
# 
# 
# code_m11.13 <- "
# data{
# int N; // number of individuals
# int K; // number of possible careers
# int career[N]; // outcome
# vector[K] career_income;
# }
# parameters{
# vector[K-1] a; // intercepts
# real<lower=0> b; // association of income with choice
# }
# model{
# vector[K] p;
# vector[K] s;
# a ~ normal( 0 , 1 );
# b ~ normal( 0 , 0.5 );
# s[1] = a[1] + b*career_income[1];
# s[2] = a[2] + b*career_income[2];
# s[3] = 0; // pivot
# p = softmax( s );
# career ~ categorical( p );
# }
# "
# dat_list <- list(N=N, K=3, career=career, career_income=income)
# m11.13 <- stan(model_code=code_m11.13)
# 
# precis(m11.13, 2)




# binomial model of overall admissions probability ------------------------
data("UCBadmit"); d <- UCBadmit; rm(UCBadmit)

m_binom <- quap(
  alist(
    admit ~ dbinom(applications, p), 
    logit(p) <- a,
    a ~ dnorm(0, 1.5)
  ),data=d
)
precis(m_binom)

# poisson
dat <- list( admit=d$admit , rej=d$reject )

m_pois <- ulam(
  alist(
    admit ~ dpois(lambda1),
    rej ~ dpois(lambda2),
    log(lambda1) <- a1,
    log(lambda2) <- a2,
    c(a1,a2) ~ dnorm(0,1.5)
  ), data=dat , chains=3 , cores=3 )

precis(m_pois)

# binom probability of admissions (posterior mean)
inv_logit(coef(m_binom))

# poisson implied probability of admissions
# p = lambda1/lambda1+lambda2 = exp(a1) / exp(a1)+exp(a2)

k <- coef(m_pois)
a1 <- k['a1']; a2 <- k['a2']
exp(a1) / (exp(a1) + exp(a2))











