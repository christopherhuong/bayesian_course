
library(rethinking)



data(rugged); d <- rugged

# log the outcome
d$log_gdp <- log(d$rgdppc_2000)
# drop missing
dd <- d[complete.cases(d$rgdppc_2000),]
# rescale
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
# divide by mean: 1=mean, 0.80=80% of mean, 1.2 = 120% of mean
dd$rugged_std <- dd$rugged / max(dd$rugged)
# divide by max: 0 = flattest, 1=most rugged

# linear model
m8.1 <- quap(alist(
  log_gdp_std ~ dnorm(mu, sigma),
  mu <- a + b * (rugged_std - mean(rugged_std)),
  a ~ dnorm(1,1),
  b ~ dnorm(0, 1),
  sigma ~ dexp(1)
),data=dd)

# plot priors
priors <- extract.prior(m8.1)

plot(NULL, xlim=c(0,1), ylim=c(0.5,1.5),
     xlab="ruggedness", ylab="log GDP")
abline(h=min(dd$log_gdp_std), lty=2)
abline(h=max(dd$log_gdp_std), lty=2)

rugged_seq <- seq(from=-0.1, to=1.1, length.out=30)
mu <- link(m8.1, post=priors, data=data.frame(rugged_std=rugged_seq))
for(i in 1:50) lines(rugged_seq, mu[i, ], col=col.alpha("black",0.3))

# better priors
m8.1 <- quap(alist(
  log_gdp_std ~ dnorm(mu, sigma),
  mu <- a + b * (rugged_std - mean(rugged_std)),
  a ~ dnorm(1, 0.1), #intercept closer to the mean
  b ~ dnorm(0, 0.3),
  sigma ~ dexp(1)
),data=dd)



precis(m8.1)

# index variable to estimate different intercept for african continent
dd$cid <- ifelse(dd$cont_africa==1, 1,2)

m8.2 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b * (rugged_std - mean(dd$rugged_std)),
    a[cid] ~ dnorm(1, 0.1),
    b ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ),data=dd)

precis(m8.2, depth=2)
# compare models with WAIC
compare(m8.1, m8.2, func=WAIC)
# using indicator variable to separately estimate an intercept 
# (mean GDP at mean ruggedness)

# posterior contrast between the two variables
post <- extract.samples(m8.2)
diff_a1_a2 <- post$a[,1] - post$a[,2]
PI(diff_a1_a2)
# different is reliably below zero

# plot m8.2 posterior predictions
rugged_seq <- seq(from=-0.1, to=1.1, length.out=30)
# compute mu over samples, fixing cid=2, then cid=1
mu.notafrica <- link(m8.2, 
                     data=data.frame(cid=2, rugged_std=rugged_seq))
mu.africa <- link(m8.2,
                  data=data.frame(cid=1, rugged_std=rugged_seq))
#summarize to means and intervals
mu.notafrica_mu <- apply( mu.notafrica , 2 , mean )
mu.notafrica_ci <- apply( mu.notafrica , 2 , PI , prob=0.97 )
mu.africa_mu <- apply( mu.africa , 2 , mean )
mu.africa_ci <- apply( mu.africa , 2 , PI , prob=0.97 )





# plot slopes for africa vs not africa
plot(dd$log_gdp_std ~ dd$rugged_std,
     col=dd$cont_africa+1, pch=16)
legend("topright", legend = c("not africa", "africa"),
       col=c("black","red"), pch = 16)

lines(x=rugged_seq, y=mu.notafrica_mu, col="black", lwd=2)
shade(mu.notafrica_ci, rugged_seq)
lines(x=rugged_seq, y=mu.africa_mu, col="red", lwd=2)
shade(mu.africa_ci, rugged_seq)


# -------------------------------------------------------------------------


# interaction term to estimate slope conditional on continent
m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm(mu, sigma),
    mu <- a[cid] + b[cid] * (rugged_std-mean(rugged_std)),
    a[cid] ~ dnorm(1, 0.1),
    b[cid] ~ dnorm(0, 0.3),
    sigma ~ dexp(1)
  ),data=dd)

precis(m8.3, depth=2)

# posterior contrast between the two variables
post <- extract.samples(m8.3)
diff_a1_a2 <- post$b[,1] - post$b[,2]
PI(diff_a1_a2)


# plot m8.3 posterior predictions
rugged_seq <- seq(from=-0.1, to=1.1, length.out=30)
# compute mu over samples, fixing cid=2, then cid=1
mu.notafrica <- link(m8.3, 
                     data=data.frame(cid=2, rugged_std=rugged_seq))
mu.africa <- link(m8.3,
                  data=data.frame(cid=1, rugged_std=rugged_seq))
#summarize to means and intervals
mu.notafrica_mu <- apply( mu.notafrica , 2 , mean )
mu.notafrica_ci <- apply( mu.notafrica , 2 , PI , prob=0.97 )
mu.africa_mu <- apply( mu.africa , 2 , mean )
mu.africa_ci <- apply( mu.africa , 2 , PI , prob=0.97 )





plot(dd$log_gdp_std ~ dd$rugged_std,
     col=dd$cont_africa+1, pch=16)
legend("topright", legend = c("not africa", "africa"),
       col=c("black","red"), pch = 16)

lines(x=rugged_seq, y=mu.notafrica_mu, col="black", lwd=2)
shade(mu.notafrica_ci, rugged_seq)
lines(x=rugged_seq, y=mu.africa_mu, col="red", lwd=2)
shade(mu.africa_ci, rugged_seq)

# compare predictive accuracy
compare(m8.1, m8.2, m8.3, func=PSIS)

# interactions are symmetrical
rugged_seq <- seq(from=-0.2,to=1.2,length.out=30)

muA <- link( m8.3 , data=data.frame(cid=1,rugged_std=rugged_seq) )
muN <- link( m8.3 , data=data.frame(cid=2,rugged_std=rugged_seq) )


delta <- muA - muN
delta_m <- apply(delta, 2, mean)
delta_ci <- apply(delta, 2, PI, prob=0.97)

plot(NULL, ylim=c(-0.3, 0.2), xlim=c(0, 1),
     xlab="ruggedness", ylab="diff GDP")
lines(x=rugged_seq, y=delta_m, lwd=2)
shade(delta_ci, rugged_seq)
abline(h=0, lty=5)
text(x=0.2,y=0.05, "Africa higher GDP")
text(x=0.2, y=-0.05, "Africa lower GDP")


# -------------------------------------------------------------------------
rm(list=ls())
library(rethinking)

data(tulips); d <- tulips; rm(tulips)
d$blooms_std <- d$blooms / max(d$blooms) #scale by maximum
d$water_cent <- d$water - mean(d$water) #center at 0
d$shade_cent <- d$shade - mean(d$shade) #center at 0


m8.4 <- quap(
  alist(
    blooms_std ~ dnorm(mu, sigma),
    mu <- a + bW*water_cent + bS*shade_cent,
    a ~ dnorm(0.5, 0.25),
    bW ~ dnorm(0, 0.25),
    bS ~ dnorm(0, 0.25),
    sigma ~ dexp(1)
  ),data=d)

# interaction term
m8.5 <- quap(
  alist(
    blooms_std ~ dnorm( mu , sigma ) ,
    mu <- a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent ,
    a ~ dnorm( 0.5 , 0.25 ) ,
    bw ~ dnorm( 0 , 0.25 ) ,
    bs ~ dnorm( 0 , 0.25 ) ,
    bws ~ dnorm( 0 , 0.25 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )

# posterior predictions of B ~ W at different values of S with no interaction term
par(mfrow=c(1,3))
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx], d$blooms_std[idx],
       xlim=c(-1,1), ylim=c(0, 1),
       xlab="water", ylab="blooms", pch=16, col=rangi2)
  mu <- link(m8.4, data=data.frame(shade_cent=s, water_cent=-1:1))
  for(i in 1:20) lines(-1:1, mu[i,], col=col.alpha("black",0.3))
}



# posterior predictions B ~ W with interaction term W * S
par(mfrow=c(1,3))
for(s in -1:1){
  idx <- which(d$shade_cent==s)
  plot(d$water_cent[idx], d$blooms_std[idx],
       xlim=c(-1,1), ylim=c(0,1),
       xlab="water",ylab="blooms", pch=16, col=rangi2)
  mu <- link(m8.5, data=data.frame(shade_cent=s, water_cent=-1:1))
  for(i in 1:20) lines(-1:1, mu[i,], col=col.alpha("black", 0.3))
}
























