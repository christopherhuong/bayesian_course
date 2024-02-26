
library(rethinking)

data("Howell1"); d <-Howell1; rm(Howell1)



mu_female <- rnorm(1e4, 178, 20)
mu_male <- rnorm(1e4, 178, 20) + rnorm(1e4,0,10)
precis(data.frame(mu_female, mu_male))
# The prior for males is wider, because it uses both parameters.



# index variables ---------------------------------------------------------
d$sex <- ifelse(d$male==1, 2, 1)
#male=2, female=1

m5.8 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a[sex],
    a[sex] ~ dnorm(178, 20),
    sigma ~ dunif(0, 50)
  ),
  data=d
)

precis(m5.8, depth=2)

#now calculate differences in M F height, using posterior samples
post <- extract.samples(m5.8, size=1e4)

post$diff_fm <- post$a[,1] - post$a[,2] #female - male
precis(post, depth=2)



# multiple categories -----------------------------------------------------
rm(list=ls())
data(milk);d<-milk;rm(milk)

levels(d$clade)
#index value for each category
d$clade_id <- as.integer(d$clade)

d$K <- standardize(d$kcal.per.g)

m5.9 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm(0, 5),
    sigma ~ dexp(1)
  )
  ,data=d
)

labels<- paste("a[", 1:4, "]: ", levels(d$clade), sep="")

plot(precis(m5.9, depth=2, pars="a"), labels=labels,
     xlab="expected kcal")


#add another categorical variable
set.seed(63)
d$house <- sample(rep(1:4, each=8), size=nrow(d))

m5.10 <- quap(
  alist(
    K ~ dnorm(mu, sigma),
    mu <- a[clade_id] + h[house],
    a[clade_id] ~ dnorm(0, 0.5),
    h[house] ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ),
  data=d
)

precis(m5.10,depth=2)














