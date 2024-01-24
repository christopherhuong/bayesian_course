

library(rstan)


schools_d <- list(J = 8,
                  y = c(28, 8, -3, 7, -1, 1, 18, 12),
                  sigma = c(15, 10, 16, 11, 9, 11, 10, 18))



fit <- stan(file = 'schools.stan', data = schools_d)


# The object fit, returned from function stan is an S4 object of class stanfit. Methods such as print, plot, and pairs are associated with the fitted result so we can use the following code to check out the results in fit. print provides a summary for the parameter of the model as well as the log-posterior with name lp__ (see the following example output). For more methods and details of class stanfit, see the help of class stanfit.
# In particular, we can use the extract function on stanfit objects to obtain the samples. extract extracts samples from the stanfit object as a list of arrays for parameters of interest, or just an array. In addition, S3 functions as.array, as.matrix, and as.data.frame are defined for stanfit objects (using help("as.array.stanfit") to check out the help document in R).

print(fit)
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = T)  #return list of arrays
mu <- la$mu

# return an array of three dimensions: iterations, chains, parameters

a <- extract(fit, permuted = F)

# use S3 fxs on stanfit objs
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
a2

x <- rnorm(100, mean=0, sd=15)
hist(x)



