
remove.packages("rstan")

if (file.exists(".RData")) file.remove(".RData")

install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
library(rstan)


example(stan_model, package = "rstan", run.dontrun = TRUE)

options(mc.cores = parallel::detectCores())








