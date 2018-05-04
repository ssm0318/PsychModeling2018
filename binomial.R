# ra_prospect_stan_singleSubj.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr), Apr 2018

rm(list=ls())  # remove all variables 

library(rstan)

# set working directory where all the files exist
setwd("~/Dropbox/5-1/ExpPsych/HW/HW3")

# source HDIofMCMC.R to calculate HDI
source("HDIofMCMC.R") 

# read the data file
T = 6  # number of trials per subject (=6)
k = c(11, 18, 15, 39, 32, 43)
n = c(20, 25, 30, 55, 45, 57)
numIter = 100             # number of iterations to find global minimum values
numPars = 1               # number of parameters

# use first subject only
dataList <- list(
  T       = T,
  k       = k,
  n       = n
)

# run!
output = stan("binomial.stan", data = dataList, pars = c("theta"),
              iter = 2000, warmup=1000, chains=2, cores=2)

# traceplot
traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

# plot posteriors 
hist(parameters$theta)

# 95% HDI of theta
HDIofMCMC(parameters$theta, credMass = 0.95)

