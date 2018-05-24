# ra_prospect_stan_singleSubj.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr), Apr 2018

rm(list=ls())  # remove all variables 

library(rstan)

# set working directory where all the files exist
setwd("~/Dropbox/5-1/ExpPsych/HW/HW4")

# source HDIofMCMC.R to calculate HDI
source("HDIofMCMC.R") 

# read the data file
dat = read.table("ra_exampleData.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
Tsubj = unname(table(dat$subjID))
T = table(dat$subjID)[1]     # number of trials per subject (=140)
numIter = 100             # number of iterations to find global minimum values
numPars = 5               # number of parameters

gain <- array(dim=c(N, T))
loss <- array(dim=c(N, T))
cert <- array(dim=c(N, T))
gamble <- array(dim=c(N, T))

for(i in 1:N) {
  tmpData = subset(dat, subjID==allSubjs[i])
  gain[i,] = tmpData$gain
  loss[i,] = abs(tmpData$loss)
  cert[i,] = tmpData$cert
  gamble[i,] = tmpData$gamble
}

dataList <- list(
  N       = N,  
  T       = T,
  Tsubj   = Tsubj,
  gain    = gain,
  loss    = loss,
  cert    = cert,
  gamble  = gamble
)

# run!
output = stan("ra_prospect_w_reparam(3).stan", data = dataList, pars = c("rho", "lambda", "tau", "mu_p", "sigma"),
              iter = 2000, warmup=1000, chains=4, cores=2)

# traceplot
traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
rho <- rstan::extract(output)$rho
lambda <- rstan::extract(output)$lambda
tau <- rstan::extract(output)$tau
mu_p1 <- rstan::extract(output)$mu_p[,1]
mu_p2 <- rstan::extract(output)$mu_p[,2]
mu_p3 <- rstan::extract(output)$mu_p[,3]
sigma1 <- rstan::extract(output)$sigma[,1]
sigma2 <- rstan::extract(output)$sigma[,2]
sigma3 <- rstan::extract(output)$sigma[,3]
rho1 <- rstan::extract(output)$rho[,1]
rho2 <- rstan::extract(output)$rho[,2]
rho3 <- rstan::extract(output)$rho[,3]
rho4 <- rstan::extract(output)$rho[,4]
rho5 <- rstan::extract(output)$rho[,5]
lambda1 <- rstan::extract(output)$lambda[,1]
lambda2 <- rstan::extract(output)$lambda[,2]
lambda3 <- rstan::extract(output)$lambda[,3]
lambda4 <- rstan::extract(output)$lambda[,4]
lambda5 <- rstan::extract(output)$lambda[,5]
tau1 <- rstan::extract(output)$tau[,1]
tau2 <- rstan::extract(output)$tau[,2]
tau3 <- rstan::extract(output)$tau[,3]
tau4 <- rstan::extract(output)$tau[,4]
tau5 <- rstan::extract(output)$tau[,5]

# plot posteriors
par(mfrow=c(3,1))
hist(mu_p1)
hist(mu_p2)
hist(mu_p3)
par(mfrow=c(3,1))
hist(sigma1)
hist(sigma2)
hist(sigma3)
par(mfrow=c(3,2))
hist(rho1)
hist(rho2)
hist(rho3)
hist(rho4)
hist(rho5)
par(mfrow=c(3,2))
hist(lambda1)
hist(lambda2)
hist(lambda3)
hist(lambda4)
hist(lambda5)
par(mfrow=c(3,2))
hist(tau1)
hist(tau2)
hist(tau3)
hist(tau4)
hist(tau5)

# 95% HDI of rho
HDIofMCMC(rho, credMass = 0.95)
HDIofMCMC(mu_p1, credMass = 0.95)
HDIofMCMC(sigma1, credMass = 0.95)
HDIofMCMC(rho1, credMass = 0.95)
HDIofMCMC(rho2, credMass = 0.95)
HDIofMCMC(rho3, credMass = 0.95)
HDIofMCMC(rho4, credMass = 0.95)
HDIofMCMC(rho5, credMass = 0.95)
HDIofMCMC(lambda, credMass = 0.95)
HDIofMCMC(mu_p2, credMass = 0.95)
HDIofMCMC(sigma2, credMass = 0.95)
HDIofMCMC(lambda1, credMass = 0.95)
HDIofMCMC(lambda2, credMass = 0.95)
HDIofMCMC(lambda3, credMass = 0.95)
HDIofMCMC(lambda4, credMass = 0.95)
HDIofMCMC(lambda5, credMass = 0.95)
HDIofMCMC(tau, credMass = 0.95)
HDIofMCMC(mu_p3, credMass = 0.95)
HDIofMCMC(sigma3, credMass = 0.95)
HDIofMCMC(tau1, credMass = 0.95)
HDIofMCMC(tau2, credMass = 0.95)
HDIofMCMC(tau3, credMass = 0.95)
HDIofMCMC(tau4, credMass = 0.95)
HDIofMCMC(tau5, credMass = 0.95)



