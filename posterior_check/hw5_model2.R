# hw5_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr), May 2018

rm(list=ls())  # remove all variables

library(rstan)

# set working directory where all the files exist
setwd("~/Dropbox/5-1/ExpPsych/HW/HW5")

# read the data file
dat = read.table("simul_data_hw5_model1.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)[1]  # number of trials per subject 

choice  <- array(-1, c(N, T))
outcome <- array(0, c(N, T))

for (i in 1:N) {
  curSubj = allSubjs[i]
  tmp     = subset(dat, subjID == curSubj)
  choice[i, 1:T] <- tmp$choice
  outcome[i, 1:T] <- tmp$outcome
}

dataList <- list(
  N       = N,
  T       = T,
  Tsubj   = rep(T, N),
  choice  = choice,
  outcome = outcome
)

# run!
output = stan("hw5_model2.stan", data = dataList, pars = c("alpha", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)
saveRDS(output, file="output2.rds")

# traceplot
traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters <- rstan::extract(output)

alpha_mean = apply(parameters$alpha, 2, mean)
alpha_sd = apply(parameters$alpha, 2, sd)
beta_mean = apply(parameters$beta, 2, mean)
beta_sd = apply(parameters$beta, 2, sd)

# posterior
par(mfrow=c(5,6), mar = c(1,1,2,2) + 1)
for(i in 1:30) {
  hist(parameters$alpha[,i], xlab="", ylab="", main=paste("alpha",i))
}
for(i in 1:30) {
  hist(parameters$beta[,i], xlab="", ylab="", main=paste("beta",i))
}

alpha_beta <- readRDS("simul_pars1.rds")

library(Hmisc)

par(mfrow=c(1,1))
plot(alpha_beta$alpha, alpha_mean, xlim=c(0,0.4), ylim=c(0,0.5), pch=19)
errbar(alpha_beta$alpha, alpha_mean, alpha_mean+alpha_sd, alpha_mean-alpha_sd, add=TRUE)
abline(a=0,b=1)

plot(alpha_beta$beta, beta_mean, xlim=c(0.5,4), ylim=c(0.5,4))
errbar(alpha_beta$beta, beta_mean, beta_mean+beta_sd, beta_mean-beta_sd, add=TRUE)
abline(a=0,b=1)
