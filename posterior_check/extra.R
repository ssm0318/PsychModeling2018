# hw5_model1_exec.R
# Programmed by Woo-Young Ahn (wahn55@snu.ac.kr), May 2018

rm(list=ls())  # remove all variables

library(rstan)

# set working directory where all the files exist
setwd("~/Dropbox/5-1/ExpPsych/HW/HW5")

# read the data file
dat = read.table("simul_data_hw5_model4.txt", header=T, sep="\t")

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
output = stan("extra.stan", data = dataList, pars = c("mu_alpha_pos", "mu_alpha_neg", "mu_beta", "alpha_pos", "alpha_neg", "beta"),
              iter = 2000, warmup=1000, chains=2, cores=2)
saveRDS(output, file="output5.rds")
output4 <- readRDS("output4.rds")
output5 <- readRDS("output5.rds")

# traceplot
traceplot(output)

# print summary
print(output)

# extract Stan fit object (parameters)
parameters4 <- rstan::extract(output4)
parameters5 <- rstan::extract(output5)

alpha_pos_mean4 = apply(parameters4$alpha_pos, 2, mean)
alpha_neg_mean4 = apply(parameters4$alpha_neg, 2, mean)
alpha_pos_sd4 = apply(parameters4$alpha_pos, 2, sd)
alpha_neg_sd4 = apply(parameters4$alpha_neg, 2, sd)
beta_mean4 = apply(parameters4$beta, 2, mean)
beta_sd4 = apply(parameters4$beta, 2, sd)
mu_alpha_pos4 <- mean(parameters4$mu_alpha_pos)
mu_alpha_neg4 <- mean(parameters4$mu_alpha_neg)
mu_beta4 <- mean(parameters4$mu_beta)

alpha_pos_mean5 = apply(parameters5$alpha_pos, 2, mean)
alpha_neg_mean5 = apply(parameters5$alpha_neg, 2, mean)
alpha_pos_sd5 = apply(parameters5$alpha_pos, 2, sd)
alpha_neg_sd5 = apply(parameters5$alpha_neg, 2, sd)
beta_mean5 = apply(parameters5$beta, 2, mean)
beta_sd5 = apply(parameters5$beta, 2, sd)
mu_alpha_pos5 <- mean(parameters5$mu_alpha_pos)
mu_alpha_neg5 <- mean(parameters5$mu_alpha_neg)
mu_beta5 <- mean(parameters5$mu_beta)

# posterior
hist(parameters4$mu_alpha_pos, breaks=20)
hist(parameters4$mu_alpha_neg, breaks=20)
hist(parameters4$mu_beta, breaks=20)

hist(parameters5$mu_alpha_pos, breaks=20)
hist(parameters5$mu_alpha_neg, breaks=20)
hist(parameters5$mu_beta, breaks=20)

par(mfrow=c(5,6), mar = c(1,1,2,2) + 1)
for(i in 1:30) {
  hist(parameters4$alpha_pos[,i], xlab="", ylab="", main=paste("pos",i))
}
for(i in 1:30) {
  hist(parameters4$alpha_neg[,i], xlab="", ylab="", main=paste("neg",i))
}
for(i in 1:30) {
  hist(parameters4$beta[,i], xlab="", ylab="", main=paste("beta",i))
}

par(mfrow=c(5,6), mar = c(1,1,2,2) + 1)
for(i in 1:30) {
  hist(parameters5$alpha_pos[,i], xlab="", ylab="", main=paste("pos",i))
}
for(i in 1:30) {
  hist(parameters5$alpha_neg[,i], xlab="", ylab="", main=paste("neg",i))
}
for(i in 1:30) {
  hist(parameters5$beta[,i], xlab="", ylab="", main=paste("beta",i))
}

alpha_beta4 <- readRDS("simul_pars3.rds")
alpha_beta5 <- readRDS("simul_pars4.rds")

library(Hmisc)

par(mfrow=c(1,1))
plot(alpha_beta4$alpha_pos, alpha_pos_mean4, xlim=c(0,0.4), ylim=c(0,0.6), pch=19)
errbar(alpha_beta4$alpha_pos, alpha_pos_mean4, alpha_pos_mean4+alpha_pos_sd4, alpha_pos_mean4-alpha_pos_sd4, add=TRUE)
abline(a=0,b=1)
abline(h=mu_alpha_pos4)
par(new=TRUE)
plot(alpha_beta5$alpha_pos, alpha_pos_mean5, xlim=c(0,0.4), ylim=c(0,0.6), pch=19, col="red")
errbar(alpha_beta5$alpha_pos, alpha_pos_mean5, alpha_pos_mean5+alpha_pos_sd5, alpha_pos_mean5-alpha_pos_sd5, add=TRUE, col="red")
abline(a=0,b=1)
abline(h=mu_alpha_pos5)

par(mfrow=c(1,1))
plot(alpha_beta4$alpha_neg, alpha_neg_mean4, xlim=c(0,0.4), ylim=c(0.25,0.5), pch=19)
errbar(alpha_beta4$alpha_neg, alpha_neg_mean4, alpha_neg_mean4+alpha_neg_sd4, alpha_neg_mean4-alpha_neg_sd4, add=TRUE)
abline(a=0,b=1)
abline(h=mu_alpha_neg4)
par(new=TRUE)
plot(alpha_beta5$alpha_neg, alpha_neg_mean5, xlim=c(0,0.4), ylim=c(0.25,0.5), pch=19, col="red")
errbar(alpha_beta5$alpha_neg, alpha_neg_mean5, alpha_neg_mean5+alpha_neg_sd5, alpha_neg_mean5-alpha_neg_sd5, add=TRUE, col="red")
abline(a=0,b=1)
abline(h=mu_alpha_neg5)

par(mfrow=c(1,1))
plot(alpha_beta4$beta, beta_mean4, xlim=c(0,4), ylim=c(0,4), pch=19)
errbar(alpha_beta4$beta, beta_mean4, beta_mean4+beta_sd4, beta_mean4-beta_sd4, add=TRUE)
abline(a=0,b=1)
abline(h=mu_beta4)
par(new=TRUE)
plot(alpha_beta5$beta, beta_mean5, xlim=c(0,4), ylim=c(0,4), pch=19, col="red")
errbar(alpha_beta5$beta, beta_mean5, beta_mean5+beta_sd5, beta_mean5-beta_sd5, add=TRUE, col="red")
abline(a=0,b=1)
abline(h=mu_beta5)
