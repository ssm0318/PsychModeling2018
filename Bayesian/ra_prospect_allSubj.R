# ra_prospect_stan_singleSubj.R
# Programmed by JaeWon Kim, May 2018

rm(list=ls())  # remove all variables 

library(rstan)
library(dplyr)

# set working directory where all the files exist
setwd("C:/Users/kleen/Dropbox/5-1/ExpPsych/HW/HW3")

# source HDIofMCMC.R to calculate HDI
source("HDIofMCMC.R") 

# read the data file
dat = read.table("ra_exampleData.txt", header=T, sep="\t")

allSubjs = unique(dat$subjID)  # all subject IDs
N = length(allSubjs)      # number of subjects
T = table(dat$subjID)[1]  # number of trials per subject (=140)
numIter = 100             # number of iterations to find global minimum values
numPars = 3               # number of parameters

# run!
bayes <- function(ID) {
  tmpData = subset(dat, subjID==ID)
  print(paste("Analyzing subject number", ID))
  dataList <- list(
    T       = T,
    gain    = tmpData$gain,
    loss    = abs(tmpData$loss),   # absolute value
    cert    = tmpData$cert,
    gamble  = tmpData$gamble
  )
  
  output = stan("ra_prospect_allSubj.stan", data = dataList, pars = c("rho", "lambda", "tau"),
                iter = 2000, warmup=1000, chains=2, cores=2)
  return (output)
}

results <- lapply(allSubjs, bayes)
res_all <- sflist2stanfit(results)
summary(res_all)

# compare MLE and Bayes. analysis results
bayes <- data.frame("method"="bayes", "rho" = c(1.020, 0.693, 0.794, 0.866, 0.962), "tau" = c(1.007, 4.742, 1.035, 3.130, 2.347), 
                    "lambda" = c(0.783, 2.475, 1.077, 0.982, 1.409))
mle <- data.frame("method"="mle", "rho" = c(1.02, 0.68, 0.79, 0.87, 0.96), "tau" = c(1.04, 4.55, 1.01, 3.13, 1.54),
                  "lambda" = c(0.79, 2.52, 1.12, 0.99, 1.38))

library(ggplot2)
ggplot() + geom_point(aes(mle$rho, bayes$rho, color="rho")) + geom_point(aes(mle$tau, bayes$tau, color="tau")) + 
  geom_point(aes(mle$lambda, bayes$lambda, color="lambda")) + geom_smooth(method="lm", se=FALSE, aes(mle$rho, bayes$rho, color="rho")) + 
  geom_smooth(method="lm", se=FALSE, aes(mle$tau, bayes$tau, color="tau")) + 
  geom_smooth(method="lm", se=FALSE, aes(mle$lambda, bayes$lambda, color="lambda")) + 
  labs(x="MLE", y="Bayesian", col="Parameters")


mod <- stan_model("ra_prospect_allSubj.stan")

posterior <- function(ID) {
  tmpData = subset(dat, subjID==ID)
  print(paste("Posterior checking subject number", ID))

  gamble <- as.numeric(select(tmpData, gamble)[,1])
  cert <- as.numeric(select(tmpData, cert)[,1])
  gain <- as.numeric(select(tmpData, gain)[,1])
  loss <- as.numeric(abs(select(tmpData, loss))[,1])
  
  fit <- sampling(mod, data=list(T, gamble, cert, gain, loss), iter = 2000, warmup = 1000)

  return (fit)
}

sample <- lapply(allSubjs, posterior)
post_all <- sflist2stanfit(sample)



posterior2 <- extract(res_all, 'y_pred')$y_pred
help(extract)

gamble <- as.numeric(select(dat, gamble)[,1])
cert <- select(dat, cert)
gain <- select(dat, gain)
loss <- select(dat, loss)
gamble[2]
str(gamble)
help(select)
str(loss)
str(T)
str(dat$gamble)
help(sampling)
head(dat)
str(gamble)
data(cars)
head(cars)
y <- select(cars, speed)
y
str(y)
??select
M <- nrow(cars)
str(M)
