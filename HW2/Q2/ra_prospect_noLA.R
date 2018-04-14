rm(list=ls())  # clear workspace
graphics.off() # close all figures
options(warn=-1) # turn off warnings

# for all subjects

mydata <- read.table("ra_exampleData.txt", header = TRUE)
subjID <- unique(mydata[, "subjID"])

N = 5  # number of subjects
Trials = 140 # number of trials per subject

param2_init <- runif(2)
param3_init <- runif(3)

param_ra_noLA_low <- c(0, 0); param_ra_noLA_up <- c(2, 5); #lower and upper bounds of ra_noLA model (0<p<2, 0<b<5)
param_ra_prospect_low <- c(0, 0, 0); param_ra_prospect_up <- c(2, 5, 10); #lower and upper bounds of ra model (0<p<2, 0<b<5, 0<d<10)

# risk aversion model with lambda
mle_ra_prospect <- function(param, df=mydata, T=Trials) {
  df <- df[which(df$subjID == subjID[loop]),]
  sum_minusLL = 0  # sum of minus log likelihood (across 140 trials). Initialize
  for (t in 1:T) {
    #evSafe: expected value of a certain (safe) option
    #evGamble: expected value of a risky option (gamble)
    #pGamble   # probability of choosing a gamble on each trial
    # free parameters: rho, tau, lambda
    
    evSafe   = df$cert[t]^param[1]
    evGamble = 0.5*(df$gain[t]^param[1] - param[3]*abs(df$loss[t])^param[1]) 
    pGamble  = 1 / (1 + exp(param[2]*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL = -log(pGamble)*df$gamble[t] - log(1-pGamble)*(1-df$gamble[t])  # LL of trial t
    sum_minusLL = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

# risk aversion model without lambda
mle_ra_noLA <- function(param, df=mydata, T=Trials) {
  df <- df[which(df$subjID == subjID[loop]),]
  sum_minusLL = 0  # sum of minus log likelihood (across 140 trials). Initialize
  for (t in 1:T) {
    #evSafe: expected value of a certain (safe) option
    #evGamble: expected value of a risky option (gamble)
    #pGamble   # probability of choosing a gamble on each trial
    # free parameters: rho, tau, lambda
    
    evSafe   = df$cert[t]^param[1]
    evGamble = 0.5*(df$gain[t]^param[1] - abs(df$loss[t])^param[1]) 
    pGamble  = 1 / (1 + exp(param[2]*(evSafe - evGamble)))
    pGamble  = pGamble * 0.9998 + 0.0001  # to make its range between 0.0001 and 0.9999
    tmp_minusLL   = -log(pGamble)*df$gamble[t] - log(1-pGamble)*(1-df$gamble[t])  # LL of trial t
    sum_minusLL   = sum_minusLL + tmp_minusLL
  }
  return(sum_minusLL)
}

# Try many different inits to escape from the local maxima
mle_iter <- function(num_par, f, low, up) {
  # Set non-interfering initial value of mle model
  mle_model <- c(NA)
  mle_model$value <- 10E10
  
  # Try many different inits to escape from the local maxima
  for(i in 1:100) {
    # Random inits generated on the spot for each iteration
    temp_mle <- optim(runif(num_par), f, method="L-BFGS-B", lower=low, upper=up)
    
    # Replace the results if the latest optimization yields better result
    if(temp_mle$value < mle_model$value) mle_model <- temp_mle
  }
  return (mle_model)
}

# Do the MLE
loop <- 1 # initialize loop
mle_model_ra_prospect <- c(NA) # empty row for matrix
for(i in 1:length(subjID)) {
  one_mle_model <- mle_iter(3, mle_ra_prospect, param_ra_prospect_low, param_ra_prospect_up)
  # save parameters and likelihood of current estimation
  one_mle_stat <- c("ra_prospect", subjID[loop], one_mle_model$par, one_mle_model$value)
  loop <- loop + 1 # loop through the next subject's data
  mle_model_ra_prospect <- rbind(mle_model_ra_prospect, one_mle_stat)
}
rownames(mle_model_ra_prospect) <- c()
mle_model_ra_prospect <- as.data.frame(mle_model_ra_prospect[-1,]) # remove empty row for matrix
colnames(mle_model_ra_prospect) <- c("model", "subjID", "rho", "tau", "lambda", "m_loglik")

loop <- 1 # initialize loop
mle_model_ra_noLA <- c(NA) # empty row for matrix
for(i in 1:length(subjID)) {
  one_mle_model <- mle_iter(2, mle_ra_noLA, param_ra_noLA_low, param_ra_noLA_up)
  # save parameters and likelihood of current estimation
  one_mle_stat <- c("ra_noLA", subjID[loop], one_mle_model$par, NA, one_mle_model$value)
  loop <- loop + 1 # loop through the next subject's data
  mle_model_ra_noLA <- rbind(mle_model_ra_noLA, one_mle_stat)
}
rownames(mle_model_ra_noLA) <- c()
mle_model_ra_noLA <- as.data.frame(mle_model_ra_noLA[-1,]) # remove empty row for matrix
colnames(mle_model_ra_noLA) <- c("model", "subjID", "rho", "tau", "lambda", "m_loglik")

# number of data points N
N = T

# Compute AIC and BIC, attach to data frame
AIC_BIC <- function(mle_summary, K) {
  # Define empty columns AIC and BIC
  mle_summary [, "AIC"] <- NA
  mle_summary [, "BIC"] <- NA
  # Loop through all given models
  for(n in 1:nrow(mle_summary)) {
    model_data <- mle_summary[n, ]
    # Compute AIC = -2*log(lik) + 2*K
    curr_AIC <- 2*as.numeric(as.character(model_data$m_loglik)) + 2*K
    # Compute BIC = -2*log(lik) + K*log(N)
    curr_BIC <- 2*as.numeric(as.character(model_data$m_loglik)) + K*log(Trials)
    mle_summary[n, "AIC"] <- round(curr_AIC, 3)
    mle_summary[n, "BIC"] <- round(curr_BIC, 3)
  }
  return (mle_summary)
}

mle_ra_prospect_summary <- AIC_BIC(mle_model_ra_prospect, 3)
mle_ra_noLA_summary <- AIC_BIC(mle_model_ra_noLA, 2)
mle_summary <- rbind(mle_ra_prospect_summary, mle_ra_noLA_summary)
mle_summary$rho <- round(as.numeric(as.character(mle_summary$rho)), 3)
mle_summary$tau <- round(as.numeric(as.character(mle_summary$tau)), 3)
mle_summary$lambda <- round(as.numeric(as.character(mle_summary$lambda)), 3)
mle_summary$m_loglik <- round(as.numeric(as.character(mle_summary$m_loglik)), 3)
mle_summary$d_AIC <- mle_summary$AIC - min(mle_summary$AIC) # calculate delta AIC
mle_summary$d_BIC <- mle_summary$BIC - min(mle_summary$BIC) #calculate delta BIC

summary(mle_summary$AIC[which(mle_summary$model == "ra_prospect")])
summary(mle_summary$AIC[which(mle_summary$model == "ra_noLA")])
summary(mle_summary$BIC[which(mle_summary$model == "ra_prospect")])
summary(mle_summary$BIC[which(mle_summary$model == "ra_noLA")])

print(mle_summary)





