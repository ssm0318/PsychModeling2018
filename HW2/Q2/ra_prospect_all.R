rm(list=ls())  # clear workspace
graphics.off() # close all figures
options(warn=-1) # turn off warnings

# for all subjects

mydata <- read.table("ra_exampleData.txt", header = TRUE)
subjID <- unique(mydata[, "subjID"])
loop <- 1

N = 5  # number of subjects
Trials = 140 # number of trials per subject

param3_init <- runif(3)

param_ra_prospect_low <- c(0, 0, 0); param_ra_prospect_up <- c(2, 5, 10); #lower and upper bounds of ra model (0<p<2, 0<b<5, 0<d<10)

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
mle_model_ra_prospect <- c(NA) # empty row for matrix
for(i in 1:length(subjID)) {
  one_mle_model <- mle_iter(3, mle_ra_prospect, param_ra_prospect_low, param_ra_prospect_up)
  # save parameters and likelihood of current estimation
  one_mle_stat <- c(one_mle_model$par, one_mle_model$value)
  loop <- loop + 1 # loop through the next subject's data
  mle_model_ra_prospect <- rbind(mle_model_ra_prospect, one_mle_stat)
}
mle_model_ra_prospect <- mle_model_ra_prospect[-1,] # remove empty row for matrix

# number of parameters K
K <- c(3)

# number of data points N
N = T

# Generate data frame
names = c("ra_prospect_1", "ra_prospect_2", "ra_prospect_3", "ra_prospect_4", "ra_prospect_5")
rownames(mle_model_ra_prospect) <- names
colnames(mle_model_ra_prospect) <- c("rho", "tau", "lambda", "-loglik")
mle_model_ra_prospect <- as.data.frame((mle_model_ra_prospect))

# Compute AIC and BIC, attach to data frame
AIC_BIC <- function(mle_summary) {
  # Define empty columns AIC and BIC
  mle_summary [, "AIC"] <- NA
  mle_summary [, "BIC"] <- NA
  # Loop through all given models
  for(n in 1:nrow(mle_summary)) {
    model_data <- mle_summary[n, ]
    # Compute AIC = -2*log(lik) + 2*K
    curr_AIC <- as.numeric(2*model_data["-loglik"] + 2*K[1])
    # Compute BIC = -2*log(lik) + K*log(N)
    curr_BIC <- as.numeric(2*model_data["-loglik"] + K[1]*log(Trials))
    mle_summary[n, "AIC"] <- round(curr_AIC, 3)
    mle_summary[n, "BIC"] <- round(curr_BIC, 3)
  }
  return (mle_summary)
}

mle_summary <- AIC_BIC(mle_model_ra_prospect)

mle_summary <- round(mle_summary, 3)

print(mle_summary)





