#######################################################
## Parm_Estm.R                                       ##
## PSYCH 7695, Spring 2017                           ##
## Maximum Likelihood Estimation (Myung, JMP, 2003)  ##
## By Yun Tang, Psychology, OSU                      ##
##                                                   ##
## Main Program                                      ##
## Code Written on 12/18/2009                        ##
##                                                   ##
## Modified by Joonsuk Park on Jan 28 2015           ##
## Modified by Jay Myung in Feb 2017                 ##
## Modified by Woo-Young Ahn in March 2018           ##
## Modified by JaeWon Kim in March 2018              ##
#######################################################

# Loading the (minus) log-likelihood functions
# Please modify the path according to the actual location of the file "MLE_LSE.R"
# e.g., setwd("/Users/youngahn/this-course/")

rm(list=ls())  # clear workspace
graphics.off() # close all figures

source("MLE_LSE.R")   # source MLE_LSE.R code

##########################
## General Setup        ##
## Data and Parameters  ##
##########################
n_total <- 50 # sample size
t_int <- c(0.5, 1, 2, 4, 8, 12, 16, 18) # time interval values
n_corr <- c(44, 34, 27, 26, 19, 17, 20, 11) # number of correct responses
p_corr <- n_corr/n_total # proportion correct

# Generate random uniform numbers between 0 and 1 to use as initials for the optim procedure
param1_init <- runif(1)
param2_init <- runif(2)
param3_init <- runif(3)

param_pow1_low <- c(0); param_pow1_up <- c(3); # lower and upper bounds of POW1 model (0<b<3)
param_pow2_low <- c(0, 0); param_pow2_up <- c(1, 3);  # lower and upper bounds of POW2 model (0<a<1, 0<b<3)
param_exp1_low <- c(0); param_exp1_up <- c(3);  # lower and upper bounds of EXP1 model (0<b<3)
param_exp2_low <- c(0, 0); param_exp2_up <- c(1, 3);  # lower and upper bounds of EXP2 model (0<a<1, 0<b<3)
param_expow_low <- c(0, 0, -Inf); param_expow_up <- c(1, Inf, 3);  # lower and upper bounds of EXPOW model (0<a<1, 0<b, c<3)
param_hyp1_low <- c(0); param_hyp1_up <- c(1);  # lower and upper bounds of HYP1 model (0<b<1)
param_hyp2_low <- c(0, 0); param_hyp2_up <- c(1, 1);  # lower and upper bounds of HYP2 model (0<a<1, 0<b<1)

##########################
## MLE                  ##
##########################

# Call general purpose optimization rountine
mle_model_pow1 <- optim(param1_init, mle_pow1, method="L-BFGS-B", lower=param_pow1_low, upper=param_pow1_up, int=t_int, n=n_total, x=n_corr)
mle_model_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
mle_model_exp1 <- optim(param1_init, mle_exp1, method="L-BFGS-B", lower=param_exp1_low, upper=param_exp1_up, int=t_int, n=n_total, x=n_corr)
mle_model_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)
mle_model_expow <- optim(param3_init, mle_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, n=n_total, x=n_corr)
mle_model_hyp1 <- optim(param1_init, mle_hyp1, method="L-BFGS-B", lower=param_hyp1_low, upper=param_hyp1_up, int=t_int, n=n_total, x=n_corr)
mle_model_hyp2 <- optim(param2_init, mle_hyp2, method="L-BFGS-B", lower=param_hyp2_low, upper=param_hyp2_up, int=t_int, n=n_total, x=n_corr)

# Try many different inits to escape from the local maxima
for(i in 1:100) {
  # Re-generate random inits. Is it the best way to do this?
  param1_init <- runif(1); param2_init <- runif(2); param3_init <- runif(3); 
  
  # Do the MLE again
  temp_pow1 <- optim(param1_init, mle_pow1, method="L-BFGS-B", lower=param_pow1_low, upper=param_pow1_up, int=t_int, n=n_total, x=n_corr)
  temp_pow2 <- optim(param2_init, mle_pow2, method="L-BFGS-B", lower=param_pow2_low, upper=param_pow2_up, int=t_int, n=n_total, x=n_corr)
  temp_exp1 <- optim(param1_init, mle_exp1, method="L-BFGS-B", lower=param_exp1_low, upper=param_exp1_up, int=t_int, n=n_total, x=n_corr)
  temp_exp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_exp2_low, upper=param_exp2_up, int=t_int, n=n_total, x=n_corr)
  temp_expow <- optim(param3_init, mle_expow, method="L-BFGS-B", lower=param_expow_low, upper=param_expow_up, int=t_int, n=n_total, x=n_corr)
  temp_hyp1 <- optim(param1_init, mle_hyp1, method="L-BFGS-B", lower=param_hyp1_low, upper=param_hyp1_up, int=t_int, n=n_total, x=n_corr)
  temp_hyp2 <- optim(param2_init, mle_exp2, method="L-BFGS-B", lower=param_hyp2_low, upper=param_hyp2_up, int=t_int, n=n_total, x=n_corr)

  # Replace the results if the latest optimization yields better result
  if(temp_pow1$value < mle_model_pow1$value) mle_model_pow1 <- temp_pow1
  if(temp_pow2$value < mle_model_pow2$value) mle_model_pow2 <- temp_pow2  
  if(temp_exp1$value < mle_model_exp1$value) mle_model_exp1 <- temp_exp1
  if(temp_exp2$value < mle_model_exp2$value) mle_model_exp2 <- temp_exp2
  if(temp_expow$value < mle_model_expow$value) mle_model_expow <- temp_expow
  if(temp_hyp1$value < mle_model_hyp1$value) mle_model_hyp1 <- temp_hyp1
  if(temp_hyp2$value < mle_model_hyp2$value) mle_model_hyp2 <- temp_hyp2
}

# Save the MLE parameter estimates
parm_pow1 <- mle_model_pow1$par
parm_pow2 <- mle_model_pow2$par
parm_exp1 <- mle_model_exp1$par
parm_exp2 <- mle_model_exp2$par
parm_expow <- mle_model_expow$par
parm_hyp1 <- mle_model_hyp1$par
parm_hyp2 <- mle_model_hyp2$par

# MLE predictions
p_prd_pow1 <- (t_int+1) ^ (-parm_pow1[1])
p_prd_pow2 <- parm_pow2[1]*(t_int+1)^(-parm_pow2[2])
p_prd_exp1 <- exp((-parm_exp1[1]) * t_int)
p_prd_exp2 <- parm_exp2[1]*exp(-parm_exp2[2]*t_int)
p_prd_expow <- parm_expow[1] * exp((-parm_expow[2])*t_int) * (t_int+1) ^ (-parm_expow[3])
p_prd_hyp1 <- 1 / (1+parm_hyp1[1]*t_int)
p_prd_hyp2 <- parm_hyp2[1] / (1+parm_hyp2[2]*t_int)

# Proportion of the explained variances for each model
r2_pow1 = 1-sum((p_corr-p_prd_pow1)^2)/sum((p_corr-mean(p_corr))^2)
r2_pow2 = 1-sum((p_corr-p_prd_pow2)^2)/sum((p_corr-mean(p_corr))^2)
r2_exp1 = 1-sum((p_corr-p_prd_exp1)^2)/sum((p_corr-mean(p_corr))^2)
r2_exp2 = 1-sum((p_corr-p_prd_exp2)^2)/sum((p_corr-mean(p_corr))^2)
r2_expow = 1-sum((p_corr-p_prd_expow)^2)/sum((p_corr-mean(p_corr))^2)
r2_hyp1 = 1-sum((p_corr-p_prd_hyp1)^2)/sum((p_corr-mean(p_corr))^2)
r2_hyp2 = 1-sum((p_corr-p_prd_hyp2)^2)/sum((p_corr-mean(p_corr))^2)

# Generate summary
minus_loglik_MLE = round(c(mle_model_pow1$value, mle_model_pow2$value, mle_model_exp1$value, mle_model_exp2$value,
                           mle_model_expow$value, mle_model_hyp1$value, mle_model_hyp2$value), 3)
r2_mle <- round(c(r2_pow1, r2_pow2, r2_exp1, r2_exp2, r2_expow, r2_hyp1, r2_hyp2), 3)
names = c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2")
pars_mle <- round(cbind(mle_model_pow1$par, mle_model_pow2$par, mle_model_exp1$par, mle_model_exp2$par,
                        mle_model_expow$par, mle_model_hyp1$par, mle_model_hyp2$par),3)
dimnames(pars_mle) = list(c('par1', 'par2', 'par3'),c('POW1', 'POW2', 'EXP1', 'EXP2', 'EXPOW', 'HYP1', 'HYP2'))

mle_summary = data.frame(Models = names, loglik = - minus_loglik_MLE, r2 = r2_mle)

# Plot the MLE results using simple R graphics. Use ggplots for fancy graphs
x <- seq(0,20, 0.05)
p_pow1 <- (x+1) ^ (-parm_pow1[1])
p_pow2 <- parm_pow2[1]*(x+1)^(-parm_pow2[2])
p_exp1 <- exp((-parm_exp1[1]) * x)
p_exp2 <- parm_exp2[1]*exp(-parm_exp2[2]*x)
p_expow <- parm_expow[1] * exp((-parm_expow[2])*x) * (x+1) ^ (-parm_expow[3])
p_hyp1 <- 1 / (1+parm_hyp1[1]*x)
p_hyp2 <- parm_hyp2[1] / (1+parm_hyp2[2]*x)
plot(x, p_pow1, ylim=c(0,1), xlab='Time t', ylab='Proportion Correct', 
     main="MLE results", type='l', lwd=2, lty='dashed')
lines(x, p_pow2, lwd=2, lty='dashed', col='red')
lines(x, p_exp1, lwd=2, lty='dashed', col='orange')
lines(x, p_exp2, lwd=2, lty='dashed', col='green')
lines(x, p_expow, lwd=2, lty='dashed', col='blue')
lines(x, p_hyp1, lwd=2, lty='dashed', col='purple')
lines(x, p_hyp2, lwd=2, lty='dashed', col='pink')
points(t_int, p_corr, pch=19, cex=1.5)
text = c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2")
ltys = c('dashed')
cols = c('black','red', 'orange', 'green', 'blue', 'purple', 'pink')
legend(14, 1, text, lwd=2, cex=1, lty=ltys, col=cols)

# print maximized likehood values
print('- MLE results ------------')
print(mle_summary,4)

# print bet-fit parameter values
print('- Best-fit parameters --------')
print(pars_mle,4)

