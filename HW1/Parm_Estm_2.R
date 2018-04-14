#######################################################
## Parm_estm_2.R                                     ##
## By JaeWon Kim, SNU                                ##
##                                                   ##
## Parameter Estimation using MLE method             ##
## Code made more reusable                           ##
## Code Written on 03/23/2018                        ##
#######################################################

# Goals of Modification: minimizing rewrite of functions / unnecessary variables

rm(list=ls())  # clear workspace
graphics.off() # close all figures
options(warn=-1) # turn off warnings

source("MLE_LSE_2.R")   # source MLE_LSE.R code

##########################
## General Setup        ##
## Data and Parameters  ##
##########################
n_total <- 50 # sample size
t_int <- c(0.5, 1, 2, 4, 8, 12, 16, 18) # time interval values
n_corr <- c(44, 34, 27, 26, 19, 17, 20, 11) # number of correct responses
p_corr <- n_corr/n_total # proportion correct

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
mle_model_pow1 <- mle_iter(1, pow1, param_pow1_low, param_pow1_up)
mle_model_pow2 <- mle_iter(2, pow2, param_pow2_low, param_pow2_up)
mle_model_exp1 <- mle_iter(1, exp1, param_exp1_low, param_exp1_up)
mle_model_exp2 <- mle_iter(2, exp2, param_exp2_low, param_exp2_up)
mle_model_expow <- mle_iter(3, expow, param_expow_low, param_expow_up)
mle_model_hyp1 <- mle_iter(1, hyp1, param_hyp1_low, param_hyp1_up)
mle_model_hyp2 <- mle_iter(2, hyp2, param_hyp2_low, param_hyp2_up)

# Make the parameter lists of the same lengths to prevent unnecessary copies of the parameters
n <- max(length(mle_model_pow1$par), length(mle_model_pow2$par), length(mle_model_exp1$par), length(mle_model_exp2$par),
         length(mle_model_expow$par), length(mle_model_hyp1$par), length(mle_model_hyp2$par))
length(mle_model_pow1$par) <- n
length(mle_model_pow2$par) <- n
length(mle_model_exp1$par) <- n
length(mle_model_exp2$par) <- n
length(mle_model_expow$par) <- n
length(mle_model_hyp1$par) <- n
length(mle_model_hyp2$par) <- n

# Save the MLE parameter estimates
names = c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2")
mle_model_par <- as.data.frame(rbind(mle_model_pow1$par, mle_model_pow2$par, mle_model_exp1$par, mle_model_exp2$par,
                   mle_model_expow$par, mle_model_hyp1$par, mle_model_hyp2$par))
rownames(mle_model_par) <- names

# Proportion of the explained variances for each model
r2 <- function(p_prd) {1-sum((p_corr-p_prd)^2)/sum((p_corr-mean(p_corr))^2)}
r2_pow1 = r2(pow1(as.numeric(mle_model_par['POW1',]), ll_flag='N'))
r2_pow2 = r2(pow2(as.numeric(mle_model_par['POW2',]), ll_flag='N'))
r2_exp1 = r2(exp1(as.numeric(mle_model_par['EXP1',]), ll_flag='N'))
r2_exp2 = r2(exp2(as.numeric(mle_model_par['EXP2',]), ll_flag='N'))
r2_expow = r2(expow(as.numeric(mle_model_par['EXPOW',]), ll_flag='N'))
r2_hyp1 = r2(hyp1(as.numeric(mle_model_par['HYP1',]), ll_flag='N'))
r2_hyp2 = r2(hyp2(as.numeric(mle_model_par['HYP2',]), ll_flag='N'))

# Generate summary
minus_loglik_MLE = round(c(mle_model_pow1$value, mle_model_pow2$value, mle_model_exp1$value, mle_model_exp2$value,
                           mle_model_expow$value, mle_model_hyp1$value, mle_model_hyp2$value), 3)
r2_mle <- round(c(r2_pow1, r2_pow2, r2_exp1, r2_exp2, r2_expow, r2_hyp1, r2_hyp2), 3)
pars_mle <- round(mle_model_par,3)
mle_summary = data.frame(Models = names, par = pars_mle, loglik = - minus_loglik_MLE, r2 = r2_mle)

# Plot the MLE results using ggplot (although there's nothing fancy about these graphs)
#install.packages("ggplot2")
library(ggplot2)
library(RColorBrewer)

# Generate data frame of data points
mydata <- data.frame(t_int, p_corr)

# Plot data points
p <- ggplot(data = data.frame(x = c(0,20)), mapping = aes(x = x)) +
  geom_point(data=mydata, aes(x=t_int, y=p_corr, size=3))

# Plot each model using stat_function
p + 
  stat_function(fun = pow1, linetype = "dashed", args=list(param=as.numeric(mle_model_par['POW1',]), ll_flag='N'), aes(color='POW1')) +
  stat_function(fun = pow2, linetype = "dashed", args=list(param=as.numeric(mle_model_par['POW2',]), ll_flag='N'), aes(color='POW2')) +
  stat_function(fun = exp1, linetype = "dashed", args=list(param=as.numeric(mle_model_par['EXP1',]), ll_flag='N'), aes(color='EXP1')) +
  stat_function(fun = exp2, linetype = "dashed", args=list(param=as.numeric(mle_model_par['EXP2',]), ll_flag='N'), aes(color='EXP2')) +
  stat_function(fun = expow, linetype = "dashed", args=list(param=as.numeric(mle_model_par['EXPOW',]), ll_flag='N'), aes(color='EXPOW')) +
  stat_function(fun = hyp1, linetype = "dashed", args=list(param=as.numeric(mle_model_par['HYP1',]), ll_flag='N'), aes(color='HYP1')) +
  stat_function(fun = hyp2, linetype = "dashed", args=list(param=as.numeric(mle_model_par['HYP2',]), ll_flag='N'), aes(color='HYP2')) +
  labs(title='MLE results', x="Time t", y="Proportion Correct") +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_hue("Model Type") + guides(size = FALSE)
 
# Print MLE results and best-fit parameters
print('- MLE results ------------')
print(mle_summary, 5)


