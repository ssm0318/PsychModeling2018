#######################################################
## Parm_Estm_LSE.R                                   ##
## By JaeWon Kim, SNU                                ##
##                                                   ##
## Code Written on 03/23/2018                        ##
#######################################################

rm(list=ls())  # clear workspace
graphics.off() # close all figures
options(warn=-1) # turn off warnings

source("LSE.R")   # source LSE.R code

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
## LSE                  ##
##########################

# Try many different inits to escape from the local maxima
lse_iter <- function(num_par, f, low, up) {
  # Set non-interfering initial value of lse model
  lse_model <- c(NA)
  lse_model$value <- -10E10
  
  # Try many different inits to escape from the local maxima
  for(i in 1:100) {
    # Random inits generated on the spot for each iteration
    temp_lse <- optim(runif(num_par), f, method="L-BFGS-B", lower=low, upper=up)
    
    # Replace the results if the latest optimization yields better result
    if(temp_lse$value > lse_model$value) lse_model <- temp_lse
  }
  return (lse_model)
}

# Do the LSE
lse_model_pow1 <- lse_iter(1, lse_pow1, param_pow1_low, param_pow1_up)
lse_model_pow2 <- lse_iter(2, lse_pow2, param_pow2_low, param_pow2_up)
lse_model_exp1 <- lse_iter(1, lse_exp1, param_exp1_low, param_exp1_up)
lse_model_exp2 <- lse_iter(2, lse_exp2, param_exp2_low, param_exp2_up)
lse_model_expow <- lse_iter(3, lse_expow, param_expow_low, param_expow_up)
lse_model_hyp1 <- lse_iter(1, lse_hyp1, param_hyp1_low, param_hyp1_up)
lse_model_hyp2 <- lse_iter(2, lse_hyp2, param_hyp2_low, param_hyp2_up)

# Make the parameter lists of the same lengths to prevent unnecessary copies of the parameters
n <- max(length(lse_model_pow1$par), length(lse_model_pow2$par), length(lse_model_exp1$par), length(lse_model_exp2$par),
         length(lse_model_expow$par), length(lse_model_hyp1$par), length(lse_model_hyp2$par))
length(lse_model_pow1$par) <- n
length(lse_model_pow2$par) <- n
length(lse_model_exp1$par) <- n
length(lse_model_exp2$par) <- n
length(lse_model_expow$par) <- n
length(lse_model_hyp1$par) <- n
length(lse_model_hyp2$par) <- n

# Save the SSE value
sse_pow1 = lse_pow1(lse_model_pow1$par)
sse_pow2 = lse_pow2(lse_model_pow2$par)
sse_exp1 = lse_exp1(lse_model_exp1$par)
sse_exp2 = lse_exp2(lse_model_exp2$par)
sse_expow = lse_expow(lse_model_expow$par)
sse_hyp1 = lse_hyp1(lse_model_hyp1$par)
sse_hyp2 = lse_hyp2(lse_model_hyp2$par)

# Proportion of the explained variances for each model
r2 <- function(sse) {
  ssto <- 0
  for(i in p_corr) {
    ssto <- ssto + (i-mean(p_corr))^2
  }
  return (1-(sse/ssto))
}
r2_pow1 = r2(sse_pow1)
r2_pow2 = r2(sse_pow2)
r2_exp1 = r2(sse_exp1)
r2_exp2 = r2(sse_exp2)
r2_expow = r2(sse_expow)
r2_hyp1 = r2(sse_hyp1)
r2_hyp2 = r2(sse_hyp2)

# Generate summary
names = c("POW1", "POW2", "EXP1", "EXP2", "EXPOW", "HYP1", "HYP2")
sse = round(c(sse_pow1, sse_pow2, sse_exp1, sse_exp2, sse_expow, sse_hyp1, sse_hyp2), 4)
r2_lse <- round(c(r2_pow1, r2_pow2, r2_exp1, r2_exp2, r2_expow, r2_hyp1, r2_hyp2), 3)
pars_lse <- round(rbind(lse_model_pow1$par, lse_model_pow2$par, lse_model_exp1$par, lse_model_exp2$par,
                        lse_model_expow$par, lse_model_hyp1$par, lse_model_hyp2$par),3)
lse_summary = data.frame(Models = names, par = pars_lse, sse = sse, r2 = r2_lse)

# Plot the LSE results using ggplot (although there's nothing fancy about these graphs)
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
  stat_function(fun = pow1, linetype = "dashed", args=list(param=lse_model_pow1$par), aes(color='POW1', size=0.3)) +
  stat_function(fun = pow2, linetype = "dashed", args=list(param=lse_model_pow2$par), aes(color='POW2', size=0.3)) +
  stat_function(fun = exp1, linetype = "dashed", args=list(param=lse_model_exp1$par), aes(color='EXP1', size=0.3)) +
  stat_function(fun = exp2, linetype = "dashed", args=list(param=lse_model_exp2$par), aes(color='EXP2', size=0.3)) +
  stat_function(fun = expow, linetype = "dashed", args=list(param=lse_model_expow$par), aes(color='EXPOW', size=0.3)) +
  stat_function(fun = hyp1, linetype = "dashed", args=list(param=lse_model_hyp1$par), aes(color='HYP1', size=0.3)) +
  stat_function(fun = hyp2, linetype = "dashed", args=list(param=lse_model_hyp2$par), aes(color='HYP2', size=0.3)) +
  labs(title='LSE results', x="Time t", y="Proportion Correct") +
  theme(plot.title = element_text(hjust = 0.5)) + scale_color_hue("Model Type") + guides(size = FALSE)

# Print LSE results and best-fit parameters
print('- LSE results ------------')
print(lse_summary, 5)


