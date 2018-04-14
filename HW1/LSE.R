#######################################################
## LSE.R                                             ##
## By JaeWon Kim, SNU                                ##
##                                                   ##
## Functions Calculating LSE                         ##   
## Code Written on 03/23/2018                        ##
#######################################################

# Least squares estimate
lse_pow1 <- function(param, x=t_int, y=p_corr) {
  sum <- 0
  for (i in c(1:length(x))) {
    dev <- ((x[i]+1) ^ (-param[1])) - y[i]
    sum <- sum + dev^2
  }
  return (sum)
}

lse_pow2 <- function(param, x=t_int, y=p_corr) {
  sum <- 0
  for (i in c(1:length(x))) {
    dev <- (param[1]*((x[i]+1) ^ (-param[2]))) - y[i]
    sum <- sum + dev^2
  }
  return (sum)
}

lse_exp1 <- function(param, x=t_int, y=p_corr) {
  sum <- 0
  for (i in c(1:length(x))) {
    dev <- (exp((-param[1]) * x[i])) - y[i]
    sum <- sum + dev^2
  }
  return (sum)
}

lse_exp2 <- function(param, x=t_int, y=p_corr) {
  sum <- 0
  for (i in c(1:length(x))) {
    dev <- (param[1] * exp((-param[2])*x[i])) - y[i]
    sum <- sum + dev^2
  }
  return (sum)
}

lse_expow <- function(param, x=t_int, y=p_corr) {
  sum <- 0
  for (i in c(1:length(x))) {
    dev <- (param[1] * exp((-param[2])*x[i]) * (x[i]+1) ^ (-param[3])) - y[i]
    sum <- sum + dev^2
  }
  return (sum)
}

lse_hyp1 <- function(param, x=t_int, y=p_corr) {
  sum <- 0
  for (i in c(1:length(x))) {
    dev <- (1 / (1+param[1]*x[i])) - y[i]
    sum <- sum + dev^2
  }
  return (sum)
}

lse_hyp2 <- function(param, x=t_int, y=p_corr) {
  sum <- 0
  for (i in c(1:length(x))) {
    dev <- (param[1] / (1+param[2]*x[i])) - y[i]
    sum <- sum + dev^2
  }
  return (sum)
}


# predicted probability by parameters
pow1 <- function(param, int=t_int) {p <- (int+1) ^ (-param[1])}
pow2 <- function(param, int=t_int) {p <- param[1] * (1+int) ^ (-param[2])}
exp1 <- function(param, int=t_int) {p <- exp((-param[1]) * int)}
exp2 <- function(param, int=t_int) {p <- param[1] * exp((-param[2])*int)}
expow <- function(param, int=t_int) {p <- param[1] * exp((-param[2])*int) * (int+1) ^ (-param[3])}
hyp1 <- function(param, int=t_int) {p <- 1 / (1+param[1]*int)}
hyp2 <- function(param, int=t_int) {p <- param[1] / (1+param[2]*int)}

