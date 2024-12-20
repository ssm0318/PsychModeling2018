#######################################################
## MLE_LSE_2.R                                       ##
## By JaeWon Kim, SNU                                ##
##                                                   ##
## Modified likelihood functions                     ##
## Code made more reusable                           ##
## Code Written on 03/23/2018                        ##
#######################################################

# minus log-likelihood of predicted probability
loglik <- function(p, n=n_total, x=n_corr) {
  # ensure 0 < p < 1
  p[p<=0] <- 10e-6
  p[p>=1] <- 1 - 10e-6
  
  # Calculate minus log-likelihood
  loglik <- (-1) * (x * log(p) + (n-x) * log(1-p))
  
  # return the summed minus log-likelihood as the function value
  sum(loglik)
}

# predicted probability by parameters
pow1 <- function(param, int=t_int, ll_flag='Y') {
  p <- (int+1) ^ (-param[1])
  if(ll_flag == 'Y') {
    return (loglik(p))
  }
  return (p)
}

pow2 <- function(param, int=t_int, ll_flag='Y') {
  p <- param[1] * (1+int) ^ (-param[2])
  if(ll_flag == 'Y') {
    return (loglik(p))
  }
  return (p)
}

exp1 <- function(param, int=t_int, ll_flag='Y') {
  p <- exp((-param[1]) * int)
  if(ll_flag == 'Y') {
    return (loglik(p))
  }
  return (p)
}

exp2 <- function(param, int=t_int, ll_flag='Y') {
  p <- param[1] * exp((-param[2])*int)
  if(ll_flag == 'Y') {
    return (loglik(p))
  }
  return (p)
}

expow <- function(param, int=t_int, ll_flag='Y') {
  p <- param[1] * exp((-param[2])*int) * (int+1) ^ (-param[3])
  if(ll_flag == 'Y') {
    return (loglik(p))
  }
  return (p)
}

hyp1 <- function(param, int=t_int, ll_flag='Y') {
  p <- 1 / (1+param[1]*int)
  if(ll_flag == 'Y') {
    return (loglik(p))
  }
  return (p)
}

hyp2 <- function(param, int=t_int, ll_flag='Y') {
  p <- param[1] / (1+param[2]*int)
  if(ll_flag == 'Y') {
    return (loglik(p))
  }
  return (p)
}




