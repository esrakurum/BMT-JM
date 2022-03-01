model <- function ()
{
  betaX0[1:nt] <- BetaBs[, 1:ncF] %*% Bs.betaX0[1:ncF]
  betaX1[1:nt] <- BetaBs[, 1:ncF] %*% Bs.betaX1[1:ncF]
  betaX2[1:nt] <- BetaBs[, 1:ncF] %*% Bs.betaX2[1:ncF]
  
  betaZ1[1:nt] <- BetaBs[, 1:ncF] %*% Bs.betaZ1[1:ncF]
  betaZ2[1:nt] <- BetaBs[, 1:ncF] %*% Bs.betaZ2[1:ncF]
  
  
  for (i in 1:dim_data.id) 
  {
  ## Longitudinal Part
    for(j in (sum(nt.subj[1:i])+1): (sum(nt.subj[1:(i+1)])))
    {
      logit(prob[j]) <- betaX0[j - sum(nt.subj[1:i])] + betaX1[j - sum(nt.subj[1:i])] * x1[i] + betaX2[j - sum(nt.subj[1:i])] * x2[i] 
      + betaZ1[j - sum(nt.subj[1:i])] * z1[i]+ betaZ2[j - sum(nt.subj[1:i])] * z2[i] + b[i] + xi[fac.id[i]]
      y[j] ~ dbern(prob[j])
    }  
    
  ## Survival part
    
    beta.hazardX0[i] <- Bs.betaX0[1:ncF] %*% GammaBs[i, 1:ncF] 
    beta.hazardX1[i] <- Bs.betaX1[1:ncF] %*% GammaBs[i, 1:ncF]
    beta.hazardX2[i] <- Bs.betaX2[1:ncF] %*% GammaBs[i, 1:ncF]
    
    beta.hazardZ1[i] <- Bs.betaZ1[1:ncF] %*% GammaBs[i, 1:ncF]
    beta.hazardZ2[i] <- Bs.betaZ2[1:ncF] %*% GammaBs[i, 1:ncF]
    
    gamma.hazardX1[i] <- Bs.gammaX1[1:ncF] %*% GammaBs[i, 1:ncF]
    gamma.hazardX2[i] <- Bs.gammaX2[1:ncF] %*% GammaBs[i, 1:ncF]
    # 
    gamma.hazardZ1[i] <- Bs.gammaZ1[1:ncF] %*% GammaBs[i, 1:ncF]
    gamma.hazardZ2[i] <- Bs.gammaZ2[1:ncF] %*% GammaBs[i, 1:ncF]
    
    prob_hazard[i] <- 1/(1+exp(-(beta.hazardX0[i] + beta.hazardX1[i]*x1[i] + beta.hazardX2[i]*x2[i] + beta.hazardZ1[i]*z1[i] 
                                 + beta.hazardZ2[i]*z2[i] + b[i] + xi[fac.id[i]])))
    
    log.h0.T[i] <-  Bs.h0[1:ncF] %*% GammaBs[i, 1:ncF]
    f.long[i] <- Bs.alpha[1:ncF] %*% GammaBs[i, 1:ncF]
    
    haz[i] <- exp(log.h0.T[i] + f.long[i] * prob_hazard[i] + gamma.hazardX1[i]*x2[i] + gamma.hazardX2[i]*x1[i] 
                  + gamma.hazardZ1[i]*z1[i] + gamma.hazardZ2[i]*z2[i])
    
    for(k in 1:K)
    {
      
      beta.survX0[i, k] <- Bs.betaX0[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF] 
      beta.survX1[i, k] <- Bs.betaX1[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF]
      beta.survX2[i, k] <- Bs.betaX2[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF]
      
      beta.survZ1[i, k] <- Bs.betaZ1[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF]
      beta.survZ2[i, k] <- Bs.betaZ2[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF]
      
      gamma.survX1[i, k] <- Bs.gammaX1[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF] 
      gamma.survX2[i, k] <- Bs.gammaX2[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF] 
      
      gamma.survZ1[i, k] <- Bs.gammaZ1[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF] 
      gamma.survZ2[i, k] <- Bs.gammaZ2[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF] 
      
      prob_surv[i, k] <- 1/(1+exp(-(beta.survX0[i, k] + beta.survX1[i, k]*x1[i] + beta.survX2[i, k]*x2[i] 
                                    + beta.survZ1[i, k]*z1[i] + beta.survZ2[i, k]*z2[i] + b[i] + xi[fac.id[i]])))
      
      log.h0.s[i, k] <- Bs.h0[1:ncF] %*% GammaBs.surv[K * (i - 1) + k, 1:ncF]
      f.longs[i, k] <- Bs.alpha[1:ncF]%*% GammaBs.surv[K * (i - 1) + k, 1:ncF]
      
      surv[i, k] <- exp(log.h0.s[i, k] + f.longs[i, k] * prob_surv[i, k]  + gamma.survX1[i, k]*x2[i] + gamma.survX2[i, k]*x1[i] 
                        + gamma.survZ1[i, k]*z1[i] + gamma.survZ2[i, k]*z2[i])
      
  
    }
    log.survival[i] <-  -T1[i]/2*inprod(wk, surv[i ,])  
    phi[i] <- 100000 - (event[i] * log(haz[i])) - log.survival[i] 
    zeros[i] ~ dpois(phi[i]) 
  }
  
  # Priors
  for (l in 1:dim_data.id)
  {
    b[l] ~ dnorm(mu0, 1/sigmasq_S)
  }
  for(k in 1:n)
  {
    xi[k] ~ dnorm(mu1, 1/sigmasq_F)
  }
  
  
  Bs.betaX0[1:ncF] ~ dmnorm(priorMean.Bs.betaX0[], Tau.betaX0 * priorTau.betaX0[, ])
  Bs.betaX1[1:ncF] ~ dmnorm(priorMean.Bs.betaX1[], Tau.betaX1 * priorTau.betaX1[, ])
  Bs.betaX2[1:ncF] ~ dmnorm(priorMean.Bs.betaX2[], Tau.betaX2 * priorTau.betaX2[, ])
  
  Bs.betaZ1[1:ncF] ~ dmnorm(priorMean.Bs.betaZ1[], Tau.betaZ1 * priorTau.betaZ1[, ])
  Bs.betaZ2[1:ncF] ~ dmnorm(priorMean.Bs.betaZ2[], Tau.betaZ2 * priorTau.betaZ2[, ])
  
  
  Tau.betaX0 ~ dgamma(1, 0.005)
  Tau.betaX1 ~ dgamma(1, 0.005)
  Tau.betaX2 ~ dgamma(1, 0.005)
  
  Tau.betaZ1 ~ dgamma(1, 0.005)
  Tau.betaZ2 ~ dgamma(1, 0.005)
  
  sigmasq_S ~ dgamma(priorA, priorB)
  
  sigmasq_F ~ dgamma(priorA, priorB)
  
  Bs.gammaX1[1:ncF] ~ dmnorm(priorMean.Bs.gammaX1[], Tau.gammaX1 * priorTau.gammaX1[, ])
  Bs.gammaX2[1:ncF] ~ dmnorm(priorMean.Bs.gammaX2[], Tau.gammaX2 * priorTau.gammaX2[, ])
   
  Tau.gammaX1 ~ dgamma(1, 0.005)
  Tau.gammaX2 ~ dgamma(1, 0.005)
   
  Bs.gammaZ1[1:ncF] ~ dmnorm(priorMean.Bs.gammaZ1[], Tau.gammaZ1 * priorTau.gammaZ1[, ])
  Bs.gammaZ2[1:ncF] ~ dmnorm(priorMean.Bs.gammaZ2[], Tau.gammaZ2 * priorTau.gammaZ2[, ])
   
  Tau.gammaZ1 ~ dgamma(1, 0.005)
  Tau.gammaZ2 ~ dgamma(1, 0.005)
  
  Bs.alpha[1:ncF] ~ dmnorm(priorMean.Bs.alpha[],  Tau.alphas * priorTau.alpha[, ])
  
  Tau.alphas ~ dgamma(0.01, 0.5) ##
  
  Bs.h0[1:ncF] ~ dmnorm(priorMean.Bs.h0[], Tau.Bs.h0 * priorTau.Bs.h0[,])
  Tau.Bs.h0 ~ dgamma(0.1, 0.005)
  
}
