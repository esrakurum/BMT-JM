
####Install packages

library(rjags)
library(R2WinBUGS)
library(parallel)

###Set the directory
setwd("....") 


###Load the data
load("sample_data.rdata")


### BMTJM - Main function to obtain the posterior samples
### n.iter = number of iterations (including burn-in)
### n.burn = number of iterations during burn-in 
### n.thin = thinning
### n.chains = number of chains that will be run in parallel

BMTJM_fit = function(data, n.iter, n.burn, n.thin, n.chains)
{  
################################################
# Prepare data:
  dat = data
  
  myEnv <- environment()
  source("prepare_data.R", local = myEnv)

################################################
# Create the txt file for the jags model

  source("jags.R")
  filename <- file.path("MJM_Bayes.txt")
  write.model(model, filename)

###############################################
# Wrapper function to run jags model in parallel

  coda.samples.wrapper <- function(x)
  {
    model.fit = jags.model(file = "MJM_Bayes.txt",
                           inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = x), 
                           data = Data, n.chains = 1)#, inits = list(b = b.random))
    update(model.fit, n.burn)
    tt = coda.samples(model.fit, parms,  n.iter = n.iter - n.burn, thin = n.thin)
  }
  
  environment(coda.samples.wrapper) <- environment()
  
  ## Number of cores =  number of chains, n.chains
  
  print("Obtaining posterior samples...")
  
  post.samples <- mclapply(1:n.chains, coda.samples.wrapper,  mc.cores = n.chains) 

  
  print("Posterior samples obtained, saving as mcmc list for easy processing")

  for(ii in 1:length(post.samples))
  { post.samples[[ii]] <- post.samples[[ii]][[1]]}
  class(post.samples) <- "mcmc.list"
  
  ##merge results from chains into one data frame.
  
  bss <- do.call(rbind, post.samples)
  n.sims <- nrow(bss)
  all.samples <- vector("list", length(parms))
  names(all.samples) <- parms
  for (p in seq_along(parms))
  {
    ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
    all.samples[[p]] <- bss[, ii]
  }
  
  list(post.samples = all.samples)
}

###Example run:
results = BMTJM_fit(data = sample_data, n.iter = 5000, n.burn = 500, n.thin = 3, n.chains = 3)

save(results, file = "results_BMTJM.rdata") 


# Example code to obtain the mean of the posterior samples
apply(results$post.samples[[1]], 2, mean)
