
####Install packages

library(rjags)
library(R2WinBUGS)
library(JMbayes2)
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
tic()
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
  names(sims.list) <- parms
  for (p in seq_along(parms))
  {
    ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
    all.samples[[p]] <- bss[, ii]
  }
  
  list(post.samples = all.samples)
}

###Example run:
results = BMTJM_fit(data = sample_data, n.iter = 100, n.burn = 5, n.thin = 2, n.chains = 3)
toc()
save(results, file = "results_BMTJM.rdata") 


# Example code to obtain the mean of the posterior samples
apply(results$post.samples[[1]], 2, mean)


Bs.betas0_est <- apply(sims.list$Bs.betaX0, 2, mean)
Bs.betas1_est <- apply(sims.list$Bs.betaX1, 2, mean)  
Bs.betas2_est <- apply(sims.list$Bs.betaX2, 2, mean)  

Bs.psi1_est <- apply(sims.list$Bs.betaZ1, 2, mean)
Bs.psi2_est <- apply(sims.list$Bs.betaZ2, 2, mean)

Bs.alpha_est <- apply(sims.list$Bs.alpha, 2, mean)
Bs.gamma_est <- apply(sims.list$Bs.h0, 2, mean)
Bs.zeta1_est <- apply(sims.list$Bs.gammaX1, 2, mean)
Bs.eta1_est <- apply(sims.list$Bs.gammaX2, 2, mean)  

Bs.zeta2_est <- apply(sims.list$Bs.gammaX2, 2, mean)
Bs.eta2_est <- apply(sims.list$Bs.gammaZ2, 2, mean)  

beta0.mat = bases %*% (Bs.betas0_est)
beta1.mat = bases %*% (Bs.betas1_est)
beta2.mat = bases %*% (Bs.betas2_est)

psi1.mat = bases %*% (Bs.psi1_est)  
psi2.mat = bases %*% (Bs.psi2_est)  

alphas.mat <- bases %*% (Bs.alpha_est)
gammas.mat <- bases %*% (Bs.gamma_est)
zeta1.mat <- bases %*% (Bs.zeta1_est)
eta1.mat <- bases %*% (Bs.eta1_est)

zeta2.mat <- bases %*% (Bs.zeta2_est)
eta2.mat <- bases %*% (Bs.eta2_est)

sigmasqb_est <- mean(sims.list$sigmasq_S)
sigmasqxi_est <- mean(sims.list$sigmasq_F)

beta0_true <- cos(times*pi*3/2) - 0.5
beta1_true <- sin(times*2*pi-1/8)
beta2_true <- -sin(times*2*pi-1/8)

psi1_true <- cos(times*pi-0.5)
psi2_true <- -cos(times*pi-0.5)


zeta1_true <- cos(times*2*pi)
zeta2_true <- -cos(times*2*pi)

eta1_true <- sin(times*0.75*pi)
eta2_true <- -sin(times*0.75*pi)

alpha_true <- sin(times*2*pi)
sigmasqb_true <- 1.43
sigmasqxi_true <- 0.30


est = gammas.mat
medianv = rowMedians(est)
####LOW PRIORITY: IF USING RASE AND PLOTTING MEDIAN IS NOT THE ACCEPTED ROUTE THAN WE'LL THINK ABOUT THIS.
##sdv = rowVars(est) ##I need to think about what would be the estimated SD if I were to plot CI for all results -- get SD from each sim and average?
###I think we should follow something similar to the data first find beta zero hat using all post samples, then
## get the variance amond those samples and then take the mean --all.beta0.mat <- bases %*% t(sims.list$Bs.betas0)
##sd.beta0.matmed <- sqrt(rowVars(all.beta0.mat))
quants = rowQuantiles(est, probs = c(0.025, 0.975))  ymin = c(min(quants[, 1])) ymax = c(max(quants[, 2]))

phi <- 1.5

pdf("baseline_sim_est.pdf")
par(cex.axis = 2.0)
par(cex.lab = 2.0)
par(mai = c(1.2, 1.2, 0.5, 0.5))
plot(times, exp(est) ,"l", lwd = 3, ylab = expression(hat(h)[0](t)), xlab = "t", xaxs = "i", yaxs = "i", ylim = c(0, 2.5),  lty = 2)
lines(times, exp(quants[,1]), lwd = 3, lty = 3)
lines(times, exp(quants[,2]), lwd = 3, lty = 3)
lines(times, exp(log(phi) + (phi - 1) * log(times)), lwd = 3, lty = 1)
dev.off()


est = beta0.mat
medianv = est
####LOW PRIORITY: IF USING RASE AND PLOTTING MEDIAN IS NOT THE ACCEPTED ROUTE THAN WE'LL THINK ABOUT THIS.
##sdv = rowVars(est) ##I need to think about what would be the estimated SD if I were to plot CI for all results -- get SD from each sim and average?
###I think we should follow something similar to the data first find beta zero hat using all post samples, then
## get the variance amond those samples and then take the mean --all.beta0.mat <- bases %*% t(sims.list$Bs.betas0)
##sd.beta0.matmed <- sqrt(rowVars(all.beta0.mat))

quants = rowQuantiles(est, probs = c(0.025, 0.975))
ymin = c(min(quants[, 1]))
ymax = c(max(quants[, 2]))


pdf("beta0_sim_est.pdf")
par(cex.axis = 2.0)
par(cex.lab = 2.0)
par(mai = c(1.2, 1.2, 0.5, 0.5))
plot(times, medianv ,"l", lwd = 3, ylab = expression(hat(beta)[X0](t)), xlab = "t", xaxs = "i", yaxs = "i", ylim = c(-2, 2), lty = 2)
lines(times, quants[,1], lwd = 3, lty = 3)
lines(times, quants[,2], lwd = 3, lty = 3)
#lines(times, meanv - 2*sdv, lwd = 3, lty = 3)
#lines(times, meanv + 2*sdv, lwd = 3, lty = 3)
lines(times, beta0_true, lty = 1, lwd = 3)
dev.off()



est = beta1.mat[,]
medianv = rowMedians(est)
##sdv = rowVars(est) ##I need to think about what would be the estimated SD if I were to plot CI for all results -- get SD from each sim and average?

quants = rowQuantiles(est, probs = c(0.025, 0.975))
ymin = c(min(quants[, 1]))
ymax = c(max(quants[, 2]))

pdf("beta1_sim_est.pdf")
par(cex.axis = 2.0)
par(cex.lab = 2.0)
par(mai = c(1.2, 1.2, 0.5, 0.5))
plot(times, medianv ,"l", lwd = 3, ylab = expression(hat(beta)[X1](t)), xlab = "t", xaxs = "i", yaxs = "i", ylim = c(ymin, ymax), lty = 2)
lines(times, quants[,1], lwd = 3, lty = 3)
lines(times, quants[,2], lwd = 3, lty = 3)
#lines(times, meanv - 2*sdv, lwd = 3, lty = 3)
#lines(times, meanv + 2*sdv, lwd = 3, lty = 3)
lines(times, beta1_true, lty = 1, lwd = 3)
dev.off()



