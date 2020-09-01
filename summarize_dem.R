
# Load packages and data --------------------------------------------------

library(tidyverse)
library(MCMCvis)
library(boot)
library(doSNOW)
library(foreach)
library(abind)

results_cjspop <- readRDS("results_cjspop.rds")
pop_station <- readRDS("pop_station.rds")

mcmc_chains <- results_cjspop$mcmc_chains
mcmc_sum <- results_cjspop$mcmc_sum
chdata <- results_cjspop$chdata
spcode <- results_cjspop$spcode


# Functions ---------------------------------------------------------------

prep.wdata <- function(chdata) {
  
  wdata <- map(chdata$weather_pca, function(x) array(x, dim = c(chdata$npop, chdata$nyear, 1)))
  wdata <- do.call(abind, list(wdata, along = 3))
  
  return(wdata)
  
}

estimate.avgR <- function(chains, wdata, time, pop) {
  
  alpha1 <- logit(chains[,"survival_juv"])
  alpha2 <- logit(chains[,"survival_ad"])
  beta <- chains[,"beta"] 
  theta <- log(chains[,"fecundity"]) 
  zeta <- chains[,"zeta"]
  
  index1 <- which(grepl("beta_l", colnames(chains)))
  index2 <- which(grepl("zeta_l", colnames(chains)))
  swchains <- chains[,index1]
  fwchains <- chains[,index2]
  
  R <- array(dim = c(pop, time))
  s_ad <- array(dim = c(pop, time))
  s_juv <- array(dim = c(pop, time))
  fecundity <- array(dim = c(pop, time))
  pb <- txtProgressBar(min = 0, max = pop*time, style = 3)
  
  a <- 0
  for (k in 1:pop) {
    
    for (t in 1:time) {
      
      ## %*% is matrix multiplication
      ## swchains is matrix of betas. Each column is a beta for a different variable.
      ## Each row is a different of the mcmc chain (posterior distribution).
      ## wdata[k,t,] is a vector where each element is a varible for population k and time t.
      s1 <- inv.logit(alpha1 + beta*(-1) + swchains %*% wdata[k,t,])
      s2 <- inv.logit(alpha2 + beta*(-1) + swchains %*% wdata[k,t,])
      fec <-  exp(theta + zeta*(-1) + fwchains %*% wdata[k,t,])
      
      R[k,t] <- mean(s2 + s1*fec) 
      s_ad[k,t] <- mean(s2)
      s_juv[k,t] <- mean(s1)
      fecundity[k,t] <- mean(fec)
      
      
      setTxtProgressBar(pb, a)
      a <- a + 1
      
    }#t
    
  }#k
  
  res <- list(R = R,
              s_ad = s_ad,
              s_juv = s_juv,
              fecundity = fecundity)
  
  return(res)
  
}

estimate.p <- function(delta, gamma1, gamma2, chdata) {
  
  pad_month <- inv.logit(gamma2 + delta*chdata$effort_cent)
  pjuv_month <- inv.logit(gamma1 + delta*chdata$effort_cent)
  
  
  npop <- length(chdata$first_pop)
  pad_year <- matrix(nrow = npop, ncol = chdata$nyear)
  pjuv_year <- matrix(nrow = npop, ncol = chdata$nyear)
  
  for (k in 1:npop) {
    for (t in 1:chdata$nyear) {
      pad_year[k,t] <- (1-(1-pad_month[k,t,1])*
                          (1-pad_month[k,t,2])*
                          (1-pad_month[k,t,3])*
                          (1-pad_month[k,t,4]))
      
      pjuv_year[k,t] <- (1-(1-pjuv_month[k,t,1])*
                           (1-pjuv_month[k,t,2])*
                           (1-pjuv_month[k,t,3])*
                           (1-pjuv_month[k,t,4]))
    }
  }
  
  results <- list()
  results$pad_year <- pad_year
  results$pjuv_year <- pjuv_year
  results$ratio <- inv.logit(gamma2)/inv.logit(gamma1)
  results$gamma1 <- inv.logit(gamma1)
  results$gamma2 <- inv.logit(gamma2)
  
  return(results)
}

estimate.N <- function(p, chdata) {
  
  pad_year <- p$pad_year
  pjuv_year <- p$pjuv_year
  Nobs_ad <- chdata$Nobs_ad
  Nobs_juv <- chdata$Nobs_juv
  
  # Population size estimated by heuristic estimator
  results <- list()
  corr_ad <- (1-pad_year)/pad_year
  corr_juv <- (1-pjuv_year)/pjuv_year
  results$corr_Nad <- Nobs_ad/pad_year + corr_ad
  results$corr_Njuv <- Nobs_juv/pjuv_year + corr_juv
  results$corr_N <- results$corr_Nad + results$corr_Njuv
  results$corr_D <- t(apply(results$corr_N, 1, function(x) x/mean(x, na.rm = T)))
  results$avgN <- apply(results$corr_N, 1, mean, na.rm = T)
  results$highN <- apply(results$corr_N, 1, function(x) quantile(x, 0.95, na.rm = T))
  results$sdN <- apply(results$corr_N, 1, sd, na.rm = T)
  results$cvN <- results$sdN/results$avgN
  results$Nobs_juv <- Nobs_juv
  results$Nobs_ad <- Nobs_ad
  
  return(results)
}

# Demography estimation ---------------------------------------------------

wdata <- map(chdata, prep.wdata)

R <- foreach (i=1:length(mcmc_chains), .packages = c("abind", "boot")) %do% {
                
  estimate.avgR(chains = mcmc_chains[[i]],
                wdata = wdata[[i]],
                time = chdata[[i]]$nyear,
                pop = chdata[[i]]$npop)
   
}

names(R) <- spcode

# Calucate number of stations for each population
# This is used for standardizing population size to per station
pops <- map(results_cjspop$chdata, function(x) as.numeric(rownames(x$Nobs_ad)))
unique_pop <- unique(pop_station$pop)
unique_pop <- unique_pop[order(unique_pop)]

nstat <- c()
for (i in 1:length(unique_pop)) { 
  z <- filter(pop_station, pop==unique_pop[i])
  nstat[i] <- length(unique(z$station))
}

nstat_sp <- map(pops, function(x) nstat[x])
names(nstat_sp) <- spcode

# Posterior distribution of average local population size per species
cl <- makeCluster(4, type = "SOCK")
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(spcode), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

N <- foreach(i=1:length(spcode), .options.snow = opts, .packages = "boot") %dopar% {
  
  chains <- results_cjspop$mcmc_chains[[i]]
  npop <- results_cjspop$chdata[[i]]$npop
  
  results <- list()
  spN <- c()
  highN <- matrix(nrow = nrow(chains), ncol = npop)
  matN <- matrix(nrow = nrow(chains), ncol = npop)
  
  for (h in 1:nrow(chains)) {
    
    p <- estimate.p(delta = chains[h,"delta"],
                    gamma1 = chains[h,"gamma[1]"],
                    gamma2 = chains[h,"gamma[2]"],
                    chdata[[i]])
    tempN <- estimate.N(p, chdata[[i]])
    
    matN[h,] <- tempN$avgN/nstat_sp[[i]]
    spN[h] <- median(tempN$avgN/nstat_sp[[i]])
    highN[h,] <- tempN$highN/nstat_sp[[i]]

  }
  
  # Average abundance of each population
  results$popN <- apply(matN, 2, mean)
  
  # 95% quantile of each population
  results$highN <- apply(highN, 2, mean)
  
  # Posterior of species level median population abundance
  results$spN <- spN

  return(results)
  
}
stopCluster(cl)
names(N) <- spcode
  
results_dem <- list()
results_dem$R <- R
results_dem$N <- N
results_dem$nstat_sp <- nstat_sp

saveRDS(results_dem, "results_dem.rds")

