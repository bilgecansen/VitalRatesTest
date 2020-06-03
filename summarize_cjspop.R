
rm(list = ls())

# Packages and data -------------------------------------------------------

library(tidyverse)
library(foreach)
library(MCMCvis)
library(boot)
library(abind)

data_ss <- readRDS("data_ss.rds")[c(1:42),]
spcode <- data_ss$spcode

chdata <- list()
for (i in 1:length(spcode)) {
  filename <- paste(spcode[i], "chdata",  sep ="_") %>%
    paste("rds", sep = ".") %>%
    paste("species", spcode[i], ., sep = "/")
  chdata[[i]] <- readRDS(filename)
}

res1 <- readRDS("results_cjspop_weather_1_14.rds")
res2 <- readRDS("results_cjspop_weather_15_21.rds")
res3 <- readRDS("results_cjspop_weather_22_28.rds")
res4 <- readRDS("results_cjspop_weather_29_35.rds")
res5 <- readRDS("results_cjspop_weather_36_42.rds")
res <- c(res1, res2, res3, res4, res5)

mcmc_sum <- map(res, MCMCsummary)
excl <- c("pi", "rho", "psi")
mcmc_chains <- map(res, function(x) MCMCchains(x, excl = excl))


# Species selection -------------------------------------------------------

# Remove species with unconverged paramaters
index <- map_lgl(mcmc_sum, function(x) any(x[,"Rhat"]>=1.1)) %>%
  which()
mcmc_sum <- mcmc_sum[-index]
mcmc_chains <- mcmc_chains[-index]
chdata <- chdata[-index]
spcode <- spcode[-index]

# Remove species with positive density dependence on R
sad <- map_dbl(mcmc_sum, function(x) x["survival_ad", "mean"])
sjuv <- map_dbl(mcmc_sum, function(x) x["survival_juv", "mean"])
fec <- map_dbl(mcmc_sum, function(x) x["fecundity", "mean"])
beta <- map_dbl(mcmc_sum, function(x) x["beta", "mean"])
zeta <- map_dbl(mcmc_sum, function(x) x["zeta", "mean"])
density <- 0:2

sad_dens <- list()
for (i in 1:length(sad)) {
  sad_dens[[i]] <- inv.logit(logit(sad[i]) + beta[i]*density)
}

sjuv_dens <- list()
for (i in 1:length(sjuv)) {
  sjuv_dens[[i]] <- inv.logit(logit(sjuv[i]) + beta[i]*density)
}

fec_dens <- list()
for (i in 1:length(fec)) {
  fec_dens[[i]] <- exp(log(fec[i]) + zeta[i]*density)
}

R_dens <- list()
for (i in 1:length(sad)) {
  R_dens[[i]] <- sad_dens[[i]] + sjuv_dens[[i]]*fec_dens[[i]]
}

R_trend <- map_dbl(R_dens, function(x) x[3]/x[1])
index2 <- which(R_trend>1)

mcmc_sum <- mcmc_sum[-index2]
mcmc_chains <- mcmc_chains[-index2]
chdata <- chdata[-index2]
spcode <- spcode[-index2]

results_cjspop <- list(mcmc_sum = mcmc_sum,
                       mcmc_chains = mcmc_chains,
                       chdata = chdata,
                       spcode = spcode)

saveRDS(results_cjspop, file = "results_cjspop.rds")
