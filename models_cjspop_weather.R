
# Run cjspop models for maps species

rm(list = ls())

# Set up initial parameters -----------------------------------------------

## Code uses parallel computing. 
## set cores to 1 if not sure about the number of cores 
cores <- 7

## Iterations of MCMC
ni <- 100000

## burn-in
nb <- 50000

## Number of chains 
## Jags.parallel function is used so each chain is run on a seperate core
## total core useage is cores*nc
nc <- 4

## Thinning
nt <- 20


# Load packages -----------------------------------------------------------

if ("pacman" %in% rownames(installed.packages())==F) 
  install.packages("pacman", repos = "http://cran.rstudio.com/")
pacman::p_load(doSNOW, foreach, rjags, R2jags, MCMCvis, magrittr)


# Load data and set up jags functions -------------------------------------

ssdata <- readRDS("data_ss.rds")

# file name argument passed down in pbs 
args <- Sys.getenv("X")
args <- unlist(strsplit(args[1],"[.]"))
n1 <- as.numeric(args[1])
n2 <- as.numeric(args[2])

spcode <- ssdata$spcode[n1:n2]

chdata <- list()
for (i in 1:length(spcode)) {
  filename <- paste(spcode[i], "chdata",  sep ="_") %>%
    paste("rds", sep = ".") %>%
    paste("species", spcode[i], ., sep = "/")
  chdata[[i]] <- readRDS(filename)
}

models  <- c()
nvar <- c()
for (i in 1:length(chdata)) {
  
  nvar[i] <- length(chdata[[i]]$weather_pca)
  models[i] <- paste("cjspop", nvar[i], sep = "") %>%
    paste(".jags", sep = "")
  
}


# ch.init function initilizes values for the latent state z, using the latent_state matrix (see
# data_generation notebook for the creation of this matrix). NA values in the latent state is given
# an initial value of 1 (Years before first or after the last capture of an individual). 
# Observations between the first and last year of captures are given NA as initial values 
# because they are provided as data below. init.func initilizes other model paramters 
# and objects usually byt random draws that correspond to their priors.

ch.init <- function(latent_state, first, last_pop, pop) {
  ch <- ifelse(is.na(latent_state),1,NA)
  for (i in 1:nrow(ch)) {
    ch[i,1:first[i]] <- NA
    if(last_pop[pop[i]]<ncol(ch)) ch[i,(last_pop[pop[i]]+1):ncol(ch)] <- NA
  }
  return(ch)    
}

jags.func <- function(chdata, model, ni, nb, nc, nt, nvar) {
  
  init.func <- function() {
    list(survival_ad = runif(1,0,1),
         beta = rnorm(1,0,1),
         beta_l = rnorm(nvar,0,1),
         gamma = rnorm(2,0,10),
         delta = rnorm(1,0,2.5),
         fecundity = runif(1,1,5),
         zeta = rnorm(1,0,1),
         zeta_l = rnorm(nvar,0,1),
         pi = runif(2,0,1),
         rho = runif(2,0,1),
         sigma_s = runif(1,0,1),
         sigma_f = runif(1,0,1),
         psi = runif(1,0,1),
         R = rep(1, nrow(chdata$stage)),
         z = ch.init(chdata$latent_state, chdata$first, 
                     chdata$last_pop, chdata$pop))
  }
  
  params <- c("survival_ad", "survival_juv", "beta_l", "beta", 
              "zeta_l", "zeta", "fecundity", "pi", "rho", "sigma_f", "psi",
              "sigma_s", "gamma", "delta")
  
  data <- list(y = chdata$ch_robust,
               z = chdata$latent_state,
               r = chdata$residents,
               effort = chdata$effort_cent,
               effort_year = chdata$effort_year,
               Nobs_ad = chdata$Nobs_ad,
               Nobs_juv1 = chdata$Nobs_juv,
               Nobs_juv2 = chdata$fec_data,
               w1 = chdata$weather_pca[[1]],
               w2 = chdata$weather_pca[[2]],
               w3 = chdata$weather_pca[[3]],
               w4 = chdata$weather_pca[[4]],
               w1_fec = chdata$weather_fec_pca[[1]],
               w2_fec = chdata$weather_fec_pca[[2]],
               w3_fec = chdata$weather_fec_pca[[3]],
               w4_fec = chdata$weather_fec_pca[[4]],
               pop = chdata$pop,
               stage = chdata$stage,
               first = chdata$first,
               first_pop = chdata$first_pop,
               last_pop = chdata$last_pop,
               first_sub = chdata$first_sub,
               pop_index = chdata$pop_index,
               year_index = chdata$year_index,
               nind = chdata$nind,
               nyear = chdata$nyear,
               nsub = 4,
               npop = chdata$npop,
               nfec = chdata$nfec)
  
  # Add additional variables if more than 4
  if (nvar > 4) {
    data$w5 <- chdata$weather_pca[[5]]
    data$w5_fec <- chdata$weather_fec_pca[[5]]
  }
  
  if (nvar > 5) {
    data$w6 <- chdata$weather_pca[[6]]
    data$w6_fec <- chdata$weather_fec_pca[[6]]
  }
  
  if (nvar > 6) {
    data$w7 <- chdata$weather_pca[[7]]
    data$w7_fec <- chdata$weather_fec_pca[[7]]
  }
  
  results <- jags.parallel(model.file = model,
                           working.directory = getwd(),
                           data = data,
                           inits = init.func,
                           parameters.to.save = params,
                           n.iter = ni,
                           n.chains = nc,
                           n.burnin = nb,
                           n.thin = nt,
                           DIC = F)
  
  return(results)
}


# Run the jags model ------------------------------------------------------

cl <- makeCluster(cores, types = "SOCK")
registerDoSNOW(cl)

system.time({
  results <- foreach(i=1:length(chdata), .packages = c("rjags", "R2jags")) %dopar% {
    jags.func(chdata[[i]], models[i], nvar = nvar[i], ni = ni, nb = nb, nc = nc, nt = nt)
  }
})

stopCluster(cl)

filename <- paste("results_cjspop_weather", n1, n2, sep = "_") %>%
  paste("rds", sep = ".")
saveRDS(results, filename)
