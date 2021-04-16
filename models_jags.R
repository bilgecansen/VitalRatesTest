
# Jags models for summarizing CJS-pop and HMSC results

# Functions ---------------------------------------------------------------

write.jags <- function(name, script) {
  
  sink(name)
  cat(script, fill = TRUE)
  sink()
  
}

lm.jags <- function(data, inits, n.chains, n.adapt, 
                    n.update, n.iter, n.thin, seed_no) {
  
  script <- "
  
    model {
      
      # Priors
      
      alpha ~ dnorm(0, pow(100, -2))
      beta ~ dnorm(0, pow(100, -2))
      sigma ~ dunif(0,1000)
    
      # Likelihood
      
      for(i in 1:n) {

        mu[i] <- alpha + beta*x[i]
        y[i] ~ dnorm(mu[i], pow(sigma, -2))
        
        ## calculate residuals
        res[i] <- (mu[i]- mean(y))^2
        error[i] <- (y[i] - mu[i])^2
        
      }
      
      R.sq <- (sum(res)/(n-1))/((sum(res)/(n-1)) + (sum(error)/(n-1)))
      

    }
  "
  
  script_name <- "lm.jags"
  write.jags(name = script_name, script = script)
  
  set.seed(seed_no)
  
  jm <- rjags::jags.model(script_name, 
                          data = data,
                          inits = inits,
                          n.chains = n.chains, 
                          n.adapt = n.adapt,
                          quiet = T)
  
  update(jm, n.iter = n.update)
  
  mcmc_samples <- coda.samples(jm, 
                               variable.names = c("alpha", "beta","sigma", "R.sq"),
                               n.iter = n.iter, 
                               n.thin = n.thin)
  
  mcmc_sum <- MCMCsummary(mcmc_samples, digits = 3)
  
  results <- list(mcmc_samples = mcmc_samples,
                  mcmc_sum = mcmc_sum)
  
  return(results)
  
}

cor.jags <- function(data, inits, n.chains, n.adapt, 
                     n.update, n.iter, n.thin, seed_no) {
  
  script <- "
  
    model {
      
      # Priors
      
      mu[1] ~ dnorm(0, pow(10, -2))
      mu[2] ~ dnorm(0, pow(10, -2))
      sigma[1] ~ dunif(0, 100)
      sigma[2] ~ dunif(0, 100)
      rho ~ dunif(-1, 1)
      
      Sigma[1,1] <- sigma[1]*sigma[1]
      Sigma[2,2] <- sigma[2]*sigma[2]
      Sigma[1,2] <- sigma[1]*sigma[2]*rho
      Sigma[2,1] <- sigma[1]*sigma[2]*rho
      Omega <- inverse(Sigma[,])
      
      # Likelihood
      
      for (i in 1:n) {
      
        y[i, 1:2] ~ dmnorm(mu[], Omega[,])
      
      }

    }
  "
  
  script_name <- "cor.jags"
  write.jags(name = script_name, script = script)
  
  set.seed(seed_no)
  
  jm <- rjags::jags.model(script_name, 
                          data = data,
                          inits = inits,
                          n.chains = n.chains, 
                          n.adapt = n.adapt)
  
  update(jm, n.iter = n.update)
  
  mcmc_samples <- coda.samples(jm, 
                               variable.names = c("mu", "sigma", "rho"),
                               n.iter = n.iter, 
                               n.thin = n.thin)
  
  mcmc_sum <- MCMCsummary(mcmc_samples, digits = 3)
  
  results <- list(mcmc_samples = mcmc_samples,
                  mcmc_sum = mcmc_sum)
  
  return(results)
  
}

beta.lm.jags <- function(data, inits, n.chains, n.adapt, 
                         n.update, n.iter, n.thin, seed_no) {
  
  script <- "
  
    model {
      
      # Priors
      
      phi ~ dunif(0,100)
      alpha ~ dnorm(0, pow(10, -2))
      beta ~ dnorm(0, pow(10, -2))
    
      # Likelihood
      
      for(i in 1:n) {

        logit(mu[i]) <- alpha + beta*x[i]
        
        a[i] <- mu[i]*phi
        b[i]  <- (1-mu[i])*phi
      
        y[i] ~ dbeta(a[i], b[i])
        
        ## calculate residuals
        res[i] <- (mu[i]- mean(y))^2
        error[i] <- (y[i] - mu[i])^2
        
      }
      
      R.sq <- (sum(res)/(n-1))/((sum(res)/(n-1)) + (sum(error)/(n-1)))
      

    }
  "
  
  script_name <- "beta_lm.jags"
  write.jags(name = script_name, script = script)
  
  set.seed(seed_no)
  
  jm <- rjags::jags.model(script_name, 
                          data = data,
                          inits = inits,
                          n.chains = n.chains, 
                          n.adapt = n.adapt,
                          quiet = T)
  
  update(jm, n.iter = n.update)
  
  mcmc_samples <- coda.samples(jm, 
                               variable.names = c("alpha", "beta", "phi", "R.sq"),
                               n.iter = n.iter, 
                               n.thin = n.thin)
  
  mcmc_sum <- MCMCsummary(mcmc_samples, digits = 3)
  
  results <- list(mcmc_samples = mcmc_samples,
                  mcmc_sum = mcmc_sum)
  
  return(results)
  
}

cor.alt.jags <- function(data, inits, n.chains, n.adapt, 
                         n.update, n.iter, n.thin) {
  
  script <- "
  
    model {
      
      # Priors
      
      mu[1] ~ dnorm(0, pow(10, -2))T(0,)
      mu[2] ~ dnorm(0, pow(10, -2))T(0,)
      sigma[1] ~ dunif(0, 100)
      sigma[2] ~ dunif(0, 100)
      rho ~ dunif(-1, 1)
      
      Sigma[1,1] <- sigma[1]*sigma[1]
      Sigma[2,2] <- sigma[2]*sigma[2]
      Sigma[1,2] <- sigma[1]*sigma[2]*rho
      Sigma[2,1] <- sigma[1]*sigma[2]*rho
      Omega <- inverse(Sigma[,])
      
      # Likelihood
      
      for (i in 1:n) {
      
        y[i, 1:2] ~ dmnorm(mu[], Omega[,])
      
      }

    }
  "
  
  script_name <- "cor_alt.jags"
  write.jags(name = script_name, script = script)
  
  
  jm <- rjags::jags.model(script_name, 
                          data = data,
                          inits = inits,
                          n.chains = n.chains, 
                          n.adapt = n.adapt)
  
  update(jm, n.iter = n.update)
  
  mcmc_samples <- coda.samples(jm, 
                               variable.names = c("mu", "sigma", "rho"),
                               n.iter = n.iter, 
                               n.thin = n.thin)
  
  mcmc_sum <- MCMCsummary(mcmc_samples, digits = 3)
  
  results <- list(mcmc_samples = mcmc_samples,
                  mcmc_sum = mcmc_sum)
  
  return(results)
  
}


