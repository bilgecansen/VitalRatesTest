

# CJS-pop with 4 variables ------------------------------------------------

sink("cjspop4.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:4) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:4) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] + 
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +  
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] +  
                      zeta_l[2]*w2_fec[l] + 
                      zeta_l[3]*w3_fec[l] +  
                      zeta_l[4]*w4_fec[l] +  
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()


# CJS-pop with 5 variables ------------------------------------------------

sink("cjspop5.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:5) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:5) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +  
                             beta_l[5]*w5[k,t] + 
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] + 
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +  
                             beta_l[5]*w5[k,t] + 
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] +  
                      zeta_l[2]*w2_fec[l] + 
                      zeta_l[3]*w3_fec[l] +  
                      zeta_l[4]*w4_fec[l] +  
                      zeta_l[5]*w5_fec[l] +  
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()


# CJS-pop with 6 variables ------------------------------------------------

sink("cjspop6.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:6) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:6) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +  
                             beta_l[5]*w5[k,t] +
                             beta_l[6]*w6[k,t] +
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] + 
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +  
                             beta_l[5]*w5[k,t] +
                             beta_l[6]*w6[k,t] +
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] +  
                      zeta_l[2]*w2_fec[l] + 
                      zeta_l[3]*w3_fec[l] +  
                      zeta_l[4]*w4_fec[l] +  
                      zeta_l[5]*w5_fec[l] +
                      zeta_l[6]*w6_fec[l] +
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()


# CJS-pop with 7 variables ------------------------------------------------

sink("cjspop7.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:7) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:7) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +  
                             beta_l[5]*w5[k,t] +
                             beta_l[6]*w6[k,t] +
                             beta_l[7]*w7[k,t] +
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] + 
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] +  
                             beta_l[2]*w2[k,t] +  
                             beta_l[3]*w3[k,t] +  
                             beta_l[4]*w4[k,t] +  
                             beta_l[5]*w5[k,t] +
                             beta_l[7]*w7[k,t] +
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] +  
                      zeta_l[2]*w2_fec[l] + 
                      zeta_l[3]*w3_fec[l] +  
                      zeta_l[4]*w4_fec[l] +  
                      zeta_l[5]*w5_fec[l] +
                      zeta_l[6]*w6_fec[l] +
                      zeta_l[7]*w7_fec[l] +
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()


# Quadratic cjs-pop with 4 variables --------------------------------------

sink("cjspop4q.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:4) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
      beta_q[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:4) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
      zeta_q[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] + zeta_q[1]*(w1_fec[l]^2) +  
                      zeta_l[2]*w2_fec[l] + zeta_q[2]*(w2_fec[l]^2) +
                      zeta_l[3]*w3_fec[l] + zeta_q[3]*(w3_fec[l]^2) + 
                      zeta_l[4]*w4_fec[l] + zeta_q[4]*(w4_fec[l]^2) + 
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()



# Quadratic cjs-pop with 5 variables --------------------------------------

sink("cjspop5q.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:5) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
      beta_q[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:5) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
      zeta_q[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +  
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +  
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +  
                             beta_l[5]*w5[k,t] + beta_q[5]*(w5[k,t]^2) +
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] + 
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +  
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +  
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +  
                             beta_l[5]*w5[k,t] + beta_q[5]*(w5[k,t]^2) +
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] + zeta_q[1]*(w1_fec[l]^2) +  
                      zeta_l[2]*w2_fec[l] + zeta_q[2]*(w2_fec[l]^2) +
                      zeta_l[3]*w3_fec[l] + zeta_q[3]*(w3_fec[l]^2) + 
                      zeta_l[4]*w4_fec[l] + zeta_q[4]*(w4_fec[l]^2) + 
                      zeta_l[5]*w5_fec[l] + zeta_q[5]*(w5_fec[l]^2) + 
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()



# Quadratic cjs-pop with 6 variables  -------------------------------------

sink("cjspop6q.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:6) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
      beta_q[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:6) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
      zeta_q[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +  
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +  
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +  
                             beta_l[5]*w5[k,t] + beta_q[5]*(w5[k,t]^2) +
                             beta_l[6]*w6[k,t] + beta_q[6]*(w6[k,t]^2) +
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] + 
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +  
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +  
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +  
                             beta_l[5]*w5[k,t] + beta_q[5]*(w5[k,t]^2) +
                             beta_l[6]*w6[k,t] + beta_q[6]*(w6[k,t]^2) +
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] + zeta_q[1]*(w1_fec[l]^2) + 
                      zeta_l[2]*w2_fec[l] + zeta_q[2]*(w2_fec[l]^2) +
                      zeta_l[3]*w3_fec[l] + zeta_q[3]*(w3_fec[l]^2) +  
                      zeta_l[4]*w4_fec[l] + zeta_q[4]*(w4_fec[l]^2) + 
                      zeta_l[5]*w5_fec[l] + zeta_q[5]*(w5_fec[l]^2) +
                      zeta_l[6]*w6_fec[l] + zeta_q[6]*(w6_fec[l]^2) +
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()



# Quadratic cjs-pop with 7 variables --------------------------------------

sink("cjspop7q.jags")
cat("
    model {
    
    # Part 1: Priors
    
    # Part 1.1: Priors for survival model
    
    survival_ad ~ dunif(0,1)
    survival_juv <- (1-survival_ad)/fecundity
    alpha[1] <- logit(survival_juv)
    alpha[2] <- logit(survival_ad)
    beta ~ dt(0,pow(10,-2),1)
    sigma_s ~ dt(0,pow(10,-2),1)T(0,)
    tau_s <- 1/(sigma_s^2)
    
    for (i in 1:7) {
    
      beta_l[i] ~ dt(0,pow(10,-2),1)
      beta_q[i] ~ dt(0,pow(10,-2),1)
    
    }
    
    for (t in 1:(nyear-1)) {
      
      epsilon[t] ~ dnorm(0, tau_s)
    
    }
    
    # Part 1.2: Priors for capture probability model
    
    gamma[1] ~ dt(0,pow(10,-2),1)
    gamma[2] ~ dt(0,pow(10,-2),1)
    delta ~ dt(0,pow(2.5,-2),1)
    
    # Part 1.3: Priors for fecundity model
    
    fecundity ~ dnorm(0,pow(10,-2))T(0,)
    theta <- log(fecundity)
    zeta ~ dt(0,pow(5,-2),1)
    sigma_f ~ dt(0,pow(5,-2),1)T(0,)
    tau_f <- 1/(sigma_f^2)
    psi ~ dunif(0,1)
    
    for (i in 1:7) {
    
      zeta_l[i] ~ dt(0,pow(5,-2),1)
      zeta_q[i] ~ dt(0,pow(5,-2),1)
    
    }
    
    for (t in 2:nyear) {
      
      omega[t] ~ dnorm(0, tau_f)
    
    }
    
    # Part 1.4: Priors for residency model
    
    pi[1] ~ dunif(0,1)
    pi[2] ~ dunif(0,1)
    rho[1] ~ dunif(0,1)
    rho[2] ~ dunif(0,1)
    
    # Part 2: Robust design mark-recapture model
    
    # Part 2.1: Survival constraints
    
    for (k in 1:npop) {
      
      for (t in first_pop[k]:(last_pop[k]-1)) {
    
        phi[1,k,t] <- ilogit(alpha[1] +
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +  
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +  
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +  
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +  
                             beta_l[5]*w5[k,t] + beta_q[5]*(w5[k,t]^2) +
                             beta_l[6]*w6[k,t] + beta_q[6]*(w6[k,t]^2) +
                             beta_l[7]*w7[k,t] + beta_q[7]*(w7[k,t]^2) +
                             epsilon[t])
    
    
        phi[2,k,t] <- ilogit(alpha[2] + 
                             beta*(density[k,t]-1) +
                             beta_l[1]*w1[k,t] + beta_q[1]*(w1[k,t]^2) +  
                             beta_l[2]*w2[k,t] + beta_q[2]*(w2[k,t]^2) +  
                             beta_l[3]*w3[k,t] + beta_q[3]*(w3[k,t]^2) +  
                             beta_l[4]*w4[k,t] + beta_q[4]*(w4[k,t]^2) +  
                             beta_l[5]*w5[k,t] + beta_q[5]*(w5[k,t]^2) +
                             beta_l[6]*w6[k,t] + beta_q[6]*(w6[k,t]^2) +
                             beta_l[7]*w7[k,t] + beta_q[7]*(w7[k,t]^2) +
                             epsilon[t])
    
      }#k
    }#t
    
    for (i in 1:nind) {
    
      # Part 2.2: Residency Process
      
      R[i] ~ dbern(pi[stage[i,first[i]]])
      r[i] ~ dbern(R[i]*rho[stage[i,first[i]]])
      
      # Part 2.3: Survival Process
      
      z[i,first[i]] ~ dbern(1)
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        z[i,t+1] ~ dbern(phi[stage[i,t],pop[i],t]*z[i,t]*R[i])
      
      }#t
      
      # Part 2.4: Observation Process
      
      y[i,first[i],first_sub[i]] ~ dbern(1)
      y_pred[i,first[i],first_sub[i]] ~ dbern(1)
      
      for (h in (first_sub[i]+1):nsub) {
      
        p[i,first[i],h] <- ilogit(gamma[stage[i,first[i]]] + delta*effort[pop[i],first[i],h])
        y[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
        y_pred[i,first[i],h] ~ dbern(z[i,first[i]]*p[i,first[i],h])
      
      }#h
      
      for (t in first[i]:(last_pop[pop[i]]-1)) {
      
        for(h in 1:nsub) {
          p[i,t+1,h] <- ilogit(gamma[stage[i,t+1]] + delta*effort[pop[i],t+1,h])
          y[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
          y_pred[i,t+1,h] ~ dbern(z[i,t+1]*p[i,t+1,h])
        
        }#h
      }#t
    }#i
    
    # Part 3: Population size and Density Estimation
    
    # Part 3.1: Probability of capture per population and year
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
      
        py_j[k,t] <- 1 - ((1-ilogit(gamma[1] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[1] + delta*effort[k,t,4])))
    
        py_a[k,t] <- 1 - ((1-ilogit(gamma[2] + delta*effort[k,t,1]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,2]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,3]))*
                          (1-ilogit(gamma[2] + delta*effort[k,t,4])))
                          
      }#t
    }#k
    
    # Part 3.2: Density estimation
    
    for (k in 1:npop) {
    
      for (t in first_pop[k]:last_pop[k]) {
    
        Nad[k,t] <- Nobs_ad[k,t]/py_a[k,t] + (1-py_a[k,t])/py_a[k,t]
        Njuv1[k,t] <- Nobs_juv1[k,t]/py_j[k,t] + (1-py_j[k,t])/py_j[k,t]
        
        ## Total population size      
        N_temp[k,t] <- Nad[k,t] + Njuv1[k,t]
        N[k,t] <- ifelse(effort_year[k,t]==0, 0, N_temp[k,t])
        
        ## Count the number of years without effort in a population. 
        ## This is used for calculating mean pop size.
        n1[k,t] <- equals(effort_year[k,t], 0)
        
        ## Population Density (Assign 1 if there is no effort)
        density[k,t] <- ifelse(effort_year[k,t]==0, 1, N[k,t]/mean.N[k])
      }#t
    
      ## Counting the number of populations with effort across study period
      n2[k] <- (last_pop[k] - first_pop[k] + 1) - sum(n1[k,first_pop[k]:last_pop[k]])
    
      ## Mean pop size of effort years
      mean.N[k] <- sum(N[k,first_pop[k]:last_pop[k]])/n2[k]
    
    }#k
    
    # Part 4: Fecundity model
    
    for (l in 1:nfec) {
    
      year[l] <- year_index[l]-1
      density_fec[l] <- density[pop_index[l],year[l]]
      
      Njuv2[l] <- exp(log(Nad[pop_index[l],year_index[l]]) + 
                      theta +
                      zeta*(density_fec[l]-1) +
                      zeta_l[1]*w1_fec[l] + zeta_q[1]*(w1_fec[l]^2) + 
                      zeta_l[2]*w2_fec[l] + zeta_q[2]*(w2_fec[l]^2) +
                      zeta_l[3]*w3_fec[l] + zeta_q[3]*(w3_fec[l]^2) +
                      zeta_l[4]*w4_fec[l] + zeta_q[4]*(w4_fec[l]^2) + 
                      zeta_l[5]*w5_fec[l] + zeta_q[5]*(w5_fec[l]^2) +
                      zeta_l[6]*w6_fec[l] + zeta_q[6]*(w6_fec[l]^2) +
                      zeta_l[7]*w7_fec[l] + zeta_q[7]*(w7_fec[l]^2) +
                      omega[year_index[l]])
      
      w[l] ~ dbern(psi)
      Nobs_juv2[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
      Nobs_juv_pred[l] ~ dpois(max(py_j[pop_index[l],year_index[l]]*Njuv2[l]*w[l], 0.00000001))
    
    }#l
    
    }", fill = TRUE)
sink()
