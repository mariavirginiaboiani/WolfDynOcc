################################################################################
##### ---------------------- ALPINE WOLF OCCUPANCY ----------------------- #####
##### ------------- EFFECT OF COV on RHO and PHI by INTERVAL ------------- #####
##### ----------- psi ~ forest, herb, pas, ung, pop, lightpoll ----------- #####
##### ----------- rh0 ~ forest, herb, pas, ung, pop, lightpoll ----------- #####
##### ----------- phi ~ forest, herb, pas, ung, pop, lightpoll ----------- #####
##### ------------------------  p ~ effort, snowfall---------------------- #####
################################################################################

## ------ IMPORT REQUIRED LIBRARIES ------

library(coda)
library(nimble)

load("inputWolfDynOcc.RData")


## ------  NIMBLE MODEL ------

modelCode <- nimbleCode({
  
  ##-- ECOLOGICAL SUB-MODEL
  incProbRJ ~ dbeta(10,10) # inclusion prob
  
  ##-- Initial occupancy probability
  ## Intercept 
  psi_alpha <- logit(mean.psi)
  mean.psi ~ dunif(0,1)
  
  ## Spatial fixed effects 
  psi_beta_forest ~ dnorm(mean = 0, sd = 0.2)
  psi_beta_herb ~ dnorm(mean = 0, sd = 0.2)
  psi_beta_ung ~ dnorm(mean = 0, sd = 0.2)
  psi_beta_pas ~ dnorm(mean = 0, sd = 0.2)
  psi_beta_pop ~ dnorm(mean = 0, sd = 0.2)
  psi_beta_lightpoll ~ dnorm(mean = 0, sd = 0.2) 
  psi_beta_iucn ~ dnorm(mean = 0, sd = 0.2) 
  ## Spatial fixed effects (RJCMC)
  z_psi_beta_forest ~ dbern(incProbRJ)
  z_psi_beta_herb ~ dbern(incProbRJ)
  z_psi_beta_ung ~ dbern(incProbRJ)
  z_psi_beta_pas ~ dbern(incProbRJ)
  z_psi_beta_pop ~ dbern(incProbRJ)
  z_psi_beta_lightpoll ~ dbern(incProbRJ)
  z_psi_beta_iucn ~ dbern(incProbRJ)
  
  ## Initial occupancy
  for(i in 1:n.sites){
    logit(psi[i]) <-  psi_alpha +
      psi_beta_forest  * z_psi_beta_forest * forest[i] +
      psi_beta_herb  * z_psi_beta_herb * herb[i] +
      psi_beta_ung * z_psi_beta_ung * ung[i] +
      psi_beta_pas * z_psi_beta_pas * pas[i] +
      psi_beta_pop * z_psi_beta_pop * pop[i] +
      psi_beta_lightpoll * z_psi_beta_lightpoll * lightpoll[i] +
      psi_beta_iucn * z_psi_beta_iucn * iucn[i] 
  }#i
  
  
  ##--- Colonization probabilities
  ## Intercept 
  rho_alpha_0 <- logit(mean.rho)
  mean.rho ~ dunif(0,1)
  
  ## Intercept (TREND)
  rho_alpha_trend ~ dnorm(mean = 0, sd = 0.2)
  ## Intercept (RJCMC)
  z_rho_alpha_trend ~ dbern(incProbRJ)
  
  ## Spatial fixed effects (INTERCEPTS) 
  rho_beta_forest_0 ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_herb_0 ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_ung_0 ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_pas_0 ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_pop_0 ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_lightpoll_0 ~ dnorm(mean = 0, sd = 0.2)
  ## Spatial fixed effects (RJCMC)
  z_rho_beta_forest_0 ~ dbern(incProbRJ)
  z_rho_beta_herb_0 ~ dbern(incProbRJ)
  z_rho_beta_ung_0 ~ dbern(incProbRJ)
  z_rho_beta_pas_0 ~ dbern(incProbRJ)
  z_rho_beta_pop_0 ~ dbern(incProbRJ)
  z_rho_beta_lightpoll_0 ~ dbern(incProbRJ)
  
  ## Spatial fixed effects (TRENDS)
  rho_beta_forest_trend ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_herb_trend ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_ung_trend ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_pas_trend ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_pop_trend ~ dnorm(mean = 0, sd = 0.2)
  rho_beta_lightpoll_trend ~ dnorm(mean = 0, sd = 0.2)
  ## Spatial fixed effects (RJCMC)
  z_rho_beta_forest_trend ~ dbern(incProbRJ)
  z_rho_beta_herb_trend ~ dbern(incProbRJ)
  z_rho_beta_ung_trend ~ dbern(incProbRJ)
  z_rho_beta_pas_trend ~ dbern(incProbRJ)
  z_rho_beta_pop_trend ~ dbern(incProbRJ)
  z_rho_beta_lightpoll_trend ~ dbern(incProbRJ)
  
  ## NO RJ-MCMC CONSTRAINTS !!!!
  
  ## Annual colonization probabilities
  for(t in 1:(n.years-1)){
    
    ## Annual intercept 
    rho_alpha[t] <- rho_alpha_0 + 
      rho_alpha_trend * z_rho_alpha_trend * t 
    
    ## Annual spatial fixed effects 
    rho_beta_forest[t] <- rho_beta_forest_0 * z_rho_beta_forest_0 + rho_beta_forest_trend * z_rho_beta_forest_trend * t
    rho_beta_herb[t] <- rho_beta_herb_0 * z_rho_beta_herb_0 + rho_beta_herb_trend * z_rho_beta_herb_trend * t
    rho_beta_ung[t] <- rho_beta_ung_0 * z_rho_beta_ung_0 + rho_beta_ung_trend * z_rho_beta_ung_trend * t
    rho_beta_pas[t] <- rho_beta_pas_0 * z_rho_beta_pas_0 + rho_beta_pas_trend * z_rho_beta_pas_trend * t
    rho_beta_pop[t] <- rho_beta_pop_0 * z_rho_beta_pop_0 + rho_beta_pop_trend * z_rho_beta_pop_trend * t
    rho_beta_lightpoll[t] <- rho_beta_lightpoll_0 * z_rho_beta_lightpoll_0 + rho_beta_lightpoll_trend * z_rho_beta_lightpoll_trend * t
    
    ## Annual colonization probabilities
    for(i in 1:n.sites){
      logit(rho_annual[i,t]) <- rho_alpha[t] +
        rho_beta_forest[t] * forest[i] +
        rho_beta_herb[t] * herb[i] +
        rho_beta_ung[t] * ung[i] +
        rho_beta_pas[t] * pas[i] +
        rho_beta_pop[t] * pop[i] +
        rho_beta_lightpoll[t] * lightpoll[i]
    }#i
  }#t
  
  
  
  ##--- Survival probabilities
  ## Intercept 
  phi_alpha_0 <- logit(mean.phi)
  mean.phi ~ dunif(0,1)
  
  ## Intercept (TRENDS)
  phi_alpha_trend ~ dnorm(mean = 0, sd = 0.2)
  ## Intercept (RJCMC)
  z_phi_alpha_trend ~ dbern(incProbRJ)
  
  ## Spatial fixed effects (INTERCEPTS) 
  phi_beta_forest_0 ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_herb_0 ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_ung_0 ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_pas_0 ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_pop_0 ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_lightpoll_0 ~ dnorm(mean = 0, sd = 0.2)
  ## Spatial fixed effects (RJCMC)
  z_phi_beta_forest_0 ~ dbern(incProbRJ)
  z_phi_beta_herb_0 ~ dbern(incProbRJ)
  z_phi_beta_ung_0 ~ dbern(incProbRJ)
  z_phi_beta_pas_0 ~ dbern(incProbRJ)
  z_phi_beta_pop_0 ~ dbern(incProbRJ)
  z_phi_beta_lightpoll_0 ~ dbern(incProbRJ)
  
  ## Spatial fixed effects (TRENDS)
  phi_beta_forest_trend ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_herb_trend ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_ung_trend ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_pas_trend ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_pop_trend ~ dnorm(mean = 0, sd = 0.2)
  phi_beta_lightpoll_trend ~ dnorm(mean = 0, sd = 0.2)
  ## Spatial fixed effects (RJCMC)
  z_phi_beta_forest_trend ~ dbern(incProbRJ)
  z_phi_beta_herb_trend ~ dbern(incProbRJ)
  z_phi_beta_ung_trend ~ dbern(incProbRJ)
  z_phi_beta_pas_trend ~ dbern(incProbRJ)
  z_phi_beta_pop_trend ~ dbern(incProbRJ)
  z_phi_beta_lightpoll_trend ~ dbern(incProbRJ)
  
  ## NO RJ-MCMC CONSTRAINTS!!!!
  
  ## Annual colonization probabilities
  for(t in 1:(n.years-1)){
    
    ## Annual intercept 
    phi_alpha[t] <- phi_alpha_0 + 
      phi_alpha_trend * z_phi_alpha_trend * t 
    
    ## Annual spatial fixed effects 
    phi_beta_forest[t] <- phi_beta_forest_0 * z_phi_beta_forest_0 + phi_beta_forest_trend * z_phi_beta_forest_trend * t
    phi_beta_herb[t] <- phi_beta_herb_0 * z_phi_beta_herb_0 + phi_beta_herb_trend * z_phi_beta_herb_trend * t
    phi_beta_ung[t] <- phi_beta_ung_0 * z_phi_beta_ung_0 + phi_beta_ung_trend * z_phi_beta_ung_trend * t
    phi_beta_pas[t] <- phi_beta_pas_0 * z_phi_beta_pas_0 + phi_beta_pas_trend * z_phi_beta_pas_trend * t
    phi_beta_pop[t] <- phi_beta_pop_0 * z_phi_beta_pop_0 + phi_beta_pop_trend * z_phi_beta_pop_trend * t
    phi_beta_lightpoll[t] <- phi_beta_lightpoll_0 * z_phi_beta_lightpoll_0 + phi_beta_lightpoll_trend * z_phi_beta_lightpoll_trend * t
    
    ## Annual colonization probabilities
    for(i in 1:n.sites){
      logit(phi_annual[i,t]) <- phi_alpha[t] +
        phi_beta_forest[t] * forest[i] +
        phi_beta_herb[t] * herb[i] +
        phi_beta_ung[t] * ung[i] +
        phi_beta_pas[t] * pas[i] +
        phi_beta_pop[t] * pop[i] +
        phi_beta_lightpoll[t] * lightpoll[i]
    }#i
  }#t
  
  
  ##-- Likelihood
  for (i in 1:n.sites){
    ## Initial occupancy
    z[i,1] ~ dbern(psi[i])
    
    ##--- Colonization/survival processes
    for (t in 1:(n.years-1)){
      ## Sample next state
      z[i,t+1] ~ dbern(z[i,t]*phi_annual[i,t] + (1-z[i,t])*rho_annual[i,t])
    }#k
  }#i
  
  
  ##-- OBSERVATION SUB-MODEL
  
  ##-- Likelihood
  for(k in 1:n.sessions){ 
    ## Intercept 
    p_alpha[k] <- logit(mean.p[k])
    mean.p[k] ~ dunif(0,1)
    
    ## Spatial fixed effects 
    p_beta_searchEffort[k] ~ dnorm(mean = 0, sd = 0.2)
    p_beta_snowfall[k] ~ dnorm(mean = 0, sd = 0.2)
    
    for(i in 1:n.sites){
      for(j in 1:n.surveys){ 
        logit(p[i,j,k]) <- p_alpha[k] +
          p_beta_searchEffort[k] * searchEffort[i,j,k] +
          p_beta_snowfall[k] * snowfall[i,j,k]
        
        y[i,j,k] ~ dbern(z[i,sess.year[k]] * p[i,j,k] * isMonitored[i,k]) ## Multiply by 0/1 if monitored this year        
      }#j
    }#k
  }#i
  
  
  ##-- Derived parameters
  for(t in 1:n.years){
    occ[t] <- sum(z[1:n.sites,t]) ## Number of "sites" occupied each year
  }#k
}) 
