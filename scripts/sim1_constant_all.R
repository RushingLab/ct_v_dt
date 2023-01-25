library(rjags)
library(expm)
load.module("bugs")
load.module("msm")
source("R/z_fn.R")
source("R/sim_dat.R")


#### Set parameters ----
nSims <- 100              # Number of simulations
N <- 100                  # Number of individuals
T <- 100                  # Length of study
h <- c(0.0025, 0.01)      # Hazard rates
eta <- 0.005              # Transition rate from state 1 to state 2    
lambda <- c(0.2, 0.1, 0)  # Detection rates 


#### MCMC settings ----
params <- c("h", "mu", "eta12")
nC <- 3
nI <- 2500
nB <- 500
nT <- 1

for(ii in 1:nSims){
  print(paste0("Simulation ", ii))
  #### Simulate data ----
  ctms_data <- sim_dat_constant(N, T, h, eta, lambda)
  
  #### Fit continuous-time model ----
  ctms_inits <- function(){list(h = runif(2, 0, 0.015),
                                eta12 = runif(1, 0, 0.25),
                                mu = runif(2, 0, 0.25))}
  print(paste0("CT ", ii))
  ctms_mod <- jagsUI::jags(data = ctms_data, inits = ctms_inits, modules = "msm",
                           parameters.to.save = params, 
                           model.file = "jags/ctms_constant_all.jags", 
                           n.chains = nC, n.iter = nI, 
                           n.burnin = nB, n.thin = nT, parallel = TRUE, verbose = FALSE)
  
  ctms_summ <- data.frame(Model = "CT", l = 1, Elapsed.time = ctms_mod$mcmc.info$elapsed.mins,
                          sim = ii,
                          Param = row.names(ctms_mod$summary)[1:5],
                          Est = ctms_mod$summary[1:5, 1],
                          LCI = ctms_mod$summary[1:5, 3],
                          UCI = ctms_mod$summary[1:5, 7],
                          width = ctms_mod$summary[1:5, 7] - ctms_mod$summary[1:5, 3],
                          n.eff = ctms_mod$summary[1:5, 9])
  
  if(ii == 1){
    ctms_Summ <- ctms_summ
  }else{
    ctms_Summ <- dplyr::bind_rows(ctms_Summ, ctms_summ)
  }
  
  #### Fit discrete-time model, l = 1 ----
  l.occasion <- 1
  conflict_rule <- "min"
  dtms_data <- ct_to_dt(det = ctms_data$det, U = ctms_data$U1, 
                        state = ctms_data$state, 
                        T, l.occasion, conflict_rule)

  dtms_inits <- function(){list(h = runif(2, 0, 0.015),
                                eta12 = runif(1, 0, 0.25),
                                mu = runif(2, 0, 0.25),
                                z = ms.init.z(dtms_data$y, dtms_data$f))}
  
  print(paste0("DT ", l.occasion, " ", ii))
  dtms_mod <- jagsUI::jags(data = dtms_data, inits = dtms_inits,
                           parameters.to.save = params, 
                           model.file = "jags/dtms_constant_all.jags", 
                           n.chains = nC, n.iter = nI, 
                           n.burnin = nB, n.thin = nT, parallel = TRUE, verbose = FALSE)
  
  dtms_summ1 <- data.frame(Model = "DT", l = l.occasion, Elapsed.time = dtms_mod$mcmc.info$elapsed.mins,
                          sim = ii,
                          Param = row.names(dtms_mod$summary)[1:5],
                          Est = dtms_mod$summary[1:5, 1],
                          LCI = dtms_mod$summary[1:5, 3],
                          UCI = dtms_mod$summary[1:5, 7],
                          width = dtms_mod$summary[1:5, 7] - dtms_mod$summary[1:5, 3],
                          n.eff = dtms_mod$summary[1:5, 9])
  
  
  if(ii == 1){
    dtms_Summ1 <- dtms_summ1
  }else{
    dtms_Summ1 <- dplyr::bind_rows(dtms_Summ1, dtms_summ1)
  }
  
  #### Fit discrete-time model, l = 2 ----
  l.occasion <- 2
  conflict_rule <- "min"
  dtms_data <- ct_to_dt(det = ctms_data$det, U = ctms_data$U1, 
                        state = ctms_data$state, 
                        T, l.occasion, conflict_rule)
  
  
  print(paste0("DT ", l.occasion, " ", ii))
  dtms_mod <- jagsUI::jags(data = dtms_data, inits = dtms_inits,
                           parameters.to.save = params, 
                           model.file = "jags/dtms_constant_all.jags", 
                           n.chains = nC, n.iter = nI, 
                           n.burnin = nB, n.thin = nT, parallel = TRUE, verbose = FALSE)
  
  dtms_summ2 <- data.frame(Model = "DT", l = l.occasion, Elapsed.time = dtms_mod$mcmc.info$elapsed.mins,
                           sim = ii,
                           Param = row.names(dtms_mod$summary)[1:5],
                           Est = dtms_mod$summary[1:5, 1],
                           LCI = dtms_mod$summary[1:5, 3],
                           UCI = dtms_mod$summary[1:5, 7],
                           width = dtms_mod$summary[1:5, 7] - dtms_mod$summary[1:5, 3],
                           n.eff = dtms_mod$summary[1:5, 9])
  
  
  if(ii == 1){
    dtms_Summ2 <- dtms_summ2
  }else{
    dtms_Summ2 <- dplyr::bind_rows(dtms_Summ2, dtms_summ2)
  }
  
  #### Fit discrete-time model, l = 5 ----
  l.occasion <- 5
  conflict_rule <- "min"
  dtms_data <- ct_to_dt(det = ctms_data$det, U = ctms_data$U1, 
                        state = ctms_data$state, 
                        T, l.occasion, conflict_rule)
  
  
  print(paste0("DT ", l.occasion, " ", ii))
  dtms_mod <- jagsUI::jags(data = dtms_data, inits = dtms_inits,
                           parameters.to.save = params, 
                           model.file = "jags/dtms_constant_all.jags", 
                           n.chains = nC, n.iter = nI, 
                           n.burnin = nB, n.thin = nT, parallel = TRUE, verbose = FALSE)
  
  dtms_summ5 <- data.frame(Model = "DT", l = l.occasion, Elapsed.time = dtms_mod$mcmc.info$elapsed.mins,
                           sim = ii,
                           Param = row.names(dtms_mod$summary)[1:5],
                           Est = dtms_mod$summary[1:5, 1],
                           LCI = dtms_mod$summary[1:5, 3],
                           UCI = dtms_mod$summary[1:5, 7],
                           width = dtms_mod$summary[1:5, 7] - dtms_mod$summary[1:5, 3],
                           n.eff = dtms_mod$summary[1:5, 9])
  
  
  if(ii == 1){
    dtms_Summ5 <- dtms_summ5
  }else{
    dtms_Summ5 <- dplyr::bind_rows(dtms_Summ5, dtms_summ5)
  }
  
  
  #### Fit discrete-time model, l = 10 ----
  l.occasion <- 10
  conflict_rule <- "min"
  dtms_data <- ct_to_dt(det = ctms_data$det, U = ctms_data$U1, 
                        state = ctms_data$state, 
                        T, l.occasion, conflict_rule)
  
  
  print(paste0("DT ", l.occasion, " ", ii))
  dtms_mod <- jagsUI::jags(data = dtms_data, inits = dtms_inits,
                           parameters.to.save = params, 
                           model.file = "jags/dtms_constant_all.jags", 
                           n.chains = nC, n.iter = nI, 
                           n.burnin = nB, n.thin = nT, parallel = TRUE, verbose = FALSE)
  
  dtms_summ10 <- data.frame(Model = "DT", l = l.occasion, Elapsed.time = dtms_mod$mcmc.info$elapsed.mins,
                           sim = ii,
                           Param = row.names(dtms_mod$summary)[1:5],
                           Est = dtms_mod$summary[1:5, 1],
                           LCI = dtms_mod$summary[1:5, 3],
                           UCI = dtms_mod$summary[1:5, 7],
                           width = dtms_mod$summary[1:5, 7] - dtms_mod$summary[1:5, 3],
                           n.eff = dtms_mod$summary[1:5, 9])
  
  
  if(ii == 1){
    dtms_Summ10 <- dtms_summ10
  }else{
    dtms_Summ10 <- dplyr::bind_rows(dtms_Summ10, dtms_summ10)
  }
  
  #### Fit discrete-time model, l = 20 ----
  l.occasion <- 20
  conflict_rule <- "min"
  dtms_data <- ct_to_dt(det = ctms_data$det, U = ctms_data$U1, 
                        state = ctms_data$state, 
                        T, l.occasion, conflict_rule)
  
  
  print(paste0("DT ", l.occasion, " ", ii))
  dtms_mod <- jagsUI::jags(data = dtms_data, inits = dtms_inits,
                           parameters.to.save = params, 
                           model.file = "jags/dtms_constant_all.jags", 
                           n.chains = nC, n.iter = nI, 
                           n.burnin = nB, n.thin = nT, parallel = TRUE, verbose = FALSE)
  
  dtms_summ20 <- data.frame(Model = "DT", l = l.occasion, Elapsed.time = dtms_mod$mcmc.info$elapsed.mins,
                           sim = ii,
                           Param = row.names(dtms_mod$summary)[1:5],
                           Est = dtms_mod$summary[1:5, 1],
                           LCI = dtms_mod$summary[1:5, 3],
                           UCI = dtms_mod$summary[1:5, 7],
                           width = dtms_mod$summary[1:5, 7] - dtms_mod$summary[1:5, 3],
                           n.eff = dtms_mod$summary[1:5, 9])
  
  
  if(ii == 1){
    dtms_Summ20 <- dtms_summ20
  }else{
    dtms_Summ20 <- dplyr::bind_rows(dtms_Summ20, dtms_summ20)
  }
  
  Summ <- dplyr::bind_rows(ctms_Summ, 
                           dtms_Summ1, 
                           dtms_Summ2, 
                           dtms_Summ5, 
                           dtms_Summ10,
                           dtms_Summ20)
  
  saveRDS(Summ, "output/sim1_constant_all.RDS")
  ggplot() + 
    geom_hline(data = true_df, aes(yintercept = value), linetype = "dashed", color = "grey50") +
    geom_boxplot(data = Summ, aes(x = as.factor(l), y = Est, color = Model)) +
    facet_wrap(~Param, scale = "free")
  
  ggplot(Summ, aes(x = as.factor(l), y = Elapsed.time, color = Model)) + geom_boxplot()
}



  true_df <- data.frame(Param = unique(Summ$Param),
                      value = c(h, lambda[1:2], eta))
library(ggplot2)

 