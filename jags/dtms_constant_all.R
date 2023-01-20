sink("jags/dtms_constant_all.jags")
cat("
model{
  # Priors
  for(i in 1:2){
    phi[i] ~ dbeta(1, 1)
    p[i] ~ dbeta(1, 1)
  }

  psi12 ~ dbeta(1, 1)
  # g ~ dgamma(1, 1)
  # k ~ dgamma(1, 1)

  # State matrix
  ps[1, 1] <- phi[1] * (1 - psi12)
  ps[1, 2] <- phi[1] * psi12
  ps[1, 3] <- 1 - phi[1]
  ps[2, 1] <- 0
  ps[2, 2] <- phi[2]
  ps[2, 3] <- 1 - phi[2]
  ps[3, 1] <- 0
  ps[3, 2] <- 0
  ps[3, 3] <- 1
      
  # Observation matrix
  po[1, 1] <- p[1]
  po[1, 2] <- 0
  po[1, 3] <- 1 - p[1]
  po[2, 1] <- 0
  po[2, 2] <- p[2]
  po[2, 3] <- 1 - p[2]
  po[3, 1] <- 0
  po[3, 2] <- 0
  po[3, 3] <- 1
  
  #for(p in 1:P){
  # Transition rates
    #eta12[p] <- k * c_p[p] ^ (k - 1)/g ^ k
    
  # Define intensity matrix
  #   Q[p, 1, 1] <- -(eta12[p] + h[1])
  #   Q[p, 1, 2] <- eta12[p]
  #   Q[p, 1, 3] <- h[1]
  #   
  #   Q[p, 2, 1] <- 0
  #   Q[p, 2, 2] <- -h[2]
  #   Q[p, 2, 3] <- h[2]
  #   
  #   Q[p, 3, 1] <- 0
  #   Q[p, 3, 2] <- 0
  #   Q[p, 3, 3] <- 0
  #   R[p, 1:3, 1:3] <- mexp(Q[p, 1:3, 1:3])
  # }
  
  # Likelihood 
  for (i in 1:N){
   # Define latent state at first capture
   z[i, f[i]] <- 1
   for (t in (f[i] + 1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i, t] ~ dcat(ps[z[i,t - 1], 1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i, t], 1:3])
      } #t
   } #i
   
   ############## Derived parameters
   # Hazard rates
   h[1] <- -log(phi[1])/l.occasion
   h[2] <- -log(phi[2])/l.occasion
   
   # Detection rates
   mu[1] <- -log(1 - p[1])/l.occasion
   mu[2] <- -log(1 - p[2])/l.occasion
   
   # Transition rate
   eta12 <- -log(1 - psi12)/l.occasion
   
}
    ", fill = TRUE)
sink()