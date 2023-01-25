sink("jags/dtms_constant_all.jags")
cat("
model{
  # Priors
  for(i in 1:2){
    h[i] ~ dgamma(1, 4)
    phi[i] <- exp(-h[i] * l.occasion)
    
    mu[i] ~ dgamma(1, 4)
    p[i] <- 1 - exp(-mu[i] * l.occasion)
  }

  eta12 ~ dgamma(1, 4)
  psi12 <- 1 - exp(-eta12 * l.occasion)
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
}
    ", fill = TRUE)
sink()
