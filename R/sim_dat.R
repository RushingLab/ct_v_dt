sim_dat_constant <- function(N, T, h, eta, lambda){
  ## Define transition rate matrix and detection matrix
  Q <- matrix(0, 3, 3)
  Q[1, 1] <- -(eta + h[1])
  Q[1, 2] <- eta
  Q[1, 3] <- h[1]
  
  Q[2, 2] <- -h[2]
  Q[2, 3] <- h[2]
  
  
  ## Simulate true states
  s <- matrix(NA, nrow = N, ncol = T)
  s[,1] <- 1
  
  for(i in 1:N){
    for(t in 2:T){
      s[i, t] <- which(rmultinom(1, 1, prob = expm(Q)[s[i, t - 1],]) == 1)
    }
  }
  
  # Change state 3 to 0 (dead)
  s[s == 3] <- 0
  
  ## Simulate encounter histories as a list because 
  ##  history length differs among individuals
  det_list <- state_list <- vector(mode = "list", length = N)
  U1 <- U2 <- vector(length = N)
  
  # How many days was each individual in each state?
  s1 <- apply(s, 1, function(x) sum(x == 1))
  s2 <- apply(s, 1, function(x) sum(x == 2))
  
  for(i in 1:N){
    U1[i] <- rpois(1, lambda[1] * s1[i]) # Number of detections
    dets1 <- sort(runif(U1[i], 0, s1[i])) # Time of detections 
    state1 <- rep(1, U1[i])
    
    U2[i] <- rpois(1, lambda[2] * s2[i]) # Number of detections
    dets2 <- sort(runif(U2[i], 0, s2[i])) + s1[i] # Time of detections 
    state2 <- rep(2, U2[i])
    
    det_list[[i]] <- c(dets1, dets2)
    state_list[[i]] <- c(state1, state2)
  }
  # Total number of detections for each individual
  U <- U1 + U2
  
  # Covert detections to matrix 
  det <- state <- matrix(0, nrow = N, ncol = max(U))
  
  for(i in 1:N){
    if(U[i] > 0){
      det[i, 1:U[i]] <- det_list[[i]]
      state[i, 1:U[i]] <- state_list[[i]]
    }
  }
  
  # Calculate time difference between detections (including last detection to end of study)
  delta <- matrix(0, nrow = N, ncol = max(U) + 1)
  delta[,1] <- det[,1]
  
  for(i in 1:N){
    if(U[i] > 1){
      for(j in 2:U[i]){
        delta[i, j] <- det[i, j] - det[i, j - 1]
      }
    }
    delta[i, U[i] + 1] <- T - max(det[i,])
  }
  
  det2 <- which(U > 1)
  det1 <- which(U == 1)
  det0 <- which(U == 0)
  
  state2 <- rbind(state[det2,], state[det1,], state[det0,])
  delta2 <-  rbind(delta[det2,], delta[det1,], delta[det0,])
  U2 <- c(U[det2], U[det1], U[det0])
  N2 <- length(det2)
  N1 <- length(det1)
  
  
  # Additional objects needed for JAGS
  zeros <- matrix(0, nrow = nrow(state), ncol = max(U2) + 1)
  
  dat <- list(N = N,
              N2 = N2,
              N1 = N1,
              T = T,
              det = det,
              delta = delta2,
              state = state2,
              U = U2,
              U1 = U,
              zeros = zeros)
  return(dat)
}