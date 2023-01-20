sim_msHawkes <- function(n.ind, n.occ, h, k, g, lambda0, alpha0, beta0){
  n.states <- length(h)
  day <- seq(from = 0, to = 10, length.out = n.occ)
  
  nu <- matrix(0, nrow = length(day), ncol = n.states)
  
  for(i in 1:(n.states - 1)){
    nu[,i] <- (k[i] * day ^ (k[i] - 1)/g[i] ^ k[i])
  }
  
  q <- array(0, dim = c(n.states + 1, n.states + 1, n.occ))
  
  for(i in 1:(n.states)){
    q[i, i + 1, ] <- nu[,i]
    q[i, n.states + 1, ] <- h[i]
    q[i,i, ] <- -apply(q[i, 1:(n.states + 1),], 2, sum)
  }
  
  s <- matrix(NA, nrow = n.ind, ncol = n.occ)
  s[,1] <- 1
  
  for(i in 1:n.ind){
    for(t in 2:n.occ){
      s[i, t] <- which(rmultinom(1, 1, prob = expm::expm(q[,,t])[s[i, t - 1],]) == 1)
    }
  }
  
  s[s == 4] <- 0
  dets <- states <- vector(mode = "list", length = dim(s)[1])
  
  for(i in 1:dim(s)[1]){
    l <- sum(s[i,] == 1)
    dets[[i]] <- hawkes::simulateHawkes(lambda0 = lambda0[1], alpha = alpha0[1], beta = beta0[1], horizon = l)[[1]]
    states[[i]] <- rep(1, length(dets[[i]]))
    for(j in 2:max(s[i, ])){
      l <- sum(s[i,] == j)
      if(l > 1){
        dets.tmp <- hawkes::simulateHawkes(lambda0 = lambda0[j], alpha = alpha0[j], beta = beta0[j], horizon = l)[[1]]
        dets.tmp2 <- dets.tmp + min(which(s[i,] == j) - 1)
        dets[[i]] <- c(dets[[i]], dets.tmp2)
        states[[i]] <- c(states[[i]], rep(j, length(dets.tmp)))
      }
    }
  }
  
  n_dets <- unlist(lapply(dets, length))
  #min(n_dets)
 
  
  ind_det <- which(n_dets > 0)
  ind_nodet <- which(n_dets == 0)
  n_ind_det <- length(ind_det)
  
  ind2 <- c(ind_det, ind_nodet)
  n_dets2 <- c(n_dets[ind_det], n_dets[ind_nodet])
  n_dets_plus_one <- n_dets2 + 1
  
  det <- state  <- matrix(0, nrow = dim(s)[1], ncol = max(n_dets2))
  for(i in 1:n.ind){
    if(n_dets2[i] > 0){
      det[i, 1:n_dets2[i]] <- dets[[ind2[i]]] #+ 1
      state[i, 1:n_dets2[i]] <- states[[ind2[i]]]
    }
  }
  
  # det <- cbind(0, det)
  d <- ceiling(det)
  # d[,1] <- 1
  
  Delta <- Delta_occ <- matrix(0, nrow = dim(s)[1], ncol = max(n_dets2) - 1)
  
  
  for(j in 2:max(n_dets2)) Delta[, j - 1] <- det[,j] - det[,j - 1]
  Delta[Delta < 0] <- 0
  for(j in 2:max(n_dets)) Delta_occ[,j-1] <- d[,j] - d[,j-1] 
  Delta_occ[Delta_occ < 0] <- 0
  
  det.Deltal <- det - d + 1
  det.Deltau <- d - det
  
  # state <- cbind(1, state)
  
  is_prev_same <- as.matrix(t(apply(state, 1, function(x) as.numeric(dplyr::lead(x) == x))))
  is_prev_same <- is_prev_same[, -max(n_dets2)]
  # is_prev_same[,1] <- 0
  # diff <- matrix(0, max(n_dets2), n_ind_det)# n.ind x max1
  
  # for(i in 1:n_ind_det){
  #     for(j in 1:n_dets2[i]){
  #       events <- (det[i, j + 1] - det[i, 1:j])
  #       events <- events[state[i, 1:j] == state[i,(j + 1)]]
  #       if(length(events) == 0) events <- 0
  #       
  #       diff[i, j, 1:length(events)] <- events
  #     }
  # }
  # 
  # length_hist <- apply(diff, c(1, 2), function(x) sum(x != 0))
  
  # length_hist[length_hist == 0] <- 1
  
  # diff[diff == 0] <- 1000
  
  data <- list(n_states = n.states, n_states_plus_1 = n.states + 1,
               n_ind = dim(s)[1], 
               # diff = diff, length_hist = length_hist,
               is_prev_same = is_prev_same,
               n_ind_det = n_ind_det,
               n_occasions = n.occasions, day = day,
               d = d, state = state, s = s, det = det,
               max_det = max(n_dets2), max_det_plus_1 = max(n_dets_plus_one) ,
               max_det_m_1 = max(n_dets2) - 1,
               n_dets = n_dets2, n_dets_plus_1 = n_dets2 + 1,
               Delta = Delta, Delta_occ = Delta_occ, 
               Delta_lower = det.Deltal, Delta_upper = det.Deltau,
               f = matrix(c(1, rep(0, n.states)), 1, n.states + 1),
               ones = matrix(rep(1, n.states + 1), n.states + 1, 1))
  return(data)
}
