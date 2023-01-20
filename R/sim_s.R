sim_s <- function(n.ind, n.occ, h, k, g, lambda){
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
  
  dets <- states <- vector(mode = "list", length = dim(s)[1])
  
  s1 <- apply(s, 1, function(x) sum(x == 1))
  s2 <- apply(s, 1, function(x) sum(x == 2))
  s3 <- apply(s, 1, function(x) sum(x == 3))

  for(i in 1:dim(s)[1]){
    det1 <- rpois(1, lambda[1] * s1[i])
    # det.tmp <- rpois(1, lambda[s[i,1]])
    # p <- sort(runif(det.tmp))
    p <- sort(runif(det1, 0, s1[i]))
    dets[[i]] <- p
    states[[i]] <- rep(1, det1)
    if(s2[i] > 0){
      det2 <- rpois(1, lambda[2] * s2[i])
      # det.tmp <- rpois(1, lambda[s[i,1]])
      # p <- sort(runif(det.tmp))
      p <- sort(runif(det2, s1[i], (s2[i] + s1[i])))
      dets[[i]] <- c(dets[[i]], p)
      states[[i]] <- c(states[[i]], rep(2, det2))
    }
    if(s3[i] > 0){
      det3 <- rpois(1, lambda[3] * s3[i])
      # det.tmp <- rpois(1, lambda[s[i,1]])
      # p <- sort(runif(det.tmp))
      p <- sort(runif(det3, (s2[i] + s1[i]), (s3[i] + s2[i] + s1[i])))
      dets[[i]] <- c(dets[[i]], p)
      states[[i]] <- c(states[[i]], rep(3, det3))
    }
    # for(t in 2:dim(s)[2]){
    #   det.tmp <- rpois(1, lambda[s[i, t]])
    #   p <- sort(runif(det.tmp))
    #   if(length(p) > 0) p <- p + (t - 1)
    #   dets[[i]] <- c(dets[[i]], p)
    #   states[[i]] <- c(states[[i]], rep(s[i,t], det.tmp))
    # }
  }
  
  n_dets <- unlist(lapply(dets, length))
  #min(n_dets)
  det <- state  <- matrix(0, nrow = dim(s)[1], ncol = max(n_dets))
  
  for(i in 1:dim(s)[1]){
    if(n_dets[i] > 0){
      det[i, 1:n_dets[i]] <- dets[[i]] #+ 1
      state[i, 1:n_dets[i]] <- states[[i]]
    }
  }
  
  d <- ceiling(det)
  
  delta <- matrix(0, nrow = n.ind, ncol = max(n_dets) + 1)
  delta[,1] <- det[,1]
  
  for(i in 1:n.ind){
    if(n_dets[i] > 1){
      for(j in 2:n_dets[i]){
        delta[i, j] <- det[i, j] - det[i, j - 1]
      }
    }
    delta[i, n_dets[i] + 1] <- n.occ - max(det[i,])
  }
  
  delta_occ <- delta_lower <-  delta_upper <- matrix(0, nrow = n.ind, ncol = max(n_dets))
  
  for(i in 1:n.ind){
    if(n_dets[i] > 1){
      delta_occ[i, 1] <- d[i, 1] - 1
      delta_lower[i, 1] <- det[i, 1] - d[i, 1] + 1
      delta_upper[i, 1] <- d[i, 1] - det[i, 1]
      for(j in 2:n_dets[i]){
        delta_occ[i, j] <- d[i, j] - d[i, j - 1]
        delta_lower[i, j] <- det[i, j] - d[i, j] + 1
        delta_upper[i, j] <- d[i, j] - det[i, j]
      }
    }
  }
  
  n_dets_plus_one <- n_dets + 1
  
  data <- list(n_ind = dim(s)[1],
               n_occasions = n.occasions, day = day,
               d = d, state = state, s = s,
               max_det = max(n_dets), max_det_plus_1 = max(n_dets_plus_one),
               n_dets = n_dets, n_dets_plus_1 = n_dets_plus_one,
               delta = delta, delta_occ = delta_occ, 
               delta_lower = delta_lower, delta_upper = delta_upper,
               f = matrix(c(1, rep(0, max(s) - 1)), 1, max(s)),
               ones = matrix(rep(1, max(s)), max(s), 1))
  return(data)
}

