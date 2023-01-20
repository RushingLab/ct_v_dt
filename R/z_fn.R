ct_to_dt <- function(det, U, state, T, l.occasion, conflict_rule){
  ## Convert to discrete encounter histories
  occasions <- seq(from = 0, to = T, by = l.occasion)
  n.occasions <- length(occasions) - 1
  N <- dim(det)[1]
  
  ch <- matrix(0, nrow = N, ncol = n.occasions)
  
  for(i in 1:N){
    if(U[i] > 0){
      dets <- data.frame(occasion = as.numeric(cut(det[i,1:U[i]], breaks = occasions)),
                         state = state[i, 1:U[i]])
      mult_det <- unique(dets$occasion[duplicated(dets$occasion)])
      for(j in 1:length(mult_det)){
        obs_states <- unique(dets$state[dets$occasion == mult_det[j]])
        conflict <- ifelse(length(obs_states) > 1, 1, 0)
        if(conflict){
          if(conflict_rule == "min") dets$state[dets$occasion == mult_det[j]] <- min(obs_states)
          if(conflict_rule == "max") dets$state[dets$occasion == mult_det[j]] <- max(obs_states)
          if(conflict_rule == "sample") dets$state[dets$occasion == mult_det[j]] <- sample(obs_states, size = 1)
        }
      }
      dets <- dets[!duplicated(dets$occasion),]
      ch[i, dets$occasion] <-  dets$state
    }
  }
  ch[, 1] <- 1
  ch[ch == 0] <- 3
  
  z.known <- ch
  z.known[z.known == 3] <- NA
  for (i in 1:N){
    m <- min(which(!is.na(z.known[i,])))
    z.known[i, m] <- NA
  }
  
  f <- rep(1, N)
  dat <- list(N = N,
              n.occasions = n.occasions,
              l.occasion = l.occasion,
              y = ch,
              f = f,
              z = z.known)
  return(dat)
}

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
  ch2 <- ch
  for (i in 1:dim(ch)[1]){
    d <- which(ch[i,] != 3)
    nd <- which(ch[i, ] == 3)
    ch2[i, nd] <- NA
    
    for(j in nd){
      m <- max(ch2[i, 1:(j-1)], na.rm = TRUE)
      ch[i, j] <- m
    }
    ch[i, d] <- NA
  }
  return(ch)
}