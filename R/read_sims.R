# Function to read in posterior means from simulations

read_sims <- function(nSims, h, g, k, lambda, dir, tidy = TRUE){

  h.hat <- matrix(NA, nrow = nSims, ncol = length(h))
  g.hat <- matrix(NA, nrow = nSims, ncol = length(g))
  k.hat <- matrix(NA, nrow = nSims, ncol = length(k))
  lambda.hat <- matrix(NA, nrow = nSims, ncol = length(lambda))
  
  for(s in 1:nSims){
    f <- here::here(paste0("output/", dir,"/sim_", s, "/summary.txt"))
    for(i in 1:ncol(h.hat)){
      temp <- grep(paste0("h\\[", i, "\\]"), readLines(f), value = TRUE)
      h.hat[s, i] <- as.numeric(substring(temp, 12, last = 20))
    }
    
    for(i in 1:ncol(g.hat)){
      temp <- grep(paste0("g\\[", i, "\\]"), readLines(f), value = TRUE)
      g.hat[s, i] <- as.numeric(substring(temp, 12, last = 20))
    }
    
    for(i in 1:ncol(k.hat)){
      temp <- grep(paste0("k\\[", i, "\\]"), readLines(f), value = TRUE)
      k.hat[s, i] <- as.numeric(substring(temp, 12, last = 20))
    }
    
    for(i in 1:ncol(lambda.hat)){
      temp <- grep(paste0("lambda\\[", i, "\\]"), readLines(f), value = TRUE)
      lambda.hat[s, i] <- as.numeric(substring(temp, 12, last = 20))
    }
  }
  
  if(tidy){
    out <- data.frame(value = c(c(h.hat), c(g.hat), c(k.hat), c(lambda.hat)),
                      parameter = c(rep("h", length(h) * nSims), 
                                    rep("g", length(g) * nSims),
                                    rep("k", length(k) * nSims),
                                    rep("lambda", length(lambda) * nSims)),
                      state = c(rep(1:length(h), each = nSims), 
                                rep(1:length(g), each = nSims),
                                rep(1:length(k), each = nSims),
                                rep(1:length(lambda), each = nSims)))
  }else{
    out <- list(h = h.hat, g = g.hat, k = k.hat, lambda = lambda.hat)
  }

  
  return(out)
}


