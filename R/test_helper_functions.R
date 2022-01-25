# helper functions for test.R


# custom testhat expect function
# mean absolute difference < tol
expect_diffLessTol <- function(object, expected, tolerance = 1e-14, 
                               label = NULL, expected.label = NULL) {
  
  act <- testthat::quasi_label(rlang::enquo(object), label, arg = "object")
  exp <- testthat::quasi_label(rlang::enquo(expected), expected.label, arg = "expected")
  
  abs_diff <- abs(unlist(act$val) - unlist(exp$val))
  
  
  
  testthat::expect(
    all(abs_diff < tolerance),
    sprintf("%s not equal to %s.\nMean absolute difference: %s \nMax. absolute difference: %s", 
            act$lab, exp$lab, mean(abs_diff), max(abs_diff))
  )
  invisible(act$val)
}


### simulate random probabilistic knowledge structure:
# |Q| = 25
# |K| = 500
# beta_q ~ unif(0, .1)
# eta_q  ~ unif(0, .1)
# P(K)   ~ unif(0, 1) -> normalized

genPKS <- function(kQ = 25, kK = 500, seed = 42, coef_gen = runif, ...){
  # cat(paste0("Generating new KS with seed: ", seed, "\n")) 
  set.seed(seed)
  # power set of Q:
  potQ <- expand.grid(rep(list(0:1), kQ), KEEP.OUT.ATTRS = F)
  # 2^ == nrow(potQ)                  # check size

  # sample Elements of potQd that are part of K, without {} and Q  
  idx <- sort(sample(2:(nrow(potQ) - 1), size = kK - 2, replace = FALSE)) 

  # add {} and Q  
  K <- potQ[c(1, idx, nrow(potQ)), ]
  rownames(K) <- 1:kK

  # Parameter beta and eta
  beta <- coef_gen(kQ, ...)
  eta  <- coef_gen(kQ, ...)

  # P(K)
  P.K <- runif(kK, 0, 1)
  P.K <- P.K / sum(P.K)
  sum(P.K)

  true_para <- list(K = as.matrix(K), 
    P.K = P.K, 
    beta = beta, 
    eta = eta, 
    nitems = kQ, 
    nstates = kK)

  return(true_para)
}



# Description: Funcion to simulate (missing) data 
#              Possible dependencies between missings and M:
#                   - MCAR: P(Mq = 1 | q in K) = P(Mq = 1 | q not in K) = Pmcar 
#                   - MNAR: P(Mq = 1 | q in K) = MUq; 
#                           P(Mq = 1 | q not in K) = MUq_


simData <- function(N, true_model, seed = 42){
  set.seed(seed)
  K    <- true_model$K
  P.K  <- true_model$P.K
  beta <- true_model$beta
  eta  <- true_model$eta

  if(is.null(true_model$mu0)){
    mu0 <- mu1 <- rep(0, length(eta)) 
  } else {
    mu0  <- true_model$mu0
    mu1  <- true_model$mu1
  }
  nitems <- ncol(K)
  nstates <- nrow(K)

  # check input
  msg <- c("Number of", "parameters is incorrect or columns of K do not represent items.")
  if(length(beta) != nitems)
    stop(paste(msg[1], "beta", msg[2]))
  if(length(eta) != nitems)
    stop(paste(msg[1], "eta", msg[2]))
  if(length(mu0) != nitems)
    stop(paste(msg[1], "mu0", msg[2]))
  if(length(mu1) != nitems)
    stop(paste(msg[1], "mu1", msg[2]))       
  if(length(P.K) != nstates)
    stop(paste(msg[1], "P.K parameters are incorrect.")) 


  ## R* = RS complete response pattern ~ BlIM

  # sample N times k in K with prob. P.K
  idx_k <- sample(1:nrow(K), size = N, prob = P.K, replace = TRUE)

  # prob. to give correct answer    
  p_q_in_R <- t(t(K[idx_k, ]) * (1 - beta) + t(1 - K[idx_k, ]) * eta)

  # apply prob. and get complete response pattern
  RS <- t(apply(p_q_in_R, 1, function(r) rbinom(length(r), 1, t(r))))

  R <- RS

  # decimation mechanism, leads to missings
  for(i in 1:nrow(R)){
    k <- as.numeric(K[idx_k[i], ])
    p_miss <- k * mu1 + (1 - k) * mu0
    idx_na <- rbinom(ncol(R), 1, p_miss)
    R[i, c(which(idx_na == 1))] <- NA
  }

  return(list(R = R, idx_k = idx_k, RS = RS, tModel = true_model))
}
