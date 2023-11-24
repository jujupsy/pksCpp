#' Basic Local Independence Models (BLIMs) using C++
#'
#' @description
#' Fits a basic local independence model (BLIM) for probabilistic 
#' knowledge structures by maximum likelihood estimation.
#'
#' @param K a state-by-problem indicator matrix representing the knowledge
#'          structure.  An element is one if the problem is contained in the 
#'          state, and else zero.
#' @param N.R a (named) vector of absolute frequencies of response patterns or 
#'            a person-by-problem (item) matrix containing the the individual answers to each item.
#'            An element is one if the problem is solved by the person, and else zero.
#'            Missing values are coded as NA.
#' @param na.use a character string indicating a method to deal with missing values (see Details).
#' @param tol tolerance, stopping criterion for iteration.
#' @param maxiter the maximum number of iterations.
#' @param fdb feedback during the optimization.
#'
#' @details See Doignon and Falmagne (1999) for details on the basic local independence
#'  model (BLIM) for probabilistic knowledge structures.
#'  
#'  The BLIM cannot handle missing data. In \code{na.use}, however, methods can 
#'  be specified how missing values are to be taken into account. Besides the 
#'  trivial methods, \code{"omit"} (data series with missing values are excluded) or 
#'  \code{"missing-as-wrong"} (missing values are coded as wrong), the missing data can be 
#'  integrated into the model estimation by specifying the missing data 
#'  generating process (decimation mechanism). This leads to two extensions of 
#'  the BLIM, the \code{"IMBLIM"} and \code{"MissBLIM"} (De Chiusole et al. 2015). 
#'  The \code{"IMBLIM"} (Ignorable Missings BLIM) assumes that the decimation mechanism 
#'  is independent of the knowledge state K. This is the case, for example, 
#'  when individuals are presented with only a random subset of items from Q. 
#'  In the \code{"MissBLIM"} (Non-ignorable missings BLIM), dependencies 
#'  between the missing values and the knowledge state K are assumed. 
#'  This means that depending on the mastered items in K, the missing data 
#'  varies. For example, if one assumes that people who have not mastered an 
#'  item q are more likely to give no answer than people who have mastered q.
#'  To capture these dependencies, additional parameters \eqn{\mu = (\mu_q)_{q \in Q}} (\code{mu1})
#'  and \eqn{\bar{\mu} = (\mu_{\bar{q}})_{q \in Q}} (\code{mu0}) are estimated.
#'  Here, \eqn{\mu_q} describes the probability that item q, if contained in K, 
#'  will not be answered, and \eqn{\mu_{\bar{q}}} describes the probability that the 
#'  item will not be answered if not contained in K.
#'  Similar to the response patterns themselves, local stochastic independence 
#'  of the missing values is assumed: It states that for a fixed knowledge 
#'  state K, not answering or answering one item has no influence on the 
#'  response behavior for the other items.
#'  
#'
#' @return An object of class \code{blim}. See  \code{?blim}.
#' 
#' @examples
#' data(DoignonFalmagne7)
#' K   <- DoignonFalmagne7$K         # knowledge structure
#' N.R <- DoignonFalmagne7$N.R       # frequencies of response patterns
#' 
#' ## Fit basic local independence model (BLIM) 
#' blimCpp(K, N.R)
#' @export
#' @useDynLib pksCpp
#' @import RcppEigen
#' @import stats
#' @import sets
#' @import pks
#' @importFrom Rcpp evalCpp
blimCpp <- function(K, N.R, na.use = c("omit", "missing-as-wrong", "IMBLIM", "MissBLIM"), 
                    tol = 1e-07 , maxiter = 10000, fdb = FALSE){
  
  K <- as.matrix(K)
  
  # N.R can be a (named) vector of absolute frequencies or a matrix of response patterns
  if(is.atomic(N.R)){ # (named) vector
    R_all <- as.binmat(N.R, uniq = FALSE)
    
  } else if(is.matrix(N.R)){ # matrix
    R_all <- N.R
    
  } else {
    stop("N.R has to be a (named) vector of absolute frequencies or matrix of repsonse pattern")
  }
  
  # input checks
  if(ncol(K) != ncol(R_all)) stop("Matrix K and R must have the same number of columns")
  
  na.use <- match.arg(na.use)
  if(na.use == "missing-as-wrong"){
    R_all[is.na(R_all)] <- 0
    
  } else if(na.use == "omit"){ 
    R_all <- na.omit(R_all)
    
  }
  
  N.R <- as.pattern(R_all, freq = TRUE)
  R   <- as.binmat(N.R, uniq = TRUE)
  
  
  N       <- sum(N.R)  # sample size
  nitems  <- ncol(K)   # number of items q in Q
  nstates <- nrow(K)   # number of states k in K
  npat    <- nrow(R)   # number of unique response pattern

  ## set initial values
  P.K  <- setNames(rep(1/nstates, nstates),  rownames(K)) # probability vector P(K)
  beta <- setNames(rep(0.1, nitems), colnames(K))         # initital parameter estimations 
  eta  <- setNames(rep(0.1, nitems), colnames(K))
  mu0  <- rep(0.1, nitems)        # mu_q_
  mu1  <- rep(0.1, nitems)        # mu_q

  .M <- (is.na(R))   # M: missing pattern: missing = 1, otherwise 0
  .R <- (R == 1)     # R: right pattern: right = 1, otherwise 0
  .R[is.na(.R)] <- 0 # 
  .W <- (R == 0)     # W: wrong pattern
  .W[is.na(.W)] <- 0
  
  iter <- 1
  maxdiff <- 2 * tol

  # convert to double for c++ function
  mode(.R) <- mode(K) <- mode(N.R) <- mode(P.K) <- mode(beta) <- 
             mode(eta) <- mode(mu0) <- mode(mu1) <- mode(.M) <- mode(.W) <- "double"
  
  # EM Loop 
  #### C++ Funktion waehlen
  if(na.use == "IMBLIM"){
    # IMBLIM
    P.M <- matrix(1, npat, nstates) * 
      (as.pattern(.M * 1, freq = TRUE)[as.pattern(.M * 1)] / N)
    mode(P.M) <- "double"  
    para <- emIMBLIMcpp(.R, .W, .M, K, N.R, P.K, beta, eta, P.M, 
                        maxiter, tol, fdb)
    
  }else if(na.use == "MissBLIM"){
    # MissBLIM
    para <- emMissBLIMcpp(.R, .W, .M, K, N.R, P.K, eta, beta, mu0, mu1, 
                         maxiter, tol, fdb)
    # stop("MissBLIM not implemented yet")
  }else if(na.use == "omit" || na.use == "missing-as-wrong"){
    # BLIM
    para <- emBLIMcpp(.R, K, N.R, P.K, beta, eta, maxiter, tol, fdb)
    
  }
  
  
  # drop dimensions 
  para$P.R <- drop(para$P.R)

  # warning if maxiter was reached
  if(para$maxiter) 
    warning(paste("Maximum number of", maxiter, "iterations was reached!"))
  if(!para$converged) 
    warning(paste("EM-algorithm did NOT converge!"))
  
  

  #######################
  # compute missing statistics for class "blim"
  # copied from blim.R

  para$npar <- 2 * nitems + nstates - 1

  ## Goodness of fit, df = number of patterns or persons
  fitted <- setNames(N*para$P.R, names(N.R))
  G2     <- 2*sum(N.R*log(N.R/fitted), na.rm=TRUE)
  df     <- min(2^nitems - 1, N) - para$npar        # number of patterns or persons
  # df     <- min(if(nitems <= zeropad) 2^nitems - 1 else npat, N) - npar
  gof    <- c(G2=G2, df=df, pval = ifelse(df > 0, 1 - pchisq(G2, df), NA))


  if (sum(para$P.R) < 1) para$P.R <- para$P.R/sum(para$P.R)     

  loglik <- sum(log(para$P.R) * N.R, na.rm=TRUE)


  ## Mean number of errors
  P.Kq <- numeric(nitems)
  for(j in seq_len(nitems))
    P.Kq[j] <- sum(para$P.K[which(K[,j] == 1)])
  nerror <- c("careless error" = sum(para$beta * P.Kq),
              "lucky guess"    = sum( para$eta * (1 - P.Kq)))


  ## Assigning state K given response R
  d.RK  <- apply(K, 1, function(k) colSums(xor(t(.R), k)))
  d.min <- apply(d.RK, 1, min, na.rm = TRUE)             # minimum discrepancy
  i.RK  <- (d.RK <= (d.min)) & !is.na(d.RK)

  ## Minimum discrepancy distribution
  disc.tab <- xtabs(N.R ~ d.min)
  disc     <- as.numeric(names(disc.tab)) %*% disc.tab / N

  para$discrepancy <-  disc
  para$disc.tab    <- disc.tab

  # add information to output
  para$method <- "blim"
  para$K <- K
  para$nitems <- nitems
  para$nstates <- nstates
  para$npatterns <- npat
  para$ntotal <- N
  para$N.R <- N.R
  para$N.RM <- N.R


  para$nerror <- nerror
  para$loglik <- loglik
  para$fitted.values <- fitted
  para$goodness.of.fit <- gof
  para$method <- "ML"
  para$na.use <- na.use
  

  # set names
  para$P.K  <- setNames(para$P.K,  as.pattern(K)) 
  para$beta <- setNames(para$beta, colnames(K))   
  para$eta  <- setNames(para$eta,  colnames(K))
  

  class(para) <- "blim"
  return(para)
}



