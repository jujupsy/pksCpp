#' Basic Local Independence Models (BLIMs) using C++
#'
#' @description
#' Fits a basic local independence model (BLIM) for probabilistic 
#' knowledge structures by maximum likelihood estimation.
#'
#' @param K a state-by-problem indicator matrix representing the knowledge
#'          structure.  An element is one if the problem is contained in the 
#'          state, and else zero.
#' @param N.R a (named) vector of absolute frequencies of response patterns.
#' @param tol tolerance, stopping criterion for iteration.
#' @param maxiter the maximum number of iterations.
#' @param fdb feedback during the optimization.
#'
#' @details See Doignon and Falmagne (1999) for details on the basic local independence
#'  model (BLIM) for probabilistic knowledge structures.
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
blimCpp <- function(K, N.R, tol = 1e-07 , maxiter = 10000, 
                            fdb = FALSE){

  model <- "blim"
  K     <- as.matrix(K)
  R     <- as.binmat(N.R, uniq = TRUE)

  # input checks
  if(ncol(K) != ncol(R)) stop("Matrix K and R must have the same number of columns")

  N       <- sum(N.R)  # sample size
  nitems  <- ncol(K)   # number of items q in Q
  nstates <- nrow(K)   # number of states k in K
  npat    <- nrow(R)   # number of unique response pattern

  ## set initial values
  P.K  <- setNames(rep(1/nstates, nstates),  rownames(K)) # probability vector P(K)
  beta <- setNames(rep(0.1, nitems), colnames(K))         # initital parameter estimations 
  eta  <- setNames(rep(0.1, nitems), colnames(K))

  R <- (R == 1) # R: right pattern

  iter <- 1
  maxdiff <- 2 * tol

  
  # convert to double for c++ function
  mode(R) <- mode(K) <- mode(N.R) <- mode(P.K) <- mode(beta) <- 
             mode(eta) <- "double"

  # EM Loop 
  para <- emBLIMcpp(R, K, N.R, P.K, beta, eta, maxiter, tol, fdb)
  
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
  # if no zero padding: normalize
  loglik <- sum(log(para$P.R) * N.R, na.rm=TRUE)


  ## Mean number of errors
  P.Kq <- numeric(nitems)
  for(j in seq_len(nitems))
    P.Kq[j] <- sum(para$P.K[which(K[,j] == 1)])
  nerror <- c("careless error" = sum(para$beta * P.Kq),
              "lucky guess"    = sum( para$eta * (1 - P.Kq)))


  ## Assigning state K given response R
  d.RK  <- apply(K, 1, function(k) colSums(xor(t(R), k)))
  d.min <- apply(d.RK, 1, min, na.rm = TRUE)             # minimum discrepancy
  i.RK  <- (d.RK <= (d.min)) & !is.na(d.RK)

  ## Minimum discrepancy distribution
  disc.tab <- xtabs(N.R ~ d.min)
  disc     <- as.numeric(names(disc.tab)) %*% disc.tab / N

  para$discrepancy <-  disc
  para$disc.tab    <- disc.tab

  # add information to output
  para$method <- model
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

  # set names
  para$P.K  <- setNames(para$P.K,  as.pattern(K)) 
  para$beta <- setNames(para$beta, colnames(K))   
  para$eta  <- setNames(para$eta,  colnames(K))

  class(para) <- "blim"
  return(para)
}



