


# test with example data [validated using blim() from pks package v 0.4-1]
# see gen_valid_data.R

val7 <- readRDS(test_path(".", "valid_data_tol7.rds"))

comp <- c("beta", "eta", "P.K", "PK.R", "P.R", "converged", "iterations")

## DoignonFalmagne7
data(DoignonFalmagne7)
context("blimCpp: DoignonFalmagne7")
#######################################################################
test_that("DoignonFalmagne7", {                                
  expect_diffLessTol(
    blimCpp(K = DoignonFalmagne7$K, N.R = DoignonFalmagne7$N.R)[c("beta", "eta", "P.K")], 
    val7$df7[c("beta", "eta", "P.K")]
  )   
})                                                             
#######################################################################


tol <- 2e-6

## chess1
data(chess)
context("blimCpp: chess dst1")
#######################################################################
test_that("Chess.dst1", {                                
  expect_diffLessTol(
    blimCpp(K = chess$dst1, N.R = chess$N.R)[c("beta", "eta", "P.K")], 
    val7$chess1[c("beta", "eta", "P.K")], tolerance = tol
  )   
})                                                             
#######################################################################


## chess3
data(chess)
context("blimCpp: chess dst3")
#######################################################################
test_that("Chess.dst3", {                                
  expect_diffLessTol(
    blimCpp(K = chess$dst3, N.R = chess$N.R)[c("beta", "eta", "P.K")], 
    val7$chess3[c("beta", "eta", "P.K")], tolerance = tol
  )     
})                                                             
#######################################################################


## chess3
data(chess)
context("blimCpp: chess dst4")
#######################################################################
test_that("Chess.dst4", {                                
  expect_diffLessTol(
    blimCpp(K = chess$dst4, N.R = chess$N.R)[c("beta", "eta", "P.K")], 
    val7$chess4[c("beta", "eta", "P.K")], tolerance = tol
  )   
})                                                             
#######################################################################




## blim <=> IMBLIM <=> MissBLIM if M = 0
comp <- c("beta", "eta", "P.K", "PK.R", "P.R", "converged", "iterations")
data(chess)
context("blim <=> IMBLIM <=> MissBLIM if M = 0")
#######################################################################
test_that("Chess.dst4", {  
  # blim <=> blim omit if M=0
  expect_diffLessTol(
    blimCpp(K = chess$dst3, N.R = chess$N.R)[comp], 
    blimCpp(K = chess$dst3, N.R = chess$N.R, na.use = "omit")[comp]
  )   

  # blim <=> blim missing-as-wrong if M=0
  expect_diffLessTol(
    blimCpp(K = chess$dst3, N.R = chess$N.R)[comp], 
    blimCpp(K = chess$dst3, N.R = chess$N.R, na.use = "missing-as-wrong")[comp]
  )   

  # blim <=> IMBLIM if M=0
  expect_diffLessTol(
    blimCpp(K = chess$dst3, N.R = chess$N.R)[comp], 
    blimCpp(K = chess$dst3, N.R = chess$N.R, na.use = "IMBLIM")[comp]
  )   

})                                                             
#######################################################################















seed <- 42



# simulate random probabilistic knowledge structure:
# 25 Items Q and 500 random knwoledge states
# beta, eta ~ unif(0, 0.1);  P(K) ~ unif(0, 1) -> normalized
true_pks <- genPKS(kQ = 10, kK = 100, seed = seed, coef_gen = runif, min = 0, max = 0.1)

res_beta <- data.frame()
res_eta  <- data.frame()
res_pk   <- data.frame()

# sim data without missings BLIM
sim1 <- simData(N = 1000, true_model = true_pks, seed = seed)

# estimate models and compare with true parameters
m1 <- blimCpp(K = sim1$tModel$K, N.R = as.pattern(sim1$R, freq = TRUE))

eta_res  <- m1$eta  - sim1$tModel$eta
beta_res <- m1$beta - sim1$tModel$beta
pk_res   <- m1$P.K  - sim1$tModel$P.K



# sim data with missings

### Missing Completely at Random (MCAR) ###
# P(Mq = 1 | q in K) = P(Mq = 1 | q not in K) = Pmcar 

# add mu0, mu1 to true model
true_pks$mu0 <- rep(0.05, true_pks$nitems)
true_pks$mu1 <- rep(0.05, true_pks$nitems)

sim1_na <- simData(N = 100, true_model = true_pks, seed = seed)

# number of missings distribution
table(rowSums(is.na(sim1_na$R)))



