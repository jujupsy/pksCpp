
# custom testhat expect function

# mean absolute difference < tol

expect_diffLessTol <- function(object, expected, tolerance = 1e-14, 
                               label = NULL, expected.label = NULL) {
  
  act <- quasi_label(rlang::enquo(object), label, arg = "object")
  exp <- quasi_label(rlang::enquo(expected), expected.label, arg = "expected")
  
  abs_diff <- abs(unlist(act$val) - unlist(exp$val))
  
  
    
    expect(
      all(abs_diff < tolerance),
      sprintf("%s not equal to %s.\nMean absolute difference: %s \nMax. absolute difference: %s", 
              act$lab, exp$lab, mean(abs_diff), max(abs_diff))
    )
    invisible(act$val)
}



if(FALSE){
#library(pks)

  data(DoignonFalmagne7)
  data(endm)
  data(chess)
  data(Taagepera)

iterm <- 100000 # maxiter
tol   <- 1e-07
data7 <- list(
  df7    = blim(DoignonFalmagne7$K, DoignonFalmagne7$N.R, method = "ML"),
  endm      = blim(endm$K, endm$N.R, method = "ML", tol = tol,
                   maxiter = iterm),
  chess1    = blim(chess$dst1, chess$N.R, method = "ML", tol = tol,
                   maxiter = iterm),
  chess3    = blim(chess$dst3, chess$N.R, method = "ML", tol = tol,  
                   maxiter = iterm),
  chess4    = blim(chess$dst4, chess$N.R, method = "ML", tol = tol,
                   maxiter = iterm),
  density97 = blim(density97$K, density97$N.R, method = "ML", tol = tol, 
                   maxiter = iterm),
  matter97  = blim(matter97$K, matter97$N.R, method = "ML", tol = tol,
                   maxiter = iterm)

)

saveRDS(data7, "valid_data_tol7.rds")



## higher tolerance
iterm <- 1000000 # maxiter
tol   <- 1e-14

data14 <- list(
  df7    = blim(DoignonFalmagne7$K, DoignonFalmagne7$N.R, method = "ML"),
  endm      = blim(endm$K, endm$N.R, method = "ML", tol = tol,
                   maxiter = iterm),
  chess1    = blim(chess$dst1, chess$N.R, method = "ML", tol = tol,
                   maxiter = iterm),
  chess3    = blim(chess$dst3, chess$N.R, method = "ML", tol = tol,  
                   maxiter = iterm),
  chess4    = blim(chess$dst4, chess$N.R, method = "ML", tol = tol,
                   maxiter = iterm),
  density97 = blim(density97$K, density97$N.R, method = "ML", tol = tol, 
                   maxiter = iterm),
  matter97  = blim(matter97$K, matter97$N.R, method = "ML", tol = tol,
                   maxiter = iterm)

)

saveRDS(data14, "valid_data_tol14.rds")
}
