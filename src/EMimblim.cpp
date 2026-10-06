
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


// [[Rcpp::export]]
Rcpp::List emIMBLIMcpp(
        const Eigen::Map<Eigen::MatrixXd> R,
        const Eigen::Map<Eigen::MatrixXd> W,
        const Eigen::Map<Eigen::MatrixXd> M,
        const Eigen::Map<Eigen::MatrixXd> K, 
        const Eigen::Map<Eigen::VectorXd> NR,
        Eigen::Map<Eigen::VectorXd> PK,
        Eigen::Map<Eigen::VectorXd> beta, 
        Eigen::Map<Eigen::VectorXd> eta,
        const Eigen::Map<Eigen::MatrixXd> PM, 
        const int maxiter = 1000, 
        const double tol  = 1e-07,
        const bool fdb    = true
        
){

// declare objects     
Eigen::MatrixXd PRMSK;
Eigen::MatrixXd PRMK;
Eigen::MatrixXd PRM;
Eigen::MatrixXd PKRM;

Eigen::MatrixXd diff(3, 1);  

Eigen::VectorXd PKold(PK.size());
Eigen::VectorXd etaold(eta.size()); 
Eigen::VectorXd betaold(beta.size());

  // generate additional matrices
  Eigen::MatrixXd notR = (1 - R.array()).matrix(); // 1 - R 
  Eigen::MatrixXd notK = (1 - K.array()).matrix(); // 1 - K
  Eigen::MatrixXd notM = (1 - M.array()).matrix(); // 1 - M
  
// log(0) and log(1 - 1) prevention limits
double eps = 1e-9;
Eigen::VectorXd lowerLimit(beta.size());
Eigen::VectorXd upperLimit(beta.size());
lowerLimit.fill(eps);
upperLimit.fill(1 - eps);

int iter = 0;
bool converged = false;

do {
   
    // E-Step
    //
    //prevent log zero
    beta = (beta.array() > (1 - eps)).select(upperLimit, 
            (beta.array() < eps).select(lowerLimit, beta));
    eta  = ( eta.array() > (1 - eps)).select(upperLimit, 
                     ( eta.array() < eps).select(lowerLimit,  eta));
 

    PRMSK = (notM.cwiseProduct(
                                       R.cwiseProduct(
                                        ((1 - beta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                       ) * K.transpose()
          +
          notM.cwiseProduct(
                                    R.cwiseProduct(
                                     (eta.array().log().matrix()
                                      ).replicate(1, M.rows()).transpose()
                                     )
                                    ) * notK.transpose()
          +
          notM.cwiseProduct(
                                    notR.cwiseProduct(
                                        (beta.array().log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * K.transpose()
          +
          notM.cwiseProduct(
                                    notR.cwiseProduct(
                                        ((1 - eta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * notK.transpose()
          ).array().exp().matrix();


    PRMK = PRMSK.cwiseProduct(PM);
    PRM = PRMK * PK;


    PKRM  = (PK * PRM.cwiseInverse().transpose()).cwiseProduct(PRMK.transpose());

    // M - Step

    // save old estimates
    PKold = PK;
    etaold = eta; 
    betaold = beta;
    
    // PK
    PK  = (PKRM * NR); 
    PK /= NR.sum();

    // beta, eta
    beta  = (
             (((PKRM.transpose() * K).cwiseProduct(W)).transpose()) * NR
             ).cwiseQuotient(
                (((PKRM.transpose() * K).cwiseProduct(notM)).transpose()) * NR);

    eta   = (
             (((PKRM.transpose() * notK).cwiseProduct(R)).transpose()) * NR
             ).cwiseQuotient(
                (((PKRM.transpose() * notK).cwiseProduct(notM)).transpose()) * NR);




    diff(0, 0) = (( PKold  -   PK).cwiseAbs().maxCoeff()) < tol;
    diff(1, 0) = (( etaold -  eta).cwiseAbs().maxCoeff()) < tol;
    diff(2, 0) = ((betaold - beta).cwiseAbs().maxCoeff()) < tol;

    iter ++;
    
    if (fdb) {
      if (iter % 50 == 0) {
        Rcpp::Rcout << ".  " << "Iteration #: " << iter << std::endl;
        Rcpp::checkUserInterrupt();
      } else {
        Rcpp::Rcout << ".";
      }
    } else if (iter % 10 == 0) {
      Rcpp::checkUserInterrupt();
    }

}
while (diff.sum() < 3 && iter < maxiter);

if(diff.sum() >= 3){
    converged = true;
}

if(fdb){
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "     **** DONE ****" << std::endl;

    if(iter >= maxiter){
        Rcpp::Rcout << "Maximum # of " << maxiter << " Iterations reached." << std::endl;
    }

    if(converged == false){
        Rcpp::Rcout << "EM-Algorithm did NOT converged!" << std::endl;
    }

}




 // compute PKR with final estimates
  beta = (beta.array() > (1 - eps)).select(upperLimit, 
                     (beta.array() < eps).select(lowerLimit, beta));
  eta  = ( eta.array() > (1 - eps)).select(upperLimit, 
                     ( eta.array() < eps).select(lowerLimit,  eta));

 
    PRMSK = (notM.cwiseProduct(
                                       R.cwiseProduct(
                                        ((1 - beta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                       ) * K.transpose()
          +
          notM.cwiseProduct(
                                    R.cwiseProduct(
                                     (eta.array().log().matrix()
                                      ).replicate(1, M.rows()).transpose()
                                     )
                                    ) * notK.transpose()
          +
          notM.cwiseProduct(
                                    notR.cwiseProduct(
                                        (beta.array().log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * K.transpose()
          +
          notM.cwiseProduct(
                                    notR.cwiseProduct(
                                        ((1 - eta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * notK.transpose()
          ).array().exp().matrix();


    PRMK = PRMSK.cwiseProduct(PM);
    PRM = PRMK * PK;


    PKRM  = (PK * PRM.cwiseInverse().transpose()).cwiseProduct(PRMK.transpose());


return Rcpp::List::create(Rcpp::Named("P.K") = PK,
                          Rcpp::Named("beta") = beta,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("iterations") = iter,
                          Rcpp::Named("maxiter") = iter >= maxiter,
                          Rcpp::Named("converged") = converged,
                          Rcpp::Named("PK.R") = PKRM,
                          Rcpp::Named("P.R") = PRM
                          );

}





