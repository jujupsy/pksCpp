
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::export]]
Rcpp::List emBLIMcpp(
  const Eigen::Map < Eigen::MatrixXd > R,
  const Eigen::Map < Eigen::MatrixXd > K,
  const Eigen::Map < Eigen::VectorXd > NR,
        Eigen::Map < Eigen::VectorXd > PK,
        Eigen::Map < Eigen::VectorXd > beta,
        Eigen::Map < Eigen::VectorXd > eta,
  const int maxiter = 100000,
  const double tol = 1e-07,
  const bool fdb = false

) {

  // declare objects
  Eigen::MatrixXd PRK(R.rows(), K.rows());
  Eigen::MatrixXd PR(R.rows(), 1);
  Eigen::MatrixXd PKR(K.rows(), R.rows());

  Eigen::MatrixXd diff(3, 1);

  Eigen::VectorXd PKold(PK.size());
  Eigen::VectorXd etaold(eta.size());
  Eigen::VectorXd betaold(beta.size());

  // generate additional matrices
  Eigen::MatrixXd notR = (1 - R.array()).matrix(); // 1 - R 
  Eigen::MatrixXd notK = (1 - K.array()).matrix(); // 1 - K
  
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

    // P(R|K)
    PRK = 
      (
       (   R * ((1 - beta.array()).log()).matrix().asDiagonal()) *    K.transpose() +
       (   R * (      eta.array( ).log()).matrix().asDiagonal()) * notK.transpose() +
       (notR * (     beta.array( ).log()).matrix().asDiagonal()) *    K.transpose() +
       (notR * ((1 -  eta.array()).log()).matrix().asDiagonal()) * notK.transpose()
      ).array().exp().matrix();

    // P(R) = sum_K(P(R, K)) = sum_K(P(R | K) * pi_K)
    PR = PRK * PK;
    
    // P(K | R)
    PKR = (PK * PR.cwiseInverse().transpose()).cwiseProduct(PRK.transpose());

    // M - Step
    //
    // save old estimates
    PKold = PK;
    etaold = eta;
    betaold = beta;

    // PK
    PK = (PKR * NR);
    PK /= NR.sum();

    // beta, eta
    beta = 
      (
        (((PKR.transpose() * K).cwiseProduct(notR)).transpose()) * NR
      ).cwiseQuotient(
        ((PKR.transpose() * K).transpose()) * NR);

    eta = 
      (
        (((PKR.transpose() * notK).cwiseProduct(R)).transpose()) * NR
      ).cwiseQuotient(
        ((PKR.transpose() * notK).transpose()) * NR);
     
    diff(0, 0) = ((PKold - PK).cwiseAbs().maxCoeff()) < tol;
    diff(1, 0) = ((etaold - eta).cwiseAbs().maxCoeff()) < tol;
    diff(2, 0) = ((betaold - beta).cwiseAbs().maxCoeff()) < tol;

    iter++;

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
  
  } while (diff.sum() < 3 && iter < maxiter);

  
  if (diff.sum() >= 3) {
    converged = true;
  }

  if (fdb) {
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "     **** DONE ****" << std::endl;

    if (iter >= maxiter) {
      Rcpp::Rcout << "Maximum # of " << maxiter << " iterations reached." << std::endl;
    }

    if (converged == false) {
      Rcpp::Rcout << "EM-Algorithm did NOT converged!" << std::endl;
    }

  }

 // compute PKR with final estimates
  beta = (beta.array() > (1 - eps)).select(upperLimit, 
                     (beta.array() < eps).select(lowerLimit, beta));
  eta  = ( eta.array() > (1 - eps)).select(upperLimit, 
                     ( eta.array() < eps).select(lowerLimit,  eta));

  PRK = 
    (
       (   R * ((1 - beta.array()).log()).matrix().asDiagonal()) *    K.transpose() +
       (   R * (      eta.array( ).log()).matrix().asDiagonal()) * notK.transpose() +
       (notR * (     beta.array( ).log()).matrix().asDiagonal()) *    K.transpose() +
       (notR * ((1 -  eta.array()).log()).matrix().asDiagonal()) * notK.transpose()
    ).array().exp().matrix();

  PR = PRK * PK;

  PKR = (PK * PR.cwiseInverse().transpose()).cwiseProduct(PRK.transpose());
  
  return Rcpp::List::create(Rcpp::Named("P.K") = PK,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("iterations") = iter,
    Rcpp::Named("maxiter") = iter >= maxiter,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("PK.R") = PKR,
    Rcpp::Named("P.R") = PR
    
  );

}
