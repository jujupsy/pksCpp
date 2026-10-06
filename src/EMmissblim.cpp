
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>


// [[Rcpp::export]]
Rcpp::List emMissBLIMcpp(
        const Eigen::Map<Eigen::MatrixXd> Rr,
        const Eigen::Map<Eigen::MatrixXd> Wr,
        const Eigen::Map<Eigen::MatrixXd> Mr,
        const Eigen::Map<Eigen::MatrixXd> Kr, 
        const Eigen::Map<Eigen::VectorXd> NRr,
        const Eigen::Map<Eigen::VectorXd> PKr,
        const Eigen::Map<Eigen::VectorXd> etar, 
        const Eigen::Map<Eigen::VectorXd> betar,
        const Eigen::Map<Eigen::VectorXd> mu0r, 
        const Eigen::Map<Eigen::VectorXd> mu1r, 
        const int maxiter = 1000, 
        const double tol  = 1e-07,
        const bool fdb    = true
        
){

//copy mapped input

Eigen::MatrixXd R    = Rr;
Eigen::MatrixXd W    = Wr;
Eigen::MatrixXd M    = Mr;
Eigen::MatrixXd K    = Kr; 
Eigen::VectorXd NR   = NRr;
Eigen::VectorXd PK   = PKr;
Eigen::VectorXd eta  = etar; 
Eigen::VectorXd beta = betar;
Eigen::VectorXd mu0  = mu0r; 
Eigen::VectorXd mu1  = mu1r;     

// declare objects including type      
Eigen::MatrixXd PMK;
Eigen::MatrixXd PRMSK;
Eigen::MatrixXd PRMK;
Eigen::MatrixXd PRM;
Eigen::MatrixXd PKRM;

Eigen::MatrixXd diff(5, 1);  

Eigen::VectorXd PKold;
Eigen::VectorXd etaold; 
Eigen::VectorXd betaold;
Eigen::VectorXd mu0old; 
Eigen::VectorXd mu1old;

int iter = 0;
double eps = 1e-9;
bool converged = false;

do {
   
    // E-Step
    //prevent log zero // maybe .array() <= eps is faster....
    for(int i = 0; i < mu0.rows(); i++){
        if(mu0(i) < eps){
            mu0(i) += eps;
        }else if (mu0(i) > (1 - eps)){
            mu0(i) -= eps; 
        }

        if(mu1(i) < eps){
            mu1(i) += eps;
        }else if (mu1(i) > (1 - eps)){
            mu1(i) -= eps; 
        }

        if(eta(i) < eps){
            eta(i) += eps;
        }else if (eta(i) > (1 - eps)){
            eta(i) -= eps; 
        }

        if(beta(i) < eps){
            beta(i) += eps;
        }else if (beta(i) > (1 - eps)){
            beta(i) -= eps; 
        }
    } 


    PMK = (M.cwiseProduct(
              (mu1.array().log().matrix()
               ).replicate(1, M.rows()).transpose()
              )* K.transpose()
          +
          M.cwiseProduct(
             (mu0.array().log().matrix()
              ).replicate(1, M.rows()).transpose()
             ) * (((1 - K.array()).matrix()).transpose())
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    ((1 - mu1.array()).log().matrix()
                                     ).replicate(1, M.rows()).transpose()
                                    ) * K.transpose()
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    ((1 - mu0.array()).log().matrix()
                                     ).replicate(1, M.rows()).transpose()
                                    ) * (
                                        ((1 - K.array()).matrix()
                                            ).transpose()
                                        )
          ).array().exp().matrix();

    PRMSK = ((1 - M.array()).matrix().cwiseProduct(
                                       R.cwiseProduct(
                                        ((1 - beta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                       ) * K.transpose()
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    R.cwiseProduct(
                                     (eta.array().log().matrix()
                                      ).replicate(1, M.rows()).transpose()
                                     )
                                    ) * (
                                         ((1 - K.array()).matrix()).transpose()
                                         )
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    (1 - R.array()).matrix().cwiseProduct(
                                        (beta.array().log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * K.transpose()
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    (1 - R.array()).matrix().cwiseProduct(
                                        ((1 - eta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * (((1 - K.array()).matrix()).transpose())
          ).array().exp().matrix();


    PRMK = PRMSK.cwiseProduct(PMK);
    PRM = PRMK * PK;


    PKRM  = (PK * PRM.cwiseInverse().transpose()).cwiseProduct(PRMK.transpose());

    // M - Step

    // save old estimates
    PKold = PK;
    etaold = eta; 
    betaold = beta;
    mu0old = mu0; 
    mu1old = mu1;

    // PK
    PK  = (PKRM * NR); 
    PK /= NR.sum();

    // beta, eta
    beta  = (
             (((PKRM.transpose() * K).cwiseProduct(W)).transpose()) * NR
             ).cwiseQuotient(
                (((PKRM.transpose() * K).cwiseProduct((1 - M.array()).matrix())).transpose()) * NR);

    eta   = (
             (((PKRM.transpose() * ((1 - K.array())).matrix()).cwiseProduct(R)).transpose()) * NR
             ).cwiseQuotient(
                (((PKRM.transpose() * ((1 - K.array()).matrix())).cwiseProduct((1 - M.array()).matrix())).transpose()) * NR);


    //mu0, mu1
    mu1   = (
             (((PKRM.transpose() * K).cwiseProduct(M)).transpose()) * NR
             ).cwiseQuotient(
                ((PKRM.transpose() * K).transpose()) * NR);

    mu0   = (
             (((PKRM.transpose() * ((1 - K.array())).matrix()).cwiseProduct(M)).transpose()) * NR
             ).cwiseQuotient(
                ((PKRM.transpose() * ((1 - K.array()).matrix())).transpose()) * NR);


    diff(0, 0) = (( PKold  -   PK).cwiseAbs().maxCoeff()) < tol;
    diff(1, 0) = (( etaold -  eta).cwiseAbs().maxCoeff()) < tol;
    diff(2, 0) = ((betaold - beta).cwiseAbs().maxCoeff()) < tol;
    diff(3, 0) = (( mu0old -  mu0).cwiseAbs().maxCoeff()) < tol; 
    diff(4, 0) = (( mu1old -  mu1).cwiseAbs().maxCoeff()) < tol;

    iter ++;
    
    if(fdb){
        if(iter % 50 == 0){
            Rcpp::Rcout << ".  " << "Iteration #: " << iter << std::endl;
            Rcpp::checkUserInterrupt();
        } else {
            Rcpp::Rcout << ".";
        }
    } else if (iter % 10 == 0) {
      Rcpp::checkUserInterrupt();
    }
}
while (diff.sum() < 5 && iter < maxiter);

if(diff.sum() >= 5){
    converged = true;
}

if(fdb){
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "     **** EM-DONE ****" << std::endl;

}


 // compute PKR with final estimates
    for(int i = 0; i < mu0.rows(); i++){
        if(mu0(i) < eps){
            mu0(i) += eps;
        }else if (mu0(i) > (1 - eps)){
            mu0(i) -= eps; 
        }

        if(mu1(i) < eps){
            mu1(i) += eps;
        }else if (mu1(i) > (1 - eps)){
            mu1(i) -= eps; 
        }

        if(eta(i) < eps){
            eta(i) += eps;
        }else if (eta(i) > (1 - eps)){
            eta(i) -= eps; 
        }

        if(beta(i) < eps){
            beta(i) += eps;
        }else if (beta(i) > (1 - eps)){
            beta(i) -= eps; 
        }
    } 


    PMK = (M.cwiseProduct(
              (mu1.array().log().matrix()
               ).replicate(1, M.rows()).transpose()
              )* K.transpose()
          +
          M.cwiseProduct(
             (mu0.array().log().matrix()
              ).replicate(1, M.rows()).transpose()
             ) * (((1 - K.array()).matrix()).transpose())
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    ((1 - mu1.array()).log().matrix()
                                     ).replicate(1, M.rows()).transpose()
                                    ) * K.transpose()
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    ((1 - mu0.array()).log().matrix()
                                     ).replicate(1, M.rows()).transpose()
                                    ) * (
                                        ((1 - K.array()).matrix()
                                            ).transpose()
                                        )
          ).array().exp().matrix();

    PRMSK = ((1 - M.array()).matrix().cwiseProduct(
                                       R.cwiseProduct(
                                        ((1 - beta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                       ) * K.transpose()
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    R.cwiseProduct(
                                     (eta.array().log().matrix()
                                      ).replicate(1, M.rows()).transpose()
                                     )
                                    ) * (
                                         ((1 - K.array()).matrix()).transpose()
                                         )
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    (1 - R.array()).matrix().cwiseProduct(
                                        (beta.array().log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * K.transpose()
          +
          (1 - M.array()).matrix().cwiseProduct(
                                    (1 - R.array()).matrix().cwiseProduct(
                                        ((1 - eta.array()).log().matrix()
                                         ).replicate(1, M.rows()).transpose()
                                        )
                                    ) * (((1 - K.array()).matrix()).transpose())
          ).array().exp().matrix();


    PRMK = PRMSK.cwiseProduct(PMK);
    PRM = PRMK * PK;


    PKRM  = (PK * PRM.cwiseInverse().transpose()).cwiseProduct(PRMK.transpose());


//Rcpp::Rcout << typeid(PK).name() << std::endl;
//Rcpp::Rcout << typeid(beta).name() << std::endl;
//Rcpp::Rcout << typeid(eta).name() << std::endl;
//Rcpp::Rcout << typeid(mu0).name() << std::endl;
//Rcpp::Rcout << typeid(mu1).name() << std::endl;
//Rcpp::Rcout << typeid(iter).name() << std::endl;
//Rcpp::Rcout << typeid(diff).name() << std::endl;


return Rcpp::List::create(Rcpp::Named("P.K") = PK,
                          Rcpp::Named("beta") = beta,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("mu0") = mu0,
                          Rcpp::Named("mu1") = mu1,
                          Rcpp::Named("iterations") = iter,
                          Rcpp::Named("maxiter") = iter >= maxiter,
                          Rcpp::Named("converged") = converged,
                          Rcpp::Named("PK.R") = PKRM,
                          Rcpp::Named("P.R") = PRM
                          );

}
