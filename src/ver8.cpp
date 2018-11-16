#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// #include "mvtnorm.h"
// #include <stdlib.h>
// #include <stdio.h>
using namespace Rcpp;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);
// const static char errorMessage0[] = "Normal Completion";                           // inform = 0
// const static char errorMessage1[] = "Completion with error > abseps";              // inform = 1
// const static char errorMessage2[] = "N greater 1000 or N < 1";                     // inform = 2
// const static char errorMessage3[] ="Covariance matrix not positive semidefinite";  // inform = 3
// const static char* errorMessage[4] = {errorMessage0, errorMessage1, errorMessage2, errorMessage3};
// const static int INFIN_BOUND_NORMAL = 2;        // (..., ...)
// const static int INFIN_BOUND_UPPER = 1;         // (..., inf)
// const static int INFIN_BOUND_LOWER = 0;         // (-inf, ..)
// const static int INFIN_BOUND_LOWER_UPPER = -1;  // (-inf, inf)


// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// double pmvnorm(int n,
//                int nu,
//                double *lower,
//                double *upper,
//                int *infin,
//                double correl,
//                double *delta, // non-central parameter
//                int maxpts,    // param
//                double abseps, // param
//                double releps, // param
//                double error,  // estimated abs. error. with 99% confidence interval
//                double value,     // results store here.
//                int inform)    // inform message goes here
// {
//   mvtdst_ (n, nu,
//            lower, upper, infin, correl, delta,
//            maxpts, abseps, releps,
//            error, value, inform);
//   printf ("error = %g, value = %g, inform = %d\n", error, value, inform);
//
//   switch (inform) {
//   case 0:
//     return value;
//   case 1:
//   case 2:
//   case 3:
//     return -1.0;
//   };
//
//   return value;
// };
//
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// double pmvnorm_P(int n,
//                  double bound,
//                  double correlationMatrix, // (2,1), (3,1), (3,2) .....
//                  double error)
// {
//   int nu_ = 0;
//   int maxpts_ = 25000;     // default in mvtnorm: 25000
//   double abseps_ = 1e-6;   // default in mvtnorm: 0.001, we make it more stringent
//   double releps_ = 0;      // default in mvtnorm: 0
//
//   // arma::vec lower = arma::zeros<arma::vec>(n);
//   // arma::vec infin = arma::zeros<arma::vec>(n);
//   // arma::vec delta = arma::zeros<arma::vec>(n);
//   double lower = new double[n];
//   int infin = new int[n];
//   double delta = new double[n];
//
//   int i = 0;
//   for (i = 0; i < n; ++i) {
//     infin[i] = 0; // (-inf, bound]
//     lower[i] = 0.0;
//     delta[i] = 0.0;
//     // infin(i) = 0; // (-inf, bound]
//     // lower(i) = 0.0;
//     // delta(i) = 0.0;
//   }
//
//   // return values
//   double value_ = 0.0;
//   int inform_ = 0.0;
//
//   double ret = pmvnorm(&n, &nu_, lower, bound, infin, correlationMatrix, delta, &maxpts_, &abseps_, &releps_, error, &value_, &inform_);
//   delete[] (lower);
//   delete[] (infin);
//   delete[] (delta);
//
//   return ret;
// }
//
// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// double pmvnorm_Q(int n,
//                  double bound,
//                  double correlationMatrix, // (2,1), (3,1), (3,2) .....
//                  double error)
// {
//   int nu_ = 0;
//   int maxpts_ = 25000;
//   double abseps_ = 1e-6;
//   double releps_ = 0;
//
//
//   double *upper = new double[n];
//   int *infin = new int[n];
//   double *delta = new double[n];
//
//   int i = 0;
//   for (i = 0; i < n; ++i) {
//     infin[i] = 1; // (-inf, bound]
//     upper[i] = 0.0;
//     delta[i] = 0.0;
//   }
//
//   // return values
//   double value_ = 0.0;
//   int inform_ = 0.0;
//
//   double ret = pmvnorm(&n, &nu_, bound, upper, infin, correlationMatrix, delta, &maxpts_, &abseps_, &releps_, error, &value_, &inform_);
//   delete[] (upper);
//   delete[] (infin);
//   delete[] (delta);
//
//   return ret;
// }


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::rowvec colmeanNA(arma::mat Y){
  int T = Y.n_cols;
  arma::rowvec colmean = arma::zeros<arma::rowvec>(T);
  for (int i=0; i < T; ++i){
    arma::vec ycol = Y.col(i);
    arma::uvec finiteind = find_finite(ycol);
    arma::vec yy = ycol(finiteind);
    colmean(i) = mean(yy);
  }
  return colmean;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
int sample_index(const int size,
                 const NumericVector prob = NumericVector::create()){
  //sample one number from 1:size
  arma::vec sequence = arma::linspace<arma::vec>(1, size, size);
  arma::vec out = Rcpp::RcppArmadillo::sample(sequence, size, false, prob);
  return out(0);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mvrnormArma(const int n,
                      const arma::vec mu,
                      const arma::mat Sigma) {
  //returns random multivariate normal vectors with mean mu and covariance Sigma
  //input : integer n for the number of vectors you'd like to draw
  //      : vector mu for the mean
  //      : matrix Sigma for the covariance - needs to be psd
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dmvnrm_arma(const arma::rowvec x,
                   const arma::rowvec mean,
                   const arma::mat sigma,
                   const bool logd = false) {
  //returns the density of a multivariate normal vector
  //input : a rowvector x whose density you'd like to know
  //      : a rowvector mean for the mean of mvn
  //      : a matrix Sigma for covariance, needs to be psd
  //      : a boolean logd, true if you like the log density
  int xdim = x.n_cols;
  if(xdim==0){return 0;}
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  arma::vec z = rooti*arma::trans(x-mean);
  out = constants - 0.5*arma::sum(z%z)+rootisum;
  if (logd == false){
    out = exp(out);
  }
  return(out);
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat em_with_zero_mean_c(arma::mat y,
                              const int maxit){
  int orig_p = y.n_cols;
  arma::vec vars = arma::zeros<arma::vec>(orig_p);
  for (int i=0; i < orig_p; ++i){
    arma::vec ycol = y.col(i);
    arma::uvec finiteind = find_finite(ycol);
    arma::vec yy = ycol(finiteind);
    vars(i) = sum((yy-mean(yy))%(yy-mean(yy)));
  }
  arma::uvec valid_ind = find(vars>1e-6);
  y = y.cols(valid_ind);
  int p = y.n_cols;
  int n = y.n_rows;
  arma::rowvec mu = arma::zeros<arma::rowvec>(p);
  arma::mat y_imputed = y;
  for (int j = 0; j < p; ++j){
    arma::uvec colind = arma::zeros<arma::uvec>(1);
    colind(0) = j;
    arma::uvec nawhere = find_nonfinite(y_imputed.col(j));
    arma::uvec nonnawhere = find_finite(y_imputed.col(j));
    arma::vec tempcolmean = mean(y_imputed(nonnawhere, colind), 0);
    y_imputed(nawhere, colind).fill(tempcolmean(0));
  }
  arma::mat oldSigma = y_imputed.t() * y_imputed / n;
  arma::mat Sigma = oldSigma;
  double diff = 1;
  int it = 1;
  while (diff>0.001 && it < maxit){
    arma::mat bias = arma::zeros<arma::mat>(p,p);
    for (int i=0; i<n; ++i){
      arma::rowvec tempdat = y.row(i);
      arma::uvec ind = find_finite(tempdat);
      arma::uvec nind = find_nonfinite(tempdat);
      if (0 < ind.size() && ind.size() < p){
        //MAKE THIS PART FASTER
        bias(nind, nind) += Sigma(nind, nind) - Sigma(nind, ind) * (Sigma(ind, ind).i()) * Sigma(ind, nind);
        arma::uvec rowind = arma::zeros<arma::uvec>(1);
        rowind(0) = i;
        arma::mat yvec = y(rowind, ind);
        //MAKE THIS PART FASTER
        y_imputed(rowind, nind) = (Sigma(nind, ind)*(Sigma(ind, ind).i())*y(rowind, ind).t()).t();
      }
    }
    Sigma = (y_imputed.t() * y_imputed + bias)/n;
    arma::mat diffmat = (Sigma-oldSigma);
    arma::mat diffsq = diffmat%diffmat;
    diff = accu(diffsq);
    oldSigma = Sigma;
    it = it + 1;
  }
  arma::mat finalSigma = arma::zeros<arma::mat>(orig_p, orig_p);
  finalSigma.submat(valid_ind, valid_ind.t()) = Sigma;
  /*
  return Rcpp::List::create(
    Rcpp::Named("Sigma") = finalSigma,
  Rcpp::Named("iteration") = it
  )*/
  return finalSigma;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::cube rinvwish_c(const int n,
                      const int v,
                      const arma::mat S){
  //draw a matrix from inverse wishart distribution with parameters S and v
  RNGScope scope;
  int p = S.n_rows;
  arma::mat L = chol(inv_sympd(S), "lower");
  arma::cube sims(p, p, n, arma::fill::zeros);
  for(int j = 0; j < n; j++){
    arma::mat A(p,p, arma::fill::zeros);
    for(int i = 0; i < p; i++){
      int df = v - (i + 1) + 1; //zero-indexing
      A(i,i) = sqrt(R::rchisq(df));
    }
    for(int row = 1; row < p; row++){
      for(int col = 0; col < row; col++){
        A(row, col) = R::rnorm(0,1);
      }
    }
    arma::mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
    sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec get_target_c(const arma::vec X,
                       const arma::mat Y,
                       const double sigmabeta,
                       const arma::mat Sigma,
                       const arma::vec gam,
                       const arma::vec beta,
                       const int aa,
                       const int bb){
  //get the target likelihood circumventing the missing value issue
  int T = Y.n_cols;
  int n = Y.n_rows;
  double L = 0;
  double B = 0;
  double G = 0;
  for (int i=0; i < n; ++i){
    arma::uvec naind = find_finite(Y.row(i).t());
    if(naind.size()>0){
      arma::uvec rowind = arma::zeros<arma::uvec>(1);
      rowind(0) = i;
      arma::rowvec Ytemp = Y(rowind, naind.t());
      L = L + dmvnrm_arma(Ytemp,
                          X(i)*beta(naind).t(),
                          Sigma(naind,naind.t()),
                          true);
    }
  }
  arma::uvec ind = find(gam==1);
  int s = ind.size();
  if(s>0){
    arma::vec ds = Sigma.diag();
    for (int j=0; j<s; ++j){
      int newind = ind(j);
      B = B + R::dnorm(beta(newind), 0, sqrt(sigmabeta), true);
    }
  }


  // double aa = 1;
  // double bb = 10;
  G = log(std::tgamma(s+aa) * std::tgamma(T+bb-s) * std::tgamma(aa+bb));
  G = G - log(std::tgamma(T+aa+bb) * std::tgamma(aa) * std::tgamma(bb));
  //
  // G = log(std::tgamma(s+1)*std::tgamma(T-s+1)/std::tgamma(T+2));
  // G = 0;
  arma::vec out = arma::zeros<arma::vec>(3);
  out(0) = L;
  out(1) = B;
  out(2) = G;
  return out;
}











// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat update_Sigma_c(const int n,
                         const int nu,
                         const arma::vec X,
                         const arma::vec beta,
                         const arma::mat Phi,
                         const arma::mat Y){
  int T = Y.n_cols;
  arma::mat X2 = arma::zeros<arma::mat>(n, 1);
  X2.col(0) = X;
  arma::mat beta2 = arma::zeros<arma::mat>(T,1);
  beta2.col(0) = beta;
  arma::mat r = Y - X2 * beta2.t();
  arma::mat emp = em_with_zero_mean_c(r,100);
  arma::cube res = rinvwish_c(1, n+nu, emp*n + Phi*nu);
  return res.slice(0);
}


// [[Rcpp::export]]
Rcpp::List update_gamma_random_c(const arma::vec X,
                                 const arma::mat Y,
                                 const arma::vec gam){
  arma::vec prob = arma::ones<arma::vec>(gam.size());
  NumericVector prob2 = wrap(prob);
  int changeind = sample_index(gam.size(), prob2) - 1;
  arma::vec newgam = gam;
  newgam(changeind) = abs(gam(changeind) - 1);
  // int change = newgam(changeind);
  return(
    Rcpp::List::create(
      Rcpp::Named("gam") = newgam,
      Rcpp::Named("changeind") = changeind)
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec betagam_accept_random_c(const arma::vec X,
                                  const arma::mat Y,
                                  const double sigmabeta1,
                                  const arma::mat inputSigma,
                                  const double Vbeta,
                                  const arma::vec gam1,
                                  const arma::vec beta1,
                                  const arma::vec gam2,
                                  const arma::vec beta2,
                                  const int changeind,
                                  const int change,
                                  const int aa,
                                  const int bb){
  double newtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam2,beta2, aa, bb));
  double oldtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam1,beta1, aa, bb));
  double proposal_iter = R::dnorm(beta1(changeind)-beta2(changeind),0,sqrt(Vbeta),true);
  double proposal_ratio = proposal_iter;
  if(change==1){
      proposal_ratio = -proposal_iter;
    }else{
      proposal_ratio = proposal_iter;
    }
  double final_ratio = newtarget-oldtarget+proposal_ratio;
  arma::vec out = arma::zeros<arma::vec>(4);
  out(0) = final_ratio;
  out(1) = newtarget;
  out(2) = oldtarget;
  out(3) = proposal_ratio;
  return(out);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_betagam_random_c(const arma::vec X,
                                   const arma::mat Y,
                                   arma::vec gam1,
                                   arma::vec beta1,
                                   const arma::mat Sigma,
                                   const double sigmabeta,
                                   const double Vbeta,
                                   const int bgiter,
                                   const double smallchange,
                                   const int aa,
                                   const int bb){
  arma::vec gam2;
  arma::vec beta2;
  arma::vec A;
  for (int i=1; i<bgiter; ++i){
    Rcpp::List temp = update_gamma_random_c(X,Y,gam1);
    arma::vec gam2 = as<arma::vec>(temp["gam"]);
    int changeind = temp["changeind"];
    int change = gam2(changeind);
    arma::vec beta2 = beta1 % gam2;
    if(change==1){
      arma::uvec ind = find(gam1==1);
      beta2(ind) = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(smallchange)));
      beta2(changeind) = beta2(changeind) + (rnorm(1, 0, sqrt(Vbeta)))[0];
    }else{
      arma::uvec ind = find(gam2==1);
      beta2(ind) = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(smallchange)));
      beta2(changeind) = 0;
    }

    A = betagam_accept_random_c(X,Y,sigmabeta,
                                Sigma,Vbeta,
                                gam1,beta1,
                                gam2,beta2,
                                changeind,change,
                                aa, bb);

    NumericVector check2 = runif(1);
    double check = check2(0);

    if(exp(A(0)) > check){
      gam1 = gam2; beta1 = beta2;
    }

  }
  return Rcpp::List::create(
    Rcpp::Named("gam")= gam1,
    Rcpp::Named("beta") = beta1,
    Rcpp::Named("target") = A
  );
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_sigmabeta_from_h_c(const double h,
                              const arma::vec gam,
                              const arma::mat Sigma,
                              const arma::vec X,
                              const int T){
  //convert h to sigmabeta conditioning on gamma and Sigma
  arma::uvec ind = find(gam == 1);
  arma::vec ds = Sigma.diag();
  int n = X.size();
  double num = h * sum(ds);
  double denom = (1-h)* sum(X%X)/n * ind.size() + 1;
  return num/denom;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_h_from_sigmabeta_c(const arma::vec X,
                              const double sigmabeta,
                              const arma::mat Sigma,
                              const arma::vec gam,
                              const int n,
                              const int T){
  //converts sigmabeta to h conditioning on gamma and Sigma
  arma::uvec ind = find(gam==1);
  arma::vec ds = Sigma.diag();
  double num = sum(X%X)/n * ind.size() * sigmabeta;
  double denom = num + sum(ds);
  return num/denom;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_h_c(const double initialh,
                      const int hiter,
                      const arma::vec gam,
                      const arma::vec beta,
                      const arma::mat Sig,
                      const arma::vec X,
                      int T){
  double h1 = initialh;
  double sigbeta1 = get_sigmabeta_from_h_c(initialh, gam, Sig, X, T);
  arma::vec lik = arma::zeros<arma::vec>(hiter);
  arma::vec ds = Sig.diag();
  for (int i=1; i<hiter; ++i){
    double h2 = h1;
    NumericVector(rr) = runif(1, -0.05, 0.05);
    double r = rr(0);
    h2 = h2 + r;
    if(h2<0.05){h2 = abs(0.1 - h2);}
    if(h2>0.9){h2 = 1.8-h2;}
    arma::uvec ind = find(gam==1);
    double sigmabeta1 = get_sigmabeta_from_h_c(h1, gam, Sig, X, T);
    double sigmabeta2 = get_sigmabeta_from_h_c(h2, gam, Sig, X, T);
    double lik1 = 0; double lik2 = 0;
    for (int j=0; j < ind.size(); ++j){
      int newind = ind(j);
      lik1 = lik1 + R::dnorm(beta(newind), 0, sqrt(sigmabeta1*ds(newind)), true);
      lik2 = lik2 + R::dnorm(beta(newind), 0, sqrt(sigmabeta2*ds(newind)), true);
    }
    double acceptanceprob = exp(lik2-lik1);
    arma::vec ee = runif(1);
    double e = ee(0);
    if(e<acceptanceprob){
      h1 = h2; sigbeta1 = sigmabeta2;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("h") = h1,
    Rcpp::Named("sigbeta") = sigbeta1
  );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List run2chains_c(const arma::vec X,
                        const arma::mat Y,
                        const Rcpp::List initial_chain1,//const Rcpp::List initial_chain2,
                        const arma::mat Phi,
                        const double sigmabeta,
                        const int niter   = 1000,
                        const int bgiter  = 500,
                        const int hiter   = 50,
                        const int burnin  = 100000,
                        const int Vbeta   = 1,
                        const int aa = 1,
                        const int bb = 10,
                        const double smallchange = 1e-2){

  //initialize if not user-defined
  int T = Y.n_cols;
  int n = Y.n_rows;
  int nu = n;

  arma::mat outbeta1 = arma::zeros<arma::mat>(T, niter);
  arma::mat outgam1 = arma::zeros<arma::mat>(T,niter);
  arma::cube outSigma1 = arma::zeros<arma::cube>(T,T,niter);
  arma::vec outsb1 = arma::zeros<arma::vec>(niter);
  arma::vec outh1 = arma::zeros<arma::vec>(niter);
  arma::mat tar1 = arma::zeros<arma::mat>(3, niter);

  outbeta1.col(0)    = as<arma::vec>(initial_chain1["beta"]);
  outgam1.col(0)     = as<arma::vec>(initial_chain1["gamma"]);
  outSigma1.slice(0) = as<arma::mat>(initial_chain1["Sigma"]);
  outsb1(0)          = initial_chain1["sigmabeta"];
  outh1(0)            = 0.5;

  for (int i=1; i<niter; ++i){
    Rcpp::List bg = update_betagam_random_c(X,
                                            Y,
                                            outgam1.col(i-1),
                                            outbeta1.col(i-1),
                                            outSigma1.slice(i-1),
                                            outsb1(i-1),
                                            Vbeta,
                                            bgiter,
                                            smallchange, aa, bb);
    outgam1.col(i)  = as<arma::vec>(bg["gam"]);
    outbeta1.col(i) = as<arma::vec>(bg["beta"]);
    outSigma1.slice(i) = update_Sigma_c(n,nu,X,outbeta1.col(i),Phi,Y);
    Rcpp::List hsig = update_h_c(outh1[i-1],
                                 hiter,
                                 outgam1.col(i),
                                 outbeta1.col(i),
                                 outSigma1.slice(i),
                                 X,
                                 T);
    outh1(i) = hsig["h"];
    outsb1(i) = hsig["sigbeta"];
    tar1.col(i) = get_target_c(X,
             Y,
             outsb1(i),
             outSigma1.slice(i),
             outgam1.col(i),
             outbeta1.col(i),
             aa,
             bb);
    // cout << tar1.col(i).t() << "\n";
    if(i%10==0){
      cout << outgam1.col(i).t() << "\n";
      cout << i << "\n";
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("gamma") = outgam1.t(),
    Rcpp::Named("beta") = outbeta1.t(),
    Rcpp::Named("Sigma") = outSigma1,
    Rcpp::Named("sigmabeta") = outsb1,
    Rcpp::Named("h") = outh1,
    Rcpp::Named("target") = tar1.t()
  );

}

