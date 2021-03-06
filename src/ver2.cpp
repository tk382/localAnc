#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <unistd.h>
using namespace Rcpp;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat em_with_zero_mean_c(arma::mat y,
                              const int maxit){
  //EM for empirical covariance matrix when y has missing values
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
  return finalSigma;
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat em_with_zero_mean_c_needs_debugging(arma::mat y,
                              const int maxit){
  //EM for empirical covariance matrix when y has missing values
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
  arma::mat Sigma = oldSigma;   double diff = 1;   int it = 1;
  while (diff>0.001 && it < maxit){
    arma::mat bias = arma::zeros<arma::mat>(p,p);
    for (int i=0; i<n; ++i){
      arma::rowvec tempdat = y.row(i);
      arma::uvec ind = find_finite(tempdat);
      arma::uvec nind = find_nonfinite(tempdat);
      if (0 < ind.size() && ind.size() < p){
        arma::mat rooti = arma::inv(trimatu(arma::chol(Sigma(ind,ind))));
        arma::mat z = Sigma(ind,nind)*rooti;
        bias(nind, nind) += Sigma(nind, nind) - z*z.t();
        arma::uvec rowind = arma::zeros<arma::uvec>(1);
        rowind(0) = i;   arma::mat yvec = y(rowind, ind);
        y_imputed(rowind, nind) = y(rowind, ind) * z;
      }
    }
    Sigma = (y_imputed.t() * y_imputed + bias)/n;
    arma::mat diffmat = (Sigma-oldSigma);
    arma::mat diffsq = diffmat%diffmat;
    diff = accu(diffsq);    oldSigma = Sigma;    it = it + 1;
  }
  arma::mat finalSigma = arma::zeros<arma::mat>(orig_p, orig_p);
  finalSigma.submat(valid_ind, valid_ind.t()) = Sigma;
  return finalSigma;
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
                   const int xdim,
                   bool logd = false) {
  //returns the density of a multivariate normal vector
  //input : a rowvector x whose density you'd like to know
  //      : a rowvector mean for the mean of mvn
  //      : a matrix Sigma for covariance, needs to be psd
  //      : a boolean logd, true if you like the log density
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  arma::vec z = rooti*arma::trans(x-mean);
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  out = constants - 0.5*arma::sum(z%z)+rootisum;
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_sigmabeta_from_h_c(const double h,
                              const arma::vec gam,
                              const arma::mat Sigma,
                              const arma::vec X,
                              const int T){
  //convert h to sigmabeta conditioning on gamma and Sigma
  int n = X.size();
  arma::vec ds = Sigma.diag();
  double num = h * sum(ds);
  arma::uvec ind = find(gam == 1);
  double denom = (1-h)*sum(ds(ind)) * sum(X%X)/n;
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
  double num = sum(X%X)/n * sum(ds(ind)) * sigmabeta;
  double denom = num + sum(ds);
  return num/denom;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec get_target_c(const arma::vec X,
                       const arma::mat Y,
                       const double sigmabeta,
                       const arma::mat Sigma,
                       const arma::vec gam,
                       const arma::vec beta,
                       double adjust){
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
                          T,
                          true);
    }
  }
  arma::uvec ind = find(gam==1);
  int s = ind.size();
  if(s>0){
    arma::vec ds = Sigma.diag();
    for (int j=0; j<s; ++j){
      int newind = ind(j);
      B = B + R::dnorm(beta(newind), 0, sqrt(sigmabeta*ds(newind)), true);
    }
  }
  G = log(std::tgamma(s+1)*std::tgamma(T-s+1)/std::tgamma(T+2));
  arma::vec out = arma::zeros<arma::vec>(3);
  //adjust for the missing values
  L = L * adjust;
  out(0) = L;
  out(1) = B;
  out(2) = G;
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec get_target_c_new(const arma::vec X,
                           arma::mat Y,
                           const double sigmabeta,
                           const arma::mat Sigma,
                           const arma::vec gam,
                           const arma::vec beta){
  //get the target likelihood circumventing the missing value issue
  int T = Y.n_cols;
  int n = Y.n_rows;
  double L = 0;
  double B = 0;
  double G = 0;
  //arma::mat tmp;
  for (int i=0; i < n; ++i){
    arma::uvec naind1 = find_nonfinite(Y.row(i).t());
    arma::uvec naind = find_finite(Y.row(i).t());
    if(naind.size()>0){
      arma::uvec rowind = arma::zeros<arma::uvec>(1);
      rowind(0) = i;
      //arma::mat tmp = mvrnormArma(1,arma::zeros<arma::vec>(naind1.size()),Sigma(naind1,naind1.t()));

      arma::mat Sigma11 = Sigma(naind1, naind1);
      arma::mat Sigma12 = Sigma(naind1, naind);
      arma::mat Sigma21 = Sigma(naind, naind1);
      arma::mat Sigma22 = Sigma(naind, naind);
      arma::vec mu1 = beta(naind1) + (Sigma12 * Sigma22.i() * (beta(naind).t()-Y.submat(rowind, naind.t())).t());

      arma::mat tmp = mvrnormArma(1,mu1,Sigma11-Sigma12*Sigma22.i()*Sigma21);
      Y.submat(rowind, naind1.t()) = tmp;
      //arma::rowvec Ytemp = Y(rowind, naind.t());
      L = L + dmvnrm_arma(Y.row(i),
                          X(i)*beta.t(),
                          Sigma,
                          T,
                          true);
    }
  }
  arma::uvec ind = find(gam==1);
  int s = ind.size();
  if(s>0){
    arma::vec ds = Sigma.diag();
    for (int j=0; j<s; ++j){
      int newind = ind(j);
      B = B + R::dnorm(beta(newind), 0, sqrt(sigmabeta*ds(newind)), true);
    }
  }
  G = log(std::tgamma(s+1)*std::tgamma(T-s+1)/std::tgamma(T+2));
  arma::vec out = arma::zeros<arma::vec>(3);
  out(0) = L;
  out(1) = B;
  out(2) = G;
  return out;
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
Rcpp::List update_h_c(const double initialh,
                      const int hiter,
                      const arma::vec gam,
                      const arma::vec beta,
                      const arma::mat Sig,
                      const arma::vec X,
                      const int T){
  double h1 = initialh;
  double sigbeta1 = get_sigmabeta_from_h_c(initialh, gam, Sig, X, T);
//  arma::vec lik = arma::zeros<arma::vec>(hiter);
  arma::vec ds = Sig.diag();
  if(all(beta==0)){
    return Rcpp::List::create(
      Rcpp::Named("h") = 0,
      Rcpp::Named("sigbeta") = 0
    );
  }
  for (int i=1; i<hiter; ++i){
    double h2 = h1;
    NumericVector(rr) = runif(1, -0.1, 0.1);
    double r = rr(0);
    h2 = h2 + r;
    if(h2<0){h2 = abs(h2);}
    if(h2>1){h2 = 2-h2;}
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
arma::vec betagam_accept_c(const arma::vec X,
                           const arma::mat Y,
                           const double sigmabeta1,
                           const arma::mat inputSigma,
                           const double Vbeta,
                           const arma::rowvec marcor,
                           const arma::vec gam1,
                           const arma::vec beta1,
                           const arma::vec gam2,
                           const arma::vec beta2,
                           const int changeind,
                           const int change,
                           const double adjust){
  //compute the target likelihood and the proposal ratio
  //to decide if you should accept the proposed beta and gamma
  double newtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam2,beta2, adjust));
  double oldtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam1,beta1, adjust));
  double proposal_ratio = R::dnorm(beta1(changeind)-beta2(changeind),0,sqrt(Vbeta),true);
  int s1 = sum(gam1==1);
  int s2 = sum(gam2==1);
  if(change==1){
    arma::uvec ind1 = find(gam1==0);
    double temp1 = marcor(changeind)/sum(marcor(ind1));
    proposal_ratio = -log(temp1)-log(s2)-proposal_ratio;
  }else{
    arma::uvec ind2 = find(gam2==0);
    double temp2 = marcor(changeind)/sum(marcor(ind2));
    proposal_ratio = log(temp2)+log(s1)+proposal_ratio;
  }
  //double final_ratio = newtarget-oldtarget+proposal_ratio;
  double final_ratio = newtarget-oldtarget+proposal_ratio;
  //cout << "gamma 1    : " << gam1.t() << "\n";
  //cout << "gamma 2    : " << gam2.t() << "\n";
  //cout << "beta  1    : " << beta1.t() << "\n";
  //cout << "beta  2    : " << beta2.t() << "\n";
  //cout << "new target : " << newtarget<< "\n";
  //cout << "old target : " << oldtarget<< "\n";
  //cout << "proposal   : " << proposal_ratio << "\n";
  arma::vec out = arma::zeros<arma::vec>(4);
  out(0) = final_ratio;
  //out(0) = newtarget-oldtarget;
  out(1) = newtarget;
  out(2) = oldtarget;
  out(3) = proposal_ratio;
  return(out);
}
//

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_gamma_c(const arma::vec X,
                          const arma::mat Y,
                          const arma::vec gam,
                          const arma::rowvec marcor){
  //update gamma once
  int changeind = 0;
  arma::vec newgam = gam;
  int T = gam.size();
  arma::vec prob = arma::zeros<arma::vec>(2);
  prob(0) = 0.5; prob(1) = 0.5;
  arma::uvec ind0 = find(gam==0);
  arma::uvec ind1 = find(gam==1);
  int s = ind1.size();
  NumericVector prob2 = wrap(prob);
  int cas = sample_index(2, prob2);
  if(s==0){
    cas = 1;
  }else if(s==T){
    cas = 2;
  }
  if (cas==1){
    arma::vec marcort = marcor.t();
    arma::vec rel_marcor = marcort(ind0);
    int add = 1;
    if(s<(T-1)){
      NumericVector marcor2 = wrap(rel_marcor);
      add = sample_index(ind0.size(), marcor2);
    }
    newgam(ind0(add-1)) = 1;
    changeind = ind0(add-1);
  }else if(cas==2){
    int remove = sample_index(s);
    remove = ind1(remove-1);
    newgam(remove) = 0;
    changeind = remove;
  }
  return(
    Rcpp::List::create(
      Rcpp::Named("gam") = newgam,
      Rcpp::Named("changeind") = changeind)
  );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_betagam_c(const arma::vec X,
                            const arma::mat Y,
                            arma::vec gam1,
                            arma::vec beta1,
                            const arma::mat Sigma,
                            const arma::rowvec marcor,
                            const double sigmabeta,
                            const double Vbeta,
                            const int bgiter,
                            const double adjust){
  //update and beta and gamma 'bgiter' times
  for (int i=1; i<bgiter; ++i){

    Rcpp::List temp = update_gamma_c(X,Y,gam1,marcor);
    arma::vec gam2 = as<arma::vec>(temp["gam"]);
    arma::vec beta2 = beta1 % gam2;
    arma::uvec ind = find(gam2==1);
    beta2(ind) = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(Vbeta)));
    int changeind = temp["changeind"];
    int change = gam2(changeind);
    arma::vec A = betagam_accept_c(X,Y,sigmabeta,
                                   Sigma,Vbeta,marcor,
                                   gam1,beta1,
                                   gam2,beta2,
                                   changeind,
                                   change,
                                   adjust);
    NumericVector check2 = runif(1);
    double check = check2(0);
    if(exp(A(0))>check){
      gam1 = gam2; beta1 = beta2;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("gam") = gam1,
    Rcpp::Named("beta") = beta1
  );
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

// // [[Rcpp::depends("RcppArmadillo")]]
// // [[Rcpp::export]]
// arma::vec betagam_accept_sw_c(const arma::vec X,
//                               const arma::mat Y,
//                               const double sigmabeta1,
//                               const arma::mat inputSigma,
//                               const double Vbeta,
//                               const arma::vec gam1,
//                               const arma::vec beta1,
//                               const arma::vec gam2,
//                               const arma::vec beta2,
//                               const int changeind,
//                               const int change,
//                               const int T,
//                               const arma::rowvec marcor){
//   double newtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam2,beta2));
//   double oldtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam1,beta1));
//   double proposal_ratio = R::dnorm(beta1(changeind)-beta2(changeind),0,sqrt(Vbeta),true);
//   arma::rowvec marcor2 = min(marcor)+max(marcor)-marcor;
//   if(change==1){
//     arma::uvec ind1 = find(gam1==0);
//     arma::uvec ind2 = find(gam2==1);
//     double tempadd = marcor(changeind)/sum(marcor(ind1));
//     double tempremove = marcor2(changeind)/sum(marcor2(ind2));
//     proposal_ratio = -log(tempadd)+log(tempremove)-proposal_ratio;
//   }else{
//     arma::uvec ind1 = find(gam1==1);
//     arma::uvec ind2 = find(gam2==0);
//     double tempadd = marcor(changeind)/sum(marcor(ind2));
//     double tempremove = marcor2(changeind) / sum(marcor2(ind1));
//     proposal_ratio = log(tempadd)-log(tempremove)+proposal_ratio;
//   }
//   double final_ratio = newtarget-oldtarget+proposal_ratio;
//   arma::vec out = arma::zeros<arma::vec>(4);
//   out(0) = final_ratio;
//   out(1) = newtarget;
//   out(2) = oldtarget;
//   out(3) = proposal_ratio;
//   return(out);
// }


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
IntegerVector csample_integer( IntegerVector x, int size, bool replace,
                               NumericVector prob = NumericVector::create()) {
  RNGScope scope;
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_gamma_sw_c(const arma::vec X,
                             const arma::mat Y,
                             const arma::vec gam,
                             const int T){
  IntegerVector where1 = csample_integer(wrap(arma::linspace(0, T-1)),round(T/2), false);
  arma::uvec where = as<arma::uvec>(where1);
  arma::vec newgam = gam;
  newgam(where) = -gam(where) + arma::ones<arma::vec>(where.size());
  return(
    Rcpp::List::create(
      Rcpp::Named("gam") = newgam,
      Rcpp::Named("changeind") = where)
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_betagam_sw_c(const arma::vec X,
                               const arma::mat Y,
                               arma::vec gam1,
                               arma::vec beta1,
                               const arma::mat Sigma,
                               const arma::rowvec marcor,
                               const double sigmabeta,
                               const double Vbeta,
                               const int bgiter,
                               const int T,
                               const double adjust){
  for (int i=1; i<bgiter; ++i){
    arma::vec gam2;
    arma::vec beta2;
    if(i%100==0){
      Rcpp::List temp = update_gamma_sw_c(X,Y,gam1,T);
      gam2 = as<arma::vec>(temp["gam"]);
      arma::uvec changeind = as<arma::uvec>(temp["changeind"]);
      beta2 = beta1 % gam2;
      arma::uvec ind = find(gam2==1);
      beta2(ind) = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(Vbeta)));
      arma::vec change = gam2(changeind);
      arma::vec tmpvar = arma::zeros<arma::vec>(T);
      double proposal_ratio = 0;
      for (int z=0; z < changeind.size(); ++z){
        proposal_ratio = proposal_ratio + R::dnorm(beta1(changeind(z)) - beta2(changeind(z)),0,sqrt(Vbeta),true);
      }
      double newtarget = sum(get_target_c(X,Y,sigmabeta,Sigma,gam2,beta2, adjust));
      double oldtarget = sum(get_target_c(X,Y,sigmabeta,Sigma,gam1,beta1, adjust));
      double A = newtarget - oldtarget + proposal_ratio;
      NumericVector check = runif(1);
      double check1 = check(0);
      if(exp(A)>check1){
        gam1 = gam2; beta1 = beta2;
      }
    }else{
      Rcpp::List temp = update_gamma_c(X,Y,gam1,marcor);
      gam2  = as<arma::vec>(temp["gam"]);
      beta2 = beta1 % gam2;
      arma::uvec ind  = find(gam2==1);
      beta2(ind)      = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(Vbeta)));
      int changeind   = temp["changeind"];
      int change      = gam2(changeind);
      arma::vec A     = betagam_accept_c(X, Y, sigmabeta,
                                         Sigma, Vbeta, marcor,
                                         gam1, beta1, gam2, beta2,
                                         changeind, change, adjust);
      NumericVector check = runif(1);
      double check2 = check(0);
      if(exp(A(0))>check2){
        gam1 = gam2; beta1 = beta2;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("gam")= gam1,
    Rcpp::Named("beta") = beta1
  );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List doMCMC_c(const arma::vec X,
                  const arma::mat Y,
                  const int n,
                  const int T,
                  const arma::mat Phi,
                  const int nu,
                  const arma::vec initialbeta,
                  const arma::vec initialgamma,
                  const arma::mat initialSigma,
                  const double initialsigmabeta,
                  const arma::rowvec marcor,
                  const double Vbeta,
                  const int niter,
                  const int bgiter,
                  const int hiter,
                  const int switer,
                  const double adjust){
  //empty arrays to save values
  arma::mat outbeta = arma::zeros<arma::mat>(T, niter);
  arma::mat outgam = arma::zeros<arma::mat>(T,niter);
  arma::cube outSigma = arma::zeros<arma::cube>(T,T,niter);
  arma::vec outsb = arma::zeros<arma::vec>(niter);
  arma::vec outh = arma::zeros<arma::vec>(niter);
  arma::mat tar = arma::zeros<arma::mat>(3, niter);
  //initialize
  outbeta.col(0) = initialbeta;
  outgam.col(0) = initialgamma;
  outSigma.slice(0) = initialSigma;
  outsb(0) = initialsigmabeta;
  tar.col(0) = arma::zeros<arma::vec>(3);
  outh(0) = get_h_from_sigmabeta_c(X,outsb(0),outSigma.slice(0),
       outgam.col(0), n, T);
  for (int i=1; i<niter; ++i){
    arma::vec gam1    = outgam.col(i-1);
    arma::vec beta1   = outbeta.col(i-1);
    arma::mat Sigma1  = outSigma.slice(i-1);
    double sigmabeta1 = outsb(i-1);
    double h1         = outh(i-1);
    Rcpp::List bg = update_betagam_c(X,Y,
                                     gam1,beta1,
                                     Sigma1, abs(marcor),
                                     sigmabeta1,Vbeta,
                                     bgiter,
                                     adjust);
    arma::vec gam2  = as<arma::vec>(bg["gam"]);
    arma::vec beta2 = as<arma::vec>(bg["beta"]);
    arma::mat Sigma2 = update_Sigma_c(n,nu,X,beta2,Phi,Y);
    Rcpp::List hsig = update_h_c(h1, hiter,
                                 gam2, beta2,
                                 Sigma2, X, T);
    outh(i) = hsig["h"];
    outsb(i) = hsig["sigbeta"];
    outgam.col(i) = gam2;
    outbeta.col(i) = beta2;
    outSigma.slice(i) = Sigma2;
    if(!arma::is_finite(outsb(i))){
      outsb(i) = 1000;
    }
    tar.col(i) = get_target_c(X,Y,outsb(i), Sigma2,gam2, beta2, adjust);
    cout << i << "\n";
  }
  return Rcpp::List::create(
    Rcpp::Named("gam") = wrap(outgam.t()),
    Rcpp::Named("beta") = wrap(outbeta.t()),
    Rcpp::Named("sigbeta") = wrap(outsb),
    Rcpp::Named("Sigma") = wrap(outSigma),
    Rcpp::Named("tar") = wrap(tar.t())
  );

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List run2chains_c(const arma::vec X,
             const arma::mat Y,
             const Rcpp::List initial_chain1,
             const Rcpp::List initial_chain2,
             const arma::mat Phi,
             const arma::rowvec marcor,
             const int niter = 1000,
             const int bgiter = 500,
             const int hiter = 50,
             const int switer = 50,
             const int burnin = 5,
             const double adjust = 1){
  //initialize if not user-defined
  int T = Y.n_cols;
  int n = Y.n_rows;
  int nu = T+5;

  //initialize Vbeta
  double Vbeta = mean(marcor%marcor) * 0.0005;

  arma::mat outbeta1 = arma::zeros<arma::mat>(T, niter);
  arma::mat outgam1 = arma::zeros<arma::mat>(T,niter);
  arma::cube outSigma1 = arma::zeros<arma::cube>(T,T,niter);
  arma::vec outsb1 = arma::zeros<arma::vec>(niter);
  arma::vec outh1 = arma::zeros<arma::vec>(niter);
  arma::mat tar1 = arma::zeros<arma::mat>(3, niter);

  arma::mat outbeta2 = arma::zeros<arma::mat>(T, niter);
  arma::mat outgam2 = arma::zeros<arma::mat>(T,niter);
  arma::cube outSigma2 = arma::zeros<arma::cube>(T,T,niter);
  arma::vec outsb2 = arma::zeros<arma::vec>(niter);
  arma::vec outh2 = arma::zeros<arma::vec>(niter);
  arma::mat tar2 = arma::zeros<arma::mat>(3, niter);

  outbeta1.col(0)    = as<arma::vec>(initial_chain1["beta"]);
  outgam1.col(0)     = as<arma::vec>(initial_chain1["gamma"]);
  outSigma1.slice(0) = as<arma::mat>(initial_chain1["Sigma"]);
  outsb1(0)          = initial_chain1["sigmabeta"];

  outbeta2.col(0)    = as<arma::vec>(initial_chain2["beta"]);
  outgam2.col(0)     = as<arma::vec>(initial_chain2["gamma"]);
  outSigma2.slice(0) = as<arma::mat>(initial_chain2["Sigma"]);
  outsb2(0)          = initial_chain2["sigmabeta"];

  for (int i=1; i<niter; ++i){
    //chain 1 update
    Rcpp::List bg = update_betagam_sw_c(X,
                                        Y,
                                        outgam1.col(i-1),
                                        outbeta1.col(i-1),
                                        outSigma1.slice(i-1),
                                        marcor,
                                        outsb1[i-1],
                                        Vbeta,
                                        bgiter,
                                        T,
                                        adjust);

    // Rcpp::List bg = update_betagam_c(X,Y,outgam1.col(i-1), outbeta1.col(i-1),
    //                       outSigma1.slice(i-1), marcor,
    //                       outsb1[i-1], Vbeta, bgiter);

    outgam1.col(i)  = as<arma::vec>(bg["gam"]);
    cout << outgam1.col(i).t() << "\n";
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
    if(!arma::is_finite(outsb1(i))){
      outsb1(i) = 1000;
    }
    tar1.col(i) = get_target_c(X,
                              Y,
                              outsb1(i),
                              outSigma1.slice(i),
                              outgam1.col(i),
                              outbeta1.col(i),
                              adjust);

    //chain 2 update
    bg = update_betagam_sw_c(X,
                             Y,
                             outgam2.col(i-1),
                             outbeta2.col(i-1),
                             outSigma2.slice(i-1),
                             marcor,
                             outsb1[i-1],
                             Vbeta,
                             bgiter,
                             T,
                             adjust);

    // bg = update_betagam_c(X,Y,outgam2.col(i-1), outbeta2.col(i-1),
    //                       outSigma2.slice(i-1), marcor,
    //                       outsb2[i-1], Vbeta, bgiter);


    outgam2.col(i)  = as<arma::vec>(bg["gam"]);
    cout << outgam2.col(i).t() << "\n";
    outbeta2.col(i) = as<arma::vec>(bg["beta"]);
    outSigma2.slice(i) = update_Sigma_c(n,nu,X,outbeta1.col(i),Phi,Y);
    hsig = update_h_c(outh1[i-1],
                                 hiter,
                                 outgam2.col(i),
                                 outbeta2.col(i),
                                 outSigma2.slice(i),
                                 X,
                                 T);
    outh2(i) = hsig["h"];
    outsb2(i) = hsig["sigbeta"];
    if(!arma::is_finite(outsb2(i))){
      outsb2(i) = 1000;
    }
    tar2.col(i) = get_target_c(X,
             Y,
             outsb2(i),
             outSigma2.slice(i),
             outgam2.col(i),
             outbeta2.col(i),
             adjust);

    //convergence criterion
    if(i>2*burnin && i%5==0){
      arma::vec rowmean1 = mean(outgam1.cols(burnin,i), 1);
      arma::vec rowmean2 = mean(outgam2.cols(burnin,i), 1);
      if(all(rowmean1<0.5) & all(rowmean2<0.5)){
        cout<< "both chains selected no variables - converged!";
        outSigma1.shed_slices(i+1, niter-1);  outSigma2.shed_slices(i+1,niter-1);
        outgam1.shed_cols(i+1,niter-1);       outgam2.shed_cols(i+1,niter-1);
        outbeta1.shed_cols(i+1, niter-1);     outbeta2.shed_cols(i+1, niter-1);
        outh1.shed_rows(i+1,niter-1);         outh2.shed_rows(i+1, niter-1);
        outsb1.shed_rows(i+1,niter-1);        outsb2.shed_rows(i+1,niter-1);
        break;
      }else{
        arma::uvec est1 = find(rowmean1 > 0.5);
        arma::uvec est2 = find(rowmean2 > 0.5);
        if(est1.size()==est2.size() && all(est1==est2)){
          arma::mat tmpbeta1 = outbeta1.cols(burnin,i);
          arma::mat tmpbeta2 = outbeta2.cols(burnin,i);
          arma::mat tmpgam1 = outgam1.cols(burnin,i);
          arma::mat tmpgam2 = outgam2.cols(burnin,i);
          double diff = 0;
          for (int k = 0; k < est1.size(); ++k){
            int kk = est1(k);
            arma::uvec ones = find(tmpgam1.row(kk)==1);
            arma::rowvec tmptmpbeta1 = tmpbeta1.row(kk);
            double beta1 = mean(tmptmpbeta1(ones));
            arma::uvec ones2 = find(tmpgam2.row(kk)==1);
            arma::rowvec tmptmpbeta2 = tmpbeta2.row(kk);
            double beta2 = mean(tmptmpbeta2(ones));
            diff = diff + (beta1-beta2)*(beta1-beta2);
          }
          if(diff/est1.size() < 1e-2){
            cout<< "beta difference is small between the two chains - converged!\n";
            outSigma1.shed_slices(i+1, niter-1);  outSigma2.shed_slices(i+1,niter-1);
            outgam1.shed_cols(i+1,niter-1);       outgam2.shed_cols(i+1,niter-1);
            outbeta1.shed_cols(i+1, niter-1);     outbeta2.shed_cols(i+1, niter-1);
            outh1.shed_rows(i+1,niter-1);         outh2.shed_rows(i+1, niter-1);
            outsb1.shed_rows(i+1,niter-1);        outsb2.shed_rows(i+1,niter-1);
            break;
          }
        }
      }
    }
    cout << i << "\n";
  }
  return Rcpp::List::create(
    Rcpp::Named("chain1") = Rcpp::List::create(
      Rcpp::Named("gamma") = outgam1.t(),
      Rcpp::Named("beta") = outbeta1.t(),
      Rcpp::Named("Sigma") = outSigma1,
      Rcpp::Named("sigmabeta") = outsb1,
      Rcpp::Named("h") = outh1,
      Rcpp::Named("target") =tar1
    ),
    Rcpp::Named("chain2") = Rcpp::List::create(
      Rcpp::Named("gamma") = outgam2.t(),
      Rcpp::Named("beta") = outbeta2.t(),
      Rcpp::Named("Sigma") = outSigma2,
      Rcpp::Named("sigmabeta") = outsb2,
      Rcpp::Named("h") = outh2,
      Rcpp::Named("target") = tar2
    )
  );
}

