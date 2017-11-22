#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List em_with_zero_mean_c(arma::mat y, int maxit){
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
        bias(nind, nind) += Sigma(nind, nind) - Sigma(nind, ind) * (Sigma(ind, ind).i()) * Sigma(ind, nind);
        arma::uvec rowind = arma::zeros<arma::uvec>(1);
        rowind(0) = i;
        arma::mat yvec = y(rowind, ind);
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
  return Rcpp::List::create(
    Rcpp::Named("Sigma") = finalSigma,
    Rcpp::Named("iteration") = it
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(Sigma);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double dmvnrm_arma(arma::rowvec x, arma::rowvec mean, arma::mat sigma, bool logd = false) {
  int xdim = x.n_cols;
  if(xdim==0){return 0;}
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  arma::vec z = rooti*arma::trans(x-mean);
  out = constants - 0.5*arma::sum(z%z)+rootisum;
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_sigmabeta_from_h_c(double h, arma::vec gam, arma::mat Sigma, arma::vec X, int T){
  int n = X.size();
  arma::vec ds = Sigma.diag();
  double num = h * sum(ds);
  arma::uvec ind = find(gam == 1);
  double denom = (1-h)*sum(ds(ind)) * sum(X%X)/n;
  return num/denom;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_h_from_sigmabeta_c(arma::vec X, double sigmabeta,
                              arma::mat Sigma, arma::vec gam,
                              int n, int T){
  arma::uvec ind = find(gam==1);
  arma::vec ds = Sigma.diag();
  double num = sum(X%X)/n * sum(ds(ind)) * sigmabeta;
  double denom = num + sum(ds);
  return num/denom;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec get_target_c(arma::vec X, arma::mat Y, double sigmabeta,
                       arma::mat Sigma, arma::vec gam, arma::vec beta){
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
      B = B + R::dnorm(beta(newind), 0, sqrt(sigmabeta*ds(newind)), true);
    }
  }
  G = log(gamma(s+1)*gamma(T-s+1)/gamma(T+2));
  arma::vec out = arma::zeros<arma::vec>(3);
  out(0) = L;
  out(1) = B;
  out(2) = G;
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
int sample_index(int size, NumericVector prob = NumericVector::create()){
  arma::vec sequence = arma::linspace<arma::vec>(1, size, size);
  arma::vec out = Rcpp::RcppArmadillo::sample(sequence, size, false, prob);
  return out(0);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_gamma_c(arma::vec X, arma::mat Y, arma::vec gam){
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
    arma::rowvec marcor = arma::zeros<arma::rowvec>(ind0.size());
    for (int t=0; t<ind0.size(); ++t){
      arma::vec temp = Y.col(ind0(t));
      arma::uvec tempind = find_finite(temp);
      marcor(t) = abs(sum(temp(tempind)%X(tempind)))/tempind.size();
    }
    int add = 1;
    if(s<(T-1)){
      NumericVector marcor2 = wrap(marcor);
      add = sample_index(ind0.size(), marcor2);
    }
    newgam(ind0(add-1)) = 1;
    changeind = ind0(add-1);
  }
  if(cas==2){
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
arma::vec betagam_accept_c(arma::vec X,
                         arma::mat Y,
                         double sigmabeta1,
                         arma::mat inputSigma,
                         double Vbeta,
                         arma::vec gam1,
                         arma::vec beta1,
                         arma::vec gam2,
                         arma::vec beta2,
                         int changeind,
                         int change){
  double newtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam2,beta2));
  double oldtarget = sum(get_target_c(X,Y,sigmabeta1,inputSigma,gam1,beta1));
  double proposal_ratio = R::dnorm(beta1(changeind)-beta2(changeind),0,sqrt(Vbeta),true);
  int T = gam1.size();
  int s1 = sum(gam1==1);
  int s2 = sum(gam2==1);
  arma::rowvec marcor = arma::zeros<arma::rowvec>(T);
  for (int t=0; t<T; ++t){
    arma::vec temp = Y.col(t);
    arma::uvec tempind = find_finite(temp);
    marcor(t) = abs(sum(temp(tempind)%X(tempind)))/tempind.size();
  }
  if(change==1){
    arma::uvec ind1 = find(gam1==0);
    double temp1 = marcor(changeind)/sum(marcor(ind1));
    proposal_ratio = -log(temp1)-log(s2)-proposal_ratio;
  }else{
    arma::uvec ind2 = find(gam2==0);
    double temp2 = marcor(changeind)/sum(marcor(ind2));
    proposal_ratio = log(temp2)+log(s1)+proposal_ratio;
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
Rcpp::List update_betagam_c(arma::vec X,
                           arma::mat Y,
                           arma::vec gam1,
                           arma::vec beta1,
                           arma::mat Sigma,
                           double sigmabeta,
                           double Vbeta,
                           int bgiter){
  int T = gam1.size();
  arma::mat outgamma = arma::zeros<arma::mat>(T,bgiter);
  arma::mat outbeta = arma::zeros<arma::mat>(T,bgiter);
  outgamma.col(0) = gam1;
  outbeta.col(0) = beta1;
  arma::vec tar = arma::zeros<arma::vec>(bgiter);
  for (int i=1; i<bgiter; ++i){
    Rcpp::List temp = update_gamma_c(X,Y,outgamma.col(i-1));
    arma::vec gam1 = outgamma.col(i-1);
    arma::vec beta1 = outbeta.col(i-1);
    arma::vec gam2 = as<arma::vec>(temp["gam"]);
    arma::vec beta2 = beta1 % gam2;
    arma::uvec ind = find(gam2==1);
    beta2(ind) = beta1(ind) + as<arma::vec>(rnorm(ind.size(), 0, sqrt(Vbeta)));
    int changeind = temp["changeind"];
    int change = gam2(changeind);
    arma::vec A = betagam_accept_c(X,Y,sigmabeta,
                                   Sigma,Vbeta,
                                   gam1,beta1,
                                   gam2,beta2,
                                   changeind,change);
    NumericVector check2 = runif(1);
    double check = check2(0);
    if(exp(A(0))>check){
      tar(i) = A(1);
      outgamma.col(i) = gam2; outbeta.col(i) = beta2;
    }else{
      tar(i) = A(2);
      outgamma.col(i) = gam1; outbeta.col(i) = beta1;
    }
  }
  return Rcpp::List::create(
      Rcpp::Named("gam") = outgamma.t(),
      Rcpp::Named("beta") = outbeta.t(),
      Rcpp::Named("tar") = tar
    );
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List update_h_c(double initialh, int hiter, arma::vec gam, arma::vec beta,
                  arma::mat Sig, arma::vec X, int T){
  arma::vec outh = arma::zeros<arma::vec>(hiter);
  outh(0) = initialh;
  arma::vec outsigbeta = arma::zeros<arma::vec>(hiter);
  outsigbeta(0) = get_sigmabeta_from_h_c(initialh, gam, Sig, X, T);
  arma::vec lik = arma::zeros<arma::vec>(hiter);
  arma::vec ds = Sig.diag();
  for (int i=1; i<hiter; ++i){
    double h1 = outh(i-1);
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
      outh(i) = h2; outsigbeta(i) = sigmabeta2; lik(i) = lik2;
    }else{
      outh(i) = h1; outsigbeta(i) = sigmabeta1; lik(i) = lik1;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("h") = outh,
    Rcpp::Named("sigbeta") = outsigbeta,
    Rcpp::Named("lik") = lik
  );
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::cube rinvwish_c(int n, int v, arma::mat S){
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
arma::mat update_Sigma_c(int n, int nu, arma::vec X, arma::vec beta, arma::mat Phi, arma::mat Y){
  int T = Y.n_cols;
  arma::mat X2 = arma::zeros<arma::mat>(n, 1);
  X2.col(0) = X;
  arma::mat beta2 = arma::zeros<arma::mat>(T,1);
  beta2.col(0) = beta;
  arma::mat r = Y - X2 * beta2.t();
  Rcpp::List li = em_with_zero_mean_c(r, 100);
  arma::mat emp = as<arma::mat>(li["Sigma"]);
  arma::cube res = rinvwish_c(1, n+nu, emp*n + Phi*nu);
  return res.slice(0);
}

