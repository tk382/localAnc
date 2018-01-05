get_target_new = function(X, Y, sigmabeta, Sigma, gam, beta, T){
  L = 0
  for (i in 1:n){
    ind = which(!is.na(Y[i,]))
    if(sum(ind)>0){
      Y[i, -ind] = mvrnormArma(1, rep(0, T-length(ind)), Sigma[-ind,-ind])
    }
    L = L+dmvnrm_arma(Y[i, ], X[i]*beta, Sigma, T, logd=TRUE)
  }
  ind = which(gam==1);
  s = length(ind);
  B = 0
  if(s>0){
    B = sum(dnorm(beta[ind], rep(0,s),
                  sqrt(sigmabeta*(diag(Sigma)[ind])),
                  log=TRUE))
  }
  G = log(beta(s+1, T-s+1))
  return(c(L,B,G))
}
