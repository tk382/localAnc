get_target_new = function(X, YY, sigmabeta, Sigma, gam, beta, T){
  L = 0
  Y = YY
  for (i in 1:n){
    ind = which(is.na(Y[i,]))
    if(length(ind)<T){
      mu1 = rep(0, length(ind))
      Sigma11 = Sigma[ind,ind]; Sigma12 = Sigma[ind, -ind]; Sigma21 = Sigma[-ind,ind]
      Sigma22 = Sigma[-ind,-ind]
      Y[i,ind] = mvrnormArma(1, mu1, Sigma11-Sigma12 %*% solve(Sigma22) %*% Sigma21)
      # Y[i, -ind] = mvrnormArma(1, rep(0, T-length(ind)), Sigma[-ind,-ind])
    }
    #print(dmvnrm_arma(Y[i, ], X[i]*beta, Sigma, T, logd=TRUE))
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
