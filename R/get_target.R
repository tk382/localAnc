get_target = function(X, Y, sigmabeta, Sigma, gam, beta, T){
  L=0
  for (i in 1:n){
    ind = which(!is.na(Y[i,]))
    if(sum(ind)>0){
      thing.to.add = dmvnrm_arma(Y[i,ind],
                                 X[i]*beta[ind],
                                 matrix(Sigma[ind,ind], nrow=length(ind)),
                                 T,
                                 logd=TRUE)/length(ind)*T
      L = L + thing.to.add
    }
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
