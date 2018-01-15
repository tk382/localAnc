get_missing_likelihood = function(X,Y,beta,Sigma,T){
  L=rep(0, n)
  for (i in 1:n){
    ind = which(!is.na(Y[i,]))
    if(sum(ind)>0){
      thing.to.add = dmvnrm_arma(Y[i,ind],
                                 X[i]*beta[ind],
                                 matrix(Sigma[ind,ind], nrow=length(ind)),
                                 T,
                                 logd=TRUE)

      L[i] = thing.to.add
    }
  }
  return(L)
  # L = rep(0,length(X))
  # for (i in 1:n){
  #   ind = which(is.na(Y[i,]))
  #   if(length(ind)<T){
  #     Sigma11 = Sigma[ind,ind]; Sigma12 = Sigma[ind, -ind]; Sigma21 = Sigma[-ind,ind]
  #     Sigma22 = Sigma[-ind,-ind]
  #     mu1 = beta[ind] + Sigma12 %*% solve(Sigma22) %*% (Y[i,-ind]-beta[-ind])
  #
  #     Y[i,ind] = mvrnormArma(1, mu1, Sigma11-Sigma12 %*% solve(Sigma22) %*% Sigma21)
  #     # Y[i, -ind] = mvrnormArma(1, rep(0, T-length(ind)), Sigma[-ind,-ind])
  #   }
  #   #print(dmvnrm_arma(Y[i, ], X[i]*beta, Sigma, T, logd=TRUE))
  #   L[i] = dmvnrm_arma(Y[i, ], X[i]*beta, Sigma, T, logd=TRUE)
  # }
  # return(L)

}
