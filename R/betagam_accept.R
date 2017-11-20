betagam_accept = function(X, Y, sigmabeta1, inputSigma, gam1, beta1, gam2, beta2, changeind, change){
  L_ratio_num = 0
  L_ratio_denom = 0
  for (i in 1:n){
    ind = !is.na(Y[i,])
    if(sum(ind)>0){
      L_ratio_num = L_ratio_num + dmvnrm_arma(Y[i,ind],
                    X[i]*beta2[ind],
                    matrix(inputSigma[ind,ind], nrow=sum(ind)),
                    logd=TRUE)
      L_ratio_denom = L_ratio_denom + dmvnrm_arma(Y[i,ind],
                      X[i]*beta1[ind],
                      matrix(inputSigma[ind,ind], nrow=sum(ind)),
                      logd=TRUE)
    }
  }
  L_ratio = L_ratio_num - L_ratio_denom

  ind1 = which(gam1==1); ind2 = which(gam2==1)
  s1 = length(ind1);     s2 = length(ind2)
  if(s1>1){
    Sigma1 = diag(diag(inputSigma)[ind1])
    B_ratio_num = dmvnrm_arma(beta1[ind1], rep(0,length(ind1)), sigmabeta1*Sigma1)
  }else if(s1==1){
    Sigma1 = matrix(diag(inputSigma)[ind1])
    B_ratio_num = dmvnrm_arma(beta1[ind1], rep(0,length(ind1)), sigmabeta1*Sigma1)
  }else{
    B_ratio_num = 0
  }
  if(s2>1){
    Sigma2 = diag(diag(inputSigma)[ind2])
    B_ratio_denom = dmvnrm_arma(beta2[ind2], rep(0,length(ind2)), sigmabeta1*Sigma2)
  }else if(s2==1){
    Sigma2 = matrix(diag(inputSigma)[ind2])
    B_ratio_denom = dmvnrm_arma(beta2[ind2], rep(0,length(ind2)), sigmabeta1*Sigma2)
  }else{
    B_ratio_denom = 0
  }
  B_ratio = B_ratio_num - B_ratio_denom


  G_ratio = log(beta(s2+1, T-s2+1)) - log(beta(s1+1, T-s1+1))

  if(change==1){
    beta1_temp = beta1[ind2]; beta2_temp = beta2[ind2]
  }else{
    beta1_temp = beta1[ind1]; beta2_temp = beta2[ind1]
  }

  proposal_ratio = sum(dnorm(beta1_temp-beta2_temp,0,Vbeta, log = TRUE))
  if(change=='add'){
    proposal_ratio = (T-s1)/(s1+1)/proposal_ratio
  }
  if(change=='remove'){
    proposal_ratio = s1 / (T-s1+1) * proposal_ratio
  }
  final_ratio = L_ratio+G_ratio+B_ratio+proposal_ratio
  return(c(L_ratio, G_ratio, B_ratio, proposal_ratio))
}
