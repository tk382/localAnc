update_h=function(initialh, hiter,gam, beta, Sig,X,T){
  h1 = initialh
  sigbeta1 = get_sigmabeta_from_h(h1, gam, Sig, X, T)
  lik = rep(0,hiter)
  if(all(beta==0)){return(list(h=0,sigbeta = 0))}
  for (i in 2:hiter){
    h2 = h1
    r = runif(1, -0.1, 0.1)
    h2 = h2+r
    if (h2 < 0){h2 = (abs(h2))}
    if (h2 > 1){h2 = 2-h2}
    ind = which(gam==1)
    sigbeta1 = get_sigmabeta_from_h_c(h1,gam,Sig,X,T)
    sigbeta2 = get_sigmabeta_from_h_c(h2,gam,Sig,X,T)
    lik1 = sum(dnorm(beta[ind], rep(0,length(ind)),
                     sqrt(sigbeta1*diag(Sig)[ind]), log=TRUE))
    lik2 = sum(dnorm(beta[ind], rep(0,length(ind)),
                     sqrt(sigbeta2*diag(Sig)[ind]), log=TRUE))
    acceptanceprob = exp(lik2-lik1)
    e = runif(1,0,1)
    if(e < acceptanceprob){
      h1 = h2
      sigbeta1 = sigbeta2
    }
  }
  return(list(h=h2, sigbeta = sigbeta2))
}
