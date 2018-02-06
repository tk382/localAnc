update_betagam_sw = function(X,
                          Y,
                          gam1,
                          beta1,
                          Sigma,
                          marcor,
                          sigmabeta,
                          Vbeta,
                          bgiter,
                          T){
  for (i in 2:bgiter){
    #small world proposal
    if(i%%100==0){ #small world proposal
      temp = update_gamma_smallworld(X, Y, gam1, T)
      gam2 = as.numeric(temp$newgamma)
      changeind = temp$changeind
      beta2 = beta1*gam2
      ind = which(gam2==1)
      beta2[ind] = beta1[ind]+rnorm(length(ind), 0, sqrt(Vbeta))
      changeind = temp$changeind
      change = gam2[changeind]
      proposal_ratio = sum(dnorm(beta1[changeind]-beta2[changeind],
                           0, sqrt(Vbeta), log = TRUE))
      newtarget = get_target(X, Y, sigmabeta, Sigma, gam2, beta2,T)$final
      oldtarget = get_target(X, Y, sigmabeta, Sigma, gam1, beta1,T)$final
      A = newtarget - oldtarget + proposal_ratio
      check = runif(1)
      if(exp(A)>check){
        gam1 = gam2; beta1 = beta2;
      }
    }else{
      temp = update_gamma(X, Y, gam1, marcor)
      gam2 = as.numeric(temp$newgamma)
      beta2 = beta1*gam2
      ind = which(gam2==1)
      beta2[ind] = beta1[ind] + rnorm(length(ind), 0, sqrt(Vbeta))
      changeind = temp$changeind
      change = gam2[changeind]
      A = betagam_accept(X,
                         Y,
                         sigmabeta,
                         Sigma,
                         Vbeta,
                         marcor,
                         gam1, beta1,
                         gam2, beta2,
                         changeind,
                         change,
                         T)
      check = runif(1,0,1)
      print(gam1); print(gam2)
      if(exp(A[1])>check){gam1 = gam2; beta1 = beta2}
    }
  }
  return(list(gam = gam1, beta = beta1))
}


update_gamma_smallworld = function(X,
                        Y,
                        gamma,
                        T){
  where = sample(1:T, size=round(T/2))
  newgamma = gamma
  newgamma[where] = -newgamma[where]+1
  return(list(newgamma=newgamma, changeind=where))
}
