update_betagam = function(X,
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
    temp = update_gamma(X, Y, gam1, marcor)
    gam2 = as.numeric(temp$newgamma);
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
    print(which(gam1!=gam2))
    print(gam2[which(gam1!=gam2)])
    print(A)
    check = runif(1,0,1)
    if(exp(A[1]) > check){
      print("make update!")
      #gam1 = gam2; beta1 = beta2;
    }
    i = i+1
  }
  return(list(gam = gam2, beta = beta2))
}

update_gamma = function(X,
                        Y,
                        gamma,
                        marcor){
  newgamma = gamma
  T = length(gamma)
  p1 = p2 = 0.5;
  ind0 = which(gamma==0)
  ind1 = which(gamma==1)
  s = length(ind1)
  case = sample(c(1,2), size=1, prob = c(p1,p2))
  if(s==0){
    ##nothing to remove
    case = 1
  }else if(s==T){
    #nothing to add
    case = 2
  }
  if (case==1){
    add = ind0[1]
    if(s<(T-1)){
      add = sample(ind0, size=1, prob=marcor[ind0])
    }
    newgamma[add] = 1
    changeind = add
  }
  if (case==2){
    remove = sample(ind1, size=1)
    newgamma[remove] = 0
    changeind = remove
  }
  return(list(newgamma=newgamma, changeind=changeind))
}

betagam_accept = function(X,
                          Y,
                          sigmabeta1,
                          inputSigma,
                          Vbeta,
                          marcor,
                          gam1,
                          beta1,
                          gam2,
                          beta2,
                          changeind,
                          change,
                          T){
  newtarget = sum(get_target_new(X, Y, sigmabeta1, inputSigma, gam2, beta2, T))
  oldtarget = sum(get_target_new(X, Y, sigmabeta1, inputSigma, gam1, beta1, T))
  proposal_ratio = dnorm(beta1[changeind]-beta2[changeind],
                         0, sqrt(Vbeta), log = TRUE)
  s1 = sum(gam1==1)
  s2 = sum(gam2==1)
  if(change==1){
    tempadd = marcor[changeind] / sum(marcor[gam1==0])
    proposal_ratio = -log(tempadd)-log(s2)-proposal_ratio
  }
  if(change==0){
    tempadd = marcor[changeind] / sum(marcor[gam2==0])
    proposal_ratio = log(tempadd) + log(s1) + proposal_ratio
  }
  final_ratio = newtarget - oldtarget + proposal_ratio
  return(c(final_ratio, newtarget, oldtarget, proposal_ratio))
}
