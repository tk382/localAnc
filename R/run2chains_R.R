run2chains_R = function(X,Y,
                        initial_chain1,
                        initial_chain2,
                        Phi,
                        marcor,
                        niter=1000,
                        bgiter=500,
                        hiter=50,
                        switer=50,
                        burnin=5){
  T = ncol(Y)
  n = nrow(Y)
  nu = T+4
  Vbeta     = mean(marcor^2) * 0.0001
  outbeta1  = outbeta2  = matrix(0,niter,T)
  outgam1   = outgam2   = matrix(0,niter,T)
  outSigma1 = outSigma2 = array(0,dim=c(T,T,niter))
  outsb1    = outsb2    = rep(0,niter)
  outh1     = outh2     = rep(0,niter)
  tar1      = tar2      = matrix(0,niter,4)

  outbeta1[1,]         = initial_chain1$beta
  outgam1[1,]          = initial_chain1$gamma
  outSigma1[,,1]       = initial_chain1$Sigma
  outsb1[1]            = initial_chain1$sigmabeta

  outbeta2[1,]         = initial_chain2$beta
  outgam2[1,]          = initial_chain2$gamma
  outSigma2[,,1]       = initial_chain2$Sigma
  outsb2[1]             = initial_chain2$sigmabeta

  for (i in 2:niter){
    ##chain 1 update
    bg = update_betagam_sw(X,
                        Y,
                        outgam1[i-1,],
                        outbeta1[i-1,],
                        outSigma1[,,i-1],
                        abs(marcor),
                        outsb1[i-1],
                        Vbeta,
                        bgiter,
                        T)
    # update betagam arguments
    # X,
    # Y,
    # gam1,
    # beta1,
    # Sigma,
    # marcor,
    # sigmabeta,
    # Vbeta,
    # bgiter,
    # T

    outgam1[i,] = bg$gam
    print(outgam1[i,])
    outbeta1[i,] = bg$beta
    outSigma1[,,i] = update_Sigma(n,nu,X,outbeta1[i,],Phi,Y)
    hsig = update_h(outh1[i-1],
                    hiter,
                    outgam1[i,],
                    outbeta1[i,],
                    outSigma1[,,i],
                    X,
                    T)
    outh1[i] = hsig$h
    outsb1[i] = hsig$sigbeta
    if(is.nan(outsb1[i])){outsb1[i]=1000}
    temptar = get_target(X, Y,
                         outsb1[i],
                         outSigma1[,,i],
                         outgam1[i,],
                         outbeta1[i,],
                         T)

    tar1[i,1] = temptar$L
    tar1[i,2] = temptar$B
    tar1[i,3] = temptar$G
    tar1[i,4] = temptar$final
    ##chain 2 update
    bg = update_betagam(X,
                        Y,
                        outgam2[i-1,],
                        outbeta2[i-1,],
                        outSigma2[,,i-1],
                        abs(marcor),
                        outsb2[i-1],
                        Vbeta,
                        bgiter,
                        T)
    outgam2[i,] = bg$gam
    #print(outgam2[i,])
    outbeta2[i,] = bg$beta
    outSigma2[,,i] = update_Sigma(n,nu,X,outbeta2[i,],Phi,Y)
    hsig = update_h(outh2[i-1],
                    hiter,
                    outgam2[i,],
                    outbeta2[i,],
                    outSigma2[,,i],
                    X,
                    T)
    outh2[i] = hsig$h
    outsb2[i] = hsig$sigbeta
    if(is.nan(outsb2[i])){outsb2[i]=1000}
    temptar = get_target(X,
                        Y,
                        outsb2[i],
                        outSigma2[,,i],
                        outgam2[i,],
                        outbeta2[i,],
                        T)
    tar2[i,1] = temptar$L
    tar2[i,2] = temptar$B
    tar2[i,3] = temptar$G
    tar2[i,4] = temptar$final
    if(i > 2*burnin && i%%10==0){
      est1 = as.numeric(colMeans(outgam1[burnin:i, ])>0.5)
      est2 = as.numeric(colMeans(outgam2[burnin:i, ])>0.5)
      if(sum(est1==est2)==T){
        tmpbeta1 = outbeta1[burnin:i, ]
        tmpbeta2 = outbeta2[burnin:i, ]
        tmpgam1  = outgam1[burnin:i, ]
        tmpgam2  = outgam2[burnin:i, ]
        tmpbeta1[tmpgam1==0] = NA
        tmpbeta2[tmpgam2==0] = NA
        betaest1 = colMeans(tmpbeta1, na.rm=TRUE)
        betaest2 = colMeans(tmpbeta2, na.rm=TRUE)
        if(sum((betaest1[est1] - betaest2[est2])^2)/length(est1) < 1e-3){
          print('beta difference is small between two chains - converged!')
          outgam1 = outgam1[1:i,];         outgam2   = outgam2[1:i,]
          outbeta1 = outbeta1[1:i,];       outbeta2  = outbeta2[1:i,]
          outSigma1 = outSigma1[,,1:i];    outSigma2 = outSigma2[,,1:i]
          outh1 = outh1[1:i];              outh2     = outh2[1:i]
          outsb1 = outsb1[1:i];            outsb2    = outsb2[1:i]
          break;
        }
      }
    }
    print(i)
  }
  return(list(
      chain1 = list(gamma = outgam1, beta = outbeta1, Sigma = outSigma1,
                    h = outh1, sigmabeta = outsb1, tar = tar1),
      chain2 = list(gamma = outgam2, beta = outbeta2, Sigma = outSigma2,
                    h = outh2, sigmabeta = outsb2, tar = tar2)
    ))

}
