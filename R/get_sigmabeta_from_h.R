get_sigmabeta_from_h = function(h,
                                gam,
                                Sigma,
                                X,
                                T){
  n = length(X)
  ds = diag(Sigma)
  num = h * sum(ds)
  denom = (1-h)*sum(ds[gam==1]) * sum(X^2)/n
  return(num/denom)
}
