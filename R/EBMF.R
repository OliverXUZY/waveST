## ===========================================================================================
## EBNM
## ===========================================================================================
### 1 normal prior (simple)
logLik = function(par, x, s){
  V = exp(par[2])
  -sum(dnorm(x, mean = par[1], sd = sqrt(s^2 + V), log = TRUE))
}

#logLik(c(0,1), rnorm(10), rep(1,10))

ebnm_normalPrior = function(x, s, init = c(0,0)){
  res = optim(par = init, fn = logLik, method = "BFGS", x = x, s = s)
  return(res$par)
}

# mu, V = cat(res[1], exp(res[2]))

#x = rnorm(10); s = rep(1,10)
#ebnm_normalPrior(x,s)
poster = function(x,s, par){
  pMu = (x/s^2 + par[1]/exp(par[2])) / (1/s^2 + 1/exp(par[2]))
  pSigma2 = 1 / (1/s^2 + 1/exp(par[2]))
  return(list(pMu = pMu, pSigma2 = pSigma2))
}

ebnm = function(x, s, init = c(.1,.1)){
  par = ebnm_normalPrior(x, s, init = init)
  posterior = poster(x,s,par)
  return(posterior)
}

## ===========================================================================================
## EBMF
## ===========================================================================================
# ----- simple one factor- rank1
# ----- Algorithm 2

EBMF_rank1 = function(Y, tol = 1e-3, maxIter = 1000){
  n = dim(Y)[1]; p = dim(Y)[2]
  
  ## init with SVD
  res = svd(Y); d = res$d[1]
  LL = matrix(res$u[,1],nrow = n) * sqrt(d)
  FF = matrix(res$v[,1],nrow=p) * sqrt(d)
  LL2 = LL^2
  FF2 = FF^2
  yhat = list(LL %*% t(FF))
  LLlist = list(LL)
  FFlist = list(FF)
  
  for(iter in 1:maxIter){
    ## update tau
    R2 = (Y - LL %*% t(FF))^2 + LL2 %*% t(FF2) - LL^2 %*% t(FF^2)
    tau = n*p/sum(R2)
    
    ## update ll
    s2_hat = 1/sum(tau*FF2)
    ll_hat = tau*Y%*%FF/sum(tau*FF2)
    ll_res = ebnm(ll_hat, sqrt(s2_hat))
    LL = ll_res$pMu; LL2 = matrix(ll_res$pSigma2 + ll_res$pMu^2,nrow = n)
    
    ## update ff
    s2_hat2 = 1/sum(tau*LL2)
    ff_hat = t(tau*t(LL)%*%Y/sum(tau*LL2))
    ff_res = ebnm(ff_hat, sqrt(s2_hat2))
    FF = ff_res$pMu; FF2 = matrix(ff_res$pSigma2 + ff_res$pMu^2,nrow=p)
    
    LLlist[[iter + 1]] = LL
    FFlist[[iter + 1]] = FF
    yhat[[iter + 1]] = LL %*% t(FF)
    
    if(sum((yhat[[iter+1]] - yhat[[iter]])^2) < tol){
    #if(sum((LLlist[[iter+1]] - LLlist[[iter]])^2)/n/p < tol){
      message("tolerance satisfied || rank 1 fitting || ", iter, "-th iteration")
      break
    }
  }
  
  return(list(LL = LL, FF = FF, LL2 = LL2, FF2 = FF2, tau = tau))
 
  
  
}
#res = EBMF_rank1(Y)

#sum((res$LL%*% t(res$FF) - LL %*% t(FF))^2)

## ===========================================================================================
# ----- Single-factor update for EBMF (rank K)
# ----- Algorithm 3
EBMF_sigFct = function(Y, LL, FF, LL2, FF2, k){ # update k-th factor
  n = dim(Y)[1]; p = dim(Y)[2]; K = dim(LL)[2]
  
  ## update tau
  R2 = (Y - LL %*% t(FF))^2 + LL2 %*% t(FF2) - LL^2 %*% t(FF^2)
  tau = n*p/sum(R2)
  
  ## Residual matrix
  Rk = Y - LL[, -k, drop = FALSE] %*% t(FF[, -k, drop = FALSE])
  
  ## update llk -- single factor, only involve FF[,k] and FF2[,k], computation same as same as rank 1
  s2_hat = 1/sum(tau*FF2[,k])
  ll_hat = tau*Rk%*%FF[,k]/sum(tau*FF2[,k])
  ll_res = ebnm(ll_hat, sqrt(s2_hat))
  lk = ll_res$pMu; lk2 = matrix(ll_res$pSigma2 + ll_res$pMu^2,nrow = n)
  
  ## update ffk
  s2_hat2 = 1/sum(tau*lk2)
  ff_hat = t(tau*t(lk)%*%Rk/sum(tau*lk2))
  ff_res = ebnm(ff_hat, sqrt(s2_hat2))
  fk = ff_res$pMu; fk2 = matrix(ff_res$pSigma2 + ff_res$pMu^2,nrow=p)
  
  return(list(lk = lk, fk = fk, lk2 = lk2, fk2 = fk2, tau = tau))
}
# res = EBMF_sigFct(Y, LL, FF, LL2, FF2, 1)


## ===========================================================================================
# ----- Greedy Algorithm for EBMF
# ----- Algorithm 4
EBMF_greedy = function(Y, tol = 1e-3, maxIter = 1000, maxK){
  n = dim(Y)[1]; p = dim(Y)[2]
  # K is estimated rank, first we do rank-1 EBMF
  # First we do rank-1 EBMF
  K = 1 
  res = EBMF_rank1(Y, tol = 1e-3)
  LL = res$LL; FF = res$FF; LL2 = res$LL2; FF2 = res$FF2; tau = res$tau
  
  #while (iter <= maxIter) {  ###  problem!
  while(K < maxK) {
    
    K = K + 1
    ## Residual matrix for deciding whether add new factor K
    RK = Y - LL %*% t(FF)
    
    ## init with SVD
    res = svd(RK)
    lk = matrix(res$u[,1], nrow = n)
    fk = matrix(res$v[,1], nrow = p)
    lk2 = lk^2
    fk2 = fk^2
    LL = cbind(LL, lk); FF = cbind(FF, fk)
    LL2 = cbind(LL2, lk2); FF2 = cbind(FF2, fk2)
    
    iter = 1
    lkfk = list(lk %*% t(fk))
    while (TRUE) {
      iter = iter + 1
      
      res = EBMF_sigFct(Y, LL, FF, LL2, FF2, K)
      LL[,K] = res$lk; FF[,K] = res$fk
      LL2[,K] = res$lk2; FF2[,K] = res$fk2
      
      lkfk[[iter]] =res$lk %*% t(res$fk)
      
      if(sum((lkfk[[length(lkfk)]] - lkfk[[length(lkfk)-1]])^2) < tol){
        message("tolerance satisfied || for ",K,"-th factor single update")
        break
      }
    }
    if(sum(abs(lk)) <= 1e-4 | sum(abs(fk)) <= 1e-4) {
      message("encounter lk or fk identically 0")
      break
    }
  }
  return(list(LL = LL, FF = FF, LL2 = LL2, FF2 = FF2, tau = tau))
}

# res = EBMF_greedy(Y,maxK = 20)

## ===========================================================================================
# ----- Backfitting algorithm for EBMF (rank K)
# ----- Algorithm 5
EBMF_bkft = function(Y, K, initMethod = c("greedy", "svd"), tol = 1e-3, maxIter = 1000){
  n = dim(Y)[1]; p = dim(Y)[2]
  
  ## init
  if(initMethod == "svd"){
    ## init with SVD
    res = svd(Y); d = res$d 
    if(K < min(n,p)) d[(K + 1):min(n,p)] <- 0
    D = sqrt(diag(d))
    LL = (matrix(res$u,nrow = n) %*% D)[, 1:K, drop = FALSE]
    FF = (matrix(res$v,nrow=p) %*% D)[, 1:K, drop = FALSE]
    LL2 = LL^2
    FF2 = FF^2
    message("SVD init complete")
  } else if (initMethod == "greedy"){
    res = EBMF_greedy(Y,maxK = K)
    LL = res$LL
    FF = res$FF
    LL2 = res$LL2
    FF2 = res$FF2
    message("greedy init complete")
  } else {
    stop("init method must be SVD or greedy")
  }
  
  yhat = list(LL %*% t(FF))
  
  for(iter in 1:maxIter) {
    for(k in 1:K) {
      res = EBMF_sigFct(Y, LL, FF, LL2, FF2, k)
      LL[,k] = res$lk; FF[,k] = res$fk
      LL2[,k] = res$lk2; FF2[,k] = res$fk2
      tau = res$tau
      message(k,"-th factor updated || ",iter, "-th iteration")
    }
    
    yhat[[iter + 1]] = LL %*% t(FF)
    
    if((sum(yhat[[iter+1]] - yhat[[iter]])^2) < tol){
      message("tolerance satisfied || ",iter, "-th iteration")
      break
    }
  }
  return(list(LL = LL, FF = FF, LL2 = LL2, FF2 = FF2, tau = tau))
}

#K = 20; initMethod = "svd"; tol = 1e-3; maxIter = 10
# EBMF_bkft(Y, K = 20, initMethod = "svd", tol = 1e-3, maxIter = 1000)






















  


