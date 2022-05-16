### helper functions
threshold_coef = function(gridValue, hard = c(TRUE,FALSE),  ...){
  genedwt = dwt.2d(gridValue, ...)
  
  ## set up layout
  lay = layout(matrix(1:4,2,2,byrow = TRUE))
  layout.show(lay)
  
  ## original
  image.plot(idwt.2d(genedwt),asp = 1)
  
  ## SURE threshold
  genedwt_sure = sure.thresh(genedwt, hard = hard)
  image.plot(idwt.2d(genedwt_sure),asp = 1)
  
  ## universal threshold with value sqrt(2logn)
  genedwt_universal = universal.thresh.dwt.2d(genedwt, hard = hard)
  image.plot(idwt.2d(genedwt_universal),asp = 1)
  
  ## hybrid threhold: SURE threshold + universal threshold with value sqrt(2logn)
  genedwt_hybrid = hybrid.thresh(genedwt)
  image.plot(idwt.2d(genedwt_hybrid),asp = 1)
  
  layout(1)
}




## This is an encapsulated function for decomposition alone
## @ raws raw data matrix with dim is D^2 by Number of genes
## @ decom_method do we decompose on raw data or wavelet?
## @ K designated number of factors
## @ threshold threshold for selecting K number of factors, not using if K is specified
## @ wf wavelet filter
## @ J wavelet level
## @ thresholdMethod: method for thresholding, c("hybrid", "manual")
## @ tau: wavelet threshold in manual thresholding, only used for manual thresholding
decomp = function(raws, decom_method = c("raw", "wave"),  K = NULL, threshold = 500, replacement = TRUE, 
                  wf = "d4", J = 5, thresholdMethod = "manual", tau = 0){
  
  ##########=================== ##########=================== ##########=================== ##########=================== 
  ##########=================== raw decomposition
  if(decom_method == "raw") {
    ##########===================  SVD
    s = svd(raws)
    if(is_null(K)){
      K = sum(s$d > threshold)
    }
    
    u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
    
    recon_SVD = u%*%D%*%t(v)
    
    SVD = list(u = u,D = D,  v = v, recon_SVD = recon_SVD)
    
    ##########==================== EBMF
    res = EBMF_bkft(as.matrix(raws), K = K, initMethod = "greedy", tol = 1e-10, maxIter = 30)
    
    recon_EBMF = res$LL%*% t(res$FF)
    
    EBMF = list(LL = res$LL, FF = res$FF, recon_EBMF = recon_EBMF)
    
    return(list(SVD = SVD, EBMF = EBMF))
  }
  
  ##########=================== ##########=================== ##########=================== ##########=================== 
  ##########=================== wavelet decomposition
  if(decom_method == "wave") {
    ## apply WaveTransCoefs and threshold to each column of raws x
    coefs_matrix = 1:dim(raws)[2] %>%
      map_dfc(~ WaveTransCoefs(raws[,.x], wf = wf, J = J, thresholdMethod = thresholdMethod, tau = tau))
    
    ##########===================  SVD
    s = svd(coefs_matrix)
    if(is_null(K)){
      K = sum(s$d > threshold)
    }
    u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
    
    SVD = u%*%D%*%t(v)
    recon_SVD = split(SVD, col(SVD)) %>%
      map(~ InvWaveTrans(.x, raws = raws)) %>%
      map_dfc(~as.vector(.)) %>% 
      as.matrix()
    
    SVD = list(u = u,D = D,  v = v, recon_SVD = recon_SVD)
    
    ##########==================== EBMF
    
    res = EBMF_bkft(as.matrix(coefs_matrix), K = K, initMethod = "greedy", tol = 1e-10, maxIter = 30)
    
    EBMF = res$LL%*% t(res$FF)
    recon_EBMF = split(EBMF, col(EBMF)) %>%
      map(~ InvWaveTrans(.x, raws = raws)) %>%
      map_dfc(~as.vector(.)) %>% 
      as.matrix()
    
    EBMF = list(LL = res$LL, FF = res$FF, recon_EBMF = recon_EBMF)
    
    return(list(SVD = SVD, EBMF = EBMF))
  }
}


## This is an encapsulated function for decomposition alone
## @ raws raw data matrix with dim is D^2 by Number of genes
## @ decom_method do we decompose on raw data or wavelet?
## @ K designated number of factors
## @ threshold threshold for selecting K number of factors, not using if K is specified
## @ wf wavelet filter
## @ J wavelet level
## @ thresholdMethod: method for thresholding, c("hybrid", "manual")
## @ tau: wavelet threshold in manual thresholding, only used for manual thresholding
decomp2 = function(raws, wavemethod = c("raw", "wave"), decom_method = c("SVD", "EBMF"),  K = NULL, threshold = 500,
                  wf = "d4", J = 5, thresholdMethod = "manual", tau = 0){
  
  raws = raws %>% as.matrix()
  ##########=================== ##########=================== ##########=================== ##########=================== 
  ##########=================== raw decomposition
  if(wavemethod == "raw") {
    s = svd(raws)
    if(is_null(K)){
      K = sum(s$d > threshold)
    }
    if(decom_method == "SVD"){
      ##########===================  SVD
      
      u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
      
      recon_SVD = u%*%D%*%t(v)
      
      SVD = list(f = u,D = D,  l = v, recon = recon_SVD)
      
      return(SVD)
    }
    
    if(decom_method == "EBMF"){
      ##########==================== EBMF
      res = EBMF_bkft(as.matrix(raws), K = K, initMethod = "greedy", tol = 1e-10, maxIter = 30)
      
      recon_EBMF = res$LL%*% t(res$FF)
      
      EBMF = list(f = res$FF, D = diag(rep(1,K)), l = res$LL, recon = recon_EBMF)
      
      return(EBMF)
    }
    
  }
  
  ##########=================== ##########=================== ##########=================== ##########=================== 
  if(wavemethod == "wave") {
    ##########=================== wavelet decomposition
    ## apply WaveTransCoefs and threshold to each column of raws x
    coefs_matrix = 1:dim(raws)[2] %>%
      map_dfc(~ WaveTransCoefs(raws[,.x], wf = wf, J = J, thresholdMethod = thresholdMethod, tau = tau))
    
    s = svd(coefs_matrix)
    if(is_null(K)){
      K = sum(s$d > threshold)
    }
    if(decom_method == "SVD") {
      ##########===================  SVD
      
      u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
      
      SVD = u%*%D%*%t(v)
      recon_SVD = split(SVD, col(SVD)) %>%
        map(~ InvWaveTrans(.x, raws = raws)) %>%
        map_dfc(~as.vector(.)) %>% 
        as.matrix()
      
      SVD = list(f = u,D = D,  l = v, recon = recon_SVD)
      return(SVD)
    }
    
    if(decom_method == "EBMF") {
      ##########==================== EBMF
      
      res = EBMF_bkft(as.matrix(coefs_matrix), K = K, initMethod = "greedy", tol = 1e-10, maxIter = 30)
      
      EBMF = res$LL%*% t(res$FF)
      recon_EBMF = split(EBMF, col(EBMF)) %>%
        map(~ InvWaveTrans(.x, raws = raws)) %>%
        map_dfc(~as.vector(.)) %>% 
        as.matrix()
      
      EBMF = list(f = res$FF, D = diag(rep(1,K)), l = res$LL, recon = recon_EBMF)
      
      return(EBMF)
    }
  }
}

