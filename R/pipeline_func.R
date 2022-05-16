
## function for generate patterns
#'Generate patterns in simulations
#'
#'Generate 9 matrix with patterns among spatial locations, this is factor genes in simulation.
#'@docType methods
#'@export
#'@param N_dup how many duplicates to generate
#'@return a list contaning nine D^2 matrix
generateCoorM = function(N_dup){
  ls = list()
  if(N_dup == 1){
    al = 0; bl = 0; cl = 0; dl = 0
  }
  al = seq(0,100,length.out = N_dup)
  bl = seq(-50,80,length.out = N_dup)
  cl = seq(-40,60,length.out = N_dup)
  dl = seq(-100,20,length.out = N_dup)
  
  for(i in 1:N_dup){
    a = al[i]; b = bl[i]; c = cl[i]; d = dl[i]; 
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[row(m) - col(m) > 0 & row(m) - col(m) <= 8] = 40 + a
    m[row(m) - col(m) > 8 & row(m) - col(m) <= 16] = 50 + b
    m[row(m) - col(m) > 16 & row(m) - col(m) <= 24] = 60 + c
    m[row(m) - col(m) > 24 & row(m) - col(m) <= 32] = 70 + d
    m[row(m) - col(m) <= 0 & row(m) - col(m) >= -8] = 30- a
    m[row(m) - col(m) < -8 & row(m) - col(m) >= -16] = 20 - b
    m[row(m) - col(m) < -16 & row(m) - col(m) >= -24] = 10 - c
    m[row(m) - col(m) < -24 & row(m) - col(m) >= -32] = 0 - d
    ls[[length(ls)+1]] = m
    #image.plot(m)
    
    ##
    # m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    # m[row(m) - col(m) > 0 & row(m) - col(m) <= 12] = 40 + a
    # m[row(m) - col(m) >12 & row(m) - col(m) <= 24] = 50 + b
    # m[row(m) - col(m) > 24 & row(m) - col(m) <= 32] = 60 + c
    # m[row(m) - col(m) <= 0 & row(m) - col(m) >= -12] = 30 - a
    # m[row(m) - col(m) < -12 & row(m) - col(m) >= -24] = 20 - b
    # m[row(m) - col(m) < -24 & row(m) - col(m) >= -32] = 10 - c
    # ls[[length(ls)+1]] = m
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[row(m) + col(m) > 0 & row(m) + col(m) <= 8] = 40  + a
    m[row(m) + col(m) > 8 & row(m) + col(m) <= 16] = 50 + b
    m[row(m) + col(m) > 16 & row(m) + col(m) <= 24] = 60 + c
    m[row(m) + col(m) > 24 & row(m) + col(m) <= 32] = 70 + d
    m[row(m) + col(m) > 32 & row(m) + col(m) <= 40] = 30 - a
    m[row(m) + col(m) > 40 & row(m) + col(m) <= -48] = 20 - b
    m[row(m) + col(m) > 48 & row(m) + col(m) <= 56] = 10 - c
    m[row(m) + col(m) > 56 & row(m) + col(m) <= 64] = 0 - d
    #image.plot(m)
    ls[[length(ls)+1]] = m
    
    # ##
    # m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    # m[row(m) + col(m) > 0 & row(m) + col(m) <= 12] = 40 + a
    # m[row(m) + col(m) >12 & row(m) + col(m) <= 24] = 50 + b
    # m[row(m) + col(m) > 24 & row(m) + col(m) <= 36] = 60 + c
    # m[row(m) + col(m) > 36 & row(m) + col(m) <= 48] = 30 - a
    # m[row(m) + col(m) > 48 & row(m) + col(m) <= 64] = 20 - b
    # #image.plot(m)
    # ls[[length(ls)+1]] = m
    
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[row(m)^2 + col(m)^2 >=1000 | row(m)^2 + col(m)^2 <= 100] = 10 + a
    m[row(m)^2 + col(m)^2 >=500 & row(m)^2 + col(m)^2 <= 1000] = 20 - a
    ls[[length(ls)+1]] = m
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[(row(m) - 32)^2 + col(m)^2 >=1000 | (row(m)- 32)^2 + col(m)^2 <= 100] = 10 
    m[(row(m) - 32)^2 + col(m)^2 >=500 & (row(m) - 32)^2 + col(m)^2 <= 1000] = 20 
    ls[[length(ls)+1]] = m
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[(row(m) - 32)^2 + (col(m) - 32)^2 >=1000 | (row(m)- 32)^2 + (col(m) - 32)^2 <= 100] = 10 
    m[(row(m) - 32)^2 + (col(m) - 32)^2 >=500 & (row(m) - 32)^2 + (col(m) - 32)^2 <= 1000] = 20 
    ls[[length(ls)+1]] = m
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[row(m)^2 + (col(m) - 32)^2 >=1000 | row(m)^2 + (col(m) - 32)^2 <= 100] = 10 
    m[row(m)^2 + (col(m) - 32)^2 >=500 & row(m)^2 + (col(m) - 32)^2 <= 1000] = 20 
    ls[[length(ls)+1]] = m
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[row(m)^2 - col(m)^2 >=1000 | row(m)^2 - col(m)^2 <= -500] = 10 + a
    m[row(m)^2 - col(m)^2 >=500 & row(m)^2 - col(m)^2 <= 1000] = 20 + b
    m[row(m)^2 - col(m)^2 >=10 & row(m)^2 - col(m)^2 <= 500] = 30 - a
    ls[[length(ls)+1]] = m
    #image.plot(m)
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[(row(m)-15)^2 + (col(m)-15)^2 >=10 & (row(m)-15)^2 + (col(m)-15)^2 <= 50] = 10 + a
    m[(row(m)-15)^2 + (col(m)-15)^2 >=50 & (row(m)-15)^2 + (col(m)-15)^2 <= 150] = 20 - a
    m[(row(m)-15)^2 + (col(m)-15)^2 >=150 & (row(m)-15)^2 + (col(m)-15)^2 <= 300] = 30 + b
    #image.plot(m,asp = 1)
    ls[[length(ls)+1]] = m
    
    ##
    m = matrix(0, nrow = 32, ncol = 32, byrow = TRUE); m
    m[3*(row(m)-14)^2 + 2*(col(m)-16)^2 >=10 & 3*(row(m)-14)^2 + 2*(col(m)-16)^2 <= 50] = 10 + b
    m[3*(row(m)-14)^2 + 2*(col(m)-16)^2 >=100 & 3*(row(m)-14)^2 + 2*(col(m)-16)^2 <= 300] = 20 + c
    m[3*(row(m)-14)^2 + 2*(col(m)-16)^2 >=300 & 3*(row(m)-14)^2 + 2*(col(m)-16)^2 <= 500] = 30 + d
    ls[[length(ls)+1]] = m
    #image.plot(m,asp = 1)
  }
  return(ls)
}



#### function for get n by p data

#' k over A method for spatially resolved transcriptomics dataset
#'
#'Collect spatially resolved transcriptomics datasets in SpatialExperiment Bioconductor format using Visium_humanDLPFC_3_13,
#'then use kOverA select genes

#'@docType methods
#'@export
#'@param k The number of elements that have to exceed A.
#'@param A The value you want to exceed.
#'@examples 
#'res = kOverA_ST(k = 3,A = 7)
### This function get n*p matrix (p genes after k over A selected) and coordinates
kOverA_ST = function(k = 5, A = 50){
#pipe1_STE = function(k = 5, A = 50){
  spe <- Visium_humanDLPFC()
  
  coor = spatialCoords(spe)
  raw = assay(spe)
  viz = as_tibble(coor)
  
  ## gene filter
  f1 = kOverA(k, A)
  ffun_combined = filterfun(f1)
  wh3 = genefilter(raw, ffun_combined)
  df = t(raw[wh3,] %>% as.matrix())
  #df = as.matrix(df)
  
  return(list(viz = viz, df = df))
}


############==============#####==============#####==============#####==============#####==============
##  refine coordinates into a grid 64*64
## This function set up group membership by range of interval in gridX
partition = function(boundary, coor1d, grid) {
  return(rep(x = boundary, times = length(coor1d[coor1d >= grid[boundary] & coor1d < grid[boundary + 1]])))
}

## set up group membership, this viz requires gene column
Gmbsp = function(viz, level = 6) {
  viz = viz[order(viz$x),]
  gridX = seq(from = range(viz$x)[1], to = (range(viz$x)[2] + 1), length.out = 2^level + 1)
  partX = mapply(partition,boundary = seq(1,2^level), MoreArgs = list(coor1d = viz$x, grid = gridX))
  
  missing = sum(unlist(lapply(partX, function(x) length(x))) == 0) # can use lengths() directly
  if( missing > 0) {
    message("there are ",missing, " groups in x-axis are empty" )
  }
  
  viz$xg = unlist(partX)
  
  viz = viz[order(viz$y),]
  gridY = seq(from = range(viz$y)[1], to = (range(viz$y)[2] + 1), length.out = 2^level + 1)
  partY = mapply(partition,boundary = seq(1,2^level), MoreArgs = list(coor1d = viz$y, grid = gridY))
  
  missing = sum(unlist(lapply(partY, function(x) length(x))) == 0)  # can use lengths() directly
  if( missing > 0) {
    message("there are ",missing, " groups in y-axis are empty" )
  }
  
  viz$yg = unlist(partY)
  
  
  return(viz)
}

coorToGrid = function(viz, level = 6){
  viz2 = Gmbsp(viz, level = level)
  viz2$x = NULL; viz2$y = NULL
  
  gridValue = aggregate(viz2, by = list(viz2$xg, viz2$yg), FUN = "mean")
  gridValue = dcast(gridValue, xg ~ yg, value.var = "gene")
  gridValue$xg = NULL
  
  gridValue[is.na(gridValue)] = 0
  
  return(as.matrix(gridValue))
}


##  refine coordinates into a grid 64*64
# res = kOverA_ST()
# viz = res$viz; df = res$df
# viz$gene = df[,1]
# gridValue = coorToGrid(viz)

### This function is the same as coorToGrid, but take gene as extra argument, 
### it takes one more step is viz$gene = gene
### this is convenient for pipeline
coorToGridGene = function(viz, gene, level = 6){
  viz$gene = gene
  viz2 = Gmbsp(viz, level = level)
  viz2$x = NULL; viz2$y = NULL
  
  gridValue = aggregate(viz2, by = list(viz2$xg, viz2$yg), FUN = "mean")
  gridValue = dcast(gridValue, xg ~ yg, value.var = "gene")
  gridValue$xg = NULL
  
  gridValue[is.na(gridValue)] = 0
  
  return(as.matrix(gridValue))
}

##  refine coordinates into a grid 64*64
# res = kOverA_ST()
# viz = res$viz; df = res$df
# gridValue = coorToGridGene(viz, df[,1])

### This function take original n*p matrix to grid by gene matrix D^2*p, apply coorToGrid to each gene of a matrix
### viz: coor matrix, n*2
### df: raw matrix n*p
gridGeneRaw = function(viz, df, level = 6){
  ## apply coorToGrid to each column of raws
  gridls = 1:dim(df)[2] %>%
    map(~ coorToGridGene(viz, gene = df[,.x], level = level))
  return(gridls)
}

# gridls = gridGeneRaw(viz, df)
# gridtoy = gridGeneRaw(viz, df[,1:5])

# lay = layout(matrix(1:6,2,3,byrow = TRUE))
# gridtoy %>% 
#   map(~.x/sqrt(sum(.x^2))) %>%
#   map(~image.plot(.x, asp = 1))
# layout(1)


# threshold_coef(viz, geneID = 1,  wf = "d4", J = 4, hard = TRUE)

## Different SNR to on simulation pipelines
## @ R 1/SNR = 1/(signal to noise ratio)
## @ truth ground truth matrix with rank K, dim is D^2 by Number of genes
## @ K number of factors
## @ wf wavelet filter
## @ J wavelet level
pipe_SNR = function(R, truth, K, wf = "d4", J = 5) {
  D2 = dim(truth)[1]; NUM_GENES = dim(truth)[2]
  raws = truth + matrix(rnorm(prod(dim(truth)), mean = 0, sd = R*sd(truth)), 
                        nrow = dim(truth)[1], ncol = dim(truth)[2])
  
  ##########==========##########==========##########==========##########==========##########==========
  ##########==========##########==========   SVD on raw
  s = svd(raws)
  u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
  
  ###   convert u to uls
  uls = split(u, col(u)) %>%
    map(~ matrix(.x, nrow = sqrt(length(.x))))
  
  # ##########==========  show all the generated pattern for u in SVD, this u corresponding to the D^2*K
  # png(file=paste("exp_on_cells/","raw_R",R,".jpeg", sep = ""), width = 880, height = 880)
  # lay = layout(matrix(1:9,3,3,byrow = TRUE))
  # #layout.show(lay)
  # #lapply(CMls, image.plot)
  # uls %>%
  #   map(~image.plot(.x, asp = 1))
  # layout(1)
  # dev.off()
  
  ##########================== sum of gradients
  gd = uls %>%
    map(~diff(as.vector(.x))) %>%
    map(~sum(.x^2))
  gd_raw = sum(unlist(gd))
  
  ##########================== recons error
  recon_raw = sum((u%*%D%*%t(v) - truth)^2)/prod(dim(raws))
  
  ##########==========##########==========##########==========##########==========##########==========
  ##########==========##########==========   SVD on wavelet coef
  ####======== This part save partial info for inverse wavelet transformation ####========
  ## gene1 is col1
  col1 = raws[,1]
  
  #WaveTransCoefs(col1)
  m = matrix(col1, nrow = sqrt(length(col1)))
  ## save coefficient list attribute info
  coef1 = dwt.2d(m, wf = wf, J = J)
  
  
  ## save coefficient long vector level info
  coef1M = coef1%>%
    map_dfr(~ tibble(coef = as.numeric(.)), .id = "level") ; coef1M
  level = coef1M$level
  ####======== This part save partial info for inverse wavelet transformation ####========
  
  ## apply WaveTransCoefs and threshold to each column of raws x
  coefs_matrix = 1:dim(raws)[2] %>%
    map_dfc(~ WaveTransCoefs(raws[,.x], wf = wf, J = J, thresholdMethod = "hybrid"))
  
  s = svd(coefs_matrix)
  u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
  
  
  ###   convert u for coef to factor for raws
  uReconsCM_ls = split(u, col(u)) %>%
    map(~ InvWaveTransCoefs(.x, coef1 = coef1, level = level))
  
  # ##########==========  show all the generated pattern for u in SVD, this u corresponding to the D^2*K
  # png(file=paste("exp_on_cells/","coef_R",R,".jpeg", sep = ""), width = 880, height = 880)
  # lay = layout(matrix(1:9,3,3,byrow = TRUE))
  # #layout.show(lay)
  # #lapply(CMls, image.plot)
  # uReconsCM_ls %>%
  #   map(~image.plot(.x, asp = 1))
  # layout(1)
  # dev.off()
  
  ##########================== sum of gradients
  gd = uReconsCM_ls %>%
    map(~diff(as.vector(.x))) %>%
    map(~sum(.x^2))
  gd_coef = sum(unlist(gd))
  
  ##########================== recons error
  SVD = u%*%D%*%t(v)
  SVD = split(SVD, col(SVD)) %>%
    map(~ InvWaveTransCoefs(.x, coef1 = coef1, level = level)) %>%
    map_dfc(~as.vector(.)) %>% 
    as.matrix()
  recon_coef = sum((SVD - truth)^2)/prod(dim(truth))
  
  return(list(gd_raw = gd_raw, gd_coef = gd_coef, recon_raw = recon_raw, recon_coef = recon_coef))
  
}
#pipe_SNR(2,truth,9)
#R = 1


## cross validation for decomposition on STE raw data pipelines
## @ raws raw data matrix with dim is D^2 by Number of genes
## @ testID ID indicating the test entry (vectorized ID, 
##      e.g.raws is a matrix: id = which(raws != 0); testID = sample(id, length(id)/k))
## @ K designated number of factors
## @ threshold threshold for selecting K number of factors, not using if K is specified
## @ replacement bool value indicating whether select test splits with replacement, if TRUE, it's
##   equivalent to select test observations with replacement k times 

pipe_decom_raw_STEdata = function(raws, testID, K = NULL, threshold = 500, replacement = TRUE) {
  # raws = raws %>% as.matrix()
  # 
  # id = which(raws != 0)
  # 
  # testID = sample(id, length(id)/k)
  
  train = raws
  train[testID] = 0

  ##########===================  SVD
  s = svd(train)
  if(is_null(K)){
    K = sum(s$d > threshold)
  }
  
  u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
  
  recon_SVD = u%*%D%*%t(v)
  error_SVD = sum((recon_SVD[testID] - raws[testID])^2)/length(testID)
  
  ##########==================== EBMF
  res = EBMF_bkft(as.matrix(train), K = K, initMethod = "greedy", tol = 1e-10, maxIter = 30)
  
  recon_EBMF = res$LL%*% t(res$FF)
  error_EBMF = sum((recon_EBMF[testID] - raws[testID])^2)/length(testID)

  return(list(error_SVD = error_SVD, error_EBMF = error_EBMF))
}
 
## hybrid thresholding
## cross validation for decomposition on wavelet coefficients of STE data pipelines using hybrid thresholding
## @ raws raw data matrix with dim is D^2 by Number of genes
## @ testID ID indicating the test entry (vectorized ID, 
##      e.g.raws is a matrix: id = which(raws != 0); testID = sample(id, length(id)/k))
## @ K designated number of factors
## @ threshold threshold for selecting K number of factors, not using if K is specified
## @ replacement bool value indicating whether select test splits with replacement, if TRUE, it's
##   equivalent to select test observations with replacement k times 
## @ wf wavelet filter
## @ J wavelet level
## @ ... pass thresholdMethod in WaveTransCoefs(), default for None
pipe_decom_wave_STEdata = function(raws, testID, K = NULL, threshold = 500, replacement = TRUE, wf = "d4", J = 5, ...) {
  # raws = raws %>% as.matrix()
  # 
  # id = which(raws != 0)
  # 
  # testID = sample(id, length(id)/k)
  
  train = raws
  train[testID] = 0
  
  ####======== This part save partial info for inverse wavelet transformation ####========
  ## gene1 is col1
  col1 = train[,1]
  
  m = matrix(col1, nrow = sqrt(length(col1)))
  ## save coefficient list attribute info
  coef1 = dwt.2d(m, wf = wf, J = J)
  
  
  ## save coefficient long vector level info
  coef1M = coef1%>%
    map_dfr(~ tibble(coef = as.numeric(.)), .id = "level") #; coef1M
  level = coef1M$level
  ####======== This part save partial info for inverse wavelet transformation ####========
  
  ## apply WaveTransCoefs and threshold to each column of raws x
  coefs_matrix = 1:dim(train)[2] %>%
    map_dfc(~ WaveTransCoefs(train[,.x], wf = wf, J = J, ...))
  
  ##########===================  SVD
  s = svd(coefs_matrix)
  if(is_null(K)){
    K = sum(s$d > threshold)
  }
  u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
  
  SVD = u%*%D%*%t(v)
  recon_SVD = split(SVD, col(SVD)) %>%
    map(~ InvWaveTransCoefs(.x, coef1 = coef1, level = level)) %>%
    map_dfc(~as.vector(.)) %>% 
    as.matrix()
  error_SVD = sum((recon_SVD[testID] - raws[testID])^2)/length(testID)
  
  ##########==================== EBMF
  
  res = EBMF_bkft(as.matrix(coefs_matrix), K = K, initMethod = "greedy", tol = 1e-10, maxIter = 30)
  
  EBMF = res$LL%*% t(res$FF)
  recon_EBMF = split(EBMF, col(EBMF)) %>%
    map(~ InvWaveTransCoefs(.x, coef1 = coef1, level = level)) %>%
    map_dfc(~as.vector(.)) %>% 
    as.matrix()
  error_EBMF = sum((recon_EBMF[testID] - raws[testID])^2)/length(testID)

  
  return(list(error_SVD = error_SVD, error_EBMF = error_EBMF))
}

#####====================== Manual thresholding
## private helper function
## for each vectorized D^2 vector, compute its threshold coefficients, return a long vectorized coefficients
## @ tau wavelet threshold in manual thresholding
WaveTransCoefsManual = function(x, wf = "d4", J = 4 ,tau) {
  ## contruct it as coor matrix image
  m = matrix(x, nrow = sqrt(length(x)))
  ## apply wavelet transformation and threshold
  coef = dwt.2d(m, wf = wf, J = J)
  coef_manual = manual.thresh.dwt.2d(coef, value = tau,  hard = TRUE)
  ## populate it into long vector whose length is the number of coefficients
  coefs = coef_manual%>%
    map_dfr(~ tibble(coef = as.numeric(.)))
  return(coefs)
}

## cross validation for decomposition on wavelet coefficients of STE data pipelines using hybrid thresholding
## @ raws raw data matrix with dim is D^2 by Number of genes
## @ testID ID indicating the test entry (vectorized ID, 
##      e.g.raws is a matrix: id = which(raws != 0); testID = sample(id, length(id)/k))
## @ K designated number of factors
## @ threshold threshold for selecting K number of factors, not using if K is specified
## @ replacement bool value indicating whether select test splits with replacement, if TRUE, it's
##   equivalent to select test observations with replacement k times 
## @ wf wavelet filter
## @ J wavelet level
## @ tau wavelet threshold in manual thresholding
pipe_decom_wave_STEdata_manual = function(raws, testID, K = NULL, threshold = 500, replacement = TRUE, wf = "d4", J = 5, tau) {
  # raws = raws %>% as.matrix()
  # 
  # id = which(raws != 0)
  # 
  # testID = sample(id, length(id)/k)
  
  train = raws
  train[testID] = 0
  
  ####======== This part save partial info for inverse wavelet transformation ####========
  ## gene1 is col1
  col1 = train[,1]
  
  m = matrix(col1, nrow = sqrt(length(col1)))
  ## save coefficient list attribute info
  coef1 = dwt.2d(m, wf = wf, J = J)
  
  
  ## save coefficient long vector level info
  coef1M = coef1%>%
    map_dfr(~ tibble(coef = as.numeric(.)), .id = "level") #; coef1M
  level = coef1M$level
  ####======== This part save partial info for inverse wavelet transformation ####========
  
  ## apply WaveTransCoefs and threshold to each column of raws x
  coefs_matrix = 1:dim(train)[2] %>%
    map_dfc(~ WaveTransCoefs(train[,.x], wf = wf, J = J, thresholdMethod = "manual", tau = tau))
    #map_dfc(~ WaveTransCoefsManual(train[,.x], wf = wf, J = J,  tau = tau))
  
  ##########===================  SVD
  s = svd(coefs_matrix)
  if(is_null(K)){
    K = sum(s$d > threshold)
  }
  u = s$u[,1:K]; D =  diag(s$d[1:K]); v = s$v[,1:K]
  
  SVD = u%*%D%*%t(v)
  recon_SVD = split(SVD, col(SVD)) %>%
    map(~ InvWaveTransCoefs(.x, coef1 = coef1, level = level)) %>%
    map_dfc(~as.vector(.)) %>% 
    as.matrix()
  error_SVD = sum((recon_SVD[testID] - raws[testID])^2)/length(testID)
  
  ##########==================== EBMF
  
  res = EBMF_bkft(as.matrix(coefs_matrix), K = K, initMethod = "greedy", tol = 1e-10, maxIter = 30)
  
  EBMF = res$LL%*% t(res$FF)
  recon_EBMF = split(EBMF, col(EBMF)) %>%
    map(~ InvWaveTransCoefs(.x, coef1 = coef1, level = level)) %>%
    map_dfc(~as.vector(.)) %>% 
    as.matrix()
  error_EBMF = sum((recon_EBMF[testID] - raws[testID])^2)/length(testID)
  
  
  return(list(error_SVD = error_SVD, error_EBMF = error_EBMF))
}


###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== ###=== 
###=== ###=== ###=== ###=== ###=== ###===###=== ###=== ###=== ###=== ###=== ###===###=== ###=== ###=== ###=== ###=== ###===###=== ###=== ###=== ###=== ###=== ###===


## error for each gene column with CV (random hold out dataset only each column)
## This function run replicates simulations, it parallel on each gene
## functionized pipeline for error
## @ raws: the raw D^2*P matrix
## @ colID indicating which column to be compute
## @ Nsim: replicates for folds (K-fold CV)
## @ ...: other arguments pass to decomp function
error_per_gene_CV = function(raws, colID, Nsim, ...){
  id = which(raws[,colID] != 0) 
  
  ## record error for gene colID
  error_SVD_raw = c(); error_EBMF_raw = c()
  error_SVD_wave = c(); error_EBMF_wave = c()
  
  for(simID in 1:Nsim){
    k = 5 ## take a portion of the data, but not really need to be k-fold, since we sample with replacement
    testID = sample(id, length(id)/k)
    train = raws
    train[testID, colID] = 0
    
    ########==============  decomp on raw
    res_raw = decomp(train, decom_method = "raw", ...)
    error_SVD = sum((res_raw$SVD$recon_SVD[testID, colID] - raws[testID, colID])^2)/length(testID)
    error_SVD_raw = append(error_SVD_raw, error_SVD)
    
    error_EBMF = sum((res_raw$EBMF$recon_EBMF[testID, colID] - raws[testID, colID])^2)/length(testID)
    error_EBMF_raw = append(error_EBMF_raw, error_EBMF)
    
    
    ########==============  decomp on wave
    res_wave = decomp(train, decom_method = "wave", ...)
    error_SVD = sum((res_wave$SVD$recon_SVD[testID, colID] - raws[testID, colID])^2)/length(testID)
    error_SVD_wave = append(error_SVD_wave, error_SVD)
    
    error_EBMF = sum((res_wave$EBMF$recon_EBMF[testID, colID] - raws[testID, colID])^2)/length(testID)
    error_EBMF_wave = append(error_EBMF_wave, error_EBMF)
    
    print("---=-====", simID)
  }
  
  df = data.frame(SVD_raw = error_SVD_raw, EBMF_raw = error_EBMF_raw, 
                  SVD_wave = error_SVD_wave, EBMF_wave = error_EBMF_wave,
                  gene = colID)
  
  return(df)
}


## recons for each gene column with CV (random hold out within each gene)
## This function use just one replicate, for plotting purposes
## functionized pipeline for plotting
## @ raws: the raw D^2*P matrix
## @ colID indicating which column to be compute
## @ ...: other arguments pass to decomp function

## return the reconstruction and testID
recons_per_gene_CV = function(raws, colID, ...){
  id = which(raws[,colID] != 0) 

  k = 5 ## take a portion of the data, but not really need to be k-fold, since we sample with replacement
  testID = sample(id, length(id)/k)
  train = raws
  train[testID, colID] = 0
  
  ########==============  decomp on raw
  res_raw = decomp(train, decom_method = "raw", ...)
  
  
  ########==============  decomp on wave
  res_wave = decomp(train, decom_method = "wave", ...)

  df = data.frame(SVD_raw = res_raw$SVD$recon_SVD[,colID], EBMF_raw = res_wave$SVD$recon_SVD[, colID], 
                  SVD_wave = res_wave$SVD$recon_SVD[, colID], EBMF_wave = res_wave$EBMF$recon_EBMF[, colID],
                  test = 0)
  df$test[testID] = 1
  
  return(df)
}

## 2/8/22
##====================##====================##====================##====================##====================##====================
##====================##====================##====================##====================##====================
###===  error for each gene column with CV (random hold out dataset whole raws)
## cross validation for decomposition (SVD_EBMF_raw_wavelet) of STE data pipelines
## The idea is same as function error_per_gene_CV
## But this function run for loop for each gene, parallel for simulation
## @ raws raw data matrix with dim is D^2 by Number of genes

## below are variables passed to decomp function
## @ K designated number of factors
## @ threshold threshold for selecting K number of factors, not using if K is specified
## @ replacement bool value indicating whether select test splits with replacement, if TRUE, it's
##   equivalent to select test observations with replacement k times 
## @ wf wavelet filter
## @ J wavelet level
## @ tau wavelet threshold in manual thresholding
pipe_error_per_gene_STEdata = function(raws, ...) {
  
  ## @ testID ID indicating the test entry (vectorized ID, 
  ##      e.g.raws is a matrix: id = which(raws != 0); testID = sample(id, length(id)/k))
  k=5
  raws = raws %>% as.matrix()
  
  ## contruct a parallel boolean matrix same as raws, indicating whether such entry is in test set
  isTest = matrix(FALSE, dim(raws)[1], dim(raws)[2]) ##
  
  id = which(raws != 0)
  testID = sample(id, length(id)/k)
  
  train = raws
  train[testID] = 0
  isTest[testID] = TRUE
  
  ## record error for gene colID
  error_SVD_raw = c(); error_EBMF_raw = c()
  error_SVD_wave = c(); error_EBMF_wave = c()
  
  ########==============  decomp on raw
  res_raw = decomp(train, decom_method = "raw", ...)
  
  ########==============  decomp on wave
  res_wave = decomp(train, decom_method = "wave", ...)
  
  
  for(colID in 1:dim(raws)[2]) {
    error_SVD = sum((res_raw$SVD$recon_SVD[isTest[,colID], colID] - raws[isTest[,colID], colID])^2)/length(isTest[,colID])
    error_SVD_raw = append(error_SVD_raw, error_SVD)
    
    error_EBMF = sum((res_raw$EBMF$recon_EBMF[isTest[,colID], colID] - raws[isTest[,colID], colID])^2)/length(isTest[,colID])
    error_EBMF_raw = append(error_EBMF_raw, error_EBMF)
    
    error_SVD = sum((res_wave$SVD$recon_SVD[isTest[,colID], colID] - raws[isTest[,colID], colID])^2)/length(isTest[,colID])
    error_SVD_wave = append(error_SVD_wave, error_SVD)
    
    error_EBMF = sum((res_wave$EBMF$recon_EBMF[isTest[,colID], colID] - raws[isTest[,colID], colID])^2)/length(isTest[,colID])
    error_EBMF_wave = append(error_EBMF_wave, error_EBMF)
  }
  
  df = data.frame(SVD_raw = error_SVD_raw, EBMF_raw = error_EBMF_raw, 
                  SVD_wave = error_SVD_wave, EBMF_wave = error_EBMF_wave,
                  gene = 1:dim(raws)[2])
  
  
  return(df)
}


## recons for each gene column with CV (random hold out for whole plot)
## This function use just one replicate, for plotting purposes
## functionized pipeline for plotting
## @ raws: the raw D^2*P matrix
## @ ...: other arguments pass to decomp function

## return the reconstruction and testID
pipe_recons_per_gene = function(raws, ...){
  id = which(raws != 0) 
  
  k = 5 ## take a portion of the data, but not really need to be k-fold, since we sample with replacement
  
  ## contruct a parallel boolean matrix same as raws, indicating whether such entry is in test set
  isTest = matrix(FALSE, dim(raws)[1], dim(raws)[2]) ##
  
  testID = sample(id, length(id)/k)
  train = raws
  train[testID] = 0
  isTest[testID] = TRUE
  
  ########==============  decomp on raw
  res_raw = decomp(train, decom_method = "raw", ...)
  
  
  ########==============  decomp on wave
  res_wave = decomp(train, decom_method = "wave", ...)

  return(list(SVD_raw = res_raw$SVD$recon_SVD, EBMF_raw = res_wave$SVD$recon_SVD,
              SVD_wave = res_wave$SVD$recon_SVD, EBMF_wave = res_wave$EBMF$recon_EBMF, 
              isTest = isTest))
}

