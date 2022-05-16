#### wavelets and shrinkage functions

#'wavelet transformation
#'
#'Performs a separable two-dimensional discrete wavelet transform (DWT) on a matrix of dyadic dimensions. Applying
#'indicated thresholds.
#'@docType methods
#'@export
#'@param x length D^2 vector, vectorized matrix for DWT
#'@param wf name of the wavelet filter to use
#'@param J depth of the wavelet basis decomposition, must be a number less than or equal to log(min{M,N},2)
#'@param thresholdMethod wavelet shrinkage method used  "hybrid" thresholding or "manual" thresholding, see "details"
#'@param tau constant threshold when using manual thresholding
#'@return a \code{tibble} with two columns, \code{level} and \code{coef}, see details
#'@details  The coefficients after wavelet thresformation and thresholding is a list structure 
#'          containing the 3J+1 sub-matrices from the decomposition. The \code{coef} is the vectorized 
#'          coefficients.
#'          
#'          The "hybrid" thresholding combined the SURE and universal thresholds. The manual thresholding set a constant threshold.
#'          see Section 3.2 in references.
#'
#'@references G. P. Nason. Wavelet methods in statistics with R. Springer, 2008.
WaveTransCoefs = function(x, wf = "d4", J = 4, thresholdMethod = NULL, tau = NULL) {
  ## contruct it as coor matrix image
  m = matrix(x, nrow = sqrt(length(x)))
  ## apply wavelet transformation and threshold
  coef = dwt.2d(m, wf = wf, J = J)
  if(!is_null(thresholdMethod)){
    if(thresholdMethod == "hybrid") coef = hybrid.thresh(coef)
    if(thresholdMethod == "manual") coef = manual.thresh.dwt.2d(coef, value = tau,  hard = TRUE)
  }
  
  ## populate it into long vector whose length is the number of coefficients
  #coefs = coef %>%
  #  map_dfr(~ tibble(coef = as.numeric(.)))
  coefs = tibble(coef = unlist(coef))
  return(coefs)
}



## for each long vectorized coefficients, compute its original vector,  return a dim DxD matrix
## @ x: length D^2 vector
## @ coef1 stores coefficient list attribute info
## @ level stores coefficient long vector level info
InvWaveTransCoefs = function(x, coef1, level) {
  coe_df = data.frame(x = x, level = level)
  coe_ls = coe_df %>%
    split(.$level) %>%
    map(~ matrix(.$x, nrow = sqrt(dim(.)[1])))
  
  # rearrange list to match the correct  order of coefficient list, and add attributes, for it fed into idwt.2d()
  coe_ls = coe_ls[names(coef1)]
  attributes(coe_ls) = attributes(coef1)
  
  reconsCM = idwt.2d(coe_ls)
  return(reconsCM)
}


#'Inverse Two-Dimensional Discrete Wavelet Transform
#'
#'Performs Inverse wavelet transformation
#'@docType methods
#'@export
#'@param x vectorized coefficients after wavelet transformation
#'@param raws original D^2 by P input
#'@return a D^2 matrix
InvWaveTrans = function(x, raws) {
  
  ####======== This part save partial info for inverse wavelet transformation ####========
  ## gene1 is col1
  raws = raws %>% as.matrix()
  col1 = raws[,1]
  
  m = matrix(col1, nrow = sqrt(length(col1)))
  ## save coefficient list attribute info
  coef1 = dwt.2d(m, wf = "d4", J = 5)
  
  
  ## save coefficient long vector level info
  coef1M = coef1 %>%
    map_dfr(~ tibble(coef = as.numeric(.)), .id = "level") #; coef1M
  level = coef1M$level
  ####======== This part save partial info for inverse wavelet transformation ####========
  
  coe_df = data.frame(x = x, level = level)
  coe_ls = coe_df %>%
    split(.$level) %>%
    map(~ matrix(.$x, nrow = sqrt(dim(.)[1])))
  
  # rearrange list to match the correct  order of coefficient list, and add attributes, for it fed into idwt.2d()
  coe_ls = coe_ls[names(coef1)]
  attributes(coe_ls) = attributes(coef1)
  
  reconsCM = idwt.2d(coe_ls)
  return(reconsCM)
}

##########===============##########===============##########===============##########===============##########===============##########===============
##########===============##########===============##########===============##########===============##########===============##########===============
##### 
manual.thresh.dwt.2d <- function(wc, max.level=4, value, hard=TRUE)
{
  wc.fine <- unlist(wc[1:3])
  factor <- median(abs(wc.fine)) / .6745
  
  wc.shrink <- wc
  
  if(hard) {
    # Hard thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * value
      wc.shrink[[i]] <- wci * (abs(wci) > unithresh)
    }
  }
  else {
    # Soft thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * value
      wc.shrink[[i]] <- sign(wci) * (abs(wci) - unithresh) * 
        (abs(wci) > unithresh)
    }
  }
  wc.shrink
}


universal.thresh.dwt.2d <- function(wc, max.level=4, hard=TRUE)
{
  n <- length(idwt.2d(wc))
  
  wc.fine <- unlist(wc[1:3])
  factor <- median(abs(wc.fine)) / .6745
  
  wc.shrink <- wc
  
  if(hard) {
    # Hard thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * sqrt(2 * log(n))
      wc.shrink[[i]] <- wci * (abs(wci) > unithresh)
    }
  }
  else {
    # Soft thresholding
    for(i in names(wc)[1:max.level]) {
      wci <- wc[[i]]
      unithresh <- factor * sqrt(2 * log(n))
      wc.shrink[[i]] <- sign(wci) * (abs(wci) - unithresh) *
        (abs(wci) > unithresh)
    }
  }
  wc.shrink
}


######============== 1. Threshold on wavelet coefficient
### below thresholding comes from sec 3.2 in G.P Nason Wavelet medthods in Statistics in R

threshold_coef = function(viz, geneID, hard = c(TRUE,FALSE),  ...){
  viz$gene = df[,geneID]
  gridValue = coorToGrid(viz)
  
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

