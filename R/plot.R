## The plot function
theme_set(theme_minimal())

#'plot function for waveST class
#'
#'visualize the top factor genes and top 5 genes with high loadings on certain factor genes
#'@docType methods
#'@export
#'@param x An object of class \code{\link{waveST}}.
#@param y the missing signature we do not use
#@param waveST An object of class \code{\link{waveST}}.
#'@param k which factor gene to visualize
#'@param wave whether use wavelet method
#'@param ... Plot parameters forwarded.
#'@return A plot object.
#'@examples 
methods::setMethod(
  "plot",
  c(x = "waveST", y = "missing"),
  #c(waveST = "waveST", k = "numeric", wave = "logical"),
  function(x, k = 1, wave = TRUE, ...) {
  
  waveST <- x
  lay = layout(matrix(1:6,2,3,byrow = TRUE))
  ########============== show on raw
  if(!wave){
    waveST@output$f[,k] %>% 
      #normalize1() %>%
      matrix(nrow = sqrt(length(.))) %>%
      image(., asp = 1, xaxt = "n", yaxt = "n")
    
    
    ids = waveST@output$l[,k] %>%
      abs() %>%
      order(decreasing = T)
    
    waveST@input[,ids[1:5]] %>%
      #map(~.x/sqrt(sum(.x^2))) %>%
      map(~ matrix(.x, nrow = sqrt(length(.x)))) %>%
      map(~image(.x, asp = 1, xaxt = "n", yaxt = "n"))
  } else if(wave){
    waveST@output$f[,k] %>% 
      InvWaveTrans(., waveST@input %>% as.matrix()) %>%
      #image.plot(., asp = 1)
      image(., asp = 1, xaxt = "n", yaxt = "n")
    
    
    ids = waveST@output$l[,k] %>%
      abs() %>%
      order(decreasing = T)
    
    waveST@input[,ids[1:5]] %>% 
      map_dfc(~ WaveTransCoefs(.x,  wf = "d4", J = 5, thresholdMethod = "manual", tau = 40)) %>% 
      map(~ InvWaveTrans(.x, raws = waveST@input %>% as.matrix())) %>%
      map_dfc(~as.vector(.)) %>%
      map(~ matrix(.x, nrow = sqrt(length(.x)))) %>%
      #map(~image.plot(.x, asp = 1))
      map(~image(.x, asp = 1, xaxt = "n", yaxt = "n"))
  }
  layout(1)
})
