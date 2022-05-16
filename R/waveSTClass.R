###Class Definition

#'S4 Class for union of data.frame, matrix, and tibble
#'
#'An S4 class for a union of data.frame, matrix, and tibble
#'@rdname dfOrMtx
#'@aliases dfOrMtx
#'@docType class
#'@exportClass dfOrMtx
setClassUnion("dfOrMtx", c("data.frame", "matrix", "tbl_df"))

#'S4 Class for a waveST
#'
#'The waveST class stores the data and additional result after computing.
#'
#'@section Slots:
#' \describe{
#'	\item{data}{original data, cell(n)-by-gene(p) matrix, can be matrix/data.frame/tibble}
#'  \item{spatial}{cell(n)-by-location(2) matrix, representing coordinates of each cell, see "details"}
#'  \item{input}{location(D^2)-by-gene(2) matrix, representing gene expression over location}
#'  \item{output}{a list containing SVD-like output: factors, diagonal matrix, loadings, fitted values, see "details"}
#' }
#'@name waveST-class
#'@rdname waveST-class
#'@aliases waveST-class
#'@docType class
#'@exportClass waveST
#'@author Zhuoyan Xu \email{zhuoyan.xu@wisc.edu}
#'@note The decomposition in this package require a waveST input.
#'@details   If \code{spatial} is not given, the class assumes the data is generated input: location(D^2)-by-gene(2).
#'           The class set \code{input} = \code{data}, this is useful when we already have location(D^2)-by-gene(2)
#'           matrix (like settings in simulations).
#'           
#'           The content of \code{output} depends on whether use wavelet, and the decomposition methods.
#'           The  \code{output} is NA once created. The values will be updated once called \code{\link{decompse}}.
#'@export
#'@examples 
#'wave = new("waveST", data = matrix(rnorm(6),2,3))
setClass("waveST", representation(data = "dfOrMtx", spatial = "dfOrMtx", input = "dfOrMtx",
                                  output = "list"), 
         prototype = list(
           spatial = matrix(),
           input = matrix(),
           output = list()
         ))


###Creation of waveST from data
#'waveST Constructor
#'
#'Create a \code{\link{waveST}} object from an \code{data.frame}, \code{matrix}, or \code{tibble}.
#'@docType methods
#'@export
#'@name waveST
#'@rdname waveST
#'@aliases waveST
#'@param data an instance of \code{data.frame}, \code{matrix}, or \code{tibble}, representing riginal data, cell(n)-by-gene(p) matrix,
#'@param spatial cell(n)-by-location(2) matrix, representing coordinates of each cell, see "details"
#'@param input location(D^2)-by-gene(2) matrix, representing gene expression over location
#'@param level the scale for the gene expression image after input generation, the size of images for input generation should be 2^level
#'@return a \code{\link{waveST}} object
#'@details  If \code{spatial} is not given, the class assumes the data is generated input: location(D^2)-by-gene(2).
#'          The class set \code{input} = \code{data}, this is useful when we already have location(D^2)-by-gene(2)
#'          matrix (like settings in simulations).
#'           
#'          The content of \code{output} depends on whether use wavelet, and the decomposition methods.
#'          The  \code{output} is NA once created. The values will be updated once called \code{\link{decompse}}.
#'@examples 
#'res = kOverA_ST(k = 3,A = 7)
#'viz = res$viz; df = res$df
#'wave = waveST(data = df, spatial = viz)
waveST <- function(data, spatial = NA, input = NA, level = 6) {
  if(sum(is.na(spatial)) == 1) {
    spatial = matrix()
    input = data
  } else{
    # if given raw data and saptial matrix, run Algorithm 1: Input Generation
    ## transfer n*p matrix to grid by gene D^2*p matrix
    
    ## gridValue is coor image, gridls is a long list with length p (num of genes),
    gridls = gridGeneRaw(spatial, data, level = level)
    
    suppressMessages(
      input <- gridls %>% 
      map_dfc(~as.vector(.))
    )
  }
  new("waveST", data = data,  spatial = spatial,  input = input)
}



## Main workflow shown in diagram

#'waveST decomposition
#'
#'Decompose input in \code{\link{waveST}} object with wavelet transformation and threshold
#'@docType methods
#'@export
#'@name decompose
#'@rdname decompose
#'@aliases decompose
#'@param waveST an \code{waveST} class
#'@param wavemethod a vector showing whether whether we use wavelet transformation, "raw" means we decompose 
#'                 directly without wavelet technique, "wave" means we apply wavelet transformation and shrinkage technique
#'@param decom_method Matrix decomposition method, "SVD" or "EBMF"
#'@param K The estimated number of factors, use "elbow" method to estimate it if not given
#'@param bar threshold for selecting number of factors, \code{K}, useless if \code{is already given}, see "details"
#'@param wf name of the wavelet filter to use
#'@param J depth of the wavelet basis decomposition, must be a number less than or equal to log(min{M,N},2)
#'@param thresholdMethod wavelet shrinkage method used  "hybrid" thresholding or "manual" thresholding, see "details"
#'@param tau constant threshold when using manual thresholding
#'@return a \code{\link{waveST}} object with \code{output} updated, see "details"
#'@details  If the number of factors \code{K} is not given, the use \code{bar} to select \code{K}, it 
#'          conduct SVD first,only keep the singular values larger than \code{bar}, the number of 
#'          kept values are selected \code{K}.
#'          
#'          The "hybrid" thresholding or "manual" thresholding, see \code{\link{WaveTransCoefs}}.
#'          
#'          The output contains:
#'          \enumerate{
#'          \item \code{f} factor matrix, representing factor genes
#'          \item \code{D} middle diagonal matrix 
#'          \item \code{l} loading matrix
#'          \item \code{recon} The fitted reconstruction matrix, estimated data
#'          }
#'          
#'          If wavelet method is used, \code{f}, \code{D}, \code{l} cooresponding to the decomposition 
#'          result for coefficient matrix. \code{recon} is still the fitted matrix of original input.
#'  
#'@examples 
#'res = kOverA_ST(k = 3,A = 7)
#'viz = res$viz; df = res$df
#'wave = waveST(data = df[,1:5], spatial = viz)
#'wave = decompose(wave, "raw", "SVD", K = 5)
decompose = function(waveST, wavemethod = c("raw", "wave"), decom_method = c("SVD", "EBMF"),  
                     K = NULL, bar = 500, wf = "d4", J = 5, thresholdMethod = "manual", tau = 0){
  stopifnot(is(waveST) == "waveST")
  
  raws = waveST@input
  if(wavemethod == "raw"){
    if(decom_method == "SVD"){
      waveST@output = decomp2(raws, wavemethod = "raw", decom_method = "SVD", K = K, threshold = bar)
    } else if(decom_method == "EBMF") {
      waveST@output = decomp2(raws, wavemethod = "raw", decom_method = "EBMF", K = K, threshold = bar)
    }
  } else if(wavemethod == "wave"){
    
    if(decom_method == "SVD"){
      waveST@output = decomp2(raws, wavemethod = "wave", decom_method = "SVD", K = K, threshold = bar,
                              wf = wf, J = J, thresholdMethod = thresholdMethod, tau = tau)
    } else if(decom_method == "EBMF") {
      waveST@output = decomp2(raws, wavemethod = "wave", decom_method = "EBMF", K = K, threshold = bar,
                              wf = wf, J = J, thresholdMethod = thresholdMethod, tau = tau)
    }
    
  }
  waveST
}

