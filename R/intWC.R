#' Wagner-Cheverud phenotypic integration index
#' 
#' intWC estimates the relative variance of the eigenvalues of a correlation 
#' matrix, corrected by sample size.
#' 
#' @param data A data frame or matrix with two or more numeric variables
#' @param correct Whether or not to correct eigenvalue variance by sample size. 
#' Default is \code{TRUE}. 
#' @param type A character string indicating if the data is data frame or matrix 
#'  with raw data (\code{'raw.data'}) or a correlation matrix (\code{'cor'}).
#' @param ... Additional arguments to be passed to function \code{cor}
#' 
#' @details The variance of the eigenvalues is estimated as a population
#' variance (i.e. over N instead of N-1), where N is the number 
#' of traits. The variance of eigenvalues ranges from zero, when all 
#' eigenvalues are equal, to a maximum of N-1, when only one 
#' eigenvalue is larger than zero. To allow comparisons among samples 
#' with different number of traits, the variance of eigenvalues is displayed 
#' as a percentage of this maximum. When \code{correct = TRUE} (default) the 
#' variance of eigenvalues is corrected by substracting the expected value of 
#' eigenvalue variance for finite sample correlation matrices for uncorrelated 
#' variables, i.e. (N-1)/M where N is the number of traits and M is 
#' the sample size.
#' 
#' @return the Wagner-Cheverud phenotypic integration index (INT). 
#'  
#' @references
#' Cheverud, J.M., Wagner, G.P., and Dow, M.M. 1989. Methods for 
#' the comparative analysis of variation patterns. Systematic Zoology 38:201-213.
#' @references
#' Wagner, G.P. 1984. On the eigenvalue distribution of genetic and 
#' phenotypic dispersion matrices: Evidence for a nonrandom 
#' organization for quantitative character variation. Journal of 
#' Mathematical Biology 21:77-95.
#'  
#' @examples
#' A <- sorbus[, 1:7]
#' intWC(A, type = "raw.data")
#' 
#' B <- cor(sorbus[, 1:7])
#' intWC(B, type = "cor")
#'  
#' @export
#' 
intWC <- function(data, correct = TRUE, type = c("cor", "raw.data"), ...){
  if(type == "cor"){
    EI <- eigen(data, symmetric = TRUE)$values
    NC <- ncol(data)
    INT <- ((sum((EI - mean(EI))^2)/NC))/(NC - 1)*100 
    return(INT)
  } else {
    if(type == "raw.data"){
      M <- cor(data, ...)
      EI <- eigen(M, symmetric = TRUE)$values
      NC <- ncol(data)
      if (correct == TRUE){
        E <- (NC - 1)/nrow(data)
        INT <- ((sum((EI - mean(EI))^2)/NC) - E)/(NC - 1)*100
        return(INT)
      } else {
        INT <- ((sum((EI - mean(EI))^2)/NC))/(NC - 1)*100 
        return(INT)
      }  
    } else {
      stop("invalid 'type' argument")
    }
  }
}