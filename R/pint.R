#' Integration pseudovalues
#' 
#' pint estimates the individual jackknife pseudovalues of the Wagner-Cheverud
#' morphological integration index (INT).
#' 
#' @param data A data frame or matrix with two or more numeric variables
#' @param correct Whether or not to correct INT by sample size. 
#' Default is \code{TRUE}. 
#' @param ... Additional arguments to be passed to function \code{cor}
#' 
#' @details The Wagner-Cheverud index of morphological integration (Wagner 1984, Cheverud 
#' et al. 1989) is first estimated using the complete sample of individuals. Then, one
#' individual is ignored in turn, and the same statistic is computed. Finally, the pseudovalues are estimated according
#' to the usual jackknife procedure (Tuckey 1958, Manly 1997). The integration pseudovalue 
#' of the i individual is defined as
#' 
#' \deqn{pINT_{i} = (N)INT - (N-1) INT_{-i}}{pINT_i = N*INT - (N-1)*INT_-i}
#' 
#' where N is the number of individuals, INT is the Wagner-Cheverud index of morphological
#' integration and -i indicates that the i individual was removed before estimating INT.
#' when \code{correct = TRUE} (default) INT is corrected by substracting to the variance
#' of eigenvalues is corrected their expected value for finite sample correlation matrices
#' with uncorrelated variables, i.e. (N-1)/M where N is the number of traits and M is 
#' the sample size.
#' The individual pINT values allow to estimate change in morphological integration due to
#' phenotypic selection, avoiding the calculation of the change in the full 
#' variance-covariance phenotypic matrix, following Benitez-Vieyra et al. (in prep).
#'  
#' @return the vector of individual jackknife pseudovalues of the Wagner-Cheverud 
#' morphological integration index (INT).
#' 
#' @seealso \code{\link{intWC}}  
#'  
#' @author Santiago Benitez-Vieyra 
#'  
#' @references
#' Benitez-Vieyra, S., Fornoni, J., Dominguez, C.A. in prep.  
#' @references
#' Cheverud, J.M., Wagner, G.P., and Dow, M.M. 1989. Methods for 
#' the comparative analysis of variation patterns. Systematic Zoology 38:201-213.
#' @references
#' Manly, B.J.S., 1997. Randomization, Bootstrap and Monte Carlo methods in Biology. 
#' CRC Press, Boca Raton.
#' @references
#' Wagner, G.P. 1984. On the eigenvalue distribution of genetic and 
#' phenotypic dispersion matrices: Evidence for a nonrandom organization for 
#' quantitative character variation. Journal of Mathematical Biology 21:77-95.
#' @references
#' Tukey, J.W. 1958. Bias and confidence in not-quite large samples. Annals of 
#' Mathematical Statistics 29:614.
#'  
#' @examples
#' pint(sorbus[, 1:7])
#'  
#' @export
#' 
pint <- function(data, correct = TRUE, ...){
  N <- nrow(data) 
  Tint <- intWC(data, correct = correct, ...)
  Lint <- numeric(N)  
  for(i in 1:N) Lint[i] <- intWC(data[-i, ], ...) 
  Dint <- N*Tint - (N-1)*Lint
  Dint  
}