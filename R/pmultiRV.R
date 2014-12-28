#' Pseudovalues of multi-set Escoufier's RV coefficient
#' 
#' pmultiRV estimates the individual jackknife pseudovalues of the multi-set RV coefficient
#' between multiple sets of variables.
#' 
#' @param data A data frame or matrix with two or more numeric variables
#' @param vars A list of vectors indicating the variables to be assigned to each module.
#' @param ... Additional arguments to be passed to function \code{cov}.
#' 
#' @details The multi-set RV coefficient (Klingenberg 2009) is first estimated using the
#' complete sample of individuals, as the average among all pair-wise Escoufier's RV 
#' coefficients between sets. Then, one individual is ignored in turn, and the same 
#' statistic is computed. Finally, the pseudovalues are estimated according to the usual
#' jackknife procedure (Tuckey 1958, Manly 1997). The pseudovalue of the multi-set RV (RVm)
#' of the i individual is defined as
#' 
#' \deqn{pRVm_{i} = (N)RVm - (N-1)RVm_{-i}}{pRVm_i = N*RVm - (N-1)*RVm_-i}
#' 
#' where N is the number of individuals, and -i indicates that the i individual was removed. before estimating INT.
#' The individual pRVm values allow to estimate change in morphological modularity due to
#' phenotypic selection, avoiding the calculation of the change in the full 
#' variance-covariance phenotypic matrix, following Benitez-Vieyra et al. (in prep).
#'  
#' @return the vector of individual jackknife pseudovalues of the multi-set RV coeffcient.
#' 
#' @seealso \code{\link{multiRV}}
#'  
#' @author Santiago Benitez-Vieyra 
#'  
#' @references
#' Benitez-Vieyra, S., Fornoni, J., Dominguez, C.A. in prep.  
#' @references
#' Klingenberg, C.P. 2009. Morphometric integration and modularity in configurations of 
#' landmarks: tools for evaluating a priori hypotheses. Evolution and Development 
#' 11:405-421.
#' @references
#' Manly, B.J.S., 1997. Randomization, Bootstrap and Monte Carlo methods in Biology. 
#' CRC Press, Boca Raton.
#' @references
#' Tukey, J.W. 1958. Bias and confidence in not-quite large samples. Annals of 
#' Mathematical Statistics 29:614.
#'  
#' @examples
#' pmultiRV(sorbus,  vars = list(c(1,2), c(3,4), c(5:7)))
#'  
#' @export
#' 
pmultiRV <- function(data, vars, ...){
  N <- nrow(data) 
  Trv <- multiRV(data = data, vars = vars, ...)
  Lrv <- numeric(N)  
  for(i in 1:N) Lrv[i] <- multiRV(data = data[-i, ], vars = vars, ...)
  Drv <- N*Trv - (N-1)*Lrv
  Drv  
}
