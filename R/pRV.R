#' Pseudovalues of Escoufier's RV coefficient
#' 
#' pRV estimates the individual jackknife pseudovalues of the Escoufier's RV coefficient
#' between two sets of variables
#' 
#' @param data A data frame or matrix with two or more numeric variables
#' @param vars1 Vector indicating the variables to be assigned to 
#' module 1.
#' @param vars2 Vector indicating the variables to be assigned to 
#' module 2.
#' @param ... Additional arguments to be passed to function \code{cov}.
#' 
#' @details The RV coefficient (Escoufier 1973, Robert and Escoufier 1978) is first 
#' estimated using the complete sample of individuals. Then, one individual is ignored 
#' in turn, and the same statistic is computed. Finally, the pseudovalues are estimated 
#' according to the usual jackknife procedure (Tuckey 1958, Manly 1997). The pseudovalue 
#' of the i individual is defined as
#' 
#' \deqn{pRV_{i} = (N)RV - (N-1)RV_{-i}}{pRV_i = N*RV - (N-1)*RV_-i}
#' 
#' where N is the number of individuals, and -i indicates that the i individual was removed. before estimating INT.
#' The individual pRV values allow to estimate change in morphological modularity due to
#' phenotypic selection, avoiding the calculation of the change in the full 
#' variance-covariance phenotypic matrix, following Benitez-Vieyra et al. (in prep).
#'  
#' @return the vector of individual jackknife pseudovalues of the RV coeffcient.
#' 
#' @seealso \code{\link{pairRV}} 
#'  
#' @author Santiago Benitez-Vieyra 
#'  
#' @references
#' Benitez-Vieyra, S., Fornoni, J., Dominguez, C.A. in prep.  
#' @references
#' Escoufier, Y. 1973. Le traitement des variables vectorielles. Biometrics 29:751-760.
#' @references
#' Manly, B.J.S., 1997. Randomization, Bootstrap and Monte Carlo methods in Biology. 
#' CRC Press, Boca Raton.
#' @references
#' Robert, P., and Escoufier, Y. 1976. A unifying tool for linear multivariate 
#' statistical analysis: the RV-coefficient. Applied Statistics 25:257-265.
#' @references
#' Tukey, J.W. 1958. Bias and confidence in not-quite large samples. Annals of 
#' Mathematical Statistics 29:614.
#'  
#' @examples
#' A <- matrix(rnorm(100), 20, 5)
#' pRV(A, vars1 = c(1:3), vars2 = c(4:5))
#'  
#' @export

pRV <- function(data, vars1, vars2, ...){
  N <- nrow(data) 
  Lrv <- numeric(N)
  Trv <- pairRV(data = data, vars1 = vars1, vars2 = vars2, ...) 
  for(i in 1:N) {
    Lrv[i] <- pairRV(data[-i, ], vars1 = vars1, vars2 = vars2, ...) 
    }
  Drv <- N*Trv - (N-1)*Lrv
  Drv  
}