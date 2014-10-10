#' Escoufier's RV coefficient
#' 
#' pairRV estimates the strength of the association between two sets of 
#' variables (modules) using the Escoufier's RV coefficient.
#' 
#' @param mat A variance-covariance matrix, as obtained using \code{cov}.
#' @param vars1 Vector indicating the variables to be assigned to 
#' module 1.
#' @param vars2 Vector indicating the variables to be assigned to 
#' module 2.
#' 
#' @details  The Escoufier's RV coefficient is deffined as
#' 
#'\deqn{RV = \frac{tr(S_{ij}S_{ji})}{\sqrt(tr(S_{i})tr(S_{j}))}}{%
#'RV = tr(S_ij*S_ji)/(tr(S_i)*tr(S_j))^1/2}
#'
#' where \eqn{S_i} and \eqn{S_j} represent the variance-covariance matrices of 
#' the i and j set of variables (modules), \eqn{S_{ij}}{S_ij} is the covariance matrix 
#' between these two sets (and \eqn{S_{ji}}{S_ji} is its transpose), and \eqn{tr} 
#' is the trace of the matrices, calculated as the sum of its diagonal elements. 
#' Therefore, RV represents the amount of covariation scaled by the amounts of variation 
#' within two groups of variables, which is analogous to the calculation of the 
#' (squared) correlation coefficient between two variables (Klingenberg 2009). 
#' RV takes values between zero, when two modules are completely separated, 
#' and one, when there are no modular structure. 
#'
#' @return The Escoufier's RV coefficient. 
#'  
#' @references
#' Escoufier, Y. 1973. Le traitement des variables vectorielles. Biometrics 29:751-760.
#' @references
#' Klingenberg, C.P. 2009. Morphometric integration and modularity in configurations of 
#' landmarks: tools for evaluating a priori hypotheses. Evolution and Development 
#' 11:405-421.
#' @references
#' Robert, P., and Escoufier, Y. 1976. A unifying tool for linear multivariate 
#' statistical analysis: the RV-coefficient. Applied Statistics 25:257-265.
#'
#' @examples
#' A <- matrix(rnorm(100), 20, 5)
#' pairRV(cov(A), vars1 = c(1:3), vars2=c(4:5))
#'  
#' @export
#' 
pairRV <-function(mat, vars1, vars2){
  S12 <- mat[vars1, vars2]
  S11 <- mat[vars1, vars1]
  S22 <- mat[vars2, vars2]
  RV<-sum(S12^2)/sqrt(sum(S11^2)*sum(S22^2))
  RV
}