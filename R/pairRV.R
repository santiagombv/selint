#'Escoufier's RV coefficient
#'
#'pairRV estimates the strength of the association between two sets of variables
#'(modules) using the Escoufier's RV coefficient.
#'
#'@param data A data frame or matrix with two or more numeric variables.
#'@param vars1 Vector indicating the variables to be assigned to module 1.
#'@param vars2 Vector indicating the variables to be assigned to module 2.
#'@param ... Additional arguments to be passed to function \code{cov}.
#'  
#'@details  The Escoufier's RV coefficient is deffined as
#'  
#'  \deqn{RV = \frac{tr(S_{ij}S_{ji})}{\sqrt(tr(S_{i})tr(S_{j}))}}{% RV = 
#'  tr(S_ij*S_ji)/(tr(S_i)*tr(S_j))^1/2}
#'  
#'  where \eqn{S_i} and \eqn{S_j} represent the variance-covariance matrices of 
#'  the i and j set of variables (modules), \eqn{S_{ij}}{S_ij} is the covariance
#'  matrix between these two sets (and \eqn{S_{ji}}{S_ji} is its transpose), and
#'  \eqn{tr} is the trace of the matrices, calculated as the sum of its diagonal
#'  elements. Therefore, RV represents the amount of covariation scaled by the 
#'  amounts of variation within two groups of variables, which is analogous to 
#'  the calculation of the (squared) correlation coefficient between two 
#'  variables (Klingenberg 2009). RV takes values between zero, when two modules
#'  are completely separated, and one, when there are no modular structure.
#'  
#'@return The Escoufier's RV coefficient.
#'  
#'@references Escoufier, Y. 1973. Le traitement des variables vectorielles. 
#'  Biometrics 29:751-760.
#'@references Klingenberg, C.P. 2009. Morphometric integration and modularity in
#'  configurations of landmarks: tools for evaluating a priori hypotheses. 
#'  Evolution and Development 11:405-421.
#'@references Robert, P., and Escoufier, Y. 1976. A unifying tool for linear 
#'  multivariate statistical analysis: the RV-coefficient. Applied Statistics 
#'  25:257-265.
#'  
#'@examples
#'pairRV(sorbus, vars1 = c(3, 4), vars2 = c(1, 2, 5:7))
#'  
#'@export
#'
pairRV <-function(data, vars1, vars2, ...){
  mat <- cov(data[, c(vars1, vars2)], ...)
  v1 <- seq(1, length(vars1))
  v2 <- seq(length(vars1)+1, length(c(vars1,vars2)))
  S12 <- mat[v1, v2]
  S11 <- mat[v1, v1]
  S22 <- mat[v2, v2]
  RV<-sum(S12^2)/sqrt(sum(S11^2)*sum(S22^2))
  RV
}