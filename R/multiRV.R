#' Multi-set Escoufier's RV coefficient
#' 
#' multiRV estimates the strength of the association between multiple sets of 
#' variables (modules) as the average of all pair-wise RV coefficients between sets.
#' 
#' @param dat A data frame or matrix with two or more numeric variables.
#' @param vars A list of vectors indicating the variables to be assigned to each module.
#' 
#' @details  This function follows the suggestion of Klingenberg (2009) to estimate
#' the strength of association among multiple sets of variables (modules), as the average 
#' among all pair-wise Escoufier's RV coefficients between sets. As the RV coeffcient, the 
#' multi-set RV takes values between zero, when modules are completely separated, and one,
#' when there are no modular structure. 
#'
#' @return The multi-set Escoufier's RV coefficient. 
#' 
#' @seealso \code{\link{pairRV}}  
#' 
#' @author Santiago Benitez-Vieyra
#'  
#' @references
#' Klingenberg, C.P. 2009. Morphometric integration and modularity in configurations of 
#' landmarks: tools for evaluating a priori hypotheses. Evolution and Development 
#' 11:405-421.
#'
#' @examples
#' A <- matrix(rnorm(200), 20, 10)
#' multiRV(A, vars = list(c(1:3), c(4:5), c(6:10)))
#'  
#' @export
#' 
multiRV <- function (data, vars){
  mat <- cov(data[, unlist(vars)])
  K <- c(1:length(vars))
  set <- t(combn(K, 2))
  rvlist <- numeric(nrow(set))
  for(i in 1:nrow(set)) rvlist[i] <- pairRV(data, vars1= vars[[set[i, 1]]], vars2=vars[[set[i, 2]]])
  res <- mean(rvlist)
  res
}
