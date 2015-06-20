#' Generation of correlation matrices.
#' 
#' \code{corsim} simulates a correlation matrix with a given set of eigenvalues.
#' 
#' @param d a numeric vector of eigenvalues.
#'  
#' @details The function uses a vector of eigenavlues to build a covariance matrix.
#' Then, applies a series of Givens rotations to transform the covariance matrix into 
#' a correlation matrix without changing the eigenvalues, following the procedure of
#' Davies and Higham (2000).
#' 
#' @return a correlation matrix.
#'
#' @author Santiago Benitez-Vieyra
#' 
#' @references Davies, P.I. and Higham, N.J. 2000. Numerically stable generation of correlation
#' matrices and their factors. BIT 40:640-651.
#' 
#' @examples
#' # simulation from a broken-stick distribution of eigenvalues
#' e.val1 <- BSeigenval(INT = 50, dim = 6, sigma=0.1, tol = 0.0001)
#' Z <- corsim(d = e.val1)
#' Z
#' 
#' eigen(Z)$values # check the eigenvalues
#' intWC(Z, type = "cor") # check the integration of the Z matrix
#'  
#' @export
#' 
corsim <- function(d){
  A <- covsim(d)
  B <- fullgivens(A)
  B
}