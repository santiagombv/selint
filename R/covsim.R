#' Generation of covariance matrices.
#' 
#' The function covsim generates a covariance matrix from a given set of eigenvalues.
#' 
#' @param d a numeric vector of eigenvalues.
#' 
#' @details The function uses a given set of eigenvalues(\eqn{\lambda}), preferably obatained from 
#' \code{BSeigenval} or \code{MODeigenval}, and a random eigenvector structure (\eqn{V})  to
#' build a covariance matrix as \eqn{V\lambda V'}. The eigenvectors are obtained from a vector of
#' random uniform numbers between -1 and 1 and using Gram-Schmidt orthogonalization, from the 
#' \code{orthonormalization} function of \code{far} package. 
#' 
#' @return a covariance matrix.
#' 
#' @author Santiago Benitez-Vieyra
#' 
#' @examples
#' A <- BSeigenval(INT = 50, dim = 8, sigma = 0.1)
#' covsim(A)
#' 
#' @importFrom far orthonormalization
#' 
#' @export
#' 
covsim <- function(d){
  if (!requireNamespace("far", quietly = TRUE)) {
    stop("far package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  D <- diag(d)
  R <- runif(length(d), -1, 1)
  V <- orthonormalization(R, basis = T)
  A <- V%*%D%*%t(V) ## covariance
  A
}