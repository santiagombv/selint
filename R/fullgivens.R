#' Apply a serie of Givens rotation to a covariance matrix. It allows to trasform a covariance
#' matrix in a correlation matrix without changing the eignevalues.
#' 
#' \code{fullgivens} apply a series of Givens rotation to a matrix in the (i, j) position.
#'  
#' @param M a matrix.
#'  
#' @details To convert the covariance matrix to a correlation matrix without changing the 
#' eigenvalues, a series of N - 1 Givens rotations (Bendel and Mickey 1978) are needed, following 
#' the algortithm of Davies and Higham (2000). The functions \code{fullgivens} to perform this 
#' serie of rotations. The code for this function is based on the SAS routine detailed in Wicklin 
#' (2013, Appendix C).
#' 
#' @return a correlation matrix.
#' 
#' @seealso \code{\link{givensrot}}  
#' 
#' @author Santiago Benitez-Vieyra, based in the SAS routine of Wicklin (2013).
#' 
#' @references Bendel, R.B. and Mickey, M.R. 1978. Population correlation matrices for sampling
#' experiments. Communications in Statistics-Simulation and Computation 7:163-182.
#' @references Davies, P.I. and Higham, N.J. 2000. Numerically stable generation of correlation
#' matrices and their factors. BIT 40:640-651.
#' @references Wicklin, R. 2013. Simulating data with SAS. Cary, NC:SAS Institute Inc.
#' 
#' @examples
#' Z <- BSeigenval(INT = 40, dim = 4, sigma = 0.1)
#' Z
#' M <- covsim(Z)
#' N <-fullgivens(M)
#' eigen(N)$values
#'  
#' @export
#'  
fullgivens <- function(M){
  X <- sum(diag(M) > 1.001 | diag(M) < 0.999) 
  while(X >1){
    d <- diag(M)
    SUP <- d[d>1]
    INF <- d[d<1]
    sup.place <- which(d == max(SUP))
    inf.place <- which(d == min(INF))
    M <- givensrot(M = M, i = sup.place, j = inf.place)
    X <- sum(diag(M) > 1.001 | diag(M) < 0.999) 
  }
  diag(M) <- 1
  M
} 