#' Apply Givens rotation to a matrix in a fixed position.
#' 
#' The function \code{givensrot} apply Givens rotation to a matrix in the (i, j) position.
#'  
#' @param M a matrix.
#' @param i indicates the row to apply Givens rotation.
#' @param j indicates the column to apply Givens rotation.
#' 
#' @details To convert the covariance matrix to a correlation matrix without changing the 
#' eigenvalues, a series of N - 1 Givens rotations (Bendel and Mickey 1978) are needed, following 
#' the algortithm of Davies and Higham (2000). The functions \code{givensrot} is used internally
#' by the function \code{fullgivens} to perform this serie of rotations. The code for this function
#' is based on the SAS routine detailed in Wicklin (2013, Appendix C).
#' 
#' @return a matrix.
#' 
#' @seealso \code{\link{fullgivens}}  
#' 
#' @author Santiago Benitez-Vieyra, based in the SAS routine of Wicklin (2013).
#' 
#' @references Bendel, R.B. and Mickey, M.R. 1978. Population correlation matrices for sampling
#' experiments. Communications in Statistics—Simulation and Computation 7:163–182.
#' @references Davies, P.I. and Higham, N.J. 2000. Numerically stable generation of correlation
#' matrices and their factors. BIT 40:640–651.
#' @references Wicklin, R. 2013. Simulating data with SAS. Cary, NC:SAS Institute Inc.
#' 
#' @export
#'  
givensrot <- function(M, i, j){
  Aii <- M[i,i]
  Aij <- M[i,j]
  Ajj <- M[j,j]
  t <- (Aij + sqrt(Aij^2 - (Aii-1)*(Ajj-1))) / (Ajj - 1)
  c <- 1/sqrt(1+t^2)
  s <- c*t
  Ai <- A[i,]
  Aj <- A[j,] 
  A[i,] = c*Ai - s*Aj
  A[j,] = s*Ai + c*Aj
  Ai <- A[,i]
  Aj <- A[,j]
  A[,i] <- c*Ai - s*Aj
  A[,j] <- s*Ai + c*Aj
  A
}