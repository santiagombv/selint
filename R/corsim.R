#' Generation of correlation matrices.
#' 
#' \code{corsim} simulates a correlation matrix with a given set of eigenvalues.
#' 
#' @param d a numeric vector of eigenvalues.
#' @param order A character string indicating the ordering method of the correlation
#' matrix. Use \code{'none'} for the original order, \code{'first'} for the order of 
#' the first eigenvector or \code{'angle'} for the angular order of eigenvectors.
#' with raw data (\code{'raw.data'}) or a correlation matrix (\code{'cor'}).
#'  
#' @details The function uses a vector of eigenvalues to build a covariance matrix.
#' Then, applies a series of Givens rotations to transform the covariance matrix into 
#' a correlation matrix without changing the eigenvalues, following the procedure of
#' Davies and Higham (2000). The argument order allows to arrange the matrix following 
#' the order of the first eigenvector (best method for correlation matrices simulated 
#' from broken-stick eigenvalue distribution) or by the angular order of eigenvectors
#' (best method for modular correlation matrices). The last method follows the procedure 
#' of Friendly (2002), but uses absolute correlations.
#' 
#' @return a correlation matrix.
#'
#' @author Santiago Benitez-Vieyra
#' 
#' @references Davies, P.I. and Higham, N.J. 2000. Numerically stable generation of correlation
#' matrices and their factors. BIT 40:640-651.
#' @references Friendly, M. 2002. Corrgrams. Exploratory Displays for Correlation Matrices. The 
#' American Statistician 56:316-324.
#' 
#' @examples
#' # simulation from a broken-stick distribution of eigenvalues
#' e.val1 <- BSeigenval(INT = 50, dim = 6, sigma=0.1, tol = 0.0001)
#' Z <- corsim(d = e.val1, order = "none")
#' Z
#' 
#' eigen(Z)$values # check the eigenvalues
#' intWC(Z, type = "cor") # check the integration of the Z matrix
#'  
#' @export
#' 
corsim <- function(d, order = c("none", "first", "angle")){
  A <- covsim(d)
  B <- fullgivens(A)
  if(order == "none"){
    B
  } else {
    if(order == "first"){
      o <- order(eigen(B)$vectors[, 1])
      B <- B[o, o]
      B
    } else {
      if(order == "angle"){
        ei <- eigen(abs(B))$vectors
        s1 <- which(ei[,1]>0)
        s2 <- which(ei[,1]<0)
        a <- numeric(length=nrow(ei))
        a[s1] <- atan(ei[s1, 2]/ei[s1, 1])
        a[s2] <- atan(ei[s2, 2]/ei[s2, 1]) + pi
        o <- order(a)
        B <- B[o, o]
        B
      }
    }
  }
  
}

