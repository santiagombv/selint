#' Generation of modular eigenvalue distribution.
#' 
#' MODeigenval generates a vector of numeric values, to be used as the 
#' eigenvalues of a correlation matrix with modular structure.  
#' 
#' @param dim number of eigenvalues.
#' @param k nunber of modules.
#' @param n
#' @param r
#' @param c
#' @param sigma the standard deviation of a normally distributed error term, which 
#' is added  to the eigenvalues.
#' @param tol if one or more of the obtained eigenvalues is negative or smaller 
#' than tol, they eare replaced by the \code{tol} value and all values are rescaled 
#' to add the value determined by \code{dim}. This procedure (eigenvalue bending)
#' assures that matrices with these eigenvalues were positive. 
#' 
#' @details In a broken-stick distribution of eigenvalues, the first (leading) 
#' eigenvalue accounts for most of the variance and the remaining eigenvalues are 
#' equal (in the extreme case where all *r* = 1, the leading eigenvalue equals 
#' the number of traits and the remaining eigenvalues are zero). This kind of 
#' eigenvalue distribution corresponds to a homogeneous matrix were all *r* are 
#' equal (Pavlicev et al. 2009) and INT is $r^{2}$. To introduce randomness around the 
#' expected homogeneous matrix, we added a normally distributed error term to the 
#' eigenvalues (with zero mean and a given standard deviation), rescaled eigenvalues
#' to ensure they totalled *N* and eliminated negative and zero eigenvalues, if they
#' existed. The variability of the correlation matrices is controlled by the standard
#' deviation of the error term. This procedure avoid to simulate uniform matrices, 
#' however, is worth to mention that when the desired integration is high, adding large
#' error terms leads to a substantially reduction in integration.
#'  
#' @return a vector of positive numeric values that would be used as eigenvalues
#' for matrix simulation
#'  
#' @references
#' Pavlicev, M., Cheverud, J.M., and Wagner, G.P. 2009. Measuring morphological 
#' integration using eigenvalue variance. Evolutionary Biology 36(1):157-170.
#'  
#' @examples
#' BSeigenval(INT = 50, dim = 8, sigma = 0.1)
#' 
#' # testing a wide range of INT
#' sim.intBS <- numeric(1000)
#' int <- seq(1, 99, length.out=1000)
#' for(i in 1:1000){
#'   dim <- 10
#'   Z <-BSeigenval(INT = int[i], dim = dim, sigma = 0.05)
#'   sim.intBS[i] <- ((sum((Z-mean(Z))^2)/dim))/(dim-1)*100
#' }
#' plot(int, sim.intBS)
#' abline(0,1, col ="red")
#'  
#' @export
#' 
MODeigenval <- function(dim, k, n, r, c, sigma, tol = 0.0001){
  values <- list(k)
  for(i in 1:k){
    L1 <- 1 + (n[i] - 1)*r[i]
    Lk <- rep((1 - r[i]), n[i]-1)
    values[[i]] <- c(L1, Lk)
  }
  d <- sort(unlist(values), decreasing = T) 
  c <- d[2] * c
  d[1] <- d[1] + c
  d[2] <- d[2] - c
  d <- d + rnorm(dim, mean = 0, sd = sigma)
  if(min(d)<tol){d <- d+abs(min(d))+tol}
  d <- sort(d*dim/sum(d), decreasing = TRUE)
  d
}
