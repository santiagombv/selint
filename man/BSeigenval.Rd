% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/BSeigenval.R
\name{BSeigenval}
\alias{BSeigenval}
\title{Generation of broken-stick eigenvalue distribution}
\usage{
BSeigenval(INT, dim, sigma, tol = 1e-04)
}
\arguments{
\item{INT}{the Wagner-Cheverud phenotypic integration index (in a 0-100 scale).}

\item{dim}{the number of eigenvalues.}

\item{sigma}{the standard deviation of a normally distributed error term, which
is added  to the eigenvalues.}

\item{tol}{if one or more of the obtained eigenvalues is negative or smaller
than tol, they eare replaced by the \code{tol} value and all eigenvalues are rescaled
to add the value determined by \code{dim}. This procedure assures that matrices with
these eigenvalues were positive.}
}
\value{
a vector of positive numeric values that would be used as eigenvalues
for matrix simulation.
}
\description{
BSeigenval generates a vector of numeric values following a broken-stick
distribution, to be used as the eigenvalues of a correlation matrix.
This allows to set the integration value of the resulting matrix.
}
\details{
In a broken-stick distribution of eigenvalues, the first (leading)
eigenvalue accounts for most of the variance and the remaining eigenvalues are
equal (in the extreme case where all r = 1, the leading eigenvalue equals
the number of traits and the remaining eigenvalues are zero). This kind of
eigenvalue distribution corresponds to a homogeneous matrix were all r are
equal (Pavlicev et al. 2009) and INT is \eqn{r^{2}}. To introduce randomness around the
expected homogeneous matrix, we added a normally distributed error term to the
eigenvalues (with zero mean and a given standard deviation), rescaled eigenvalues
to ensure they totalled N and eliminated negative and zero eigenvalues, if they
existed. The variability of the correlation matrices is controlled by the standard
deviation of the error term. This procedure avoid to simulate uniform matrices,
however, is worth to mention that when the desired integration is high, adding large
error terms leads to a substantially reduction in integration.
}
\examples{
BSeigenval(INT = 50, dim = 8, sigma = 0.1)

# testing a wide range of INT
sim.intBS <- numeric(1000)
int <- seq(1, 99, length.out=1000)
for(i in 1:1000){
  dim <- 10
  Z <-BSeigenval(INT = int[i], dim = dim, sigma = 0.05)
  sim.intBS[i] <- ((sum((Z-mean(Z))^2)/dim))/(dim-1)*100
}
plot(int, sim.intBS)
abline(0,1, col ="red")
}
\author{
Santiago Benitez-Vieyra
}
\references{
Pavlicev, M., Cheverud, J.M., and Wagner, G.P. 2009. Measuring morphological
integration using eigenvalue variance. Evolutionary Biology 36(1):157-170.
}
\seealso{
\code{\link{MODeigenval}}
}

