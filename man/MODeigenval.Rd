% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/MODeigenval.R
\name{MODeigenval}
\alias{MODeigenval}
\title{Generation of modular eigenvalue distribution.}
\usage{
MODeigenval(dim, n, r, c, sigma, tol = 1e-04)
}
\arguments{
\item{dim}{number of eigenvalues.}

\item{n}{a numeric vector of length 2 indicating the number of variables to be
assigned to each module}

\item{r}{a numeric vector of length 2 indicating which proportion of the
variance within a module have to be assigned to the leading eigenvalue. If this
proportion approaches 1, the leading eigenvalue will account for all the variance.
On the other hand, when this proportion approaches \eqn{1/n_i}, variance will be equally
distributed among eigenvalues.}

\item{c}{a constant to be subtracted from one leading eigenvalue and added to the
other. When c equals \eqn{\lambda_{1,2}} the resulting eigenvalue structure is the typical
from a homogeneous correlation matrix, and when c equals zero the resulting eigenvalue
structure from a matrix with homogeneous and independent modules. c ranges from zero to
\eqn{n_{2}*r_{2}}.}

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
MODeigenval generates a vector of numeric values, to be used as the
eigenvalues of a correlation matrix with two modules.
}
\details{
Modular correlation matrices with completely independent, non-overlapping
modules (i.e. when among-module correlations, \eqn{r_{a}}, are zero) and with homogeneous
modules (i.e. with homogeneous within-module correlations, \eqn{r_{w}}), behave as the sum
of isolated homogeneous matrices (Pavlicev et al. 2009). Each module has one leading eigenvalue,
which value depends on the within-module correlation and on the number of traits included in the
module,
\deqn{\lambda_{1,i} = 1 + (N_{i} - 1)r_{wi}}
while the remaining eigenvalues equal \eqn{1 - r_{wi}}. Thus, the whole matrix has as many
leading eigenvalues as modules. If matrices were simulated from this kind of eigenvalue
distribution, the among-module correlations will be zero. To obtain matrices with two
modules and among-module correlations different from zero, we introduced c, a constant
to be subtracted from one leading eigenvalue and added to the other leading eigenvalue.
Supposing that c is subtracted from \eqn{\lambda_{1,1}} and added to \eqn{\lambda_{1,2}}, its
easy to see that when c equals \eqn{\lambda_{1,2}} the resulting eigenvalue structure is the
typical from a homogeneous correlation matrix, and when c equals zero the resulting
eigenvalue structure from a matrix with homogeneous and independent modules. Intermediate
values of c may result in different outcomes, depending on the chosen value: matrices with
among-module correlations greater than zero or matrices with modules of different sizes.
}
\examples{
MODeigenval(dim = 6, c = 0.1, sigma = 0.01, n = c(3,3), r =c(0.8, 0.7))
}
\author{
Santiago Benitez-Vieyra
}
\references{
Pavlicev, M., Cheverud, J.M., and Wagner, G.P. 2009. Measuring morphological
integration using eigenvalue variance. Evolutionary Biology 36(1):157-170.
}
\seealso{
\code{\link{BSeigenval}}
}

