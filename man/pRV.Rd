\name{pRV}
\alias{pRV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pRV(data, vars1, vars2, method = c("cov", "cor"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{vars1}{
%%     ~~Describe \code{vars1} here~~
}
  \item{vars2}{
%%     ~~Describe \code{vars2} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, vars1, vars2, method = c("cov", "cor"), ...) 
{
    if (method == "cov") 
        mat <- cov(data[, c(vars1, vars2)], ...)
    else mat <- cor(data[, c(vars1, vars2)], ...)
    N <- nrow(data)
    Trv <- pairRV(mat = mat, vars1 = vars1, vars2 = vars2)
    Lrv <- numeric(N)
    for (i in 1:N) {
        if (method == "cov") 
            Lrv[i] <- pairRV(mat = cov(data[-i, ], ...), vars1 = vars1, 
                vars2 = vars2)
        else Lrv[i] <- pairRV(mat = cor(data[-i, ], ...), vars1 = vars1, 
            vars2 = vars2)
    }
    Drv <- N * Trv - (N - 1) * Lrv
    Drv
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
