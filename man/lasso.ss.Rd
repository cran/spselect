\name{lasso.ss}
\alias{lasso.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial scale lasso
}
\description{
This function fits a spatial scale (SS) lasso model.
}
\usage{
lasso.ss(y, X, ss, a.lst, S.v, C.v, col.plot, verbose=TRUE, plot=TRUE)
}

\arguments{
  \item{y}{
  A numeric response vector
}
  \item{X}{
  A data frame of numeric variables
}
  \item{ss}{
  A vector of names to identify the different levels of covariates available as potential candidates for model input
}
  \item{a.lst}{
  A list of identity matrices, where each column indicates a particular level or spatial scale for a specified covariate (e.g., ss1_x2)
}
  \item{S.v}{
  A vector of positive integers, where each number denotes the number of spatial scales associated with a particular covariate
}
  \item{C.v}{
  A vector, where all values are initialized to 0
}
  \item{col.plot}{
  A vector of colors (corresponding to each SS) used in the coefficient path plot
}
  \item{verbose}{
  If TRUE, details are printed as the algorithm progresses
}
  \item{plot}{
  If TRUE, a coefficient path plot is generated
}
}
\details{
This function estimates coefficients using the SS lasso modeling approach.
The function also provides summary details and plots a coefficient path plot.
}
\value{
 A list with the following items:

  \item{beta }{Regression coefficient estimates from all set of model solutions}
  \item{beta.aic }{Regression coefficient estimates from final model}
  \item{ind.v }{Vector of indices to denote the corresponding columns of X associated with each active predictor }
  \item{aic.v }{Vector of Akaike information criterion (AIC) values}
  \item{stack.ss }{Vector of indices to indicate the level at which each covariate enters the model}

}
\references{
Grant LP, Gennings C, Wheeler, DC. (2015). Selecting spatial scale of covariates in regression models of environmental exposures.
Cancer Informatics, 14(S2), 81-96. doi: 10.4137/CIN.S17302
}
\author{
Lauren Grant, David Wheeler
}


\examples{
data(y)
data(X)

names.X <- colnames(X)

ss <- c("ind", "ss1", "ss2")

a.lst <- list(NULL)
a.lst[[1]] <- 1
dim(a.lst[[1]]) <- c(1,1)
dimnames(a.lst[[1]]) <- list(NULL, names.X[1])

a.lst[[2]] <- diag(2)
dimnames(a.lst[[2]]) <- list(NULL, names.X[c(2,3)])

a.lst[[3]] <- diag(2)
dimnames(a.lst[[3]]) <- list(NULL, names.X[c(4,5)])

S.v <- c(1,2,2)
C.v <- rep(0,length(a.lst))

mod_lasso.ss <- lasso.ss(y, X, ss, a.lst, S.v, C.v, c("black", "red", "green"))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

