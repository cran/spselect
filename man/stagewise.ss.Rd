\name{stagewise.ss}
\alias{stagewise.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial scale incremental forward stagewise regression
}
\description{
This function fits a spatial scale (SS) incremental forward stagewise regression model.
}
\usage{
stagewise.ss(y, X, X.3D, ss, increment, tolerance, col.plot, verbose=TRUE, plot=TRUE)
}

\arguments{
  \item{y}{
  A numeric response vector
}
  \item{X}{
  A data frame of numeric variables
}
  \item{X.3D}{
  A 3-D or stacked array of numeric variables, where each stack represents a particular level of covariates (i.e., individual- and area-level
  variables at more than one spatial scale). In cases where values are only present for a covariate at certain levels, that covariate is assigned
  missing values at all other levels.
}
  \item{ss}{
  A vector of names to identify the different levels of covariates available as potential candidates for model input
}
  \item{increment}{
  A positive step size
}
  \item{tolerance}{
  A small, positive value used as a stopping criterion when none of the predictors are correlated with the residuals. The algorithm stops if the
  overall maximum correlation is less than a specified tolerance.
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
This function estimates coefficients using the SS forward stagewise regression approach.
The function also provides summary details and plots a coefficient path plot.
}
\value{
 A list with the following items:

  \item{beta.final }{Regression coefficient estimates from final model}
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
data(X.3D)
ss <- c("ind", "ss1", "ss2")
mod_forward.stage.ss_0.1 <- stagewise.ss(y, X ,X.3D, ss, 0.1, 0.1, c("black", "red", "green"))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}