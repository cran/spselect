\name{stepwise.ss}
\alias{stepwise.ss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial scale forward stepwise regression
}
\description{
This function fits a spatial scale (SS) forward stepwise regression model.
}
\usage{
stepwise.ss(y, X.3D, y.name, ss, epsilon, verbose=TRUE)
}

\arguments{
  \item{y}{
  A numeric response vector
}
  \item{X.3D}{
  A 3-D or stacked array of numeric variables, where each stack represents a particular level of covariates (i.e., individual- and area-level
  variables at more than one spatial scale). In cases where values are only present for a covariate at certain levels, that covariate is assigned
  missing values at all other levels.                                
}
  \item{y.name}{
  A name for y
}
  \item{ss}{
  A vector of names to identify the different levels of covariates available as potential candidates for model input
}
  \item{epsilon}{
  A positive value used as a stopping criterion when there is inadequate improvement in the model's performance. The algorithm stops if the
  difference in the Akaike information criterion (AIC) between the current model and the proposed model is less than epsilon.  
}
  \item{verbose}{
  If TRUE, details are printed as the algorithm progresses
}
}
\details{
This function estimates coefficients using the SS forward stepwise regression approach.
The function also estimates the model fit and provides summary details.
}
\value{
 A list with the following items:

  \item{beta.final }{Regression coefficient estimates from final model}
  \item{aic.final }{AIC for final model}
  \item{summary.final }{Summary output of final model}
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
data(X.3D)
y.name <- "y"
ss <- c("ind", "ss1", "ss2")
mod_forward.step.ss_1 <- stepwise.ss(y, X.3D, y.name, ss, 1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

