\name{X.3D}
\alias{X.3D}
\docType{data}
\title{
Input data X.3D
}
\description{
Simulated input data (in stacked array format)
}
\usage{data(X.3D)}
\format{
  A 20x3x3 stacked array with the following stacks:

  \describe{
    \item{\code{ind}}{a numeric array containing an individual-level variable (x1)}
    \item{\code{ss1}}{a numeric array containing area-level variables (x2, x3) available at ss1}
    \item{\code{ss2}}{a numeric array containing area-level variables (x2, x3) available at ss2}
   }
}
\details{
The data consist of simulated variables, including an individual-level covariate (x1) and two area-level covariates (x2, x3) available at two
different spatial scales (ss1, ss2). The data are in the form of a 3-D or stacked array, where each stack represents a particular level of
covariates, including spatial scale. The first stack contains the individual-level variable; the second and third stacks contain the area-level
variables at the ss1 and ss2 levels, respectively. Note that in cases where values are only present for a covariate at certain levels, that
covariate is assigned missing values at all other levels.
}
\references{
Grant LP, Gennings C, Wheeler, DC. (2015). Selecting spatial scale of covariates in regression models of environmental exposures.
Cancer Informatics, 14(S2), 81-96. doi: 10.4137/CIN.S17302
} 
\examples{
data(X.3D)
}
\keyword{datasets}
