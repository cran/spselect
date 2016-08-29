\name{spselect-package}
\alias{spselect-package}
\alias{spselect}
\docType{package}
\title{
Selecting spatial scale of area-level covariates in regression models
}
\description{
Fits spatial scale (SS) forward stepwise regression, SS incremental forward stagewise regression, SS least angle regression (LARS), and SS lasso
models.  All area-level covariates are considered at all available scales to enter a model, but the SS algorithms are constrained to select each area-level covariate at a single spatial scale.
}
\details{
\tabular{ll}{
Package: \tab spselect\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2016-08-29\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}

}
\author{
Lauren Grant, David Wheeler

Maintainer: Lauren Grant <pacele@vcu.edu>

}
\references{
Grant LP, Gennings C, Wheeler, DC. (2015). Selecting spatial scale of covariates in regression models of environmental exposures.
Cancer Informatics, 14(S2), 81-96. doi: 10.4137/CIN.S17302
}

\keyword{ package }

\examples{
data(y)
data(X.3D)
y.name <- "y"
ss <- c("ind", "ss1", "ss2")
mod_forward.step.ss_1 <- stepwise.ss(y, X.3D, y.name, ss, 1)
}
