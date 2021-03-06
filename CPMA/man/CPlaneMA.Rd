% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CPlaneMA.R
\name{CPlaneMA}
\alias{CPlaneMA}
\title{Change Plane Model Average method for subgroup identification.}
\usage{
CPlaneMA(
  data.tr,
  data.te = NULL,
  yind,
  xind,
  zind,
  ini.theta = matrix(1, 1, length(zind)),
  tol = 0.001,
  K = 10,
  subm.vol = rep(1, length(xind) + 1)
)
}
\arguments{
\item{data.tr}{a matrix or data frame of the training dataset.}

\item{data.te}{a matrix or data frame of the testing dataset.}

\item{yind}{the column number or column name of the response variable.}

\item{xind}{a numeric or character vector, containing the column numbers or column names of the predictors (covariates).}

\item{zind}{a numeric or character vector, containing the column numbers or column names of the threshold variables.}

\item{ini.theta}{a numeric matrix, with each row containing the initial values of change plane parameters for one
submodel. Note that the initial values (the rows) should be arranged to match the order of submodels.
If there is no enough initial value entered, i.e. the row number of ini.theta is smaller than the
number of submodels, the provided initial values will be repeated for use. For better performance
of this function, it is highly recommended that the user provide an input based on preknowledge
or a reasonable guess of the change plane structure. Default is ini.theta=matrix(1, 1, length(zind)).}

\item{tol}{error tolerance in parameter estimation and optimization. Default is tol=1e-3.}

\item{K}{maximum time of iterations in parameter estimation and optimization. Dault is K=10.}

\item{subm.vol}{a numeric vector, containing the number of predictors (starting from the constant term) whose
varying covariate effect should be considered in the corresponding submodel.
Default is subm.vol=rep(1, length(xind)+1).}
}
\value{
A list consisting of the following components:
\item{train.res}{estimation results for the training dataset, including estimated response value (Yhat), subgroup id (subgroup), mean squared error (mse), estimated change plane parameter (theta), thresholds on the linear combination of threshold variables (threshold) and the segment regression coefficients (coefficient).}
\item{test.res}{estimation results for the testing dataset, including estimated response value (Yhat), subgroup id (subgroup) and mean squared error (mse).}
\item{submodel.res}{estimation results for each of the submodels.}
}
\description{
While CPlane() assumes that subgroups are charaterised by parallel change planes, this function
admits change planes that are no necessarily parallel and yields multiple vectors of change plane
parameters through model averaging.
}
\details{
The methodology of this function consists of two levels. In the first level, a number of individual
change plane regression submodels are fitted to model the varying covariate effect of some of the
given predictors (the constant term may also be included) in the training dataset. In the second
level, their model averaging weights are estimated so that a weighted ensemble of these submodels
can be used to further approximate the true model. Since the change plane parameters yielded by
different submodels are typically not the same, the averaged full model admits change planes that
are no necessarily parallel. The structure of submodels can be specified by the user through the
input parameter subm.vol.
}
\examples{
Theta <- matrix(0, nrow = 6, ncol = 5)
Theta[1,] <- rep(1, 5)
Theta[2,] <- c(sqrt(0.5), -sqrt(0.5), 0, 0, 0)
Theta[3,] <- c(0.75, 0, -0.5, 0, sqrt(1-(0.75)^2-(-0.5)^2))
Theta[4,] <- c(sqrt(0.5), -sqrt(0.5), 0, 0, 0)
Theta[5,] <- c(0.75, 0, -0.5, 0, sqrt(1-(0.75)^2-(-0.5)^2))
Theta[6,] <- rep(1, 5)
out <- CPlaneMA(data.tr=CPdata[1:350,], data.te=CPdata[-(1:350),], yind=1, xind=2:6, zind=2:6,
ini.theta=Theta, subm.vol=rep(1, 6))
}
\references{
Li, J., Y. Li, B. Jin, and M. R. Kosorok (2021). Multithreshold change plane model: estimation theory and applications in subgroup identification. Statistics in Medicine 40(15), 3440-3459.
}
