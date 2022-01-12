
#' @title Change Plane method for subgroup identification.
#' @description
#' Idenfiy subgroups marked by multiple parallel change planes.
#' @details
#' This function identifies subgroups and makes predictions by fitting a multithreshold change plane regression
#' model to the given training dataset. Subgroups are thus charaterised by a linear combination of the specified
#' threshold variables. The change plane parameters, thresholds and regression coefficients are estimated by a
#' two stage approach proposed by Li et al (2018). Note that this function should only be used when there are
#' multiple threshold variables so that a linear combination of them makes sense, i.e. length(zind) > 1. Otherwise,
#' one should turn to the CPoint() function instead.
#' @param data.tr a matrix or data frame of the training dataset.
#' @param data.te a matrix or data frame of the testing dataset.
#' @param yind the column number or column name of the response variable.
#' @param xind a numeric or character vector, containing the column numbers or column names of the predictors (covariates).
#' @param zind a numeric or character vector, containing the column numbers or column names of the threshold variables.
#' @param ini.theta
#' a numeric vector, containing the initial values of change plane parameters (the linear combination coefficients
#' of the threshold variables). For better performance of this function, it is highly recommended that the user
#' provide an input based on preknowledge or a reasonable guess of the change plane structure.
#' Default is ini.theta=rep(1, length(zind)).
#' @param tol error tolerance in parameter estimation and optimization. Default is tol=1e-3.
#' @param K maximum time of iterations in parameter estimation and optimization. Dault is K=10.
#'
#' @return A list consisting of the following components:
#' \item{train.res}{estimation results for the training dataset, including estimated response value (Yhat), subgroup id (subgroup), mean squared error (mse), estimated change plane parameter (theta), thresholds on the linear combination of threshold variables (threshold) and the segment regression coefficients (coefficient).}
#' \item{test.res}{estimation results for the testing dataset, including estimated response value (Yhat), subgroup id (subgroup) and mean squared error (mse).}
#' @references Li, J., Y. Li, B. Jin, and M. R. Kosorok (2021). Multithreshold change plane model: estimation theory and applications in subgroup identification. Statistics in Medicine 40(15), 3440-3459.
#' @importFrom stats coef lm pnorm
#' @importFrom ncvreg ncvreg
#' @importFrom BB BBoptim
#' @importFrom nloptr nloptr
#' @import grpreg
#' @export
#'
#' @examples
#' out <- CPlane(data.tr=CPdata[1:350,], data.te=CPdata[-(1:350),], yind = 1, xind = 2:6, zind = 2:6,
#' ini.theta = c(sqrt(0.5), -sqrt(0.5), 0, 0, 0))
CPlane <- function(data.tr, data.te=NULL, yind, xind, zind,
                   ini.theta=rep(1, length(zind)), tol=1e-3, K=10) {

  data.tr <- as.matrix(data.tr)
  Y.tr <- data.tr[, yind]
  X.tr <- data.tr[, xind]
  Z.tr <- data.tr[, zind]

  n <- dim(X.tr)[1]
  p <- dim(X.tr)[2]
  q <- length(zind)

  ini.theta <- ini.theta/norm(ini.theta, '2')

  obj.fit <- subm.fun(X = X.tr, Y = Y.tr, Z = Z.tr, l = 0, len = p+1,
                      ini.theta = ini.theta, tol = tol, K = K)
  train.est <- pointest(X = X.tr, Y = Y.tr, Z = Z.tr, l = 0, len = p+1,
                        s = obj.fit$s.est, obj = obj.fit$sm.obj)

  subgroup <- rep(1, n)
  if (obj.fit$s.est >= 1){
    for (j in 1:obj.fit$s.est) {
      id <-  which( Z.tr %*% obj.fit$sm.obj$theta.hat > obj.fit$sm.obj$a.hat[j])
      subgroup[id] <- subgroup[id]+1
    }
  }
  delta <- matrix(obj.fit$sm.obj$delta.hat, nrow = obj.fit$s.est, ncol = p+1)
  coeff <- list(beta=obj.fit$sm.obj$beta.hat, delta=delta)
  train.res <- list(Yhat=train.est$Yhat, subgroup=as.factor(subgroup), mse=train.est$mse,
                    theta=obj.fit$sm.obj$theta.hat, threshold=obj.fit$sm.obj$a.hat, coefficient=coeff)

  test.res <- NULL
  if (!is.null(data.te)) {

    data.te <- as.matrix(data.te)
    Y.te <- data.te[, yind]
    X.te <- data.te[, xind]
    Z.te <- data.te[, zind]

    test.est <- pointest(X = X.te, Y = Y.te, Z = Z.te, l = 0, len = p+1,
                         s = obj.fit$s.est, obj = obj.fit$sm.obj)
    subgroup.te <- rep(1, length(Y.te))
    if (obj.fit$s.est >= 1){
      for (j in 1:obj.fit$s.est) {
        id <-  which( Z.te %*% obj.fit$sm.obj$theta.hat > obj.fit$sm.obj$a.hat[j])
        subgroup.te[id] <- subgroup.te[id]+1
      }
    }

    test.res <- list(Yhat=test.est$Yhat, subgroup=as.factor(subgroup.te), mse=test.est$mse)

  }

  return(list(train.res=train.res, test.res=test.res))

}

