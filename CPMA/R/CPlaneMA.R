
#' @title Change Plane Model Average method for subgroup identification.
#' @description
#' While CPlane() assumes that subgroups are charaterised by parallel change planes, this function
#' admits change planes that are no necessarily parallel and yields multiple vectors of change plane
#' parameters through model averaging.
#' @details
#' The methodology of this function consists of two levels. In the first level, a number of individual
#' change plane regression submodels are fitted to model the varying covariate effect of some of the
#' given predictors (the constant term may also be included) in the training dataset. In the second
#' level, their model averaging weights are estimated so that a weighted ensemble of these submodels
#' can be used to further approximate the true model. Since the change plane parameters yielded by
#' different submodels are typically not the same, the averaged full model admits change planes that
#' are no necessarily parallel. The structure of submodels can be specified by the user through the
#' input parameter subm.vol.
#' @param data.tr a matrix or data frame of the training dataset.
#' @param data.te a matrix or data frame of the testing dataset.
#' @param yind the column number or column name of the response variable.
#' @param xind a numeric or character vector, containing the column numbers or column names of the predictors (covariates).
#' @param zind a numeric or character vector, containing the column numbers or column names of the threshold variables.
#' @param ini.theta
#' a numeric matrix, with each row containing the initial values of change plane parameters for one
#' submodel. Note that the initial values (the rows) should be arranged to match the order of submodels.
#' If there is no enough initial value entered, i.e. the row number of ini.theta is smaller than the
#' number of submodels, the provided initial values will be repeated for use. For better performance
#' of this function, it is highly recommended that the user provide an input based on preknowledge
#' or a reasonable guess of the change plane structure. Default is ini.theta=matrix(1, 1, length(zind)).
#' @param tol error tolerance in parameter estimation and optimization. Default is tol=1e-3.
#' @param K maximum time of iterations in parameter estimation and optimization. Dault is K=10.
#' @param subm.vol
#' a numeric vector, containing the number of predictors (starting from the constant term) whose
#' varying covariate effect should be considered in the corresponding submodel.
#' Default is subm.vol=rep(1, length(xind)+1).
#'
#' @return A list consisting of the following components:
#' \item{train.res}{estimation results for the training dataset, including estimated response value (Yhat), subgroup id (subgroup), mean squared error (mse), estimated change plane parameter (theta), thresholds on the linear combination of threshold variables (threshold) and the segment regression coefficients (coefficient).}
#' \item{test.res}{estimation results for the testing dataset, including estimated response value (Yhat), subgroup id (subgroup) and mean squared error (mse).}
#' \item{submodel.res}{estimation results for each of the submodels.}
#' @references Li, J., Y. Li, B. Jin, and M. R. Kosorok (2021). Multithreshold change plane model: estimation theory and applications in subgroup identification. Statistics in Medicine 40(15), 3440-3459.
#' @importFrom plus plus
#' @importFrom stats coef lm pnorm
#' @importFrom ncvreg ncvreg
#' @importFrom BB BBoptim
#' @importFrom nloptr nloptr
#' @import grpreg
#' @export
#'
#' @examples
#' Theta <- matrix(0, nrow = 6, ncol = 5)
#' Theta[1,] <- rep(1, 5)
#' Theta[2,] <- c(sqrt(0.5), -sqrt(0.5), 0, 0, 0)
#' Theta[3,] <- c(0.75, 0, -0.5, 0, sqrt(1-(0.75)^2-(-0.5)^2))
#' Theta[4,] <- c(sqrt(0.5), -sqrt(0.5), 0, 0, 0)
#' Theta[5,] <- c(0.75, 0, -0.5, 0, sqrt(1-(0.75)^2-(-0.5)^2))
#' Theta[6,] <- rep(1, 5)
#' out <- CPlaneMA(data.tr=CPdata[1:350,], data.te=CPdata[-(1:350),], yind=1, xind=2:6, zind=2:6,
#' ini.theta=Theta, subm.vol=rep(1, 6))
CPlaneMA <- function(data.tr, data.te=NULL, yind, xind, zind,
                   ini.theta=matrix(1, 1, length(zind)), tol=1e-3, K=10,
                   subm.vol=rep(1, length(xind)+1)) {

  data.tr <- as.matrix(data.tr)
  Y.tr <- data.tr[, yind]
  X.tr <- data.tr[, xind]
  Z.tr <- data.tr[, zind]

  n <- dim(X.tr)[1]
  p <- dim(X.tr)[2]
  q <- length(zind)
  N <- length(subm.vol)

  ini.theta <- t(apply(ini.theta, 1, function(x) x/norm(x, '2')))
  if (dim(ini.theta)[1] < N) ini.theta <- matrix(rep(t(ini.theta), length.out=N*q),
                                                   nrow = N, ncol = q)
  varind <- c(0, cumsum(subm.vol[1:(N-1)]))

  V=NULL
  V.te=NULL
  result=NULL

  for (l in 1:N) {

    obj.fit <- subm.fun(X = X.tr, Y = Y.tr, Z = Z.tr, l = varind[l], len = subm.vol[l],
                        ini.theta = ini.theta[l,], tol = tol, K = K)
    train.est <- pointest(X = X.tr, Y = Y.tr, Z = Z.tr, l = varind[l], len = subm.vol[l],
                          s = obj.fit$s.est, obj = obj.fit$sm.obj)

    subgroup <- rep(1, n)
    if (obj.fit$s.est >= 1){
      for (j in 1:obj.fit$s.est) {
        id <-  which( Z.tr %*% obj.fit$sm.obj$theta.hat > obj.fit$sm.obj$a.hat[j] )
        subgroup[id] <- subgroup[id]+1
      }
      theta <- obj.fit$sm.obj$theta.hat
      threshold <- obj.fit$sm.obj$a.hat
      delta <- matrix(obj.fit$sm.obj$delta.hat, nrow = obj.fit$s.est, ncol = subm.vol[l])
      coeff <- list(beta=obj.fit$sm.obj$beta.hat, delta=delta)
    } else {
      theta <-'none'
      threshold <- 'none'
      coeff <- obj.fit$sm.obj$coefficients
    }

    train.res <- list(Yhat=train.est$Yhat, subgroup=as.factor(subgroup), mse=train.est$mse,
                      theta=theta, threshold=threshold, coefficient=coeff)
    V <- cbind(V, train.est$Yhat)

    test.res <- NULL
    if (!is.null(data.te)) {

      data.te <- as.matrix(data.te)
      Y.te <- data.te[, yind]
      X.te <- data.te[, xind]
      Z.te <- data.te[, zind]

      test.est <- pointest(X = X.te, Y = Y.te, Z = Z.te, l = varind[l], len = subm.vol[l],
                           s = obj.fit$s.est, obj = obj.fit$sm.obj)
      subgroup.te <- rep(1, length(Y.te))
      if (obj.fit$s.est >= 1){
        for (j in 1:obj.fit$s.est) {
          id <-  which( Z.te %*% obj.fit$sm.obj$theta.hat > obj.fit$sm.obj$a.hat[j] )
          subgroup.te[id] <- subgroup.te[id]+1
        }
      }

      test.res <- list(Yhat=test.est$Yhat, subgroup=as.factor(subgroup.te), mse=test.est$mse)
      V.te <- cbind(V.te, test.est$Yhat)

    }

    result[[l]] <- list(train.res=train.res, test.res=test.res)

    #print(l)

  }

  object <- plus(V, Y.tr, method = 'scad', gamma=3.1, intercept = F, normalize = F, eps = 1e-30)

  bic=log(n)*object$dim+n*log(as.vector((1-object$r.square)*sum(Y.tr^2))/length(Y.tr))
  beta.min=c(apply(object$beta,1,min))
  iid=which(beta.min>=0)
  bic[-iid]=Inf
  step.bic=which.min(bic)

  maver.weight <- coef(object,lam=object$lam[step.bic])

  Yhat <- V %*% maver.weight
  sigmaep <- mean((Y.tr-Yhat)^2)

  subgr=rep('',n)
  for (l in 1:N) {
    if (maver.weight[l]>0) {
      subgr=paste(subgr, result[[l]]$train.res$subgroup, sep = '')
    }
  }
  coeff.hat=NULL
  a.hat=NULL
  for (l in 1:N){
    if (maver.weight[l]>0) {
      coeff.hat[[paste('M',l,sep='')]] <- result[[l]]$train.res$coefficient
      a.hat[[paste('M',l,sep='')]] <- result[[l]]$train.res$threshold
    } else {
      coeff.hat[[paste('M',l,sep='')]] <- 'ineffective'
      a.hat[[paste('M',l,sep='')]] <- 'ineffective'
    }
  }

  train.res <- list(weight=maver.weight, Yhat=Yhat, subgroup=as.factor(subgr), threshold=a.hat, coefficient=coeff.hat, mse=sigmaep)

  test.res <- NULL

  if (!is.null(data.te)) {
    Y.te <- data.te[, yind]
    Yhat.te <- V.te %*% maver.weight
    sigmaep.te <- mean((Y.te-Yhat.te)^2)
    subgr.te=rep('', dim(data.te)[1])
    for (l in 1:N) {
      if (maver.weight[l]>0) {
        subgr.te=paste(subgr.te, result[[l]]$test.res$subgroup, sep = '')
      }
    }
  }

  test.res <- list(Yhat=Yhat.te, subgroup=as.factor(subgr.te), mse=sigmaep.te)

  return(list(train.res=train.res, test.res=test.res, submodel.res=result))

}
