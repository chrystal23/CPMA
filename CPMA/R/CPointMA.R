
#' @title Change Point Model Average method for subgroup identification.
#' @description 
#' While CPoint() identifies subgroups based on a single threshold variable, this function 
#' detects change points each threshold variable separately and then assembles the subgrouping 
#' results through model averaging.
#' @details 
#' The methodology of this function consists of two levels. An individual change point (threshold regression) 
#' submodel is fitted based on each of the given threshold variables in the training dataset. 
#' Then their model averaging weights are estimated so that a weighted ensemble of these submodels 
#' can be used to further approximate the true model. Note that this function should only be used
#' when there are multiple threshold variables, i.e. length(zind) > 1. Otherwise, one should turn to
#' the CPoint() function instead.
#' @param data.tr a matrix or data frame of the training dataset.
#' @param data.te a matrix or data frame of the testing dataset.
#' @param yind the column number or column name of the response variable.
#' @param xind a numeric or character vector, containing the column numbers or column names of the predictors (covariates).
#' @param zind a numeric or character vector, containing the column numbers or column names of the threshold variables.
#' @param c a numeric vector to determine the initial segment length in the splitting stage of TSMCD for submodel estimation, i.e. the tentative choices of initial segment length m = c * sqrt(n), where n is the sample size of the training dataset. Default is c = seq(0.5, 1.5, 0.1)
#' @param penalty the penalty to be used in the model average step, including 'SCAD', 'MCP' and 'LASSO'. Default is 'SCAD'.
#'
#' @return A list consisting of the following components:
#' \item{train.res}{estimation results for the training dataset, including estimated response value (Yhat), subgroup id (subgroup), mean squared error (mse), detected change points (threshold) and the segment regression coefficients (coefficient).}
#' \item{test.res}{estimation results for the testing dataset, including estimated response value (Yhat), subgroup id (subgroup) and mean squared error (mse).}
#' \item{submodel.res}{estimation results for each of the submodels.}
#' @references Li, J. and B. Jin (2018). Multi-threshold accelerated failure time model. The Annals of Statistics 46(6A), 2657-2682.
#' @importFrom plus plus
#' @export
#'
#' @examples out <- CPointMA(data.tr=CPdata[1:350,], data.te=CPdata[-(1:350),], yind=1, xind=2:6, zind=2:6)
CPointMA <- function(data.tr, data.te=NULL, yind, xind, zind, 
                   c=seq(0.5,1.5,0.1), penalty=c('SCAD', 'MCP', 'LASSO')) {
  
  data.tr <- as.matrix(data.tr)
  n <- dim(data.tr)[1]
  p <- length(xind)
  q <- length(zind)
  
  result=NULL
  V=NULL
  V.te=NULL
  
  for (l in 1:q) {
  
    Z.tr <- data.tr[, zind[l]]
    ord <- order(Z.tr)
    XY.tr <- data.tr[ord,] # re-sort by the thresholding variable
    Y.tr <- XY.tr[, yind]
    X.tr <- XY.tr[, xind]
    
    m <- ceiling(c*sqrt(n))
    c <- c[which(m>p+1)]
    bicy <- c
    tsmc <- NULL
    
    for(i in 1:length(c)) {
      tsm <- TSMCD(Y.tr,X.tr,c[i])
      bicy[i] <- log(n)*((length(tsm$cp)+1)*(p+1))+n*log(tsm$sigmaep0) 
      tsmc[[i]] <- tsm
    }
    
    tsmcd <- tsmc[[which(bicy==min(bicy))[1]]] # choose the optimal results by BIC
    
    Yhat <- tsmcd$Yhat
    Yhat[ord] <- tsmcd$Yhat
    thre <- (XY.tr[, zind[l]])[tsmcd$cp]
    coeff <- matrix(tsmcd$coef, nrow = length(thre)+1, ncol = p+1, byrow = T)
    subgroup <- rep(1, n)
    if (length(thre) >= 1) {
      for (i in 1:length(thre)) {
        id <- which(Z.tr>thre[i])
        subgroup[id] <- subgroup[id]+1
      }
    }
    
    train.res <- list(Yhat=Yhat, subgroup=as.factor(subgroup), threshold=thre, coefficient=coeff, mse=tsmcd$sigmaep0)
    V <- cbind(V, Yhat)
    
    test.res <- NULL
    
    if (!is.null(data.te)) {
      data.te <- as.matrix(data.te)
      Y.te <- data.te[, yind]
      X.te <- data.te[, xind]
      Z.te <- data.te[, zind[l]]
      tsm.te=TSMCD.pred(Y.te, X.te, Z.te, thre, coeff)
      subgroup.te <- rep(1, length(Y.te))
      if (length(thre) >= 1) {
        for (i in 1:length(thre)) {
          id <- which(Z.te>thre[i])
          subgroup.te[id] <- subgroup.te[id]+1
        }
      }
      test.res <- list(Yhat=tsm.te$Yhat, subgroup=as.factor(subgroup.te), mse=tsm.te$sigmaep0)
      V.te <- cbind(V.te, tsm.te$Yhat)
    }
    
    result[[l]] <- list(train.res=train.res, test.res=test.res)
  
  }
  
  YY <- data.tr[, yind]
  object <- plus(V,YY,method = 'scad', gamma=2.4, intercept = F, normalize = F, eps=1e-30)
  
  bic=log(n)*object$dim+n*log(as.vector((1-object$r.square)*sum(YY^2))/length(YY))
  beta.min=c(apply(object$beta,1,min))
  iid=which(beta.min>=0)
  bic[-iid]=Inf
  step.bic=which.min(bic)  
  
  maver.weight <- coef(object,lam=object$lam[step.bic])
  
  Yhat <- V %*% maver.weight
  sigmaep <- mean((YY-Yhat)^2)
  
  subgr=rep('',n)
  for (l in 1:q) {
    if (maver.weight[l]>0) {
      subgr=paste(subgr, result[[l]]$train.res$subgroup, sep = '')
    }
  }
  coeff.hat=NULL
  a.hat=NULL
  for (l in 1:q){
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
    YY.te <- data.te[, yind]
    Yhat.te <- V.te %*% maver.weight
    sigmaep.te <- mean((YY.te-Yhat.te)^2)
    subgr.te=rep('', dim(data.te)[1])
    for (l in 1:q) {
      if (maver.weight[l]>0) {
        subgr.te=paste(subgr.te, result[[l]]$test.res$subgroup, sep = '')
      }
    }
  }
  
  test.res <- list(Yhat=Yhat.te, subgroup=as.factor(subgr.te), mse=sigmaep.te)
  
  return(list(train.res=train.res, test.res=test.res, submodel.res=result))
  
}
