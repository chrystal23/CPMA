
#' @title Change Point method for subgroup identification.
#' @description 
#' Subgroup identification by change point detection on a (single) specified threshold variable Z.
#' @details 
#' This function identifies subgroups and makes predictions by fitting a change point (threshold regression) 
#' model to the given training dataset. The change points (thresholds) and regression coefficients 
#' are estimated by a two stage multiple change-points detection (TSMCD) method proposed by Li and Jin (2018).
#' @param data.tr a matrix or data frame of the training dataset.
#' @param data.te a matrix or data frame of the testing dataset.
#' @param yind the column number or column name of the response variable.
#' @param xind a numeric or character vector, containing the column numbers or column names of the predictors (covariates).
#' @param zind the column number or column name of the threshold variable.
#' @param c a numeric vector to determine the initial segment length in the splitting stage of TSMCD, i.e. the tentative choices of initial segment length m = c * sqrt(n), where n is the sample size of the training dataset. Default is c = seq(0.5, 1.5, 0.1)
#' 
#' @return A list consisting of the following components:
#' \item{train.res}{estimation results for the training dataset, including estimated response value (Yhat), subgroup id (subgroup), mean squared error (mse), detected change points (threshold) and the segment regression coefficients (coefficient).}
#' \item{test.res}{estimation results for the testing dataset, including estimated response value (Yhat), subgroup id (subgroup) and mean squared error (mse).}
#' @references Li, J. and B. Jin (2018). Multi-threshold accelerated failure time model. The Annals of Statistics 46(6A), 2657-2682.
#' @importFrom plus plus
#' @export
#'
#' @examples out <- CPoint(data.tr=CPdata[1:350,], data.te=CPdata[-(1:350),], yind=1, xind=2:6, zind=2)
CPoint <- function(data.tr, data.te=NULL, yind, xind, zind,
                   c=seq(0.5,1.5,0.1)) {

  data.tr <- as.matrix(data.tr)
  Z.tr <- data.tr[, zind]

  ord <- order(Z.tr)
  data.tr <- data.tr[ord,] # re-sort by the thresholding variable

  Y.tr <- data.tr[, yind]
  X.tr <- data.tr[, xind]
  n <- dim(X.tr)[1]
  p <- dim(X.tr)[2]

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
  thre <- (data.tr[, zind])[tsmcd$cp]
  coeff <- matrix(tsmcd$coef, nrow = length(thre)+1, ncol = p+1, byrow = T)
  subgroup <- rep(1, n)
  if (length(thre) >= 1) {
    for (i in 1:length(thre)) {
      id <- which(Z.tr>thre[i])
      subgroup[id] <- subgroup[id]+1
    }
  }

  train.res <- list(Yhat=Yhat, subgroup=as.factor(subgroup), mse=tsmcd$sigmaep0, threshold=thre, coefficient=coeff)

  test.res <- NULL

  if (!is.null(data.te)) {

    data.te <- as.matrix(data.te)
    Y.te <- data.te[, yind]
    X.te <- data.te[, xind]
    Z.te <- data.te[, zind]

    tsm.te=TSMCD.pred(Y.te, X.te, Z.te, thre, coeff)

    subgroup.te <- rep(1, length(Y.te))
    if (length(thre) >= 1) {
      for (i in 1:length(thre)) {
        id <- which(Z.te>thre[i])
        subgroup.te[id] <- subgroup.te[id]+1
      }
    }

    test.res <- list(Yhat=tsm.te$Yhat, subgroup=as.factor(subgroup.te), mse=tsm.te$sigmaep0)

  }

  return(list(train.res=train.res, test.res=test.res))

}
