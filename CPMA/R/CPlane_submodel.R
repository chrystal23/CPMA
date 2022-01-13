

##
split.fun <- function(X, Y, l, len=1, c) {

  n=length(Y)
  p=dim(X)[2]
  m=ceiling(sqrt(n)*c)
  q=floor(n/m)

  ###########transform###################
  K_temp=matrix(0, nrow = q, ncol=q-1, byrow=TRUE)

  X_temp=cbind(1,X)

  Y_temp=c(Y)

  for(i in 2:q)
    K_temp[i,1:(i-1)]=rep(1,i-1)


  x=NULL
  y=NULL
  x[[1]]=as.matrix(X_temp[1:((n-(q-1)*m)),])
  y[[1]]=Y_temp[1:((n-(q-1)*m))]

  for(i in 2:q) {
    x[[i]] = as.matrix(X_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m)),])
    y[[i]] = Y_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
  }
  X_temp1 <- lapply(1:length(x), function(j, mat, list)
    kronecker(K_temp[j, , drop=F], x[[j]][,(l+1):(l+len),drop=F]), mat = K_temp, list = x)
  Xn=do.call("rbind",X_temp1)
  Xn=cbind(X_temp,Xn)

  Yn=NULL
  for(i in 1:q) {
    Yn=c(Yn,y[[i]])
  }

  group <- c(rep(1,p+1), kronecker(c(2:q),rep(1,len)))
  colnames(Xn) <- group
  fit <- grpreg::grpreg(Xn[,-1], Yn, group[-1], penalty = "grMCP")
  mcp.coef <- grpreg::select(fit,"BIC" )$beta
  mcp.lamb <- grpreg::select(fit,"BIC" )$lambda

  mcp.coef.l=c(mcp.coef[(l+1):(l+len)],mcp.coef[-(1:(p+1))])
  mcp.coef.v.m <- abs(matrix(c(mcp.coef.l), q, len, byrow = T))
  mcp.coef.m <- c(apply(mcp.coef.v.m, 1, max))
  mcp.cp <- which(mcp.coef.m!=0)

  if (length(mcp.cp) > 1) {
    for (i in 2:length(mcp.cp))
    {
      if (mcp.cp[i] - mcp.cp[i - 1] == 1)
        mcp.cp[i] = 0
    }
  }

  mcp.cp1 <- mcp.cp[mcp.cp > 1 & mcp.cp < q]

  ## consistent change location interval
  s.est <- length(mcp.cp1)
  if(s.est > 1){

    ci <- NULL
    ci[[1]] <- c( (n-(q - mcp.cp1[1]+2)*m) : (n-(q - mcp.cp1[1])*m) )
    numb <- n-(q-1)*m + (mcp.cp1[1] - 2)*m
    for(i in 2:s.est) {
      ci[[i]] = c( (n-(q - mcp.cp1[i]+2)*m) : (n-(q - mcp.cp1[i])*m) )
      numb <- c( numb, n-(q-1)*m+(mcp.cp1[i]-2)*m )
    }

    xx_temp <- NULL
    xx_temp[[1]] <- X_temp
    for(i in 1:s.est) {
      xx_temp[[i+1]] <- Xn[,which(colnames(Xn)==mcp.cp1[i])]
    }
    xx <- do.call(cbind, xx_temp)
    rss <- sum(lm(Yn~xx-1)$res^2)
    BIC <- log(n)*(p+1+s.est*len)+n*log(rss/n)


  } else if (s.est == 1) {

    ci <- c( (n-(q - mcp.cp1+2)*m) : (n-(q - mcp.cp1)*m) )
    numb <- n-(q-1)*m + (mcp.cp1 - 2)*m

    xx_temp <- NULL
    xx_temp[[1]] <- X_temp
    for(i in 1:s.est) {
      xx_temp[[i+1]] <- Xn[,which(colnames(Xn)==mcp.cp1[i])]
    }
    xx <- do.call(cbind, xx_temp)
    rss <- sum(lm(Yn~xx-1)$res^2)
    BIC <- log(n)*(p+1+s.est*len)+n*log(rss/n)

  }  else {

    rss <- sum(lm(Y~ X_temp -1 )$res^2)
    BIC <- log(n)*(p+1)+n*log(rss/n)
    mcp.cp1 = mcp.cp1
    ci = 1:n
    numb = c()

  }


  return(list(BIC = BIC, mcp.cp1 = mcp.cp1, numb = numb) )


}

##

obj.f <- function(X, Y, Z, l, len=1, a, theta, beta, delta, h) {

  X <- as.matrix(X)
  n <- dim(X)[1]
  X1 <- cbind(1,X)
  s <- length(a)
  delta <- matrix(delta, len, s, byrow = F)
  lt <- length(theta)+1
  tmp <- NULL
  for (i in 1:s) {
    tmp[[i]] <- sapply(1:n, function(j, X1, Z, l, a, theta, delta, h) {
      (X1[j, (l+1):(l+len)]%*%delta[,i]) * pnorm( drop( c(-1, Z[j,]) %*% c(a[i], theta) ) / h)
    }, X1=X1, Z=Z, l=l, a = a, theta = theta, delta = delta, h = h)
  }

  sm <- Reduce("+", tmp)

  da <- cbind(sm, Y, X)
  Q <- apply(da, 1, function(x) {
    x[2] - c(1,x[-c(1:2)])%*%beta - x[1]
  })

  return(sum(Q^2)/n)

}


## given beta, theta, estimate a

obj.f_a <- function(a, X, Y, Z, l, len=1, theta, beta, delta, h) {

  X <- as.matrix(X)
  n <- dim(X)[1]
  X1 <- cbind(1,X)
  s <- length(a)
  delta <- matrix(delta, len, s, byrow = F)
  lt <- length(theta)+1
  tmp <- NULL
  for (i in 1:s) {
    tmp[[i]] <- sapply(1:n, function(j, X1, Z, l, a, theta, delta, h) {
      (X1[j, (l+1):(l+len)]%*%delta[,i]) * pnorm( drop( c(-1, Z[j,]) %*% c(a[i], theta) ) / h)
    }, X1=X1, Z=Z, l=l, a = a, theta = theta, delta = delta, h = h)
  }

  sm <- Reduce("+", tmp)

  da <- cbind(sm, Y, X)
  Q <- apply(da, 1, function(x) {
    x[2] - c(1,x[-c(1:2)])%*%beta - x[1]
  })

  return(sum(Q^2)/n)

}

## given beta, a, estimate theta

obj.f_theta <- function(theta, X, Y, Z, l, len=1, a, beta, delta, h, lamb) {

  X <- as.matrix(X)
  n <- dim(X)[1]
  X1 <- cbind(1,X)
  s <- length(a)
  delta <- matrix(delta, len, s, byrow = F)
  lt <- length(theta)+1
  tmp <- NULL
  for (i in 1:s) {
    tmp[[i]] <- sapply(1:n, function(j, X1, Z, l, a, theta, delta, h) {
      (X1[j, (l+1):(l+len)]%*%delta[,i]) * pnorm( drop( c(-1, Z[j,]) %*% c(a[i], theta) ) / h)
    }, X1=X1, Z=Z, l=l, a = a, theta = theta, delta = delta, h = h)
  }

  sm <- Reduce("+", tmp)

  da <- cbind(sm, Y, X)
  Q <- apply(da, 1, function(x) {
    x[2] - c(1,x[-c(1:2)])%*%beta - x[1]
  })

  gamma <- 2.4
  ttheta <- abs(theta)
  penal <- ttheta

  id1 <- which(ttheta <= lamb)
  id2 <- which(ttheta > lamb & ttheta <= gamma*lamb)
  id3 <- which(ttheta > gamma*lamb)

  penal[id1] <- lamb * ttheta[id1]
  penal[id2] <- - (ttheta[id2]^2 - 2*gamma*lamb*ttheta[id2] + lamb^2) / (2*(gamma-1))
  penal[id3] <- (gamma+1) * lamb^2 / 2

  return(sum(Q^2)/n+sum(penal))

}



## given theta, a, estimate beta
lm.f_beta <- function(X, Y, Z, l, len=1, theta, a, h) {

  X <- as.matrix(X)
  n <- dim(X)[1]
  X1 <- cbind(1,X)
  d <- ncol(X1)
  s <- length(a)
  lt <- length(theta)+1

  XX <- X1
  for(i in 1:s) {
    XX <- cbind( XX, t(sapply(1:n, function(j, X1, Z, theta, a, h){
      X1[j,] * pnorm(drop( c(-1, Z[j,]) %*% c(a, theta) ) / h)
      }, X1=X1, Z=Z, theta= theta, a = a[i], h = h ))[,(l+1):(l+len),drop=F] )
  }

  object <- ncvreg(XX[,-1], c(Y), penalty = 'MCP')
  bic <-  n*log(object$loss/n) + colSums(object$beta!=0) * log(n) ## BIC
  lamb=object$lambda[which.min(bic)]
  psi.hat <- coef(object, lam = lamb)

  beta.hat <- psi.hat[1:d]
  delta.hat <- psi.hat[-c(1:d)]

  return(list(beta.hat = beta.hat, delta.hat = delta.hat, lamb = lamb))

}



smooth <- function(X, Y, Z, l, len=1, s.est, ini.a = c(0, rep(0, length.out = s.est-1)),
                   ini.theta = c(1, 0, 0, 0, 0),
                   tol= 1e-4, K = 10) {

  X <- as.matrix(X)
  Y <- c(Y)
  Z <- as.matrix(Z)
  n <- dim(X)[1]
  p <- dim(X)[2]

  k <- 1; d <- p+1
  a.hat0 <- ini.a
  theta.hat0 <- ini.theta
  h <- n^(-1)*log(n)
  lt <- length(theta.hat0) + 1
  beta.trace <- matrix(0, d, K)
  delta.trace <- matrix(0, s.est*len, K)
  theta.trace <- matrix(0, lt-1, K)
  a.trace <- matrix(0, s.est, K)

  while (k <= K) {

    ## given theta, a, estimate beta
    ord <- order(Z %*% theta.hat0)
    Zo <- Z[ord, ]
    Xo <- X[ord, ]
    Yo <- Y[ord]
    fit <- lm.f_beta(X=Xo, Y=Yo, Z=Zo, l=l, len=len, theta = theta.hat0, a = a.hat0, h = h)
    beta.hat <- fit$beta.hat
    delta.hat <- fit$delta.hat
    lamb <- fit$lamb

    ## given beta, theta, estimate a
    fit.a <- BBoptim(par = a.hat0, fn = obj.f_a, lower = rep(-2,s.est), upper = rep(2,s.est),
                     X=Xo, Y=Yo, Z=Zo, l=l, len=len, theta = theta.hat0, beta = beta.hat, delta = delta.hat, h = h )
    a.hat <- fit.a$par
    #cat('a=', a.hat, '\n')

    ## given beta, a, estimate theta
    fit.theta <- nloptr(x0 = theta.hat0, eval_f = obj.f_theta, lb = rep(-1,lt-1), ub = rep(1,lt-1),
                        opts = list('algorithm'='NLOPT_LN_BOBYQA', 'maxeval'=1000, 'xtol_rel'=1.0e-6, 'print_level'=0),
                        X=Xo, Y=Yo, Z=Zo, l=l, len=len, a = a.hat, beta = beta.hat, delta = delta.hat, h = h, lamb = lamb )
    theta.hat <- fit.theta$solution
    theta.hat <- theta.hat/norm(theta.hat, '2')
    theta.hat[is.nan(theta.hat)] <- 0

    #fit.theta <- BBoptim(par = theta.hat0, fn = obj.f_theta, lower = rep(-1,lt-1), upper = rep(1,lt-1),
    #                     l=l, a = a.hat, beta = beta.hat, delta = delta.hat, ZZ = ZZ, h = h, lamb = lamb )
    #theta.hat <- fit.theta$par/norm(fit.theta$par, "2")

    beta.trace[,k] <- beta.hat
    delta.trace[,k] <- delta.hat
    theta.trace[,k] <- theta.hat
    a.trace[,k] <- a.hat
    k <- k + 1
    if ( all( abs(a.hat0 - a.hat) < tol, abs(theta.hat0 - theta.hat) < tol)  ) {
      print("ok!")
      break
    }
    theta.hat0 <- theta.hat
    a.hat0 <- a.hat

  }

  ## given theta, a, estimate beta
  ord=order(Z %*% theta.hat)
  Zo <- Z[ord, ]
  Xo <- X[ord, ]
  Yo <- Y[ord]
  fit <- lm.f_beta(X=Xo, Y=Yo, Z=Zo, l=l, len=len, theta = theta.hat, a = a.hat, h = h)
  beta.hat <- fit$beta.hat
  delta.hat <- fit$delta.hat

  mse <- obj.f(X=Xo, Y=Yo, Z=Zo, l=l, len=len, theta = theta.hat, a = a.hat,
                      beta = beta.hat, delta = delta.hat, h = h)


  return(obj <- list(beta.hat = beta.hat, delta.hat = delta.hat, theta.hat = theta.hat,
                     a.hat = a.hat, mse = mse))

}


##
subm.fun <- function(X, Y, Z, l, len=1, ini.theta, tol= 1e-4, K= 10){

  X <- as.matrix(X)
  Y <- c(Y)
  Z <- as.matrix(Z)
  n <- dim(X)[1]
  p <- dim(X)[2]

  theta.hat0 <- ini.theta
  lt <- length(theta.hat0) + 1
  theta.trace <- matrix(0, lt-1, K)
  s.trace <- matrix(0, 1, K)
  k <- 1

  while (k <= K){

    ord <- order(Z %*% theta.hat0)
    Zo <- Z[ord, ]
    Xo <- X[ord, ]
    Yo <- Y[ord]

    ## choose the best split c
    c <- seq(0.5,1.5,0.1)
    split.fit <- NULL
    for(j in 1:length(c)) {
      split.fit[[j]] <- split.fun(Xo, Yo, l=l, len=len, c=c[j])
    }
    BIC <- sapply(split.fit, function(x){x$BIC}, simplify = T)
    mcp.cp1 <- sapply(split.fit, function(x){x$mcp.cp1}, simplify = F)
    mcp.cp1 <- mcp.cp1[[which.min(BIC)]]
    s.est <- length(mcp.cp1)
    numb <- sapply(split.fit, function(x){x$numb}, simplify = F)
    numb <- numb[[which.min(BIC)]]
    a.hat0 <- (Zo %*% theta.hat0)[numb]

    if(s.est >= 1) {

      ############### step 2
      ## smoothing estimate with initial (theta, a, beta, delta)
      #cat('ahat0=', a.hat0, '\n')
      sm.obj <- smooth(X=Xo, Y=Yo, Z=Zo, l=l, len=len, s.est=s.est, ini.a=a.hat0,
                       ini.theta=theta.hat0, tol=tol/10, K=K)
      theta.hat <- sm.obj$theta.hat

      theta.trace[,k] <- theta.hat
      s.trace[k] <- s.est
      k <- k + 1
      if ( all(theta.hat0 - theta.hat < tol) ) {
        print("Finish!")
        break
      }
      theta.hat0 <- theta.hat

    } else {
      sm.obj <- lm(Y~X)
      break
    }

  }

  return(list(s.est = s.est, sm.obj = sm.obj, s.trace = s.trace, theta.trace = theta.trace))

}




pointest <- function(X, Y, Z, l, len=1, s, obj){

  X <- as.matrix(X)
  Y <- c(Y)
  Z <- as.matrix(Z)
  n <- dim(X)[1]
  p <- dim(X)[2]

  X1 <- cbind(1, X)
  d <- ncol(X1)

  if (s>=1){

    beta <- obj$beta.hat
    delta <- obj$delta.hat
    delta <- matrix(delta, len, s, byrow = F)
    theta <- obj$theta.hat
    a <- obj$a.hat

    tmp <- NULL

    for (i in 1:s) {

      ind <- which(Z %*% theta > a[i])
      ktemp <- matrix(0, nrow=n, ncol=1)
      ktemp[ind,] <- 1
      tmp[[i]] <- X1[,(l+1):(l+len),drop=F] %*% delta[,i,drop=F] * ktemp

    }

    sm <- Reduce("+", tmp)

    Yhat <- X1 %*% beta + sm

  } else Yhat <- X1 %*% obj$coefficients

  mse <- mean((Y-Yhat)^2)

  return(list(Yhat=Yhat, mse=mse))

}




