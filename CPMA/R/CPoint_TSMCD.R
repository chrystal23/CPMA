


cumcp<-function(Y, X) {

  p=dim(X)[2]
  n=dim(X)[1]
  cp=rep(0,n)

  if(n>2*p){

    for(i in (p+1):(n-p)) {
      lmr1 = lm(Y[1:i]~as.matrix(X[1:i,]))
      lmr2 = lm(Y[(i+1):n]~as.matrix(X[(i+1):n,]))
      cp[i] = i*sum((lmr1$residuals)^2) + (n-i)*sum((lmr2$residuals)^2)
    }

    cpm=min(cp[which(cp>0)])
    return(which(cp==cpm))

  }
  else
    return(0)

}




# the main funtion (The input Y,X has to be in order based on the thresholding variable)
TSMCD <- function(Y, X, c) {

  n=length(Y)
  p=dim(X)[2]
  m=ceiling(c*sqrt(n))
  q=floor(n/m)

  x=NULL
  y=NULL

  X_temp=cbind(1,X)
  Y_temp=c(Y)

  x[[1]] = sqrt(n-(q-1)*m) * as.matrix(X_temp[1:(n-(q-1)*m),])
  y[[1]] = sqrt(n-(q-1)*m) * Y_temp[1:(n-(q-1)*m)]

  for(i in 2:q) {
    x[[i]] = sqrt(m) * as.matrix(X_temp[(n-(q-i+1)*m+1):(n-(q-i)*m),])
    y[[i]] = sqrt(m) * Y_temp[(n-(q-i+1)*m+1):(n-(q-i)*m)]
  }

  K_temp=matrix(0, nrow = q, ncol=q, byrow=TRUE)
  for(i in 1:q)
    K_temp[i,1:i]=rep(1,i) # a lower triangular matrix (for the matrix representation of least square)

  X_temp1 <- lapply(1:length(x), function(j, mat, list)
    kronecker(K_temp[j, , drop = FALSE], x[[j]]), mat=K_temp, list=x)
  Xn=do.call("rbind",X_temp1)

  Yn=NULL
  for(i in 1:q) {
    Yn=c(Yn,y[[i]])
  }

  group=kronecker(c(1:q),rep(1,p+1))
  fit <- grpreg(Xn[,-1], Yn, group[-1], penalty="grSCAD")
  mcp.coef <-grpreg::select(fit,"BIC")$beta

  mcp.coef.s=sum(abs(mcp.coef))
  mcp.coef.v.m=abs(matrix(c(mcp.coef),q,(p+1),byrow=T))
  mcp.coef.m=c(apply(mcp.coef.v.m,1,max))
  mcp.cp=which(mcp.coef.m!=0)

  if(q>10) {

    if(length(mcp.cp)>1)
    {
      for(i in 2:length(mcp.cp))
      {
        if(mcp.cp[i]-mcp.cp[i-1]==1)
          mcp.cp[i]=0
      }
    }  # pick out the possible range of change points

    mcp.cp1=mcp.cp[mcp.cp>1&mcp.cp<q]

    d1=length(mcp.cp1)

    if(d1==0) {
      mcpcss.cp=NULL
      adcpcss.cp=mcpcss.cp
    }

    if(d1>=1) {
      mcpcss.cp=NULL

      for(i in 1:d1) {
        t=mcp.cp1[i]
        id=c((n-(q-t+2)*m+1):(n-(q-t-1)*m))
        cp.t=cumcp(Y[id],X[id,])[[1]]
        if(cp.t>0)
          mcpcss.cp=c(mcpcss.cp,cp.t+n-(q-t+2)*m+1-1)
      } # get the list of change points (mcpcss.cp) by minimizing Q

      id2=c(mcpcss.cp,n)
      id3=c(0,mcpcss.cp)

      tt2=which(id2-id3<p+2)
      tt2[which(tt2==length(id2))]=length(id2)-1

      if(length(tt2)>0) {
        mcpcss.cp=mcpcss.cp[-tt2]
      } # remove some change points so that the number of uncensored obervations in each segment is not smaller than p+2

      tt=which(abs(diff(mcpcss.cp))<2*m)
      if(length(tt)>0)
        mcpcss.cp=mcpcss.cp[-tt]
    }

  }


  if(q<=10) {

    mcp.cp1=mcp.cp[mcp.cp>1&mcp.cp<q]

    d1=length(mcp.cp1)

    if(d1==0) {
      mcpcss.cp=NULL
      adcpcss.cp=mcpcss.cp
    }

    if(d1>=1) {
      mcpcss.cp=NULL

      for(i in 1:d1) {
        t=mcp.cp1[i]
        if(t==2) {
          id=c(1:(n-(q-t-1)*m))
          cp.t=cumcp(Y[id],X[id,])[[1]]
          if(cp.t>0)
            mcpcss.cp=c(mcpcss.cp,cp.t)
        }
        if(t==q) {
          id=c((n-(q-t+2)*m+1):n)
          cp.t=cumcp(Y[id],X[id,])[[1]]
          if(cp.t>0)
            mcpcss.cp=c(mcpcss.cp,cp.t+n-(q-t+2)*m+1-1)
        }
        if(t<q&t>2) {
          id=c((n-(q-t+2)*m+1):(n-(q-t-1)*m))
          cp.t=cumcp(Y[id],X[id,])[[1]]
          if(cp.t>0)
            mcpcss.cp=c(mcpcss.cp,cp.t+n-(q-t+2)*m+1-1)
        }
      }

      id2=c(mcpcss.cp,n)
      id3=c(0,mcpcss.cp)

      tt2=which(id2-id3<p+2)
      tt2[which(tt2==length(id2))]=length(id2)-1

      if(length(tt2)>0) {
        mcpcss.cp=mcpcss.cp[-tt2]
      }

      tt=which(abs(diff(mcpcss.cp))<2*m)
      if(length(tt)>0)
        mcpcss.cp=mcpcss.cp[-tt]

    }

  }


  # with the estimated change points (mcpcss.cp), use weighted least squares to compute the final estimates of the coefficients

  if(length(mcpcss.cp)>=1 && mcpcss.cp[1]!=0) {

    id2=c(mcpcss.cp,n)
    x=NULL
    y=NULL
    x0=NULL
    y0=NULL

    X_temp=cbind(1,X)
    Y_temp=c(Y)

    x[[1]] = sqrt(mcpcss.cp[1]) * as.matrix(X_temp[1:mcpcss.cp[1],])
    y[[1]] = sqrt(mcpcss.cp[1]) * Y_temp[1:mcpcss.cp[1]]
    x0[[1]] = as.matrix(X_temp[1:mcpcss.cp[1],])
    y0[[1]] = Y_temp[1:mcpcss.cp[1]]

    for(i in 2:length(id2)) {
      x[[i]] = sqrt(id2[i]-id2[i-1]) * as.matrix(X_temp[(id2[i-1]+1):id2[i],])
      y[[i]] = sqrt(id2[i]-id2[i-1]) * Y_temp[(id2[i-1]+1):id2[i]]
      x0[[i]] = as.matrix(X_temp[(id2[i-1]+1):id2[i],])
      y0[[i]] = Y_temp[(id2[i-1]+1):id2[i]]
    }

    K_temp=matrix(0, nrow = length(id2), ncol=length(id2), byrow=TRUE)
    for(i in 1:length(id2))
      K_temp[i,1:i]=rep(1,i)

    X_temp1 <-lapply(1:length(x), function(j, mat, list)
      kronecker(K_temp[j,,drop=FALSE], x[[j]]), mat=K_temp, list=x)
    Xn=do.call("rbind",X_temp1)

    X0_temp1 <-lapply(1:length(x0), function(j, mat, list)
      kronecker(K_temp[j,,drop=FALSE], x0[[j]]), mat=K_temp, list=x0)
    X0=do.call("rbind",X0_temp1)

    Yn=NULL
    for(i in 1:length(id2)) {
      Yn=c(Yn,y[[i]])
    }

    Y0=NULL
    for(i in 1:length(id2)) {
      Y0=c(Y0,y0[[i]])
    }


    object <- plus(Xn,Yn,method="scad",gamma=2.4,intercept = F,normalize =F,eps = 1e-30)
    #step 1 estimate coef using BIC.
    bic=log(dim(Xn)[1])*object$dim+dim(Xn)[1]*log(as.vector((1-object$r.square)*sum(Yn^2))/length(Yn))
    step.bic=which.min(bic)

    mcp.coef <- coef(object,lam=object$lam[step.bic])

    ep=(Yn-Xn%*%mcp.coef)
    sigmaep=sum(ep^2)/length(Yn)

    ep0=(Y0-X0%*%mcp.coef)
    sigmaep0=sum(ep0^2)/length(Y0)

    Yhat=X0%*%mcp.coef

    mcp.coef.v.m=matrix(c(mcp.coef),length(id2),p+1,byrow=T)
    mcp.coef.m=c(apply(abs(mcp.coef.v.m),1,max))
    mcp.m.id=which(mcp.coef.m==0)
    if(length(mcp.m.id)>0) {
      mcp.coef <-c(t(mcp.coef.v.m[-mcp.m.id,]))
      mcpcss.cp=mcpcss.cp[-c(mcp.m.id-1)] # remove the segments (change points) with coefficient 0
    }
  }


  if(mcpcss.cp[1]==0 || length(mcpcss.cp)==0) {

    id2=c(mcpcss.cp,n)
    X_temp=cbind(1,X)
    Y_temp=c(Y)

    Xn = sqrt(n) * as.matrix(X_temp)
    Yn = sqrt(n) * Y_temp
    X0 = as.matrix(X_temp)
    Y0 = Y_temp

    object <- plus(Xn,Yn,method="scad",gamma = 2.4,intercept = F,normalize =F,eps = 1e-30)
    #step 1 estimate  coef using BIC.
    bic=log(dim(Xn)[1])*object$dim+dim(Xn)[1]*log(as.vector((1-object$r.square)*sum(Yn^2))/length(Yn))
    step.bic=which.min(bic)

    mcp.coef <- coef(object, lam=object$lam[step.bic])

    ep=(Yn-Xn%*%mcp.coef)
    sigmaep=sum(ep^2)/length(Y)

    ep0=(Y0-X0%*%mcp.coef)
    sigmaep0=sum(ep0^2)/length(Y0)

    Yhat=X0%*%mcp.coef

    mcp.coef.v.m=matrix(c(mcp.coef),length(id2),p+1,byrow=T)
    mcp.coef.m=c(apply(abs(mcp.coef.v.m),1,max))
    mcp.m.id=which(mcp.coef.m==0)
    if(length(mcp.m.id)>0) {
      mcp.coef <-c(t(mcp.coef.v.m[-mcp.m.id,]))
      mcpcss.cp=mcpcss.cp[-c(mcp.m.id-1)]
    }

  }

  return(list(cp=mcpcss.cp, coef=mcp.coef, sigmaep0=sigmaep0, Yhat=Yhat))
}







TSMCD.pred<- function(Y, X, Z, thre, mcp.coef){

  Y <- c(Y)
  X1 <- cbind(1, as.matrix(X))
  Z <- c(Z)

  if (length(thre)>=1){

    tmp <- NULL

    for (i in 1:length(thre)) {

      ind <- which(Z>thre[i])
      ktemp = matrix(0, nrow=length(Y), ncol=1)
      ktemp[ind,] = 1
      tmp[[i]] <- X1 %*% mcp.coef[i+1,] * ktemp

    }

    sm <- Reduce("+", tmp)

    Yhat <- X1 %*% mcp.coef[1,] + sm

  } else {

    Yhat <- X1 %*% t(mcp.coef)

  }

  ep0 <- Y-Yhat
  sigmaep0 <- sum(ep0^2)/length(Y)

  return(list(Yhat=Yhat, sigmaep0=sigmaep0))

}




