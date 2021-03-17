
#' Fit a Monotone Smoothing Spline
#' @export
#' @param x numeric variable
#' @param y numeric variable
#' @import "NlcOptim"
#' @import "kit"
#' @import "graphics"
#' @import "stats"





snorm<-function(x,y){

  # library(NlcOptim)
  # library(pracma)
  # library(nleqslv)

  N<-20

  opt.rss<-rep(0,N)

  wgt.mat<-list()

  d.mat<-list()

  k.mat<-list()



  for(i in 1:N){

    set.seed(N)

    k<-sort(runif(i,0,10))

    w<-runif(i,0,1)

    sig<-0.3

    x.mat<-matrix(rep(0,length(k)*length(x)),nrow = length(x),ncol = length(k))

    for(j in 1:length(k)){

      phi<-1-pnorm(x,k[j],sig)

      x.mat[,j]<-phi

    }

    rss<-function(w) {

      yhat<-x.mat %*% w

      rss<-sum((y-yhat)^2)

      return(rss)

    }

    confun<-function(w) {

      f<-sum(w)-1

      return(list(ceq=f,c=NULL))

    }

    ub<-rep(1,ncol(x.mat))
    lb<-rep(0.000001,ncol(x.mat))

    wgt<-solnl(X=w,objfun = rss,confun = confun,lb=lb,ub=ub)$par

    resid<-solnl(X=w,objfun = rss,confun = confun,lb=lb,ub=ub)$fn

    ######
    opt.rss[i]<-resid

    wgt.mat[[i]]<-wgt

    d.mat[[i]]<-x.mat

    k.mat[[i]]<-k
  }



  ######

  best.opt<-min(opt.rss)

  index<-which.min(opt.rss)

  best.xmat<-d.mat[[index]]

  best.wgt<-wgt.mat[[index]]

  best.k<-k.mat[[index]]

  ### build data frame


  wgtk<-data.frame(best.k,best.wgt)

  new.wgtk<-subset(wgtk,wgtk$best.wgt>0.0001)

  best.wgtk<-new.wgtk[order(new.wgtk$best.wgt,decreasing = T),]

  #### pick maximum weights

  ##library(kit)


  id<-topn(best.wgtk$best.wgt,n=length(best.wgtk$best.wgt))

  new.wgtk<-data.frame(best.wgtk,id)

  for(i in 1:length(new.wgtk$best.wgt)){

    k<-new.wgtk$best.k[id<=i]

    w<-new.wgtk$best.wgt[id<=i]

    sig<-0.3

    x.mat<-matrix(rep(0,length(k)*length(x)),nrow = length(x),ncol = length(k))

    for(j in 1:length(k)){

      phi<-1-pnorm(x,k[j],sig)

      x.mat[,j]<-phi

    }

    yhat<-x.mat %*% w

    rss<-sum((y-yhat)^2)

    if(rss<=best.opt)

      best.opt<-rss
      best.k<-k
      best.xmat<-x.mat
      best.wgt<-w


  }

  ##### prediciton

  ypre<-best.xmat %*% best.wgt

  #plot(x,y)

  lines(x,ypre,col='red')

  data.frame(optimal.k=best.k,estimated.weights=best.wgt)


}




