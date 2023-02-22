#' @title confidence intervals for fitted values of (Generalized) Linear Model
#' @description determines confidence intervals for fitted values of (Generalized) Linear Model to be added to plot
#' @param model.res model fitted using lm, glm, or glm.nb.
#' @param resol numeric scalar; number values of covariates (between their min and max) for which
#' confidence intervals are determined (default: 100).
#' @param confidence numeric; level of confidence (default: 0.95)
#' @param use character (scalar or vector), names of predictors in model for which confidence
#' intervals are to be determined. Note that only covariates need to be named here (factors are
#' included by default); default: "all".
#' @param set.offs.to value which the offset should have (if any was included) when determining
#' the fitted value and its CI. When left at its default the offset will be set to its average.
#' @return A data frame with columns for the predictors in the model, the fitted values and their
#' confidence intervals (column headers: "lwr", "upr")
#' @details The function is largely a wrapper for predict.lm and predict.glm, respectively. But it
#' it adds additional stuff aiming at making plotting CI for a fitted model easier.
#' @examples
#' ##create some data:
#' n=100
#' xdata=data.frame(
#'   xcov=runif(n=n, min=-3, max=3),
#'   ycov=rnorm(n=n),
#'   xfac=as.factor(sample(x=letters[1:2], size=n, replace=T))
#' )
#' xdata$resp=-1*xdata$xcov^2*ifelse(xdata$xfac=="b", 1, 0)
#' xdata$resp[xdata$xfac=="a"]=0.5*xdata$xcov[xdata$xfac=="a"]
#' xdata$resp=xdata$resp+rnorm(n)
#' full=lm(resp~xfac*(xcov+I(xcov^2))+ycov, data=xdata)
#' ci=CIglm(model.res=full, use="xcov")
#' par(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.3, 0), las=1, tcl=-0.15)
#' plot(x=xdata$xcov, y=xdata$resp, pch=as.numeric(xdata$xfac), xlab="predictor", ylab="response")
#' lines(x=ci$xcov[ci$xfac=="a"], y=ci$fit[ci$xfac=="a"], lty=2, lwd=2)
#' lines(x=ci$xcov[ci$xfac=="b"], y=ci$fit[ci$xfac=="b"], lty=3, lwd=2)
#' lines(x=ci$xcov[ci$xfac=="a"], y=ci$lwr[ci$xfac=="a"], lty=2)
#' lines(x=ci$xcov[ci$xfac=="a"], y=ci$upr[ci$xfac=="a"], lty=2)
#' lines(x=ci$xcov[ci$xfac=="b"], y=ci$lwr[ci$xfac=="b"], lty=3)
#' lines(x=ci$xcov[ci$xfac=="b"], y=ci$upr[ci$xfac=="b"], lty=3)
#' @export

ci.glm<-function(model.res, resol=100, level=0.95, use="all", set.offs.to=NULL){
  if(class(model.res)[1]=="lm"){
    xfamily="gaussian"
  }else if(class(model.res)[1]!="negbin"){
    xfamily=as.character(model.res$call)[3]
  }else if(class(model.res)[1]=="negbin"){
    xfamily="negbin"
  }else{
    stop("unsupported model")
  }
  model.terms=names(model.res$model)[-1]
  if(length(model.res$offset)>0){model.terms=model.terms[substr(x=model.terms, start=1, stop=6)!="offset"]}
  if(use[1]=="all"){use=model.terms}
  #create new data to be used to determine fitted values:
  new.data=vector("list", length(model.terms))
  use=model.terms%in%use
  for(i in 1:length(model.terms)){
    if(is.factor(model.res$model[, model.terms[i]]) &	 use[i]){
      new.data[[i]]=levels(model.res$model[, model.terms[i]])
    }else if(!is.factor(model.res$model[, model.terms[i]]) & use[i]){
      new.data[[i]]=seq(from=min(model.res$model[, model.terms[i]]), to=max(model.res$model[, model.terms[i]]), length.out=resol)
    }else if(is.factor(model.res$model[, model.terms[i]]) &	!use[i]){
      new.data[[i]]=levels(model.res$model[, model.terms[i]])
    }else{
      new.data[[i]]=0
    }
  }
  new.data=data.frame(expand.grid(new.data))
  names(new.data)=model.terms
  #deal with potential squared terms
  new.data=new.data[, !names(new.data)%in%paste("I(", model.terms[use], "^2)", sep=""), drop=F]

  if(length(model.res$offset)>0){
    new.data=data.frame(new.data, mean(model.res$offset))
    xx=names(model.res$model)[substr(x=names(model.res$model), start=1, stop=6)=="offset"]
    xx=gsub(xx, pattern="offset(log(", replacement="", fixed=T)
    xx=gsub(xx, pattern="offset(", replacement="", fixed=T)
    xx=gsub(xx, pattern="))", replacement="", fixed=T)
    xx=gsub(xx, pattern=")", replacement="", fixed=T)
    names(new.data)[ncol(new.data)]=xx
    if(!is.null(set.offs.to)){
      new.data[, ncol(new.data)]=set.offs.to
    }
  }
  if(xfamily!="gaussian"){
    ci.plot=predict.glm(object=model.res, newdata=new.data, type="link", se.fit=T)
    ci.plot=data.frame(
      orig=ci.plot$fit,
      lwr=ci.plot$fit-ci.plot$se.fit*abs(qt(p=(1-level)/2, df=model.res$df.residual)),
      upr=ci.plot$fit+ci.plot$se.fit*abs(qt(p=(1-level)/2, df=model.res$df.residual))
    )
  }else{
    ci.plot=predict.lm(object=model.res, newdata=new.data, interval="confidence")
  }
  if(xfamily=="binomial"){ci.plot=exp(ci.plot)/(1+exp(ci.plot))}
  if(xfamily=="poisson" | xfamily=="negbin"){ci.plot=exp(ci.plot)}
  return(data.frame(new.data, ci.plot))
}
