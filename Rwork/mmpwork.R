
#############################################
# marginal model plots
#############################################
marginal.model.plot <- function(...){mmp(...)}
mmp <- function(object, ...){UseMethod("mmp")}

mmp.lm <-
function (object, u, mean = TRUE, sd = FALSE,
    xlab = deparse(substitute(u)), degree = 1, span = 2/3, key="topleft", 
    lineColors = c("blue","red"), ...)
{
    if (missing(u)) {
        xlab <- "Fitted values"
        u <- fitted(update(object,na.action=na.exclude))
    }
    na.cases <- attr(object$model,"na.action")
    if(length(na.cases)>0) u <- u[-na.cases]
    zpred <- function(...){pmax(predict(...),0)}
    plot(u, object$model[, 1], xlab = xlab, ylab = colnames(object$model[1]),
        ...)
    if(!is.null(key))legend(key,c("Data","Model"),lty=c(1,2),col=lineColors)
#    mtext("Data",col=lineColors[1],adj=0,cex=.7)
#    mtext("Model",col=lineColors[2],adj=1,cex=.7)
    loess.y <- loess(object$model[, 1] ~ u, degree = degree,
        span = span)
    loess.yhat <- loess(predict(object) ~ u, degree = degree,
        span = span)
    new <- seq(min(u), max(u), length = 200)
    if (mean == TRUE) {
        lines(new, predict(loess.y, data.frame(u = new)), lty = 1,
            col = lineColors[1])
        lines(new, predict(loess.yhat, data.frame(u = new)),
            lty = 2, col = lineColors[2])
    }
    if (sd == TRUE) {
        loess.y.var <- loess(residuals(loess.y)^2 ~ u, degree = degree,
            span = span)
        lines(new, predict(loess.y, data.frame(u = new)) +
            sqrt(zpred(loess.y.var,data.frame(u = new))), lty = 1, col = lineColors[1])
        lines(new, predict(loess.y, data.frame(u = new)) - sqrt(zpred(loess.y.var,
            data.frame(u = new))), lty = 1, col = lineColors[1])
        loess.yhat.var <- loess(residuals(loess.yhat)^2 ~ u,
            degree = degree, span = span)
        s2 <- summary(object)$sigma^2
        lines(new, predict(loess.yhat, data.frame(u = new)) +
            sqrt(s2 + zpred(loess.yhat.var, data.frame(u = new))),
            lty = 2, col = lineColors[2])
        lines(new, predict(loess.yhat, data.frame(u = new)) -
            sqrt(s2 + zpred(loess.yhat.var, data.frame(u = new))),
            lty = 2, col = lineColors[2])
    }
}
  
mmp.glm <- function (object, u, mean = TRUE, sd = FALSE, 
    xlab = deparse(substitute(u)), degree = 1, span = 2/3, key="topleft",
    lineColors = c("blue", "red"), ...) 
{
    if (missing(u)) {
        xlab <- "Linear Predictor"
        u <- fitted(update(object,na.action=na.exclude))
    }    
    na.cases <- attr(object$model,"na.action")
    if(length(na.cases)>0) u <- u[-na.cases]
    fr.mmp <- function(family, x) {
        if (family == "binomial") 
            pmax(0, pmin(1, x))
        else if (family == "poisson") 
            pmax(0, x)
        else if (family == "gamma") 
            pmax(0, x)
        else x
    }
    response <- object$model[, 1]
    fam <- object$family$family
    if (is.matrix(response)) 
        response <- response[, 1]/apply(response, 1, sum)
    plot(u, response, xlab = xlab, ylab = colnames(object$model[1]), 
        ...)
    if(!is.null(key))legend(key,c("Data","Model"),lty=c(1,2),col=lineColors,...)
    loess.y <- loess(response ~ u, degree = degree, span = span)
    loess.yhat <- loess(predict(object, type = "response") ~ 
        u, degree = degree, span = span)
    new <- seq(min(u), max(u), length = 200)
    pred.loess.y <- fr.mmp(fam, predict(loess.y, data.frame(u = new)))
    pred.loess.yhat <- fr.mmp(fam, predict(loess.yhat, data.frame(u = new)))
    if (mean == TRUE) {
        lines(new, pred.loess.y, lty = 1, col = lineColors[1])
        lines(new, pred.loess.yhat, lty = 2, col = lineColors[2])
    }
    if (sd == TRUE) {
        loess.y.var <- loess(residuals(loess.y)^2 ~ u, degree = degree, 
            span = span)
        pred.loess.y.var <- pmax(0, predict(loess.y.var, data.frame(u = new)))
        lines(new, fr.mmp(fam, pred.loess.y + sqrt(pred.loess.y.var)), 
            lty = 1, col = lineColors[1])
        lines(new, fr.mmp(fam, pred.loess.y - sqrt(pred.loess.y.var)), 
            lty = 1, col = lineColors[1])
        loess.yhat.var <- loess(residuals(loess.yhat)^2 ~ u, 
            degree = degree, span = span)
        pred.loess.yhat.var <- pmax(0, predict(loess.yhat.var, 
            data.frame(u = new)))
        varfun <- summary(object)$dispersion * object$family$variance(predict(object, 
            type = "response"))/if (!is.null(object$prior.weights)) 
            object$prior.weights
        else 1
        loess.varfun <- loess(varfun ~ u, degree = degree, span = span)
        pred.loess.varfun <- pmax(0, predict(loess.varfun, data.frame(u = new)))
        sd.smooth <- sqrt(pred.loess.yhat.var + pred.loess.varfun)
        lines(new, fr.mmp(fam, pred.loess.yhat + sd.smooth), 
            lty = 2, col = lineColors[2])
        lines(new, fr.mmp(fam, pred.loess.yhat - sd.smooth), 
            lty = 2, col = lineColors[2])
    }
}

 
mmps <- function(object,vars=~.,fitted=TRUE,layout=NULL,ask,...){
  vars <- update(object,vars,na.action=NULL,method="model.frame")
  dataClasses <- attr(attr(vars,"terms"),"dataClasses")[-1]
  terms <- names(dataClasses)[dataClasses == "numeric"]
  nt <- length(terms)+fitted
  if(is.null(layout)){
   layout <- switch(min(nt,9),
                           c(1,1),c(1,2),c(2,2),c(2,2),c(3,2),c(3,2),
                           c(3,3),c(3,3),c(3,3))}
  op<-par(no.readonly=TRUE)
  ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
  on.exit(par(op))
  if(prod(layout) > 1)
    par(mfrow=layout,ask=ask,mai=c(.6,.6,.1,.1),mgp=c(2,1,0),cex.lab=1.0,cex=0.7)
    else par(mfrow=layout,ask=ask)
  for (term in terms){
    j <- match(term,names(vars))
    mmp(object,vars[,j],xlab=term,...)}
  if(fitted==TRUE) mmp(object,...)
  }
  
data(sleep1,package="alr3")
s1 <- lm(TS ~log(BodyWt)+log(BrainWt)+log(GP)+log(Life)+P+factor(D),sleep1)
s2 <- glm(TS ~log(BodyWt)+log(BrainWt)+log(GP)+log(Life)+P+factor(D),data=sleep1)
sleep1$D1 <- ifelse(sleep1$D<=2,1,2)
mmp(s1)
mmp(s2)
mmps(s1,pch=sleep1$D1)
mmps(s1,pch=sleep1$D1,col=sleep1$D1)
mmps(s1,~sqrt(BodyWt))
s4 <- update(s1,subset=-12)
s5 <- update(s2,subset=-12)
mmps(s1,~.+sqrt(BodyWt),sd=TRUE)
mmps(s2,~.+sqrt(BodyWt),sd=TRUE)
mmps(s4,~.+sqrt(BodyWt),sd=TRUE)
mmps(s5,~.+sqrt(BodyWt),sd=TRUE)


plot(1:10, (-4:5)^2, main="Parabola Points", xlab="xlab")
mtext("10 of them",col="red")
for(s in 1:4)
    mtext(paste("mtext(..., line= -1, {side, col, font} = ",s,
          ", cex = ", (1+s)/2, ")"), line = -1,
          side=s, col=s, font=s, cex= (1+s)/2)
mtext("mtext(..., line= -2)", line = -2)
mtext("mtext(..., line= -2, adj = 0)", line = -2, adj =0)
