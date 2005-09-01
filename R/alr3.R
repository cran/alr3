#  R version of the alr3 package

#######################################################################
# Chapter 2
#######################################################################
# method to return sigma hat
   
sigma.hat <- function(object){UseMethod("sigma.hat")}
sigma.hat.default <- function(object){summary(object)$sigma}
sigma.hat.glm <- function(object){sqrt(summary(object)$dispersion)}
### sd function for splus only; idential to R function

#######################################################################
#  confidence intervals 
# This could probably be rewritted as intervals.lm to work with the
# intervals function in nlme

#intervals.default <- function(object,...){conf.intervals(object,...)}

conf.intervals <- function(object,level=.95,f=qnorm((1-level)/2))
   {UseMethod("conf.intervals")}
conf.intervals.default <- function(object,level=.95,f=qnorm((1-level)/2)) {
  ans <- cbind(coef(object),
   coef(object) + f * sqrt(diag(vcov(object))),
   coef(object) - f * sqrt(diag(vcov(object))))
  dimnames(ans)[[2]] <- c("Coef est", "Lower", "Upper")
  ans}
conf.intervals.lm <- function(object,level=.95,f=qt((1-level)/2,object$df.residual)){
 conf.intervals.default(object,alpha,f)}
 
#######################################################################
# Chapter 4
#######################################################################

####  Case resampling bootstrap, based on work by Lexin Li
####  Includes error handling

next.boot <- function(object,sample){UseMethod("next.boot")}
next.boot.default <- function(object,sample){ 
   assign("boot.sample",sample,inherits=TRUE)
# next line assures resampling only rows in the original subset 9/1/2005
   update(object,data=model.frame(object),subset=boot.sample)}
next.boot.nls <- function(object,sample){
# modify to assure resampling only rows in the original subset 9/1/2005
   update(object,subset=sample,start=coef(object),
    data=data.frame(update(object,model=TRUE)$model))}

boot.case <- function(object, f=coef, B=999){UseMethod("boot.case")}
boot.case.default <- function (object, f=coef, B = 999) 
{
    n <- length(resid(object))
    options(show.error.messages = FALSE)
    on.exit(options(show.error.messages = TRUE))
    coef.boot <- NULL
    count.error <- 0
    i <- 0
    while (i < B) {
        obj.nl <- try(next.boot(object, sample=sample(1:n, replace = TRUE)))
        if (is.null(class(obj.nl))) {
            count.error <- 0
            i <- i + 1
            coef.boot <- rbind(coef.boot, f(obj.nl))
        }
        else {
            if (class(obj.nl)[1] != "try-error") {
                count.error <- 0
                i <- i + 1
                coef.boot <- rbind(coef.boot, f(obj.nl))
            }
            else {
                count.error <- count.error + 1
            }
        }
        if (count.error >= 25) {
            options(show.error.messages = TRUE)
            stop("25 consecutive bootstraps did not converge.  Bailing out.")}
    } 
    return(coef.boot)
}


###########################################################################
# Chapter 5  Pure error anova
###########################################################################
pure.error.anova <- function(mod) {
 if (!(class(mod) == "lm")) stop("Pure error for lm objects only")
 if (is.null(mod$model)) mod <- update(mod,model=TRUE)
 #save.seed <- .Random.seed  # keep the current random number seed
 p <- dim(m1$model)[2] -1
 mod$model$Lack.of.Fit <- 
   factor(random.lin.comb(model.matrix(mod),101319853))
 #set.seed(save.seed) # restore random number seed
 if (length(levels(mod$model$Lack.of.Fit)) == length(mod$model$Lack.of.Fit)) 
  anova(mod) else {
  anova(update(mod, ~.+Lack.of.Fit,data=mod$model,singular.ok=TRUE))}
 }
 
random.lin.comb <- function(X,seed=NULL) {UseMethod("random.lin.comb")}

random.lin.comb.default <- function(X,seed=NULL) {
 if(!is.null(seed)) set.seed(seed)
 std <- function(x){ 
    s <- sd(x)
    if( s > 0) (x-mean(x))/s else x}
 apply(X,2,std)%*% as.vector(2*rnorm(dim(X)[2])-1)
 }
 
random.lin.comb.lm <- function(X,...) {
 if(is.null(X$model)) X <- update(X,model=TRUE)
 random.lin.comb(X$model[,-1],...)}


###########################################################################
# Chapter 6  Delta Method
###########################################################################

delta.method <- function(object,g,parameter.prefix="b",print=TRUE)
  {UseMethod("delta.method")}

### lm, glm, and others with unnamed parameters:
delta.method.default<-function(object, g, parameter.prefix="b",print=TRUE)
{
   if(!is.character(g)) stop("function argument must be a string")
   g<-parse(text=g)
# variance of estimated parameters
   V <- vcov(object)        # vcov is in R base, in alr3 for S-Plus
   para.val<-coef(object)   # values of coef estimates
   q<-length(para.val)      # number of coef estimates
   para.name<-paste(parameter.prefix,0:(q-1),sep="") #coef names
   compute.delta.method(V,g,para.val,para.name,print)}
   
# nls has named parameters so parameter.prefix is ignored
delta.method.nls<-function(object, g, parameter.prefix=NULL,print, ...)
{
   if(!is.character(g)) stop("function argument must be a string")
   g<-parse(text=g)
   V<-vcov(object)  # variance of estimated parameters
   compute.delta.method(V,g,coef(object),names(coef(object)),
       print=print)}
   
# computes g evaluated at the data, and t(g')Vt(g'), the estimated
# standard error
compute.delta.method <- function(Var,g,values,para.name,print=TRUE){
   q <- length(values)
   for(i in 1:q) {assign(para.name[i], values[i])}
   est<-eval(g) 
   names(est)<-NULL
# derivative of function g of parameters
   gd<-NULL
   for(i in 1:q) {gd<-c(gd, eval(D(g, para.name[i])))}
# compute se
   se.est<-as.vector(sqrt(t(gd) %*% Var %*% gd))
# print
   if(print){
     cat("Functions of parameters:  ")
     print.default(g)
     cat("Estimate =", est, "with se =", se.est, "\n")}
# output
   ans<-list(estimate=est, se=se.est, func=g)
   return(invisible(ans)) 
}


##################################################################
# pod models; Cook and Weisberg (2004), American Statistician
##################################################################

# This code was written by Lexin Li.  It was modified and simplified by S. Weisberg
# to work with standard R stuff in nls.

pod <- function(x, ...) {UseMethod("pod")}

pod.lm <- function(x,group,mean.function,control,...) {
 call <- x$call
 call[[1]] <- as.name("pod")
 if(!missing(group)) call$group <- as.name(as.character(substitute(group)))
 if(!missing(mean.function)) call$mean.function <- 
     as.name(as.character(substitute(mean.function)))
 if(!missing(control)) call$control <- 
     as.name(as.character(substitute(control)))
 eval(call)}

pod.formula <-
function (formula, data=sys.parent(), group, subset, weights=NULL, na.action, 
    mean.function=c("pod","common","parallel","general"),
    singular.ok = FALSE, contrasts = NULL, offset, control=nls.control(), ...) 
{ 
    mod <- match.arg(mean.function)
    call <- match.call(expand.dots=FALSE)
    call$... <- NULL
    subset <- if (missing(subset)) NULL else call$subset
    g <- factor(eval(substitute(group),data))
    gname <- substitute(group)  
    gname <- as.character(if(length(gname) == 1) gname else gname[3])
    if (mod != "pod") { 
      if (mod == "common") {
       call$group <- call$mean.function <- NULL
       call[[1]] <- as.name("lm")
       ans <- eval(call,parent.frame())} else
     if (mod == "parallel") {
       assign(gname,g)
       call$formula <- update(as.formula(call$formula),
                       as.formula(paste("~.+",gname,sep="")))
                        #as.formula("~.+g"))
       call$group <- call$mean.function <- NULL
       call[[1]] <- as.name("lm")
       ans <- eval(call)} else
     if (mod == "general") {
       assign(gname,g)
       call$formula <- update(as.formula(call$formula),
                        as.formula(paste("~(.)*",gname,sep="")))
       call$group <- call$mean.function <- NULL
       call[[1]] <- as.name("lm")
       ans <- eval(call)}
     ans$pod.mean.function <- mod
     ans$group <- structure(if (!is.null(subset)) g[eval(subset)] else g,  
                            name=gname)
     class(ans) <- c("pod.lm","lm")
     ans} else
    { # POD model
    l<-nlevels(g) 
    cl <- call
    cl$group <- cl$mean.function <- NULL
    cl[[1]] <- as.name("lm")
    fit1 <- eval(cl)
    coef1 <- coef(fit1)
    p1 <- predict(fit1) - coef(fit1)[1]
# If subset is not NULL, find the subscripts of the cases used:
    if(!is.null(call$subset)){
       rows <- row.names(model.matrix(update(fit1,subset=NULL))) %in% 
               row.names(model.matrix(fit1))
       } else {
       rows <- rep(TRUE,dim(model.matrix(fit1))[1])}
       bign <- length(rows)
       ans <- rep(NA,bign)
       ans[rows] <- p1
       p1 <- ans
# update the formula in cl to fit parallel within group regression
    cl[[2]] <- update(as.formula(formula),as.formula("~ p1*g"))
# update cl$data
    if (is.null(cl$data)){
       cl$data <- data.frame(p1=p1,g=g)} else {
       cl$data <- data.frame(data,p1=p1,g=g)}
# fit parallel lines model  
    fit2 <- eval(cl)
    coef2 <- coef(fit2)
# group.data.subset excludes the subset
    group.data.subset <- 
              model.matrix(fit2$terms,fit2$model,contrasts)[,2+1:(l-1)]
    group.data.subset <- data.frame(group=group.data.subset)
# again, fix for non-null subset
    if(!is.null(subset)){
       group.data <- data.frame(matrix(NA,nrow=bign,
          ncol=dim(group.data.subset)[2]))
       group.data[rows,] <- group.data.subset} else {
       group.data <- group.data.subset}    
# construct new data with X, y, and G
    formula<-as.formula(formula)
    y.name<-as.character(attr(fit1$terms,"predvars")[2])
    x.name<-attr(fit1$terms,"term.labels")
    g.name<-paste(gname,levels(g)[2:l],sep="")
    names(group.data) <- g.name
    ncols <- length(x.name)+length(y.name)+length(g.name)
    p<-length(x.name)
# generate model
    para.name<-c("eta0", "eta1")
    form1<-paste("eta1 *", x.name[1])
    for(i in 2:p) {
       eta<-paste("eta", i, sep="")
       form1<-paste(form1, "+", eta, "*", x.name[i])
       para.name<-c(para.name, eta)  
    } 
   form2<-paste("eta0 +", form1)
   for(i in 1:(l-1)) {
      G<-g.name[i]
      th0<-paste("th0", i+1, sep="")
      th1<-paste("th1", i+1, sep="")
      form2<-paste(form2, "+", G, "* (", th0, "+", th1, "* (", form1, "))") 
      para.name<-c(para.name, c(th0, th1))
   } 
# New June 6, 2005
   if (is.null(weights))
      {form<-as.formula(paste(y.name, "~", form2))} else
      {wts <- substitute(weights)
       form<-as.formula(paste("sqrt(",wts,")*",y.name,"~",
             "sqrt(",wts,")*(",form2,")", sep=""))}
# End change
   start <- c(coef2[1],coef2[2]*coef1[-1],coef2[3:(l+1)],
                coef2[(l+2):(2*l)]/coef2[2])
   names(start) <- c(paste("eta",0:(length(x.name)),sep=""),
                     paste("th0",2:l,sep=""),paste("th1",2:l,sep=""))
   if (is.null(call$data)){
        for (j in 1:length(g.name)) assign(g.name[j],group.data[,j])
        obj.nl<-nls(form, start=start, subset=rows,
            control=control, ...)} else {
        obj.nl<-nls(form, cbind(data,group.data), start=start, 
            subset=rows, control=control, ...)}  
# evaluate linear part, but only for the subset of cases used
   eta.name<-para.name[2:(p+1)]
   eta.value<-coef(obj.nl)[2:(p+1)]
   for(i in 1:p)
      assign(eta.name[i], eta.value[i])
   form1.expr<-parse(text=form1)
   envr <-
     if (is.null(call$data) & is.null(subset)){group.data} else{
      if(is.null(call$data)) group.data.subset else {
       if(is.null(subset)) as.data.frame(cbind(data,group.data)) else{
        as.data.frame(cbind(data[rows,],group.data.subset))}}}
   if (!is.null(subset)) envr <- as.data.frame(envr[rows,])
   linear.part<-eval(form1.expr, envr)
   ans <- NULL
   ans$nls.fit <- obj.nl
   ans$linear.part <- linear.part
   ans$group <- if (!is.null(subset)) g[rows] else g
   ans$call <- call
   class(ans)<-c("pod")
   ans
}}

podnls.fit <-       function(x){ x$nls.fit }
print.pod<-        function(x, ...){ print(podnls.fit(x)) }
summary.pod <-     function(object, ...){ summary(podnls.fit(object),...) }
coef.pod <-        function(object, ...){ coef(podnls.fit(object),...)}
deviance.pod <-    function(object, ...){ deviance(podnls.fit(object),...)} 
vcov.pod <-        function(object, ...){ vcov(podnls.fit(object),...)}
residuals.pod <-   function(object,...){ residuals(podnls.fit(object))}
formula.pod<-      function(x, ...){      formula(x$call)}
fitted.pod <-      function(object, ...){ fitted(podnls.fit(object),...)}
podresponse<-     function(object, ...){ residuals(object)+fitted(object)}
df.residual.pod <- function(object,...){ length(resid(object))-length(coef(object))}
predict.pod     <- function(object,...){ predict(podnls.fit(object))}

# Plot one dimensional models by group
plot.pod.lm <- function(x, colors=rainbow(nlevels(x$group)),
      pch=1:nlevels(x$group),key=FALSE,identify=FALSE,
      xlab="Linear Predictor", ylab=as.character(c(formula(x)[[2]])),...) {
  mean.function <- x$pod.mean.function
  if(mean.function == "general") stop("No 2D plot for the general pod model")
  g1 <- x$group
  g1.name <- attr(x$group,"name")
  levels(g1) <- 1:nlevels(x$group)
  g1 <- as.numeric(as.character(g1))
  gloc <- match("group",names(x$call))
  gname <- as.character(x$call[[gloc]])
# common regressions
  if(mean.function == "common"){
    plot(predict(x),x$model[,1],ylab=ylab,pch=pch[g1],col=colors[g1],
    xlab=paste(xlab, ", ignore groups", sep=""),...)
    abline(lm(x$model[,1]~predict(x)))}
# parallel regressions
  if(mean.function == "parallel"){
      c2 <- coef(x)
      xp <-0
      for (j in 2:(length(c2)-nlevels(x$group)+1)) 
         xp <- xp +c2[j] * model.matrix(x)[,j]
      plot(xp,x$model[,1],pch=pch[g1],col=colors[g1],
       ylab=paste(ylab, ", Groups = ", g1.name, sep=""),
       xlab=paste(xlab,  ", parallel mean function", sep=""), ...)
      for (j in 1:nlevels(x$group))
       abline(if(j==1) c2[1] else c2[1]+c2[length(c2)-nlevels(x$group)+j],1,
           lty=j,col=colors[j])}
# key
  if (class(key) == "logical") {
   if (key == TRUE) {
      print("Click mouse on plot to locate the key, or press Escape")
      loc <-locator(n=1) 
      legend(loc[1],loc[2], legend = as.character(levels(x$group)),
           lty=1:nlevels(x$group),col=colors[1:nlevels(x$group)],
           pch=pch[1:nlevels(x$group)])}}
   else { 
      loc <- key
      legend(loc[1],loc[2], legend = as.character(levels(x$group)),
           lty=1:nlevels(x$group),col=colors[1:nlevels(x$group)],
           pch=pch[1:nlevels(x$group)])}
# identify
  if(identify == TRUE){
      identify(xp,x$model[,1],row.names(x$model))}
  invisible()}

plot.pod<-function(x, colors=rainbow(nlevels(x$group)),
  pch=1:nlevels(x$group),key=FALSE,identify=FALSE,
  xlab="Linear Predictor", ylab=as.character(c(formula(x)[[2]])),...)
{ 
   yp<-podresponse(x)
   group <- x$group
   gloc <- match("group",names(x$call))
   xp<-x$linear.part
   g1 <- x$group
   levels(g1) <- 1:nlevels(group)
   g1 <- as.numeric(as.character(g1))
   gname <- as.character(x$call[[gloc]])
   plot(xp,yp,pch=pch[g1],col=colors[g1],
    ylab=paste(ylab, ", Groups = ", as.character(x$call$group), sep=""),
    xlab=paste(xlab,  ", pod mean function",sep=""), ...)
      for (j in 1:nlevels(group))
        {abline(lm(yp~xp,subset=g1==j),lty=j,col=colors[j])}
   if (class(key) == "logical") {
   if (key == TRUE) {
      print("Click mouse on plot to locate the key, or press Escape")
      loc <-locator(n=1) 
      legend(loc[1],loc[2], legend = as.character(levels(group)),
           lty=1:nlevels(group),col=colors[1:nlevels(group)],
           pch=pch[1:nlevels(group)])}} else
    { 
      loc <- key
      legend(loc[1],loc[2], legend = as.character(levels(group)),
           lty=1:nlevels(group),col=colors[1:nlevels(group)],
           pch=pch[1:nlevels(group)])}
   invisible()
}

anova.pod <- 
function (object, scale = 0, test = "F", ...) 
{
    m1 <- update(object,mean.function="common")
    m2 <- update(object,mean.function="parallel")
    m4 <- update(object,mean.function="general")
    objects <- list(m1,m2,object,m4)
    resdf  <- as.numeric(lapply(objects, df.residual))
    resdev <- as.numeric(lapply(objects, deviance))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, 
        -diff(resdev)))
    variables <- c("1: common","2: parallel","3: pod","4: pod + 2fi")
    dimnames(table) <- list(variables, c("Res.Df", "RSS", "Df", 
        "Sum of Sq"))
    title <- paste("POD Analysis of Variance Table for ",
            deparse(formula(objects[[1]])[[2]]),", grouped by ",
            as.character(object$call$group),"\n" ,sep="")
    topnote <- c(paste("1: ", deparse(formula(objects[[1]])),sep=""),
                 paste("2: ", deparse(formula(objects[[2]])),sep=""),
                 paste("3: ", deparse(formula(object$nls.fit)),sep=""),
                 paste("4: ", deparse(formula(objects[[4]])),sep=""))
    if (!is.null(test)) {
        bigmodel <- order(resdf)[1]
        scale <- if (scale > 0) 
            scale
        else resdev[bigmodel]/resdf[bigmodel]
        table <- stat.anova(table = table, test = test, scale = scale, 
            df.scale = resdf[bigmodel], n = length(objects[bigmodel$residuals]))
    }
    structure(table, heading = c(title,topnote), class = c("anova", 
        "data.frame"))
}


#######################################################################
# Chapter 7, Transformations
#######################################################################
# The function powtran computes power transformations for 3 families:

powtran <- function(U,lambda,family,modified){UseMethod("powtran")}

powtran.default <- function(U, lambda, family="box.cox",modified=TRUE){
 if (any(class(U) == "factor")) stop("Attempt to transform a factor")
 f <- parse(text=paste("tran.family.",family,sep=""))
 eval(f)(U,lambda,modified=modified)}
 
powtran.matrix <- function(U, lambda, family="power",modified=FALSE){
 Z <- U
 if(dim(U)[2] != length(lambda))
  stop("The length of lambda must equal the number of colums of data")
 for (j in 1:dim(Z)[2]) 
      Z[,j]<-powtran(U[,j],lambda[j],family,modified=modified)
 if (!is.null(colnames(Z))){
       names <- colnames(Z)
       for (j in 1:(length(names))){
         names[j] <- if(abs(lambda[j]) < 1.e-6)
                       paste("log",names[j],sep="") else
                       if (abs(lambda[j]-1)< .001) names[j] else
                       if (lambda[j] > 0)
                          paste(names[j],round(lambda[j],2),sep=".") else
                          paste(names[j],".minus",round(-lambda[j],2),sep="")} 
     colnames(Z) <- names}
 Z}
 
powtran.bctrans <- function(U, lambda=coef(U), family=U$family,modified=FALSE){
 ans <- powtran(U$X,lambda,family,modified)
 if(!is.null(U$Y)) ans <- cbind(ans,U$Y)
 as.data.frame(ans)}
    
powtran.data.frame <- function(U, lambda, family="power",modified=FALSE){
    powtran.matrix(U,lambda,family,modified)}

tran.family <- function(U,lambda,modified) {UseMethod("tran.family")}

tran.family.default <- function(U,lambda,modified)
                                {tran.family.box.cox(U,lambda,modified)}

# if modified is FALSE, return an unnormalized transfomation of U with one parameter
# lambda.  If modified is TRUE, modify by multiplying by the inverse of the Jacobian
# to assure that the scale is the same for each lambda.
tran.family.box.cox <- function(U,lambda,modified=TRUE) {
  z <- if (abs(lambda) <= 1.e-6) log(U) else ((U^lambda) - 1)/lambda
  if (modified == TRUE) {
    z * (exp(mean(log(U),na.omit=TRUE)))^(1-lambda)} else z
  } 


tran.family.boxcox <- function(U,...) {tran.family.box.cox(U,...)}

# The Yeo-Johnson (2000).  {\it Biometrika}, 87, 954-959.

tran.family.yeo.johnson <- function(U,lambda,modified=TRUE) {
  nonnegs <- U >= 0
  z <- rep(NA,length(U))
  z[nonnegs] <- tran.family.box.cox(U[nonnegs]+1,lambda,modified=FALSE)
  z[!nonnegs] <- tran.family.box.cox(-U[!nonnegs]+1,2-lambda,modified=FALSE)
  if (modified == TRUE)
        z * (exp(mean(log((1 + abs(U))^(2 * nonnegs - 1)),na.rm=TRUE)))^(1 -
            lambda)
    else z
  }

tran.family.YJ <- function(U,...){tran.family.yeo.johnson(U,...)}

# modified arg is ignored in this function.
tran.family.power <- function(U,lambda,modified=FALSE){
  if (abs(lambda) <= 1.e-6) log(U) else U^(lambda)}
tran.family.basic <- function(U,lambda,modified=FALSE){
  if (abs(lambda) <= 1.e-6) log(U) else U^(lambda)}

####################################################################
inv.tran.plot<- function(x,y,lambda=c(-1,0,1),lty=1:(1+length(lambda)),
        col=rainbow(length(lambda)+1),xlab=deparse(substitute(x)),
        ylab=deparse(substitute(y)),family="box.cox",optimal=TRUE,
        key=FALSE,...){
 if (is.factor(x)) stop("Predictor variable may not be a factor")
 if (is.factor(y)) stop("Response variable may not be a factor")
 if (optimal){opt <- inv.tran.estimate(x,y,family=family)
              lam <- c(opt$lambda,lambda)} else lam <- lambda
 plot(x,y,xlab=xlab,ylab=ylab,...) 
 rss <- NULL
 new <- seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=100)
 for (j in 1:length(lam)){
     m1 <- lm(y~powtran(x,lam[j],family=family,modified=FALSE),na.action=na.omit)
     rss <- c(rss,sum(m1$residuals^2))
     lines(new,predict(m1,data.frame(x=new)),lty=lty[j],col=col[j])}
 if (class(key) == "logical") {
    if (key == TRUE) {
      print("Click mouse on plot to locate the key, or press Escape")
      loc <-locator(n=1) 
      legend(loc[1],loc[2], legend = as.character(round(lam,2)),lty=lty,col=col,cex=.75)}}
    else { 
      loc <- key
      legend(loc[1],loc[2], legend = as.character(round(lam,2)),lty=lty,col=col,cex=.75)}
 data.frame(lambda=lam,RSS=rss)
}

inv.tran.estimate <- function(x,y,family="box.cox",...){
  if (is.factor(x)) stop("Predictor variable may not be a factor")
  if (is.factor(y)) stop("Response variable may not be a factor")
   f <- function(lambda,x,y,family){
         predict(lm(y~powtran(x,lambda,family=family,modified=TRUE),
                    na.action=na.omit))}   
  lhat <- optimize(f = function(lambda) sum((y-f(lambda,x,y,family))^2),
                 interval=c(-10,10))
  g <- lm(y~powtran(x,lhat$minimum),na.action=na.omit)             
  n1 <- nls( y ~ b0 + b1*powtran(x,lam,family=family,modified=TRUE),
        start=list(b0 = coef(g)[1], b1=coef(g)[2], lam=lhat$minimum),
        na.action=na.omit,...)
  s1 <- summary(n1)
  list(lambda=s1$parameters[3,1], se=s1$parameters[3,2], 
        RSS = s1$df[2] * s1$sigma^2)}       
     
inverse.response.plot <- function(m,lambda=c(0,1),maxiter=100,
    xlab=NULL, ...) 
    UseMethod("inverse.response.plot")

inverse.response.plot.default <- function(m,lambda=c(0,1),
   maxiter=100,xlab=NULL,...) {
   call <- match.call(expand.dots=FALSE)
   if (!is.null(call$xlab)) {
      xlab <- call$xlab
      call$xlab <- NULL} else
      xlab <- names(m$model)[1]
   yhat <- predict(m)
   y <- m$model[,1]
   inv.tran.plot(y,yhat,lambda=lambda,xlab=xlab,...)
}
   
# The following is not used
#inverse.response.plot.formula <- 
#     function(m,data= NULL, subset, na.action=na.omit, ...)
#{
#  formula <- m
#  m1 <- match.call(expand.dots=FALSE)
#  m1[[1]] <- as.name("lm")
#  inverse.response.plot(eval(m1), ...)}

inv.res.plot <- function(m, ...) UseMethod("inverse.response.plot")
  
##########################################################################
#  This is a modification of the function box.cox.powers      
# last modified 15 April 03 by J. Fox
# Here are the differences:
#    1.  Allows X to have missing values.
#    2.  The family argument allows use of other than box.cox.  only
#        yeo.johnson is implemented, but other families can be added by
#        writing a method tran.family.nameoffamily(X,lambda,modified=TRUE).
#    3.  Uses powtran to get the power transformations.
#    4.  New plot method uses pairs to draw a scatterplot matrix of transformations.
#    5.  New 
# (with bug fixes by S. Weisberg)
# reorganized and expanded by S. Weisberg 4/18/04

bctrans<-
function (formula, data = NULL, subset, na.action=na.omit, ...) 
{
    mf <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(mf$data, parent.frame()))) 
        mf$data <- as.data.frame(data)
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mf <- na.action(mf) 
    y <- model.response(mf, "numeric")
    if (!is.null(y)){
      y <- matrix(y,ncol=1)
      colnames(y) <- as.character(attr(attr(mf,"terms"),"variables"))[2]
      mf <- mf[,-1]} 
    if (any(sapply(1:dim(mf)[2],function(j) is.factor(mf[[j]]))))
       stop("At least one specified term is a factor.")
    bctrans1(mf, Y=y, ..., call=match.call(expand.dots=TRUE))
}

bctrans1 <-function(X, Y=NULL, start=NULL, family="box.cox",call=NULL,...){
    modified.power<-function(x, lambda, family, gm){
        powtran(x,lambda,family=family,modified=FALSE) * gm^(1-lambda)}
    get.gm <- function(X,family){
        gm <- powtran(X,0,family=family,modified=TRUE)/
              powtran(X,0,family=family,modified=FALSE)
        gm[!is.nan(gm)][1]}
    neg.kernel.profile.logL<-function(X, lambda, family, gm){
        for (j in 1:ncol(X)){
            X[,j]<-modified.power(X[,j],lambda[j],family, gm[j])
            }
        (nrow(X)/2)*log(((nrow(X)-1)/nrow(X))*det(var(X)))
        }
    univ.neg.kernel.logL <- function(x, lambda, family, gm){
        x <- modified.power(x, lambda, family, gm)
        (length(x)/2)*log(((length(x)-1)/length(x))*var(x))
        }
    result<- list()
    X<-as.matrix(X) 
    result$X <- X
    result$Y <- Y
    nc <- ncol(X) 
    if(any(X<=0) & family == "box.cox") 
      stop("All values must be > 0; use family=\"yeo.johnson\"")
    gm<-apply(X, 2, function(X){ get.gm(X,family=family)})
    if (is.null(start)) {
        start <- rep(1, nc)
        for (j in 1:nc){
            res<- suppressWarnings(optimize(
                f = function(lambda) univ.neg.kernel.logL(x=X[,j], 
                                       lambda=lambda, family=family, gm=gm[j]),
                lower=-50, upper=+50))
            start[j] <- res$minimum
            }
        }  
    res<-optim(start, neg.kernel.profile.logL, hessian=TRUE, 
                         method="L-BFGS-B", X=X, family=family, gm=gm, ...)
    result$start<-start
    result$family<-family
    result$optim <- res
    if(res$convergence != 0) 
        warning(paste("Convergence failure: return code =",res$convergence))
    class(result)<-"bctrans"
    result$call <- if(is.null(call)) match.call(expand.dots=TRUE) else call
    result$call$call <- NULL
    result
    }

    
lrt.bctrans <-function(object, lrt=NULL, ones=TRUE, zeroes=TRUE){
    modified.power<-function(x, lambda, family){
        powtran(x,lambda,family=family,modified=TRUE)}
    neg.kernel.profile.logL<-function(X, lambda, family){
        for (j in 1:ncol(X)){
            X[,j]<-modified.power(X[,j],lambda[j],family)
            }
        (nrow(X)/2)*log(((nrow(X)-1)/nrow(X))*det(var(X)))
        }
    univ.neg.kernel.logL <- function(x, lambda, family){
        x <- modified.power(x, lambda, family)
        (length(x)/2)*log(((length(x)-1)/length(x))*var(x))
        }
    result <- NULL
    rnames <- NULL
    df <- length(object$optim$par)
    if (zeroes==TRUE){
      LR<- 2*(neg.kernel.profile.logL(object$X,rep(0,df),object$family)-object$optim$value)
      result <- rbind(result,c(LR,df,1-pchisq(LR,df)))
      rnames <- c(rnames,"LR test, all lambda equal 0")}
    if (ones==TRUE){
      LR<-2*(neg.kernel.profile.logL(object$X,rep(1,df),object$family)-object$optim$value)
      result <- rbind(result,c(LR,df,1-pchisq(LR,df)))
      rnames <- c(rnames,"LR test, all lambda equal 1")}
    if (!is.null(lrt)) {
        for (i in 1:length(lrt)){
            if (length(lrt[[i]]) != df) 
                stop(paste("hypothesis", i, "that powers =", lrt[[i]], "does not have", df, "lrt"))
            LR<-2*(neg.kernel.profile.logL(object$X,lrt[[i]],object$family)-object$optim$value)
            result <- rbind(result,c(LR,df,1-pchisq(LR,df)))
            rnames <- c(rnames,
                         paste("LR test, lambda =",paste(lrt[[i]],collapse=" ")))}
        }
    rownames(result) <- rnames
    colnames(result) <- c("LRT", "df", "p-value")
    data.frame(result)
    }
      
plot.bctrans <- function(x,y=x$Y,lambda=coef(x),family=x$family,
                 plot=pairs,...){ 
    Z <- powtran(x,lambda,family)
    plot(Z,...)}
    
print.bctrans<-function(x, ...) {
   cat("Estimated transformation parameters \n")
   print(coef(x))
   invisible(x)}
      
summary.bctrans<-function(object, digits=4,...){
    one<-1==length(object$optim$lambda)
    cat(paste(object$family, 
       (if(one) "Transformation to Normality" else 
                "Transformations to Multinormality"),"\n\n"))
    lambda<-coef(object)
    stderr<-sqrt(diag(vcov(object)))
    df<-length(lambda) 
    result<-cbind(lambda,stderr,lambda/stderr,(lambda-1)/stderr)
    rownames(result)<-colnames(object$X)
    colnames(result)<-c("Est.Power","Std.Err.","Wald(Power=0)","Wald(Power=1)")
    if (one)rownames(result)<-""
    print(round(result,digits))
    print(lrt.bctrans(object,...))
    invisible(object)
    }

coef.bctrans <- function(object,...){
    lambda <- object$optim$par
    names(lambda) <- colnames(object$X)
    lambda
    }
    
vcov.bctrans <- function(object,...){solve(object$optim$hessian)}

    
################################################################################
#    Chapter 8
################################################################################
# Test for curvature in a residual plot (Chapter 8)
# Residual plots, and curvature tests
  
resid.curv.test <- function(m,varname) {
 if(varname == "tukey") tukey.nonadd.test(m) else {
  if(is.na(match(varname, attr(m$terms,"term.labels"))))
     stop(paste(varname,"is not a term in the mean function")) else {
     form <- update.formula(formula(m),
           as.formula(paste("~ I(",varname,"^2) +.",sep="")))
     mup <- update(m,form)
     if(mup$rank > m$rank) summary(mup)$coef[2,3:4]  else c(NA,NA)
     }}}
     
tukey.nonadd.test <- function(m){ 
# yhat^2 must be of the length of the data without deleting any cases.
  horiz <- predict(m,model.frame(update(m,subset=NULL)))
  assign("zZ5481f",horiz^2,pos=sys.frame()) # create variable in system frame
  mup <- update(m,~zZ5481f+.)               # to make it visible to update
  rm(zZ5481f,envir=sys.frame())             # delete from system frame
  if(mup$rank > m$rank){
   ans <- summary(mup)$coef[2,3:4] 
   ans[2] <- pnorm(-abs(ans[1]))*2
   ans} else c(NA,NA)
   }
  
resplot <- function(m,varname="tukey",type="pearson",
                    plot=TRUE,add.quadratic=TRUE,...){
 string.capitalize <- function(string) {
     paste(toupper(substring(string,1,1)),substring(string,2),sep="")}
 col <- match(varname,names(m$model))
 if(is.na(col) && varname != "tukey")
   stop(paste(varname,"is not a term in the mean function"))
 horiz <- if(varname == "tukey") predict(m) else m$model[[col]]
 
 lab <- if(varname == "tukey") {"Fitted values"} else varname
 if(plot==TRUE){ 
  plot(horiz,residuals(m,type=type),xlab=lab,
            ylab=paste(string.capitalize(type),"Residuals",...))
  abline(h=0,lty=2)
  if(class(horiz) != "factor") {
    if(add.quadratic==TRUE){
        new <- seq(min(horiz),max(horiz),length=200)
        lm2 <- lm(residuals(m,type=type)~poly(horiz,2))
        lines(new,predict(lm2,list(horiz=new)),lty=3)
        }
    resid.curv.test(m,varname)} else c(NA,NA)}}

residual.plots <- function(m, ...){UseMethod("residual.plots")}
residual.plots.lm <- function(m,tukey=TRUE,exclude=NULL,plot=TRUE,
     layout=NULL,ask,...){
  term.labels <- attr(m$terms,"term.labels") # this is a list
  term.classes <- attr(m$terms,"dataClasses")
  nt <- length(term.labels) # number of terms excluding intercept
  if(is.null(layout)){
   layout <- switch(min(nt+1,9),c(1,1),c(1,2),c(2,2),c(2,2),c(3,2),c(3,2),
                                c(3,3),c(3,3),c(3,3))}
  nr <- 0
  op<-par(no.readonly=TRUE)
  ask <- if(missing(ask) || is.null(ask))
    prod(layout)<nt-length(exclude) else ask
  on.exit(par(op))
  if(prod(layout) > 1)
    par(mfrow=layout,mai=c(.6,.6,.1,.1),mgp=c(2,1,0),
        cex.lab=1.0,cex=0.7,ask=ask) 
    else par(mfrow=layout,ask=ask)   
  ans <- NULL
  for (j in 1:nt){ 
   if(is.na(match(j,exclude))){ 
     nr <- nr+1
     ans <- rbind(ans,resplot(m,term.labels[j],plot=plot,...))
     row.names(ans)[nr] <- term.labels[j]
    }}
  # Tukey's test
  if (tukey == TRUE){
   ans <- rbind(ans,resplot(m,"tukey",plot=plot,...))
   row.names(ans)[nr+1] <- "Tukey test"
   ans[nr+1,2] <- 2*pnorm(abs(ans[nr+1,1]),lower.tail=FALSE)}
  dimnames(ans)[[2]] <- c("Test stat", "Pr(>|t|)")
  ans}

#############################################
# marginal model plots
#############################################
marginal.model.plot <- function(...){mmp(...)}
mmp <- function(object, ...){UseMethod("mmp")}

mmp.lm <- function(object,u=predict(object),mean=TRUE,sd=FALSE,
         label=deparse(substitute(u)),degree=1,span=2/3,
         colors=c("blue","red"),...){
  if(label=="predict(object)"){label <- "Fitted values"}
  plot(u,object$model[,1],xlab=label,ylab=colnames(object$model[1]),...)
  loess.y <- loess(object$model[,1]~u, degree=degree, span=span)
  loess.yhat <- loess(predict(object) ~ u, degree=degree,span=span)
  new <- seq(min(u),max(u),length=200)
  if(mean==TRUE) {
   lines(new,predict(loess.y,data.frame(u=new)),lty=1,col=colors[1])
   lines(new,predict(loess.yhat,data.frame(u=new)),lty=2,col=colors[2])}
  if(sd==TRUE) { # add \pm sd lines
   loess.y.var <- loess(residuals(loess.y)^2~u,degree=degree,span=span)
   lines(new, predict(loess.y,data.frame(u=new)) + 
              sqrt(predict(loess.y.var,data.frame(u=new))), lty=1,col=colors[1])
   lines(new, predict(loess.y,data.frame(u=new)) - 
              sqrt(predict(loess.y.var,data.frame(u=new))), lty=1,col=colors[1]) 
   loess.yhat.var <- loess(residuals(loess.yhat)^2~u,
              degree=degree,span=span)
   s2 <- summary(object)$sigma^2
   lines(new, predict(loess.yhat,data.frame(u=new)) + 
              sqrt(s2+predict(loess.yhat.var,data.frame(u=new))), lty=2,col=colors[2])
   lines(new, predict(loess.yhat,data.frame(u=new)) - 
              sqrt(s2+predict(loess.yhat.var,data.frame(u=new))), lty=2,col=colors[2])} 
  }
  
mmp.glm <- function(object,u=predict(object),mean=TRUE,sd=FALSE,
         label=deparse(substitute(u)),degree=1,span=2/3,
         colors=c("blue","red"),...){
  fr.mmp <- function(family,x) {
   if(family == "binomial") pmax(0, pmin(1, x)) else
    if(family == "poisson") pmax(0, x) else
     if(family == "gamma") pmax(0, x) else x}
  if(label=="predict(object)"){label <- "Linear Predictor"}
  response <- object$model[,1]
  fam <- object$family$family
  if (is.matrix(response)) response <- response[,1]/apply(response,1,sum)
  plot(u,response,xlab=label,ylab=colnames(object$model[1]),...)
  loess.y <- loess(response~u, degree=degree, span=span)
  loess.yhat <- loess(predict(object,type="response") ~ u, degree=degree,span=span)
  new <- seq(min(u),max(u),length=200)
  pred.loess.y <- fr.mmp(fam,predict(loess.y,data.frame(u=new)))
  pred.loess.yhat <- fr.mmp(fam,predict(loess.yhat,data.frame(u=new)))
  if(mean==TRUE) {
   lines(new,pred.loess.y,   lty=1,col=colors[1])
   lines(new,pred.loess.yhat,lty=2,col=colors[2])}
  if(sd==TRUE) { # add \pm sd lines
   loess.y.var <- loess(residuals(loess.y)^2~u,degree=degree,span=span)
   pred.loess.y.var <- pmax(0,predict(loess.y.var,data.frame(u=new)))
   lines(new, fr.mmp(fam,pred.loess.y + sqrt(pred.loess.y.var)), lty=1,col=colors[1])
   lines(new, fr.mmp(fam,pred.loess.y - sqrt(pred.loess.y.var)), lty=1,col=colors[1])  
   loess.yhat.var <- loess(residuals(loess.yhat)^2~u,degree=degree,span=span)
   pred.loess.yhat.var <- pmax(0,predict(loess.yhat.var,data.frame(u=new))) 
   varfun <- summary(object)$dispersion * 
               object$family$variance(predict(object,type="response"))/
               if(!is.null(object$prior.weights)) object$prior.weights else 1  
   loess.varfun <- loess(varfun~u,degree=degree,span=span)
   pred.loess.varfun <- pmax(0,predict(loess.varfun,data.frame(u=new)))
   sd.smooth <- sqrt(pred.loess.yhat.var + pred.loess.varfun)
   lines(new, fr.mmp(fam,pred.loess.yhat + sd.smooth), lty=2,col=colors[2])
   lines(new, fr.mmp(fam,pred.loess.yhat - sd.smooth), lty=2,col=colors[2])} 
  }
  
mmps <- function(object,exclude=NULL,layout=NULL,
          ask,...){
  predvars <- attr(object$terms,"predvars") # this is a list
  nt <- length(predvars)-1 # number of terms excluding intercept
  factors <- which(sapply(2:dim(object$model)[2],function(j) 
                    is.factor(object$model[[j]])))
  if (length(factors)>0)warning("Factors were skipped")
  exclude <- c(exclude,factors)
  if(is.null(layout)){
   layout <- switch(min(nt+1-length(exclude),9),
                           c(1,1),c(1,2),c(2,2),c(2,2),c(3,2),c(3,2),
                           c(3,3),c(3,3),c(3,3))}
  op<-par(no.readonly=TRUE)
  ask <- if(missing(ask) || is.null(ask))
    prod(layout)<nt-length(exclude) else ask
  on.exit(par(op))
  if(prod(layout) > 1)
    par(mfrow=layout,mai=c(.6,.6,.1,.1),mgp=c(2,1,0),cex.lab=1.0,cex=0.7,ask=ask) 
    else par(mfrow=layout,ask=ask) 
  for (j in 3:(nt+1)){
      if(is.na(match(j-2,exclude))){ 
       mmp(object,object$model[,j-1],label=deparse(predvars[[j]]),...)}}
  mmp(object,,...)
  }
 


################################################################
# Chapter 9
################################################################

############################################
# influence index plot 

inf.index <- function(m,cooks.distance,rstudent,outlier.t.test,leverages,...)
         {UseMethod("inf.index")}

inf.index.lm <- function(m,cooks.distance=TRUE,rstudent=TRUE,
                           outlier.t.test=TRUE,leverages=TRUE,...) {
sel <- which(c(cooks.distance,rstudent,outlier.t.test,leverages))
names <- c("Cook's distance", "Studentized residuals", 
           "Bonferroni p-value", "Leverage")
# check for row.names, and use them if they are numeric.
xaxis <- sapply(row.names(m$model),function(x) eval(parse(text=x)))
if (any (sapply(xaxis,is.character))) xaxis <- 1:length(xaxis)
op <- par(mfrow=c(length(sel),1),mai=c(.6,.6,.2,.2),mgp=c(2,1,0))
on.exit(par(op))
for (j in sel){
 y <- switch(j,cooks.distance(m),rstudent(m),
          outlier.t.test(m,order=FALSE,bound=1.01)$Bonf.pvals,
          hatvalues(m))
 plot(xaxis,y,type="b",xlab="Index",ylab=names[j],...)}
invisible()
}

outlier.t.test <- function(m,order=TRUE,bound=1) {UseMethod("outlier.t.test")}
outlier.t.test.default <- function(m,order=TRUE,bound=1){
  res <- rstudent(m)
  df <- m$df.residual - 1
  pvals <- pmin(length(res)*2*(1-pt(abs(res),df)),1)
  or <- if(order==TRUE) order(pvals) else seq(length(res))
  d<-data.frame(tvalue=res[or],Bonf.pvals=pvals[or])
  d[d$Bonf.pvals<bound,]
  }
