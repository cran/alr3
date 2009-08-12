delta.method <- function(object,...)
  {UseMethod("delta.method")}
  
delta.method.default <- function(object,g,var,...){
   if(!is.character(g)) stop("The argument 'g' must be a character string")
   coefest <- object
   charSubst <- function(x,oldchar) { 
      sapply(x,function(x) paste(unlist(strsplit(x,oldchar)),collapse="_"))}
   func <- g  
# remove characters from g and para names that would cause D to fail 
   g <- charSubst(g,oldchar=":")
   for (c in c(":","\\^"," ")){
          names(coefest) <- charSubst(names(coefest),oldchar=c)}
   g <- parse(text=g)
   q <- length(coefest)
   for(i in 1:q) {assign(names(coefest)[i], coefest[i])}
   est<-eval(g) 
   names(est)<-NULL
# derivative of function g of parameters
   gd<-NULL
   for(i in 1:q) {gd<-c(gd, eval(D(g, names(coefest)[i])))}
# compute se          
   se.est<-as.vector(sqrt(t(gd) %*% var %*% gd))
# output
   data.frame(Estimate=est, SE=se.est, row.names=c(func))
}

### lm, glm, and others with unnamed parameters:
delta.method.lm<-function(object,g,var=vcov,parameterPrefix="b",...)
{
   para <- coef(object)
   var <- if(is.function(var)) var(object) else var
   names(para) <- if ("(Intercept)" %in% names(para))
    paste(parameterPrefix,0:(length(para)-1),sep="") else
    paste(parameterPrefix,1:length(para),sep="")
   delta.method.default(para,g,var)}
   
# nls has named parameters so parameterPrefix is ignored
delta.method.nls<-function(object, g, var=vcov,...)
{
   var <- if(is.function(var)) var(object)
   delta.method.default(coef(object),g,var)}
  
delta.method.drc <-  
function(object,g,var=vcov,...)
{  var <- if(is.function(var)) var(object)
   delta.method.default(coef(object),g,var)
   }
   
delta.method.lmList <- function(object,g,var=vcov,parameterPrefix="b",...){
  out <- NULL
  for (j in 1:length(object)){
   ans <- delta.method(object[[j]],g,var=var,parameterPrefix="b")
   rownames(ans)[1]<-paste(names(object)[j],rownames(ans)[1])
   out <- rbind(out,ans)
   }
   out} 
   
delta.method.nlsList <- function(object,g,var=vcov,...){
  out <- NULL
  for (j in 1:length(object)){
   ans <- delta.method(object[[j]],g,var=vcov)
   rownames(ans)[1]<-paste(names(object)[j],rownames(ans)[1])
   out <- rbind(out,ans)
   }
   out} 
   
delta.method.lme <-  function(object,g,var=vcov,parameterPrefix="b",...)
{  para <- fixef(object)
   var <- if(is.function(var)) var(object)
   names(para) <- if ("(Intercept)" %in% names(para))
    paste(parameterPrefix,0:(length(para)-1),sep="") else
    paste(parameterPrefix,1:length(para),sep="")
   delta.method.default(para,g,var)}
   
delta.method.nlme <-  
function(object,g,var=vcov,...)
{  para <- fixef(object)
   var <- if(is.function(var)) var(object)
   delta.method.default(para,g,var)}
   
delta.method.mer <-  
function(object,g,var=vcov,parameterPrefix="b",...)
{  para <- fixef(object)
   var <- if(is.function(var)) var(object)
   names(para) <- if ("(Intercept)" %in% names(para))
    paste(parameterPrefix,0:(length(para)-1),sep="") else
    paste(parameterPrefix,1:length(para),sep="")
   delta.method.default(para,g,var)}
   
delta.method.multinom<-function(object,g,var=vcov,parameterPrefix="b",...)
{
   out <- NULL
   coefs <- coef(object)
   nc <- dim(coefs)[2]
  for (i in 1:dim(coefs)[1]){
   para <- coefs[i,]
   names(para) <- if ("(Intercept)" %in% names(para))
    paste(parameterPrefix,0:(length(para)-1),sep="") else
    paste(parameterPrefix,1:length(para),sep="")
   ans <- delta.method.default(para,g,var(object)[(i-1)+1:nc,(i-1)+1:nc])
   rownames(ans)[1]<-paste(rownames(coefs)[i],rownames(ans)[1])
   out <- rbind(out,ans)
   }
   out}
   
delta.method.polr<-function(object,g,var=vcov,parameterPrefix="b",...)
{
   para <- coef(object)
   var <- if(is.function(var)) var(object)[1:length(para),1:length(para)]
   names(para) <- if ("(Intercept)" %in% names(para))
    paste(parameterPrefix,0:(length(para)-1),sep="") else
    paste(parameterPrefix,1:length(para),sep="")
   delta.method.default(para,g,var)}
   
