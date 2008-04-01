# residual plot tests
data(forbes)
X <- X1 <- forbes$X1 <-forbes$Temp
X1[3] <- forbes$X1[3] <- NA
Y <- Y1 <- forbes$Y1 <- forbes$Pressure
Y1[4] <- forbes$Y1[4] <- NA

m1 <- lm(log(Pressure)~log(Temp),forbes)
m2 <- lm(log(Y)~log(X))
m3 <- lm(log(Y1)~log(X1),forbes)
m4 <- lm(log(Y1)~log(X1))
m5 <- update(m3,subset=-10)
m6 <- update(m4,subset=-10)
tukey.nonadd.test(m1)
tukey.nonadd.test(m2)
tukey.nonadd.test(m3)
tukey.nonadd.test(m4)
tukey.nonadd.test(m5)
tukey.nonadd.test(m6)

data(sleep1)
s1 <- lm(TS ~log(BodyWt)+log(BrainWt)+log(GP)+log(Life)+P+factor(D),sleep1)
residual.plots(s1)

tukey.nonadd.test <- function(m){
  envir <- environment(formula(m))
	dd <- eval(m$call$data, envir)
	subs <- eval(m$call$subset, envir)
	wgts <- eval(m$call$weights, envir)
	naa <- m$call$na.action
  preds.sq <- fitted(update(m,na.action=na.exclude))^2
  if(!is.null(dd)) dd <- data.frame(dd, preds.sq=predict(m, dd)^2)
	uf <- update.formula(formula(m$terms), ~ . + preds.sq)
	environment(uf) <- environment(NULL)
  if(is.null(dd) & !is.null(subs)) { 
    out <- update(m,subset=NULL,method="model.frame")
    out <- rep(NA,dim(out)[1]+length(attr(out,"na.action")))
    out[subs] <- preds.sq
    preds.sq <- out} 
	mup <- if(is.null(naa) & is.null(dd))
      lm(uf, subset=subs, weights=wgts)
   else if(is.null(naa)) lm(uf, data=dd, subset=subs, weights=wgts)
   else if(is.null(dd))  lm(uf, subset=subs, weights=wgts)
   else lm(uf, data=dd, subset=subs, weights=wgts,na.action=naa)
  if(mup$rank > m$rank){
   ans <- summary(mup)$coef[,3:4]
   ans <- ans[match("preds.sq",rownames(ans)),]
   ans[2] <- pnorm(-abs(ans[1]))*2
   ans} else c(NA,NA)
   }


lm1 <-
function (formula, data, subset, weights, na.action, method = "qr",
    model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
    contrasts = NULL, offset, ...)
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame")
        return(mf)
}