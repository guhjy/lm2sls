fit2sls <- function(y, X, Z, wt=NULL, singular.ok=FALSE, qr=TRUE, ...){

  # handle NAs
  
  na.action <- which(rowSums(is.na(cbind(y, X, Z, wt))) > 0)
  if (length(na.action) > 0) {
    class(na.action) <- "exclude" # to be used for fits and residuals
    warning(length(na.action), " cases with missing values deleted")
    y <- y[-na.action]
    X <- X[-na.action, , drop=FALSE]
    Z <- Z[-na.action, , drop=FALSE]
    if (!is.null(wt)) wt <- wt[-na.action]
  }
  else na.action <- NULL
  
  # stage 1 regression
  
  stage1 <- suppressWarnings(lsfit(Z, X, wt=wt, intercept=FALSE))
  rk.1 <- stage1$qr$rank
  if (rk.1 < ncol(Z)){
    msg <- paste0("rank of stage-1 model matrix = ", rk.1, 
                  " < number of stage-1 coefficients = ", ncol(Z))
    if (singular.ok) warning(msg) else stop(msg)
  }
  
  # stage 2 regression
  
  X.fit <- X - stage1$residuals
  colnames(X.fit) <- colnames(X)
  stage2 <- suppressWarnings(lsfit(X.fit, y, wt=wt, intercept=FALSE))
  rk.2 <- stage2$qr$rank
  if (rk.2 < ncol(X.fit)){
    msg <- paste0("rank of stage-2 model matrix = ", rk.2, 
                  " < number of stage-2 coefficients = ", ncol(X.fit))
    if (singular.ok) warning(msg) else stop(msg)
  }
  fitted <- X %*% stage2$coef
  residuals <- y - fitted
  p <- ncol(X)
  n <- nrow(X)
  df.res <- n - p
  sigma2 <- (if (is.null(wt)) sum(residuals^2) else sum(wt*residuals^2))/df.res
  vcov <- sigma2*chol2inv(stage2$qr$qr)
  rownames(vcov) <- colnames(vcov) <- names(stage2$coef)
  
  list(
    n              = n,
    p              = p,
    q              = ncol(Z),
    qr             = if (qr) stage2$qr else NULL,
    rank           = rk.2,
    qr.1           = if (qr) stage1$qr else NULL,
    rank.1         = rk.1,
    coefficients   = stage2$coef,
    coefficients.1 = stage1$coef,
    vcov           = vcov,
    df.residual    = df.res,
    sigma          = sqrt(sigma2),
    residuals      = naresid(na.action, as.vector(residuals)),
    fitted         = napredict(na.action, as.vector(fitted)),
    residuals.1    = naresid(na.action, stage1$residuals),
    fitted.1       = napredict(na.action, X.fit),
    residuals.2    = naresid(na.action, stage2$residuals),
    fitted.2       = napredict(na.action, y - stage2$residuals)
  )
}

lm2sls <- function (formula, instruments=rhs(formula), data, subset, weights, 
                    na.action=getOption("na.action"), contrasts = NULL, 
                    singular.ok=FALSE, model=TRUE, x=TRUE, y=TRUE, qr=TRUE, ...){
  rhs <- function(formula) if (length(formula) == 3) formula[-2] else formula
  combineFormulas <- function(formula1, formula2){
    rhs <- as.character(formula2)[length(formula2)]
    formula2 <- paste("~ . +", rhs)
    update(formula1, formula2)
  }
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) m$data <- as.data.frame(data)
  response.name <- deparse(formula[[2]])
  m$formula <- combineFormulas(formula, instruments)
  m$instruments <- m$contrasts <- m$singular.ok <- NULL
  m[[1]] <- as.name("model.frame")
  mt.model <- mt.instruments <- m
  mt.model$formula <- formula
  mt.instruments$formula <- instruments
  mf <- eval(m, parent.frame())
  mt.model <- attr(eval(mt.model, parent.frame()), "terms")
  mt.instruments <- attr(eval(mt.instruments, parent.frame()), "terms")
  na.act <- attr(mf, "na.action")
  w <- as.vector(model.weights(mf))
  wt.var <- if (!is.null(w)) deparse(substitute(weights)) else NULL
  Z <- model.matrix(instruments, data = mf, contrasts)
  y. <- mf[, response.name]
  X <- model.matrix(formula, data = mf, contrasts)
  result <- fit2sls(y., X, Z, w, singular.ok=singular.ok, qr=qr)
  result <- c(result, list(
    response.name            = response.name,
    formula                  = formula,
    instruments              = instruments,
    model.matrix             = if (x) X else NULL,
    y                        = if (y) y. else NULL,
    model.matrix.instruments = if (x) Z else NULL,
    weights                  = w,
    wt.var                   = wt.var,
    na.action                = na.act,
    call                     = cl,
    contrasts                = attr(X, "contrasts"),
    contrasts.instruments    = attr(Z, "contrasts"),
    xlevels                  = .getXlevels(mt.model, mf),
    xlevels.instruments      = .getXlevels(mt.instruments, mf),
    terms                    = mt.model,
    terms.instruments        = mt.instruments,
    model                    = if (model) mf else NULL
  ))
  class(result) <- c("2sls", "lm")
  result
}

model.matrix.2sls <- function(object, type=c("model", "instruments", "stage2"), ...){
  type <- match.arg(type)
  switch(type,
         model = object$model.matrix,
         instruments = object$model.matrix.instruments,
         stage2 = object$fitted.1)
}

avPlot.2sls <- function(model, ...){
  model$model.matrix <- model.matrix(model, type="stage2")
  NextMethod()
}

vcov.2sls <- function(object, ...) object$vcov

residuals.2sls <- function(object, type=c("model", "stage1", "stage2", "working", 
                                          "deviance", "pearson", "partial"), ...){
  type <- match.arg(type)
  w <- object$weights
  if (is.null(w)) w <- 1
  res <- switch(type,
         working  =,
         model    = object$residuals,
         deviance =,
         pearson  = sqrt(w)*object$residuals,
         stage1   = object$residuals.1,
         stage2   = object$residuals.2,
         partial  = object$residuals + predict(object, type = "terms"))
  naresid(object$na.action, res)
}

fitted.2sls <- function(object, type=c("model", "stage1", "stage2"), ...){
  type <- match.arg(type)
  switch(type,
         model  = napredict(object$na.action, object$fitted),
         stage1 = napredict(object$na.action, object$fitted.1),
         stage2 = napredict(object$na.action, object$fitted.2))
}

print.2sls <- function (x, digits = getOption("digits") - 2, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(coef(x), digits=digits)
  cat("\nResidual standard deviation =", 
      paste0(format(x$sigma, digits=digits), ",  R-squared analog ="), Rsq(x))
  invisible(x)
}

summary.2sls <- function (object, digits = getOption("digits") - 2, vcov.=vcov, ...){
  df <- df.residual(object)
  std.errors <- if (is.function(vcov.)) {
    sqrt(diag(V <- vcov.(object)))
  } else {
    if (!is.matrix(vcov.)) stop("vcov. must be a function or a matrix")
    if (!(ncol(vcov.) == nrow(vcov.)) || ncol(vcov.) != length(coef(object))) 
      stop("vcov. matrix is of the wrong dimensions")
    sqrt(diag(V <- vcov.))
  }
  b <- coef(object)
  t <- b/std.errors
  p <- 2 * (1 - pt(abs(t), df))
  table <- cbind(b, std.errors, t, p)
  rownames(table) <- names(b)
  colnames(table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  intercept <- which(names(b) == "(Intercept)")
  b <- b[-intercept]
  V <- V[-intercept, -intercept]
  F <-  (b %*% solve(V) %*% b) / length(b)
  pval <- pf(F, length(b), df, lower.tail=FALSE)
  result <- list(call=object$call,
                 residuals = summary(residuals(object)), 
                 coefficients = table, digits = digits, s = object$s, 
                 df = df, dfn=length(b), na.action=object$na.action, r2=Rsq(object), 
                 r2adj=Rsq(object, adjust=TRUE), F=F, pval=pval)
  class(result) <- "summary.2sls"
  result
}

print.summary.2sls <- function (x, ...) {
  digits <- x$digits
  cat("Call:\n")
  print(x$call)
  cat("\nResiduals:\n")
  print(round(x$residuals, digits))
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits)
  cat("\nResidual standard deviation =", round(x$s, digits), 
            "on", x$df, "degrees of freedom")
  cat("\nR-squared analog =", paste0(format(x$r2, digits=digits), 
                                     ",  Adjusted R-squared ="), 
      format(x$r2adj, digits=digits))
  cat("\nF =", paste0(format(x$F, digits=digits), " on ", x$dfn, " and ", x$df), 
      "degrees of freedom, p", if (x$pval < .Machine$double.eps) "" else "=", 
      format.pval(x$pval, digits=digits), "\n\n")
  if (!is.null(x$na.action)){
    cat(paste(length(x$na.action), "cases deleted due to NAs\n\n"))
  }
  invisible(x)
}

anova.2sls <- function(object, model.2, s2, dfe, ...){
  if (!inherits(model.2, "2sls")) stop('requires two models of class 2sls')
  s2.1 <- object$s^2
  dfe.1 <- df.residual(object)
  s2.2 <- model.2$s^2
  dfe.2 <- df.residual(model.2)
  SS.1 <- s2.1 * dfe.1
  SS.2 <- s2.2 * dfe.2
  SS <- abs(SS.1 - SS.2)
  Df <- abs(dfe.2 - dfe.1)
  if (missing(s2)){
    s2 <- if (dfe.1 > dfe.2) s2.2 else s2.1
    f <- (SS/Df) / s2
    RSS <- c(SS.1, SS.2)
    Res.Df <- c(dfe.1, dfe.2)
    SS <- c(NA, SS)
    P <- c(NA, 1 - pf(f, Df, min(dfe.1, dfe.2)))
    Df <- c(NA, Df)
    f <- c(NA, f)
    rows <- c("Model 1", "Model 2")
  }
  else{
    f <- (SS/Df) / s2
    RSS <- c(SS.1, SS.2, s2*dfe)
    Res.Df <- c(dfe.1, dfe.2, dfe)
    SS <- c(NA, SS, NA)
    P <- c(NA, 1 - pf(f, Df, dfe), NA)
    Df <- c(NA, Df, NA)
    f <- c(NA, f, NA)
    rows <- c("Model 1", "Model 2", "Error")
  }
  table <- data.frame(Res.Df, RSS, Df, SS, f, P)
  head.1 <- paste("Model 1: ",format(object$formula), "  Instruments:", 
                  format(object$instruments))
  head.2 <- paste("Model 2: ",format(model.2$formula), "  Instruments:", 
                  format(model.2$instruments))
  names(table) <- c("Res.Df", "RSS", "Df", "Sum of Sq", "F", "Pr(>F)")
  row.names(table) <- rows
  structure(table, heading = c("Analysis of Variance", "", head.1, head.2, ""), 
            class = c("anova", "data.frame"))
}

update.2sls <- function (object, formula., instruments., ..., evaluate=TRUE){
  # adapted from stats::update.default()
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) 
    call$formula <- update(formula(object), formula.)
  if (!missing(instruments.)) 
    call$instruments <- update(object$instruments, instruments.)
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) 
    eval(call, parent.frame())
  else call
}

bread.2sls <- function(x, ...){
  x$vcov*x$n/x$sigma^2
}

diagprod <- function(d, X){
  # equivalent to diag(d) %*% X
  if (!is.vector(d)) stop("d is not a vector")
  if (!is.matrix(X)) stop("X is not a matrix")
  if (length(d) != nrow(X)) stop("d and X not conformable")
  d*X
}

estfun.2sls <- function (x, ...) {
  if (x$rank < x$p) stop("second stage model matrix is of deficient rank")
  w <- x$weights
  if (is.null(w)) w <- 1
  diagprod(w*residuals(x), model.matrix(x, type="stage2"))
}


Rsq <- function(model, ...){
  UseMethod("Rsq")
}

Rsq.default <- function(model, adjusted=FALSE, ...){
  SSE <- sum(residuals(model)^2, na.rm=TRUE)
  y <- na.rm(model.response(model.frame(model)))
  SST <- sum((y - mean(y))^2)
  if (adjusted) {
    1 - (SSE/df.residual(model))/(SST/(model$n - 1))
  }
  else {
    1 - SSE/SST
  }
}
