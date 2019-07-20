influence.2sls <- function(model, sigma. = n <= 1e3, type=c("stage2", "both"), ...){

  if (!is.null(model$weights)) stop("weights not supported") #TODO: support weights?
  
  type <- match.arg(type)
  
  Z <- model$model.matrix.instruments
  X <- model$model.matrix
  X.fit <- model$fitted.1
  y <- model$y
  b <- model$coefficients
  res <- model$residuals
  sigma2 <- model$sigma^2
  hatvalues <-  na.omit(hatvalues(model, type=type))
  
  na.action <- model$na.action
  
  rnames <- rownames(X)
  cnames <- colnames(X)
  
  names(hatvalues) <- rnames
  
  rss <- sum(res^2)
  ZtZinv <- solve(crossprod(Z)) #TODO: avoid matrix inversions?
  XtZ <- crossprod(X, Z)
  A <- XtZ %*% ZtZinv %*% t(XtZ)
  Ainv <- solve(A)
  pi <- ZtZinv %*% crossprod(Z, y)
  r <- XtZ %*% ZtZinv %*% t(Z)
  XfXfinv <- solve(crossprod(X.fit))
  
  n <- model$n
  p <- model$p
  dfbeta <- matrix(0, n, p)
  rownames(dfbeta) <- rnames
  colnames(dfbeta) <- cnames
  dffits <- cookd <- rep(0, n)
  sigma <- rep(sqrt(sigma2), n)
  names(dffits) <- names(sigma) <- names(cookd) <- rnames
  
  for (i in 1:n){ #TODO: move this loop to cpp code?
    c <- as.vector(Z[i, ] %*% ZtZinv %*% Z[i, ])
    Xmr <- X[i, ] - r[, i]
    XiAinvXi <- as.vector(X[i, ] %*% Ainv %*% X[i, ])
    XmrAinvXi <- as.vector(Xmr %*% Ainv %*% X[i, ])
    XmrAinvXmr <- as.vector(Xmr %*% Ainv %*% Xmr)
    delta <- 1 - XiAinvXi + XmrAinvXi^2/(1 - c + XmrAinvXmr)
    denom <- (1 - c + XmrAinvXmr)*delta
    h <- Xmr * (1 - XiAinvXi) / denom + X[i, ] * XmrAinvXi / denom
    j <- Xmr * XmrAinvXi / denom - X[i, ]/delta
    g <- h * as.vector((y[i] - Z[i, ] %*% pi) - (y[i] - r[, i] %*% b)) +
      (h + j)*as.vector(y[i] - X[i, ] %*% b)
    dfbeta[i, ] <- - Ainv %*% g
    if (sigma.){
      ss <- rss + as.vector(g %*% Ainv %*% crossprod(X[-i, ]) %*% Ainv %*% g) -
        2 * as.vector(g %*% Ainv %*% t(X[-i, ]) %*% (y[-i] - X[-i, ] %*% b))  - 
        as.vector(y[i] - X[i, ] %*% b)^2
      sigma[i] <- sqrt(ss/(n - p - 1))
    }
    dffits[i] <- as.vector(X[i, ] %*% dfbeta[i, ])/
      (sigma[i] * as.vector(sqrt(X[i, ] %*% XfXfinv %*% X[i, ])))
  }
  
  rstudent <- res/(sigma * sqrt(1 - hatvalues))
  cookd <- (sigma^2/sigma2)*dffits^2/p
  
  result <- list(dfbeta = naresid(na.action, dfbeta), 
                 sigma = naresid(na.action, sigma), 
                 dffits = naresid(na.action, dffits), 
                 cookd = naresid(na.action, cookd), 
                 hatvalues = naresid(na.action, hatvalues), 
                 rstudent = naresid(na.action, rstudent))
  class(result) <- "influence.2sls"
  result
}

rstudent.influence.2sls <- function(model, ...) model$rstudent

rstudent.2sls <- function(model, ...) influence(model)$rstudent

hatvalues.influence.2sls <- function(model, ...) model$hatvalues

cooks.distance.influence.2sls <- function(model, ...) model$cookd

cooks.distance.2sls <- function(model, ...) influence(model)$cookd

dfbeta.influence.2sls <- function(model, ...) model$dfbeta

dfbeta.2sls <- function(model, ...) influence(model)$dfbeta

hatvalues.2sls <- function(model, type=c("stage2", "both"), ...){
  type <- match.arg(type)
  hatvalues <- if (type == "stage2") NextMethod() else {
    hat.2 <- lm.influence(model)$hat
    model[c("qr", "rank", "residuals", "coefficients")] <- 
      list(model$qr.1, model$rank.1, model$residuals.1, model$coefficients.1)
    hat.1 <- lm.influence(model)$hat
    sqrt(hat.1*hat.2)
  }
  na.action <- model$na.action
  if(class(na.action) == "exclude") hatvalues[na.action]  <- NA
  hatvalues
}
