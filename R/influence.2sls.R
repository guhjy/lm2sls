
#' Deletion Diagnostic Methods for \code{"2sls"} Objects
#'
#' @aliases 2SLS_Diagnostics
#' @param model A \code{"2sls"} or \code{"influence.2sls"} object.
#' @param sigma. If \code{TRUE} (the default for 1000 or fewer cases), the deleted value
#' of the residual standard deviation is computed for each case; if \code{FALSE}, the
#' overall residual standard deviation is used to compute other deletion diagnostics.
#' @param type If \code{"stage2"} (the default), hatvalues are for the second stage regression;
#' if \code{"both"}, the hatvalues are the geometric mean of the casewise hatvalues for the
#' two stages; if \code{"maximum"}, the hatvalues are the larger of the casewise
#' hatvalues for the two stages. In computing the geometric mean or casewise maximum hatvalues,
#' the hatvalues for each stage are first divided by their average (number of coefficients in
#' stage regression/number of cases); the geometric mean or casewise maximum values are then
#' multiplied by the average hatvalue from the second stage.
#' @param ... arguments to be passed down.
#'
#' @return In the case of \code{influence.2sls}, an object of class \code{"influence.2sls"}
#' with the following components:
#' \describe{
#' \item{\code{coefficients}}{the estimated regression coefficients}
#' \item{\code{model}}{the model matrix}
#' \item{\code{dfbeta}}{influence on coefficients}
#' \item{\code{sigma}}{deleted values of the residual standard deviation}
#' \item{\code{dffits}}{overall influence on the regression coefficients}
#' \item{\code{cookd}}{Cook's distances}
#' \item{\code{hatvalues}}{hatvalues}
#' \item{\code{rstudent}}{Studentized residuals}
#' \item{\code{df.residual}}{residual degrees of freedom}
#' }
#' In the case of other methods, such as \code{rstudent.2sls} or
#' \code{rstudent.influence.2sls}, the corresponding diagnostic statistics.
#'
#' @description Methods for computing deletion diagnostics for 2SLS regression.
#' It's generally more efficient to compute the diagnostics via the \code{influence}
#' method and then to extract the various specific diagnostics with the methods for
#' \code{"influence.2sls"} objects. Other diagnostics for linear models, such as
#' added-variable plots (\code{\link[car]{avPlots}}) and component-plus-residual
#' plots (\code{\link[car]{crPlots}}), also work, as do effect plots
#' (e.g., \code{\link[effects]{predictorEffects}}) with residuals (see the examples below).
#' The pointwise confidence envelope for the \code{qqPlot} methods assumes an independent random sample
#' from the t distribution with degrees of freedom equal to the residual degrees of
#' freedom for the model and so are approximate, because the studentized residuals aren't
#' independent.
#'
#' @importFrom stats influence
#' @export
#' @seealso \code{\link{lm2sls}}, \link{2SLS_Methods}, \code{\link[car]{avPlots}},
#'   \code{\link[car]{crPlots}}, \code{\link[effects]{predictorEffects}},
#'   \code{\link[car]{qqPlot}}, \code{\link[car]{influencePlot}},
#'   \code{\link[car]{infIndexPlot}}
#' @examples
#' kmenta.eq1 <- lm2sls(Q ~ P + D, ~ D + F + A, data=Kmenta)
#' car::avPlots(kmenta.eq1)
#' car::crPlots(kmenta.eq1)
#' car::influencePlot(kmenta.eq1)
#' car::influenceIndexPlot(kmenta.eq1)
#' car::qqPlot(kmenta.eq1)
#' if (require(effects)){
#'   plot(effects::predictorEffects(kmenta.eq1, residuals=TRUE))
#' }
influence.2sls <- function(model, sigma. = n <= 1e3, type=c("stage2", "both"), ...){

  type <- match.arg(type)

  Z <- model$model.matrix.instruments
  X <- model$model.matrix
  X.fit <- model$fitted.1
  y <- model$y
  b <- model$coefficients
  res <- na.remove(model$residuals)
  sigma2 <- model$sigma^2
  hatvalues <-  na.remove(hatvalues(model, type=type))

  na.action <- model$na.action

  rnames <- rownames(X)
  cnames <- colnames(X)

  names(hatvalues) <- rnames

  w <- na.remove(model$weights)
  if (!is.null(w)){
    w <- sqrt(w)
    X <- diagprod(w, X)
    Z <- diagprod(w, Z)
    X.fit <- diagprod(w, X.fit)
    y <- w*y
  }
  else w <- 1

  rss <- sum((w*res)^2)
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


  rstudent <- w*res/(sigma * sqrt(1 - hatvalues))

  cookd <- (sigma^2/sigma2)*dffits^2/p

  result <- list(model = model.matrix(model),
                 coefficients=coef(model),
                 dfbeta = naresid(na.action, dfbeta),
                 sigma = naresid(na.action, sigma),
                 dffits = naresid(na.action, dffits),
                 cookd = naresid(na.action, cookd),
                 hatvalues = naresid(na.action, hatvalues),
                 rstudent = naresid(na.action, rstudent),
                 df.residual = df.residual(model))
  class(result) <- "influence.2sls"
  result
}

#' @rdname influence.2sls
#' @importFrom stats rstudent
#' @export
rstudent.2sls <- function(model, ...) {
  influence(model)$rstudent
}

#' @rdname influence.2sls
#' @importFrom stats cooks.distance
#' @method cooks.distance 2sls
#' @export
cooks.distance.2sls <- function(model, ...) {
  influence(model)$cookd
}

#' @rdname influence.2sls
#' @export
dfbeta.influence.2sls <- function(model, ...) {
  model$dfbeta
}

#' @rdname influence.2sls
#' @importFrom stats dfbeta
#' @export
dfbeta.2sls <- function(model, ...) {
  influence(model)$dfbeta
}

#' @rdname influence.2sls
#' @importFrom stats hatvalues lm.influence
#' @export
hatvalues.2sls <- function(model, type=c("stage2", "both", "maximum"), ...){
  type <- match.arg(type)
  hatvalues <- if (type == "stage2") NextMethod() else {
    n <- model$n
    p <- model$p
    q <- model$q
    mean1 <- q/n
    mean2 <- p/n
    hat.2 <- lm.influence(model)$hat/mean2
    model[c("qr", "rank", "residuals", "coefficients")] <-
      list(model$qr.1, model$rank.1, model$residuals.1, model$coefficients.1)
    hat.1 <- lm.influence(model)$hat/mean1
    hat <- if (type == "both") {
      sqrt(hat.1*hat.2)
    } else {
      pmax(hat.1, hat.2)
    }
    mean2*hat
  }
  na.action <- model$na.action
  if(class(na.action) == "exclude") hatvalues[na.action]  <- NA
  hatvalues
}

#' @rdname influence.2sls
#' @method rstudent influence.2sls
#' @export
rstudent.influence.2sls <- function(model, ...) {
  model$rstudent
}

#' @rdname influence.2sls
#' @method hatvalues influence.2sls
#' @export
hatvalues.influence.2sls <- function(model, ...) {
  model$hatvalues
}

#' @rdname influence.2sls
#' @export
#' @method cooks.distance influence.2sls
cooks.distance.influence.2sls <- {
  function(model, ...) model$cookd
}

#' @rdname influence.2sls
#' @importFrom car qqPlot
#' @importFrom graphics par
#' @export
qqPlot.2sls <- function(x,
                        ylab=paste("Studentized Residuals(",deparse(substitute(x)), ")", sep=""),
                        distribution=c("t", "norm"), ...){
  distribution <- match.arg(distribution)
  rstudent <- rstudent(x)
  if (distribution == "t"){
    car::qqPlot(rstudent, ylab=ylab, distribution="t", df=df.residual(x), ...)
  } else {
    car::qqPlot(rstudent, ylab=ylab, distribution="norm", ...)
  }
}

#' @rdname influence.2sls
#' @method qqPlot influence.2sls
#' @param x A \code{"2sls"} or \code{"influence.2sls"} object.
#' @param distribution \code{"t"} (the default) or \code{"norm"}.
#' @param ylab The vertical axis label.
#' @export
qqPlot.influence.2sls <- function(x,
                                  ylab=paste("Studentized Residuals(",deparse(substitute(x)), ")", sep=""),
                                  distribution=c("t", "norm"), ...){
  distribution <- match.arg(distribution)
  rstudent <- rstudent(x)
  if (distribution == "t"){
    car::qqPlot(rstudent, ylab=ylab, distribution="t", df=df.residual(x), ...)
  } else {
    car::qqPlot(rstudent, ylab=ylab, ...)
  }
}

#' @rdname influence.2sls
#' @method influencePlot 2sls
#' @importFrom car influencePlot
#' @export
influencePlot.2sls <- function(model, ...){
  influencePlot(influence(model), ...)
}

#' @rdname influence.2sls
#' @method influencePlot influence.2sls
#' @export
influencePlot.influence.2sls <- function(model, ...){
  if (length(class(model)) == 1) {
    class(model) <- c(class(model), "lm")
    influencePlot(model)
  }
  else NextMethod()
}

#' @rdname influence.2sls
#' @method infIndexPlot 2sls
#' @export
infIndexPlot.2sls <- function(model, ...){
  infIndexPlot(influence(model), ...)
}

#' @rdname influence.2sls
#' @method infIndexPlot influence.2sls
#' @importFrom car infIndexPlot
#' @export
infIndexPlot.influence.2sls <- function(model, ...){
  if (length(class(model)) == 1) {
    class(model) <- c(class(model), "lm")
    infIndexPlot(model, ...)
  }
  else NextMethod()
}

#' @rdname influence.2sls
#' @method model.matrix influence.2sls
#' @param object An \code{"influence.2sls"} object.
#' @export
model.matrix.influence.2sls <- function(object, ...){
  object$model
}
