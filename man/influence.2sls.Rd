% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/influence.2sls.R
\name{influence.2sls}
\alias{influence.2sls}
\alias{2SLS_Diagnostics}
\alias{rstudent.2sls}
\alias{cooks.distance.2sls}
\alias{dfbeta.influence.2sls}
\alias{dfbeta.2sls}
\alias{hatvalues.2sls}
\alias{rstudent.influence.2sls}
\alias{hatvalues.influence.2sls}
\alias{cooks.distance.influence.2sls}
\alias{qqPlot.2sls}
\alias{qqPlot.influence.2sls}
\alias{influencePlot.2sls}
\alias{influencePlot.influence.2sls}
\alias{infIndexPlot.2sls}
\alias{infIndexPlot.influence.2sls}
\alias{model.matrix.influence.2sls}
\title{Deletion Diagnostic Methods for \code{"2sls"} Objects}
\usage{
\method{influence}{2sls}(model, sigma. = n <= 1000, type = c("stage2",
  "both"), ...)

\method{rstudent}{2sls}(model, ...)

\method{cooks.distance}{2sls}(model, ...)

\method{dfbeta}{influence.2sls}(model, ...)

\method{dfbeta}{2sls}(model, ...)

\method{hatvalues}{2sls}(model, type = c("stage2", "both", "maximum"),
  ...)

\method{rstudent}{influence.2sls}(model, ...)

\method{hatvalues}{influence.2sls}(model, ...)

\method{cooks.distance}{influence.2sls}(model, ...)

\method{qqPlot}{2sls}(x, ylab = paste("Studentized Residuals(",
  deparse(substitute(x)), ")", sep = ""), distribution = c("t", "norm"),
  ...)

\method{qqPlot}{influence.2sls}(x, ylab = paste("Studentized Residuals(",
  deparse(substitute(x)), ")", sep = ""), distribution = c("t", "norm"),
  ...)

\method{influencePlot}{2sls}(model, ...)

\method{influencePlot}{influence.2sls}(model, ...)

\method{infIndexPlot}{2sls}(model, ...)

\method{infIndexPlot}{influence.2sls}(model, ...)

\method{model.matrix}{influence.2sls}(object, ...)
}
\arguments{
\item{model}{A \code{"2sls"} or \code{"influence.2sls"} object.}

\item{sigma.}{If \code{TRUE} (the default for 1000 or fewer cases), the deleted value
of the residual standard deviation is computed for each case; if \code{FALSE}, the
overall residual standard deviation is used to compute other deletion diagnostics.}

\item{type}{If \code{"stage2"} (the default), hatvalues are for the second stage regression;
if \code{"both"}, the hatvalues are the geometric mean of the casewise hatvalues for the
two stages; if \code{"maximum"}, the hatvalues are the larger of the casewise
hatvalues for the two stages. In computing the geometric mean or casewise maximum hatvalues,
the hatvalues for each stage are first divided by their average (number of coefficients in
stage regression/number of cases); the geometric mean or casewise maximum values are then
multiplied by the average hatvalue from the second stage.}

\item{...}{arguments to be passed down.}

\item{x}{A \code{"2sls"} or \code{"influence.2sls"} object.}

\item{ylab}{The vertical axis label.}

\item{distribution}{\code{"t"} (the default) or \code{"norm"}.}

\item{object}{An \code{"influence.2sls"} object.}
}
\value{
In the case of \code{influence.2sls}, an object of class \code{"influence.2sls"}
with the following components:
\describe{
\item{\code{coefficients}}{the estimated regression coefficients}
\item{\code{model}}{the model matrix}
\item{\code{dfbeta}}{influence on coefficients}
\item{\code{sigma}}{deleted values of the residual standard deviation}
\item{\code{dffits}}{overall influence on the regression coefficients}
\item{\code{cookd}}{Cook's distances}
\item{\code{hatvalues}}{hatvalues}
\item{\code{rstudent}}{Studentized residuals}
\item{\code{df.residual}}{residual degrees of freedom}
}
In the case of other methods, such as \code{rstudent.2sls} or
\code{rstudent.influence.2sls}, the corresponding diagnostic statistics.
}
\description{
Methods for computing deletion diagnostics for 2SLS regression.
It's generally more efficient to compute the diagnostics via the \code{influence}
method and then to extract the various specific diagnostics with the methods for
\code{"influence.2sls"} objects. Other diagnostics for linear models, such as
added-variable plots (\code{\link[car]{avPlots}}) and component-plus-residual
plots (\code{\link[car]{crPlots}}), also work, as do effect plots
(e.g., \code{\link[effects]{predictorEffects}}) with residuals (see the examples below).
The pointwise confidence envelope for the \code{qqPlot} methods assumes an independent random sample
from the t distribution with degrees of freedom equal to the residual degrees of
freedom for the model and so are approximate, because the studentized residuals aren't
independent.
}
\examples{
kmenta.eq1 <- lm2sls(Q ~ P + D, ~ D + F + A, data=Kmenta)
car::avPlots(kmenta.eq1)
car::crPlots(kmenta.eq1)
car::influencePlot(kmenta.eq1)
car::influenceIndexPlot(kmenta.eq1)
car::qqPlot(kmenta.eq1)
if (require(effects)){
  plot(effects::predictorEffects(kmenta.eq1, residuals=TRUE))
}
}
\seealso{
\code{\link{lm2sls}}, \link{2SLS_Methods}, \code{\link[car]{avPlots}},
  \code{\link[car]{crPlots}}, \code{\link[effects]{predictorEffects}},
  \code{\link[car]{qqPlot}}, \code{\link[car]{influencePlot}},
  \code{\link[car]{infIndexPlot}}
}
