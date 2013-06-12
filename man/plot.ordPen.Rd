\name{plot.ordPen}
\alias{plot.ordPen}

\title{Plot method for ordPen objects}
\description{
  Takes a fitted \code{ordPen} object and plots estimated dummy coefficients
  of ordinal predictors for different \code{lambda} values.
}
\usage{
\method{plot}{ordPen}(x, whichlam = NULL, whichx = NULL, 
  type = NULL, xlab = NULL, ylab = NULL, main = NULL, 
  xlim = NULL, ylim = NULL, col = NULL, ...)
}

\arguments{
  \item{x}{an \code{ordPen} object.}
  \item{whichlam}{a vector of indices of \code{lambda} values corresponding
    to \code{object$lambda} for which plotting is done; if \code{NULL}, all
    values from \code{object$lambda} are considered.}
  \item{whichx}{a vector of indices indicating the ordinal predictors whose
    dummy coefficients are plotted; e.g., set \code{whichx=2}, if you just want 
    the plot for the second smooth term.}
  \item{type}{1-character string giving the type of plot desired, see
  \code{\link{plot.default}}.}
  \item{xlab}{a label for the x axis; if supplied then this will be used as the
    x label for all plots.}
  \item{ylab}{a label for the y axis; if supplied then this will be used as the
    y label for all plots.}
  \item{main}{a main title for the plot(s); if supplied then this will be used
    as the title for all plots.}
  \item{xlim}{the x limits; if supplied then this pair of numbers are used
    as the x limits for each plot.}
  \item{ylim}{the y limits; if supplied then this pair of numbers are used
    as the y limits for each plot.}
  \item{col}{the plotting color; can be a vector of the same length as
    \code{whichlam} specifying different colors for different \code{lambda}
    values. Default is shades of gray: the higher \code{lambda}
    the darker.}
  \item{...}{additional graphical parameters (see \code{\link{plot.default}},
    or  \code{\link{par}}).}
}
\value{
  The function simply generates plots.
}

\author{Jan Gertheiss}

\seealso{\code{\link{ordSelect}}, \code{\link{ordSmooth}}}

\examples{
# see for example
help(ordSelect)
}

\keyword{methods}
