\name{ordPens-package}
\alias{ordPens-package}
\alias{ordPens}
\docType{package}
\title{
Selection and/or Smoothing of Ordinal Predictors
}
\description{
Selection and/or smoothing of ordinally scaled independent variables using 
a Group Lasso or generalized Ridge penalty.
}
\details{
\tabular{ll}{
Package: \tab ordPens\cr
Type: \tab Package\cr
Version: \tab 0.1-6\cr
Date: \tab 2011-05-27\cr
Depends: \tab grplasso\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

Smoothing and selection of ordinal predictors is done by the function 
\code{\link{ordSelect}}; smoothing only, by \code{\link{ordSmooth}}.
}
\author{
Jan Gertheiss

Maintainer: Jan Gertheiss <jan.gertheiss@stat.uni-muenchen.de>
}
\references{
Gertheiss, J. and G. Tutz (2009). \emph{Penalized regression with ordinal 
predictors}. International Statistical Review, 77, 345-365.

Gertheiss, J., S. Hogger, C. Oberhauser and G. Tutz (2011). \emph{Selection
of ordinally scaled independent variables with applications to international
classification of functioning core sets}.
Journal of the Royal Statistical Society C (Applied Statistics), 60, 377-395.

Meier, L., S. van de Geer and P. Buehlmann (2008). \emph{The group lasso for 
logistic regression}. Journal of the Royal Statistical Society B, 70, 53-71.

Yuan, M. and Y. Lin (2006). \emph{Model selection and estimation in regression 
with grouped variables}. Journal of the Royal Statistical Society B, 68, 49-67.
}

\keyword{ package }

\seealso{\code{\link{ordSelect}}, \code{\link{ordSmooth}}}


\examples{
# smooth modeling of a random dataset
set.seed(123)

# generate (ordinal) predictors
x1 <- sample(1:8,100,replace=TRUE)
x2 <- sample(1:6,100,replace=TRUE)
x3 <- sample(1:7,100,replace=TRUE)

# the response
y <- -1 + log(x1) + sin(3*(x2-1)/pi) + rnorm(100)

# x matrix
x <- cbind(x1,x2,x3)

# lambda values
lambda <- c(1000,500,200,100,50,30,20,10,1)

# smooth modeling
o1 <- ordSmooth(x = x, y = y, lambda = lambda)

# results
round(o1$coef,digits=3)
plot(o1)

# If for a certain plot the x-axis should be annotated in a different way,
# this can (for example) be done as follows:
plot(o1, whichx = 1, xlim = c(0,9), xaxt = "n")
axis(side = 1, at = c(1,8), labels = c("no agreement","total agreement"))
}
