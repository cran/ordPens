\name{ordPens-package}
\alias{ordPens-package}
\alias{ordPens}
\docType{package}
\title{
Selection and/or Smoothing of Ordinal Predictors
}
\description{
Selection and/or smoothing of ordinally scaled independent variables using 
a group lasso or generalized ridge penalty.
}
\details{
\tabular{ll}{
Package: \tab ordPens\cr
Type: \tab Package\cr
Version: \tab 0.2-1\cr
Date: \tab 2013-06-12\cr
Depends: \tab grplasso, mgcv, RLRsim\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

Smoothing and selection of ordinal predictors is done by the function 
\code{\link{ordSelect}}; smoothing only, by \code{\link{ordSmooth}}. For
ANOVA with ordinal factors, use \code{\link{ordAOV}}.
}
\author{
Jan Gertheiss

Maintainer: Jan Gertheiss, \email{jgerthe@uni-goettingen.de}
}
\references{
Gertheiss, J., S. Hogger, C. Oberhauser and G. Tutz (2011). \emph{Selection
of ordinally scaled independent variables with applications to international
classification of functioning core sets}.
Journal of the Royal Statistical Society C (Applied Statistics), 60, 377-395.

Gertheiss, J. and F. Oehrlein (2011). \emph{Testing Relevance and Linearity of
Ordinal Predictors}, Electronic Journal of Statistics, 5, 1935-1959.

Gertheiss, J. and G. Tutz (2009). \emph{Penalized regression with ordinal 
predictors}. International Statistical Review, 77, 345-365.
}

\keyword{ package }

\seealso{\code{\link{ordSelect}}, \code{\link{ordSmooth}},
\code{\link{ordAOV}}}


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
