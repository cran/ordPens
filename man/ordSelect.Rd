\name{ordSelect}
\alias{ordSelect}
\alias{ordSelect.default}
\title{Selection and smoothing of dummy coefficients of ordinal predictors}
\description{Fits dummy coefficients of ordinally scaled independent variables
with a group lasso penalty on differences of adjacent dummy coefficients.}

\usage{
ordSelect(x, ...)

\method{ordSelect}{default}(x, y, u = NULL, z = NULL, 
  offset = rep(0,length(y)), lambda, nu = 1, zeta = 1,  
  model = c("linear", "logit", "poisson"), penscale = sqrt, 
  scalex = TRUE, scalez = TRUE, scaleu = TRUE,
  nonpenx = NULL, nonpenz = NULL, nonpenu = NULL, 
  intercept = TRUE, eps = 1e-3, ...)
}
\arguments{
  \item{x}{the matrix of ordinal predictors, with each column corresponding to
    one predictor and containing numeric values from \{1,2,...\}; for each 
    covariate, category 1 is taken as reference category with zero dummy 
    coefficient.}
  \item{y}{the response vector.}
  \item{u}{a matrix (or \code{data.frame}) of additional categorical (nominal) 
    predictors, with each column corresponding to one (additional) predictor and
    containing numeric values from \{1,2,...\}; corresponding dummy coefficients
    will be regularized using a groupwise simple Ridge penalty, and for each 
    covariate category 1 is taken as reference category.}
  \item{z}{a matrix (or \code{data.frame}) of additional metric predictors, with 
    each column corresponding to one (additional) predictor; corresponding
    coefficients will be regularized using a simple Lasso penalty.}
  \item{offset}{vector of offset values.}
  \item{lambda}{vector of penalty parameters (in decreasing order).
    Optimization starts with the first component. See details below.}
  \item{nu}{additional tuning parameter to control the strength of the penalty
    imposed on dummy coefficients corresponding to \code{u}. See details below.}
  \item{zeta}{additional tuning parameter  to control the strength of the penalty
    imposed on coefficients corresponding to \code{z}. See details below.}
  \item{model}{the model which is to be fitted. Possible choices are "linear"
    (default), "logit" or "poisson". See details below.}
  \item{penscale}{rescaling function to adjust the value of the penalty
    parameter to the degrees of freedom of the parameter group. See the
    references below.}
  \item{scalex, scaleu, scalez}{logical. Should (split/dummy-coded) design
    matrices corresponding to \code{x}, \code{u}, resp. \code{z} be scaled
    to have unit variance over columns? See details below.}
  \item{nonpenx, nonpenu, nonpenz}{vectors of indices indicating columns of
    \code{x}, \code{u}, resp. \code{z} whose regression coefficients are not
    penalized.}
  \item{intercept}{logical. Should a (non-penalized) intercept be included
    in the model? Default is TRUE.}
  \item{eps}{a (small) constant to be added to the columnwise standard
    deviations when scaling design matrices, to control the effect of very small
    stds. See details below.}
  \item{...}{additional arguments.}
}
\details{
  The method assumes that categorical covariates (contained in \code{x} and 
  \code{u}) take values 1,2,...,max, where max denotes the (columnwise) highest 
  level observed in the data. If any level between 1 and max is not observed, 
  a corresponding (dummy) coefficient is fitted anyway. If any level > max is 
  not observed but possible in principle, and a corresponding coefficient is to 
  be fitted, the easiest way is to add a corresponding row to \code{x} (and 
  \code{u},\code{z}) with corresponding \code{y} value being \code{NA}.

  
  In order to adjust the strength of penalty that is imposed on nominal or 
  metric covariates, \code{nu}, resp. \code{zeta} can be chosen. The penalty 
  parameter used for nominal variables (contained in \code{u}) is 
  \code{lambda*nu}, for metric predictors (contained in \code{z}) it is 
  \code{lambda*zeta}. 
  
  
  If a linear regression model is fitted, response vector \code{y} may contain 
  any numeric values; if a logit model is fitted, \code{y} has to be 0/1 coded;
  if a poisson model is fitted, \code{y} has to contain count data. 
  
  
  If \code{scalex}, \code{scaleu} or \code{scalez} are \code{TRUE}, design 
  matrices constructed from \code{x}, \code{u}, resp. \code{z} are scaled to have 
  unit variance over columns. In the case of \code{x}, the design matrix is 
  split-coded, in the case of \code{u}, it is dummy-coded, and in case of 
  \code{z}, it is just \code{z}. If a certain \code{x}- or \code{u}- category, 
  however, is observed only a few times, variances may become very small and
  scaling has enormous effects on the result and may cause numerical problems.
  Hence a small constant \code{eps} can be added to each standard deviation 
  when used for scaling. 
  
}

\value{An \code{ordPen} object, which is a list containing:
  \item{fitted}{the matrix of fitted response values of the training data.
    Columns correspond to different \code{lambda} values.}
  \item{coefficients}{the matrix of fitted coefficients with respect to 
    dummy-coded (ordinal or nominal) categorical input variables (including the 
    reference category) as well as metric predictors. Columns correspond to 
    different lambda values.}
  \item{model}{the type of the fitted model: "linear", "logit", or "poisson".}
  \item{lambda}{the used lambda values.}
  \item{xlevels}{a vector given the number of levels of the ordinal predictors.}
  \item{ulevels}{a vector given the number of levels of the nominal predictors.}
  \item{zcovars}{the number of metric covariates.}
}
\references{
  Meier, L., S. van de Geer and P. Buehlmann (2008). \emph{The
  group lasso for logistic regression}. Journal of the Royal
  Statistical Society B, 70, 53 - 71.
  
  Gertheiss, J., S. Hogger, C. Oberhauser and G. Tutz (2011). \emph{Selection
  of ordinally scaled independent variables with applications to international
  classification of functioning core sets}.
  Journal of the Royal Statistical Society C (Applied Statistics), to appear.
 }

\author{Jan Gertheiss, \email{jan.gertheiss@stat.uni-muenchen.de}}


\seealso{\code{\link{plot.ordPen}}, \code{\link{predict.ordPen}}, 
\code{\link{ICFCoreSetCWP}}}


\examples{
# smoothing and selection of ordinal covariates using a random dataset
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

# smoothing and selection
o2 <- ordSelect(x = x, y = y, lambda = lambda)

# results
round(o2$coef,digits=3)
plot(o2)

# If for a certain plot the x-axis should be annotated in a different way,
# this can (for example) be done as follows:
plot(o2, whichx = 1, xlim = c(0,9), xaxt = "n")
axis(side = 1, at = c(1,8), labels = c("no agreement","total agreement"))
}

\keyword{models}
\keyword{regression}
