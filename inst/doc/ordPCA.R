## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",  
  out.extra = 'style="border:0; margin: auto"', 
fig.width=5,
fig.height=5,
out.width="400px",  
out.height="400px"  
)

## ----setup, results="hide", warning=FALSE,message=FALSE-----------------------
library(ordPens)

## -----------------------------------------------------------------------------
library(psy) 
data(ehd)
H <- ehd + 1
head(H)

## -----------------------------------------------------------------------------
ehd_pca1 <- ordPCA(H, p = 5, lambda = 0.5, maxit = 100, crit = 1e-7,
                   qstart = NULL, Ks = apply(H, 2, max), constr = rep(TRUE, ncol(H)),
                   CV = FALSE, k = 5, CVfit = FALSE)

## -----------------------------------------------------------------------------
summary(ehd_pca1$pca)

## -----------------------------------------------------------------------------
ehd_pca1$qs

## ---- fig.align="center"------------------------------------------------------
plot(1:5, ehd_pca1$qs[[2]], type = "b", xlab = "category", ylab = "quantification", 
     col = 2, main = colnames(H)[2], bty = "n")

## -----------------------------------------------------------------------------
ehd_pca2 <- ordPCA(H, p = 5, lambda = c(5, 0.5, 0), maxit = 100, crit = 1e-7,
                   qstart = NULL, Ks = apply(H, 2, max), constr = rep(TRUE, ncol(H)),
                   CV = FALSE)

## -----------------------------------------------------------------------------
summary(ehd_pca2$pca[[1]])

## -----------------------------------------------------------------------------
ehd_pca2$qs[[2]]

## ---- fig.align="center"------------------------------------------------------
plot(ehd_pca2$qs[[2]][,3], type = "b", xlab = "category", ylab = "quantification", col = 1, 
     ylim = range(ehd_pca2$qs[[2]]), main = colnames(H)[2], bty = "n")
lines(ehd_pca2$qs[[2]][,2], type = "b", col = 2, lty = 2, pch = 2, lwd=2)
lines(ehd_pca2$qs[[2]][,1], type = "b", col = 3, lty = 3, pch = 3, lwd=2)

## ----test_a, figures-side, fig.show="hold", fig.width = 3, fig.height = 3, out.width="31%",out.height="30%"----
par(mar = c(4.1, 4.1, 3.1, 1.1))

for(j in c(1, 9, 12, 13, 15, 19)){ 
  plot(ehd_pca2$qs[[j]][,3], type = "b", main = colnames(H)[j], xlab = "category", 
       ylab = "quantification", lwd = 2, bty = "n") 
  lines(ehd_pca2$qs[[j]][,2], type = "b", col = 2, lty = 2, pch = 2, lwd = 2)
  lines(ehd_pca2$qs[[j]][,1], type = "b", col = 3, lty = 3, pch = 3, lwd = 2)
} 

## -----------------------------------------------------------------------------
ehd_pca3 <- ordPCA(H, p = 5, lambda = c(5, 0.5, 0.001), maxit = 100, crit = 1e-7,
                   qstart = NULL, Ks = apply(H, 2, max), constr = rep(TRUE, ncol(H)),
                   CV = TRUE, k = 5)

ehd_pca3$VAFtest

## -----------------------------------------------------------------------------
lambda <- 10^seq(4, -4, by = -0.1)
set.seed(456)
ehd_pca4 <- ordPCA(H, p = 5, lambda = lambda, maxit = 100, crit = 1e-7,
                   qstart = NULL, Ks = apply(H, 2, max), constr = rep(TRUE, ncol(H)),
                   CV = TRUE, k = 5, CVfit = FALSE)

ehd_pca4$VAFtest

## ----test, figures-side, fig.show="hold", fig.width = 5, fig.height = 5, out.width="45%",out.height="50%"----
plot(log10(lambda), apply(ehd_pca4$VAFtrain,2,mean), type = "l",
     xlab = expression(log[10](lambda)), ylab = "proportion of variance explained",
     main = "training data", cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.4)
plot(log10(lambda), apply(ehd_pca4$VAFtest,2,mean), type = "l",
     xlab = expression(log[10](lambda)), ylab = "proportion of variance explained",
     main = "validation data", cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.4)
abline(v = log10(lambda)[which.max(apply(ehd_pca4$VAFtest,2,mean))])

## -----------------------------------------------------------------------------
data(ICFCoreSetCWP) 
H <- ICFCoreSetCWP[, 1:67] + matrix(c(rep(1,50), rep(5,16), 1), 
                                    nrow(ICFCoreSetCWP), 67, byrow = TRUE)
head(H)
xnames <- colnames(H)

## ----test_environment, figures-side, fig.show="hold", fig.width = 5, fig.height = 5, out.width="45%",out.height="50%"----
icf_pca1 <- ordPCA(H, p = 2, lambda = c(5, 0.5, 0.001), maxit = 100, crit = 1e-7, qstart = NULL, 
                   Ks = c(rep(5,50), rep(9,16), 5), 
                   constr = c(rep(TRUE,50), rep(FALSE,16), TRUE), 
                   CV = FALSE, k = 5, CVfit = FALSE) 
 
icf_pca1C <- ordPCA(H, p = 2, lambda = c(5, 0.5, 0.001), maxit = 100, crit = 1e-7, qstart = NULL, 
                    Ks = c(rep(5,50), rep(9,16), 5), constr = rep(TRUE, ncol(H)), 
                    CV = FALSE, k = 5, CVfit = FALSE) 

plot(icf_pca1$qs[[51]][,3], type = "b", xlab = "category", ylab = "quantification", col = 1, 
     ylim = range(icf_pca1$qs[[51]]), main = xnames[51], bty = "n", xaxt = "n") 
lines(icf_pca1$qs[[51]][,2], type = "b", col = 2, lty = 2, pch = 2, lwd=2)
lines(icf_pca1$qs[[51]][,1], type = "b", col = 3, lty = 3, pch = 3, lwd=2)
axis(1, at = 1:length(icf_pca1$qs[[51]][,1]), labels = -4:4)   

plot(icf_pca1C$qs[[51]][,3], type = "b", xlab = "category", ylab = "quantification", col = 1, 
     ylim = range(icf_pca1C$qs[[51]]), main = xnames[51], bty = "n", xaxt = "n")  
lines(icf_pca1C$qs[[51]][,2], type = "b", col = 2, lty = 2, pch = 2, lwd=2)
lines(icf_pca1C$qs[[51]][,1], type = "b", col = 3, lty = 3, pch = 3, lwd=2)
axis(1, at = 1:length(icf_pca1C$qs[[51]][,1]), labels = -4:4)   

## ----test_pca, figures-side, fig.show="hold", fig.width = 3, fig.height = 3, out.width="31%",out.height="30%"----
wx <- which(xnames=="b265"|xnames=="d450"|xnames=="d455"|xnames=="e1101"|xnames=="e460"
            |xnames=="e325") 
xmain <- c()
xmain[wx] <- list("Touch function",
                  "Walking",
                  "Moving around",
                  "Drugs",
                  "Societal attitudes",
                   paste(c("Acquaintances,colleagues,","peers,community members")))

par(mar = c(4.1, 4.1, 3.1, 1.1))
for (j in wx){
plot(icf_pca1$qs[[j]][,3], type = "b", main = xmain[j], xlab = "category", bty = "n", 
     ylim = range(icf_pca1$qs[[j]]), ylab = "quantification", xaxt = "n", cex.main= 1)  
lines(icf_pca1$qs[[j]][,2], type = "b", col = 2, lty = 2, pch = 2, lwd=2)
lines(icf_pca1$qs[[j]][,1], type = "b", col = 3, lty = 3, pch = 3, lwd=2)
axis(1, at = 1:length(icf_pca1$qs[[j]][,1]), 
     labels = ((1:length(icf_pca1$qs[[j]][,1])) - c(rep(1,50), rep(5,16), 1)[j]))   
}

