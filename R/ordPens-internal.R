coding <- function(x, constant=TRUE, splitcod=TRUE)
  {
  # Converts ordinal/categorical data into split-coding, or dummy-coding.

  # Input:
  # x: data matrix with ordinal/categorical data
  # constant: should a constant added in the regression

  # Output:
  # Split/dummy-coded data-matrix
  # Reference-Category = 1
  n <- nrow(x)
  kx <- apply(x,2,max) - 1
  xds <- matrix(0,n,sum(kx))

  # Loop to run through the columns of x
  for (j in 1:ncol(x))
		{
      j1col <- ifelse(j>1,sum(kx[1:(j-1)]),0)
  		# Loop to run through the rows of x
  		for (i in 1:n)
  			{
          if (x[i,j] > 1)
            {
              if (splitcod)
                xds[i,j1col + 1:(x[i,j]-1)] <- 1
              else
                xds[i,j1col + (x[i,j]-1)] <- 1
            }
  		  }
    }
  ## Output
  if (constant)
    return(cbind(1,xds))
  else
    return(xds)
 	}



genRidge <- function(x, y, offset, omega, lambda, model, delta=1e-6, maxit=25)
  {
    coefs <- matrix(0,ncol(x),length(lambda))
    fits <- matrix(NA,nrow(x),length(lambda))
    l <- 1
    for (lam in lambda)
      {
        if (model == "linear")
          {
            yw <- y - offset
            chollam <- chol(t(x)%*%x + lam*omega)
            coefs[,l] <- backsolve(chollam,
            backsolve(chollam, t(x)%*%yw, transpose=TRUE))
            fits[,l] <- x%*%coefs[,l] + offset
          }
        else
          {
            ## penalized logistic/poisson regression
            conv <- FALSE
    
            # start values
            if (l > 1)
              {
                b.start <- coefs[,l-1]
              }
            else
              {
                if (model == "logit")
                  {
                    gy <- rep(log(mean(y)/(1-mean(y))),length(y)) - offset
                  }
                else
                  {
                    gy <- rep(log(mean(y)),length(y)) - offset
                  }
                chollam <- chol(t(x)%*%x + lam*omega)
                b.start <- backsolve(chollam,
                backsolve(chollam, t(x)%*%gy, transpose=TRUE))
              }
    
            # fisher scoring
            b.old <- b.start
            i <- 1
            while(!conv)
              {
                eta <- x%*%b.old + offset
                if (model == "logit")
                  {
                    mu <- plogis(eta)
                    sigma <- as.numeric(mu*(1-mu))
                  }
                else
                  {
                    mu <- exp(eta)
                    sigma <- as.numeric(mu)
                  }
                score <- t(x)%*%(y-mu)
                fisher <- t(x)%*%diag(sigma)%*%x
                choli <- chol(fisher + lam*omega)
                b.new <- b.old + backsolve(choli,
                backsolve(choli, (score - lam*omega%*%b.old), transpose=TRUE))
  
                if(sqrt(sum((b.new-b.old)^2)/sum(b.old^2)) < delta | i>=maxit)
                  {
                    # check the stop criterion
                    conv <- TRUE
                  }
                b.old <- b.new
                i <- i+1
              }  # end while
            coefs[,l] <- b.old   
            fits[,l] <- mu       
          }
        l <- l+1
      }
    rownames(fits) <- NULL
    colnames(fits) <- lambda
    rownames(coefs) <- NULL
    colnames(coefs) <- lambda
    
    return(list(fitted = fits, coefficients = coefs))   
  }