# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# Solves the weighted conditional logit with Box-Cox transforms with U(Cost), U(Cost, Duration) and U(Cost, Duration, Length)
# using the mnLogit package
#
library(mlogit)
library(mnlogit)

# Solves the conditional logit with Box-Cox transforms with U(Cost), U(Cost, Duration) or U(Cost, Duration, Length)
solveBoxCoxLogit <- function(x, lambdas, nbCores = 1) {
  nbLambdas <- length(lambdas)
  
  # Replace null quantities with a small one
  smallQty <- .Machine$double.xmin
  x$qty.1[x$qty.1 == 0] <- smallQty
  x$qty.2[x$qty.2 == 0] <- smallQty
  x$qty.3[x$qty.3 == 0] <- smallQty
  
  # Cost is always used as independent variable
  
  # Replace missing costs with an high value
  highValue <- max(x$cost.1, x$cost.2, x$cost.3, na.rm = TRUE) * 1000
  x$cost.1[is.na(x$cost.1)] <- highValue
  x$cost.2[is.na(x$cost.2)] <- highValue
  x$cost.3[is.na(x$cost.3)] <- highValue
  
  # Box-Cox transform the variable
  lambdaCost <- lambdas[1]
  
  if (lambdaCost != 0) {
    x$cost.1 <- (x$cost.1^lambdaCost - 1) / lambdaCost
    x$cost.2 <- (x$cost.2^lambdaCost - 1) / lambdaCost
    x$cost.3 <- (x$cost.3^lambdaCost - 1) / lambdaCost
  }
  else {
    x$cost.1 <- log(x$cost.1)
    x$cost.2 <- log(x$cost.2)
    x$cost.3 <- log(x$cost.3)
  }
  
  # Formula for the conditional multinomial logit
  f <- mFormula(mode ~ cost | 1 | 1)
  
  
  # Bivariate and trivariate cases : duration is also used as independent variable
  if (nbLambdas > 1) {
    # Replace missing durations with an high value
    highValue <- max(x$duration.1, x$duration.2, x$duration.3, na.rm = TRUE) * 1000
    x$duration.1[is.na(x$duration.1)] <- highValue
    x$duration.2[is.na(x$duration.2)] <- highValue
    x$duration.3[is.na(x$duration.3)] <- highValue
    
    # Box-Cox transform the variable
    lambdaDuration <- lambdas[2]
    if (lambdaDuration != 0) {
      x$duration.1 <- (x$duration.1^lambdaDuration - 1) / lambdaDuration
      x$duration.2 <- (x$duration.2^lambdaDuration - 1) / lambdaDuration
      x$duration.3 <- (x$duration.3^lambdaDuration - 1) / lambdaDuration
    }
    else {
      x$duration.1 <- log(x$duration.1)
      x$duration.2 <- log(x$duration.2)
      x$duration.3 <- log(x$duration.3)
    }
    
    # Formula for the conditional multinomial logit
    f <- mFormula(mode ~ cost + duration | 1 | 1)
  }
  
  # Trivariate case : add Length to the explanatory variables
  if (nbLambdas > 2) {
    # Replace missing lengths with an high value
    highValue <- max(x$length.1, x$length.2, x$length.3, na.rm = TRUE) * 1000
    x$length.1[is.na(x$length.1)] <- highValue
    x$length.2[is.na(x$length.2)] <- highValue
    x$length.3[is.na(x$length.3)] <- highValue
    
    lambdaLength <- lambdas[3]
    if (lambdaLength != 0) {
      x$length.1 <- (x$length.1^lambdaLength - 1) / lambdaLength
      x$length.2 <- (x$length.2^lambdaLength - 1) / lambdaLength
      x$length.3 <- (x$length.3^lambdaLength - 1) / lambdaLength
    }
    else {
      x$length.1 <- log(x$length.1)
      x$length.2 <- log(x$length.2)
      x$length.3 <- log(x$length.3)
    }
    
    # Formula for the conditional multinomial logit
    f <- mFormula(mode ~ cost + duration + length | 1 | 1)
  }
  
  # Total transported quantity
  x$totQty <- x$qty.1 + x$qty.2 + x$qty.3
  
  # Create wideData data, with one record per mode for each OD pair (see mlogit documentation)
  wideData <- data.frame()
  for (mode in 1:3) {
    wd <- data.frame(mode = integer(nrow(x)))
    wd$mode <- mode
    
    wd$cost.1 <- x$cost.1
    wd$cost.2 <- x$cost.2
    wd$cost.3 <- x$cost.3
    
    if (nbLambdas > 1) {
      wd$duration.1 <- x$duration.1
      wd$duration.2 <- x$duration.2
      wd$duration.3 <- x$duration.3
    }
    
    if (nbLambdas > 2) {
      wd$length.1 <- x$length.1
      wd$length.2 <- x$length.2
      wd$length.3 <- x$length.3
    }
    
    wd$qty <- x$qty.1
    if (mode == 2) {
      wd$qty <- x$qty.2
    } else if (mode == 3) {
      wd$qty <- x$qty.3
    }
    
    wideData <- rbind(wideData, wd)
  }
  
  # Transform into "long" format data (see mlogit documentation)
  n <- as.integer(3 * nbLambdas + 1)
  longData <-
    mlogit.data(wideData,
                choice = "mode",
                shape = "wide",
                varying = 2:4
    ) # First column is "mode", variables are in columns 2 to 4.
  
  # mlogit version 1.1 returns a dfidx object, but mnlogit only accepts data frames for now
  # Be sure to have a dataframe
  longData = as.data.frame(longData)
  
  # Solve the model, using the mnlogit package, faster (parallelized) that mlogit
  model <- mnlogit(
    f,
    choiceVar = "alt",
    longData,
    weights = wideData$qty, # This is a weighted logit
    na.rm = FALSE,
    ncores = nbCores
  )
  
  return(model)
}
