# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# A few convenience functions used in more than one script in this project
#


# Test if all the estimators are of the expected sign (must be negative)
allSignsAreExpected <- function(model) {
  c <- coef(model)
  correctSign <- TRUE
  # Browse de coefficients names (see output of "summary(model)")
  for (j in 1:length(c)) {
    name <- names(c[j])
    if (substring(name, 1, 1) != "(") {
      # "(Intercept)" must not be tested
      if (c[name] > 0) {
        correctSign <- FALSE
        break
      }
    }
  }

  return(correctSign)
}

# Get Pr(>|t|) of model for all coefficients
allCoefsAreSignificant <- function(model, nbStars, withSignificantIntercepts) {
  b <- coef(model, order = TRUE)
  std.err2 <- sqrt(diag(vcov(model)))
  std.err <- b
  std.err[names(std.err2)] <- std.err2
  
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  
  # Test or not the signif.level of the 2 intercepts
  startIdx <- 3
  if (withSignificantIntercepts) {
    startIdx <- 1
  }
  startIdx = 1
  for (j in startIdx:length(p)) {
    if (nbStars == 3 && p[j] > 0.001) {
      return(FALSE)
    }
    if (nbStars == 2 && p[j] > 0.01) {
      return(FALSE)
    }
    if (nbStars == 1 && p[j] > 0.05) {
      return(FALSE)
    }
    if (nbStars == 0 && p[j] > 0.1) {
      return(FALSE)
    }
    if (nbStars == 0 && p[j] > 1) {
      return(FALSE)
    }
  }
  return(TRUE)
}

# Get Pr(>|t|) of model for a given coefficient
getStars <- function(model, coefName) {
  b <- coef(model, order = TRUE)
  std.err2 <- sqrt(diag(vcov(model)))
  std.err <- b
  std.err[names(std.err2)] <- std.err2
  
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  
  pp <- p[coefName]
  if (pp <= 0.001) {
    return("***")
  }
  if (pp <= 0.01) {
    return("**")
  }
  if (pp <= 0.05) {
    return("*")
  }
  if (pp <= 0.1) {
    return(".")
  }
  if (pp <= 1) {
    return(" ")
  }
}

# Returns TRUE if signs are expected and  all the estimators are enough significant
isValid <- function(solution) {
  
  # When looking in stored results (HeuristicVsBruteForce.R)
  if("keep" %in% colnames(solution)) {
    return (solution$keep) 
  }
  
  # When used in "real" conditions (Heuristic.R)
  if (solution$error == "") {
    return(TRUE)
  }
  return(FALSE)
}

# Returns TRUE if all the signs are expected
hasExpectedSigns <- function(solution) {
  if (solution$error == "") {
    return(TRUE)
  }
  return(FALSE)
}


# Returns the lambdas of a given solution
getLambdas <- function(solution, nbVariables) {
  lambdas <- c()
  for (i in 1:nbVariables) {
    if (i == 1) s <- "solution$lambda.cost"
    if (i == 2) s <- "solution$lambda.duration"
    if (i == 3) s <- "solution$lambda.length"
    value <<- eval(parse(text = s))
    lambdas <- c(lambdas, value)
  }
  return(lambdas)
}

# Draw a random combination of lambda's
randomDrawLambdas <- function(nbLambdas, range, granularity) {
  lambdas <- c()
  for (j in 1:nbLambdas) {
    z <- sample(1:nbSteps, 1)
    lambda <- -range - granularity + (z * granularity)
    lambdas <- c(lambdas, round(lambda, 1))
  }
  return(lambdas)
}

# Identify, in the brute force results, the solutions with a given signif level.
# One can decide to test or not the signif level of the intercepts.
# (This is a piece of ugly code that could be improved)
#
# bfSolution : dataframe that contains the brute force solutions
# nbVariables : 1, 2 or 3
# minSigLevel : minimum signif. level to retain for the estimators (# " " = 1, ." = 0, "*" = 1, "**" = 2, "***" = 3)
# withSignificantIntercepts : if TRUE, the signif. levels of the intercepts must also be larger or equal that minSigLevel
markValidSolutions <- function(bfSolutions, nbVariables, minSigLevel, withSignificantIntercepts) {
  bfSolutions$sig1b <- -1
  bfSolutions$sig1b[bfSolutions$sig.const.iww == "."] <- 0
  bfSolutions$sig1b[bfSolutions$sig.const.iww == "*"] <- 1
  bfSolutions$sig1b[bfSolutions$sig.const.iww == "**"] <- 2
  bfSolutions$sig1b[bfSolutions$sig.const.iww == "***"] <- 3
  
  bfSolutions$sig2b <- -1
  bfSolutions$sig2b[bfSolutions$sig.const.rail == "."] <- 0
  bfSolutions$sig2b[bfSolutions$sig.const.rail == "*"] <- 1
  bfSolutions$sig2b[bfSolutions$sig.const.rail == "**"] <- 2
  bfSolutions$sig2b[bfSolutions$sig.const.rail == "***"] <- 3
  
  
  bfSolutions$sig3b <- -1
  bfSolutions$sig3b[bfSolutions$sig.cost == "."] <- 0
  bfSolutions$sig3b[bfSolutions$sig.cost == "*"] <- 1
  bfSolutions$sig3b[bfSolutions$sig.cost == "**"] <- 2
  bfSolutions$sig3b[bfSolutions$sig.cost == "***"] <- 3
  
  if (nbVariables > 1) {
    bfSolutions$sig4b <- -1
    bfSolutions$sig4b[bfSolutions$sig.duration == "."] <- 0
    bfSolutions$sig4b[bfSolutions$sig.duration == "*"] <- 1
    bfSolutions$sig4b[bfSolutions$sig.duration == "**"] <- 2
    bfSolutions$sig4b[bfSolutions$sig.duration == "***"] <- 3
  }
  
  if (nbVariables > 2) {
    bfSolutions$sig5b <- -1
    bfSolutions$sig5b[bfSolutions$sig.length == "."] <- 0
    bfSolutions$sig5b[bfSolutions$sig.length == "*"] <- 1
    bfSolutions$sig5b[bfSolutions$sig.length == "**"] <- 2
    bfSolutions$sig5b[bfSolutions$sig.length == "***"] <- 3
  }
  
  bfSolutions$keep <- FALSE
  bfSolutions$best <- FALSE
  
  if (nbVariables == 1) {
    if (withSignificantIntercepts) {
      bfSolutions$keep[bfSolutions$cost < 0 & bfSolutions$sig1b >= minSignifLevel & bfSolutions$sig2b >= minSignifLevel & bfSolutions$sig3b >= minSignifLevel] <- TRUE
    } else {
      bfSolutions$keep[bfSolutions$cost < 0 & bfSolutions$sig3b >= minSignifLevel] <- TRUE
    }
  }
  
  
  if (nbVariables == 2) {
    if (withSignificantIntercepts) {
      bfSolutions$keep[bfSolutions$cost < 0 & bfSolutions$duration < 0 & bfSolutions$sig1b >= minSignifLevel & bfSolutions$sig2b >= minSignifLevel & bfSolutions$sig3b >= minSignifLevel & bfSolutions$sig4b >= minSignifLevel] <- TRUE
    } else {
      bfSolutions$keep[bfSolutions$cost < 0 & bfSolutions$duration < 0 & bfSolutions$sig3b >= minSignifLevel & bfSolutions$sig4b >= minSignifLevel] <- TRUE
    }
  }
  
  if (nbVariables == 3) {
    if (withSignificantIntercepts) {
      bfSolutions$keep[bfSolutions$cost < 0 & bfSolutions$duration < 0 & bfSolutions$length < 0 & bfSolutions$sig1b >= minSignifLevel & bfSolutions$sig2b >= minSignifLevel & bfSolutions$sig3b >= minSignifLevel & bfSolutions$sig4b >= minSignifLevel & bfSolutions$sig5b >= minSignifLevel] <- TRUE
    } else {
      bfSolutions$keep[bfSolutions$cost < 0 & bfSolutions$duration < 0 & bfSolutions$length < 0 &bfSolutions$sig3b >= minSignifLevel & bfSolutions$sig4b >= minSignifLevel & bfSolutions$sig5b >= minSignifLevel] <- TRUE
    }
  }
  
  bfSolutions$sig1b <- NULL
  bfSolutions$sig2b <- NULL
  bfSolutions$sig3b <- NULL
  bfSolutions$sig4b <- NULL
  bfSolutions$sig5b <- NULL
  
  return (bfSolutions)
}
