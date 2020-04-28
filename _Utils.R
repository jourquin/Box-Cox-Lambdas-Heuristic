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
allCoefsAreSignificant <- function(object, nbStars = 3) {
  b <- coef(object, order = TRUE)
  std.err2 <- sqrt(diag(vcov(object)))
  std.err <- b
  std.err[names(std.err2)] <- std.err2

  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))

  for (j in 1:length(p)) {
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
getStars <- function(object, coefName) {
  b <- coef(object, order = TRUE)
  std.err2 <- sqrt(diag(vcov(object)))
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

# Returns TRUE if signs are expected and  all the estimators are enough significant (no error)
isValid <- function(solution) {
  if (solution$error == "") {
    return(TRUE)
  }
  return(FALSE)
}

# Returns TRUE if all the signs are expected (some estimators may have lower significance)
hasExpectedSigns <- function(solution) {
  if (solution$error == "" | solution$error == "Low signif.") {
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
