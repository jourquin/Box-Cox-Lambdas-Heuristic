# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# A few convenience functions used in more than one script in this project
#


# Test if all the estimators are of the expected sign (must be negative)
allSignsAreExpected <- function(model) {
  c = coef(model)
  correctSign = TRUE
  # Browse de coefficients names (see output of "summary(model)")
  for (j in 1:length(c)) {
    name = names(c[j])
    if (substring(name, 1, 1) != "(") {
      # "(Intercept)" must not be tested
      if (c[name] > 0) {
        correctSign = FALSE
        break
      }
    }
  }
  return(correctSign)
}

# Get Pr(>|t|) of model for all coefficients
allCoefsAreSignificant <- function(object, nbStars=3)
{
  b <- coef(object, order=TRUE)
  std.err2 <- sqrt(diag(vcov(object)))
  std.err <- b
  std.err[names(std.err2)] <- std.err2
  
  z <- b / std.err
  p <- 2 * ( 1 - pnorm(abs(z)))
  
  for (j in 1:length(p)) {
    if (nbStars == 3 && p[j] > 0.001 ) return(FALSE)
    if (nbStars == 2 && p[j] > 0.01 ) return(FALSE)
    if (nbStars == 1 && p[j] > 0.05 ) return(FALSE)
    if (nbStars == 0 && p[j] > 0.1 ) return(FALSE)
    if (nbStars == 0 && p[j] > 1 ) return(FALSE)
  }
  return(TRUE)
}

# Get Pr(>|t|) of model for a given coefficient
getStars <- function(object, coefName)
{
  b <- coef(object, order=TRUE)
  std.err2 <- sqrt(diag(vcov(object)))
  std.err <- b
  std.err[names(std.err2)] <- std.err2
  
  z <- b / std.err
  p <- 2 * ( 1 - pnorm(abs(z)))
  
  pp = p[coefName]
  if (pp <= 0.001 ) return("***")
  if (pp <= 0.01 ) return("**")
  if (pp <= 0.05 ) return("*")
  if (pp <= 0.1 ) return(".")
  if (pp <= 1 ) return(" ")
}

# Returns TRUE if signs are expected and  all the estimators are enough significant (no error)
isValid <- function(solution) {
  if (solution$error == "") return(TRUE)
  return(FALSE)
}

# Returns TRUE if all the signs are expected (some estimators may have lower significance)
hasExpectedSigns <- function(solution) {
  if (solution$error == "" | solution$error == "Low signif.") return(TRUE)
  return(FALSE)
}


# Create a dataframe that will contain the solutions for all the possible (combinations of) lambdas.
# The table contains a row for each lambda's combination. Each row contains the characteristics of 
# a solution (LL, values of the estimated parameters and their signif. level)
createSolutionsTable <- function(nbVariables, range, granularity) {
  
  # Number of steps in the range of values to test
  nbSteps = 1 + 2*range / granularity # Number of steps in the range of values to test
  nbCombinations = nbSteps^nbVariables
  
  # Fill a matrix all possible lambda combinations
  m = matrix(ncol = nbVariables, nrow = nbCombinations)
  mode(m) = "numeric"
  for (col in 1:nbVariables) {
    nbOccurences = nbCombinations / nbSteps^col
    row = 1
    value = -range
    while (row <= nbCombinations) {
      for (j in 1: nbOccurences) {
        m[row, col] <- value
        row = row + 1
      }
      
      # Next value
      value = round(value + granularity, 1)
      if (value > range) {
        value = -range
      }
    }
  }
  
  # Transform to data frame
  solutionsTable <- data.frame(m)
  
  names(solutionsTable)[1] = "lambda.cost"
  if (nbVariables > 1) names(solutionsTable)[2] = "lambda.duration"
  if (nbVariables > 2) names(solutionsTable)[3] = "lambda.length"
  
  solutionsTable$LL = NA
  solutionsTable$error = ""
  
  s = "const.iww"
  solutionsTable[, s] = NA
  s = "const.rail"
  solutionsTable[, s] = NA
  
  solutionsTable$cost = NA
  if (nbVariables > 1) solutionsTable$duration = NA
  if (nbVariables > 2) solutionsTable$length = NA
  
  s = "sig.const.iww"
  solutionsTable[, s] = ""
  s = "sig.const.rail"
  solutionsTable[, s] = ""
  
  s = "sig.cost"
  solutionsTable[, s] = ""
  
  if (nbVariables > 1) {
    s = "sig.duration"
    solutionsTable[, s] = ""
  }
  
  if (nbVariables > 2) {
    s = "sig.length"
    solutionsTable[, s] = ""
  }
  
  
  # Add two first column to store the group id and a flag to mark the best solution
  solutionsTable = cbind(data.frame(matrix(FALSE, nrow = nrow(solutionsTable), ncol = 1)), solutionsTable)
  solutionsTable = cbind(data.frame(matrix(NA, nrow = nrow(solutionsTable), ncol = 1)), solutionsTable)
  names(solutionsTable)[1] = "group"
  names(solutionsTable)[2] = "best"
  
  return(solutionsTable)
}

# Compute index of a given lambda combination in the solutions table
getIdx <- function(lambdas, range, granularity) {
  l = length(lambdas)
  
  # Replace lambda values with their index in the range
  for (i in 1:l) {
    lambdas[i] = ((range + lambdas[i]) / granularity) + 1
  }
  
  # Compute the index of this combination into the solutions table 
  idx = 0
  z = l-1
  for (i in 1:l) {
    if (i == l) {
      idx = idx + lambdas[i] 
    } else {
      idx = idx + (lambdas[i]-1)*nbSteps^z   
    }
    z = z-1
  }
  
  return (round(idx,0))
}

# Returns the lambdas of a given solution
getLambdas <- function(solution, nbVariables) {
  lambdas = c()
  for (i in 1:nbVariables) {
    if (i == 1) s = "solution$lambda.cost"
    if (i == 2) s = "solution$lambda.duration"
    if (i == 3) s = "solution$lambda.length"
    value <<- eval(parse(text = s))
    lambdas = c(lambdas,value)
  } 
  return(lambdas)
}

# Draw a random combination of lambda's
randomDrawLambdas <- function(nbLambdas, range, granularity) {
  lambdas <- c()
  for (j in 1:nbLambdas) {
    z = sample(1:nbSteps, 1)
    lambda = -range - granularity + (z*granularity)
    lambdas = c(lambdas, round(lambda, 1))
  }
  return(lambdas)
}
