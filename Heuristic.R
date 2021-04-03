# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# R script that implements a weighted conditional multinomial logit with a Box-Cox transform for 1, 2 or 3 explanatory
# variables when the following input data is available for an origin-destination pair:
# - Origin and a destination ID's (OD relation)
# - Group of commodities (the logit is solved separately for each group)
# - Unit cost for each mode on the OD relation
# - Distance for each mode on the OD relation
# - Length for each mode on the OD relation
# - Observed demand (tons per year) for each mode between O and D
#
# The provided input dataset is related to freight transport, with 3 modes (1: road, 2: inland waterways, 3: rail). Not all modes
# are available between all OD pairs (NA values in the dataframe). Unit costs are expressed in euros/ton, transit times
# (durations) in hours and lengths in kilometers.
#
# The script implementents a meta heuristic able to quickly identify good (combination of) lambda(s) to use for the Box-Cox
# transforms in a given range and for a given step (granularity).
#
# Output: identifed lambda's, log-likelihood, values of the estimated parameters and their level of significance of the logit model.
#

# Clear memory if wanted
rm(list = ls(all.names = TRUE))

# Set working dir to location of this script
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
options(width = 200)

# Input data
modelName <- "europeL2-Input.Rda"

# Group of commodities for which the logit must be estimated
group <- 0

# Number of independant variables that must be Box-Cox transformed
nbVariables <- 2

# Min significant level to accept for all the estimators (valid solutions)
minSigLevel <- 1

# If TRUE, the signif. level of the intercepts is also tested
withSignificantIntercepts <- FALSE

# Number of random combinations to draw during the first step of the heuristic (for 1, 2 or 3 lambdas)
nbRandomCombinations <- c(1, 5, 15)

# Number of heuristic repetitions for one group (shotgun hill climbing). Depends on nbVariables
nbShotguns <- c(1, 5, 10)

# Range to look in (from -range to +range)
range <- 2

# Step to use between -range and +range
granularity <- 0.1

# Successive steps to use to reach the final granularity
steps <- c(0.4, 0.2, 0.1)

# A random seed can be set here (for replicability). If seed = -1, no seed is set
seed <- -1
#seed <- 2020

# Maximum number of Logit computations before stop if no convergence to a solution
maxLogitComputations <- 500 * nbShotguns

# Maximum number of solution checks (retrieval of stored results) before stop if no convergence
maxChecks <- 2000 * nbShotguns


######################################################################
# Nothing should be changed beyond this line                     #####
######################################################################


# import the source code of the solver
source("_BoxCoxLogitSolver.R", local = TRUE)

# import some convenience functions
source("_Utils.R", local = TRUE)

# Scope level variables
nbLogitComputations <- 0
nbModelChecks <- 0
bestSolution <- NULL
solutionsToExplore <- list()
storedSolutions <- new.env(hash=TRUE) 

# Solve the logit for a given combination of lambda's, that will be stored in row idx of the solution table
solveLogit <- function(lambdas) {
  # Prepare empty dataframe
  colNames <- c(
    "group",
    "lambda.cost",
    "lambda.duration",
    "lambda.length",
    "LL",
    "const.iww",
    "sig.const.iww",
    "const.rail",
    "sig.const.rail",
    "cost",
    "sig.cost",
    "duration",
    "sig.duration",
    "length",
    "sig.length",
    "error"
  )
  solution <- data.frame(matrix(ncol = length(colNames), nrow = 1))
  solution <- setNames(solution, colNames)

  # As some lambda's can lead to numerical singularity, intercept error.
  res <- try(
    {
      # Solve logit (inputDataForGroup and nbCores are global variables, as R doesn't pass parameters as references...)
      model <- solveBoxCoxLogit(inputDataForGroup, lambdas, nbCores)
      #print(summary(model))

      # Retrieve and store results
      solution[1, "group"] <- group

      solution[1, "LL"] <- model$logLik

      s <- "(Intercept):2"
      solution[1, "const.iww"] <- coef(model)[s]
      solution[1, "sig.const.iww"] <- getStars(model, s)

      s <- "(Intercept):3"
      solution[1, "const.rail"] <- coef(model)[s]
      solution[1, "sig.const.rail"] <- getStars(model, s)

      s <- "cost"
      solution[1, "cost"] <- coef(model)[s]
      solution[1, "sig.cost"] <- getStars(model, s)
      solution[1, "lambda.cost"] <- lambdas[1]

      if (nbVariables > 1) {
        s <- "duration"
        solution[1, "duration"] <- coef(model)[s]
        solution[1, "sig.duration"] <- getStars(model, s)
        solution[1, "lambda.duration"] <- lambdas[2]
      }

      if (nbVariables > 2) {
        s <- "length"
        solution[1, "length"] <- coef(model)[s]
        solution[1, "sig.length"] <- getStars(model, s)
        solution[1, "lambda.length"] <- lambdas[3]
      }

      # Is this a valid solution ?
      if (allSignsAreExpected(model)) {
        if (allCoefsAreSignificant(model, minSigLevel, withSignificantIntercepts)) {
          error <- ""
        } else {
          error <- "Low signif."
        }
      } else {
        error <- "Wrong Sign(s)"
      }
    },
    silent = TRUE
  )
  if (inherits(res, "try-error")) {
    error <- "Singularity"
  }

  solution[1, "error"] <- error

  return(solution)
}


# Get the solution for a given lambda combination. Compute the model if not yet done.
retrieveOrComputeLogit <- function(lambdas) {
  solution <- storedSolutions[[getKey(lambdas)]]
  if (is.null(solution)) {
    cat("C")
    nbLogitComputations <<- nbLogitComputations + 1
    solution <- solveLogit(lambdas)
    storedSolutions[[getKey(lambdas)]] <<- solution
  } else {
    cat(".")
    nbModelChecks <<- nbModelChecks + 1
  }

  return(solution)
}


# Explore the lambda combinations arround a given solution
exploreAround <- function(solution, stepSize, dimLevel = 1) {
  # Test a step backward in the current dimension, but remain in range
  p <- getLambdas(solution, nbVariables)
  if (p[dimLevel] - stepSize >= -range) {
    p[dimLevel] <- round(p[dimLevel] - stepSize, 1)
  } else {
    p[dimLevel] <- -range
  }
  newSolution <- retrieveOrComputeLogit(p)

  # Is this the best solution found so far?
  if (isValid(newSolution) & newSolution$LL > bestSolution$LL) {
    bestSolution <<- newSolution
  }

  # Add this solution to the list to explore later if it has better LL and expected signs
  if (hasExpectedSigns(newSolution) &
    newSolution$LL > solution$LL) {
    solutionsToExplore[[length(solutionsToExplore) + 1]] <<- newSolution
  }

  # Test a step forward in the current dimension, but remain in range
  p <- getLambdas(solution, nbVariables)
  if (p[dimLevel] + stepSize <= range) {
    p[dimLevel] <- round(p[dimLevel] + stepSize, 1)
  } else {
    p[dimLevel] <- range
  }
  newSolution <- retrieveOrComputeLogit(p)

  # Is this the best solution found so far?
  if (isValid(newSolution) & newSolution$LL > bestSolution$LL) {
    bestSolution <<- newSolution
  }

  # Add this solution to the list to explore later if it has better LL and expected signs
  if (hasExpectedSigns(newSolution) &
    newSolution$LL > solution$LL) {
    solutionsToExplore[[length(solutionsToExplore) + 1]] <<- newSolution
  }

  # Recursive entrance to next dimension
  if (dimLevel < length(p)) {
    exploreAround(solution, stepSize, dimLevel + 1)
  }
}


#########################################################################################
# MAIN ENTRY
#########################################################################################

# Load input data needed to solve the logit for the group of commodities
inputData <- NULL
load(modelName)
inputDataForGroup <- inputData[inputData$grp == group, ]

# Number of steps in the range of values to testwarni wa
nbSteps <- 1 + 2 * range / granularity # Number of steps in the range of values to test

# Use parallel computing in mnLogit. Using 2 cores seems to be the most efficient
#nbCores <- parallel:::detectCores()
nbCores <- 2

# Force seed for replicability
if (seed != -1) {
  set.seed(seed)
}

cat(paste("Solving for group", group, "\n"))

for (shotgun in 1:nbShotguns[nbVariables]) {
  ###################################################################################
  # 1. Randomly choose a series of lambdas and keep the one that gives the highest LL
  ###################################################################################

  # Loop until the requested random combinations is reached
  nbCombinationsDrawn <- 0
  repeat {
    # Randomly choose a combination of lambdas
    lambdas <- randomDrawLambdas(nbVariables, range, granularity)

    # Compute the logit for this combination
    solution <- retrieveOrComputeLogit(lambdas)

    # If solution for this combination is valid
    if (isValid(solution)) {
      # Keep this solution if better∑
      if (is.null(bestSolution)) {
        bestSolution <- solution
      } else if (solution$LL > bestSolution$LL) {
        bestSolution <- solution
      }

      # Break loop if number of valid drawns is reached
      nbCombinationsDrawn <- nbCombinationsDrawn + 1
      if (nbCombinationsDrawn == nbRandomCombinations[nbVariables]) {
        break
      }
    }

    # Break if max attempts has been reached
    if (nbLogitComputations >= maxLogitComputations ||
      nbModelChecks >= maxChecks) {
      break
    }
  }
  cat(" - ")

  ###################################################################################
  # 2. Explore around this combination and try to increase the LL
  ###################################################################################

  if (nbLogitComputations < maxLogitComputations &&
    nbModelChecks < maxChecks) {
    # The list of combinations to test will dynamically grow/shrink
    solutionsToExplore[[length(solutionsToExplore) + 1]] <- bestSolution

    # Repeat the procedure until the list contains no points anymore
    for (j in 1:length(steps)) {
      currentStep <- steps[j]
      repeat {
        # Stop heuristic if one of the thresholds is reached
        if (nbLogitComputations >= maxLogitComputations ||
          nbModelChecks >= maxChecks) {
          break
        }

        # Explore solutions arround current combination
        exploreAround(solutionsToExplore[[1]], currentStep)

        # Break if only one (and last) solution remains in the list
        if (length(solutionsToExplore) == 1) {
          break
        }

        # Remove the current solution from the list
        solutionsToExplore <- solutionsToExplore[-1]
      }
    }
  }
}


# Print solution found by the heuristic
if (nbVariables < 3) {
  bestSolution$lambda.length <- NULL
  bestSolution$length <- NULL
  bestSolution$sig.length <- NULL
}

if (nbVariables < 2) {
  bestSolution$lambda.duration <- NULL
  bestSolution$duration <- NULL
  bestSolution$sig.duration <- NULL
}

cat(
  paste(
    "\n",
    nbLogitComputations,
    " Logit computations were needed instead of ",
    nbSteps^nbVariables,
    "\n",
    sep = ""
  )
)
print(bestSolution, row.names = FALSE)
