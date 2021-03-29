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
# The script computes a logit model for each possible (combination of) lambda(s) to use for the Box-Cox transform in a given
# range and for a given step (granularity).
#
# For each model, the log-likelihood, the values of the estimated parameters and their level of significance are stored, along
# with an error code, if any.
#

# Clear memory if wanted
rm(list = ls(all.names = TRUE))

# Model to solve (input data to use)
modelName <- "europeL2"

# Number of variables that will be Box-Cox transformed
nbVariables <- 2

# Range and step size (granularity) for the Box-Cox lambda values to test
range <- 2.0
granularity <- 0.1

######################################################################
#                      Main entry point                          #####
######################################################################

# Change working directory to the location of this script
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# import the source code of the solver
source("_BoxCoxLogitSolver.R")

# import some convenience functions
source("_Utils.R")

# Change the output width
options(width = 200)

# Create a dataframe that will contain the solutions for all the possible (combinations of) lambdas
# and a given group.
# The table contains a row for each lambda's combination. Each row contains the characteristics of
# a solution (LL, values of the estimated parameters and their signif. level)
createSolutionsTable <- function(group, nbVariables, range, granularity) {

  # Number of steps in the range of values to test
  nbSteps <- 1 + 2 * range / granularity # Number of steps in the range of values to test
  nbCombinations <- nbSteps^nbVariables

  # Fill a matrix all possible lambda combinations
  m <- matrix(ncol = nbVariables, nrow = nbCombinations)
  mode(m) <- "numeric"
  for (col in 1:nbVariables) {
    nbOccurences <- nbCombinations / nbSteps^col
    row <- 1
    value <- -range
    while (row <= nbCombinations) {
      for (j in 1:nbOccurences) {
        m[row, col] <- value
        row <- row + 1
      }

      # Next value
      value <- round(value + granularity, 1)
      if (value > range) {
        value <- -range
      }
    }
  }

  # Transform to data frame
  solutionsTable <- data.frame(m)

  names(solutionsTable)[1] <- "lambda.cost"
  if (nbVariables > 1) names(solutionsTable)[2] <- "lambda.duration"
  if (nbVariables > 2) names(solutionsTable)[3] <- "lambda.length"

  solutionsTable$LL <- NA
  solutionsTable$error <- ""

  s <- "const.iww"
  solutionsTable[, s] <- NA
  s <- "const.rail"
  solutionsTable[, s] <- NA

  solutionsTable$cost <- NA
  if (nbVariables > 1) solutionsTable$duration <- NA
  if (nbVariables > 2) solutionsTable$length <- NA

  s <- "sig.const.iww"
  solutionsTable[, s] <- ""
  s <- "sig.const.rail"
  solutionsTable[, s] <- ""

  s <- "sig.cost"
  solutionsTable[, s] <- ""

  if (nbVariables > 1) {
    s <- "sig.duration"
    solutionsTable[, s] <- ""
  }

  if (nbVariables > 2) {
    s <- "sig.length"
    solutionsTable[, s] <- ""
  }


  # Add two first column to store the group id and a flag to mark the best solution
  solutionsTable <- cbind(data.frame(matrix(FALSE, nrow = nrow(solutionsTable), ncol = 1)), solutionsTable)
  solutionsTable <- cbind(data.frame(matrix(NA, nrow = nrow(solutionsTable), ncol = 1)), solutionsTable)
  names(solutionsTable)[1] <- "group"
  names(solutionsTable)[2] <- "best"

  # Set the group of commodities that must be solved
  solutionsTable$group <- group
  
  return(solutionsTable)
}

# Dataset that contains the explanatory variables
inputData <- NULL
inputFile <- paste(modelName, "-Input.Rda", sep = "")
load(inputFile)

# Dataframe that will contain all the computed models
bfSolutions <- data.frame()

# Get a vector with the available groups
groups <- as.numeric(levels(as.factor(inputData$grp)))

# Force a subset of groups
#groups <- c(0)

# Use parallel computing in mnLogit. using 2 cores seems to be the most efficient
#nbCores <- parallel:::detectCores()
nbCores <- 2

# Solve for each group
for (i in 1:length(groups)) {
  # Get the data for current group
  group <- groups[i]
  inputDataForGroup <- inputData[inputData$grp == group, ]

  # Dataframe that will contain all the computed models  for this group
  solutionsForGroup <- createSolutionsTable(group, nbVariables, range, granularity)
  
  # Solve for each (combination of) lambda(s)
  n <- nrow(solutionsForGroup)
  
  for (r in 1:n) {
    # Get lambda combination for this row
    cat(paste("Solving logit for group", group, "and lambda(s) "))
    lambdas <- c()
    for (j in 1:nbVariables) {
      # Lambdas start at column 3 in table
      cat(paste(solutionsForGroup[r, j + 2], " "))
      lambdas <- c(lambdas, solutionsForGroup[r, j + 2])
    }
    cat("\n")

    # As some lambda's can lead to numerical singularity, intercept error.
    res <- try(
      {
        # Solve logit
        model <- solveBoxCoxLogit(inputDataForGroup, lambdas, nbCores)
        
        # Retrieve and store results
        solutionsForGroup[r, "LL"] <- model$logLik

        s <- "(Intercept):2"
        solutionsForGroup[r, "const.iww"] <- coef(model)[s]
        solutionsForGroup[r, "sig.const.iww"] <- getStars(model, s)

        s <- "(Intercept):3"
        solutionsForGroup[r, "const.rail"] <- coef(model)[s]
        solutionsForGroup[r, "sig.const.rail"] <- getStars(model, s)

        solutionsForGroup[r, "cost"] <- coef(model)["cost"]
        solutionsForGroup[r, "sig.cost"] <- getStars(model, "cost")

        if (nbVariables > 1) {
          solutionsForGroup[r, "duration"] <- coef(model)["duration"]
          solutionsForGroup[r, "sig.duration"] <- getStars(model, "duration")
        }

        if (nbVariables > 2) {
          solutionsForGroup[r, "length"] <- coef(model)["length"]
          solutionsForGroup[r, "sig.length"] <- getStars(model, "length")
        }

        # Is this a valid solution ?
        if (allSignsAreExpected(model)) {
          error <- ""
        } else {
          error <- "Wrong Sign(s)"
        }

        solutionsForGroup[r, "error"] <- error
      },
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      solutionsForGroup[r, "error"] <- "Singularity"
    }
  }
  
  # Identify best valid solution
  tmp <- subset(solutionsForGroup, error == "")
  bestLL <- max(tmp$LL)
  solutionsForGroup[which(solutionsForGroup$LL == bestLL), "best"] <- TRUE
  
  # Merge with solutions for other groups
  bfSolutions <- rbind(bfSolutions, solutionsForGroup)
  
  # Next group
  cat("\n")
}

# Print best solutions
bestSolutions <- bfSolutions[bfSolutions$best == TRUE, ]
bestSolutions$best <- NULL
bestSolutions$error <- NULL
cat("Best solutions (correct signs but regardless of signif. levels of estimators):\n")
print(bestSolutions, row.names = FALSE)

# Remove the "best" column before saving
bfSolutions$best <- NULL

# Save all the results (further used for plotting and comparison with heuristic solutions)
z <- paste("0", 10 * granularity, sep = "")
outputFile <- paste(modelName, "-", z, "-BruteForce", nbVariables, ".Rda", sep = "")
save(bfSolutions, file = outputFile)
