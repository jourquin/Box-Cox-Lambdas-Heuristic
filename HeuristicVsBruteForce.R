# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# R script that compares the results of the heuristic to the optimal solutions obtained using the brute force approach.
# The output table gives, for each group of commodities :
# - The Log-Likelihood of the two solutions (exact and heuristic)
# - The number of logit computations needed to obtain both solutions ("nb.exact" and "nb.heur")
# - The number of "checks" (retrieval of an already computed solution) needed by the heuristic.
# - The number of "valid" solutions identified by the brute force approach
# - The rank of the solution found by the heuristic in all the possible "valid" solutions
# - The relative difference between the highest and lowest Log-Likelihood among all the "valid" solutions
# - The relative difference between the Log-Likelihood of the exact solution and the one of the heuristic
#

# Clear memory if wanted
rm(list = ls(all.names = TRUE))

# Number of explanatory variables (1, 2 or 3)
nbVariables = 2

# Successive steps to use to reach the final granularity
steps = c(0.4, 0.2, 0.1)

# Number of random combinations to draw during the first step of the heuristic (for 1, 2 or 3 lambdas)
nbRandomCombinations = c(1, 5, 15)

# Number of heuristic repetitions for one group (shotgun hill climbing)
nbShotguns = 1

# A random seed can be set here (for replicability). If seed = -1, no seed is set
#seed = -1
seed = 2020

# Maximum number of Logit computations before stop if no convergence to a solution
maxLogitComputations = 500*nbShotguns

# Maximum number of solution checks (retrieval of stored results) before stop if no convergence
maxChecks = 2000*nbShotguns


# Use the previously solved solutions using the brute force method
modelName = "europeL2-01"

######################################################################
# Nothing should be changed beyond this line                     #####
######################################################################

# import some convenience functions
source("_Utils.R")

# This script uses the already computed solutions by the brute force algorithm
# instead of using a real solver
solveLogit <- function(idx, lambdas) {
  nbLogitComputations <<- nbLogitComputations + 1
  
  # Just retrieve an already computed solution using a brute force approach
  solution = bfSolutionsForGroup[idx,]
  
  # Copy the results
  for (i in 1:ncol(solution)) {
    s = solution[1, i]
    if (nlevels(s) > 0) s = paste(s) # paste needed, as a factor level is returned otherwise
    solutionsForGroup[idx, i] <<- s 
  }
  solutionsForGroup$isSolved[idx] <<- TRUE
}


# Get the solution for a given lambda combination. Compute the model if not yet done.
retrieveOrComputeLogit <- function(lambdas) {
  
  # The index can be computed from the lambda combination. 
  idx = getIdx(lambdas, range, granularity)
  
  # Solve it not yet done
  if (!solutionsForGroup$isSolved[idx]) {
    solveLogit(idx, lambdas)
  } else {
    nbModelChecks <<- nbModelChecks + 1
  }
  
  return(solutionsForGroup[idx, ])
}


# Explore the lambda combinations arround a given solution
exploreAround <- function(solution, stepSize, dimLevel = 1) {
  
  # Test a step backward in the current dimension, but remain in range
  p = getLambdas(solution, nbVariables)
  if (p[dimLevel] - stepSize >= -range) {
    p[dimLevel] = round(p[dimLevel] - stepSize, 1)
  } else {
    p[dimLevel] = -range
  }
  newSolution = retrieveOrComputeLogit(p)
  
  # Is this the best solution found so far?
  if (isValid(newSolution) & newSolution$LL > bestSolution$LL) {
    bestSolution <<- newSolution
  }
  
  # Add this solution to the list to explore later if it has better LL and expected signs
  if (hasExpectedSigns(newSolution) & newSolution$LL > solution$LL) {
    solutionsToExplore[[length(solutionsToExplore)+1]] <<- newSolution
  }
  
  # Test a step forward in the current dimension, but remain in range
  p = getLambdas(solution, nbVariables)
  if (p[dimLevel] + stepSize <= range) {
    p[dimLevel] = round(p[dimLevel] + stepSize, 1)
  } else {
    p[dimLevel] = range
  }
  newSolution = retrieveOrComputeLogit(p)
  
  # Is this the best solution found so far?
  if (isValid(newSolution) & newSolution$LL > bestSolution$LL) {
    bestSolution <<- newSolution
  }
  
  # Add this solution to the list to explore later if it has better LL and expected signs
  if (hasExpectedSigns(newSolution) & newSolution$LL > solution$LL) {
    solutionsToExplore[[length(solutionsToExplore)+1]] <<- newSolution
  }
  
  # Recursive entrance to next dimension
  if (dimLevel < length(p)) {
    exploreAround(solution, stepSize, dimLevel + 1)
  } 
}

#########################################################################################
# MAIN ENTRY 
#########################################################################################

# Set working dir to location of this script
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
options(width=250)

# Load brute force results (bfSolutions)
load(paste(modelName, "-BruteForce", nbVariables, ".Rda", sep=""))

# Get the groups of commodities present in the input data
groups = as.numeric(levels(as.factor(bfSolutions$group)))

# Get the range of lambdas
range = max(bfSolutions$lambda.cost)

# Get the granularity
granularity = steps[length(steps)]

# Number of steps in the range of values to test
nbSteps = 1 + 2*range / granularity # Number of steps in the range of values to test

# Force a subset of groups
#groups = c(2)

# Initialize a dataframe that will contain the results (exact vs heuristic)
n = nbVariables*2 + 10
output <- data.frame(matrix(ncol = n, nrow = length(groups)))
s = c("group", "LL.exact", "LL.heur", "nb.exact", "nb.heur", "nb.checks", "nb.solutions", "rank", "max.diff.LL", "diff.LL")
for (k in 1:nbVariables) {
  if (k == 1) s = append(s, "lambda.cost.exact")
  if (k == 2) s = append(s, "lambda.duration.exact")
  if (k == 3) s = append(s, "lambda.length.exact")
}

for (k in 1:nbVariables) {
  if (k == 1) s = append(s, "lambda.cost.heur")
  if (k == 2) s = append(s, "lambda.duration.heur")
  if (k == 3) s = append(s, "lambda.length.heur")
}

colnames(output) <- s

# Dataframe that will contain all the computed solutions by the heuristic
heuristicPath = data.frame()

# Solve the problem separately for each group of commodities
for (g in 1:length(groups)) {
  group = groups[g]
  bfSolutionsForGroup = bfSolutions[bfSolutions$group == group, ]
  
  # Force seed for replicability
  if (seed != -1) {
    set.seed(seed + group)
  }
  
  cat(paste("Solving for group", group, "\n"))
  
  # Track the already tested lambda combinations in a table
  solutionsForGroup = createSolutionsTable(nbVariables, range, granularity)
  solutionsForGroup$group = group
  solutionsForGroup$isSolved = FALSE
  
  # Counters
  nbLogitComputations = 0
  nbModelChecks = 0
  
  # Initial values
  bestSolution = NULL
  solutionsToExplore = list()
  
  for (shotgun in 1:nbShotguns) {
    ###################################################################################
    # 1. Randomly choose a series of lambdas and keep the one that gives the highest LL
    ###################################################################################
    
    # Loop until the requested random combinations is reached
    nbCombinationsDrawn = 0
    repeat {
      # Randomly choose a combination of lambdas
      lambdas = randomDrawLambdas(nbVariables, range, granularity)
      
      # Compute the logit for this combination
      solution = retrieveOrComputeLogit(lambdas)
      
      # If solution for this combination is valid
      if (isValid(solution)) {
        # Keep this solution if better
        if (is.null(bestSolution)) {
          bestSolution = solution
        } else if (solution$LL > bestSolution$LL) {
          bestSolution = solution
        }
        
        # Break loop if number of valid drawns is reached
        nbCombinationsDrawn = nbCombinationsDrawn + 1
        if (nbCombinationsDrawn == nbRandomCombinations[nbVariables]) break
      }
      
      # Break if max attempts has been reached
      if (nbLogitComputations >= maxLogitComputations || nbModelChecks >= maxChecks) break
    }
    
    
    ###################################################################################
    # 2. Explore around this combination and try to increase the LL
    ###################################################################################
   
    if (nbLogitComputations < maxLogitComputations && nbModelChecks < maxChecks) {
      
      # The list of combinations to test will dynamically grow/shrink
      solutionsToExplore[[length(solutionsToExplore)+1]] = bestSolution
      
      # Repeat the procedure until the list contains no points anymore
      for (j in 1:length(steps)) {
        currentStep = steps[j]
        repeat {
          # Stop heuristic if one of the thresholds is reached
          if (nbLogitComputations >= maxLogitComputations || nbModelChecks >= maxChecks) break
          
          # Explore solutions arround current combination
          exploreAround(solutionsToExplore[[1]], currentStep)
          
          # Break if only one (and last) solution remains in the list
          if (length(solutionsToExplore) == 1) {
            break
          }
          
          # Remove the current solution from the list
          solutionsToExplore = solutionsToExplore[-1]
        }
      }
    }
  }
  
  # Mark solution found by the heuristic
  lambdas = getLambdas(bestSolution, nbVariables)
  solutionsForGroup[getIdx(lambdas, range, granularity), "best"] = TRUE

  
  # Rank of this solution
  LLs = sort(bfSolutionsForGroup$LL[bfSolutionsForGroup$error == ""], decreasing = TRUE)
  for (j in 1:length(LLs)) {
    if (LLs[j] == bestSolution$LL) rank = j
  }
  
  # Difference between best and worse LL of the valid solutions
  maxDiffLL = round((LLs[1] - LLs[length(LLs)]) / LLs[1], 2)
  diffLL = round((LLs[1] - bestSolution$LL) / LLs[1], 2)
  
  # Target LL to reach
  tmp = bfSolutionsForGroup[bfSolutionsForGroup$error == "", ]
  veryBestLL = max(tmp$LL) 
  
  # Fill the row of the output dataframe corresponding to the current group
  output[g, "group"] = group
  output[g, "nb.exact"] = (1 + 2*range/granularity)^nbVariables
  output[g, "nb.heur"] = nbLogitComputations
  output[g, "nb.checks"] = nbModelChecks
  output[g, "nb.solutions"] = length(LLs)
  output[g, "rank"] = rank
  output[g, "max.diff.LL"] = paste((maxDiffLL*100),"%", sep="")
  output[g, "diff.LL"] = paste((diffLL*100),"%", sep="")
  
  idx1 = which(bfSolutionsForGroup$LL == veryBestLL)
  idx2 = which(bfSolutionsForGroup$LL == bestSolution$LL)
  
  output[g, "LL.exact"] = bfSolutionsForGroup$LL[idx1]
  output[g, "LL.heur"] = bfSolutionsForGroup$LL[idx2]
  
  output[g, "lambda.cost.exact"] = bfSolutionsForGroup$lambda.cost[idx1]
  output[g, "lambda.cost.heur"] = bfSolutionsForGroup$lambda.cost[idx2]
  
  if (nbVariables > 1) {
    output[g, "lambda.duration.exact"] = bfSolutionsForGroup$lambda.duration[idx1]
    output[g, "lambda.duration.heur"] = bfSolutionsForGroup$lambda.duration[idx2]
  }
  
  if (nbVariables > 2) {
    output[g, "lambda.length.exact"] = bfSolutionsForGroup$lambda.length[idx1]
    output[g, "lambda.length.heur"] = bfSolutionsForGroup$lambda.length[idx2]
  }
  
  # Add the solutions computed for this group to the solutions
  heuristicPath = rbind(heuristicPath, solutionsForGroup)
}

# Save the solutions computing by the heuristic (further used for plotting)
save(heuristicPath,file=paste(modelName, "-Heuristic", nbVariables, ".Rda", sep=""))

# Print the output table
cat("\nComparison between exact (brute force) and heuristic solutions:\n")
print(output, row.names = FALSE)
