# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# R script that plots the Log-Likelihood of the logits computed for all the possible combinations of lambdas (trivariate case).
#
# The results are taken from the output of the brute force algorithm. Only the "valid" solutions are plot.
#
# The exact best solution is identified by a green dot.
#
# If wanted, the unconstrained max Log-Likelihood can be plot (in red). It only appears if it is different from the constrained
# best solution (green dot).
#
# The solutions tested by the heuristic can also be plotted, using the output of the "HeuristicVsBruteForce.R" script. The final
# solution of the heuristic is plot in blue if it is different from the (exact) best solution identified by the brute force
# approach.
#
# Predefined plot settings are provided for each group. These are specific for the provided dataset in this project.
#
# The resulting plot is interactive: the mouse can be used to zoom and modify the perpective.
#

library(rgl)
library(car)

# Clear memory if wanted
rm(list = ls(all.names = TRUE))

# Change working directory to the location of this script
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Group of commodities to plot
group <- 9


# Minimum signif. level for the estimators: " " = -1, ." = 0, "*" = 1, "**" = 2, "***" = 3
minSignifLevel <- 1

# If TRUE, the signif. level of the intercepts is also tested
withSignificantIntercepts <- FALSE

# Plot unconstrained max LL (red dot)
plotUnconstrainedMax <- TRUE

# Limit to 0.5 step. If TRUE, plotPathToSolution must be set to FALSE
limitGranularity <- TRUE

# Plot all the solutions tested by the heuristic. The best one will be in blue if different from the exact bes solution
# If set to TRUE, limitGranularity must be FALSE
plotPathToSolution <- FALSE

# input dataset
modelName <- "europeL2-01"

# Can be adjusted to improve output (depends on previous settings)
sphereSize <- 2
zoomLevel <- 0.12
theta <- 0
emphaseLevel <- 1


######################################################################
# Nothing should be changed beyond this line                     #####
######################################################################

# Size of the plot area and position of the title
windowSize <- 800
titleLine <- -50

# import some convenience functions
source("_Utils.R")

# Number of variables in the utility function of the logit model
nbVariables <- 3

# Get the data for the current group from the brute force results
bfSolutions <- NULL
load(paste(modelName, "-BruteForce3.Rda", sep = ""))

if (limitGranularity) {
  # Limit to 0.5 granularity
  bfSolutions$l.c <- abs(bfSolutions$lambda.cost) *2
  bfSolutions$l.d <- abs(bfSolutions$lambda.duration) *2
  bfSolutions$l.l <- abs(bfSolutions$lambda.length) *2
  bfSolutions <- bfSolutions[bfSolutions$l.c%%1 == 0 & bfSolutions$l.d%%1 == 0 & bfSolutions$l.l%%1 == 0, ]
  bfSolutions$l.c <- NULL
  bfSolutions$l.d <- NULL
  bfSolutions$l.l <- NULL
}

data <- bfSolutions[bfSolutions$group == group, ]

# Mark the solutions to retain
data <- markValidSolutions(data, nbVariables, minSigLevel, withSignificantIntercepts)

# Mark constrained max
data$best[data$LL == max(data$LL[data$keep])] <- TRUE

unconstrainedMaxLL <- max(data$LL, na.rm = TRUE)

# Get all the valid solutions, with or without the best unconstrained one
if (plotUnconstrainedMax) {
  # Keep valid LL's and the max LL without constraints
  LambdaCost <- data$lambda.cost[data$keep | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
  LambdaDuration <- data$lambda.duration[data$keep | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
  LambdaLength <- data$lambda.length[data$keep | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
  LogLik <- data$LL[data$keep | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
  best <- data$best[data$keep | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
} else {
  # Keep valid LL's
  LambdaCost <- data$lambda.cost[data$keep]
  LambdaDuration <- data$lambda.duration[data$keep]
  LambdaLength <- data$lambda.length[data$keep]
  LogLik <- data$LL[data$keep]
  best <- data$best[data$keep]
}

# Default color and size of the dots
sizes <- (LogLik - min(LogLik)) / min(LogLik)
colors <- rep("gray95", length(sizes))
emphased <- rep(FALSE, length(sizes))

# Force to use the -2, 2 range in all dimensions using invisible spheres
LambdaCost <- c(LambdaCost, -2)
LambdaDuration <- c(LambdaDuration, -2)
LambdaLength <- c(LambdaLength, -2)
sizes <- c(sizes, 0)
colors <- c(colors, "white")

LambdaCost <- c(LambdaCost, 2.0)
LambdaDuration <- c(LambdaDuration, 2.0)
LambdaLength <- c(LambdaLength, 2.0)
sizes <- c(sizes, 0)
colors <- c(colors, "white")

# Find best LL under constraints
constrainedMaxIdx <- which(best == TRUE)
constrainedMaxLL <- LogLik[constrainedMaxIdx]

# Identify path to best LL if wanted
if (plotPathToSolution) {
  # Load the solutions tested by the heuristic for the group
  inputFile <- paste(modelName, "-Heuristic3.Rda", sep = "")
  load(inputFile)
  solutions <- heuristicPath[heuristicPath$group == group, ]

  lambdaCost2 <- solutions$lambda.cost
  lambdaDuration2 <- solutions$lambda.duration
  lambdaLength2 <- solutions$lambda.length
  LogLik2 <- solutions$LL

  best2 <- solutions$best
  bestLL <- LogLik2[which(best2 == TRUE)]
  bestLLIdx <- which(LogLik == bestLL)
  
  # Mark all the tested solutions
  for (i in 1:length(LogLik2)) {
    idx <- which(LogLik == LogLik2[i])
    if (length(idx) == 0) {
      # When a non valid solutions is tested
      LambdaCost <- c(LambdaCost, lambdaCost2[i])
      LambdaDuration <- c(LambdaDuration, lambdaDuration2[i])
      LambdaLength <- c(LambdaLength, lambdaLength2[i])
      LogLik <- c(LogLik, LogLik2[i])
      colors <- c(colors, "black")
      emphased <- c(emphased, FALSE)
    } else {
      colors[idx] <- "black"
    }
  }

  # Recomput sizes of nodes
  LogLik[is.na(LogLik)] <- min(LogLik, na.rm = TRUE)
  sizes <- (LogLik - min(LogLik, na.rm = TRUE)) / min(LogLik, na.rm = TRUE)
  
  # Identify solution of heuristic
  colors[bestLLIdx] <- "blue"
  if (!emphased[bestLLIdx]) {
    sizes[bestLLIdx] <- emphaseLevel * sizes[bestLLIdx]
    emphased[bestLLIdx] <- TRUE
  }
}

# Identify max LL without constraint
if (plotUnconstrainedMax) {
  unconstrainedMaxIdx <- which(LogLik == unconstrainedMaxLL)
  colors[unconstrainedMaxIdx] <- "red"
  if (!emphased[unconstrainedMaxIdx]) {
    sizes[unconstrainedMaxIdx] <- emphaseLevel * sizes[unconstrainedMaxIdx]
    emphased[unconstrainedMaxIdx] <- TRUE
  }
}

# identify best valid solution
colors[constrainedMaxIdx] <- "green"
if (!emphased[constrainedMaxIdx]) {
  sizes[constrainedMaxIdx] <- emphaseLevel * sizes[constrainedMaxIdx]
}

# Display
open3d()
par3d(windowRect = c(0, 0, windowSize, windowSize))

# Labels with greek lambda
l1 <- parse(text = paste(intToUtf8(955), "[C]^", group, sep = ""))
l2 <- parse(text = paste(intToUtf8(955), "[T]^", group, sep = ""))
l3 <- parse(text = paste(intToUtf8(955), "[L]^", group, sep = ""))

# Hide the default axis of the scatter3d plot
axisCol <- rgb(255, 255, 255, max = 255, alpha = 0, names = "transparentWhite")

# Plot and adapt perspective to settings
scatter3d(x = LambdaCost, y = LambdaDuration, z = LambdaLength, fogtype = "exp2", surface = FALSE, sphere.size = sphereSize, radius = sizes, point.col = colors, axis.col = c(axisCol, axisCol, axisCol), xlab = l1, ylab = l2, zlab = l3)
rgl.viewpoint(theta = 0, phi = -90, fov = 170, zoom = zoomLevel, interactive = TRUE)
U <- par3d("userMatrix")
par3d(userMatrix = rotate3d(U, theta, 0, 0, 1))

# Add title to the bottom
title <- paste("NST-R", group)
bgplot3d({
  plot.new()
  title(main = title, line = titleLine, cex.main = 2)
})

# Add the central axis with their labels
axisCol <- "gray20"
labelCol <- "blue"
rgl::axis3d("x", pos = 0.5, lwd = 2, col = axisCol, labels = c("-2", "0", "2"), tick = FALSE, nticks = 2, cex = 1.2)
rgl::axis3d("y", pos = 0.5, lwd = 2, col = axisCol, labels = c("-2", "2"), tick = FALSE, nticks = 1, cex = 1.2)
rgl::axis3d("z", pos = 0.5, lwd = 2, col = axisCol, labels = c("-2", "2"), tick = FALSE, nticks = 1, cex = 1.3)
rgl::mtext3d(text = l1, edge = "x", pos = 0.5, col = labelCol, at = 1.1, cex = 1.5)
rgl::mtext3d(text = l2, edge = "y", pos = 0.5, col = labelCol, at = 1.1, cex = 1.5)
rgl::mtext3d(text = l3, edge = "z", pos = 0.5, col = labelCol, at = 1.1, cex = 1.5)

# Print LL's of exact and heuristic solutions
cat(paste("Exact max LL = ", constrainedMaxLL, "\n"))
if (plotPathToSolution) {
  cat(paste("Heur. max LL = ", bestLL, "\n"))
}

# Save the plot to a PNG file if wanted
# rgl.snapshot(filename = "plot3.png")

