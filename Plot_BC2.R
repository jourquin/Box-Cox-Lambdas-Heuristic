# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# R script that plots the Log-Likelihood of the logits computed for all the possible combinations of lambdas (bivariate case).
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

# Group of commodities to plot
group <- 5

# Plot unconstrained max LL (red dot)
plotUnconstrainedMax <- TRUE

# Plot all the solutions tested by the heuristic. The best one will be in blue if different from the exact bes solution
plotPathToSolution <- TRUE

# input dataset
modelName <- "europeL2-01"

######################################################################
# Nothing should be changed beyond this line                     #####
######################################################################

# Size of the plot area and position of the title
windowSize <- 800
titleLine <- -50

# Number of variables in the utility function of the logit model
nbVariables <- 2

# Get the data for the current group from the brute force results
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
load(paste(modelName, "-BruteForce2.Rda", sep = ""))
data <- bfSolutions[bfSolutions$group == group, ]

# Specific settings (radius, zoom and emphase levels for each group)
source("_PlotSettings.R")

# Get unconstrained max Log Likelihood
unconstrainedMaxLL <- max(data$LL, na.rm = TRUE)

# Get all the valid solutions, with or without the best unconstrained one
if (plotUnconstrainedMax) {
  LambdaCost <- data$lambda.cost[data$error == "" | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
  LambdaDuration <- data$lambda.duration[data$error == "" | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
  LogLik <- data$LL[data$error == "" | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
  best <- data$best[data$error == "" | (data$LL == unconstrainedMaxLL & !is.na(data$LL))]
} else {
  LambdaCost <- data$lambda.cost[data$error == ""]
  LambdaDuration <- data$lambda.duration[data$error == ""]
  LogLik <- data$LL[data$error == ""]
  best <- data$best[data$error == ""]
}

# Default color and size of the dots
colors <- rep("gray95", length(LogLik))
sizes <- (LogLik - min(LogLik)) / min(LogLik)
emphased <- rep(FALSE, length(LogLik))

# Force to use the -2, 2 range in all dimensions using invisible spheres
LambdaCost <- c(LambdaCost, -2)
LambdaDuration <- c(LambdaDuration, -2)
colors <- c(colors, "white")
LogLik <- c(LogLik, min(LogLik))
sizes <- c(sizes, 0)

LambdaCost <- c(LambdaCost, 2)
LambdaDuration <- c(LambdaDuration, 2)
colors <- c(colors, "white")
LogLik <- c(LogLik, min(LogLik))
sizes <- c(sizes, 0)

# Find best LL under constraints
constrainedMaxIdx <- which(best == TRUE)
constrainedMaxLL <- LogLik[constrainedMaxIdx]

# Identify path to best LL if wanted
if (plotPathToSolution) {
  # Load the solutions tested by the heuristic for the group
  inputFile <- paste(modelName, "-Heuristic2.Rda", sep = "")
  load(inputFile)
  solutions <- heuristicPath[heuristicPath$group == group, ]

  lambdaCost2 <- solutions$lambda.cost
  lambdaDuration2 <- solutions$lambda.duration
  LogLik2 <- solutions$LL

  best2 <- solutions$best
  bestLL <- LogLik2[which(best2 == TRUE)]
  bestLLIdx <- which(LogLik == bestLL)

  # Mark all the tested solutions
  for (i in 1:length(LogLik2)) {
    idx <- which(LogLik == LogLik2[i])
    if (length(idx) == 0) {
      # When a non valid solutions is tested, add it to the solutions to plot
      LambdaCost <- c(LambdaCost, lambdaCost2[i])
      LambdaDuration <- c(LambdaDuration, lambdaDuration2[i])
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


# Identify the unconstrained max using a larger red dot
if (plotUnconstrainedMax == TRUE) {
  unconstrainedMaxIdx <- which(LogLik == unconstrainedMaxLL)
  colors[unconstrainedMaxIdx] <- "red"
  if (!emphased[unconstrainedMaxIdx]) {
    sizes[unconstrainedMaxIdx] <- emphaseLevel * sizes[unconstrainedMaxIdx]
    emphased[unconstrainedMaxIdx] <- TRUE
  }
}

# Identify the best solution under constraints using a larger green dot
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
l3 <- "Log Likelihood"

# Hide the default axis of the scatter3d plot
axisCol <- rgb(255, 255, 255, max = 255, alpha = 0, names = "transparentWhite")

# Plot and adapt perspective to settings
scatter3d(x = LambdaCost, y = LambdaDuration, z = LogLik, fogtype = "exp2", surface = FALSE, sphere.size = sphereSize, radius = sizes, point.col = colors, axis.col = c(axisCol, axisCol, axisCol), xlab = l1, ylab = l2, zlab = l3)
rgl.viewpoint(theta = 180, phi = 90, fov = 45, zoom = zoomLevel, interactive = TRUE)
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
rgl::axis3d("x", pos = c(0, 0.5, 0), lwd = 2, col = axisCol, labels = c("-2", "0", "2"), tick = FALSE, nticks = 2, cex = 1.2)
rgl::axis3d("y", pos = c(0.5, 0, 0), lwd = 2, col = axisCol, labels = c("-2", "2"), tick = FALSE, nticks = 1, cex = 1.2)
rgl::axis3d("z", pos = 0.5, lwd = 2, col = axisCol, tick = FALSE, labels = "", cex = 1.2)
rgl::mtext3d(text = l1, edge = "x", pos = c(0, 0.5, 0), col = labelCol, at = 1.0, cex = 1.5)
rgl::mtext3d(text = l2, edge = "y", pos = c(0.5, 0, 0), col = labelCol, at = 1.0, cex = 1.5)
rgl::mtext3d(text = l3, edge = "z", pos = 0.5, col = labelCol, at = 1.1, cex = 1.2)

# Print LL's of exact and heuristic solutions
cat(paste("Exact max LL = ", constrainedMaxLL, "\n"))
if (plotPathToSolution) {
  cat(paste("Heur. max LL = ", bestLL, "\n"))
}

# Save the plot to a PNG file if wanted
# rgl.snapshot(filename = "plot2.png")
