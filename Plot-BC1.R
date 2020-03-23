# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# R script that plots the Log-Likelihood of the logits computed for all the possible lambdas (univariate case).
# The results are taken from the output of the brute force algorithm.
# The best solution is identified by a green dot. Only the "valid" solutions are plot.
#
# The solutions tested by the heuristic can also be plotted, using the output of the "HeuristicVsBruteForce.R" script. The final
# solution of the heuristic is plot in blue if it is different from the (exact) best solution identified by the brute force
# approach.
#
# While being a little bit tricky, a 3D plot is used in order to obtain similar visual outputs as for the
# bivariate and trivariate cases.
#
# Predefined plot settings are provided for each group. These are specific for the provided dataset in this project.
#

library(rgl)
library(car)

# Clear memory if wanted
rm(list = ls(all.names = TRUE))

# Group to plot
group = 7

# Plot unconstrained max LL (red dot, but not used in the univariate plot as the heuristic always finds the best solution)
plotUnconstrainedMax = FALSE

# Plot all the solutions tested by the heuristic. The best one will be in blue if different from the exact bes solution
plotPathToSolution = FALSE

# input dataset
modelName = "europeL2-01"

######################################################################
# Nothing should be changed beyond this line                     #####
######################################################################

# Size of the plot area and position of the title
windowSize = 800
titleLine = -47

# Number of variables in the utility function of the logit model
nbVariables = 1

# Get the data for the current group from the brute force results
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
load(paste(modelName,"-BruteForce1.Rda", sep=""))
data = bfSolutions[bfSolutions$group == group, ]

# Specific settings (radius, zoom and emphase levels for each group)
source("_PlotSettings.R")

# Get the lambda and LL of the best valid solution 
LambdaCost = data$lambda.cost[data$error == ""]
LogLik = data$LL[data$error == ""]

# Set default color and sizes of th dots
colors <- rep("gray95", length(LogLik))
sizes =  (LogLik - min(LogLik)) /min(LogLik)
emphased <-rep(FALSE, length(LogLik))



# Identify path to best LL if wanted
if (plotPathToSolution) {
  # Load the solutions tested by the heuristic for the group
  inputFile = paste(modelName, "-Heuristic1.Rda", sep="")
  load(inputFile)
  solutions = heuristicPath[heuristicPath$group == group, ]
  
  lambdaCost2 = solutions$lambda.cost[solutions$isSolved == TRUE]
  LogLik2 = solutions$LL[solutions$isSolved == TRUE]
  
  best2 = solutions$best[solutions$isSolved == TRUE]
  bestLL = LogLik2[which(best2 == TRUE)]
  bestLLIdx = which(LogLik == bestLL)
  
  # Mark all the tested solutions
  for (i in 1:length(LogLik2)) {
    idx = which(LogLik == LogLik2[i]) 
    if (length(idx) == 0) {
      # When a non valid solutions is tested, add it to the solutions to plot
      LambdaCost = c(LambdaCost, lambdaCost2[i])
      LogLik = c(LogLik, LogLik2[i])
      colors = c(colors, "black")
      emphased = c(emphased, FALSE)
    } else  {
      colors[idx] = "black"
    }
  }
  
  # Recomput sizes of nodes
  LogLik[is.na(LogLik)] = min(LogLik, na.rm=TRUE)
  sizes =  (LogLik - min(LogLik, na.rm=TRUE)) /min(LogLik, na.rm=TRUE)
  
  # Identify solution of heuristic
  colors[bestLLIdx] = "blue"
  if (!emphased[bestLLIdx]) {
    sizes[bestLLIdx] = emphaseLevel * sizes[bestLLIdx]
    emphased[bestLLIdx] = TRUE  
  }
} 

# Identify the best solution
bestLLIdx = which(LogLik == max(LogLik)) 
colors[bestLLIdx] = "green"
if (!emphased[bestLLIdx]) {
  sizes[bestLLIdx] = emphaseLevel*sizes[bestLLIdx] 
}


# Normalize the LL
minLL = min(LogLik)
for (i in 1:length(LogLik)) {
  LogLik[i] = LogLik[i] - minLL
}

# Force a 3D display and -2, 2 range with 2 dummy spheres that will not be displayed
zz = rep(-2, length(LogLik))

LambdaCost = c(LambdaCost, -2)
zz = c(zz, -0.1)
LogLik = c(LogLik, min(LogLik))
sizes = c(sizes, 0)

LambdaCost = c(LambdaCost, 2)
zz = c(zz, -0.1)
LogLik = c(LogLik, min(LogLik))
sizes = c(sizes, 0)

# Title and axis labels with greek lambda
l1 = parse(text=paste(intToUtf8(955), "[C]^", group, sep=""))
l2 = "Log Likelihood"
title = paste("NST-R", group)

# Display
open3d()
par3d(windowRect=c(0,0,windowSize,windowSize))

# Plot and adapt perspective to settings
tWhite <- rgb(255, 255, 255, max = 255, alpha = 0, names = "transparentWhite")
scatter3d(x=LambdaCost, y=LogLik, z=zz, surface = FALSE, fogtype = "none", sphere.size = sphereSize, radius = sizes, point.col = colors, axis.col =c("gray20",tWhite,tWhite), xlab="", zlab = "", ylab="")
rgl.viewpoint( theta = 0, phi = 0, fov = 0, zoom = zoomLevel, interactive = TRUE )

axisCol = "gray20"
labelCol = "blue"

# Add title and axis with labels
bgplot3d({
  plot.new()
  title(main = title, line = titleLine, cex.main=1.5)
})
rgl::axis3d("y", pos = c(0.5, 0.0, 0), lwd=2, col = axisCol, tick=FALSE, nticks=0, labels = c(""))
rgl::mtext3d(text=l1, edge="x", pos=c(0.5, 0, 0), col=labelCol, at=1.15, cex=1.5)
rgl::mtext3d(text=l2, edge="y", pos=c(0.5,0, 0), col=labelCol, at=1.1, cex=1.2)

# Save the plot to a PNG file if wanted
#rgl.snapshot(filename = "plot1.png")
#rgl.close()
