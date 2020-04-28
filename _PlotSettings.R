# Copyright (c) 2020 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# Specific sphere sizes, emphase and zoom levels for the differents groups of commodities. These values are specific to the
# dataset provided in this project.
#

if (nbVariables == 1) {
  if (group == 0) {
    sphereSize <- 0.6
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 1) {
    sphereSize <- 0.6
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 2) {
    sphereSize <- 0.6
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 3) {
    sphereSize <- 0.6
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 4) {
    sphereSize <- 0.6
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 5) {
    sphereSize <- 0.6
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 6) {
    sphereSize <- 0.6
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 7) {
    sphereSize <- 0.5
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 8) {
    sphereSize <- 0.55
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }

  if (group == 9) {
    sphereSize <- 0.5
    emphaseLevel <- 1
    zoomLevel <- 0.75
  }
}

if (nbVariables == 2) {
  # rotate aroud Z axis - specific to each group
  if (group == 0) {
    theta <- -1.6
    sphereSize <- 2
    emphaseLevel <- 1.1
    zoomLevel <- 0.7
  }

  if (group == 1) {
    theta <- -0.8
    sphereSize <- 2
    emphaseLevel <- 1.1
    zoomLevel <- 0.75
  }

  if (group == 2) {
    theta <- -1.2
    sphereSize <- 1.3
    emphaseLevel <- 1.1
    zoomLevel <- 0.7
  }

  if (group == 3) {
    theta <- 1.9
    sphereSize <- 1.8
    emphaseLevel <- 1.2
    zoomLevel <- 0.75
  }

  if (group == 4) {
    theta <- 2.2
    sphereSize <- 1.8
    emphaseLevel <- 1.2
    zoomLevel <- 0.7
  }

  if (group == 5) {
    theta <- 1.6
    sphereSize <- 1.5
    emphaseLevel <- 1.2
    zoomLevel <- 0.75
  }

  if (group == 6) {
    theta <- 2.2
    sphereSize <- 2
    emphaseLevel <- 1.1
    zoomLevel <- 0.7
  }

  if (group == 7) {
    theta <- 2.0
    sphereSize <- 1.4
    sphereSize <- 1.5
    emphaseLevel <- 1.2
    zoomLevel <- 0.75
  }

  if (group == 8) {
    theta <- -0.8
    sphereSize <- 1.5
    sphereSize <- 1.6
    emphaseLevel <- 1.2
    zoomLevel <- 0.75
  }

  if (group == 9) {
    theta <- 1.9
    sphereSize <- 1.5
    emphaseLevel <- 1.2
    zoomLevel <- 0.75
  }
}

if (nbVariables == 3) {
  if (group == 0) {
    sphereSize <- 3.5
    theta <- -2.9
    theta <- 1.9
    zoomLevel <- 0.09
    emphaseLevel <- 1.5
  }

  if (group == 1) {
    sphereSize <- 3.8
    theta <- 2.0
    zoomLevel <- 0.085
    emphaseLevel <- 1.5
  }

  if (group == 2) {
    sphereSize <- 2.7
    theta <- -0.5
    zoomLevel <- 0.08
    emphaseLevel <- 1.2
  }

  if (group == 3) {
    sphereSize <- 3
    theta <- -2.8
    zoomLevel <- 0.08
    emphaseLevel <- 1.5
  }

  if (group == 4) {
    sphereSize <- 2.6
    theta <- 2.2
    zoomLevel <- 0.075
    emphaseLevel <- 1.2
  }

  if (group == 5) {
    sphereSize <- 3.8
    theta <- 0.9
    zoomLevel <- 0.14
    emphaseLevel <- 1.75
  }

  if (group == 6) {
    sphereSize <- 3.7
    theta <- 2.9
    theta <- 2.65
    zoomLevel <- 0.09
    emphaseLevel <- 2
  }

  if (group == 7) {
    sphereSize <- 3.5
    theta <- 2.5
    zoomLevel <- 0.085
    emphaseLevel <- 1.2
  }

  if (group == 8) {
    sphereSize <- 4
    theta <- 1.7
    zoomLevel <- 0.1
    emphaseLevel <- 1.6
  }

  if (group == 9) {
    sphereSize <- 3.3
    theta <- -1.5
    theta <- 2.8
    zoomLevel <- 0.11
    emphaseLevel <- 1.4
  }
}
