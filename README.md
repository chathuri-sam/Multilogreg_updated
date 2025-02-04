# **Multilogreg Updated**

## Overview

This repository contains an updated version of the Multilogreg package, incorporating fixes and enhancements for improved functionality. The key updates include:

Fixed spatstat package dependencies to ensure compatibility with recent versions.
Updated the simulation code to support:
 - Multiple covariates (previously restricted to a single covariate).
 - Arbitrary spatial windows (previously limited to a unit square).

---

## Installation

To install the updated package, run the following command in R:

# Install dependencies
install.packages(c("devtools", "spatstat"))

# Install from GitHub
devtools::install_github("chathuri-sam/Multilogreg_updated")

## Acknowledgments

This package builds on the original [Multilogreg package](https://github.com/kristianhessellund/Multilogreg.git). We acknowledge the contributions of the original authors and maintainers.
