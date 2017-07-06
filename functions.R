# Universal functions for sunitinib project
# ------------------------------------------------------------------------------
# Package libraries
  library(MASS)	# mvrnorm function
  library(ggplot2)  # Plotting
  library(grid)   #Plotting
  library(plyr)  # Split and rearrange data, ddply function
  library(dplyr)  # New plyr
  library(mrgsolve) # Metrum differential equation solver for pharmacometrics
  library(reshape)	# melt function
  library(survival)	# Kaplan-Meier plots

# ------------------------------------------------------------------------------
# Statistical functions
# Simulation seed for reproducible numbers
  set.seed(123456)
# 95% prediction interval functions
  CI95lo <- function(x) quantile(x,probs = 0.025)
  CI95hi <- function(x) quantile(x,probs = 0.975)
# Interquartile range
  CI50lo <- function(x) quantile(x,probs = 0.25)
  CI50hi <- function(x) quantile(x,probs = 0.75)

# ------------------------------------------------------------------------------
# Data re-arrangement functions
# Take the first line for each individual
  headperID <- function(x) {
    y <- head(x,1)
    return(y)
  }

# Assign dead or alive
  alive.function <- function(x) {
    x$ALIVE <- NA
    for (i in 1:nrow(x)) {
      if (x$time[i] == 0) {
        x$ALIVE[i] <- 1
      } else {
        prev.ALIVE <- x$ALIVE[i-1]
        if (prev.ALIVE == 0) {
          x$ALIVE[i] <- 0
        } else {
          Sdrop <- (x$SURD[i-1]-x$SURD[i])/x$SURD[i-1]
          Sevent <- (x$SUR[i-1]-x$SUR[i])/x$SUR[i-1]
          if (Sdrop > x$RDROP[i]) {
            x$ALIVE[i] <- 0
          } else {
            if (Sevent < x$REVENT[i]) x$ALIVE[i] <- 1
            if (Sevent > x$REVENT[i]) x$ALIVE[i] <- 0
          }
        }
      }
    }
    x
  }

# Calculate proportion of individuals alive at each time-point
  pro.alive.function <- function(x) {
    n <- length(x$ALIVE)	# Total number of individuals
    n.alive <- length(x$ALIVE[x$ALIVE == 1])
    pro.alive <- n.alive/n
    result <- c("n" = n,"n.alive" = n.alive,"pro.alive" = pro.alive)
    return(result)
  }

# ------------------------------------------------------------------------------
# Plotting functions
# Define a custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
