# Universal functions for sunitinib project
# ------------------------------------------------------------------------------
# Package libraries
  library(MASS)	# mvrnorm function
  library(MBESS)	# cor2cov function
  library(ggplot2)  # Plotting
  library(grid)   # Plotting
  library(GGally)	# For plotting plotmatrix
  library(plyr)  # Split and rearrange data, ddply function
  library(dplyr)  # New plyr
  library(mrgsolve) # Metrum differential equation solver for pharmacometrics
  library(reshape)	# melt function
  library(survival)	# Kaplan-Meier plots
  library(stringr)	# Split character strings

# ------------------------------------------------------------------------------
# Statistical functions
# Simulation seed for reproducible numbers
  set.seed(123456)
# 95% prediction interval functions
  CI95lo <- function(x) quantile(x,probs = 0.025,na.rm = TRUE)
  CI95hi <- function(x) quantile(x,probs = 0.975,na.rm = TRUE)
# 90% prediction interval functions
  CI90lo <- function(x) quantile(x,probs = 0.05,na.rm = TRUE)
  CI90hi <- function(x) quantile(x,probs = 0.95,na.rm = TRUE)
# Summary function (median and confidence intervals)
  summary.function <- function(x) {
    n <- length(x)
    lo <- CI95lo(x)
    hi <- CI95hi(x)
    med <- median(x)
    result <- c(n,med,lo,hi)
    names(result) <- c("n","med","lo95","hi95")
    result
  }
# Graded summary function
  graded.summary <- function(x) {
    n <- length(x)
    CI90lo <- quantile(x,probs = 0.05,na.rm = TRUE)
    CI80lo <- quantile(x,probs = 0.1,na.rm = TRUE)
    CI60lo <- quantile(x,probs = 0.2,na.rm = TRUE)
    CI40lo <- quantile(x,probs = 0.3,na.rm = TRUE)
    CI20lo <- quantile(x,probs = 0.4,na.rm = TRUE)
    med <- median(x)
    CI20hi <- quantile(x,probs = 0.6,na.rm = TRUE)
    CI40hi <- quantile(x,probs = 0.7,na.rm = TRUE)
    CI60hi <- quantile(x,probs = 0.8,na.rm = TRUE)
    CI80hi <- quantile(x,probs = 0.9,na.rm = TRUE)
    CI90hi <- quantile(x,probs = 0.95,na.rm = TRUE)
    result <- c(n,med,CI90lo,CI90hi,CI80lo,CI80hi,CI60lo,CI60hi,CI40lo,CI40hi,
      CI20lo,CI20hi)
    names(result) <- c("n","med","CI90lo","CI90hi","CI80lo","CI80hi","CI60lo",
      "CI60hi","CI40lo","CI40hi","CI20lo","CI20hi")
    result
  }
# Summary count function
  summary.count.function <- function(x) {
    total.n <- ntotal
    n <- length(x)
    result <- n/total.n
    names(result) <- "pro"
    result
  }
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

# ------------------------------------------------------------------------------
# Pharmacokinetic functions not written in mrgsolve code
# Calculate the AUC in the previous 24-hours
  auc24.function <- function(df) {
    x <- df$AUC
    y <- x
    first.row <- 24/pk.increment+1
    for (i in first.row:length(x)) {
      y[i] <- x[i]-x[i-(first.row-1)]
    }
    df$AUC24 <- y
    df
  }

# ------------------------------------------------------------------------------
# Pharmacodynamic functions not written in mrgsovle code
# Account for drop out during tumour simulations
  tumour.dropout <- function(df) {
  # Drop out model parameters
    TINT <- -3.49 # Intercept of drop out
    TDP <- 1.12 # Parameter related to occurrence of disease progression (a 20% increase in tumour size since nadir, yes/no)
    TSLD <- 0.00105 # Parameter related to tumour size at drop out (sum of longest diamaters, mm)
    TTIME <- 0.00707/7/24 # Parameter related to time since start of study, hours^-1
  # Identify the nadir of the tumour size
    tummin <- min(df$TUMOUR)	# Individual's minimum tumour size
    t.tummin <- df$time[df$TUMOUR == tummin]	# Time of minimum tumour size
    df$DP <- 0	# Initially no disease progression is considered until nadir is observed
    for (i in 1:nrow(df)) {
      if (df$time[i] < t.tummin) {
        df$DP[i] <- 0
      }
      if (df$time[i] >= t.tummin) {
        curr.tum <- df$TUMOUR[i]
        if (curr.tum > tummin*1.2) {
          dp <- 1
        } else {
          dp <- 0
        }
        df$DP[i] <- dp
      }
      tumdrop <- TINT+TDP*df$DP[i]+TSLD*df$TUMOUR[i]+TTIME*df$time[i]
      Plogitdrop <- exp(tumdrop)/(1+exp(tumdrop))
      if (i == 1) {
        df$Pdrop[i] <- 0
      } else {
        df$Pdrop[i] <- 1-(1-Plogitdrop)^((df$time[i]-df$time[i-1])/7/24)
      }
    }
    df$Tdrop <- 0
    for (i in 2:nrow(df)) {
      prev.Tdrop <- df$Tdrop[i-1]
      if (prev.Tdrop == 1) {
        Tdrop <- 1
      } else {
        Tdrop <- sample(c(1,0),size = 1,prob = c(df$Pdrop[i],1-df$Pdrop[i]))
      }
      df$Tdrop[i] <- Tdrop
    }
    df
  }

# Simulate fatigue grade based on probabilities
  simulate.FAT.grade <- function(df) {
    df$FAT <- NA
    scales <- c(0,1,2,3)
    for (i in 1:nrow(df)) {
      if (i == 1) {
        FAT <- df$FATBASE[i]
      } else {
        prev.FAT <- df$FAT[i-1]
        if (prev.FAT == 0) {
          fat.probs <- c(df$FATPROB00[i],
            df$FATPROB01[i],df$FATPROB02[i],df$FATPROB03[i])
        }
        if (prev.FAT == 1) {
          fat.probs <- c(df$FATPROB10[i],
            df$FATPROB11[i],df$FATPROB12[i],df$FATPROB13[i])
        }
        if (prev.FAT == 2) {
          fat.probs <- c(df$FATPROB20[i],
            df$FATPROB21[i],df$FATPROB22[i],df$FATPROB23[i])
        }
        if (prev.FAT == 3) {
          fat.probs <- c(df$FATPROB30[i],
            df$FATPROB31[i],df$FATPROB32[i],df$FATPROB33[i])
        }
        FAT <- sample(scales,size = 1,prob = fat.probs)
      }
      df$FAT[i] <- FAT
    }
    df
  }

# Simulate hand-foot syndrome grade based on probabilities
  simulate.HFS.grade <- function(df) {
    df$HFS <- NA
    scales <- c(0,1,2,3)
    for (i in 1:nrow(df)) {
      if (i == 1) {
        HFS <- df$HFSBASE[i]
      } else {
        prev.HFS <- df$HFS[i-1]
        if (prev.HFS == 0) {
          hfs.probs <- c(df$HFSPROB00[i],
            df$HFSPROB01[i],df$HFSPROB02[i],df$HFSPROB03[i])
        }
        if (prev.HFS == 1) {
          hfs.probs <- c(df$HFSPROB10[i],
            df$HFSPROB11[i],df$HFSPROB12[i],df$HFSPROB13[i])
        }
        if (prev.HFS == 2) {
          hfs.probs <- c(df$HFSPROB20[i],
            df$HFSPROB21[i],df$HFSPROB22[i],df$HFSPROB23[i])
        }
        if (prev.HFS == 3) {
          hfs.probs <- c(df$HFSPROB30[i],
            df$HFSPROB31[i],df$HFSPROB32[i],df$HFSPROB33[i])
        }
        HFS <- sample(scales,size = 1,prob = hfs.probs)
      }
      df$HFS[i] <- HFS
    }
    df
  }

# ------------------------------------------------------------------------------
# Functions for overall survival
# Determine whether individual is dead/alive or dropped out
  alive.function <- function(df) {
    df$status <- 1
    for (i in 2:nrow(df)) {
      prev.status <- df$status[i-1]
      if (prev.status == 0) {
        df$status[i] <- 0
      } else {
        Ssurv <- (df$SURV[i-1]-df$SURV[i])/df$SURV[i-1]
        Sdrop <- (df$DROP[i-1]-df$DROP[i])/df$DROP[i-1]
        surv.status <- sample(c(0,1),size = 1,prob = c(Ssurv,1-Ssurv))
        drop.status <- sample(c(0,1),size = 1,prob = c(Sdrop,1-Sdrop))
        if (drop.status == 0) {
          df$status[i] <- drop.status
        } else {
          df$status[i] <- surv.status
        }
      }
    }
    df
  }

# Calculate proportion of individuals alive at each time-point
  pro.alive.function <- function(df) {
    n <- length(df$status)	# Total number of individuals
    n.alive <- length(df$status[df$status == 1])
    pro.alive <- n.alive/n
    result <- c("n" = n,"n.alive" = n.alive,"pro.alive" = pro.alive)
    return(result)
  }

# For each individual calculate their stop time (death, study end, censored)
# Required for Kaplan-Meier plot function
  stop.time.function <- function(df) {
    if (tail(df$status,1) == 1) {
      stop.time <- max(pd.times)
      event <- 0
    } else {
      stop.time <- head(df$time[df$status == 0],1)
      event <- 1
    }
    df$stop <- stop.time
    df$event <- event
    df
  }

# ------------------------------------------------------------------------------
# Plotting functions
# Define a custom ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
