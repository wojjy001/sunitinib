# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J,
# Karlsson MO, Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT
# as Predictors of Tumor Dynamics and Overall Survival Following Sunitinib
# Treatment in GIST. CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9
# ------------------------------------------------------------------------------
# Source model code
  source("hansson_2013_biomarker_model.R")

# ------------------------------------------------------------------------------
# Generate andom effects
  bio.ETA.matrix <- mvrnorm(nsim,
    mu = rep(0,times = dim(bio.OMEGA)[1]),bio.OMEGA) %>%
    as.data.frame

# ------------------------------------------------------------------------------
# Input data frame for simulatione
  ID.data <- data.frame(
    ID = ID.seq,
    DOSE = 0,
    bio.ETA.matrix
  )
  names(ID.data)[c(3:dim(ID.data)[2])] <- c("ETAMRT23","ETABM0","ETADPSLO","ETABM02","ETABM03","ETABM0S","ETAMRTS",
  "ETADPSLOS","ETAIC50","ETAIC502","ETAIC503","ETAIC50S")
  input.data <- lapply(ID.data,rep.int,times = length(bio.times)) %>% as.data.frame
  input.data <- input.data[with(input.data,order(input.data$ID)),]
  input.data$time <- bio.times
  input.data$DOSE[input.data$time %in% on.times] <- DOSE
  input.data$CL <- PK.CL
  input.data$cmt <- 1

# ------------------------------------------------------------------------------
# Simulate
  bio.data <- bio.mod %>% mrgsim(data = input.data) %>% as.data.frame

# ------------------------------------------------------------------------------
# Rearrange data frame for plotting
  melt.bio.data <- melt(bio.data,id = c("ID","time","DOSE","CL"),
    measure = c("IPRE_VEGF","IPRE_VEGFR2","IPRE_VEGFR3","IPRE_SKIT"))

# Plot
  plotobj <- NULL
  plotobj <- ggplot(melt.bio.data)
  plotobj <- plotobj + stat_summary(aes(x = time/7/24,y = value),
    geom = "line",fun.y = median,colour = "red")
  plotobj <- plotobj + stat_summary(aes(x = time/7/24,y = value),
    geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
    fill = "red",alpha = 0.3)
  plotobj <- plotobj + scale_y_log10("Biomarker (pg/mL)\n",
    breaks = c(0.1,1,10,100),labels = c(0.1,1,10,100),lim = c(1,NA))
  plotobj <- plotobj + scale_x_continuous("\nTime (weeks)",
    breaks = seq(from = 0,to = 18,by = 2))
  plotobj <- plotobj + facet_wrap(~variable)
  plotobj
