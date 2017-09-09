# Read in simulated pharmacokinetic data for a specific study
# Simulate pharmacodynamics
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define output directory and read simulated pharmacokinetic data
# Output is saved into a totally different folder from scripts because my code
# is tracked by git and stored online.  Keeping all output (i.e., .csv and .png)
# files in the git repository slows down synchronisation
# PD output is saved in the same folder as PK output
  project.name <- "ExposureResponse"
  output.dir <- paste0("/Volumes/Prosecutor/sunitinib_nogit/",project.name)
  study.name <- "nodrug"
  study.dir <- paste0(output.dir,"/",study.name)
  setwd(study.dir)
  pk.data <- read.csv(file = paste0(study.name,"_pk_data.csv"))	# Read in PK data

# ------------------------------------------------------------------------------
# Source and run files required prior to PD simulation
# Define simulation file directory for the project
  file.dir <- "/Volumes/Prosecutor/sunitinib/Project/SimulationFiles/"
  project.dir <- paste0(file.dir,project.name)
  setwd(project.dir)
# Functions file (package libraries, summary functions, plotting)
  source("functions.R")
# Population pharmacodynamic model of biomarkers, adverse effects and overall
# survival due to sunitinib exposure
  source("sunitinib_pd_model.R")
# Define PK and PD simulation times
  source("times.R")

# ------------------------------------------------------------------------------
# PD simulation
  setwd(study.dir)
# Pull population characteristics from pk.data
  nid <- length(unique(pk.data$ID))	# Number of unique individuals
  nsim <- length(unique(pk.data$SIM))	# Number of times individuals with x
  # characteristics will be simulated
  ntotal <- nid*nsim	# Total number of individuals to be simulated
  SIM.seq <- rep(1:nsim,times = nid) %>% sort	# Sequence of simulation numbers
  ID.seq <- rep(1:nid,times = nsim)	# ID sequence
# Simulate pharmacodynamic random effect parameters
# Covariance terms, therefore need to use multivariate random number generator
  if (ntotal > 1) {
    pd.ETA.matrix <- mvrnorm(ntotal,
      mu = rep(0,times = dim(pd.OMEGA)[1]),pd.OMEGA) %>%
      as.data.frame
    names(pd.ETA.matrix) <- c("ETAVEGFR3BASE","ETAVEGFR3MRT","ETAVEGFR3I50",
      "ETASKITBASE","ETASKITMRT","ETASKITI50","ETASKITSLP","ETAKG","ETAKRSKIT",
      "ETAKRD","ETAOBASE","ETAANCBASE","ETAANCMTT","ETAANCEMAX","ETAANCE50",
      "ETABPBASE","ETABPSLP","ETABPMRT","ETAHFS0","ETAHFS1","ETAHFS2","ETAFAT0",
      "ETAFAT1","ETAFAT2","ETAFAT3")
  } else {
    pd.ETA.matrix <- data.frame(ETAVEGFR3BASE = 0,ETAVEGFR3MRT = 0,
      ETASKITBASE = 0,ETASKITMRT = 0,ETASKITSLP = 0,ETAVEGFR3I50 = 0,
      ETASKITI50 = 0,ETAKG = 0,ETAKRSKIT = 0,ETAKRD = 0,ETAOBASE = 0,
      ETAANCBASE = 0,ETAANCMTT = 0,ETAANCEMAX = 0,ETAANCE50 = 0,ETABPBASE = 0,
      ETABPSLP = 0,ETABPMRT = 0,ETAHFS0 = 0,ETAHFS1 = 0,ETAHFS2 = 0,ETAFAT0 = 0,
      ETAFAT1 = 0,ETAFAT2 = 0,ETAFAT3 = 0)
  }
# Input data frame for simulation
  pkpd.data <- pk.data[c("SIM","ID","time","cyc","amt","WT","OBASE","HFSBASE",
    "FATBASE","IPREP","IPREM","IPRE","AUC24")]
  pkpd.data <- pkpd.data[pkpd.data$time %in% pd.times,]
  pd.ID.data <- data.frame(SIM = SIM.seq,ID = ID.seq,pd.ETA.matrix)
  input.pd.data <- merge(pkpd.data,pd.ID.data,by = c("SIM","ID"),all = T)
  input.pd.data$cmt <- 0
  input.pd.data$evid <- 0
  input.pd.data <- input.pd.data[with(input.pd.data,order(input.pd.data$SIM,
    input.pd.data$ID,input.pd.data$time)),]
# For individual, replicate the cycle 2 portion of pk.data to extend for the
# full PD duration
  cycle1.data <- input.pd.data[input.pd.data$cyc == 1,]
  cycle2.data <- input.pd.data[input.pd.data$cyc == 2 &
    input.pd.data$time < max(pk.times),]
  end.cycle2.data <- input.pd.data[input.pd.data$time == max(pk.times),]
  rep.factor <- pd.ncycles-2
  extend.data <- lapply(cycle2.data,rep.int,times = rep.factor) %>%
    as.data.frame
  input.pd.data <- rbind(cycle1.data,cycle2.data,extend.data,end.cycle2.data)
  input.pd.data <- input.pd.data[with(input.pd.data,order(input.pd.data$SIM,
    input.pd.data$ID)),]
  input.pd.data$time <- pd.times
  input.pd.data$cyc <- pd.cyc
# Generate random numbers of each individual at each time-point to calculate
# status of alive status and dropout status
  input.pd.data$RSURV <- runif(ntotal*length(pd.times),min = 0,max = 1)
  input.pd.data$RDROP <- runif(ntotal*length(pd.times),min = 0,max = 1)
# Simulate
  pd.data <- pd.mod %>% mrgsim(data = input.pd.data,
    carry.out = c("SIM","cyc","amt","IPREP","IPREM","IPRE","AUC24")) %>%
    as.data.frame
# Use hand-foot syndrome state probabilities to simulate HFS profile
  pd.data <- ddply(pd.data, .(SIM,ID), simulate.HFS.grade, .progress = "text")
# Use fatigue state probabilities to simulate FAT profile
  pd.data <- ddply(pd.data, .(SIM,ID), simulate.FAT.grade, .progress = "text")
# Determine overall survival for the population over the study period
  pd.data <- ddply(pd.data, .(SIM,ID), alive.function, .progress = "text")
# Clean up and save the PD simulated data
  output.pd.data <- pd.data[c("SIM","ID","time","cyc","amt","IPREP",
    "IPREM","IPRE","AUC24","TUMOUR","ANC","BP","WT","OBASE","HFSBASE","FATBASE",
    "IPRE_VEGFR3","IPRE_SKIT","HFS","FAT","status")]
  write.csv(output.pd.data,file = paste0(study.name,"_pd_data.csv"),
    quote = FALSE,row.names = FALSE)

# ------------------------------------------------------------------------------
# Plot results
# Overall survival
  # For overall survival plots, calculate the proportion of individuals alive
  # at each time-point
  # Kaplan Meier Plot for survival confidence intervals
  # All individuals have the same start time, i.e., time == 0
    pd.data$start <- 0
    pd.data <- ddply(pd.data, .(SIM,ID), stop.time.function) # For each
    # individual calculate their stop time
    km.data <- ddply(pd.data, .(SIM,ID), headperID)	# First line per individual
    km.data <- km.data[c("SIM","ID","start","stop","event")]
    S <- Surv(time = km.data$start,time2 = km.data$stop,event = km.data$event)
    result <- survfit(formula = S ~ 1,data = km.data)
    cols <- lapply(2:12, function(x) summary(result)[x])
    surv.data <- do.call(data.frame, cols)
  # Plot overall survival
    plotobj1 <- NULL
    plotobj1 <- ggplot()
    plotobj1 <- plotobj1 + geom_line(aes(x = time/24/7,y = surv),
      colour = "brown3",data = surv.data)
    plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper),data = surv.data,fill = "brown3",alpha = 0.3)
    plotobj1 <- plotobj1 + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj1 <- plotobj1 + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(surv.data$time)/24/7,by = 12),
      lim = c(0,max(surv.data$time/24/7)))
    # print(plotobj1)
    ggsave(plot = plotobj1,
      filename = paste0(study.name,"_overallsurvival.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
# Summarise and plot sVEGFR-3 over time
  summary.VEGFR3 <- ddply(pd.data, .(time),
    function(pd.data) graded.summary(pd.data$IPRE_VEGFR3))
# Plot median and confidence intervals and facet
  plotobj2 <- NULL
  plotobj2 <- ggplot(summary.VEGFR3[summary.VEGFR3$time <= 50*7*24,])
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
    ymax = CI90hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
    ymax = CI80hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
    ymax = CI60hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
    ymax = CI40hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
    ymax = CI20hi),fill = "skyblue4",alpha = 0.2)
  plotobj2 <- plotobj2 + geom_line(aes(x = time/24/7,y = med),
    colour = "skyblue4")
  plotobj2 <- plotobj2 + scale_y_continuous("sVEGFR-3 Concentration (pg/mL)",
    lim = c(8,12))
  plotobj2 <- plotobj2 + scale_x_continuous("Time (weeks)")
  # print(plotobj2)
  ggsave(plot = plotobj2,filename = paste0(study.name,"_sVEGFR3.png"),
    width = 30,height = 15,unit = "cm",dpi = 300)
# Summarise and plot sKIT over time
  summary.SKIT <- ddply(pd.data, .(time),
    function(pd.data) graded.summary(pd.data$IPRE_SKIT))
# Plot median and confidence intervals and facet
  plotobj3 <- NULL
  plotobj3 <- ggplot(summary.SKIT[summary.SKIT$time <= 50*7*24,])
  plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
    ymax = CI90hi),fill = "skyblue4",alpha = 0.2)
  plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
    ymax = CI80hi),fill = "skyblue4",alpha = 0.2)
  plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
    ymax = CI60hi),fill = "skyblue4",alpha = 0.2)
  plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
    ymax = CI40hi),fill = "skyblue4",alpha = 0.2)
  plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
    ymax = CI20hi),fill = "skyblue4",alpha = 0.2)
  plotobj3 <- plotobj3 + geom_line(aes(x = time/24/7,y = med),
    colour = "skyblue4")
  plotobj3 <- plotobj3 + scale_y_continuous("sKIT Concentration (pg/mL)",
    lim = c(8,14))
  plotobj3 <- plotobj3 + scale_x_continuous("Time (weeks)")
  # print(plotobj3)
  ggsave(plot = plotobj3,filename = paste0(study.name,"_sKIT.png"),
    width = 30,height = 15,unit = "cm",dpi = 300)
# Summarise and plot tumour size over time
  summary.TUMOUR <- ddply(pd.data, .(time),
    function(pd.data) graded.summary(pd.data$TUMOUR))
# Plot median and confidence intervals and facet
  plotobj4 <- NULL
  plotobj4 <- ggplot(summary.TUMOUR[summary.TUMOUR$time <= 50*7*24,])
  plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
    ymax = CI90hi),fill = "skyblue4",alpha = 0.2)
  plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
    ymax = CI80hi),fill = "skyblue4",alpha = 0.2)
  plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
    ymax = CI60hi),fill = "skyblue4",alpha = 0.2)
  plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
    ymax = CI40hi),fill = "skyblue4",alpha = 0.2)
  plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
    ymax = CI20hi),fill = "skyblue4",alpha = 0.2)
  plotobj4 <- plotobj4 + geom_line(aes(x = time/24/7,y = med),
    colour = "skyblue4")
  plotobj4 <- plotobj4 + scale_y_continuous("Sum of Longest Diameters (mm)")
  plotobj4 <- plotobj4 + scale_x_continuous("Time (weeks)")
  # print(plotobj4)
  ggsave(plot = plotobj4,filename = paste0(study.name,"_tumour.png"),
    width = 30,height = 15,unit = "cm",dpi = 300)
# Summarise and plot absolute neutrophil count over time
  summary.ANC <- ddply(pd.data, .(time),
    function(pd.data) graded.summary(pd.data$ANC))
# Plot median and confidence intervals and facet
  plotobj5 <- NULL
  plotobj5 <- ggplot(summary.ANC[summary.ANC$time <= 50*7*24,])
  plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
    ymax = CI90hi),fill = "skyblue4",alpha = 0.2)
  plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
    ymax = CI80hi),fill = "skyblue4",alpha = 0.2)
  plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
    ymax = CI60hi),fill = "skyblue4",alpha = 0.2)
  plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
    ymax = CI40hi),fill = "skyblue4",alpha = 0.2)
  plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
    ymax = CI20hi),fill = "skyblue4",alpha = 0.2)
  plotobj5 <- plotobj5 + geom_line(aes(x = time/24/7,y = med),
    colour = "skyblue4")
  plotobj5 <- plotobj5 + scale_y_log10("Absolute Neutrophil Count (x10^9)")
  plotobj5 <- plotobj5 + scale_x_continuous("Time (weeks)")
  # print(plotobj5)
  ggsave(plot = plotobj5,filename = paste0(study.name,"_anc.png"),
    width = 30,height = 15,unit = "cm",dpi = 300)
# Summarise and plot diastolic blood pressure over time
  summary.BP <- ddply(pd.data, .(time),
    function(pd.data) graded.summary(pd.data$BP))
# Plot median and confidence intervals and facet
  plotobj6 <- NULL
  plotobj6 <- ggplot(summary.BP[summary.BP$time <= 50*7*24,])
  plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
    ymax = CI90hi),fill = "skyblue4",alpha = 0.2)
  plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
    ymax = CI80hi),fill = "skyblue4",alpha = 0.2)
  plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
    ymax = CI60hi),fill = "skyblue4",alpha = 0.2)
  plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
    ymax = CI40hi),fill = "skyblue4",alpha = 0.2)
  plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
    ymax = CI20hi),fill = "skyblue4",alpha = 0.2)
  plotobj6 <- plotobj6 + geom_line(aes(x = time/24/7,y = med),
    colour = "skyblue4")
  plotobj6 <- plotobj6 + scale_y_continuous("Diastolic Blood Pressure (mmHg)")
  plotobj6 <- plotobj6 + scale_x_continuous("Time (weeks)")
  # print(plotobj6)
  ggsave(plot = plotobj6,filename = paste0(study.name,"_bp.png"),
    width = 30,height = 15,unit = "cm",dpi = 300)
# Summarise and plot probability of hand-foot syndrome grade over time
  summary.HFS <- ddply(pd.data, .(time,HFS),
    function(pd.data) summary.count.function(pd.data$HFS))
  summary.HFS$Grade <- as.factor(summary.HFS$HFS)
  plotobj7 <- NULL
  plotobj7 <- ggplot(summary.HFS)
  plotobj7 <- plotobj7 + geom_line(aes(x = time/24/7,y = pro,colour = Grade))
  plotobj7 <- plotobj7 + scale_y_continuous("Probability of Hand-Foot Syndrome Grade",
    lim = c(0,1),
    breaks = seq(from = 0,to = 1,by = 0.2),
    labels = seq(from = 0,to = 1,by = 0.2))
  plotobj7 <- plotobj7 + scale_x_continuous("Time (weeks)")
  # print(plotobj7)
  ggsave(plot = plotobj7,filename = paste0(study.name,"_HFS.png"),
    width = 30,height = 15,unit = "cm",dpi = 300)
# Summarise and plot probability of fatigue grade over time
  summary.FAT <- ddply(pd.data, .(time,FAT),
    function(pd.data) summary.count.function(pd.data$FAT))
  summary.FAT$Grade <- as.factor(summary.FAT$FAT)
  plotobj8 <- NULL
  plotobj8 <- ggplot(summary.FAT)
  plotobj8 <- plotobj8 + geom_line(aes(x = time/24/7,y = pro,colour = Grade))
  plotobj8 <- plotobj8 + scale_y_continuous("Probability of Fatigue Grade",
    lim = c(0,1),
    breaks = seq(from = 0,to = 1,by = 0.2),
    labels = seq(from = 0,to = 1,by = 0.2))
  plotobj8 <- plotobj8 + scale_x_continuous("Time (weeks)")
  # print(plotobj8)
  ggsave(plot = plotobj8,filename = paste0(study.name,"_FAT.png"),
    width = 30,height = 15,unit = "cm",dpi = 300)
