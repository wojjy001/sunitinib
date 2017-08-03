# Read in simulated study data and perform statistical summaries or recreate
# plots
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define global directory for project
  global.dir <- "/Volumes/Prosecutor/sunitinib/Project/"
# Source required files
  source(paste0(global.dir,"functions.R"))	# Universal functions file
# Set directory where simulation output files are saved
  sim.dir <- paste0(global.dir,"Output")
# Set working directory
  setwd(sim.dir)
# Read in simulation files in the output directory
# List the folders in the output directory
  folder.list <- list.dirs(sim.dir)
  folder.list <- folder.list[2:length(folder.list)]
  for (i in 1:length(folder.list)) {
    folder.list[i] <- str_split(folder.list[i],pattern = sim.dir)
    folder.list[i] <- folder.list[[i]][2]
  }
  folder.list <- unlist(folder.list)
  # folder.list <- "/bayes_tdm06"
# Read in .csv files from each folder
  read.pd.data <- function(x) {
    pd.data <- read.csv(file = paste0(sim.dir,x,"/",x,"_pd_data.csv"))
    pd.data$study <- x
    pd.data
  }
  pd.data <- ldply(folder.list,read.pd.data,.progress = "text")
# For biomarker and adverse effect simulations, only plot the first 50 weeks
  early.data <- pd.data[pd.data$time <= 50*24*7,]

# ------------------------------------------------------------------------------
# For overall survival plots, calculate the proportion of individuals alive
# at each time-point
  survival.data <- ddply(pd.data, .(study,time), pro.alive.function)

# ------------------------------------------------------------------------------
# Plot overall survival over time
  plot.os.function <- function(x) {
    plotobj1 <- NULL
    plotobj1 <- ggplot(x)
    plotobj1 <- plotobj1 + geom_step(aes(x = time/24/7,y = pro.alive),
      colour = "red")
    plotobj1 <- plotobj1 + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj1 <- plotobj1 + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(survival.data$time)/24/7,by = 24),
      lim = c(0,max(survival.data$time/24/7)))
    print(plotobj1)

    ggsave(plot = plotobj1,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_overallsurvival.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(survival.data, .(study), plot.os.function)

# ------------------------------------------------------------------------------
# Summarise and plot sVEGFR-3 over time
  summary.VEGFR3 <- ddply(early.data, .(study,time),
    function(early.data) graded.summary(early.data$IPRE_VEGFR3))
  plot.vegfr3.function <- function(x) {
    plotobj2 <- NULL
    plotobj2 <- ggplot(x)
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi),fill = "blue",alpha = 0.1)
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi),fill = "blue",alpha = 0.1)
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi),fill = "blue",alpha = 0.1)
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi),fill = "blue",alpha = 0.1)
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi),fill = "blue",alpha = 0.1)
    plotobj2 <- plotobj2 + geom_line(aes(x = time/24/7,y = med),colour = "blue")
    plotobj2 <- plotobj2 + scale_y_continuous("sVEGFR-3 Concentration (pg/mL)",
      lim = c(8,12))
    plotobj2 <- plotobj2 + scale_x_continuous("Time (weeks)")
    print(plotobj2)

    ggsave(plot = plotobj2,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_sVEGFR3.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(summary.VEGFR3, .(study), plot.vegfr3.function)

# Summarise and plot sKIT over time
  summary.SKIT <- ddply(early.data, .(study,time),
    function(early.data) graded.summary(early.data$IPRE_SKIT))
  plot.skit.function <- function(x) {
    plotobj3 <- NULL
    plotobj3 <- ggplot(x)
    plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi),fill = "red",alpha = 0.1)
    plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi),fill = "red",alpha = 0.1)
    plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi),fill = "red",alpha = 0.1)
    plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi),fill = "red",alpha = 0.1)
    plotobj3 <- plotobj3 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi),fill = "red",alpha = 0.1)
    plotobj3 <- plotobj3 + geom_line(aes(x = time/24/7,y = med),colour = "red")
    plotobj3 <- plotobj3 + scale_y_continuous("sKIT Concentration (pg/mL)",
      lim = c(8,14))
    plotobj3 <- plotobj3 + scale_x_continuous("Time (weeks)")
    print(plotobj3)

    ggsave(plot = plotobj3,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_sKIT.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(summary.SKIT, .(study), plot.skit.function)

# Sumamrise and plot tumour size over time
  summary.TUMOUR <- ddply(early.data, .(study,time),
    function(early.data) graded.summary(early.data$TUMOUR))
  plot.tumour.function <- function(x) {
    plotobj4 <- NULL
    plotobj4 <- ggplot(x)
    plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi),fill = "darkgreen",alpha = 0.1)
    plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi),fill = "darkgreen",alpha = 0.1)
    plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi),fill = "darkgreen",alpha = 0.1)
    plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi),fill = "darkgreen",alpha = 0.1)
    plotobj4 <- plotobj4 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi),fill = "darkgreen",alpha = 0.1)
    plotobj4 <- plotobj4 + geom_line(aes(x = time/24/7,y = med),colour = "darkgreen")
    plotobj4 <- plotobj4 + scale_y_continuous("Sum of Longest Diameters (mm)")
    plotobj4 <- plotobj4 + scale_x_continuous("Time (weeks)")
    print(plotobj4)

    ggsave(plot = plotobj4,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_tumour.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(summary.TUMOUR, .(study), plot.tumour.function)

# Summarise and plot absolute neutrophil count over time
  summary.ANC <- ddply(early.data, .(study,time),
    function(early.data) graded.summary(early.data$ANC))
  plot.anc.function <- function(x) {
    plotobj5 <- NULL
    plotobj5 <- ggplot(x)
    plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi),fill = "purple",alpha = 0.1)
    plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi),fill = "purple",alpha = 0.1)
    plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi),fill = "purple",alpha = 0.1)
    plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi),fill = "purple",alpha = 0.1)
    plotobj5 <- plotobj5 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi),fill = "purple",alpha = 0.1)
    plotobj5 <- plotobj5 + geom_line(aes(x = time/24/7,y = med),colour = "purple")
    plotobj5 <- plotobj5 + scale_y_log10("Absolute Neutrophil Count (x10^9 cells)",
      breaks = c(0.1,0.3,1,3,5,10),labels = c(0.1,0.3,1,3,5,10),
      lim = c(0.1,15))
    plotobj5 <- plotobj5 + scale_x_continuous("Time (weeks)")
    print(plotobj5)

    ggsave(plot = plotobj5,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_ANC.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(summary.ANC, .(study), plot.anc.function)

# Summarise and plot diastolic blood pressure over time
  summary.BP <- ddply(early.data, .(study,time),
    function(early.data) graded.summary(early.data$BP))
  plot.bp.function <- function(x) {
    plotobj6 <- NULL
    plotobj6 <- ggplot(x)
    plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI90lo,
      ymax = CI90hi),fill = "orange",alpha = 0.1)
    plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI80lo,
      ymax = CI80hi),fill = "orange",alpha = 0.1)
    plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI60lo,
      ymax = CI60hi),fill = "orange",alpha = 0.1)
    plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI40lo,
      ymax = CI40hi),fill = "orange",alpha = 0.1)
    plotobj6 <- plotobj6 + geom_ribbon(aes(x = time/24/7,ymin = CI20lo,
      ymax = CI20hi),fill = "orange",alpha = 0.1)
    plotobj6 <- plotobj6 + geom_line(aes(x = time/24/7,y = med),colour = "orange")
    plotobj6 <- plotobj6 + scale_y_continuous("Diastolic Blood Pressure (mmHg)",
      lim = c(50,120))
    plotobj6 <- plotobj6 + scale_x_continuous("Time (weeks)")
    print(plotobj6)

    ggsave(plot = plotobj6,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_BP.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(summary.BP, .(study), plot.bp.function)

# Summarise and plot probability of hand-foot syndrome grade over time
  nid <- length(unique(early.data$ID))
  nsim <- length(unique(early.data$SIM))
  ntotal <- nid*nsim
  summary.HFS <- ddply(early.data, .(study,time,HFS),
    function(early.data) summary.count.function(early.data$HFS))
  summary.HFS$HFS <- as.factor(summary.HFS$HFS)
  plot.hfs.function <- function(x) {
    plotobj7 <- NULL
    plotobj7 <- ggplot(x)
    plotobj7 <- plotobj7 + geom_line(aes(x = time/24/7,y = pro,colour = HFS))
    plotobj7 <- plotobj7 + scale_y_continuous("Probability of Hand-Foot Syndrome Grade",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj7 <- plotobj7 + scale_x_continuous("Time (weeks)")
    print(plotobj7)

    ggsave(plot = plotobj7,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_HFS.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(summary.HFS, .(study), plot.hfs.function)

# Summarise and plot probability of fatigue grade over time
  summary.FAT <- ddply(early.data, .(study,time,FAT),
    function(early.data) summary.count.function(early.data$FAT))
  summary.FAT$FAT <- as.factor(summary.FAT$FAT)
  plot.fat.function <- function(x) {
    plotobj8 <- NULL
    plotobj8 <- ggplot(x)
    plotobj8 <- plotobj8 + geom_line(aes(x = time/24/7,y = pro,colour = FAT))
    plotobj8 <- plotobj8 + scale_y_continuous("Probability of Fatigue Grade",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj8 <- plotobj8 + scale_x_continuous("Time (weeks)")
    print(plotobj8)

    ggsave(plot = plotobj8,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_FAT.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(summary.FAT, .(study), plot.fat.function)

# ------------------------------------------------------------------------------
# Subset individuals into "above" or "below" target trough by day 28
  trough.target <- function(pd.data) {
    target.lower <- 0.05	# mg/L
    target.upper <- 0.1	# mg/L
    t.attain <- 0
    if (pd.data$IPRE[pd.data$time == 28*24] >= target.lower &
      pd.data$IPRE[pd.data$time == 28*24] < target.upper) t.attain <- 1
    if (pd.data$IPRE[pd.data$time == 28*24] >= target.upper) t.attain <- 2
    pd.data$trough.target <- t.attain
    pd.data
  }
  trough.data <- ddply(pd.data, .(study,ID), trough.target)
  trough.survival <- ddply(trough.data, .(study,time,trough.target),
    pro.alive.function)
  trough.survival$trough.target <- as.factor(trough.survival$trough.target)
  plot.trough.os <- function(x) {
    plotobj9 <- NULL
    plotobj9 <- ggplot(x)
    plotobj9 <- plotobj9 + geom_step(aes(x = time/24/7,y = pro.alive,
      colour = trough.target))
    plotobj9 <- plotobj9 + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj9 <- plotobj9 + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(survival.data$time)/24/7,by = 24),
      lim = c(0,max(survival.data$time/24/7)))
    print(plotobj9)

    ggsave(plot = plotobj9,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_os_trough.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(trough.survival, .(study), plot.trough.os)

# Subset individuals into "above" or "below" target AUC by day 28
  auc.target <- function(pd.data) {
    target.lower <- 1.4	# mg*h/L
    target.upper <- 2.6
    t.attain <- 0
    if (pd.data$AUC24[pd.data$time == 28*24] >= target.lower &
      pd.data$AUC24[pd.data$time == 28*24] < target.upper) t.attain <- 1
    if (pd.data$AUC24[pd.data$time == 28*24] >= target.upper) t.attain <- 2
    pd.data$auc.target <- t.attain
    pd.data
  }
  auc.data <- ddply(pd.data, .(study,ID), auc.target)
  auc.survival <- ddply(auc.data, .(study,time,auc.target),
    pro.alive.function)
  auc.survival$auc.target <- as.factor(auc.survival$auc.target)
  plot.auc.os <- function(x) {
    plotobj10 <- NULL
    plotobj10 <- ggplot(x)
    plotobj10 <- plotobj10 + geom_step(aes(x = time/24/7,y = pro.alive,
      colour = auc.target))
    plotobj10 <- plotobj10 + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj10 <- plotobj10 + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(survival.data$time)/24/7,by = 24),
      lim = c(0,max(survival.data$time/24/7)))
    print(plotobj10)

    ggsave(plot = plotobj10,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_os_auc.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }
  ddply(auc.survival, .(study), plot.auc.os)
