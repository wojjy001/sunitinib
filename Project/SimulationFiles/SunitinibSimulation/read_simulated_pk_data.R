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
  # folder.list <- list.dirs(sim.dir)
  # folder.list <- folder.list[2:length(folder.list)]
  # for (i in 1:length(folder.list)) {
  #   folder.list[i] <- str_split(folder.list[i],pattern = sim.dir)
  #   folder.list[i] <- folder.list[[i]][2]
  # }
  # folder.list <- unlist(folder.list)
  folder.list <- "/standard_dose02"
# Read in .csv files from each folder
  read.pk.data <- function(x) {
    pk.data <- read.csv(file = paste0(sim.dir,x,"/",x,"_pk_data.csv"))
    pk.data$study <- x
    pk.data
  }
  pk.data <- ldply(folder.list,read.pk.data,.progress = "text")
  early.data <- pk.data[pk.data$cyc <= 2,]

# ------------------------------------------------------------------------------
# Recreate IPRE versus time plots so that y-axes are the same
# Summarise the data for median line and confidence intervals
  summary.IPRE <- ddply(early.data, .(study,time),
    function(early.data) graded.summary(early.data$IPRE))

  summarise.IPRE <- function(x) {
    plotobj <- NULL
    plotobj <- ggplot(x)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI90lo,
      ymax = CI90hi),fill = "red",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI80lo,
      ymax = CI80hi),fill = "red",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI60lo,
      ymax = CI60hi),fill = "red",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI40lo,
      ymax = CI40hi),fill = "red",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI20lo,
      ymax = CI20hi),fill = "red",alpha = 0.1)
    plotobj <- plotobj + geom_line(aes(x = time/24,y = med),colour = "red")
    plotobj <- plotobj + geom_hline(aes(yintercept = 0.05),
      linetype = "dashed")
    plotobj <- plotobj + geom_hline(aes(yintercept = 0.1),
      linetype = "dashed")
    plotobj <- plotobj + scale_y_continuous("Total Sunitinib Concentration (mg/L)",
      breaks = seq(from = 0,to = 0.15,by = 0.025),
      labels = seq(from = 0,to = 0.15,by = 0.025),
      lim = c(0,0.15))
    plotobj <- plotobj + scale_x_continuous("Time (days)")
    print(plotobj)

    ggsave(plot = plotobj,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_IPREvtime.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }

  ddply(summary.IPRE, .(study), summarise.IPRE)

# ------------------------------------------------------------------------------
# Recreate AUC24 versus time plots so that y-axes are the same
# Summarise the data for median line and confidence intervals
  summary.AUC24 <- ddply(early.data, .(study,time),
    function(early.data) graded.summary(early.data$AUC24))

  summarise.AUC24 <- function(x) {
    plotobj <- NULL
    plotobj <- ggplot(x)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI90lo,
      ymax = CI90hi),fill = "blue",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI80lo,
      ymax = CI80hi),fill = "blue",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI60lo,
      ymax = CI60hi),fill = "blue",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI40lo,
      ymax = CI40hi),fill = "blue",alpha = 0.1)
    plotobj <- plotobj + geom_ribbon(aes(x = time/24,ymin = CI20lo,
      ymax = CI20hi),fill = "blue",alpha = 0.1)
    plotobj <- plotobj + geom_line(aes(x = time/24,y = med),colour = "blue")
    plotobj <- plotobj + geom_hline(aes(yintercept = 1.4),
      linetype = "dashed")
    plotobj <- plotobj + geom_hline(aes(yintercept = 2.6),
      linetype = "dashed")
    plotobj <- plotobj + scale_y_continuous("Total Sunitinib 24-hour AUC (mg*h/L)",
      breaks = seq(from = 0,to = 4,by = 0.5),
      labels = seq(from = 0,to = 4,by = 0.5),
      lim = c(0,4))
    plotobj <- plotobj + scale_x_continuous("Time (days)")
    print(plotobj)

    ggsave(plot = plotobj,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_AUC24vtime.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }

  ddply(summary.AUC24, .(study), summarise.AUC24)

# ------------------------------------------------------------------------------
# Plot trough versus AUC24 at the end of the on-periods for each study
  trough.AUC24 <- function(x) {
    plotobj <- NULL
    plotobj <- ggplot(x[x$time %in% c((42+27)*24),])
    plotobj <- plotobj + geom_point(aes(x = IPRE,y = AUC24),
      colour = "blue",shape = 1,alpha = 0.3)
    plotobj <- plotobj + geom_hline(aes(yintercept = 1.4),
      linetype = "dashed")
    plotobj <- plotobj + geom_hline(aes(yintercept = 2.6),
      linetype = "dashed")
    plotobj <- plotobj + geom_vline(aes(xintercept = 0.05),
      linetype = "dashed")
    plotobj <- plotobj + geom_vline(aes(xintercept = 0.1),
      linetype = "dashed")
    plotobj <- plotobj + scale_y_continuous("Total Sunitinib 24-hour AUC (mg*h/L)",
      breaks = seq(from = 0,to = 6.5,by = 1),
      labels = seq(from = 0,to = 6.5,by = 1),
      lim = c(0,6.5))
    plotobj <- plotobj + scale_x_continuous("Total Sunitinib Concentration (mg/L)",
      breaks = seq(from = 0,to = 0.3,by = 0.05),
      labels = seq(from = 0,to = 0.3,by = 0.05),
      lim = c(0,0.3))
    print(plotobj)

    ggsave(plot = plotobj,filename = paste0(sim.dir,x$study[1],x$study[1],
      "_AUC24vIPRE.png"),
      width = 20,height = 15,unit = "cm",dpi = 300)
  }

  ddply(early.data, .(study), trough.AUC24)
