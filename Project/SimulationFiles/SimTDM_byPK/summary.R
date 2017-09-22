# Read in simulated pharmacokinetic pharmacodynamic data
# Summarise overall survival
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define simulation file directory for the project
  file.dir <- "/Volumes/Prosecutor/sunitinib/Project/SimulationFiles/"
  project.name <- "SimTDM_byPK"
  project.dir <- paste0(file.dir,project.name)
  setwd(project.dir)
# Functions file (package libraries, summary functions, plotting)
  source("functions.R")
# Times file (particularly the "on.times" and "pd.times")
  source("times.R")

# ------------------------------------------------------------------------------
# Define output directory and read simulated pharmacodynamic data
# Output is saved into a totally different folder from scripts because my code
# is tracked by git and stored online.  Keeping all output (i.e., .csv and .png)
# files in the git repository slows down synchronisation
# PD output is saved in the same folder as PK output
  output.dir <- paste0("/Volumes/Prosecutor/sunitinib_nogit/",project.name)
# Read in the PD data for standard
  # Everyone receives 50 mg
  standard.dir <- paste0(output.dir,"/standard")
  setwd(standard.dir)
  standard.pk <- read.csv(file = "standard_pk_data.csv"))
  standard.pd <- read.csv(file = "standard_pd_data.csv")
# Read in the PD data for target_auc_exact
  target.auc.exact.dir <- paste0(output.dir,"/target_auc_exact")
  setwd(target.auc.exact.dir)
  target.auc.exact.pk <- read.csv(file = "target_auc_exact_pk_data.csv"))
  target.auc.exact.pd <- read.csv(file = "target_auc_exact_pd_data.csv")
# Read in the PD data for target_auc_round
  # Doses capped between 12.5 and 87.5 mg
  target.auc.round.dir <- paste0(output.dir,"/target_auc_round")
  setwd(target.auc.round.dir)
  target.auc.round.pk <- read.csv(file = "target_auc_round_pk_data.csv")
  target.auc.round.pd <- read.csv(file = "target_auc_round_pd_data.csv")
# Read in the PD data for target_trough_exact
  target.trough.exact.dir <- paste0(output.dir,"/target_trough_exact")
  setwd(target.trough.exact.dir)
  target.trough.exact.pk <- read.csv(file = "target_trough_exact_pk_data.csv")
  target.trough.exact.pd <- read.csv(file = "target_trough_exact_pd_data.csv")
# Read in the PD data for target_trough_round
  # Doses capped between 12.5 and 87.5 mg
  target.trough.round.dir <- paste0(output.dir,"/target_trough_round")
  setwd(target.trough.round.dir)
  target.trough.round.pk <- read.csv(file = "target_trough_round_pk_data.csv")
  target.trough.round.pd <- read.csv(file = "target_trough_round_pd_data.csv")

# ------------------------------------------------------------------------------
# Create a folder to save summary results
  save.name <- "summary"
  save.dir <- paste0(output.dir,"/",save.name)
  dir.create(save.dir)	# Create the folder if not done so
  setwd(save.dir)	# Set to working directory

# ------------------------------------------------------------------------------
# Overall survival
  # For overall survival plots, calculate the proportion of individuals alive
  # at each time-point
  # Kaplan Meier Plot for survival confidence intervals
  # All individuals have the same start time, i.e., time == 0
    km.function <- function(pd.data) {
      pd.data$start <- 0
      pd.data <- ddply(pd.data, .(SIM,ID), stop.time.function) # For each
      # individual calculate their stop time
      km.data <- ddply(pd.data, .(SIM,ID), headperID)	# First line per individual
      km.data <- km.data[c("SIM","ID","start","stop","event")]
      S <- Surv(time = km.data$start,time2 = km.data$stop,event = km.data$event)
      result <- survfit(formula = S ~ 1,data = km.data)
      cols <- lapply(2:12, function(x) summary(result)[x])
      surv.data <- do.call(data.frame, cols)
    }
# Apply function to all data frames)
  standard.surv <- km.function(standard.pd)
  target.auc.exact.surv <- km.function(target.auc.exact.pd)
  target.auc.round.surv <- km.function(target.auc.round.pd)
  target.trough.exact.surv <- km.function(target.trough.exact.pd)
  target.trough.round.surv <- km.function(target.trough.round.pd)

# ------------------------------------------------------------------------------
# Plot results
  # Exact doses TDM versus standard dose
    plotobj1 <- NULL
    plotobj1 <- ggplot()
    # plotobj1 <- plotobj1 + ggtitle("Exact Dose to Achieve Target Administered")
    # plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    #   ymax = upper,fill = "Standard Dosing (50 mg)"),
    #   data = standard.surv,alpha = 0.3)
    plotobj1 <- plotobj1 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Standard Dosing (50 mg)"),
      data = standard.surv)
    # plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    #   ymax = upper,fill = "Target Trough (50 ng/mL)"),
    #   data = target.trough.exact.surv,
    #   alpha = 0.3)
    plotobj1 <- plotobj1 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.surv)
    # plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    #   ymax = upper,fill = "Target AUC (1.5 mg*h/L)"),
    #   data = target.auc.exact.surv,
    #   alpha = 0.3)
    plotobj1 <- plotobj1 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.surv)
    plotobj1 <- plotobj1 + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj1 <- plotobj1 + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(standard.surv$time)/24/7,by = 30),
      lim = c(0,max(standard.surv$time/24/7)))
    plotobj1 <- plotobj1 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy")
    plotobj1 <- plotobj1 + theme(legend.position = "bottom")
    print(plotobj1)
    # ggsave(plot = plotobj1,
    #   filename = paste0("exactdose_overallsurvival.png"),
    #   height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses TDM versus standard dose
    plotobj2 <- NULL
    plotobj2 <- ggplot()
    # plotobj2 <- plotobj2 + ggtitle("Dose Rounded to Nearest 12.5 mg")
    # plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    #   ymax = upper,fill = "Standard Dosing (50 mg)"),
    #   data = standard.surv,alpha = 0.3)
    plotobj2 <- plotobj2 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Standard Dosing (50 mg)"),
      data = standard.surv)
    # plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    #   ymax = upper,fill = "Target Trough (50 ng/mL)"),
    #   data = target.trough.round.surv,
    #   alpha = 0.3)
    plotobj2 <- plotobj2 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Target Trough (50 ng/mL)"),
      data = target.trough.round.surv)
    # plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    #   ymax = upper,fill = "Target AUC (1.5 mg*h/L)"),
    #   data = target.auc.round.surv,
    #   alpha = 0.3)
    plotobj2 <- plotobj2 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.surv)
    plotobj2 <- plotobj2 + scale_y_continuous("Probability of Survival",
      lim = c(0,1),
      breaks = seq(from = 0,to = 1,by = 0.2),
      labels = seq(from = 0,to = 1,by = 0.2))
    plotobj2 <- plotobj2 + scale_x_continuous("Time (weeks)",
      breaks = seq(from = 0,to = max(standard.surv$time)/24/7,by = 30),
      lim = c(0,max(standard.surv$time/24/7)))
    plotobj2 <- plotobj2 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy")
    plotobj2 <- plotobj2 + theme(legend.position = "bottom")
    print(plotobj2)
    # ggsave(plot = plotobj2,
    #   filename = paste0("rounddose_overallsurvival.png"),
    #   height = 15,width = 25,units = "cm",dpi = 300)

# Combined overall survival plots
  plotobj3 <- plotobj1 + theme(legend.position = "none")
  plotobj3 <- plotobj3 + ggtitle("(e)")
  plotobj4 <- plotobj2 + theme(legend.position = "none")
  plotobj4 <- plotobj4 + ggtitle("(f)")
  plotobj5 <- grid.arrange(grobs = list(plotobj3,plotobj4),
    widths = c(5,5),ncol = 2)
  ggsave(plot = plotobj5,filename = "overallsurvival.png",
    width = 17.4,height = 7,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot concentrations
  # Exact doses
    plotobj6 <- NULL
    plotobj6 <- ggplot()
  # Standard
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pk,
      geom = "line",fun.y = median)
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pk,
      geom = "line",fun.y = median)
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj6 <- plotobj6 + geom_hline(aes(yintercept = 50),
      linetype = "dashed")
    plotobj6 <- plotobj6 + scale_y_continuous(
      "Combined Sunitinib Concentrations (ng/mL)",
      breaks = seq(0,300,25))
    plotobj6 <- plotobj6 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj6 <- plotobj6 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    print(plotobj6)
    ggsave(plot = plotobj6,
      filename = paste0("exactdose_IPRE.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj7 <- NULL
    plotobj7 <- ggplot()
  # Standard
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pk,
      geom = "line",fun.y = median)
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pk,
      geom = "line",fun.y = median)
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj7 <- plotobj7 + geom_hline(aes(yintercept = 50),
      linetype = "dashed")
    plotobj7 <- plotobj7 + scale_y_continuous(
      "Combined Sunitinib Concentrations (ng/mL)",
      breaks = seq(0,300,25))
    plotobj7 <- plotobj7 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj7 <- plotobj7 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    print(plotobj7)
    ggsave(plot = plotobj7,
      filename = paste0("rounddose_IPRE.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# Combined concentration plots
  plotobj8 <- plotobj6 + theme(legend.position = "none")
  plotobj8 <- plotobj8 + ggtitle("(a)")
  plotobj9 <- plotobj7 + theme(legend.position = "none")
  plotobj9 <- plotobj9 + ggtitle("(b)")
  plotobj10 <- grid.arrange(grobs = list(plotobj8,plotobj9),
    widths = c(5,5),ncol = 2)
  ggsave(plot = plotobj10,filename = "IPRE.png",
    width = 17.4,height = 7,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot AUC
  # Exact doses
    plotobj11 <- NULL
    plotobj11 <- ggplot()
  # Standard
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pk,
      geom = "line",fun.y = median)
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pk,
      geom = "line",fun.y = median)
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj11 <- plotobj11 + geom_hline(aes(yintercept = 1.5),
      linetype = "dashed")
    plotobj11 <- plotobj11 + scale_y_continuous(
      "Combined Sunitinib 24-hour AUC (mg*h/L)",
      breaks = seq(0,10,0.5))
    plotobj11 <- plotobj11 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj11 <- plotobj11 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    print(plotobj11)
    ggsave(plot = plotobj11,
      filename = paste0("exactdose_AUC24.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj12 <- NULL
    plotobj12 <- ggplot()
  # Standard
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pk,
      geom = "line",fun.y = median)
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pk,
      geom = "line",fun.y = median)
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj12 <- plotobj12 + geom_hline(aes(yintercept = 1.5),
      linetype = "dashed")
    plotobj12 <- plotobj12 + scale_y_continuous(
      "Combined Sunitinib 24-hour AUC (mg*h/L)",
      breaks = seq(0,10,0.5))
    plotobj12 <- plotobj12 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj12 <- plotobj12 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    print(plotobj12)
    ggsave(plot = plotobj12,
      filename = paste0("rounddose_AUC24.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# Combined concentration plots
  plotobj13 <- plotobj11 + theme(legend.position = "none")
  plotobj13 <- plotobj13 + ggtitle("(c)")
  plotobj14 <- plotobj12 + theme(legend.position = "none")
  plotobj14 <- plotobj14 + ggtitle("(d)")
  plotobj15 <- grid.arrange(grobs = list(plotobj13,plotobj14),
    widths = c(5,5),ncol = 2)
  ggsave(plot = plotobj15,filename = "AUC24.png",
    width = 17.4,height = 7,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Combine ALL plots
  plotobj16 <- grid.arrange(grobs = list(plotobj10,plotobj15,plotobj5),
    ncol = 1)
  ggsave(plot = plotobj16,filename = "Standard_v_TDM.pdf",
    width = 17.4,height = 21,units = "cm",dpi = 300)
  ggsave(plot = plotobj16,filename = "Standard_v_TDM.png",
    width = 17.4,height = 21,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Combine plots of ONLY rounded doses
# Not much difference between exact and rounded doses
# plotobj7, plotobj12, plotobj2
  plotobj17 <- plotobj7 + theme(legend.position = "none")
  plotobj17 <- plotobj17 + ggtitle("(a)")
  plotobj18 <- plotobj12 + theme(legend.position = "none")
  plotobj18 <- plotobj18 + ggtitle("(b)")
  plotobj19 <- plotobj2 + ggtitle("(c)")
  plotobj19 <- plotobj19 + guides(col = guide_legend(ncol = 1,byrow = TRUE))
  plotobj20 <- grid.arrange(plotobj17,plotobj18,plotobj19,ncol = 1,
    heights = c(5,5,5.5))
  ggsave(plot = plotobj20,filename = "Standard_v_TDMround.pdf",
    width = 8.4,height = 21,units = "cm",dpi = 300)
  ggsave(plot = plotobj20,filename = "Standard_v_TDMround.png",
    width = 8.4,height = 21,units = "cm",dpi = 300)
