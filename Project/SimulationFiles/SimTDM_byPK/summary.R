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
  standard.pk <- read.csv(file = "standard_pk_data.csv")
  standard.pd <- read.csv(file = "standard_pd_data.csv")
# Read in the PD data for target_auc_exact
  target.auc.exact.dir <- paste0(output.dir,"/target_auc_exact")
  setwd(target.auc.exact.dir)
  target.auc.exact.pk <- read.csv(file = "target_auc_exact_pk_data.csv")
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
    plotobj1 <- plotobj1 + ggtitle("(c)")
    plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Standard Dosing (50 mg)"),
      data = standard.surv,alpha = 0.3)
    plotobj1 <- plotobj1 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Standard Dosing (50 mg)"),
      data = standard.surv)
    plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.surv,
      alpha = 0.3)
    plotobj1 <- plotobj1 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.surv)
    plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.surv,
      alpha = 0.3)
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
    plotobj1 <- plotobj1 + theme(legend.position = "none")
    print(plotobj1)
    ggsave(plot = plotobj1,
      filename = paste0("exactdose_overallsurvival.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses TDM versus standard dose
    plotobj2 <- NULL
    plotobj2 <- ggplot()
    plotobj2 <- plotobj2 + ggtitle("(d)")
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Standard Dosing (50 mg)"),
      data = standard.surv,alpha = 0.3)
    plotobj2 <- plotobj2 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Standard Dosing (50 mg)"),
      data = standard.surv)
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.surv,
      alpha = 0.3)
    plotobj2 <- plotobj2 + geom_step(aes(x = time/24/7,y = surv,
      colour = "Target Trough (50 ng/mL)"),
      data = target.trough.round.surv)
    plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.surv,
      alpha = 0.3)
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
    plotobj2 <- plotobj2 + theme(legend.position = "none")
    print(plotobj2)
    ggsave(plot = plotobj2,
      filename = paste0("rounddose_overallsurvival.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot concentrations
  # Exact doses
    plotobj3 <- NULL
    plotobj3 <- ggplot()
    plotobj3 <- plotobj3 + ggtitle("(a)")
  # Standard
    plotobj3 <- plotobj3 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj3 <- plotobj3 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj3 <- plotobj3 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pk,
      geom = "line",fun.y = median)
    plotobj3 <- plotobj3 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj3 <- plotobj3 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pk,
      geom = "line",fun.y = median)
    plotobj3 <- plotobj3 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj3 <- plotobj3 + geom_hline(aes(yintercept = 50),
      linetype = "dashed")
    plotobj3 <- plotobj3 + scale_y_continuous(
      "Trough Concentration (ng/mL)",
      breaks = seq(0,300,25))
    plotobj3 <- plotobj3 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj3 <- plotobj3 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj3 <- plotobj3 + theme(legend.position = "none")
    print(plotobj3)
    ggsave(plot = plotobj3,
      filename = paste0("exactdose_IPRE.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj4 <- NULL
    plotobj4 <- ggplot()
    plotobj4 <- plotobj4 + ggtitle("(c)")
  # Standard
    plotobj4 <- plotobj4 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj4 <- plotobj4 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj4 <- plotobj4 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pk,
      geom = "line",fun.y = median)
    plotobj4 <- plotobj4 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj4 <- plotobj4 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pk,
      geom = "line",fun.y = median)
    plotobj4 <- plotobj4 + stat_summary(aes(x = time/24/7,y = IPRE,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj4 <- plotobj4 + geom_hline(aes(yintercept = 50),
      linetype = "dashed")
    plotobj4 <- plotobj4 + scale_y_continuous(
      "Trough Concentration (ng/mL)",
      breaks = seq(0,300,25))
    plotobj4 <- plotobj4 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj4 <- plotobj4 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj4 <- plotobj4 + theme(legend.position = "none")
    print(plotobj4)
    ggsave(plot = plotobj4,
      filename = paste0("rounddose_IPRE.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot AUC
  # Exact doses
    plotobj5 <- NULL
    plotobj5 <- ggplot()
    plotobj5 <- plotobj5 + ggtitle("(b)")
  # Standard
    plotobj5 <- plotobj5 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj5 <- plotobj5 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj5 <- plotobj5 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pk,
      geom = "line",fun.y = median)
    plotobj5 <- plotobj5 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj5 <- plotobj5 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pk,
      geom = "line",fun.y = median)
    plotobj5 <- plotobj5 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj5 <- plotobj5 + geom_hline(aes(yintercept = 1.5),
      linetype = "dashed")
    plotobj5 <- plotobj5 + scale_y_continuous(
      "24-hour AUC (mg*h/L)",
      breaks = seq(0,10,0.5))
    plotobj5 <- plotobj5 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj5 <- plotobj5 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj5 <- plotobj5 + theme(legend.position = "none")
    print(plotobj5)
    ggsave(plot = plotobj5,
      filename = paste0("exactdose_AUC24.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj6 <- NULL
    plotobj6 <- ggplot()
    plotobj6 <- plotobj6 + ggtitle("(d)")
  # Standard
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)"),data = standard.pk,
      geom = "line",fun.y = median)
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pk,
      geom = "line",fun.y = median)
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pk,
      geom = "line",fun.y = median)
    plotobj6 <- plotobj6 + stat_summary(aes(x = time/24/7,y = AUC24,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pk,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj6 <- plotobj6 + geom_hline(aes(yintercept = 1.5),
      linetype = "dashed")
    plotobj6 <- plotobj6 + scale_y_continuous(
      "24-hour AUC (mg*h/L)",
      breaks = seq(0,10,0.5))
    plotobj6 <- plotobj6 + scale_x_continuous("Time (weeks)",
      lim = c(0,12),breaks = seq(0,12,3))
    plotobj6 <- plotobj6 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj6 <- plotobj6 + theme(legend.position = "none")
    print(plotobj6)
    ggsave(plot = plotobj6,
      filename = paste0("rounddose_AUC24.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot sVEGFR-3
  # Exact doses
    plotobj7 <- NULL
    plotobj7 <- ggplot()
    plotobj7 <- plotobj7 + ggtitle("(e)")
  # Standard
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pd,
      geom = "line",fun.y = median)
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pd,
      geom = "line",fun.y = median)
    plotobj7 <- plotobj7 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj7 <- plotobj7 + scale_y_continuous(
      "sVEGFR-3 Concentration (pg/mL)",
      breaks = seq(9,12,0.5),labels = seq(9,12,0.5))
    plotobj7 <- plotobj7 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj7 <- plotobj7 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj7 <- plotobj7 + theme(legend.position = "none")
    print(plotobj7)
    ggsave(plot = plotobj7,
      filename = paste0("exactdose_SVEGFR3.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj8 <- NULL
    plotobj8 <- ggplot()
    plotobj8 <- plotobj8 + ggtitle("(f)")
  # Standard
    plotobj8 <- plotobj8 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj8 <- plotobj8 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj8 <- plotobj8 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pd,
      geom = "line",fun.y = median)
    plotobj8 <- plotobj8 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj8 <- plotobj8 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pd,
      geom = "line",fun.y = median)
    plotobj8 <- plotobj8 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR3,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj8 <- plotobj8 + scale_y_continuous(
      "sVEGFR-3 Concentration (pg/mL)",
      breaks = seq(9,12,0.5),labels = seq(9,12,0.5))
    plotobj8 <- plotobj8 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj8 <- plotobj8 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj8 <- plotobj8 + theme(legend.position = "none")
    print(plotobj8)
    ggsave(plot = plotobj8,
      filename = paste0("rounddose_SVEGFR3.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot sKIT
  # Exact doses
    plotobj9 <- NULL
    plotobj9 <- ggplot()
    plotobj9 <- plotobj9 + ggtitle("(g)")
  # Standard
    plotobj9 <- plotobj9 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj9 <- plotobj9 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj9 <- plotobj9 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pd,
      geom = "line",fun.y = median)
    plotobj9 <- plotobj9 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj9 <- plotobj9 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pd,
      geom = "line",fun.y = median)
    plotobj9 <- plotobj9 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj9 <- plotobj9 + scale_y_continuous(
      "sKIT Concentration (pg/mL)",
      breaks = seq(9,13,0.5),labels = seq(9,13,0.5))
    plotobj9 <- plotobj9 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj9 <- plotobj9 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj9 <- plotobj9 + theme(legend.position = "none")
    print(plotobj9)
    ggsave(plot = plotobj9,
      filename = paste0("exactdose_SKIT.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj10 <- NULL
    plotobj10 <- ggplot()
    plotobj10 <- plotobj10 + ggtitle("(h)")
  # Standard
    plotobj10 <- plotobj10 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj10 <- plotobj10 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj10 <- plotobj10 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pd,
      geom = "line",fun.y = median)
    plotobj10 <- plotobj10 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj10 <- plotobj10 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pd,
      geom = "line",fun.y = median)
    plotobj10 <- plotobj10 + stat_summary(aes(x = time/24/7,y = IPRE_SKIT,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj10 <- plotobj10 + scale_y_continuous(
      "sKIT Concentration (pg/mL)",
      breaks = seq(9,13,0.5),labels = seq(9,13,0.5))
    plotobj10 <- plotobj10 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj10 <- plotobj10 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj10 <- plotobj10 + theme(legend.position = "none")
    print(plotobj10)
    ggsave(plot = plotobj10,
      filename = paste0("rounddose_SKIT.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot dBP
  # Exact doses
    plotobj11 <- NULL
    plotobj11 <- ggplot()
    plotobj11 <- plotobj11 + ggtitle("(c)")
  # Standard
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pd,
      geom = "line",fun.y = median)
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pd,
      geom = "line",fun.y = median)
    plotobj11 <- plotobj11 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj11 <- plotobj11 + scale_y_continuous(
      "Diastolic Blood Pressure (mmHg)",
      breaks = seq(40,240,20),labels = seq(40,240,20))
    plotobj11 <- plotobj11 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj11 <- plotobj11 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj11 <- plotobj11 + theme(legend.position = "none")
    print(plotobj11)
    ggsave(plot = plotobj11,
      filename = paste0("exactdose_dBP.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj12 <- NULL
    plotobj12 <- ggplot()
    plotobj12 <- plotobj12 + ggtitle("(d)")
  # Standard
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pd,
      geom = "line",fun.y = median)
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pd,
      geom = "line",fun.y = median)
    plotobj12 <- plotobj12 + stat_summary(aes(x = time/24/7,y = BP,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj12 <- plotobj12 + scale_y_continuous(
      "Diastolic Blood Pressure (mmHg)",
      breaks = seq(40,240,20),labels = seq(40,240,20))
    plotobj12 <- plotobj12 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj12 <- plotobj12 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj12 <- plotobj12 + theme(legend.position = "none")
    print(plotobj12)
    ggsave(plot = plotobj12,
      filename = paste0("rounddose_dBP.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot ANC
  # Exact doses
    plotobj13 <- NULL
    plotobj13 <- ggplot()
    plotobj13 <- plotobj13 + ggtitle("(a)")
  # Standard
    plotobj13 <- plotobj13 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj13 <- plotobj13 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj13 <- plotobj13 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pd,
      geom = "line",fun.y = median)
    plotobj13 <- plotobj13 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj13 <- plotobj13 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pd,
      geom = "line",fun.y = median)
    plotobj13 <- plotobj13 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj13 <- plotobj13 + scale_y_log10(
      "Absolute Neutrophil Count (x10^9)",
      breaks = seq(1,10,1),labels = seq(1,10,1))
    plotobj13 <- plotobj13 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj13 <- plotobj13 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj13 <- plotobj13 + theme(legend.position = "none")
    print(plotobj13)
    ggsave(plot = plotobj13,
      filename = paste0("exactdose_ANC.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj14 <- NULL
    plotobj14 <- ggplot()
    plotobj14 <- plotobj14 + ggtitle("(b)")
  # Standard
    plotobj14 <- plotobj14 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj14 <- plotobj14 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj14 <- plotobj14 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pd,
      geom = "line",fun.y = median)
    plotobj14 <- plotobj14 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj14 <- plotobj14 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pd,
      geom = "line",fun.y = median)
    plotobj14 <- plotobj14 + stat_summary(aes(x = time/24/7,y = ANC,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj14 <- plotobj14 + scale_y_log10(
      "Absolute Neutrophil Count (x10^9)",
      breaks = seq(1,10,1),labels = seq(1,10,1))
    plotobj14 <- plotobj14 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj14 <- plotobj14 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj14 <- plotobj14 + theme(legend.position = "none")
    print(plotobj14)
    ggsave(plot = plotobj14,
      filename = paste0("rounddose_ANC.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot TUMOUR
  # Exact doses
    plotobj15 <- NULL
    plotobj15 <- ggplot()
    plotobj15 <- plotobj15 + ggtitle("(a)")
  # Standard
    plotobj15 <- plotobj15 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj15 <- plotobj15 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj15 <- plotobj15 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pd,
      geom = "line",fun.y = median)
    plotobj15 <- plotobj15 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj15 <- plotobj15 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pd,
      geom = "line",fun.y = median)
    plotobj15 <- plotobj15 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj15 <- plotobj15 + scale_y_continuous(
      "Tumour Size (Sum of Longest Diameters, mm)",
      breaks = seq(0,1000,100),labels = seq(0,1000,100))
    plotobj15 <- plotobj15 + scale_x_continuous("Time (weeks)",
      lim = c(0,48),breaks = seq(0,48,3))
    plotobj15 <- plotobj15 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj15 <- plotobj15 + theme(legend.position = "none")
    print(plotobj15)
    ggsave(plot = plotobj15,
      filename = paste0("exactdose_TUMOUR.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj16 <- NULL
    plotobj16 <- ggplot()
    plotobj16 <- plotobj16 + ggtitle("(b)")
  # Standard
    plotobj16 <- plotobj16 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj16 <- plotobj16 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj16 <- plotobj16 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pd,
      geom = "line",fun.y = median)
    plotobj16 <- plotobj16 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj16 <- plotobj16 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pd,
      geom = "line",fun.y = median)
    plotobj16 <- plotobj16 + stat_summary(aes(x = time/24/7,y = TUMOUR,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj16 <- plotobj16 + scale_y_continuous(
      "Tumour Size (Sum of Longest Diameters, mm)",
      breaks = seq(0,1000,100),labels = seq(0,1000,100))
    plotobj16 <- plotobj16 + scale_x_continuous("Time (weeks)",
      lim = c(0,48),breaks = seq(0,48,3))
    plotobj16 <- plotobj16 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj16 <- plotobj16 + theme(legend.position = "none")
    print(plotobj16)
    ggsave(plot = plotobj16,
      filename = paste0("rounddose_TUMOUR.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot VEGF
  # Exact doses
    plotobj17 <- NULL
    plotobj17 <- ggplot()
    plotobj17 <- plotobj17 + ggtitle("(a)")
  # Standard
    plotobj17 <- plotobj17 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj17 <- plotobj17 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj17 <- plotobj17 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pd,
      geom = "line",fun.y = median)
    plotobj17 <- plotobj17 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj17 <- plotobj17 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pd,
      geom = "line",fun.y = median)
    plotobj17 <- plotobj17 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj17 <- plotobj17 + scale_y_continuous(
      "VEGF Concentration (pg/mL)",
      breaks = seq(1,20,0.5),labels = seq(1,20,0.5))
    plotobj17 <- plotobj17 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj17 <- plotobj17 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj17 <- plotobj17 + theme(legend.position = "none")
    print(plotobj17)
    ggsave(plot = plotobj17,
      filename = paste0("exactdose_VEGF.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj18 <- NULL
    plotobj18 <- ggplot()
    plotobj18 <- plotobj18 + ggtitle("(b)")
  # Standard
    plotobj18 <- plotobj18 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj18 <- plotobj18 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj18 <- plotobj18 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pd,
      geom = "line",fun.y = median)
    plotobj18 <- plotobj18 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj18 <- plotobj18 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pd,
      geom = "line",fun.y = median)
    plotobj18 <- plotobj18 + stat_summary(aes(x = time/24/7,y = IPRE_VEGF,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj18 <- plotobj18 + scale_y_continuous(
      "VEGF Concentration (pg/mL)",
      breaks = seq(1,20,0.5),labels = seq(1,20,0.5))
    plotobj18 <- plotobj18 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj18 <- plotobj18 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj18 <- plotobj18 + theme(legend.position = "none")
    print(plotobj18)
    ggsave(plot = plotobj18,
      filename = paste0("rounddose_VEGF.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Plot VEGFR2
  # Exact doses
    plotobj19 <- NULL
    plotobj19 <- ggplot()
    plotobj19 <- plotobj19 + ggtitle("(c)")
  # Standard
    plotobj19 <- plotobj19 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj19 <- plotobj19 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj19 <- plotobj19 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.pd,
      geom = "line",fun.y = median)
    plotobj19 <- plotobj19 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj19 <- plotobj19 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.pd,
      geom = "line",fun.y = median)
    plotobj19 <- plotobj19 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj19 <- plotobj19 + scale_y_continuous(
      "sVEGFR-2 Concentration (pg/mL)",
      breaks = seq(1,20,0.5),labels = seq(1,20,0.5))
    plotobj19 <- plotobj19 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj19 <- plotobj19 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj19 <- plotobj19 + theme(legend.position = "none")
    print(plotobj19)
    ggsave(plot = plotobj19,
      filename = paste0("exactdose_VEGFR2.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

  # Rounded doses
    plotobj20 <- NULL
    plotobj20 <- ggplot()
    plotobj20 <- plotobj20 + ggtitle("(d)")
  # Standard
    plotobj20 <- plotobj20 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Standard Dosing (50 mg)"),data = standard.pd,
      geom = "line",fun.y = median)
    plotobj20 <- plotobj20 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Standard Dosing (50 mg)",fill = "Standard Dosing (50 mg)"),
      data = standard.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target Trough
    plotobj20 <- plotobj20 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.round.pd,
      geom = "line",fun.y = median)
    plotobj20 <- plotobj20 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target Trough (50 ng/mL)",fill = "Target Trough (50 ng/mL)"),
      data = target.trough.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
  # Target AUC
    plotobj20 <- plotobj20 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.pd,
      geom = "line",fun.y = median)
    plotobj20 <- plotobj20 + stat_summary(aes(x = time/24/7,y = IPRE_VEGFR2,
      colour = "Target AUC (1.5 mg*h/L)",fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.round.pd,
      geom = "ribbon",fun.ymin = CI95lo,fun.ymax = CI95hi,
      linetype = "dotted",alpha = 0.1)
    plotobj20 <- plotobj20 + scale_y_continuous(
      "sVEGFR-2 Concentration (pg/mL)",
      breaks = seq(1,20,0.5),labels = seq(1,20,0.5))
    plotobj20 <- plotobj20 + scale_x_continuous("Time (weeks)",
      lim = c(0,18),breaks = seq(0,18,3))
    plotobj20 <- plotobj20 + labs(fill = "Dosing Strategy",
      colour = "Dosing Strategy",linetype = "Dosing Strategy")
    plotobj20 <- plotobj20 + theme(legend.position = "none")
    print(plotobj20)
    ggsave(plot = plotobj20,
      filename = paste0("rounddose_VEGFR2.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Hand Foot Syndrome
  ntotal <- 1000
  count.hfs.standard.pd <- ddply(standard.pd, .(time,HFS),
    function(standard.pd) summary.count.function(standard.pd$HFS))
  count.hfs.target.trough.exact.pd <- ddply(target.trough.exact.pd, .(time,HFS),
    function(target.trough.exact.pd) summary.count.function(target.trough.exact.pd$HFS))
  count.hfs.target.trough.round.pd <- ddply(target.trough.round.pd, .(time,HFS),
    function(target.trough.round.pd) summary.count.function(target.trough.round.pd$HFS))
  count.hfs.target.auc.exact.pd <- ddply(target.auc.exact.pd, .(time,HFS),
    function(target.auc.exact.pd) summary.count.function(target.auc.exact.pd$HFS))
  count.hfs.target.auc.round.pd <- ddply(target.auc.round.pd, .(time,HFS),
    function(target.auc.round.pd) summary.count.function(target.auc.round.pd$HFS))

# Exact dose
  plotobj21 <- NULL
  plotobj21 <- ggplot()
  plotobj21 <- plotobj21 + ggtitle("(a)")
# Standard dosing
  plotobj21 <- plotobj21 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Standard Dosing (50 mg)"),data = count.hfs.standard.pd)
# Target Trough (50 ng/mL)
  plotobj21 <- plotobj21 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target Trough (50 ng/mL)"),data = count.hfs.target.trough.exact.pd)
# Target AUC (1.5 mg*h/L)
  plotobj21 <- plotobj21 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target AUC (1.5 mg*h/L)"),data = count.hfs.target.auc.exact.pd)
  plotobj21 <- plotobj21 + scale_y_continuous(
    "Probability of Grade",
    lim = c(0,1),breaks = seq(0,1,0.5),labels = seq(0,1,0.5))
  plotobj21 <- plotobj21 + scale_x_continuous("Time (weeks)",
    lim = c(0,48),breaks = seq(0,48,6))
  plotobj21 <- plotobj21 + facet_wrap(~HFS)
  plotobj21 <- plotobj21 + labs(fill = "Dosing Strategy",
    colour = "Dosing Strategy",linetype = "Dosing Strategy")
  plotobj21 <- plotobj21 + theme(legend.position = "none")
  print(plotobj21)
  ggsave(plot = plotobj21,
    filename = paste0("exactdose_HFS.png"),
    height = 15,width = 25,units = "cm",dpi = 300)

# round dose
  plotobj22 <- NULL
  plotobj22 <- ggplot()
  plotobj22 <- plotobj22 + ggtitle("(b)")
# Standard dosing
  plotobj22 <- plotobj22 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Standard Dosing (50 mg)"),data = count.hfs.standard.pd)
# Target Trough (50 ng/mL)
  plotobj22 <- plotobj22 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target Trough (50 ng/mL)"),data = count.hfs.target.trough.round.pd)
# Target AUC (1.5 mg*h/L)
  plotobj22 <- plotobj22 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target AUC (1.5 mg*h/L)"),data = count.hfs.target.auc.round.pd)
  plotobj22 <- plotobj22 + scale_y_continuous(
    "Probability of Grade",
    lim = c(0,1),breaks = seq(0,1,0.5),labels = seq(0,1,0.5))
  plotobj22 <- plotobj22 + scale_x_continuous("Time (weeks)",
    lim = c(0,48),breaks = seq(0,48,6))
  plotobj22 <- plotobj22 + facet_wrap(~HFS)
  plotobj22 <- plotobj22 + labs(fill = "Dosing Strategy",
    colour = "Dosing Strategy",linetype = "Dosing Strategy")
  plotobj22 <- plotobj22 + theme(legend.position = "none")
  print(plotobj22)
  ggsave(plot = plotobj22,
    filename = paste0("rounddose_HFS.png"),
    height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Fatigue
  count.fat.standard.pd <- ddply(standard.pd, .(time,FAT),
    function(standard.pd) summary.count.function(standard.pd$FAT))
  count.fat.target.trough.exact.pd <- ddply(target.trough.exact.pd, .(time,FAT),
    function(target.trough.exact.pd) summary.count.function(target.trough.exact.pd$FAT))
  count.fat.target.trough.round.pd <- ddply(target.trough.round.pd, .(time,FAT),
    function(target.trough.round.pd) summary.count.function(target.trough.round.pd$FAT))
  count.fat.target.auc.exact.pd <- ddply(target.auc.exact.pd, .(time,FAT),
    function(target.auc.exact.pd) summary.count.function(target.auc.exact.pd$FAT))
  count.fat.target.auc.round.pd <- ddply(target.auc.round.pd, .(time,FAT),
    function(target.auc.round.pd) summary.count.function(target.auc.round.pd$FAT))

  # Exact dose
  plotobj23 <- NULL
  plotobj23 <- ggplot()
  plotobj23 <- plotobj23 + ggtitle("(c)")
  # Standard dosing
  plotobj23 <- plotobj23 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Standard Dosing (50 mg)"),data = count.fat.standard.pd)
  # Target Trough (50 ng/mL)
  plotobj23 <- plotobj23 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target Trough (50 ng/mL)"),data = count.fat.target.trough.exact.pd)
  # Target AUC (1.5 mg*h/L)
  plotobj23 <- plotobj23 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target AUC (1.5 mg*h/L)"),data = count.fat.target.auc.exact.pd)
  plotobj23 <- plotobj23 + scale_y_continuous(
    "Probability of Grade",
    lim = c(0,1),breaks = seq(0,1,0.5),labels = seq(0,1,0.5))
  plotobj23 <- plotobj23 + scale_x_continuous("Time (weeks)",
    lim = c(0,48),breaks = seq(0,48,6))
  plotobj23 <- plotobj23 + facet_wrap(~FAT)
  plotobj23 <- plotobj23 + labs(fill = "Dosing Strategy",
    colour = "Dosing Strategy",linetype = "Dosing Strategy")
  plotobj23 <- plotobj23 + theme(legend.position = "none")
  print(plotobj23)
  ggsave(plot = plotobj23,
    filename = paste0("exactdose_FAT.png"),
    height = 15,width = 25,units = "cm",dpi = 300)

  # round dose
  plotobj24 <- NULL
  plotobj24 <- ggplot()
  plotobj24 <- plotobj24 + ggtitle("(d)")
  # Standard dosing
  plotobj24 <- plotobj24 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Standard Dosing (50 mg)"),data = count.fat.standard.pd)
  # Target Trough (50 ng/mL)
  plotobj24 <- plotobj24 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target Trough (50 ng/mL)"),data = count.fat.target.trough.round.pd)
  # Target AUC (1.5 mg*h/L)
  plotobj24 <- plotobj24 + geom_line(aes(x = time/24/7,y = pro,
    colour = "Target AUC (1.5 mg*h/L)"),data = count.fat.target.auc.round.pd)
  plotobj24 <- plotobj24 + scale_y_continuous(
    "Probability of Grade",
    lim = c(0,1),breaks = seq(0,1,0.5),labels = seq(0,1,0.5))
  plotobj24 <- plotobj24 + scale_x_continuous("Time (weeks)",
    lim = c(0,48),breaks = seq(0,48,6))
  plotobj24 <- plotobj24 + facet_wrap(~FAT)
  plotobj24 <- plotobj24 + labs(fill = "Dosing Strategy",
    colour = "Dosing Strategy",linetype = "Dosing Strategy")
  plotobj24 <- plotobj24 + theme(legend.position = "none")
  print(plotobj24)
  ggsave(plot = plotobj24,
    filename = paste0("rounddose_FAT.png"),
    height = 15,width = 25,units = "cm",dpi = 300)

# ------------------------------------------------------------------------------
# Combine plots together
# PK
  plotobj25 <- grid.arrange(grobs = list(plotobj3,plotobj5,plotobj4,plotobj6),
    ncol = 2)
  ggsave(plot = plotobj25,filename = "IPRE_AUC24.png",
    dpi = 300,units = "cm",width = 17.4,height = 15)
  ggsave(plot = plotobj25,filename = "IPRE_AUC24.pdf",
    dpi = 300,units = "cm",width = 17.4,height = 15)
# PD Biomarkers
  plotobj26 <- grid.arrange(grobs = list(plotobj17,plotobj18,plotobj19,
    plotobj20,plotobj7,plotobj8,plotobj9,plotobj10),ncol = 2)
  ggsave(plot = plotobj26,filename = "biomarkers.png",
    dpi = 300,units = "cm",width = 17.4,height = 21)
  ggsave(plot = plotobj26,filename = "biomarkers.pdf",
    dpi = 300,units = "cm",width = 17.4,height = 21)
# PD Adverse effects - continuous
  plotobj27 <- grid.arrange(grobs = list(plotobj13,plotobj14,plotobj11,
    plotobj12),ncol = 2)
  ggsave(plot = plotobj27,filename = "adverse_effects_cont.png",
    dpi = 300,units = "cm",width = 17.4,height = 15)
  ggsave(plot = plotobj27,filename = "adverse_effects_cont.pdf",
    dpi = 300,units = "cm",width = 17.4,height = 15)
# PD Adverse effects - categorical
  plotobj28 <- grid.arrange(grobs = list(plotobj21,plotobj22,plotobj23,
    plotobj24),ncol = 2)
  ggsave(plot = plotobj28,filename = "adverse_effects_cat.png",
    dpi = 300,units = "cm",width = 17.4,height = 15)
  ggsave(plot = plotobj28,filename = "adverse_effects_cat.pdf",
    dpi = 300,units = "cm",width = 17.4,height = 15)
# PD End-points
  plotobj29 <- grid.arrange(grobs = list(plotobj15,plotobj16,plotobj1,plotobj2),
    ncol = 2)
  ggsave(plot = plotobj29,filename = "end_points.png",
    dpi = 300,units = "cm",width = 17.4,height = 15)
  ggsave(plot = plotobj29,filename = "end_points.pdf",
    dpi = 300,units = "cm",width = 17.4,height = 15)
