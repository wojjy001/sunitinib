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
  standard.pd <- read.csv(file = "standard_pd_data.csv")
# Read in the PD data for target_auc_exact
  target.auc.exact.dir <- paste0(output.dir,"/target_auc_exact")
  setwd(target.auc.exact.dir)
  target.auc.exact.pd <- read.csv(file = "target_auc_exact_pd_data.csv")
# Read in the PD data for target_auc_round
  # Doses capped between 12.5 and 87.5 mg
  target.auc.round.dir <- paste0(output.dir,"/target_auc_round")
  setwd(target.auc.round.dir)
  target.auc.round.pd <- read.csv(file = "target_auc_round_pd_data.csv")
# Read in the PD data for target_trough_exact
  target.trough.exact.dir <- paste0(output.dir,"/target_trough_exact")
  setwd(target.trough.exact.dir)
  target.trough.exact.pd <- read.csv(file = "target_trough_exact_pd_data.csv")
# Read in the PD data for target_trough_round
  # Doses capped between 12.5 and 87.5 mg
  target.trough.round.dir <- paste0(output.dir,"/target_trough_round")
  setwd(target.trough.round.dir)
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
# Apply function to all data frames
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
    plotobj1 <- plotobj1 + ggtitle("Exact Dose to Achieve Target Administered")
    plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Standard Dosing (50 mg)"),
      data = standard.surv,alpha = 0.3)
    plotobj1 <- plotobj1 + geom_line(aes(x = time/24/7,y = surv,
      colour = "Standard Dosing (50 mg)"),data = standard.surv)
    plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Target Trough (50 ng/mL)"),
      data = target.trough.exact.surv,
      alpha = 0.3)
    plotobj1 <- plotobj1 + geom_line(aes(x = time/24/7,y = surv,
      colour = "Target Trough (50 ng/mL)"),data = target.trough.exact.surv)
    plotobj1 <- plotobj1 + geom_ribbon(aes(x = time/24/7,ymin = lower,
      ymax = upper,fill = "Target AUC (1.5 mg*h/L)"),
      data = target.auc.exact.surv,
      alpha = 0.3)
    plotobj1 <- plotobj1 + geom_line(aes(x = time/24/7,y = surv,
      colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.exact.surv)
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
    ggsave(plot = plotobj1,
      filename = paste0("exactdose_overallsurvival.png"),
      height = 15,width = 25,units = "cm",dpi = 300)

# Rounded doses TDM versus standard dose
  plotobj2 <- NULL
  plotobj2 <- ggplot()
  plotobj2 <- plotobj2 + ggtitle("Dose Rounded to Nearest 12.5 mg")
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    ymax = upper,fill = "Standard Dosing (50 mg)"),
    data = standard.surv,alpha = 0.3)
  plotobj2 <- plotobj2 + geom_line(aes(x = time/24/7,y = surv,
    colour = "Standard Dosing (50 mg)"),data = standard.surv)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    ymax = upper,fill = "Target Trough (50 ng/mL)"),
    data = target.trough.round.surv,
    alpha = 0.3)
  plotobj2 <- plotobj2 + geom_line(aes(x = time/24/7,y = surv,
    colour = "Target Trough (50 ng/mL)"),data = target.trough.round.surv)
  plotobj2 <- plotobj2 + geom_ribbon(aes(x = time/24/7,ymin = lower,
    ymax = upper,fill = "Target AUC (1.5 mg*h/L)"),
    data = target.auc.round.surv,
    alpha = 0.3)
  plotobj2 <- plotobj2 + geom_line(aes(x = time/24/7,y = surv,
    colour = "Target AUC (1.5 mg*h/L)"),data = target.auc.round.surv)
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
  ggsave(plot = plotobj2,
    filename = paste0("rounddose_overallsurvival.png"),
    height = 15,width = 25,units = "cm",dpi = 300)

  plotobj3 <- plotobj1 + theme(legend.position = "none")
  plotobj4 <- grid.arrange(grobs = list(plotobj3,plotobj2),
    heights = c(5,6))
  ggsave(plot = plotobj4,filename = "overallsurvival.png",
    width = 11.2,height = 12.6,units = "in",dpi = 300)
