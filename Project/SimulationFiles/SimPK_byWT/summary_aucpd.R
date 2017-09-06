# Read in simulated pharmacokinetic pharmacodynamic data
# Summarise exposure-response relationships
# ------------------------------------------------------------------------------
# Remove objects in workspace
  rm(list = ls(all = TRUE))
# Define simulation file directory for the project
  file.dir <- "/Volumes/Prosecutor/sunitinib/Project/SimulationFiles/"
  project.name <- "SimPK_byWT"
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
# Read in the PK and PD data for standard
  # Everyone receives 50 mg
  standard.dir <- paste0(output.dir,"/standard")
  setwd(standard.dir)
  standard.pk <- read.csv(file = "standard_pk_data.csv")
  standard.pd <- read.csv(file = "standard_pd_data.csv")
# Read in the PK and PD data for mgkg_002_exact
  # 0.714 mg/kg dosing, doses not rounded to the nearest 12.5 mg
  # Doses capped between 12.5 and 87.5 mg
  mgkg.exact.dir <- paste0(output.dir,"/mgkg_002_exact")
  setwd(mgkg.exact.dir)
  mgkg.exact.pk <- read.csv(file = "mgkg_002_exact_pk_data.csv")
  mgkg.exact.pd <- read.csv(file = "mgkg_002_exact_pd_data.csv")

# Subset data for just week 50
  assess.time <- 50*7*24
  assess.standard.pd <- standard.pd[standard.pd$time == assess.time,]
  assess.mgkg.exact.pd <- mgkg.exact.pd[mgkg.exact.pd$time == assess.time,]

# ------------------------------------------------------------------------------
# Create a folder to save summary results
  save.name <- "summary_aucpd"
  save.dir <- paste0(output.dir,"/",save.name)
  dir.create(save.dir)	# Create the folder if not done so
  setwd(save.dir)	# Set to working directory

# ------------------------------------------------------------------------------
# Plot AUC24 versus PD outcomes
# Standard
  melt.standard.pd <- melt(data = assess.standard.pd,
    id.vars = c("SIM","ID","WT","AUC24","IPRE"),
    measure.vars = c("TUMOUR","ANC","BP","IPRE_VEGFR3",
    "IPRE_SKIT","HFS","FAT")
  )
# Perform linear regression of PD outcomes versus AUC24
  lm.tumour.auc24 <- lm(TUMOUR ~ AUC24,data = assess.standard.pd)
  lm.anc.auc24 <- lm(ANC ~ AUC24,data = assess.standard.pd)
  lm.bp.auc24 <- lm(BP ~ AUC24,data = assess.standard.pd)
  lm.vegfr3.auc24 <- lm(IPRE_VEGFR3 ~ AUC24,data = assess.standard.pd)
  lm.skit.auc24 <- lm(IPRE_SKIT ~ AUC24,data = assess.standard.pd)
  R2 <- data.frame(
    variable = c("TUMOUR","ANC","BP","IPRE_VEGFR3","IPRE_SKIT"),
    r2 = c(
      summary(lm.tumour.auc24)[["adj.r.squared"]],
      summary(lm.anc.auc24)[["adj.r.squared"]],
      summary(lm.bp.auc24)[["adj.r.squared"]],
      summary(lm.vegfr3.auc24)[["adj.r.squared"]],
      summary(lm.skit.auc24)[["adj.r.squared"]]
    )
  )
  R2$r2 <- round(R2$r2,digits = 3)
  R2$r2 <- paste0("R-squared: ",R2$r2)
  melt.standard.pd <- merge(melt.standard.pd,R2,by = c("variable"),all = T)
# Plot results
# Tumour
  plotobj1 <- NULL
  plotobj1 <- ggplot(melt.standard.pd[melt.standard.pd$variable == "TUMOUR",])
  plotobj1 <- plotobj1 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj1 <- plotobj1 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj1 <- plotobj1 + scale_y_continuous(
    "Tumour Size (Sum of Longest Diameters, mm)")
  plotobj1 <- plotobj1 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj1 <- plotobj1 + facet_wrap(~r2)
  # print(plotobj1)
# Absolute Neutrophil Count
  plotobj2 <- NULL
  plotobj2 <- ggplot(melt.standard.pd[melt.standard.pd$variable == "ANC",])
  plotobj2 <- plotobj2 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj2 <- plotobj2 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj2 <- plotobj2 + scale_y_continuous(
    "Absolute Neutrophil Count (x10^9)\n")
  plotobj2 <- plotobj2 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj2 <- plotobj2 + facet_wrap(~r2)
  # print(plotobj2)
# Diastolic Blood Pressure
  plotobj3 <- NULL
  plotobj3 <- ggplot(melt.standard.pd[melt.standard.pd$variable == "BP",])
  plotobj3 <- plotobj3 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj3 <- plotobj3 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj3 <- plotobj3 + scale_y_continuous(
    "Diastolic Blood Pressure (mmHg)")
  plotobj3 <- plotobj3 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj3 <- plotobj3 + facet_wrap(~r2)
  # print(plotobj3)
# sVEGFR-3
  plotobj4 <- NULL
  plotobj4 <- ggplot(
    melt.standard.pd[melt.standard.pd$variable == "IPRE_VEGFR3",])
  plotobj4 <- plotobj4 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj4 <- plotobj4 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj4 <- plotobj4 + scale_y_continuous(
    "sVEGFR-3 Concentration (pg/mL)\n")
  plotobj4 <- plotobj4 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj4 <- plotobj4 + facet_wrap(~r2)
  # print(plotobj4)
# sVEGFR-3
  plotobj5 <- NULL
  plotobj5 <- ggplot(
    melt.standard.pd[melt.standard.pd$variable == "IPRE_SKIT",])
  plotobj5 <- plotobj5 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj5 <- plotobj5 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj5 <- plotobj5 + scale_y_continuous(
    "sKIT Concentration (pg/mL)")
  plotobj5 <- plotobj5 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj5 <- plotobj5 + facet_wrap(~r2)
  # print(plotobj5)
  plotobj.auc24.standard <- grid.arrange(plotobj1,plotobj2,plotobj3,plotobj4,
    plotobj5,ncol = 3)
  ggsave(plot = plotobj.auc24.standard,filename = "auc24_standard.png")

# mgkg.exact
  melt.mgkg.exact.pd <- melt(data = assess.mgkg.exact.pd,
    id.vars = c("SIM","ID","WT","AUC24","IPRE"),
    measure.vars = c("TUMOUR","ANC","BP","IPRE_VEGFR3",
    "IPRE_SKIT","HFS","FAT")
  )
# Perform linear regression of PD outcomes versus AUC24
  lm.tumour.auc24 <- lm(TUMOUR ~ AUC24,data = assess.mgkg.exact.pd)
  lm.anc.auc24 <- lm(ANC ~ AUC24,data = assess.mgkg.exact.pd)
  lm.bp.auc24 <- lm(BP ~ AUC24,data = assess.mgkg.exact.pd)
  lm.vegfr3.auc24 <- lm(IPRE_VEGFR3 ~ AUC24,data = assess.mgkg.exact.pd)
  lm.skit.auc24 <- lm(IPRE_SKIT ~ AUC24,data = assess.mgkg.exact.pd)
  R2 <- data.frame(
    variable = c("TUMOUR","ANC","BP","IPRE_VEGFR3","IPRE_SKIT"),
    r2 = c(
      summary(lm.tumour.auc24)[["adj.r.squared"]],
      summary(lm.anc.auc24)[["adj.r.squared"]],
      summary(lm.bp.auc24)[["adj.r.squared"]],
      summary(lm.vegfr3.auc24)[["adj.r.squared"]],
      summary(lm.skit.auc24)[["adj.r.squared"]]
    )
  )
  R2$r2 <- round(R2$r2,digits = 3)
  R2$r2 <- paste0("R-squared: ",R2$r2)
  melt.mgkg.exact.pd <- merge(melt.mgkg.exact.pd,R2,by = c("variable"),all = T)
# Plot results
# Tumour
  plotobj6 <- NULL
  plotobj6 <- ggplot(melt.mgkg.exact.pd[melt.mgkg.exact.pd$variable == "TUMOUR",])
  plotobj6 <- plotobj6 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj6 <- plotobj6 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj6 <- plotobj6 + scale_y_continuous(
    "Tumour Size (Sum of Longest Diameters, mm)")
  plotobj6 <- plotobj6 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj6 <- plotobj6 + facet_wrap(~r2)
  # print(plotobj6)
# Absolute Neutrophil Count
  plotobj7 <- NULL
  plotobj7 <- ggplot(melt.mgkg.exact.pd[melt.mgkg.exact.pd$variable == "ANC",])
  plotobj7 <- plotobj7 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj7 <- plotobj7 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj7 <- plotobj7 + scale_y_continuous(
    "Absolute Neutrophil Count (x10^9)\n")
  plotobj7 <- plotobj7 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj7 <- plotobj7 + facet_wrap(~r2)
  # print(plotobj7)
# Diastolic Blood Pressure
  plotobj8 <- NULL
  plotobj8 <- ggplot(melt.mgkg.exact.pd[melt.mgkg.exact.pd$variable == "BP",])
  plotobj8 <- plotobj8 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj8 <- plotobj8 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj8 <- plotobj8 + scale_y_continuous(
    "Diastolic Blood Pressure (mmHg)")
  plotobj8 <- plotobj8 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj8 <- plotobj8 + facet_wrap(~r2)
  # print(plotobj8)
# sVEGFR-3
  plotobj9 <- NULL
  plotobj9 <- ggplot(
    melt.mgkg.exact.pd[melt.mgkg.exact.pd$variable == "IPRE_VEGFR3",])
  plotobj9 <- plotobj9 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj9 <- plotobj9 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj9 <- plotobj9 + scale_y_continuous(
    "sVEGFR-3 Concentration (pg/mL)\n")
  plotobj9 <- plotobj9 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj9 <- plotobj9 + facet_wrap(~r2)
  # print(plotobj9)
# sVEGFR-3
  plotobj10 <- NULL
  plotobj10 <- ggplot(
    melt.mgkg.exact.pd[melt.mgkg.exact.pd$variable == "IPRE_SKIT",])
  plotobj10 <- plotobj10 + geom_point(aes(x = AUC24,y = value),
    shape = 1,colour = "skyblue4")
  plotobj10 <- plotobj10 + geom_smooth(aes(x = AUC24,y = value),
    method = "lm",se = F,colour = "brown3")
  plotobj10 <- plotobj10 + scale_y_continuous(
    "sKIT Concentration (pg/mL)")
  plotobj10 <- plotobj10 + scale_x_continuous(
    "24-hour Total Sunitinib AUC (mg*h/L)")
  plotobj10 <- plotobj10 + facet_wrap(~r2)
  # print(plotobj10)
  plotobj.auc24.mgkg.exact <- grid.arrange(plotobj6,plotobj7,plotobj8,plotobj9,
    plotobj10,ncol = 3)
  ggsave(plot = plotobj.auc24.mgkg.exact,filename = "auc24_mgkg_exact.png")
