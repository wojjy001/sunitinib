# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J,
# Karlsson MO, Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT
# as Predictors of Tumor Dynamics and Overall Survival Following Sunitinib
# Treatment in GIST. CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9
# ------------------------------------------------------------------------------
# Define model
bsur.code <- '
$SET			atol = 1e-8, rtol = 1e-8
          maxsteps = 100000000

$CMT			// List of compartments
          // Initial conditions specified in $MAIN
          VEGFR3,
          SURVIVAL,
          DROPOUT

$PARAM		// Population parameters
          // SURVIVAL
          LAMBH = 0.00596,	// Scale parameter in the Weibull density function for the survival model
          ALPHH = 1.226,	// Shape parameter in the Weibull density function for the survival model
          VEGFR3HAZ = -3.77,	// Parameter relating VEGFR3 to the hazard
          OBASEHAZ = -0.002371,	// Parameter relating baseline tumour size to the hazard

          // DROP OUT
          LAMBD = 0.00196,	// Scale parameter in the Weibull density function for the drop out model
          ALPHD = 1.273,	// Shape parameter in the Weibull density function for the drop out model

          // Individual biomarker parameters
          // VEGFR3
          IMAX = 1,
          BM03 = 44085,
          IC503 = 0.94008,
          MRT3 = 392.95,

          // Covariates
          OBASE = 200,	// Observed tumour size at baseline

          // Pharmacokinetic parameters
          DOSE = 50,
          CL = 34.5

$MAIN			// Individual parameter values
          double iIC503 = IC503*24*7;
          double iMRT3 = MRT3/24/7;

          // Compartment initial conditions
          VEGFR3_0 = BM03;
          SURVIVAL_0 = 0;
          DROPOUT_0 = 0;

          // Rate constants
          double KOUT3 = 1/iMRT3;

          // Sunitinib exposure
          double AUC = (DOSE/CL)*24*7;

$ODE			// VEGFR3
          double EFF3 = IMAX*AUC/(iIC503+AUC);
          double KIN3 = BM03*KOUT3;
          dxdt_VEGFR3 = KIN3*(1-EFF3)-KOUT3*VEGFR3;
          double RBM03 = (VEGFR3-BM03)/BM03;	// Relative change in VEGFR3 from baseline

          // SURVIVAL
          dxdt_SURVIVAL = 0;
          if (SOLVERTIME > 4) dxdt_SURVIVAL = pow(LAMBH*ALPHH*SOLVERTIME,ALPHH-1)*exp(-VEGFR3HAZ*RBM03-OBASEHAZ*OBASE);	// Patients who died before the first follow-up (week 4) were not included in the study

          // DROPOUT
          dxdt_DROPOUT = 0;
          if (SOLVERTIME > 4) dxdt_DROPOUT = pow(LAMBD*ALPHD*SOLVERTIME,ALPHD-1);

$TABLE		double iVEGFR3 = VEGFR3;	// VEGFR3 time-course
          double iRBM03	= (VEGFR3-BM03)/BM03;	// Relative change in VEGFR3 from baseline
          double CHZ = SURVIVAL;	// Cumulative hazard for survival
          double CHZD = DROPOUT;	// Cumulative hazard for drop out
          double SUR = exp(-CHZ);	// Probability of survival
          double SURD = exp(-CHZD);	// Probability of drop out
          double HAZNOW = 0;
          if (TIME > 4) HAZNOW = pow(LAMBH*ALPHH*TIME,ALPHH-1)*exp(-VEGFR3HAZ*iRBM03-OBASEHAZ*OBASE);
          double HAZDNOW = 0;
          if (TIME > 4) HAZDNOW = pow(LAMBD*ALPHD*TIME,ALPHD-1);
          double PEVENT = SUR*HAZNOW;	// Probability of event at time = TIME
          double PDROP = SURD*HAZDNOW;	// Probability of drop out at time = TIME

$CAPTURE	DOSE CL CHZ CHZD SUR SURD HAZNOW HAZDNOW PEVENT PDROP
'
# Compile the model code
  bsur.mod <- mcode("popSUR",bsur.code)
