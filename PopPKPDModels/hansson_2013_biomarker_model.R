# Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE, French J, Karlsson MO,
# Friberg LE. PKPD Modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT as Predictors
# of Tumor Dynamics and Overall Survival Following Sunitinib Treatment in GIST.
# CPT: pharmacometrics & systems pharmacology. 2013;2(11):1-9.
# ------------------------------------------------------------------------------
# Define model
code <- '
$SET		atol = 1e-8, rtol = 1e-8
        maxsteps = 100000

$CMT		// List of compartments
        // Initial conditions specified in $MAIN
        VEGF,
        VEGFR2,
        VEGFR3,
        SKIT

$PARAM	// Population parameters
        // Common for some biomarkers
        POPIMAX = 1,	// All 4
        POPIC50 = 1,	// All 4
        POPDPSLO = 0.0261/1000,	// Slope for the disease progression model (VEGF and SKIT), months^-1

        // VEGF
        POPBM0 = 59.8,	// Baseline VEGF, pg/mL
        POPMRT = 3.75,	// Mean residence time, days
        POPHILL = 3.31,	// Hill coefficient

        // VEGFR2
        POPBM02 = 8660,	// Baseline VEGFR2, pg/mL
        POPMRT2 = 23.1,	// Mean residence time, days
        POPHILL2 = 1.54,	// Hill coefficient

        // VEGFR3
        POPBM03 = 63900,	// Baseline VEGFR3, pg/mL
        POPMRT3 = 16.7,	// Mean residence time, days

        // SKIT
        POPBM0S = 39200,	// Baseline SKIT, pg/mL
        POPMRTS = 101,	// Mean residence time, days

        // Pre-defined between subject variability
        ETAMRT23 = 0,
        ETABM0 = 0,
        ETAIC50 = 0,
        ETADPSLO = 0,
        ETABM02 = 0,
        ETAIC502 = 0,
        ETABM03 = 0,
        ETAIC503 = 0,
        ETABM0S = 0,
        ETAMRTS = 0,
        ETAIC50S = 0,
        ETADPSLOS = 0,

        // Pharmacokinetic parameters
        DOSE = 50,
        CL = 34.5

$OMEGA	// Between-subject variability
        @name diag @annotated
        // Common for some biomarkers
        BSVMRT23 : 0.06 : ETA on mean residence time for VEGF, VEGFR2 and VEGFR3

        // VEGF
        BSVBM0: 0.252 : ETA on baseline VEGF
        BSVDPSLO: 2.95 : ETA on slope for the disease progression model

        // VEGFR2
        BSVBM02: 0.0369 : ETA on baseline VEGFR2

        // VEGFR3
        BSVBM03: 0.186 : ETA on baseline VEGFR3

        // SKIT
        BSVBM0S: 0.254 : ETA on baseline SKIT
        BSVMRTS: 0.0753: ETA on mean residence time for SKIT
        BSVDPSLOS: 3.01 : ETA on slope for the disease progression model

$OMEGA	// Between-subject variability
        @name block @annotated @block
        BSVIC50: 0.253 : ETA on VEGF IC50
        BSVIC502: 0.198 0.189 : ETA on VEGFR2 IC50
        BSVIC503: 0.238 0.252 0.398 : ETA on VEGFR3 IC50
        BSVIC50S: 0.218 0.297 0.936 5.77 : ETA on SKIT IC50

$MAIN		// Individual parameter values
        // VEGF
        double BM0 = POPBM0*exp(ETABM0);
        double MRT = POPMRT*exp(ETAMRT23);
        double IMAX = POPIMAX;
        double IC50 = POPIC50*exp(ETAIC50);
        double HILL = POPHILL;
        double DPSLO = POPDPSLO*exp(ETADPSLO);

        // VEGFR2
        double BM02 = POPBM02*exp(ETABM02);
        double MRT2 = POPMRT2*exp(ETAMRT23);
        double IMAX2 = POPIMAX;
        double IC502 = POPIC50*exp(ETAIC502);
        double HILL2 = POPHILL2;

        // VEGFR3
        double BM03 = POPBM03*exp(ETABM03);
        double MRT3 = POPMRT3*exp(ETAMRT23);
        double IMAX3 = POPIMAX;
        double IC503 = POPIC50*exp(ETAIC503);

        // SKIT
        double BM0S = POPBM0S*exp(ETABM0S);
        double MRTS = POPMRTS*exp(ETAMRTS);
        double IMAXS = POPIMAX;
        double IC50S = POPIC50*exp(ETAIC50S);
        double DPSLOS = POPDPSLO*exp(ETADPSLOS);

        // Compartment initial conditions
        VEGF_0 = BM0;
        VEGFR2_0 = BM02;
        VEGFR3_0 = BM03;
        SKIT_0 = BM0S;

        // Rate constants
        double KOUT = 1/MRT;
        double KOUT2 = 1/MRT2;
        double KOUT3 = 1/MRT3;
        double KOUTS = 1/MRTS;
        double KIN2 = BM02*KOUT2;
        double KIN3 = BM03*KOUT3;

        // Sunitinib exposure
        double AUC = DOSE/CL;

        // Linear disease progression
        double DP = BM0*(1+DPSLO*TIME);
        double DPS = BM0S*(1+DPSLOS*TIME);

$ODE		// Differential equations
        // Pharmacodynamics - sunitinib effect
        double KIN = DP*KOUT;
        double KINS = DPS*KOUTS;

        double EFF = IMAX*pow(AUC,HILL)/(pow(IC50,HILL)+pow(AUC,HILL));
        double EFF2 = IMAX2*pow(AUC,HILL2)/(pow(IC502,HILL2)+pow(AUC,HILL2));
        double EFF3 = IMAX3*AUC/(IC503+AUC);
        double EFFS = IMAXS*AUC/(IC50S+AUC);

        dxdt_VEGF = KIN-KOUT*(1-EFF)*VEGF;
        dxdt_VEGFR2 = KIN2*(1-EFF2)-KOUT2*VEGFR2;
        dxdt_VEGFR3 = KIN3*(1-EFF3)-KOUT3*VEGFR3;
        dxdt_SKIT = KINS*(1-EFFS)-KOUTS*SKIT;

$TABLE 	double IPRE_VEGF = log(VEGF);
        double IPRE_VEGFR2 = log(VEGFR2);
        double IPRE_VEGFR3 = log(VEGFR3);
        double IPRE_SKIT = log(SKIT);

$CAPTURE DOSE CL IPRE_VEGF IPRE_VEGFR2 IPRE_VEGFR3 IPRE_SKIT
'
# Compile the model code
  mod <- mcode("popVEGF",code)

# ------------------------------------------------------------------------------
# Pull out variance parameters from model file
  OMEGA <- as.matrix(omat(mod))
