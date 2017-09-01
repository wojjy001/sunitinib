#################################
### Directories to be changed ###
#################################
1. Open "run_pk.R"
2. Line 8; set "file.dir" to the directory where the "SimTDM_byPK" folder is
saved (not including "SimTDM_byPK" in the name)
3. Line 16; set "output.dir" to the same directory.  All simulated output will
be saved into the folder where the .R files are saved.
4. Line 37; "study.name" can only be one of the 6 studies that I have defined:
  - nodrug
  - standard
  - target_trough_exact
  - target_trough_round
  - target_auc_exact
  - target_auc_round
Each of these "studies" have their own .R script in "SimTDM_byPK".  Changing
"study.name" allows "run_pk.R" to source the correct study script
5. Open "run_pd.R"
6. Line 12; change "output.dir" to the same as Step 3
7. Line 13; change "study.name" to the same as Step 4
8. Line 21; change "file.dir" to the same as Step 2

###########################
### Simulation Sequence ###
###########################
1. Need to run "run_pk.R" first
2. Run "run_pd.R" second. You can run the PD simulation at a separate instance
if you had previously run and saved the PK simulation
