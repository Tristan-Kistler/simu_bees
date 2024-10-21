# .libPaths("/bao/lib_R3.6.2") # Server
.libPaths("/bao/lib_R3.5.0") # Server
libpath=paste("/opt/intel/compilers_and_libraries_2016.1.150/linux/compiler/lib/intel64", Sys.getenv("LD_LIBRARY_PATH"), sep=":") # Server
Sys.setenv(LD_LIBRARY_PATH=libpath) # Server
rm(list=ls()) # Cleaning the environment


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
t1 <- proc.time()#variable useful to get running time


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
arg_shell <- commandArgs(trailingOnly = TRUE)
print(arg_shell)
WD_PATH_SCRIPTS <- arg_shell[1] # folder where the scripts are stored
WD_PATH_RUN <- arg_shell[2] # user should indicate in WD_PATH_RUN the path to the folder that will serve as a working directory for Rstudio
WD_PATH_Group_id <- arg_shell[3] # folder for results of repetited runs of a same scenario
Rep <- arg_shell[4]
Rep <- as.numeric(Rep)

setwd(WD_PATH_RUN)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#charging packages
library(tidyr)
library(dplyr)
options(dplyr.summarise.inform=F)
library(purrr) # For map2 function
library(data.table)
library(mvtnorm)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
source(paste0(WD_PATH_Group_id, "/par.R"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
#files necessary for Blup evaluation by BLUPF90 will be stored in specific subfolders
if(SEL_BLUP){
  WD_PATH_BLUP <- paste0(WD_PATH_RUN, "/Data_BLUP") #we create a path to a Output subfolder
  
  if(dir.exists(file.path(WD_PATH_BLUP))){#if it exist, nothing is done
    
  }else{
    dir.create(WD_PATH_BLUP)#else we create the subfolder
  }
}


## ---- eval=T------------------------------------------------------------------------------------------------------------------------------------------------------
if (FIXED_SEED != 0){
  set.seed(FIXED_SEED)
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Will be overwritten by sourcing the environment of burnin
# source(paste0(WD_PATH_SCRIPTS, "/breeding_functions.R"))
# source(paste0(WD_PATH_SCRIPTS, "/matrix_functions.R"))
# source(paste0(WD_PATH_SCRIPTS, "/variable_functions.R"))
# source(paste0(WD_PATH_SCRIPTS, "/output_functions.R"))

# To run Pim's programs, we wills source them.
# To have Pim's programs running in a separate environment, the trick I found was to wrap inside a function Pim's program.
# When sourcing the file with Pim's program inside a function, the master_script (simulation script) learn the function containing Pim's program.
# This can then be called by the simulation program to run Pim's program.
# Because Pim's program works inside that functions, all variables it creates are private to the function, so they don't pollute the simulation environment.
if (SEL_BLUP){if(pedigree_Pim){source(paste0(PATH_PIM_PEDIGREE_SCRIPT), chdir = T, echo = F, local=F)}}
if (SEL_BLUP){if(pedigree_Pim){source(paste0(PATH_PIM_AINV_SCRIPT), chdir = T, echo = F, local=F)}}
## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# Loading burn-in environment
load(paste0("Burn_in_", LAST_YEAR_BURN_IN, ".RData"))

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
year <- year+1


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
if(Q_AND_TQ_SAME_APIARIES == TRUE){
  if(KEEPER_CONNEXION == "MINIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
    Q_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
  }
  if(KEEPER_CONNEXION == "MAXIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
    Q_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
  }
  if(KEEPER_CONNEXION == "ANY"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
    Q_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=Q_n, indiv_type="Q", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_KEEPERS, Nb_sister_groups_per_keeper=NB_SISTER_GROUPS_PER_KEEPER, Nb_keepers_per_sister_group=NB_KEEPERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_SISTERS_PER_GROUP_PER_KEEPER, NB_BQ, Performance_year=year)
  }
}else{ # If Qs and TQs are on different apiaries we distribute Qs to NB_BREEDERS
  if(BREEDER_CONNEXION == "MINIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
    Q_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_BREEDERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
  }
  if(BREEDER_CONNEXION == "MAXIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
    Q_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_BREEDERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
  }
  if(BREEDER_CONNEXION == "ANY"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
    Q_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=Q_n, indiv_type="Q", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_BREEDERS, Nb_sister_groups_per_keeper=NB_Q_SISTER_GROUPS_PER_BREEDER, Nb_keepers_per_sister_group=NB_BREEDERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_Q_SISTERS_PER_GROUP_PER_BREEDER, NB_BQ, Performance_year=year)
  }
}
#We transpose this distribution to actually survived queens:
WQ_n$Beekeeper_year <- Q_n$Beekeeper_year[Q_W1_n]


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Nb_keepers should be NB_BREEDERS if Qs and TQs are on separate apiaries
if(Q_AND_TQ_SAME_APIARIES == TRUE){
  Beekeepers_fixed_effects_n <- Get_beekeepers_fixed_effects(Nb_keepers=NB_KEEPERS, indiv_type="Q", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES)
}else{
  # If Qs and TQs are or different apiaries:
  Beekeepers_fixed_effects_n <- Get_beekeepers_fixed_effects(Nb_keepers=NB_BREEDERS, indiv_type="Q", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES)
}

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("Beekeepers_fixed_effects_", year)
  assign(Name_data_frame, Beekeepers_fixed_effects_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
birth_year_WQ <- year-1
WQ_n <- Add_colony_Performance(WQ_or_WTQ_n=WQ_n, indiv_type="Q", Beekeepers_fixed_effects_n, VAR_E)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", birth_year_WQ)
  assign(Name_data_frame,WQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Output_Q <- Update_Output_Q_from_WQ_n(WQ_n, Output_Q, pedigree, year, Rep)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(Q_AND_TQ_SAME_APIARIES == TRUE){
  if(KEEPER_CONNEXION == "MINIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
    TQ_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
  }
  if(KEEPER_CONNEXION == "MAXIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
    TQ_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
  }
  if(KEEPER_CONNEXION == "ANY"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
    TQ_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=TQ_n, indiv_type="TQ", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_KEEPERS, Nb_sister_groups_per_keeper=NB_SISTER_GROUPS_PER_KEEPER, Nb_keepers_per_sister_group=NB_KEEPERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_SISTERS_PER_GROUP_PER_KEEPER, NB_BQ, Performance_year=year)
  }
}else{ # If Qs and TQs are on different apiaries we distribute Qs to NB_TESTERS
  if(TESTER_CONNEXION == "MINIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
    TQ_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_TESTERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
  }
  if(TESTER_CONNEXION == "MAXIMAL"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
    TQ_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_TESTERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
  }
  if(TESTER_CONNEXION == "ANY"){
    # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
    TQ_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=TQ_n, indiv_type="TQ", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_TESTERS, Nb_sister_groups_per_keeper=NB_TQ_SISTER_GROUPS_PER_TESTER, Nb_keepers_per_sister_group=NB_TESTERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_TQ_SISTERS_PER_GROUP_PER_TESTER, NB_BQ, Performance_year=year)
  }
}
#We transpose this distribution to actually survived queens:
WTQ_n$Beekeeper_year <- TQ_n$Beekeeper_year[TQ_W1_n]


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# If TQs are on different apiaries than Qs, we use the effect of Testers to fill the Beekeepers_fixed_effects_n table, and draw new effects for the beekeeper_x_year effect
if(Q_AND_TQ_SAME_APIARIES == FALSE){
  # If Qs and TQs are or different apiaries:
  Beekeepers_fixed_effects_n <- Get_beekeepers_fixed_effects(Nb_keepers=NB_TESTERS, indiv_type="TQ", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES)
}

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("Beekeepers_fixed_effects_", year)
  assign(Name_data_frame, Beekeepers_fixed_effects_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WTQ_n <- Add_colony_Performance(WQ_or_WTQ_n=WTQ_n, indiv_type="TQ", Beekeepers_fixed_effects_n, VAR_E)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WTQ_", birth_year_WQ)
  assign(Name_data_frame,WTQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Output_TQ <- Update_Output_TQ(WTQ_n, Output_TQ, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_or_TQ_and_Q_W2_n <- Winter_loss_Q(Q_or_TQ_n=WQ_n, NB_W_Q_or_TQ=NB_W2_Q)

W2_Q_n <- Q_or_TQ_and_Q_W2_n[[1]]#WQ stands for Wintered Queens
Q_W2_n <- Q_or_TQ_and_Q_W2_n[[2]]#contains part of the IDs of surviving queens

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("W2_Q_", (year-1))
  assign(Name_data_frame, W2_Q_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_or_TQ_and_Q_W2_n <- Winter_loss_Q(Q_or_TQ_n=WTQ_n, NB_W_Q_or_TQ=NB_W2_TQ)

W2_TQ_n <- Q_or_TQ_and_Q_W2_n[[1]]#WQ stands for Wintered Queens
TQ_W2_n <- Q_or_TQ_and_Q_W2_n[[2]]#contains part of the IDs of surviving queens

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("W2_TQ_", (year-1))
  assign(Name_data_frame,W2_TQ_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
year <- year + 1


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(SEL_BLUP){
  perf_Bras <- Create_Perf_Bras(Output_Q, Output_TQ, erase_Perf_T1_Q=ERASE_PERF_T1_Q_BASE, erase_Perf_T2_Q=ERASE_PERF_T2_Q_BASE, erase_Perf_T1_TQ=ERASE_PERF_T1_TQ_BASE, erase_Perf_T2_TQ=ERASE_PERF_T2_TQ_BASE,
   year_erase_Perf_T1_Q=YEAR_ERASE_PERF_T1_Q_BASE, year_erase_Perf_T2_Q=YEAR_ERASE_PERF_T2_Q_BASE, year_erase_Perf_T1_TQ=YEAR_ERASE_PERF_T1_TQ_BASE, year_erase_Perf_T2_TQ=YEAR_ERASE_PERF_T2_TQ_BASE)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(SEL_BLUP){
  if(pedigree_Pim){ # if an external call to pedigree.R will be done
    #ped_blup
    fwrite(pedigree, file=paste0(WD_PATH_BLUP, "/input-pedigree.txt"), sep = " ", dec = ".", row.names = F, col.names = F)
    
    # pedigree.R also needs the parameters NB_D and NB_Q_PER_PS to calculate BLUP. We wright these in one dq.txt file to pass these parameters
    steer <- c(NS_Q_open, NB_D_FECUNDATING, 0, 
               0, 0, 0, 0, 0, 0, 0, 0, 0) # only three first cols to be filled, rest is to have correct dimensions
    fwrite(as.list(steer), paste0(WD_PATH_BLUP, "/steer.txt"), col.names = FALSE, row.names = FALSE, sep=" ")
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(SEL_BLUP){
  if(pedigree_Pim){ # if an external call to pedigree.R is done
    Pim_pedigree(WD_PATH_BLUP)
  }else{
    write_AINV_input(pedigree, WD_PATH_BLUP)
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(SEL_BLUP){
  if(pedigree_Pim){ # if an external call to pedigree.R is done
    Pim_AINV(WD_PATH_BLUP)
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
perf_Bras_joinID <- perf_Bras_RenumIDsAINV(WD_PATH_BLUP, perf_Bras)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Writing perf file with recoded IDs
fwrite(perf_Bras_joinID, paste0(WD_PATH_BLUP, "/perf_Bras.txt"), col.names = F, sep = " ")


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## VE1_ini <- 1
## VA1_D_ini <- 10
## VA1_M_ini <- 1
## VE2_ini <- 1
## VA2_D_ini <- 1000
## VA2_M_ini <- 1


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(SEL_BLUP){
  system(command = paste("cd", WD_PATH_BLUP, ";ulimit -s unlimited;sh",PATH_AIREML_SCRIPT, WD_PATH_BLUP, SINGLE_TRAIT, VA1_D_ini, VA1_M_ini, VA2_D_ini, VA2_M_ini, VE1_ini, VE2_ini, MODEL_TYPE), ignore.stdout = TRUE, ignore.stderr = TRUE)
}
# VAR_A


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(SEL_BLUP){
  CandidateBQ_n_DPQ_n <- Retrieve_BLUP_sol(WD_PATH_BLUP, perf_Bras_joinID, SINGLE_TRAIT, CandidateBQ_n=W2_Q_n, CandidateDPQ_n=WTQ_n)
  
  W2_Q_n <- CandidateBQ_n_DPQ_n[[1]]
  WTQ_n <- CandidateBQ_n_DPQ_n[[2]] #first winter surviving TQs
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# if the estimation failed to reach the convergence criterion, we stop it
if(SEL_BLUP){
  if(STOP_IF_CONV_FAILED){
    conv_crit_failed <- unique(as.numeric(W2_Q_n$Conv_crit))
    round_failed <- unique(as.numeric(W2_Q_n$Conv_round))
    if(conv_crit_failed > 1e-11 | conv_crit_failed == 0){
      write(
        paste0("BLUP estimation did not reach the convergence criterion in scenario Base at year ",
               year, ",\nat round ", round_failed, " reaching a convergence criterion of ", conv_crit_failed),
        file=paste0(WD_PATH_RUN, "/BLUP_not_converged.txt"))
      stop("BLUP estimation did not reach the convergence criterion, conv_crit was: ", conv_crit_failed, " at year: ", year)
    }
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
if(SEL_BLUP){
  # We also report EBVs on second wintered WTQ_n (family sel is done before first mortality, but DPQ sel is done after)
  W2_TQ_n <- Join_EBVs_Q_n(Table_to_update_n=W2_TQ_n, CandidateDPQ_n=WTQ_n, SINGLE_TRAIT)
  
  # INITIALIZING: AS OUTPUT_Q DOES NOT HAVE EBVs Cols yet, we do a join without out selecting these columns as done in function (Join_EBVs_to_Output) in later years
  # We uptdate Output_Q and Output_TQ with the EBVs as predicted before selection
  if(SINGLE_TRAIT == FALSE){
    Output_Q <- Output_Q %>% left_join(W2_Q_n %>% select(ID,
                                                         EBV_d1_Q, EBV_m1_Q, EBV_d2_Q, EBV_m2_Q,
                                                         EBV_d1_W, EBV_m1_W, EBV_d2_W, EBV_m2_W,
                                                         SE_EBV_d1_Q, SE_EBV_m1_Q, SE_EBV_d2_Q, SE_EBV_m2_Q,
                                                         SE_EBV_d1_W, SE_EBV_m1_W, SE_EBV_d2_W, SE_EBV_m2_W,
                                                         E_Beekeeper_x_year_effect_T1, E_Beekeeper_x_year_effect_T2,
                                                         SE_Beekeeper_x_year_effect_T1, SE_Beekeeper_x_year_effect_T2,
                                                         Conv_crit, Conv_round,
                                                         E_VAR_A1_D, E_VAR_A1_M, E_VAR_A2_D, E_VAR_A2_M, 
                                                         E_COV_A1_D_M, E_COV_A2_D_M,
                                                         E_COV_A1_A2_D, E_COV_A1_A2_M,
                                                         E_COV_A1_D_A2_M, E_COV_A1_M_A2_D,
                                                         E_VAR_E1, E_VAR_E2, E_COV_E1_E2),
                                       by="ID")
    Output_TQ <- Output_TQ %>% left_join(WTQ_n %>% select(ID,
                                                          EBV_d1_Q, EBV_m1_Q, EBV_d2_Q, EBV_m2_Q,
                                                          EBV_d1_W, EBV_m1_W, EBV_d2_W, EBV_m2_W,
                                                          SE_EBV_d1_Q, SE_EBV_m1_Q, SE_EBV_d2_Q, SE_EBV_m2_Q,
                                                          SE_EBV_d1_W, SE_EBV_m1_W, SE_EBV_d2_W, SE_EBV_m2_W,
                                                          E_Beekeeper_x_year_effect_T1, E_Beekeeper_x_year_effect_T2,
                                                          SE_Beekeeper_x_year_effect_T1, SE_Beekeeper_x_year_effect_T2,
                                                          Conv_crit, Conv_round,
                                                          E_VAR_A1_D, E_VAR_A1_M, E_VAR_A2_D, E_VAR_A2_M, 
                                                          E_COV_A1_D_M, E_COV_A2_D_M,
                                                          E_COV_A1_A2_D, E_COV_A1_A2_M,
                                                          E_COV_A1_D_A2_M, E_COV_A1_M_A2_D,
                                                          E_VAR_E1, E_VAR_E2, E_COV_E1_E2),
                                         by="ID")
    
  }else{ # if single trait sim
    
    Output_Q <- Output_Q %>% left_join(W2_Q_n %>% select(ID,
                                                         EBV_d1_Q, EBV_m1_Q,
                                                         EBV_d1_W, EBV_m1_W,
                                                         SE_EBV_d1_Q, SE_EBV_m1_Q,
                                                         SE_EBV_d1_W, SE_EBV_m1_W,
                                                         E_Beekeeper_x_year_effect_T1,
                                                         SE_Beekeeper_x_year_effect_T1,
                                                         Conv_crit, Conv_round,
                                                         E_VAR_A1_D, E_VAR_A1_M, E_COV_A1_D_M, 
                                                         E_VAR_E1),
                                       by="ID")
    Output_TQ <- Output_TQ %>% left_join(WTQ_n %>% select(ID,
                                                          EBV_d1_Q, EBV_m1_Q,
                                                          EBV_d1_W, EBV_m1_W,
                                                          SE_EBV_d1_Q, SE_EBV_m1_Q,
                                                          SE_EBV_d1_W, SE_EBV_m1_W,
                                                          E_Beekeeper_x_year_effect_T1,
                                                          SE_Beekeeper_x_year_effect_T1,
                                                          Conv_crit, Conv_round,
                                                          E_VAR_A1_D, E_VAR_A1_M, E_COV_A1_D_M, 
                                                          E_VAR_E1),
                                         by="ID")
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
BQ_n <- Select_BQ(W2_Q_n, NB_BQ, Q_SELECTION_CRITERION)

birth_year_BQ <- (year-1)
if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("BQ_", birth_year_BQ)
  assign(Name_data_frame,BQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WQ_n_save_one <- WQ_n #Q_n_save_one will store queend of the previous year


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_n <- Create_indiv_n_blanck(BQ_n, indiv_type="Q", Nb_indiv=NB_Q, Nb_indiv_per_BQ=NB_Q_PER_BQ, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_n <- Fill_indiv_n(indiv_n=Q_n, indiv_type="Q", D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ=NB_Q_PER_BQ, NB_BQ, year)


if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("Q_", year)
  assign(Name_data_frame,Q_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=Q_n, NB_W_Q_or_TQ=NB_W1_Q)

WQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
Q_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains part of the IDs of surviving queens


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, NB_D_FECUNDATING)
}else{
  pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(NB_W1_Q+NB_W1_TQ+NB_D) )


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=(NB_W1_Q+NB_W1_TQ+NB_D), NB_W1_indiv=NB_W1_Q)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WQ_n, A, ped_A)  

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", year)
  assign(Name_data_frame,WQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# WTQ_n_save_two <- WTQ_n_save_one
# WTQ_n_save_one <- WTQ_n

WTQ_n_save_one <- WTQ_n


## --------------------------------------------------------------------------------------------------------------------------------------------------------
TQ_n <- Create_indiv_n_blanck(BQ_n, indiv_type="TQ", Nb_indiv=NB_TQ, Nb_indiv_per_BQ=NB_TQ_PER_BQ, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
TQ_n <- Fill_indiv_n(indiv_n=TQ_n, indiv_type="TQ", D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ=NB_TQ_PER_BQ, NB_BQ, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=TQ_n, NB_W_Q_or_TQ=NB_W1_TQ)

WTQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
TQ_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains part of the IDs of surviving queens


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# We add the unknown mates of (open mated) TQs when UNKNOWN_1b == 'O'
if(pedigree_Pim == F){
  if(UNKNOWN_1b == '0'){ # then we add dummy open mating sires for open mated queens
    pedigree <- Add_unknown_PGM_mated_to_pedigree(pedigree, indiv_n=WTQ_n, NS_Q_or_TQ_open=NS_TQ, NB_D_FECUNDATING) 
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, NB_D_FECUNDATING)
}else{
  pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WTQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(2*NB_W1_Q+NB_W1_TQ+NB_D))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=0, NB_W1_indiv=NB_W1_TQ)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WTQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WTQ_n, A, ped_A) # Saving with yearly name is done after second overwintering. Could be changed


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==TRUE){
  Q_of_PS_n <- Select_PS_from_WTQ(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_BQ, Criterion=PGM_SELECTION_CRITERION)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("Q_of_PS_", year-2)
    assign(Name_data_frame, Q_of_PS_n)
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==TRUE){
  DPQ_n <- Create_DPQ_n_blanck(W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==TRUE){
  DPQ_n <- Fill_DPQ_n(Q_of_PS_n, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS, DPQ_n, Criterion=DPQ_SELECTION_CRITERION)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("DPQ_", year-2)
    assign(Name_data_frame,DPQ_n)
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==FALSE){
  DPQ_n <- Select_DPQs(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_TQ, NB_DPQ, Criterion=DPQ_SELECTION_CRITERION)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("DPQ_", year-2)
    assign(Name_data_frame, DPQ_n)
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
PS_or_TQ_to_WQ_n_and_WQ_n <- Assign_PGM_or_TQ_to_WQ(DPQ_n, WQ_n)

PS_or_TQ_to_WQ_n <- PS_or_TQ_to_WQ_n_and_WQ_n[[1]]
WQ_n <- PS_or_TQ_to_WQ_n_and_WQ_n[[2]]

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", year)
  assign(Name_data_frame, WQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- Add_PGM_mated_to_pedigree_Pim(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE)
}else{
  pedigree <- Add_PGM_mated_to_pedigree(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim == FALSE){
  pedigree <- order_off_give_seq2(pedigree, birth_year_Q_offspring=year, NB_Wg_offspring=NB_W1_Q+NB_W1_TQ, ID_1b=UNKNOWN_1b)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
D_n <- Produce_and_assign_D(Q_to_mate=WQ_n, DPQ_n, NB_D, NB_D_FECUNDATING, VAR_A_meiose, year)

if(MATING_SCALE=="Q_x_PS"){
  WQ_n <- select(WQ_n, -ID_PGM)#We have to withdraw column "ID_PGM" because it is not present in WQ_n when MATING_SCALE==Q_x_DPQ. Instead, there will be a column ID_TQ. This difference then creates a table WQ_n with a column ID_PGM which gets transfered to Output_Q_n in the process of updating Output_Q with function Update_Output_Q_from_WQ_n, which then creates multiple columns of ID_PGM... in table Output_Q. By erasing this ID_PGM column from WQ_n just after it was used by Produce_and_assign_D we keep all further functions working in both MATING_SCALE parameters
}

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("D_", year)
  assign(Name_data_frame,D_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_D_n(WTQ_with_DPQs=WTQ_n_save_one, D_n,
                          A,
                          DPQ_index=NB_W1_Q, 
                          NB_D)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_D_n(A, ped_A, NB_D)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WQ_n <- Produce_Wgroup(W_indiv_n=WQ_n, Ds_n=D_n, indiv_type="Q")

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", year)
  assign(Name_data_frame,WQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
TD_n <- Produce_indiv_randomly(Nb_indiv=NB_TD, indiv_type="TD", Q_to_mate=WTQ_n, Var_A=(0.5*VAR_A), year)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("TD_", year)
  assign(Name_data_frame,TD_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WTQ_n <- Produce_Wgroup(W_indiv_n=WTQ_n, Ds_n=TD_n, indiv_type="TQ")#WTQ_n should be the table of TQ born in current year that survived winter


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- Add_unknown_PGM_mated_to_pedigree_Pim(pedigree, indiv_n=WTQ_n, ID_1b=UNKNOWN_1b)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
for (year in seq(from=year+1, to=(LAST_YEAR-1), by=2)){ # LAST_YEAR-1 because for each loop here, 2 consecutive years are simulated
  tn_start <- proc.time()
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(Q_AND_TQ_SAME_APIARIES == TRUE){
    if(KEEPER_CONNEXION == "MINIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
      Q_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
    }
    if(KEEPER_CONNEXION == "MAXIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
      Q_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
    }
    if(KEEPER_CONNEXION == "ANY"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
      Q_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=Q_n, indiv_type="Q", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_KEEPERS, Nb_sister_groups_per_keeper=NB_SISTER_GROUPS_PER_KEEPER, Nb_keepers_per_sister_group=NB_KEEPERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_SISTERS_PER_GROUP_PER_KEEPER, NB_BQ, Performance_year=year)
    }
  }else{ # If Qs and TQs are on different apiaries we distribute Qs to NB_BREEDERS
    if(BREEDER_CONNEXION == "MINIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
      Q_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_BREEDERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
    }
    if(BREEDER_CONNEXION == "MAXIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
      Q_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=Q_n, Nb_keepers=NB_BREEDERS, Indiv_per_line=NB_Q_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year)
    }
    if(BREEDER_CONNEXION == "ANY"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_BREEDERS
      Q_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=Q_n, indiv_type="Q", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_BREEDERS, Nb_sister_groups_per_keeper=NB_Q_SISTER_GROUPS_PER_BREEDER, Nb_keepers_per_sister_group=NB_BREEDERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_Q_SISTERS_PER_GROUP_PER_BREEDER, NB_BQ, Performance_year=year)
    }
  }
  #We transpose this distribution to actually survived queens:
  WQ_n$Beekeeper_year <- Q_n$Beekeeper_year[Q_W1_n]
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Nb_keepers should be NB_BREEDERS if Qs and TQs are on separate apiaries
  if(Q_AND_TQ_SAME_APIARIES == TRUE){
    Beekeepers_fixed_effects_n <- Get_beekeepers_fixed_effects(Nb_keepers=NB_KEEPERS, indiv_type="Q", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES)
  }else{
    # If Qs and TQs are or different apiaries:
    Beekeepers_fixed_effects_n <- Get_beekeepers_fixed_effects(Nb_keepers=NB_BREEDERS, indiv_type="Q", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES)
  }
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("Beekeepers_fixed_effects_", year)
    assign(Name_data_frame, Beekeepers_fixed_effects_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  birth_year_WQ <- year-1
  WQ_n <- Add_colony_Performance(WQ_or_WTQ_n=WQ_n, indiv_type="Q", Beekeepers_fixed_effects_n, VAR_E)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", birth_year_WQ)
    assign(Name_data_frame,WQ_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Output_Q <- Update_Output_Q_from_WQ_n(WQ_n, Output_Q, pedigree, year, Rep)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(Q_AND_TQ_SAME_APIARIES == TRUE){
    if(KEEPER_CONNEXION == "MINIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
      TQ_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
    }
    if(KEEPER_CONNEXION == "MAXIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
      TQ_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
    }
    if(KEEPER_CONNEXION == "ANY"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
      TQ_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=TQ_n, indiv_type="TQ", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_KEEPERS, Nb_sister_groups_per_keeper=NB_SISTER_GROUPS_PER_KEEPER, Nb_keepers_per_sister_group=NB_KEEPERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_SISTERS_PER_GROUP_PER_KEEPER, NB_BQ, Performance_year=year)
    }
  }else{ # If Qs and TQs are on different apiaries we distribute Qs to NB_TESTERS
    if(TESTER_CONNEXION == "MINIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
      TQ_n <- Assign_indiv_keepers_min_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_TESTERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
    }
    if(TESTER_CONNEXION == "MAXIMAL"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
      TQ_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=TQ_n, Nb_keepers=NB_TESTERS, Indiv_per_line=NB_TQ_PER_BQ, Q_AND_TQ_SAME_APIARIES, indiv_type="TQ", Performance_year=year)
    }
    if(TESTER_CONNEXION == "ANY"){
      # If Qs should be on same apiaries than TQs, Nb_keepers=NB_KEEPERS, otherwise Nb_keepers=NB_TESTERS
      TQ_n <- Assign_indiv_keepers_any_connection(Q_or_TQ_n=TQ_n, indiv_type="TQ", Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_TESTERS, Nb_sister_groups_per_keeper=NB_TQ_SISTER_GROUPS_PER_TESTER, Nb_keepers_per_sister_group=NB_TESTERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_TQ_SISTERS_PER_GROUP_PER_TESTER, NB_BQ, Performance_year=year)
    }
  }
  #We transpose this distribution to actually survived queens:
  WTQ_n$Beekeeper_year <- TQ_n$Beekeeper_year[TQ_W1_n]
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # If TQs are on different apiaries than Qs, we use the effect of Testers to fill the Beekeepers_fixed_effects_n table, and draw new effects for the beekeeper_x_year effect
  if(Q_AND_TQ_SAME_APIARIES == FALSE){
    # If Qs and TQs are or different apiaries:
    Beekeepers_fixed_effects_n <- Get_beekeepers_fixed_effects(Nb_keepers=NB_TESTERS, indiv_type="TQ", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES)
  }
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("Beekeepers_fixed_effects_", year)
    assign(Name_data_frame, Beekeepers_fixed_effects_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  WTQ_n <- Add_colony_Performance(WQ_or_WTQ_n=WTQ_n, indiv_type="TQ", Beekeepers_fixed_effects_n, VAR_E)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WTQ_", birth_year_WQ)
    assign(Name_data_frame,WTQ_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Output_TQ <- Update_Output_TQ(WTQ_n, Output_TQ, year)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W2_n <- Winter_loss_Q(Q_or_TQ_n=WQ_n, NB_W_Q_or_TQ=NB_W2_Q)
  
  W2_Q_n <- Q_or_TQ_and_Q_W2_n[[1]]#WQ stands for Wintered Queens
  Q_W2_n <- Q_or_TQ_and_Q_W2_n[[2]]#contains part of the IDs of surviving queens
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("W2_Q_", (year-1))
    assign(Name_data_frame, W2_Q_n)# The tbl is named by marking the birth year of queens
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W2_n <- Winter_loss_Q(Q_or_TQ_n=WTQ_n, NB_W_Q_or_TQ=NB_W2_TQ)
  
  W2_TQ_n <- Q_or_TQ_and_Q_W2_n[[1]]#WQ stands for Wintered Queens
  TQ_W2_n <- Q_or_TQ_and_Q_W2_n[[2]]#contains part of the IDs of surviving queens
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("W2_TQ_", (year-1))
    assign(Name_data_frame,W2_TQ_n)# The tbl is named by marking the birth year of queens
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  year <- year + 1
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(SEL_BLUP){
    perf_Bras <- Create_Perf_Bras(Output_Q, Output_TQ, erase_Perf_T1_Q=ERASE_PERF_T1_Q_BASE, erase_Perf_T2_Q=ERASE_PERF_T2_Q_BASE, erase_Perf_T1_TQ=ERASE_PERF_T1_TQ_BASE, erase_Perf_T2_TQ=ERASE_PERF_T2_TQ_BASE,
     year_erase_Perf_T1_Q=YEAR_ERASE_PERF_T1_Q_BASE, year_erase_Perf_T2_Q=YEAR_ERASE_PERF_T2_Q_BASE, year_erase_Perf_T1_TQ=YEAR_ERASE_PERF_T1_TQ_BASE, year_erase_Perf_T2_TQ=YEAR_ERASE_PERF_T2_TQ_BASE)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(SEL_BLUP){
    if(pedigree_Pim){ # if an external call to pedigree.R will be done
    #ped_blup
    fwrite(pedigree, file=paste0(WD_PATH_BLUP, "/input-pedigree.txt"), sep = " ", dec = ".", row.names = F, col.names = F)
    }
  }
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(SEL_BLUP){
    if(pedigree_Pim){ # if an external call to pedigree.R is done
      Pim_pedigree(WD_PATH_BLUP)
    }else{
      write_AINV_input(pedigree, WD_PATH_BLUP)
    }
  }
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------
  if(SEL_BLUP){
    if(pedigree_Pim){ # if an external call to pedigree.R is done
      Pim_AINV(WD_PATH_BLUP)
    }
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  perf_Bras_joinID <- perf_Bras_RenumIDsAINV(WD_PATH_BLUP, perf_Bras)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Writing perf file with recoded IDs
  fwrite(perf_Bras_joinID, paste0(WD_PATH_BLUP, "/perf_Bras.txt"), col.names = F, sep = " ")
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(SEL_BLUP){
    system(command = paste("cd", WD_PATH_BLUP, ";ulimit -s unlimited;sh",PATH_AIREML_SCRIPT, WD_PATH_BLUP, SINGLE_TRAIT, VA1_D_ini, VA1_M_ini, VA2_D_ini, VA2_M_ini, VE1_ini, VE2_ini, MODEL_TYPE), ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(SEL_BLUP){
    CandidateBQ_n_DPQ_n <- Retrieve_BLUP_sol(WD_PATH_BLUP, perf_Bras_joinID, SINGLE_TRAIT, CandidateBQ_n=W2_Q_n, CandidateDPQ_n=WTQ_n)
    
    W2_Q_n <- CandidateBQ_n_DPQ_n[[1]]
    WTQ_n <- CandidateBQ_n_DPQ_n[[2]]
  }
  
  if(SEL_BLUP){
    if(STOP_IF_CONV_FAILED){
      conv_crit_failed <- unique(as.numeric(W2_Q_n$Conv_crit))
      round_failed <- unique(as.numeric(W2_Q_n$Conv_round))
      if(conv_crit_failed > 1e-11 | conv_crit_failed == 0){
        write(
          paste0("BLUP estimation did not reach the convergence criterion in scenario Base at year ",
                 year, ",\nat round ", round_failed, " reaching a convergence criterion of ", conv_crit_failed),
          file=paste0(WD_PATH_RUN, "/BLUP_not_converged.txt"))
        stop("BLUP estimation did not reach the convergence criterion, conv_crit was: ", conv_crit_failed, " at year: ", year)
      }
    }
  }
  
  if(SEL_BLUP){
    # We also report EBVs on second wintered WTQ_n (family sel is done before first mortality, but DPQ sel is done after)
    W2_TQ_n <- Join_EBVs_Q_n(Table_to_update_n=W2_TQ_n, CandidateDPQ_n=WTQ_n, SINGLE_TRAIT)
    # We uptdate Output_Q and Output_TQ witht the EBVs as predicted before selection
    Output_Q <- Join_EBVs_to_Output(Output_Q_or_TQ=Output_Q, Candidates_n=W2_Q_n, SINGLE_TRAIT)
    
    Output_TQ <- Join_EBVs_to_Output(Output_Q_or_TQ=Output_TQ, Candidates_n=WTQ_n, SINGLE_TRAIT)
  }
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  BQ_n <- Select_BQ(W2_Q_n, NB_BQ, Q_SELECTION_CRITERION)
  
  birth_year_BQ <- (year-1)
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("BQ_", birth_year_BQ)
    assign(Name_data_frame,BQ_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  WQ_n_save_one <- WQ_n #Q_n_save_one will store queend of the previous year
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Q_n <- Create_indiv_n_blanck(BQ_n, indiv_type="Q", Nb_indiv=NB_Q, Nb_indiv_per_BQ=NB_Q_PER_BQ, year)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Q_n <- Fill_indiv_n(indiv_n=Q_n, indiv_type="Q", D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ=NB_Q_PER_BQ, NB_BQ, year)
  
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("Q_", year)
    assign(Name_data_frame,Q_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=Q_n, NB_W_Q_or_TQ=NB_W1_Q)
  
  WQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
  Q_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains part of the IDs of surviving queens
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(pedigree_Pim){
    pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, NB_D_FECUNDATING)
  }else{
    pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ped_A <- Create_ped_A_Q_n(W_indiv_n=WQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(NB_W1_Q+NB_W1_TQ+NB_D) )
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=(NB_W1_Q+NB_W1_TQ+NB_D), NB_W1_indiv=NB_W1_Q)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  WQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WQ_n, A, ped_A)  
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", year)
    assign(Name_data_frame,WQ_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # WTQ_n_save_two <- WTQ_n_save_one
  # WTQ_n_save_one <- WTQ_n
  
  WTQ_n_save_one <- WTQ_n
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  TQ_n <- Create_indiv_n_blanck(BQ_n, indiv_type="TQ", Nb_indiv=NB_TQ, Nb_indiv_per_BQ=NB_TQ_PER_BQ, year)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  TQ_n <- Fill_indiv_n(indiv_n=TQ_n, indiv_type="TQ", D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ=NB_TQ_PER_BQ, NB_BQ, year)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=TQ_n, NB_W_Q_or_TQ=NB_W1_TQ)
  
  WTQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
  TQ_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains part of the IDs of surviving queens
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # We add the unknown mates of (open mated) TQs when UNKNOWN_1b == 'O'
  if(pedigree_Pim == F){
    if(UNKNOWN_1b == '0'){ # then we add dummy open mating sires for open mated queens
      pedigree <- Add_unknown_PGM_mated_to_pedigree(pedigree, indiv_n=WTQ_n, NS_Q_or_TQ_open=NS_TQ, NB_D_FECUNDATING) 
    }
  }
  if(pedigree_Pim){
    pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, NB_D_FECUNDATING)
  }else{
    pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ped_A <- Create_ped_A_Q_n(W_indiv_n=WTQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(2*NB_W1_Q+NB_W1_TQ+NB_D))
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=0, NB_W1_indiv=NB_W1_TQ)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  WTQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WTQ_n, A, ped_A) # Saving with yearly name is done after second overwintering. Could be changed
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==TRUE){
    Q_of_PS_n <- Select_PS_from_WTQ(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_BQ, Criterion=PGM_SELECTION_CRITERION)
    
    if(SAVE_ANNUAL_TABLES){
      Name_data_frame <- paste0("Q_of_PS_", year-2)
      assign(Name_data_frame, Q_of_PS_n)
    }
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==TRUE){
    DPQ_n <- Create_DPQ_n_blanck(W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==TRUE){
    DPQ_n <- Fill_DPQ_n(Q_of_PS_n, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS, DPQ_n, Criterion=DPQ_SELECTION_CRITERION)
    
    if(SAVE_ANNUAL_TABLES){
      Name_data_frame <- paste0("DPQ_", year-2)
      assign(Name_data_frame,DPQ_n)
    }
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==FALSE){
    DPQ_n <- Select_DPQs(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_TQ, NB_DPQ, Criterion=DPQ_SELECTION_CRITERION)
    
    if(SAVE_ANNUAL_TABLES){
      Name_data_frame <- paste0("DPQ_", year-2)
      assign(Name_data_frame, DPQ_n)
    }
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  PS_or_TQ_to_WQ_n_and_WQ_n <- Assign_PGM_or_TQ_to_WQ(DPQ_n, WQ_n)
  
  PS_or_TQ_to_WQ_n <- PS_or_TQ_to_WQ_n_and_WQ_n[[1]]
  WQ_n <- PS_or_TQ_to_WQ_n_and_WQ_n[[2]]
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", year)
    assign(Name_data_frame, WQ_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(pedigree_Pim){
    pedigree <- Add_PGM_mated_to_pedigree_Pim(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE)
  }else{
    pedigree <- Add_PGM_mated_to_pedigree(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(pedigree_Pim == FALSE){
    pedigree <- order_off_give_seq2(pedigree, birth_year_Q_offspring=year, NB_Wg_offspring=NB_W1_Q+NB_W1_TQ, ID_1b=UNKNOWN_1b)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  D_n <- Produce_and_assign_D(Q_to_mate=WQ_n, DPQ_n, NB_D, NB_D_FECUNDATING, VAR_A_meiose, year)
  
  if(MATING_SCALE=="Q_x_PS"){
    WQ_n <- select(WQ_n, -ID_PGM)#We have to withdraw column "ID_PGM" because it is not present in WQ_n when MATING_SCALE==Q_x_DPQ. Instead, there will be a column ID_TQ. This difference then creates a table WQ_n with a column ID_PGM which gets transfered to Output_Q_n in the process of updating Output_Q with function Update_Output_Q_from_WQ_n, which then creates multiple columns of ID_PGM... in table Output_Q. By erasing this ID_PGM column from WQ_n just after it was used by Produce_and_assign_D we keep all further functions working in both MATING_SCALE parameters
  }
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("D_", year)
    assign(Name_data_frame,D_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ped_A <- Create_ped_A_D_n(WTQ_with_DPQs=WTQ_n_save_one, D_n,
                            A,
                            DPQ_index=NB_W1_Q, 
                            NB_D)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  A <- Extend_A_with_D_n(A, ped_A, NB_D)
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  WQ_n <- Produce_Wgroup(W_indiv_n=WQ_n, Ds_n=D_n, indiv_type="Q")
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", year)
    assign(Name_data_frame,WQ_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  TD_n <- Produce_indiv_randomly(Nb_indiv=NB_TD, indiv_type="TD", Q_to_mate=WTQ_n, Var_A=(0.5*VAR_A), year)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("TD_", year)
    assign(Name_data_frame,TD_n)
  }
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  WTQ_n <- Produce_Wgroup(W_indiv_n=WTQ_n, Ds_n=TD_n, indiv_type="TQ")#WTQ_n should be the table of TQ born in current year that survived winter
  
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(pedigree_Pim){
    pedigree <- Add_unknown_PGM_mated_to_pedigree_Pim(pedigree, indiv_n=WTQ_n, ID_1b=UNKNOWN_1b)
  }
  
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Output_Q
# initialization
Output_Q$was_selectedBQ <- FALSE
Output_Q$was_selectedBQ[which(Output_Q$ID %in% unique(Output_Q$ID_BQ))] <- TRUE

## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Output_TQ: only works for Q_x_DPQ mating
# initialization
Output_TQ$was_selectedTQ <- FALSE
Output_TQ$was_selectedTQ[which(Output_TQ$ID %in% unique(Output_Q$ID_TQ))] <- TRUE



## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Output_Q
# initialization
Output_Q$Scheme_scenario[is.na(Output_Q$Scheme_scenario)] <- "BaseScen"

## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Output_TQ: only works for Q_x_DPQ mating
# initialization
Output_TQ$Scheme_scenario[is.na(Output_TQ$Scheme_scenario)] <- "BaseScen"


## ---- eval=T---------------------------------------------------------------------------------------------------------------------------------------------
WD_PATH_OUTPUT <- paste0(WD_PATH_RUN, "/Output") #we create a path to a Output subfolder
if(dir.exists(file.path(WD_PATH_OUTPUT))){#if it exist, nothing is done
}else{
  dir.create(WD_PATH_OUTPUT)#else we create the subfolder
}

date_1 <- format(Sys.time())
date_2 <- gsub(" ", "_", date_1)
date_3 <- gsub(":", "_", date_2)

fwrite(Output_Q, file=paste0(WD_PATH_OUTPUT,"/Output_Q_BaseScen.csv"), sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)

fwrite(Output_TQ, file=paste0(WD_PATH_OUTPUT,"/Output_TQ_BaseScen.csv"), sep = ";", dec = ".", row.names = FALSE, col.names = TRUE)

if(PRODUCE_GRAPHS){
## --------------------------------------------------------------------------------------------------------------------------------------------------------
  library(ggplot2)
  library(ggpubr)
  library(officer)
  library(rvg)
  
  WD_PATH_GRAPHS <- paste0(WD_PATH_RUN, "/Graphs")
  if(!(dir.exists(file.path(WD_PATH_GRAPHS)))){#if it does not exist already, we create an output folder
    dir.create(WD_PATH_GRAPHS)
  }
  
  # Creating a graph file in pptx
  if(!file.exists(file.path(paste0(WD_PATH_GRAPHS, "/graphs.pptx")))){ # if the file does not exist already
    read_pptx() %>% print(paste0(WD_PATH_GRAPHS, "/graphs.pptx")) # we create an empty pptx file
  }
  
  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  PATH_AIREML_SCRIPT_script <- gsub(".*?/", "", PATH_AIREML_SCRIPT)
  
  
  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Phenotypic trend
  
  # BV_d_Q per year
  g <- ggplot(Output_Q, aes(x=Q_birth_year, y=Perf_T1+Perf_T2)) +
    # Adding a boxplot for each year, dodge between was_selectedBQ values
    # coloring boxplot outliers in white
    geom_boxplot(aes(group=interaction(as.character.numeric_version(Output_Q$Q_birth_year), was_selectedBQ),
                     color=was_selectedBQ, fill=was_selectedBQ), alpha=0.5, 
                 position=position_dodge(width=0.9), outlier.alpha = 0) +
    geom_point(aes(color=was_selectedBQ),
               position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.9),
               size=0.5,
               alpha=0.2) + 
    # jitter.width = 0.15, dodge.width = 0.4
    scale_color_manual(values=c("#BBBBBB", "#EE6677"),
                       labels=c("unselected candidate", "selected candidate (BQ)"),
                       name="Was selected as BQ") +
    scale_fill_manual(values=c("#BBBBBB", "#EE6677"),
                      labels=c("unselected candidate", "selected candidate (BQ)"),
                      name="Was selected as BQ") +
    # scale_alpha_manual(labels="Was selected as BQ") +
    
    geom_smooth(data=Output_Q %>% filter(Q_birth_year <= LAST_YEAR_BURN_IN), method=lm, se=FALSE, color="#BBBBBB", size=0.5) + # c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year >= LAST_YEAR_BURN_IN), method=lm, se=FALSE, color="#BBBBBB", size=0.5) + # c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year <= LAST_YEAR_BURN_IN), aes(y=BV_d1_Q+BV_m1_Q), method=lm, se=FALSE, color="#CCBB44", linetype = "dashed", size=0.5) + #c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year >= LAST_YEAR_BURN_IN), aes(y=BV_d1_Q+BV_m1_Q), method=lm, se=FALSE, color="#CCBB44", linetype = "dashed", size=0.5) + #c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year <= LAST_YEAR_BURN_IN), aes(y=BV_d2_Q+BV_m2_Q), method=lm, se=FALSE, color="#66CCEE", linetype = "dashed", size=0.5) +
    geom_smooth(data=Output_Q %>% filter(Q_birth_year >= LAST_YEAR_BURN_IN), aes(y=BV_d2_Q+BV_m2_Q), method=lm, se=FALSE, color="#66CCEE", linetype = "dashed", size=0.5) +
    
    geom_vline(xintercept = LAST_YEAR_BURN_IN, linetype="dotted", color = "#4477AA") + 
    
    labs(x="Q birth year", 
         title=paste0("BaseScen\n", NB_Q, "Q, ", "Within_mat", MAT_LINE_SELECTION, ", Within_pat", PAT_LINE_SELECTION, "\n",
                      "VapiaryT1: ", VAR_BEEKEEPER_x_YEAR_T1, ", VapiaryT2: ", VAR_BEEKEEPER_x_YEAR_T2, ", VresT1: ", VAR_E1, ", VresT2: ", VAR_E2, "\n",
                      PATH_AIREML_SCRIPT_script)) +
    # stat_regline_equation(color="#BBBBBB", size = 5) +
    theme_bw()
  
  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Creating a vectorized form of the graph
  fig_vg <- dml(ggobj = g)
  # Adding the graph on a new slide of "graphs_out.pptx"
  read_pptx(paste0(WD_PATH_GRAPHS, "/graphs.pptx")) %>% 
    add_slide(layout = "Title Only", master = "Office Theme") %>% 
    ph_with(fig_vg, location = ph_location(left = 0, top = 1, height = 4, width = 10)) %>% 
    print(paste0(WD_PATH_GRAPHS, "/graphs.pptx"))
  

  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Genetic trend
  # BV_d_Q per year
  g <- ggplot(Output_Q, aes(x=Q_birth_year, y=BV_d1_Q+BV_m1_Q+BV_d2_Q+BV_m2_Q)) +
    # Adding a boxplot for each year, dodge between was_selectedBQ values
    # coloring boxplot outliers in white
    
    # interaction: one group per year and was_selectedBQ level
    geom_boxplot(aes(group=interaction(as.character.numeric_version(Output_Q$Q_birth_year), was_selectedBQ),
                     color=was_selectedBQ, fill=was_selectedBQ), alpha=0.5,
                 position=position_dodge(width=0.9), outlier.alpha = 0) +
    geom_point(aes(
      color=was_selectedBQ),
      position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.9),
      size=0.5,
      alpha=0.2) + 
    # jitter.width = 0.15, dodge.width = 0.4
    scale_color_manual(values=c("#BBBBBB", "#EE6677"),
                       labels=c("unselected candidate", "selected candidate (BQ)"),
                       name="Was selected as BQ") +
    scale_fill_manual(values=c("#BBBBBB", "#EE6677"),
                      labels=c("unselected candidate", "selected candidate (BQ)"),
                      name="Was selected as BQ") +
    
    geom_smooth(data=Output_Q %>% filter(Q_birth_year <= LAST_YEAR_BURN_IN), method=lm, se=FALSE, color="#BBBBBB", size=0.5) + # c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year >= LAST_YEAR_BURN_IN), method=lm, se=FALSE, color="#BBBBBB", size=0.5) + # c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year <= LAST_YEAR_BURN_IN), aes(y=BV_d1_Q+BV_m1_Q), method=lm, se=FALSE, color="#CCBB44", linetype = "dashed", size=0.5) + #c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year >= LAST_YEAR_BURN_IN), aes(y=BV_d1_Q+BV_m1_Q), method=lm, se=FALSE, color="#CCBB44", linetype = "dashed", size=0.5) + #c("#66CCEE", "#EE6677", "#CCBB44")
    geom_smooth(data=Output_Q %>% filter(Q_birth_year <= LAST_YEAR_BURN_IN), aes(y=BV_d2_Q+BV_m2_Q), method=lm, se=FALSE, color="#66CCEE", linetype = "dashed", size=0.5) +
    geom_smooth(data=Output_Q %>% filter(Q_birth_year >= LAST_YEAR_BURN_IN), aes(y=BV_d2_Q+BV_m2_Q), method=lm, se=FALSE, color="#66CCEE", linetype = "dashed", size=0.5) +
    
    geom_vline(xintercept = LAST_YEAR_BURN_IN, linetype="dotted", color = "#4477AA") + 
    
    labs(x="Q birth year", 
         title=paste0("BaseScen\n", NB_Q, "Q, ", "Within_mat", MAT_LINE_SELECTION, ", Within_pat", PAT_LINE_SELECTION, "\n",
                      "VapiaryT1: ", VAR_BEEKEEPER_x_YEAR_T1, ", VapiaryT2: ", VAR_BEEKEEPER_x_YEAR_T2, ", VresT1: ", VAR_E1, ", VresT2: ", VAR_E2, "\n",
                      PATH_AIREML_SCRIPT_script)) +
    # stat_regline_equation(color="#BBBBBB", size = 5) +
    
    # stat_regline_equation(color="#BBBBBB", size = 5) +
    theme_bw()
  
  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Creating a vectorized form of the graph
  fig_vg <- dml(ggobj = g)
  # Adding the graph on a new slide of "graphs_out.pptx"
  read_pptx(paste0(WD_PATH_GRAPHS, "/graphs.pptx")) %>% 
    add_slide(layout = "Title Only", master = "Office Theme") %>% 
    ph_with(fig_vg, location = ph_location(left = 0, top = 1, height = 4, width = 10)) %>% 
    print(paste0(WD_PATH_GRAPHS, "/graphs.pptx"))
  
  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Genetic trend
  # BV_d_Q per year
  g <- ggplot(Output_Q, aes(x=Q_birth_year, y=F_ID)) +
    # Adding a boxplot for each year, dodge between was_selectedBQ values
    # coloring boxplot outliers in white
    
    # interaction: one group per year and was_selectedBQ level
    geom_boxplot(aes(group=interaction(as.character.numeric_version(Output_Q$Q_birth_year), was_selectedBQ),
                     color=was_selectedBQ, fill=was_selectedBQ), alpha=0.5,
                 position=position_dodge(width=0.9), outlier.alpha = 0) +
    geom_point(aes(
      color=was_selectedBQ),
      position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.9),
      size=0.5,
      alpha=0.2) + 
    # jitter.width = 0.15, dodge.width = 0.4
    scale_color_manual(values=c("#BBBBBB", "#EE6677"),
                       labels=c("unselected candidate", "selected candidate (BQ)"),
                       name="Was selected as BQ") +
    scale_fill_manual(values=c("#BBBBBB", "#EE6677"),
                      labels=c("unselected candidate", "selected candidate (BQ)"),
                      name="Was selected as BQ") +
    
    geom_smooth(data=Output_Q %>% filter(Q_birth_year >= LAST_YEAR_BURN_IN), method=lm, se=FALSE, color="#BBBBBB", size=0.5) + # c("#66CCEE", "#EE6677", "#CCBB44")
    
    geom_vline(xintercept = LAST_YEAR_BURN_IN, linetype="dotted", color = "#4477AA") + 
    
    labs(x="Q birth year", 
         title=paste0("BaseScen\n", NB_Q, "Q, ", "Within_mat", MAT_LINE_SELECTION, ", Within_pat", PAT_LINE_SELECTION, "\n",
                      "VapiaryT1: ", VAR_BEEKEEPER_x_YEAR_T1, ", VapiaryT2: ", VAR_BEEKEEPER_x_YEAR_T2, ", VresT1: ", VAR_E1, ", VresT2: ", VAR_E2, "\n",
                      PATH_AIREML_SCRIPT_script)) +
    # stat_regline_equation(color="#BBBBBB", size = 5) +
    
    stat_regline_equation(color="#BBBBBB", size = 5) +
    theme_bw()
  
  ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Creating a vectorized form of the graph
  fig_vg <- dml(ggobj = g)
  # Adding the graph on a new slide of "graphs_out.pptx"
  read_pptx(paste0(WD_PATH_GRAPHS, "/graphs.pptx")) %>% 
    add_slide(layout = "Title Only", master = "Office Theme") %>% 
    ph_with(fig_vg, location = ph_location(left = 0, top = 1, height = 4, width = 10)) %>% 
    print(paste0(WD_PATH_GRAPHS, "/graphs.pptx"))
}

## --------------------------------------------------------------------------------------------------------------------------------------------------------
# We save the R objects of the burn-in
# save.image(paste0("BaseScen_", year, ".RData"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------