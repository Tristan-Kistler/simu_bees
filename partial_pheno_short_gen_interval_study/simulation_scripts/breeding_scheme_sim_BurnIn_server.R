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
WD_PATH_Group_id <- arg_shell[3]  # folder for results of repetited runs of a same scenario
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
source(paste0(WD_PATH_SCRIPTS, "/breeding_functions.R"))
source(paste0(WD_PATH_SCRIPTS, "/matrix_functions.R"))
source(paste0(WD_PATH_SCRIPTS, "/variable_functions.R"))
source(paste0(WD_PATH_SCRIPTS, "/output_functions.R"))
# To run Pim's programs, we wills source them.
# To have Pim's programs running in a separate environment, the trick I found was to wrap it inside a function Pim's program.
# When sourcing the file with Pim's program inside a function, the master_script (simulation script) learn the function containing Pim's program.
# This can then be called by the simulation program to run Pim's program.
# Because Pim's program works inside that functions, all variables it creates are private to the function, so they don't pollute the simulation environment.

# if (SEL_BLUP){source(paste0(PATH_PIM_PEDIGREE_SCRIPT), chdir = T, echo = F, local=F)}

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
year <- 1


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rep=1 #this parameter is fixed to 1 in this script where only one iteration of a simulation is performed, but will be set as a variable when repetition are wanted


## --------------------------------------------------------------------------------------------------------------------------------------------------------
BQ_n <- Produce_indiv_randomly(Nb_indiv=NB_BQ, indiv_type="Q", Q_to_mate, Var_A=VAR_A, year=year)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("BQ_", year)
  assign(Name_data_frame, BQ_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=BQ_n, NB_W_Q_or_TQ=NB_BQ)

WBQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
BQ_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains the numerical part of the IDs of surviving queens

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WBQ_", year)
  assign(Name_data_frame, WBQ_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
#Creation of a Identity matrix as initial coancestry table with the base queens
#For this generation, all individuals are considered non inbred and unrelated (base population)
A <- diag(x=1, ncol = NB_BQ, nrow = NB_BQ)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
D_n <- Produce_indiv_randomly(Nb_indiv=NB_BQ*NB_D_FECUNDATING, indiv_type="D", Q_to_mate=WBQ_n, Var_A=VAR_A_HAPLOIDS, year)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("D_", year)
  assign(Name_data_frame,D_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_random_D_n(A, Nb_D=NB_BQ*NB_D_FECUNDATING)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- create_pedigree_Pim_female(WQ_n=WBQ_n, year, NS_Q_or_TQ=NS_Q_open, NB_D_FECUNDATING, Id_1b=UNKNOWN_1b)
}else{
  pedigree <- create_pedigree_female(WQ_n=WBQ_n, year, NS_Q_or_TQ=NS_Q_open, NB_D_FECUNDATING, Id_1b=UNKNOWN_1b)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WBQ_n <- Produce_Wgroup(W_indiv_n=WBQ_n, Ds_n=D_n, indiv_type="Q")

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WBQ_", year)
  assign(Name_data_frame, WBQ_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
year <- 2


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Because the first year we only have BQ, we will ditribute them to only one beekeeper to still use the same function event hough it is thought for distributing Qs, not BQs, which differ in numbers
BQ_n <- Assign_indiv_keepers_max_connexion(Q_or_TQ_n=BQ_n, Nb_keepers=1, Indiv_per_line=1, Q_AND_TQ_SAME_APIARIES, indiv_type="Q", Performance_year=year) # Indiv_per_line is one because here we only have BQs

#We transpose this distribution to actually survived queens:
WBQ_n$Beekeeper_year <- BQ_n$Beekeeper_year[BQ_W1_n]


## --------------------------------------------------------------------------------------------------------------------------------------------------------
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
birth_year_WBQ <- year-1
WBQ_n <- Add_colony_Performance(WQ_or_WTQ_n=WBQ_n, indiv_type="Q", Beekeepers_fixed_effects_n, VAR_E)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WBQ_", birth_year_WBQ)
  assign(Name_data_frame, WBQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Output_Q <- Create_Output_Q(WQ_n=WBQ_n, NB_W1_Q=NB_BQ, Rep)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
#In the first selection, foundress BQ are selected and each initiate a breeding line
birth_year_BQ <- (year-1)
BQ_n <- Select_BQ_foundress(WQ_n=WBQ_n, NB_BQ, Q_SELECTION_CRITERION=Q_SELECTION_CRITERION_INI)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("BQ_", birth_year_BQ) 
  assign(Name_data_frame,BQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_n <- Create_indiv_n_blanck(BQ_n=BQ_n, indiv_type="Q", Nb_indiv=NB_Q, Nb_indiv_per_BQ=NB_Q_PER_BQ, year=year)



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

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", year)
  assign(Name_data_frame,WQ_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim == F){
  if(UNKNOWN_1b == 'OM'){
    WQ_n <- Add_unknown_PGM_mated_to_indiv_n(indiv_n=WQ_n, ID_1b=UNKNOWN_1b) # TO SUPPRESS?
  }else if (UNKNOWN_1b == '0'){ # then we add dummy open mating sires for open mated queens
    pedigree <- Add_unknown_PGM_mated_to_pedigree(pedigree, indiv_n=WQ_n, NS_Q_or_TQ_open=NS_Q_open, NB_D_FECUNDATING) 
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q_open, NB_D_FECUNDATING)
}else{
  pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q_open, birth_year_dams=(year-1), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WQ_n, Dams=WBQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_BQ, index_indiv=(NB_BQ+NB_D_FECUNDATING*NB_BQ))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=0, NB_W1_indiv=NB_W1_Q)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WQ_n, A, ped_A)  

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", year)
  assign(Name_data_frame,WQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
TQ_n <- Create_indiv_n_blanck(BQ_n, indiv_type="TQ", Nb_indiv=NB_TQ, Nb_indiv_per_BQ=NB_TQ_PER_BQ, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
TQ_n <- Fill_indiv_n(indiv_n=TQ_n, indiv_type="TQ", D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ=NB_TQ_PER_BQ, NB_BQ, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=TQ_n, NB_W_Q_or_TQ=NB_W1_TQ)

WTQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
TQ_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains part of the IDs of surviving queens

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WTQ_", year)
  assign(Name_data_frame,WTQ_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# We add the unknown mates of (open mated) TQs when UNKNOWN_1b == 'O'
if(pedigree_Pim == F){
  if(UNKNOWN_1b == '0'){ # then we add dummy open mating sires for open mated queens
    pedigree <- Add_unknown_PGM_mated_to_pedigree(pedigree, indiv_n=WTQ_n, NS_Q_or_TQ_open=NS_TQ, NB_D_FECUNDATING) 
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# We then add the offspring (TQ) queens
if(pedigree_Pim){
  pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, NB_D_FECUNDATING)
}else{
  pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, birth_year_dams=(year-1), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim == FALSE){
  pedigree <- order_off_give_seq2(pedigree, birth_year_Q_offspring=year, NB_Wg_offspring=NB_W1_Q+NB_W1_TQ, ID_1b=UNKNOWN_1b)
}



## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- new_pedigree_PS(pedigree, BQ_n, NB_BQ, NB_W1_Q, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WTQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_BQ, index_indiv=(NB_BQ+NB_BQ*NB_D_FECUNDATING+NB_W1_Q))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=0, NB_W1_indiv=NB_W1_TQ)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WTQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WTQ_n, A, ped_A)  

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WTQ_", year)
  assign(Name_data_frame,WTQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
D_n <- Produce_indiv_randomly(Nb_indiv=NB_D, indiv_type="D", Q_to_mate=WQ_n, Var_A=(0.5*VAR_A), year)
if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("D_", year)
  assign(Name_data_frame,D_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
#adding the fathers (drones) to A
A <- Extend_A_with_random_D_n(A, Nb_D=NB_D)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WQ_n <- Produce_Wgroup(W_indiv_n=WQ_n, Ds_n=D_n, indiv_type="Q")


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- Add_unknown_PGM_mated_to_pedigree_Pim(pedigree, indiv_n=WQ_n, ID_1b=UNKNOWN_1b)
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


## --------------------------------------------------------------------------------------------------------------------------------------------------------
year <- 3


## --------------------------------------------------------------------------------------------------------------------------------------------------------
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
# SHOULD I PUT AN IF FOR THE CASE THAT Qs AND TQs ARE ON SEPARATE APIARIES?
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
Output_TQ <- Create_Output_TQ(WTQ_n, Rep)


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
  assign(Name_data_frame, W2_TQ_n)# The tbl is named by marking the birth year of queens
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
year <- 4


## --------------------------------------------------------------------------------------------------------------------------------------------------------
birth_year_BQ <- (year-1)
BQ_n <- Select_BQ(W2_Q_n, NB_BQ, Q_SELECTION_CRITERION_INI)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("BQ_", birth_year_BQ)
  assign(Name_data_frame,BQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Output_Q <- Add_BQ_sel_dif_to_Output_Q(BQ_n, Output_Q)


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


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- Fill_ped_PS_of_Q(pedigree, New_virgin_Qs=WQ_n)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(NB_Q_TQ_AND_D) )


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=(NB_BQ+NB_BQ*NB_D_FECUNDATING), NB_W1_indiv=NB_W1_Q)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WQ_n, A, ped_A)  

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", year)
  assign(Name_data_frame,WQ_n)
}


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- Fill_ped_F_Q(pedigree, A, ped_A)


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


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# We add the unknown mates of (open mated) TQs when UNKNOWN_1b == 'O'
if(pedigree_Pim == F){
  if(UNKNOWN_1b == '0'){ # then we add dummy open mating sires for open mated queens
    pedigree <- Add_unknown_PGM_mated_to_pedigree(pedigree, indiv_n=WTQ_n, NS_Q_or_TQ_open=NS_TQ, NB_D_FECUNDATING) 
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, NB_D_FECUNDATING)
}else{
  pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
}


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- Fill_ped_PS_of_Q(pedigree, New_virgin_Qs=WTQ_n)


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- new_pedigree_PS(pedigree, BQ_n, NB_BQ, NB_W1_Q, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WTQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(2*NB_W1_Q+NB_W1_TQ+NB_D))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=0, NB_W1_indiv=NB_W1_TQ)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WTQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WTQ_n, A, ped_A)  # Saving with yearly name is done after second overwintering. Could be changed


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==TRUE){
  Q_of_PS_n <- Select_PS_from_WTQ(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_BQ, Criterion=PGM_SELECTION_CRITERION_INI)
  
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
  DPQ_n <- Fill_DPQ_n(Q_of_PS_n, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS, DPQ_n, Criterion=DPQ_SELECTION_CRITERION_INI)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("DPQ_", year-2)
    assign(Name_data_frame,DPQ_n)
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==FALSE){
  DPQ_n <- Select_DPQs(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_TQ, NB_DPQ, Criterion=DPQ_SELECTION_CRITERION_INI)
  
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

# if(MATING_SCALE=="Q_x_PS"){
# WQ_n <- select(WQ_n, -ID_PGM)#We have to withdraw column "ID_PGM" because it is not present in WQ_n when MATING_SCALE==Q_x_DPQ. Instead, there will be a column ID_TQ. This difference then creates a table WQ_n with a column ID_PGM which gets transfered to Output_Q_n in the process of updating Output_Q with function Update_Output_Q_from_WQ_n, which then creates multiple columns of ID_PGM... in table Output_Q. By erasing this ID_PGM column from WQ_n just after it was used by Produce_and_assign_D we keep all further functions working in both MATING_SCALE parameters
# }

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


## --------------------------------------------------------------------------------------------------------------------------------------------------------
year <- 5


## --------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## if(MAT_LINE_SELECTION){
## #ID_Line are factors. For now, Output_Q$ID_Line only contains NA (logical type vector), so we have to transform it to factor column, with adequate factor levels
## Output_Q$ID_Line <- factor(Output_Q$ID_Line, levels = unique(WQ_n$ID_Line))
## }


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
year <- 6


## --------------------------------------------------------------------------------------------------------------------------------------------------------
birth_year_BQ <- (year-1)
BQ_n <- Select_BQ(W2_Q_n, NB_BQ, Q_SELECTION_CRITERION_INI)

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("BQ_", birth_year_BQ)
  assign(Name_data_frame,BQ_n)
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# Output_Q <- Add_BQ_sel_dif_to_Output_Q(BQ_n, Output_Q)


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


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- Fill_ped_PS_of_Q(pedigree, New_virgin_Qs=WQ_n)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(NB_Q_TQ_AND_D) )


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=(NB_Q_TQ_AND_D), NB_W1_indiv=NB_W1_Q)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WQ_n, A, ped_A)  

if(SAVE_ANNUAL_TABLES){
  Name_data_frame <- paste0("WQ_", year)
  assign(Name_data_frame,WQ_n)
}


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- Fill_ped_F_Q(pedigree, A, ped_A)


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


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# We add the unknown mates of (open mated) TQs when UNKNOWN_1b == 'O'
if(pedigree_Pim == F){
  if(UNKNOWN_1b == '0'){ # then we add dummy open mating sires for open mated queens
    pedigree <- Add_unknown_PGM_mated_to_pedigree(pedigree, indiv_n=WTQ_n, NS_Q_or_TQ_open=NS_TQ, NB_D_FECUNDATING) 
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(pedigree_Pim){
  pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, NB_D_FECUNDATING)
}else{
  pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
}


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- Fill_ped_PS_of_Q(pedigree, New_virgin_Qs=WTQ_n)


## ---- eval=F---------------------------------------------------------------------------------------------------------------------------------------------
## pedigree <- new_pedigree_PS(pedigree, BQ_n, NB_BQ, NB_W1_Q, year)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
ped_A <- Create_ped_A_Q_n(W_indiv_n=WTQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(2*NB_W1_Q+NB_W1_TQ+NB_D))


## --------------------------------------------------------------------------------------------------------------------------------------------------------
A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=0, NB_W1_indiv=NB_W1_TQ)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
WTQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WTQ_n, A, ped_A)  # Saving with yearly name is done after second overwintering. Could be changed


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==TRUE){
  Q_of_PS_n <- Select_PS_from_WTQ(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_BQ, Criterion=PGM_SELECTION_CRITERION_INI)
  
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
  DPQ_n <- Fill_DPQ_n(Q_of_PS_n, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS, DPQ_n, Criterion=DPQ_SELECTION_CRITERION_INI)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("DPQ_", year-2)
    assign(Name_data_frame,DPQ_n)
  }
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
if(PAT_LINE_SELECTION==FALSE){
  DPQ_n <- Select_DPQs(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_TQ, NB_DPQ, Criterion=DPQ_SELECTION_CRITERION_INI)
  
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


## --------------------------------------------------------------------------------------------------------------------------------------------------------
#Boucle sur NB_SELECTION_YEARS_AFTER_INITIAL_POP-1 années supplémentaires
for (year in seq(from=year+1, to=(LAST_YEAR_BURN_IN), by=2)){
  tn_start <- proc.time()
  
  ## -------------------------------------------------------------------------------------------------------------------
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
  
  ## -------------------------------------------------------------------------------------------------------------------
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
  
  ## -------------------------------------------------------------------------------------------------------------------
  birth_year_WQ <- year-1
  WQ_n <- Add_colony_Performance(WQ_or_WTQ_n=WQ_n, indiv_type="Q", Beekeepers_fixed_effects_n, VAR_E)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", birth_year_WQ)
    assign(Name_data_frame,WQ_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  # if(MAT_LINE_SELECTION){
  # #ID_Line are factors. For now, Output_Q$ID_Line only contains NA (logical type vector), so we have to transform it to factor column, with adequate factor levels
  # Output_Q$ID_Line <- factor(Output_Q$ID_Line, levels = unique(WQ_n$ID_Line))
  # }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Output_Q <- Update_Output_Q_from_WQ_n(WQ_n, Output_Q, pedigree, year, Rep)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
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
  
  
  # If TQs are on different apiaries than Qs, we use the effect of Testers to fill the Beekeepers_fixed_effects_n table, and draw new effects for the beekeeper_x_year effect
  if(Q_AND_TQ_SAME_APIARIES == FALSE){
    # If Qs and TQs are or different apiaries:
    Beekeepers_fixed_effects_n <- Get_beekeepers_fixed_effects(Nb_keepers=NB_TESTERS, indiv_type="TQ", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES)
  }
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("Beekeepers_fixed_effects_", year)
    assign(Name_data_frame, Beekeepers_fixed_effects_n)
  }
  
  ## -------------------------------------------------------------------------------------------------------------------
  WTQ_n <- Add_colony_Performance(WQ_or_WTQ_n=WTQ_n, indiv_type="TQ", Beekeepers_fixed_effects_n, VAR_E)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WTQ_", birth_year_WQ)
    assign(Name_data_frame,WTQ_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Output_TQ <- Update_Output_TQ(WTQ_n, Output_TQ, year)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W2_n <- Winter_loss_Q(Q_or_TQ_n=WQ_n, NB_W_Q_or_TQ=NB_W2_Q)
  
  W2_Q_n <- Q_or_TQ_and_Q_W2_n[[1]]#WQ stands for Wintered Queens
  Q_W2_n <- Q_or_TQ_and_Q_W2_n[[2]]#contains part of the IDs of surviving queens
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("W2_Q_", (year-1))
    assign(Name_data_frame, W2_Q_n)# The tbl is named by marking the birth year of queens
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W2_n <- Winter_loss_Q(Q_or_TQ_n=WTQ_n, NB_W_Q_or_TQ=NB_W2_TQ)
  
  W2_TQ_n <- Q_or_TQ_and_Q_W2_n[[1]]#WQ stands for Wintered Queens
  TQ_W2_n <- Q_or_TQ_and_Q_W2_n[[2]]#contains part of the IDs of surviving queens
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("W2_TQ_", (year-1))
    assign(Name_data_frame,W2_TQ_n)# The tbl is named by marking the birth year of queens
  }
  
  # We save W2_TQ_n, to not replace it with next gen TQs, useful after burn in in AltScenario for selecting TQs of year n-3
  W2_TQ_n_save_one <- W2_TQ_n
  
  ## -------------------------------------------------------------------------------------------------------------------
  year <- year+1
  
  ## -------------------------------------------------------------------------------------------------------------------
  birth_year_BQ <- (year-1)
  BQ_n <- Select_BQ(W2_Q_n, NB_BQ, Q_SELECTION_CRITERION_INI)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("BQ_", birth_year_BQ)
    assign(Name_data_frame,BQ_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  # Output_Q <- Add_BQ_sel_dif_to_Output_Q(BQ_n, Output_Q)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  WQ_n_save_one <- WQ_n #Q_n_save_one will store queend of the previous year
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Q_n <- Create_indiv_n_blanck(BQ_n, indiv_type="Q", Nb_indiv=NB_Q, Nb_indiv_per_BQ=NB_Q_PER_BQ, year)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Q_n <- Fill_indiv_n(indiv_n=Q_n, indiv_type="Q", D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ=NB_Q_PER_BQ, NB_BQ, year)
  
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("Q_", year)
    assign(Name_data_frame,Q_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=Q_n, NB_W_Q_or_TQ=NB_W1_Q)
  
  WQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
  Q_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains part of the IDs of surviving queens
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  if(pedigree_Pim){
    pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, NB_D_FECUNDATING)
  }else{
    pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  # pedigree <- Fill_ped_PS_of_Q(pedigree, New_virgin_Qs=WQ_n)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  ped_A <- Create_ped_A_Q_n(W_indiv_n=WQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(NB_Q_TQ_AND_D) )
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=(NB_Q_TQ_AND_D), NB_W1_indiv=NB_W1_Q)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  WQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WQ_n, A, ped_A)  
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", year)
    assign(Name_data_frame,WQ_n)
  }
  
  ## -------------------------------------------------------------------------------------------------------------------
  # pedigree <- Fill_ped_F_Q(pedigree, A, ped_A)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  # WTQ_n_save_two <- WTQ_n_save_one
  # WTQ_n_save_one <- WTQ_n
  
  WTQ_n_save_one <- WTQ_n
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  TQ_n <- Create_indiv_n_blanck(BQ_n, indiv_type="TQ", Nb_indiv=NB_TQ, Nb_indiv_per_BQ=NB_TQ_PER_BQ, year)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  TQ_n <- Fill_indiv_n(indiv_n=TQ_n, indiv_type="TQ", D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ=NB_TQ_PER_BQ, NB_BQ, year)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  Q_or_TQ_and_Q_W1_n <- Winter_loss_Q(Q_or_TQ_n=TQ_n, NB_W_Q_or_TQ=NB_W1_TQ)
  
  WTQ_n <- Q_or_TQ_and_Q_W1_n[[1]]#WQ stands for Wintered Queens
  TQ_W1_n <- Q_or_TQ_and_Q_W1_n[[2]]#contains part of the IDs of surviving queens
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  # We add the unknown mates of (open mated) TQs when UNKNOWN_1b == 'O'
  if(pedigree_Pim == F){
    if(UNKNOWN_1b == '0'){ # then we add dummy open mating sires for open mated queens
      pedigree <- Add_unknown_PGM_mated_to_pedigree(pedigree, indiv_n=WTQ_n, NS_Q_or_TQ_open=NS_TQ, NB_D_FECUNDATING) 
    }
  }
  # pedigree <- new_pedigree_PS(pedigree, BQ_n, NB_BQ, NB_W1_Q, year)
  if(pedigree_Pim){
    pedigree <- new_pedigree_Pim_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, NB_D_FECUNDATING)
  }else{
    pedigree <- new_pedigree_female(pedigree, New_virgin_Qs=WTQ_n, NS_Q_or_TQ=NS_TQ, birth_year_dams=(year-2), NB_D_FECUNDATING, ID_1b=UNKNOWN_1b)
  }
  
  ## -------------------------------------------------------------------------------------------------------------------
  # pedigree <- Fill_ped_PS_of_Q(pedigree, New_virgin_Qs=WTQ_n)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  ped_A <- Create_ped_A_Q_n(W_indiv_n=WTQ_n, Dams=BQ_n, Hap_Sires=D_n, index_dams=0, index_D=NB_W1_Q+NB_W1_TQ, index_indiv=(2*NB_W1_Q+NB_W1_TQ+NB_D))
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  A <- Extend_A_with_Q_n(A, ped_A, NB_indiv_to_truncate=0, NB_W1_indiv=NB_W1_TQ)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  WTQ_n <- Fill_WQ_or_WTQ_n_with_F(WQ_or_WTQ_n=WTQ_n, A, ped_A) # Saving with yearly name is done after second overwintering. Could be changed
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==TRUE){
    Q_of_PS_n <- Select_PS_from_WTQ(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_BQ, Criterion=PGM_SELECTION_CRITERION_INI)
    
    if(SAVE_ANNUAL_TABLES){
      Name_data_frame <- paste0("Q_of_PS_", year-2)
      assign(Name_data_frame, Q_of_PS_n)
    }
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==TRUE){
    DPQ_n <- Create_DPQ_n_blanck(W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==TRUE){
    DPQ_n <- Fill_DPQ_n(Q_of_PS_n, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS, DPQ_n, Criterion=DPQ_SELECTION_CRITERION_INI)
    
    if(SAVE_ANNUAL_TABLES){
      Name_data_frame <- paste0("DPQ_", year-2)
      assign(Name_data_frame,DPQ_n)
    }
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  if(PAT_LINE_SELECTION==FALSE){
    DPQ_n <- Select_DPQs(W1_CandidateDPQ=WTQ_n_save_one, W2_CandidateDPQ=W2_TQ_n, NB_TQ, NB_DPQ, Criterion=DPQ_SELECTION_CRITERION_INI)
    
    if(SAVE_ANNUAL_TABLES){
      Name_data_frame <- paste0("DPQ_", year-2)
      assign(Name_data_frame, DPQ_n)
    }
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  PS_or_TQ_to_WQ_n_and_WQ_n <- Assign_PGM_or_TQ_to_WQ(DPQ_n, WQ_n)
  
  PS_or_TQ_to_WQ_n <- PS_or_TQ_to_WQ_n_and_WQ_n[[1]]
  WQ_n <- PS_or_TQ_to_WQ_n_and_WQ_n[[2]]
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", year)
    assign(Name_data_frame, WQ_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  if(pedigree_Pim){
    pedigree <- Add_PGM_mated_to_pedigree_Pim(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE)
  }else{
    pedigree <- Add_PGM_mated_to_pedigree(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE)
  }
  
  
  if(pedigree_Pim == FALSE){
    pedigree <- order_off_give_seq2(pedigree, birth_year_Q_offspring=year, NB_Wg_offspring=NB_W1_Q+NB_W1_TQ, ID_1b=UNKNOWN_1b)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  D_n <- Produce_and_assign_D(Q_to_mate=WQ_n, DPQ_n, NB_D, NB_D_FECUNDATING, VAR_A_meiose, year)
  
  if(MATING_SCALE=="Q_x_PS"){
    WQ_n <- select(WQ_n, -ID_PGM)#We have to withdraw column "ID_PGM" because it is not present in WQ_n when MATING_SCALE==Q_x_DPQ. Instead, there will be a column ID_TQ. This difference then creates a table WQ_n with a column ID_PGM which gets transfered to Output_Q_n in the process of updating Output_Q with function Update_Output_Q_from_WQ_n, which then creates multiple columns of ID_PGM... in table Output_Q. By erasing this ID_PGM column from WQ_n just after it was used by Produce_and_assign_D we keep all further functions working in both MATING_SCALE parameters
  }
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("D_", year)
    assign(Name_data_frame,D_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  ped_A <- Create_ped_A_D_n(WTQ_with_DPQs=WTQ_n_save_one, D_n,
                            A,
                            DPQ_index=NB_W1_Q, 
                            NB_D)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  A <- Extend_A_with_D_n(A, ped_A, NB_D)
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  WQ_n <- Produce_Wgroup(W_indiv_n=WQ_n, Ds_n=D_n, indiv_type="Q")
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("WQ_", year)
    assign(Name_data_frame,WQ_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  TD_n <- Produce_indiv_randomly(Nb_indiv=NB_TD, indiv_type="TD", Q_to_mate=WTQ_n, Var_A=(0.5*VAR_A), year)
  
  if(SAVE_ANNUAL_TABLES){
    Name_data_frame <- paste0("TD_", year)
    assign(Name_data_frame,TD_n)
  }
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  WTQ_n <- Produce_Wgroup(W_indiv_n=WTQ_n, Ds_n=TD_n, indiv_type="TQ")#WTQ_n should be the table of TQ born in current year that survived winter
  
  
  ## -------------------------------------------------------------------------------------------------------------------
  if(pedigree_Pim){
    pedigree <- Add_unknown_PGM_mated_to_pedigree_Pim(pedigree, indiv_n=WTQ_n, ID_1b=UNKNOWN_1b)
  }
  
  ## -------------------------------------------------------------------------------------------------------------------
  print(paste0("Selection year n.", year, " has been completed" ))
  tn_finish <- proc.time()
  running_time_n <- tn_finish-tn_start
  print(paste0("elapsed time in minutes to simulate year n.", year, ": ", running_time_n[3]/60))
  
  Output_Q$Scheme_scenario <- "BurnIn"
  Output_TQ$Scheme_scenario <- "BurnIn"
}


## --------------------------------------------------------------------------------------------------------------------------------------------------------
t3 <- proc.time()
running_time <- t3-t1
paste0("elapsed time in minutes for whole burn-in simulation: ", running_time[3]/60)


## --------------------------------------------------------------------------------------------------------------------------------------------------------
# We save the R objects of the burn-in
save.image(paste0("Burn_in_", year, ".RData"))


## --------------------------------------------------------------------------------------------------------------------------------------------------------