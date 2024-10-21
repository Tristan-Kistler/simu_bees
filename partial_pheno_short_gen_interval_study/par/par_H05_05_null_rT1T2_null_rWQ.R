#############################################################################
#    These are the parameters names and values for the Simu_bees program    #
#############################################################################

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#General parameters

# LAST_YEAR_BURN_IN <- 10 # must be >/= 3. Represents the number of actual selection years performed, not counting the 4 first years simulated in the program used to generate the initial population
# LAST_YEAR <- 31 # must be >/= LAST_YEAR_BURN_IN+5. Represents the number of actual selection years performed, not counting the 4 first years simulated in the program used to generate the initial population
# REP <- 1 # number of repetition of the whole selection process 

LAST_YEAR_BURN_IN <- 10 # must be >/= 3. Represents the number of actual selection years performed, not counting the 4 first years simulated in the program used to generate the initial population
LAST_YEAR <- 32 # must be >/= LAST_YEAR_BURN_IN+3. Represents the number of actual selection years performed, not counting the 4 first years simulated in the program used to generate the initial population
REP <- 1 # used in interactive mode to give a numerical ID to the rep

NB_Q <- 432 # 48 Determines the total number of Queens produced each year, before eventual winter losses
NB_BQ <- 24 # Determines the number of Qs selected as mother queen for generating Qs of the next generation. In case MAT_LINE_SELECTION==TRUE it also determines the number of maternal lines of the breeding program.
NB_TQ <- 432 # Number of test queens. This are special queens in the GPGR selection group, representing daughters of selected BQ fecundated by random drones, used to test performances of selected BQs and to produce drones (for TQs selected as DPQs).
NB_DPQ <- 24 # total number of DPQs used to produce drones. together with NB_DPQ_PER_PS it determines the number of PGMs used for Drone production (NB_PGM = NB_DPQ/NB_DPQ_PER_PS). NB_PGM must be </= to NB_BQ.

# NB_Q <- 96 # 48 Determines the total number of Queens produced each year, before eventual winter losses
# NB_BQ <- 12 # Determines the number of Qs selected as mother queen for generating Qs of the next generation. In case MAT_LINE_SELECTION==TRUE it also determines the number of maternal lines of the breeding program.
# NB_TQ <- 96 # Number of test queens. This are special queens in the GPGR selection group, representing daughters of selected BQ fecundated by random drones, used to test performances of selected BQs and to produce drones (for TQs selected as DPQs).
# NB_DPQ_PER_PS <- 3 # determines how many sister TQs will be used as DPQs for each selected line (for each Pseudo sire). Commonly reffered to as q in the litterature (for example in "Genetic evaluation in the honey bee considering queen and worker effects – A BLUP-Animal Model approach*.", Bienfeld et al., 2006.
# NB_DPQ <- 24 # total number of DPQs used to produce drones. together with NB_DPQ_PER_PS it determines the number of PGMs used for Drone production (NB_PGM = NB_DPQ/NB_DPQ_PER_PS). NB_PGM must be </= to NB_BQ.

Q_WINTER1_LOSS_RATE <- 0.25 #Annual wintering loss rate for first overwintering. See restriction below
Q_WINTER2_LOSS_RATE <- 0.25 #Annual wintering loss rate for second overwintering. Must be between :0</=Q_WINTER1_LOSS_RATE+Q_WINTER2_LOSS_RATE</=1

# Q_WINTER1_LOSS_RATE <- 0 #Annual wintering loss rate for first overwintering. See restriction below
# Q_WINTER2_LOSS_RATE <- 0 #Annual wintering loss rate for second overwintering. Must be between :0</=Q_WINTER1_LOSS_RATE+Q_WINTER2_LOSS_RATE</=1

NB_D_FECUNDATING <- 8 #Commonly reffered to as d in the litterature (for example in "Genetic evaluation in the honey bee considering queen and worker effects – A BLUP-Animal Model approach*.", Bienfeld et al., 2006.

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Simulation Options: breeding parameters
SEL_BLUP <- TRUE #if TRUE: Selection on BQ will be based en performance BLUPs calculated by BLUPF90
#if FALSE, BQs will be directly selected on the performance of Qs the year before BQs are selected
MAT_LINE_SELECTION <- F # FALSE TRUE if TRUE, maternal line selection with one replacement per sib group, if FALSE, mass selection on maternal path
PAT_LINE_SELECTION <- F #if TRUE, paternal line selection will be used to select PS for matings, if FALSE, DPQs are selected among all DPQs without restrictions on their ascent
MATING_SCALE <- "Q_x_DPQ" #if "Q_x_PS", Qs are fecundated by drones produced by randomly selected DPQs constituing one PS (PAT_LINE_SELECTION needs to be TRUE), if Q_x_DPQ, Qs are fecundated by drones produced by a single DPQ
if(PAT_LINE_SELECTION == TRUE){
  NB_DPQ_PER_PS <- 3 # determines how many sister TQs will be used as DPQs for each selected line (for each Pseudo sire). Commonly reffered to as q in the litterature (for example in "Genetic evaluation in the honey bee considering queen and worker effects – A BLUP-Animal Model approach*.", Bienfeld et al., 2006.
}
PAIRING <- "SEMI_RANDOM" #Indicates how DPQs or PSs producing Ds for Qs are distributed to the Qs. If RANDOM is chosen, this distribution is fully random, if SEMI_RANDOM it is random but equilibrating ditribution of DPQs or PSs to Qs. If MAT_LINE_SELECTION and PAT_LINE_SELECTION are set to TRUE, a 3rd PAIRING parameter is possible: AVOID_INLINE, in which the distribution of DPQs or PSs to Qs is semi random and avoids inline (same line maternally and paternally) matings.

# APAIRY GENETIC CONNEXION
Q_AND_TQ_SAME_APIARIES <- TRUE #if TRUE, Qs and TQs are distributed on the same apiaries, if FALSE, Qs go to breeders and TQs to testers, with different apiary effects
if(Q_AND_TQ_SAME_APIARIES == TRUE){
  NB_KEEPERS <- 18 #Number of beekeepers in which Qs and TQs are distributed
  KEEPER_CONNEXION <- "ANY"# if ANY: NB_BQ >= NB_SISTER_GROUPS_PER_KEEPER, & NB_SISTER_GROUPS_PER_KEEPER * NB_KEEPERS/NB_BQ >= 2, & (NB_BQ*NB_Q_PER_BQ/NB_KEEPERS)/NB_SISTER_GROUPS_PER_KEEPER must be integer, & NB_KEEPERS*NB_SISTER_GROUPS_PER_KEEPER/NB_BQ must be integer, & NB_SISTER_GROUPS_PER_KEEPER <= NB_BQ*NB_Q_PER_BQ/NB_KEEPERS
  if(KEEPER_CONNEXION == "ANY"){
    NB_SISTER_GROUPS_PER_KEEPER <- 4 # Number of sister groups per keeper. Must be a divisor of apiary size (apiary size=NB_Q/NB_KEEPERS and NB_TQ/NB_KEEPERS). Allocation to keepers is done separately for TQs and Qs, so apiary size for each must respect the condition
  }
}else{
  BREEDER_CONNEXION <- "ANY" #if MAXIMAL : NB_BREEDERS </= NB_BQ & NB_Q_PER_BQ must be a multiple of NB_BREEDER
  # if MINIMAL : NB_BREEDER == NB_BQ & NB_BQ > 1 & NB_Q_PER_BQ must be a multiple of 2
  # if ANY : NB_Q_PER_BQ must be a multiple of 4 & NB_BREEDER == NB_BQ & NB_BQ > 3. If NB_BQ </= NB_Q_PER_BQ, it equals MAXIMAL connection
  NB_BREEDERS <- 3 # please read BREEDER_CONNEXION restrictions
  if(BREEDER_CONNEXION == "ANY"){
    NB_Q_SISTER_GROUPS_PER_KEEPER <- 6 # Number of sister groups per keeper. Must be a divisor of apiary size (apiary size=NB_Q/NB_KEEPERS)
  }
  TESTER_CONNEXION <- "ANY" #if MAXIMAL : NB_TESTERS </= NB_BQ & NB_TQ_PER_BQ must be a multiple of NB_TESTER
  # if MINIMAL : NB_TESTER == NB_BQ & NB_BQ > 1 & NB_TQ_PER_BQ must be a multiple of 2
  # if ANY : NB_Q_PER_BQ must be a multiple of 4 & NB_BREEDER == NB_BQ & NB_BQ > 3. If NB_BQ </= NB_Q_PER_BQ, it equals MAXIMAL connection
  NB_TESTERS <- 3 # please read Tester_CONNEXION restrictions
  if(TESTER_CONNEXION == "ANY"){
    NB_TQ_SISTER_GROUPS_PER_KEEPER <- 6 # Number of sister groups per keeper. Must be a divisor of apiary size (apiary size=NB_TQ/NB_KEEPERS)
  }
}

# Selection criterion for initialization (burn in, before BLUP estimation)
Q_SELECTION_CRITERION_INI <- "0.5*(Perf_T1/sd(Perf_T1))+0.5*(Perf_T2/sd(Perf_T2))" # 0.5*(Perf_T1/sd(Perf_T1))+0.5*(Perf_T2/sd(Perf_T2)) Perf_T1/sd(Perf_T1) + Perf_T2/sd(Perf_T2) "Perf_T1" "Perf_T2" "Random"  etc. or simply Random
PGM_SELECTION_CRITERION_INI <- "0.5*(Perf_T1/sd(Perf_T1))+0.5*(Perf_T2/sd(Perf_T2))" # "Perf_T1" "Perf_T2" "Random"  etc. or simply Random to select 4as
DPQ_SELECTION_CRITERION_INI <- "0.5*(Perf_T1/sd(Perf_T1))+0.5*(Perf_T2/sd(Perf_T2))" # "Perf_T1" "Perf_T2" "Random"  etc. or simply Random
# Selection criterion for selection (after burn in, with possibly BLUP estimation)
Q_SELECTION_CRITERION <- "0.5*(EBV_d1_W+EBV_m1_W)+0.5*(EBV_d2_W+EBV_m2_W)" #IF SINGLE TRAIT, CANNOT USE EBVs OF T2 0.5*(EBV_d1_W+EBV_m1_W)+0.5*(EBV_d2_W+EBV_m2_W) "EBV_d1_W+EBV_m1_W+EBV_d2_W+EBV_m2_W" "Perf" "Random" calculus using Blup_d_col, Blup_m_col, Perf etc. or simply Random
PGM_SELECTION_CRITERION <- "0.5*(EBV_d1_W+EBV_m1_W)+0.5*(EBV_d2_W+EBV_m2_W)" # Perf_T1 + Perf_T2 calculus using EBV_d1_Q, EBV_m1_Q, EBV_d2_Q,EBV_m2_Q, Perf etc. or simply Random to select 4as
DPQ_SELECTION_CRITERION <- "0.5*(EBV_d1_Q+EBV_m1_Q)+0.5*(EBV_d2_Q+EBV_m2_Q)" #Perf_T2 calculus using EBV_d1_Q, EBV_m1_Q, EBV_d2_Q,EBV_m2_Q, Perf etc. or simply Random

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#BLUP parameters
NS_TQ <- 100
NS_Q_open <- 100 # for input-pedigree.txt (Pim's program) for open-mated Qs (initial Qs)
Assumed_relationship_founders <- 0 # 0 if unrelated, 1 if equilibrium
Assumed_dist_Poisson_EqualContrib_Doff_DPQoff <- 0 # 0 if Poisson, 1 if equal contribution of Ds and DPQs to a Q's offspring (used to calculate p2 when building A)
UNKNOWN_1b <- "0" # "OM" or "0"

# Getting NS for queens
if(MATING_SCALE == "Q_x_DPQ"){
  NS_Q <- 1
}else if(MATING_SCALE == "Q_x_PS"){ # only possible if PAT_LINE_SELECTION == TRUE
  NS_Q <- NB_DPQ_PER_PS
}
# Program paths
# PATH_TO_pedigree_R <- "/abeilles/Tristan/Scripts/EAAP2022/Scripts_Brascamp_vf/For_Simu_bees/pedigree-18_BeeMuSe-adjust_simu.r"
PATH_AIREML_SCRIPT <- "/abeilles/Tristan/Scripts/EAAP2022/Scripts_AIREML/AIREML_Tristan_Simu_2Traits.sh"
# PATH_AIREML_SCRIPT <- "/abeilles/Tristan/Scripts/EAAP2022/Scripts_AIREML/AIREML_Tristan_Simu_2TraitsCov0.sh"
# PATH_AIREML_SCRIPT <- "/abeilles/Tristan/Scripts/EAAP2022/Scripts_AIREML/AIREML_Tristan_Simu_2TraitsCov0Round1.sh"
STOP_IF_CONV_FAILED <- T
pedigree_Pim <- T
if(pedigree_Pim){
  PATH_PIM_PEDIGREE_SCRIPT <- "/MODIFY_for_path_to_first_program_inv_A/simu_function_Pim_pedigreev20.R"
  PATH_PIM_AINV_SCRIPT <- "/MODIFY_for_path_to_second_program_inv_A/simu_function_AMD-AINV-20.R"
}

# Model type: CM, QM or WM
MODEL_TYPE="C" # C Q W

# Perf to erase when producing perf file
## Base scenario
ERASE_PERF_T1_Q_BASE=F
ERASE_PERF_T2_Q_BASE=F
ERASE_PERF_T1_TQ_BASE=F
ERASE_PERF_T2_TQ_BASE=F
YEAR_ERASE_PERF_T1_Q_BASE='year' # only used if ERASE_PERF_T1_Q=T
YEAR_ERASE_PERF_T2_Q_BASE='year'
YEAR_ERASE_PERF_T1_TQ_BASE='year'
YEAR_ERASE_PERF_T2_TQ_BASE='year'
## Alt scenario
ERASE_PERF_T1_Q_ALT=F
ERASE_PERF_T2_Q_ALT=T
ERASE_PERF_T1_TQ_ALT=F
ERASE_PERF_T2_TQ_ALT=F
YEAR_ERASE_PERF_T1_Q_ALT='year' # only used if ERASE_PERF_T1_Q=T
YEAR_ERASE_PERF_T2_Q_ALT='c((LAST_YEAR_BURN_IN+1):year)'
YEAR_ERASE_PERF_T1_TQ_ALT='year'
YEAR_ERASE_PERF_T2_TQ_ALT='year'

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Simulation options: running parameters
PRODUCE_OUTPUT_TQ <- TRUE #Creates a table with all TQs created in the whole simulation
# CALCULATE_COANCESTRY <- TRUE #calculate a coancestry table between all individuals (Qs and Ds) adjusted for haplo-diploidy 
# CALCULATE_FULL_A_MATRIX <- FALSE #if TRUE: calculates the whole coancestry table, containing all individuals simulated during the program run. It can get extremely big and require a lot of process time and RAM in andvanced selection years if the simulated selection population is big. If FALSE, it keeps a reduced coancestry table sufficient to calculate, with the same precision as if TRUE, coancestry coefficients. This coancestry keeps coancestry coeficients between individuals born in n-2 to n, each year.
SAVE_ANNUAL_TABLES <- F #debug mode. It savec annual tables giving names mentioning the birth year of individuals in those tables
WRITE_REP_TABLES <- FALSE #Writes output files for each rep, in case all reps are not finished, and the program aborted, reps simulated can be save this ways.WRITE_RAW_TABLES
# WRITE_RAW_TABLES <- TRUE
FIXED_SEED <- 0 #if non zero, a seed will be set to obtain repetitive runs (debug) if 0, each run will be different
PRODUCE_GRAPHS <- F #if TRUE, produces graphs in a pptx file

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Genetic parameters
#we charge genetic variance parameters, ie a genetic variance-covariance matrix
# Genetic variances for each of both traits
VAR_A1_D <- 20
VAR_A1_M <- 10
VAR_A2_D <- 20 # 0 20 Put to 0 all T2 variances for single trait simulation
VAR_A2_M <- 10 # 0 10

# VAR_A1_D <- 100
# VAR_A1_M <- 30
# VAR_A2_D <- 0 # 0 20 Put to 0 all T2 variances for single trait simulation
# VAR_A2_M <- 0 # 0 10

# # Genetic covariance among each of both traits
# COV_A1_D_M <- -7.07
# COV_A2_D_M <- -7.07 #direct-maternal genetic covariance
# # Genetic covariance between each of both traits
# COV_A1_A2_D <- -10
# COV_A1_A2_M <- -5
# COV_A1_D_A2_M <- 3.535 #direct-maternal genetic covariance between traits
# COV_A1_M_A2_D <- 3.535 #direct-maternal genetic covariance between traits

# Genetic covariance among each of both traits
COV_A1_D_M <- 0
COV_A2_D_M <- 0 #direct-maternal genetic covariance
# Genetic covariance between each of both traits
COV_A1_A2_D <- 0 # 0 -6 -12
COV_A1_A2_M <- 0 # 0 -3 -6
COV_A1_D_A2_M <- 0 # direct-maternal genetic covariance between traits
COV_A1_M_A2_D <- 0 # direct-maternal genetic covariance between traits


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Environmental effects
#we charge environmental variance parameters
VAR_E1 <- 30 #residual effect of the first trait
VAR_E2 <- 30 #residual effect of the second trait
COV_E1_E2 <- 0 #residual covariance between traits

# VAR_E1 <- 10 #residual effect of the first trait
# VAR_E2 <- 0 #residual effect of the second trait
# COV_E1_E2 <- 0 #residual covariance between traits

# Fixed effects
VAR_BEEKEEPER_x_YEAR_T1 <- 16.5 # effect of the year*beekeeper on trait 1, one drawing per simulated year with phenotyping, same effect for all colonies managed by a same beekeeper a same year
VAR_BEEKEEPER_x_YEAR_T2 <- 33 # effect of the year*beekeeper on trait 1, one drawing per simulated year with phenotyping, same effect for all colonies managed by a same beekeeper a same year

# VAR_BEEKEEPER_x_YEAR_T1 <- 5 # effect of the year*beekeeper on trait 1, one drawing per simulated year with phenotyping, same effect for all colonies managed by a same beekeeper a same year
# VAR_BEEKEEPER_x_YEAR_T2 <- 0 # effect of the year*beekeeper on trait 1, one drawing per simulated year with phenotyping, same effect for all colonies managed by a same beekeeper a same year






#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Calculating secondary parameters from the previously user-specified ones °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

##General parameters
NB_Q_PER_BQ <- ceiling(NB_Q/NB_BQ) # we take the ceiling of this value if the user entered a number of Q not divisible by NB_BQ
NB_TQ_PER_BQ <- ceiling(NB_TQ/NB_BQ)
if(PAT_LINE_SELECTION == TRUE){
  NB_PGM <- NB_DPQ/NB_DPQ_PER_PS
}
NB_W1_Q <- floor(NB_Q*(1-Q_WINTER1_LOSS_RATE))#number of queens surviving their first winter
NB_W1_TQ <- floor(NB_TQ*(1-Q_WINTER1_LOSS_RATE))

NB_D <- NB_W1_Q*NB_D_FECUNDATING
NB_TD <- NB_W1_TQ*NB_D_FECUNDATING
NB_Q_AND_D <- NB_W1_Q + NB_D
NB_Q_TQ_AND_D <- NB_Q_AND_D + NB_W1_TQ

NB_W2_Q <- floor(NB_W1_Q*(1-Q_WINTER2_LOSS_RATE))
NB_W2_TQ <- floor(NB_W1_TQ*(1-Q_WINTER2_LOSS_RATE))
NB_W2_BQ <- floor(NB_BQ*(1-Q_WINTER2_LOSS_RATE))

# Genetic conncetion parameters
if(Q_AND_TQ_SAME_APIARIES == TRUE){
  if(KEEPER_CONNEXION == "ANY"){
    NB_KEEPERS_PER_SISTER_GROUP <- NB_SISTER_GROUPS_PER_KEEPER*NB_KEEPERS/NB_BQ
    APIARY_SIZE <- NB_Q/NB_KEEPERS
    NB_SISTERS_PER_GROUP_PER_KEEPER <- APIARY_SIZE/NB_SISTER_GROUPS_PER_KEEPER
  }
}else{ # if Q and TQ are not in the same apiaries
  if(BREEDER_CONNEXION == "ANY"){
    # For Qs
    NB_BREEDERS_PER_SISTER_GROUP <- NB_Q_SISTER_GROUPS_PER_BREEDER*NB_BREEDERS/NB_BQ
    Q_APIARY_SIZE <- NB_Q/NB_BREEDERS
    NB_Q_SISTERS_PER_GROUP_PER_BREEDER <- Q_APIARY_SIZE/NB_Q_SISTER_GROUPS_PER_BREEDER
  }
  if(TESTER_CONNEXION == "ANY"){
    # For TQs
    NB_TESTERS_PER_SISTER_GROUP <- NB_TQ_SISTER_GROUPS_PER_TESTER*NB_TESTERS/NB_BQ
    TQ_APIARY_SIZE <- NB_TQ/NB_TESTERS
    NB_TQ_SISTERS_PER_GROUP_PER_TESTER <- TQ_APIARY_SIZE/NB_TQ_SISTER_GROUPS_PER_TESTER
  }
} 

########
# Checks
########
# NB_Q_PER_BREEDER <- (NB_Q_PER_BQ*NB_BQ)/NB_BREEDERS
# NB_TQ_PER_TESTER <- (NB_TQ_PER_BQ*NB_BQ)/NB_TESTERS

# if(NB_Q_PER_BREEDER != ceiling(NB_Q_PER_BREEDER)){
#   print("The number of breeders does not divide the number of queens, please change one of those parameters")
#   stopifnot(NB_Q_PER_BREEDER == ceiling(NB_Q_PER_BREEDER))
# }

# if(NB_TQ_PER_TESTER != ceiling(NB_TQ_PER_TESTER)){
#   print("The number of Testers does not divide the number of testor queens, please change one of those parameters")
#   stopifnot(NB_TQ_PER_TESTER==ceiling(NB_TQ_PER_TESTER) )
# }

# # Checking for impossible combinations of parameters
# if(BREEDER_CONNEXION == "MAXIMAL"){
#   stopifnot(NB_BREEDERS <= NB_BQ)
#   stopifnot(NB_Q_PER_BQ %% NB_BREEDERS == 0)
# }else if(BREEDER_CONNEXION == "MINIMAL"){
#   stopifnot(NB_BREEDERS == NB_BQ)
#   stopifnot(NB_BQ > 1)
#   stopifnot(NB_Q_PER_BQ %% 2 == 0)
# }

# if ANY: NB_BQ >= NB_SISTER_GROUPS_PER_KEEPER, & NB_SISTER_GROUPS_PER_KEEPER * NB_KEEPERS/NB_BQ >= 2, & (NB_BQ*NB_Q_PER_BQ/NB_KEEPERS)/NB_SISTER_GROUPS_PER_KEEPER must be integer, & NB_KEEPERS*NB_SISTER_GROUPS_PER_KEEPER/NB_BQ must be integer, & NB_SISTER_GROUPS_PER_KEEPER <= NB_BQ*NB_Q_PER_BQ/NB_KEEPERS
if(Q_AND_TQ_SAME_APIARIES == TRUE){
  if(KEEPER_CONNEXION == "ANY"){
    if(NB_BQ < NB_SISTER_GROUPS_PER_KEEPER){
      print(paste0("NB_SISTER_GROUPS_PER_KEEPER: ", NB_SISTER_GROUPS_PER_KEEPER, " should be larger or equal to NB_BQ: ", NB_BQ, ". Change one of those parameters"))
      stop()
    }
    if(NB_SISTER_GROUPS_PER_KEEPER*NB_KEEPERS/NB_BQ < 2){
      print(paste0("NB_SISTER_GROUPS_PER_KEEPER * NB_KEEPERS / NB_BQ: ", NB_SISTER_GROUPS_PER_KEEPER * NB_KEEPERS / NB_BQ, " should be greater or equal to 2. Change one of those parameters"))
      stop()
    }
    if(NB_SISTER_GROUPS_PER_KEEPER > NB_BQ*NB_Q_PER_BQ/NB_KEEPERS){
      print(paste0("NB_SISTER_GROUPS_PER_KEEPER: ", NB_SISTER_GROUPS_PER_KEEPER, " should be smaller or equal to NB_BQ*NB_Q_PER_BQ/NB_KEEPERS: ", NB_BQ*NB_Q_PER_BQ/NB_KEEPERS, ". Change one of those parameters"))
      stop()
    }
    if(floor((NB_BQ*NB_Q_PER_BQ/NB_KEEPERS)/NB_SISTER_GROUPS_PER_KEEPER) != ((NB_BQ*NB_Q_PER_BQ/NB_KEEPERS)/NB_SISTER_GROUPS_PER_KEEPER)){ # if it's not integer
      print(paste0("(NB_BQ*NB_Q_PER_BQ/NB_KEEPERS)/NB_SISTER_GROUPS_PER_KEEPER: ", (NB_BQ*NB_Q_PER_BQ/NB_KEEPERS)/NB_SISTER_GROUPS_PER_KEEPER, " should be an integer. Change one of those parameters"))
      stop()
    }
    if(floor(NB_KEEPERS*NB_SISTER_GROUPS_PER_KEEPER/NB_BQ) != (NB_KEEPERS*NB_SISTER_GROUPS_PER_KEEPER/NB_BQ)){
      print(paste0("NB_KEEPERS*NB_SISTER_GROUPS_PER_KEEPER/NB_BQ: ", NB_KEEPERS*NB_SISTER_GROUPS_PER_KEEPER/NB_BQ, " should be an integer. Change one of those parameters"))
      stop()
    }
  }
}else{ # if Q and TQ are not in the same apiaries
  if(BREEDER_CONNEXION == "ANY"){
    # For Qs
    if(NB_BQ < NB_Q_SISTER_GROUPS_PER_BREEDER){
      print(paste0("NB_Q_SISTER_GROUPS_PER_BREEDER: ", NB_Q_SISTER_GROUPS_PER_BREEDER, " should be larger or equal to NB_BQ: ", NB_BQ, ". Change one of those parameters"))
      stop()
    }
    if(NB_Q_SISTER_GROUPS_PER_BREEDER*NB_BREEDERS/NB_BQ < 2){
      print(paste0("NB_Q_SISTER_GROUPS_PER_BREEDER * NB_BREEDERS / NB_BQ: ", NB_Q_SISTER_GROUPS_PER_BREEDER * NB_BREEDERS / NB_BQ, " should be greater or equal to 2. Change one of those parameters"))
      stop()
    }
    if(NB_Q_SISTER_GROUPS_PER_BREEDER > NB_BQ*NB_Q_PER_BQ/NB_BREEDERS){
      print(paste0("NB_Q_SISTER_GROUPS_PER_BREEDER: ", NB_Q_SISTER_GROUPS_PER_BREEDER, " should be smaller or equal to NB_BQ*NB_Q_PER_BQ/NB_BREEDERS: ", NB_BQ*NB_Q_PER_BQ/NB_BREEDERS, ". Change one of those parameters"))
      stop()
    }
    if(floor((NB_BQ*NB_Q_PER_BQ/NB_BREEDERS)/NB_Q_SISTER_GROUPS_PER_BREEDER) != ((NB_BQ*NB_Q_PER_BQ/NB_BREEDERS)/NB_Q_SISTER_GROUPS_PER_BREEDER)){
      print(paste0("(NB_BQ*NB_Q_PER_BQ/NB_BREEDERS)/NB_Q_SISTER_GROUPS_PER_BREEDER: ", (NB_BQ*NB_Q_PER_BQ/NB_BREEDERS)/NB_Q_SISTER_GROUPS_PER_BREEDER, " should be an integer. Change one of those parameters"))
      stop()
    }
    if(floor(NB_BREEDERS*NB_Q_SISTER_GROUPS_PER_BREEDER/NB_BQ) != (NB_BREEDERS*NB_Q_SISTER_GROUPS_PER_BREEDER/NB_BQ)){
      print(paste0("NB_BREEDERS*NB_Q_SISTER_GROUPS_PER_BREEDER/NB_BQ: ", NB_BREEDERS*NB_Q_SISTER_GROUPS_PER_BREEDER/NB_BQ, " should be an integer. Change one of those parameters"))
      stop()
    }
  }
  if(TESTER_CONNEXION == "ANY"){
    # For TQs
    if(NB_BQ < NB_TQ_SISTER_GROUPS_PER_TESTER){
      print(paste0("NB_TQ_SISTER_GROUPS_PER_TESTER: ", NB_TQ_SISTER_GROUPS_PER_TESTER, " should be larger or equal to NB_BQ: ", NB_BQ, ". Change one of those parameters"))
      stop()
    }
    if(NB_TQ_SISTER_GROUPS_PER_TESTER*NB_TESTERS/NB_BQ < 2){
      print(paste0("NB_TQ_SISTER_GROUPS_PER_TESTER * NB_TESTERS / NB_BQ: ", NB_TQ_SISTER_GROUPS_PER_TESTER * NB_TESTERS / NB_BQ, " should be greater or equal to 2. Change one of those parameters"))
      stop()
    }
    if(NB_TQ_SISTER_GROUPS_PER_TESTER > NB_BQ*NB_Q_PER_BQ/NB_TESTERS){
      print(paste0("NB_TQ_SISTER_GROUPS_PER_TESTER: ", NB_TQ_SISTER_GROUPS_PER_TESTER, " should be smaller or equal to NB_BQ*NB_Q_PER_BQ/NB_TESTERS: ", NB_BQ*NB_Q_PER_BQ/NB_TESTERS, ". Change one of those parameters"))
      stop()
    }
    if(floor((NB_BQ*NB_Q_PER_BQ/NB_TESTERS)/NB_TQ_SISTER_GROUPS_PER_TESTER) != ((NB_BQ*NB_Q_PER_BQ/NB_TESTERS)/NB_TQ_SISTER_GROUPS_PER_TESTER)){
      print(paste0("(NB_BQ*NB_Q_PER_BQ/NB_TESTERS)/NB_TQ_SISTER_GROUPS_PER_TESTER: ", (NB_BQ*NB_Q_PER_BQ/NB_TESTERS)/NB_TQ_SISTER_GROUPS_PER_TESTER, " should be an integer. Change one of those parameters"))
      stop()
    }
    if(floor(NB_TESTERS*NB_TQ_SISTER_GROUPS_PER_TESTER/NB_BQ) != (NB_TESTERS*NB_TQ_SISTER_GROUPS_PER_TESTER/NB_BQ)){
      print(paste0("NB_TESTERS*NB_TQ_SISTER_GROUPS_PER_TESTER/NB_BQ: ", NB_TESTERS*NB_TQ_SISTER_GROUPS_PER_TESTER/NB_BQ, " should be an integer. Change one of those parameters"))
      stop()
    }
  }
} 

#Genetic parameters
#we then create and fill the covariance matrix of the genetic effects 
VAR_A_up <- matrix(0, 4, 4)
VAR_A_up[1,1] <- VAR_A1_M
VAR_A_up[1,2] <- COV_A1_D_M 
VAR_A_up[1,3] <- COV_A1_A2_M
VAR_A_up[1,4] <- COV_A1_M_A2_D
VAR_A_up[2,2] <- VAR_A1_D
VAR_A_up[2,3] <- COV_A1_D_A2_M
VAR_A_up[2,4] <- COV_A1_A2_D
VAR_A_up[3,3] <- VAR_A2_M
VAR_A_up[3,4] <- COV_A2_D_M
VAR_A_up[4,4] <- VAR_A2_D

VAR_A_up_t <- t(VAR_A_up)
diag(VAR_A_up_t) <- 0
VAR_A <- VAR_A_up + VAR_A_up_t

VAR_A_meiose <- 0.25*VAR_A
VAR_A_HAPLOIDS <- 0.5*VAR_A

#Environmental parameters
VAR_E <- matrix(0, 2, 2)
VAR_E[1,1] <- VAR_E1
VAR_E[1,2] <- COV_E1_E2 
VAR_E[2,1] <- COV_E1_E2
VAR_E[2,2] <- VAR_E2


SIGMA_E1 <- sqrt(VAR_E1)
SIGMA_E2 <- sqrt(VAR_E2)

SIGMA_BEEKEEPER_x_YEAR_T1 <- sqrt(VAR_BEEKEEPER_x_YEAR_T1)
SIGMA_BEEKEEPER_x_YEAR_T2 <- sqrt(VAR_BEEKEEPER_x_YEAR_T2)

## AIREMLF90 input parameters
# Starting values for ReML parameter estimation
VA1_D_ini=0.9*VAR_A1_D
VA1_M_ini=0.9*VAR_A1_M
VA2_D_ini=0.9*VAR_A2_D
VA2_M_ini=0.9*VAR_A2_M
VE1_ini=1.1*VAR_E1
VE2_ini=1.1*VAR_E2

# VA1_D_ini=VAR_A1_D
# VA1_M_ini=VAR_A1_M
# VA2_D_ini=VAR_A2_D
# VA2_M_ini=VAR_A2_M
# VE1_ini=VAR_E1
# VE2_ini=VAR_E2

# If all genetic and non genetic variances of trait 2 are 0, we are in a single trait simulation
if( (VAR_A2_D==0)&(VAR_A2_M==0)&(COV_A2_D_M==0)&(COV_A1_A2_D==0)&(COV_A1_A2_M==0)&(COV_A1_D_A2_M==0)&(COV_A1_M_A2_D==0)&(VAR_E2==0)&(COV_E1_E2==0)&(VAR_BEEKEEPER_x_YEAR_T2==0) ){
  SINGLE_TRAIT <- TRUE
}else{
  SINGLE_TRAIT=FALSE
}