#!/bin/bash
#$ -o /travail/tkistler/output_error_Brascamp/output_AIREML.txt
#$ -e /travail/tkistler/output_error_Brascamp/error_AIREML.txt
#$ -M tristan.kistler@inrae.fr
#$ -m bea
#$ -q workq
#$ -l h_vmem=2G
#############################################################################################
##          Script to run AIREML estimating and not estimating                             ##
##      correlation between direct and maternal effects on simulation sets,                ##
##      with an A-1 matrix obtained with the Brascamp prg suite,                           ##
##      used with different values of ND and NS, defining a Brascamp_group                 ##
#############################################################################################

prog=/CHANGE_TO_AIREMLF90_PATH/BLUPF90_prog/
conv_crit=1e-11
maxrounds=1000

WD_PATH_BLUP=${1}
echo "WD_PATH_BLUP is: ${WD_PATH_BLUP}"
SINGLE_TRAIT=${2} # T: single trait, F: 2 traits
echo "SINGLE_TRAIT is: ${SINGLE_TRAIT}"

VA1_D_ini=${3}
echo "VA1_D_ini is: ${VA1_D_ini}"
VA1_M_ini=${4}
echo "VA1_M_ini is: ${VA1_M_ini}"
VA2_D_ini=${5}
echo "VA2_D_ini is: ${VA2_D_ini}"
VA2_M_ini=${6}
echo "VA2_M_ini is: ${VA2_M_ini}"
VE1_ini=${7}
echo "VE1_ini is: ${VE1_ini}"
VE2_ini=${8}
echo "VE2_ini is: ${VE2_ini}"

MODEL_TYPE=${9}
echo "MODEL_TYPE is: ${MODEL_TYPE}"

inv_A_file_name=AINV_f.giv
PERF_FILE=perf_Bras.txt


#----------------------------------------------------------------------
# Creating dir
cd ${WD_PATH_BLUP}
echo "The root directory is ${WD_PATH_BLUP}"

# inverse of the relationship matrix
echo "The inverse of the relationship matrix file name is ${inv_A_file_name}"

# Info columns of perf : 
# seq_Wgroup seq_Q gen_mean_fix Beekeeper Perf_T1 Perf_T2 ID_Wgroup ID 
# seq_Wgroup: 1 
# seq_Q: 2
# gen_mean_fix: 3
# Beekeeper: 4
# Perf_T1: 5
# Perf_T2: 6
# ID_Wgroup: 7
# ID: 8

#Initialisation of the genetic and phenotypic variances

# We read the number of levels of random effects (in last line of inv_A_file_name)
nb_animal_levels=`awk 'NR==1{max = $1 + 0; next} {if ($1 > max) max = $1;} END {print max}' ${inv_A_file_name}`
echo "there are ${nb_animal_levels} animal levels"

# We read the number of non-genetic effects' levels
fixed_GC_levels=`awk 'NR==1{max = $4 + 0; next} {if ($4 > max) max = $4;} END {print max}' ${PERF_FILE}`
echo "There are ${fixed_GC_levels} fixed_GC_levels fixed effect levels"


if [ ${SINGLE_TRAIT} = FALSE ]; then
# Colony model with 2 traits

# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
# Colony model
# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
if [ ${MODEL_TYPE} = C ]; then
cat >renf90_T1_T2_CM_corr1.par <<EOF
DATAFILE
${PERF_FILE} 
NUMBER_OF_TRAITS
  2       
NUMBER_OF_EFFECTS
  4
OBSERVATION(S)
  5 6
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
  1 1     ${nb_animal_levels} cross # worker group num ID (position in AINV)
  2 2     ${nb_animal_levels} cross # queen num ID (position in AINV) 
  3 3     1 cross
  4 4     ${fixed_GC_levels} cross
RANDOM_RESIDUAL VALUES
  ${VE1_ini} 0
  0 ${VE2_ini}  
RANDOM_GROUP
   1 2  # queen and worker genetic effects
RANDOM_TYPE
 user_file
FILE
 ${inv_A_file_name}                                      
(CO)VARIANCES # if cov value is different than 0.0 it will be estimated
 ${VA1_D_ini} 0.1 0.1 0.1
 0.1 ${VA2_D_ini} 0.1 0.1
 0.1 0.1 ${VA1_M_ini} 0.1
 0.1 0.1 0.1 ${VA2_M_ini}
OPTION alpha_size 50
OPTION missing -999
OPTION conv_crit ${conv_crit}
OPTION maxrounds ${maxrounds}
OPTION sol se

EOF

# _ _ _ _ _ _ _ _ _ _ RUN AIREML _ _ _ _ _ _ _ _ _ _ 
echo renf90_T1_T2_CM_corr1.par | ${prog}/airemlf90 | tee aireml_T1_T2_CM_corr1.lst

cp solutions aireml_T1_T2_CM_corr1.sol
cp airemlf90.log aireml_T1_T2_CM_corr1.log

echo "Estimation finished : solutions are stored in aireml_T1_T2_CM_corr1.sol and variance components in file aireml_T1_T2_CM_corr1.log"

#---------------------------------------------------------------


sleep 3 # we pause the script to let the copies finish, otherwise, there could be mix-ups between the generic files produced by AIREMLF90 (with generic names) relative to the first run or the second one
fi


# Queen model with 2 traits

# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
# Queen model
# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
if [ ${MODEL_TYPE} = Q ]; then
echo "MODEL_TYPE is: ${MODEL_TYPE}"

cat >renf90_T1_T2_QM_corr1.par <<EOF
DATAFILE
${PERF_FILE} 
NUMBER_OF_TRAITS
  2       
NUMBER_OF_EFFECTS
  3
OBSERVATION(S)
  5 6
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
  2 2     ${nb_animal_levels} cross # queen num ID (position in AINV)
  3 3     1 cross
  4 4     ${fixed_GC_levels} cross
RANDOM_RESIDUAL VALUES
  ${VE1_ini}  0.1
  0.1 ${VE2_ini}  
RANDOM_GROUP
   1   # queen genetic effects
RANDOM_TYPE
 user_file
FILE
 ${inv_A_file_name}                                      
(CO)VARIANCES # if cov value is different than 0.0 it will be estimated
 ${VA1_ini} 0.1
 0.1 ${VA2_ini} 
OPTION alpha_size 50
OPTION missing -999
OPTION conv_crit ${conv_crit}
OPTION maxrounds ${maxrounds}
OPTION sol se

EOF

# _ _ _ _ _ _ _ _ _ _ RUN AIREML _ _ _ _ _ _ _ _ _ _ 
echo renf90_T1_T2_QM_corr1.par | ${prog}/airemlf90 | tee aireml_T1_T2_QM_corr1.lst

cp solutions aireml_T1_T2_QM_corr1.sol
cp airemlf90.log aireml_T1_T2_QM_corr1.log

echo "Estimation finished : solutions are stored in aireml_T1_T2_QM_corr1.sol and variance components in file aireml_T1_T2_QM_corr1.log"

#---------------------------------------------------------------


sleep 3 # we pause the script to let the copies finish, otherwise, there could be mix-ups between the generic files produced by AIREMLF90 (with generic names) relative to the first run or the second one
fi


# Worker model with 2 traits

# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
# Worker model
# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
if [ ${MODEL_TYPE} = W ]; then
cat >renf90_T1_T2_WM_corr1.par <<EOF
DATAFILE
${PERF_FILE} 
NUMBER_OF_TRAITS
  2       
NUMBER_OF_EFFECTS
  3
OBSERVATION(S)
  5 6
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
  1 1     ${nb_animal_levels} cross # worker group num ID (position in AINV) 
  3 3     1 cross
  4 4     ${fixed_GC_levels} cross
RANDOM_RESIDUAL VALUES
  ${VE1_ini}  0.1
  0.1 ${VE2_ini}  
RANDOM_GROUP
   1 # worker and queen genetic effects
RANDOM_TYPE
 user_file
FILE
 ${inv_A_file_name}                                      
(CO)VARIANCES # if cov value is different than 0.0 it will be estimated
 ${VA1_ini} 0.1
 0.1 ${VA2_ini} 
OPTION alpha_size 50
OPTION missing -999
OPTION conv_crit ${conv_crit}
OPTION maxrounds ${maxrounds}
OPTION sol se

EOF

# _ _ _ _ _ _ _ _ _ _ RUN AIREML _ _ _ _ _ _ _ _ _ _ 
echo renf90_T1_T2_WM_corr1.par | ${prog}/airemlf90 | tee aireml_T1_T2_WM_corr1.lst

cp solutions aireml_T1_T2_WM_corr1.sol
cp airemlf90.log aireml_T1_T2_WM_corr1.log

echo "Estimation finished : solutions are stored in aireml_T1_T2_WM_corr1.sol and variance components in file aireml_T1_T2_WM_corr1.log"

#---------------------------------------------------------------


sleep 3 # we pause the script to let the copies finish, otherwise, there could be mix-ups between the generic files produced by AIREMLF90 (with generic names) relative to the first run or the second one
fi

fi # end of 2-traits estimation




















if [ ${SINGLE_TRAIT} = TRUE ]; then
# Colony model with 2 traits

# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
# Colony model
# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
if [ ${MODEL_TYPE} = C ]; then
cat >renf90_T1_CM_corr1.par <<EOF
DATAFILE
${PERF_FILE} 
NUMBER_OF_TRAITS
  1       
NUMBER_OF_EFFECTS
  4
OBSERVATION(S)
  5
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
  1     ${nb_animal_levels} cross # worker group num ID (position in AINV)
  2     ${nb_animal_levels} cross # queen num ID (position in AINV)
  3     1 cross
  4     ${fixed_GC_levels} cross
RANDOM_RESIDUAL VALUES
  ${VE1_ini} 
RANDOM_GROUP
   1 2  # queen and worker genetic effects
RANDOM_TYPE
 user_file
FILE
 ${inv_A_file_name}                                      
(CO)VARIANCES # if cov value is different than 0.0 it will be estimated
 ${VA1_D_ini} 0.1
 0.1 ${VA1_M_ini}

OPTION alpha_size 50
OPTION missing -999
OPTION conv_crit ${conv_crit}
OPTION maxrounds ${maxrounds}
OPTION sol se

EOF

# _ _ _ _ _ _ _ _ _ _ RUN AIREML _ _ _ _ _ _ _ _ _ _ 
echo renf90_T1_CM_corr1.par | ${prog}/airemlf90 | tee aireml_T1_CM_corr1.lst

cp solutions aireml_T1_CM_corr1.sol
cp airemlf90.log aireml_T1_CM_corr1.log

echo "Estimation finished : solutions are stored in aireml_T1_CM_corr1.sol and variance components in file aireml_T1_CM_corr1.log"

#---------------------------------------------------------------


sleep 3 # we pause the script to let the copies finish, otherwise, there could be mix-ups between the generic files produced by AIREMLF90 (with generic names) relative to the first run or the second one
fi


# Queen model with 2 traits

# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
# Queen model
# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
if [ ${MODEL_TYPE} = Q ]; then
echo "MODEL_TYPE is: ${MODEL_TYPE}"

cat >renf90_T1_QM_corr1.par <<EOF
DATAFILE
${PERF_FILE} 
NUMBER_OF_TRAITS
  1       
NUMBER_OF_EFFECTS
  3
OBSERVATION(S)
  5
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
  2     ${nb_animal_levels} cross # queen num ID (position in AINV)
  3     1 cross
  4     ${fixed_GC_levels} cross
RANDOM_RESIDUAL VALUES
  ${VE1_ini} 
RANDOM_GROUP
   1   # queen genetic effects
RANDOM_TYPE
 user_file
FILE
 ${inv_A_file_name}                                      
(CO)VARIANCES # if cov value is different than 0.0 it will be estimated
 ${VA1_ini}
OPTION alpha_size 50
OPTION missing -999
OPTION conv_crit ${conv_crit}
OPTION maxrounds ${maxrounds}
OPTION sol se

EOF

# _ _ _ _ _ _ _ _ _ _ RUN AIREML _ _ _ _ _ _ _ _ _ _ 
echo renf90_T1_QM_corr1.par | ${prog}/airemlf90 | tee aireml_T1_QM_corr1.lst

cp solutions aireml_T1_QM_corr1.sol
cp airemlf90.log aireml_T1_QM_corr1.log

echo "Estimation finished : solutions are stored in aireml_T1_QM_corr1.sol and variance components in file aireml_T1_QM_corr1.log"

#---------------------------------------------------------------


sleep 3 # we pause the script to let the copies finish, otherwise, there could be mix-ups between the generic files produced by AIREMLF90 (with generic names) relative to the first run or the second one
fi


# Worker model with 2 traits

# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
# Worker model
# _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _ __ _ _ _ _ _
if [ ${MODEL_TYPE} = W ]; then
cat >renf90_T1_WM_corr1.par <<EOF
DATAFILE
${PERF_FILE} 
NUMBER_OF_TRAITS
  1       
NUMBER_OF_EFFECTS
  3
OBSERVATION(S)
  5
WEIGHT(S)
 
EFFECTS: POSITIONS_IN_DATAFILE NUMBER_OF_LEVELS TYPE_OF_EFFECT[EFFECT NESTED]
  1     ${nb_animal_levels} cross # worker group num ID (position in AINV)  
  3     1 cross
  4     ${fixed_GC_levels} cross
RANDOM_RESIDUAL VALUES
  ${VE1_ini}  
RANDOM_GROUP
   1 # worker and queen genetic effects
RANDOM_TYPE
 user_file
FILE
 ${inv_A_file_name}                                      
(CO)VARIANCES # if cov value is different than 0.0 it will be estimated
 ${VA1_ini}
OPTION alpha_size 50
OPTION missing -999
OPTION conv_crit ${conv_crit}
OPTION maxrounds ${maxrounds}
OPTION sol se

EOF

# _ _ _ _ _ _ _ _ _ _ RUN AIREML _ _ _ _ _ _ _ _ _ _ 
echo renf90_T1_WM_corr1.par | ${prog}/airemlf90 | tee aireml_T1_WM_corr1.lst

cp solutions aireml_T1_WM_corr1.sol
cp airemlf90.log aireml_T1_WM_corr1.log

echo "Estimation finished : solutions are stored in aireml_T1_WM_corr1.sol and variance components in file aireml_T1_WM_corr1.log"

#---------------------------------------------------------------


sleep 3 # we pause the script to let the copies finish, otherwise, there could be mix-ups between the generic files produced by AIREMLF90 (with generic names) relative to the first run or the second one
fi

fi # end of single trait estimation