#######################################################################################
#    These are the functions which are needed for the calculation                     #
# of the inbreeding coefficients of queens via a haplo-diploid relationship matrix A  #
#                       of the simu_bees program                                      #
#                                                                                     #
#######################################################################################

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Create_ped_A_Q_n                                       #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# This function initiates ped_A for the first Q or TQ descendants of BQs

#New ultimate flexible version
Create_ped_A_Q_n <- function(W_indiv_n, Dams=WQ_n_save_one, Hap_Sires=D_n,
                             index_dams=NB_Q_TQ_AND_D, 
                             index_D=NB_Q_TQ_AND_D + NB_W1_Q + NB_W1_TQ,
                             index_indiv=(NB_Q_TQ_AND_D*2)){
  #index_dams: numerical ID in A from which dams numbers will start (first dam is at index_dam +1)
  #index_indiv: numerical ID in A from which new entries (queens) numbers will start (first queen is at index_indiv +1)
  #index_D: numerical ID in A from which drones (sires of new entries) numbers will start (first drone is at index_D +1)
  
  NB_W1_indiv <- dim(W_indiv_n)[1]
  #=W_indiv_n should contain queens of the current year and Hap_Sires drones that fecundated the BQ, so drones from current year -1. Dams should contain the Q that became BQ and mothered the Qs of W_indiv_n
  #The cohortes in A always appear in this order : 
  #first year : Qs then Ds, second year : Qs, TQs, Ds and TDs. Following years repeat this order.
  
  ped_A <- data.frame(matrix(ncol = 3, nrow = NB_W1_indiv))
  colnames(ped_A) <- c("ID", "ID_BQ", "ID_D")
  
  ped_A$ID <- seq(from=1+index_indiv, to=NB_W1_indiv+index_indiv, by=1) #giving Ids to the newly born queens
  
  for (i in 1:NB_W1_indiv){
    #filling the table giving the ID numbers of the mothers:
    ped_A$ID_BQ[i] <-  which(Dams$ID %in% W_indiv_n$ID_BQ[i]) + index_dams #numerical ID of the mother of a queen. Dams should contain queens of the current year-1.
    
    #filling the table giving the ID numbers of the fathers:
    ped_A$ID_D[i] <- which(Hap_Sires$ID %in% W_indiv_n$ID_D[i]) + index_D
  }
  return(ped_A)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Extend_A_with_Q_n                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Constant-size A matrix ###############################################################
Extend_A_with_Q_n <- function(A, ped_A, NB_indiv_to_truncate=NB_Q_TQ_AND_D, NB_W1_indiv){
  #first year : Qs then Ds, second year : Qs, TQs, Ds and TDs. Following years repeat this order.
  
  # ATTENTION
  # NB_W1_indiv is also the length of ped_A: could be infered instead of entered as argument
  
  #We keep only individuals of the 2 last years in this reduced A matrix
  #Once 3 generations are represented in A, older individuals are withdrawn to keep A at a constant size
  # with only individuals relevant for inbreeding coefficient computation of new individuals
  if (NB_indiv_to_truncate > 0){
    A <- A[-c(1:NB_indiv_to_truncate), -c(1:NB_indiv_to_truncate)]
  }
  
  index_indiv = dim(A)[1]
  # index_indiv=: numerical ID in A from which new entries (queens) numbers will start (first queen is at index_indiv +1)
  # : size of the matrix after truncation
  
  #Adding new columns to A
  A <- cbind(A, matrix(0L, ncol = NB_W1_indiv, nrow = index_indiv))
  #Adding new rows to A
  A <- rbind(A, matrix(0L, ncol = (index_indiv + NB_W1_indiv), nrow = NB_W1_indiv))
  
  #Filling A with the 2*coancestry coefficients
  #We will need again the constants Q_index to navigate between the position of the parents in A and in the yearly produced (of constant size) ped_A.
  js <- c(1:index_indiv)#for columns of the matrix containing before-born individuals
  New_Q <- index_indiv+1
  #we proceed as indicated in Fernando and Grossman, 1989: A[i,j]=0.5A[mother_index,j] + A[father_index,j]
  
  for(i in (New_Q):(index_indiv + NB_W1_indiv)){
    k <- i-index_indiv
    A[i, i] <- 1 + A[ped_A$ID_BQ[k], ped_A$ID_D[k] ]#We then add A[mother_index, father_index] to A[i,i]
    
    A[i, js] <- 0.5*A[ped_A$ID_BQ[k], js] + A[ped_A$ID_D[k], js]
  }
  
  A[(js), (New_Q):(index_indiv + NB_W1_indiv)] <-t(A[(New_Q):(index_indiv + NB_W1_indiv), js])#transposing the submatrix and filling the symmetrical submatrix
  
  #we will now calculate coancestry between newly produced queens(triangular matrix)
  for(i in (index_indiv+2):(index_indiv + NB_W1_indiv)){#for the newly added rows +1
    i_1 <- i-1
    k <- i-index_indiv
    
    A[i, (New_Q):(i_1)] <- 0.5*A[ped_A$ID_BQ[k], (New_Q):(i_1)] + A[ped_A$ID_D[k], (New_Q):(i_1)] #we sum the term relative to the mother and to the father as indexed in ped_A
    
    A[(New_Q):(i_1), i] <- t( A[i, (New_Q):(i_1)] )#we fill the upper-right triangle matrix by symmetry
    
  } # end for i
  return(A)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Create_ped_A_D_n                                       #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Create_ped_A_D_n <- function(WTQ_with_DPQs=WTQ_n_save_two, D_n,
                             A,
                             DPQ_index=NB_W1_Q, 
                             NB_D){
  #WTQ_with_DPQs should contain TQs of the current year minus 2 and D_n drones of the current year
  D_index <- nrow(A)
  
  ped_A <- data.frame(matrix(ncol = 2, nrow = NB_D))
  colnames(ped_A) <- c("ID", "ID_TQ")
  
  ped_A$ID <- seq(from=1+D_index, to=NB_D+D_index, by=1)#giving Ids to the newly born drones
  
  # The cohortes in A always appear in this order : first year : Qs then Ds, second year : Qs, TQs. Following years repeat this order.
  
  for (i in unique(D_n$ID_TQ) ){#for each ID_TQ in D_n
    which_D <- which(D_n$ID_TQ == i)#line of D_n where we will write numerical IDs of the DPQ that produced them
    TQ_num_ID <- which(WTQ_with_DPQs$ID %in% D_n$ID_TQ[which_D]) # yearly TQ num ID
    ped_A$ID_TQ[which_D] <- TQ_num_ID + DPQ_index # TQ num ID in A: Adding DPQ index
  }#end of for loop
  
  return(ped_A)
}#end of function

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Extend_A_with_random_D_n                               #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Extend_A_with_random_D_n <- function(A, Nb_D){
  
  index_indiv = dim(A)[1]
  #adding the columns and lines to A for random drones
  A <- cbind(A, matrix(0L, ncol = Nb_D, nrow = index_indiv))
  
  A <- rbind(A, matrix(0L, ncol = index_indiv + Nb_D, nrow = Nb_D))
  
  #filling A with 0.5 in diagonal for males
  for(i in (index_indiv + 1):(index_indiv + Nb_D)){#for all drones of this year
    A[i,i] <- 0.5
  }
  return(A)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Extend_A_with_D_n                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Extend_A_with_D_n <- function(A, ped_A, NB_D){
  
  index_indiv = dim(A)[1]
  # index_indiv=: numerical ID in A from which new entries (drones) numbers will start (first D is at index_indiv +1)
  # : size of the matrix
  
  #Adding the new columns to A
  A <- cbind(A, matrix(0L, ncol = NB_D, nrow = index_indiv))
  #Adding the new rows to A
  A <- rbind(A, matrix(0L, ncol = (index_indiv + NB_D) , nrow = NB_D) )
  
  #Filling A with the 2*coancestry coefficients
  #We will need again the constants index_indiv to nagigate between the position of the parents in A and in the yearly produces (of constant size) ped_A.
  js <- c(1:index_indiv)#for columns of the matrix containing before-born individuals
  New_D <- index_indiv+1
  
  #we proceed as indicated in Fernando and Grossman, 1989:
  #We first add 0.5 for A[i,i]
  for(i in (New_D):(index_indiv + NB_D)){
    A[i, i] <- 0.5
    
    #A[i,j]=0.5A[mother_index,j]
    k <- i-index_indiv
    A[i, js] <- 0.5*A[ ped_A$ID_TQ[k] , js] #we sum the term relative to the mother and to the father as indexed in ped_A
  }
  A[(js), (New_D):(index_indiv + NB_D)] <-t(A[(New_D):(index_indiv + NB_D), js])#transposing the submatrix and filling the symmetrical submatrix
  
  
  #we will now calculate coancestry between newly produced drones(triangular matrix)
  for(i in (index_indiv+2):(index_indiv + NB_D)){#for the newly added rows +1
    i_1 <- i-1
    k <- i-index_indiv
    A[i, (New_D):(i_1)] <- 0.5*A[ped_A$ID_TQ[k], (New_D):(i_1)]#we sum the term relative to the mother and to the father as indexed in ped_A
    
    A[(New_D):(i_1), i] <- t( A[i, (New_D):(i_1)] )#we fill the upper-right triangle matrix by symmetry
  } # end for i
  
  return(A)
}#end of function

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Fill_WQ_or_WTQ_n_with_F                                   #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Fill_WQ_or_WTQ_n_with_F <- function(WQ_or_WTQ_n, A, ped_A) {
  
  WQ_or_WTQ_n$F_ID <- diag(A[c(ped_A$ID), c(ped_A$ID)]) - 1
  return (WQ_or_WTQ_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                                         pedigree
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                           create_pedigree_Pim_female                                #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# This function is for use of pedigree-X.R of Pim, that uses as input input-pedigree.txt, created with this function below.
create_pedigree_Pim_female <- function(WQ_n=WQ_n, year, NS_Q_or_TQ, NB_D_FECUNDATING, Id_1b=UNKNOWN_1b){
  
  NB_W1_Q <- length(WQ_n$ID)
  
  pedigree <- data.frame(matrix(ncol = 18, nrow = NB_W1_Q))
  colnames(pedigree) <- c("Performance_year", "ID_Breeder", "ID_test_location",
                          "Q_birth_year", "ID", "BQ_birth_year", "ID_BQ", "PGM_mated_birth_year", "ID_PGM_mated",
                          "C10", "C11", "C12", "C13", "NS", "ND", "Q1b_birth_year", "ID_1b", "extra_digit")
  
  pedigree <- pedigree %>% mutate(
    Performance_year = year+1,
    ID_Breeder = "B1",
    ID_test_location = paste0(ID_Breeder, "_", year+1),
    Q_birth_year = year,
    ID = WQ_n$ID,
    BQ_birth_year = 0,
    ID_BQ = "0", #0 codes for missing value in BLUPF90, which will use this table
    PGM_mated_birth_year = 0,
    ID_PGM_mated = "0",
    C10 = "0",
    C11 = "0",
    C12 = "0",
    C13 = "0",
    NS = NS_Q_or_TQ,
    ND = NB_D_FECUNDATING,
    Q1b_birth_year = 0,
    ID_1b = Id_1b,
    extra_digit = "0"
  )
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                           Create_pedigree_female                                #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# This function is for creating directly the files otherwise produced by pedigree-X.R of Pim, using input-pedigree.txt
create_pedigree_female <- function(WQ_n=WQ_n, year, NS_Q_or_TQ, NB_D_FECUNDATING, Id_1b=UNKNOWN_1b){
  
  NB_W1_Q <- length(WQ_n$ID)
  
  # If the unknown sires are entered as '0', we will have a dummy sire for each mating
  if(Id_1b == "0"){
    index_0sires <- NB_W1_Q
  }else if(Id_1b == "OM"){
    index_0sires <- 1
  }
  
  pedigree <- data.frame(matrix(ncol = 22, nrow = NB_W1_Q))
  colnames(pedigree) <- c("ID_seq1", "Sex", "Mate",
                          "Mating_type", "yob_mod", "seq_dam", "seq_sire", 
                          "c8", "c9", "c10", "Test_station", "c12", "c13", "c14", "c15",
                          "seq_1", "seq_2", "yob", "NS", "ND", "ID", "is_SS")
  
  pedigree_sire <- pedigree[1:index_0sires, ]
  pedigree_sire <- pedigree_sire %>% mutate(
    ID_seq1 = seq(1:index_0sires), # we start with unknown sires
    Sex = 2,
    Mate = 0,
    Mating_type = 0, # don't know what it's supposed to code for
    # if Id_1b == "0", yob_mod = 1, else yob_mod = 0
    yob_mod = ifelse(Id_1b == "0", 1, 0),
    seq_dam = 0,
    seq_sire = 0,
    
    c8 = 0,
    c9 = 0,
    c10 = 0,
    Test_station = 0,
    c12 = 0,
    c13 = 0,
    c14 = 0,
    c15 = 0,
    seq_1 = seq(1:index_0sires),
    seq_2 = seq(1:index_0sires),
    yob = ifelse(Id_1b == "0", 1, 0),
    NS = NS_Q_or_TQ,
    ND = NB_D_FECUNDATING
  )
  
  # Id_1b is not in pedigree, so we cannot make a case_when mutation base on its value, so we add the following after:
  if(Id_1b == "0"){
    pedigree_sire$ID <- gsub("Q", "OpSire_Q", WQ_n$ID)
  }else{
    pedigree_sire$ID <- Id_1b
  }
  # We add column is_SS at the end to have it as last column, although its content is independent of 1d_1b
  pedigree_sire$is_SS <- NA
  
  
  # For Qs
  pedigree <- pedigree %>% mutate(
    ID_seq1 = seq(index_0sires+1, index_0sires+NB_W1_Q),
    Sex = 1,
    # Mate = 1
    Mating_type = 0, # don't know what it's supposed to code for
    yob_mod = 1,
    seq_dam = 0,
    seq_sire = 0,
    
    c8 = 0,
    c9 = 0,
    c10 = 0,
    Test_station = 0,
    c12 = 0,
    c13 = 0,
    c14 = 0,
    c15 = 0,
    seq_1 = seq(index_0sires+1, index_0sires+NB_W1_Q),
    seq_2 = seq(index_0sires+1, index_0sires+NB_W1_Q),
    yob = 1,
    NS = NS_Q_or_TQ,
    ND = NB_D_FECUNDATING,
    ID = WQ_n$ID,
    is_SS = NA,
  )
  
  if(Id_1b == "0"){
    pedigree$Mate <- seq(1:index_0sires) # contrary to Pim, we put the correct mate seq here so that new_pedigree female can use the Mate col content to obtain the Q offspring correct sire (Pim has 0 here)
  }else{ # if open mate is 'OM'
    pedigree$Mate <- 1
  }
  
  # For Wgroup of Qs
  pedigree_Wgroup <- pedigree %>% mutate(
    ID_seq1 = seq(index_0sires+NB_W1_Q+1, (index_0sires+NB_W1_Q+NB_W1_Q)),
    Sex = 3,
    Mate = 0,
    yob_mod = yob_mod + 1,
    seq_dam = pedigree$seq_2[1:NB_W1_Q], # beware that pedigree is only the dam pedigree, so with NB_W1_Q lines
    seq_sire = seq(1:index_0sires), # sires are in first position
    
    seq_1 = seq(index_0sires+NB_W1_Q+1, (index_0sires+NB_W1_Q+NB_W1_Q)), # sires, Qs, and then start the Wgroups
    seq_2 = seq((index_0sires+NB_W1_Q+1), (index_0sires+NB_W1_Q+NB_W1_Q)),
    yob = 1,
    NS = 0,
    ND = 0,
    ID = gsub("Q", "WgQ", pedigree$ID),
    is_SS = NA,
    ID_BQ = pedigree$ID
  )
  
  # row_binding the two dataframes
  pedigree <- bind_rows(pedigree_sire, pedigree, pedigree_Wgroup)
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              new_pedigree_Pim_female                                    #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
new_pedigree_Pim_female <- function(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, NB_D_FECUNDATING){
  
  Q_index <- length(pedigree$ID) #index of the last queen in the pedigree before addition of the new ones
  NB_W1_Q <- length(New_virgin_Qs$ID) #number of queens born this year
  
  WQ_n_ped <- New_virgin_Qs %>% #we create a table with the same infos as New_virgin_Qs, but with the same number of rows as pedigree
    select(ID, ID_BQ) %>% 
    mutate(
      Performance_year =  as.numeric(gsub(".*_", "", New_virgin_Qs$ID))+1,
      ID_Breeder = "B1",
      ID_test_location = paste0(ID_Breeder, "_", as.numeric(gsub(".*_", "", New_virgin_Qs$ID))+1),
      Q_birth_year = as.numeric(gsub(".*_", "", New_virgin_Qs$ID)),
      # ID = "0",
      BQ_birth_year = as.numeric(gsub(".*_", "", New_virgin_Qs$ID_BQ)), # We get BQs' birth year from their name
      # ID_BQ = "0",
      # PGM_mated_birth_year = 0,
      # ID_PGM_mated = "0",
      C10 = "0",
      C11 = "0",
      C12 = "0",
      C13 = "0",
      NS = NS_Q_or_TQ,
      ND = NB_D_FECUNDATING,
      # Q1b_birth_year = 0,
      # ID_1b = "0",
      extra_digit = "0"
    )
  
  pedigree <- pedigree %>% #we add to the current pedigree table infos gathered from New_virgin_Qs
    bind_rows(WQ_n_ped)  #F_Q will be filled with the inbreeding coefficient calculated after calculation in A
  
  # for (i in (Q_index+1):(Q_index+NB_W1_Q)){#for every queen born this year
  #   pedigree$ID_PGM[i] <- pedigree %>% filter(pedigree$ID %in% pedigree$ID_BQ[i]) %>% pull(ID_PGM_mated)#ID_PGM contains the ID_PGM_mated corresponding to the BQ. It is , genetically, the mother of the pseudo sire of a newly born queen
  # }
  
  return(pedigree)
}
# pedigree$Q1b_birth_year[newly_mated_Qs] <- gsub(".*_", "", PS_or_TQ_to_WQ_n$ID_TQ)


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              new_pedigree_female                                    #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
new_pedigree_female <- function(pedigree, New_virgin_Qs=WQ_n, NS_Q_or_TQ=NS_Q, birth_year_dams=(year-2),
                                NB_D_FECUNDATING, ID_1b=UNKNOWN_1b){
  
  Q_index <- length(pedigree$ID) #index of the last queen in the pedigree before addition of the new ones
  NB_W1_Q <- length(New_virgin_Qs$ID) #number of queens born this year
  
  New_virgin_Qs_ped <- New_virgin_Qs %>% #we create a table with the same infos as New_virgin_Qs, but with the same number of rows as pedigree
    select(ID, ID_BQ) %>% 
    mutate(
      ID_seq1 = seq((Q_index+1), (Q_index+NB_W1_Q)),
      Sex = 1,
      # Mate = 1, # default value: mate is OM
      Mating_type = 0, # don't know what it's supposed to code for
      yob_mod = as.numeric(gsub(".*_", "", New_virgin_Qs$ID)),
      # seq_dam = 0,
      # seq_sire = 0, #defaults value
      
      c8 = 0,
      c9 = 0,
      c10 = 0,
      Test_station = 0,
      c12 = 0,
      c13 = 0,
      c14 = 0,
      c15 = 0,
      seq_1 = seq((Q_index+1), (Q_index+NB_W1_Q)),
      seq_2 = seq((Q_index+1), (Q_index+NB_W1_Q)),
      yob = as.numeric(gsub(".*_", "", New_virgin_Qs$ID)),
      NS = NS_Q_or_TQ,
      ND = NB_D_FECUNDATING,
      is_SS = NA,
    )
  # Adding the new Qs' mate
  if(ID_1b != '0'){
    New_virgin_Qs_ped$Mate <- 1 # default value: mate is OM. Will be updated if Qs were controlled mated
  }else if(ID_1b == '0'){
    ##########################################################################################################################################################
    # PROBLEM: the seq_2 at this point are not unique (Wg have temporary ones), seq_2 of 0sires is not the definitive one
    # 0sires (open mating different from OM) and controlled mates will be added Add_PGM* functions and order_off_give_seq2
    # For now, they are NA, needed for proper updating later
  }
  ##########################################################################################################################################################
  
  # Adding the new Qs' seq_dam and seq_sire (which has already been ordered and given a definitive seq_2 by order function)
  # seq_dam
  pedigree_dams <- pedigree %>% filter(yob_mod==birth_year_dams) %>% select(ID, seq_2, Mate) %>% rename(seq_dam=seq_2, seq_sire=Mate)
  New_virgin_Qs_ped <- New_virgin_Qs_ped %>% left_join(pedigree_dams, by=c("ID_BQ"="ID"))
  
  # Binding the Q_offspring set to the pedigree
  pedigree <- pedigree %>% #we add to the current pedigree table infos gathered from New_virgin_Qs
    bind_rows(New_virgin_Qs_ped) # %>%
  
  
  # Adding the pedigree of new queens' worker groups (their offspring)
  pedigree <- pedigree %>% bind_rows(
    New_virgin_Qs_ped %>%
      mutate(
        ID_seq1 = seq((Q_index+NB_W1_Q+1), (Q_index+NB_W1_Q+NB_W1_Q)), # ped length + Qs have just been added
        Sex = 3,
        Mate = 0,
        Mating_type = 0, # don't know what it's supposed to code for
        yob_mod = as.numeric(gsub(".*_", "", New_virgin_Qs$ID)) + 1,
        seq_dam = seq_1, # as their dams, offpring Qs, still need to be ordered and obtain a seq_2, we give them their dam's seq_1 for now
        # seq_sire = Mates_list, # should be mate of dam, can be updated when controlled mating happens
        seq_sire = 1, # should be mate of dam, but we don't know it yet, except for Wg of open mated Qs when OM is the dummy open mating sire
        
        c8 = 0,
        c9 = 0,
        c10 = 0,
        Test_station = 0,
        c12 = 0,
        c13 = 0,
        c14 = 0,
        c15 = 0,
        seq_1 = seq((Q_index+NB_W1_Q+1), (Q_index+NB_W1_Q+NB_W1_Q)),
        # seq_2 = seq((Q_index+1), (Q_index+NB_W1_Q)), # will be given by function 3)	order_off_give_seq2
        yob = as.numeric(gsub(".*_", "", New_virgin_Qs$ID)),
        NS = 0,
        ND = 0,
        # The ID start with 'Q_' or 'TQ_'. In both case, we change if for 'Wg_Q' or 'Wg_TQ'
        # ID = gsub("Q", "W_Q", New_virgin_Qs_ped$ID),
        ID = gsub("(Q|TQ)", "Wg\\1", New_virgin_Qs$ID), # We match either Q or TQ, and append before Wg, depending on what was matched in the firt pattern
        ID_BQ = New_virgin_Qs$ID,
        is_SS = NA,
      )
  )
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              order_off_give_seq2                                    #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# orders dam's offspring (their Wg, offspring Qs and TQs), gives them a seq_2 base on the new order, and updates the seq_dam seq_sire of the offspring Qs' Wgs
# After seq2s are updated, we can update blocks table: first col is the first seq2 of a block, second column is last seq2 of that block

order_off_give_seq2 <- function(pedigree, birth_year_Q_offspring=year, NB_Wg_offspring=NB_W1_Q+NB_W1_TQ, ID_1b=UNKNOWN_1b){
  # because a seq_2 has to be given to Wgs, Qs, and TQs, so that if they have the same dam, they are together, we order them base on dam seq 
  # and give them a seq2 once all offspring of a dam have been born and added to the pedigree
  
  ped_length <- nrow(pedigree)
  
  ##########################################################
  # Sorting 
  ##########################################################
  # ordering dam's offspring (Wgs, Qs, and TQs) based on dam's seq (could be done only to newly added indiv to not order all the dataset for efficiency)
  # Difficulty:
  # unknown sires (when UNKNOWN_1b == '0') have 0 as seq_dam, so we first sort by yob_mod (and not yob so that Wgroup cluster together with offspring Qs born the year after) to keep them at the right year position
  # 0s sires have a seq_dam == 0, and offspring of a queen can be born in different years.
  # We need sorting of yob_mod to get 0 sires correctly place relative to others, and then offspring of a queen (Wgroups and Qs) sorted on seq_dam.
  # to do so, we will modify a yob_mod col to have all offspring of a queen with a same value
  pedigree <- pedigree %>%
    # We create a new column for the adjusted yob_mod based on the seq_dam grouping
    group_by(seq_dam) %>%
    mutate(yob_mod_adj = if_else(seq_dam != 0, min(yob_mod), yob_mod)) %>% # we put all offspring of a queen at a same yob_mod_adj to cluster them together in blocks
    ungroup() %>%
    # We sort the data frame by the new yob_mod_adj and seq_dam
    arrange(yob_mod_adj, seq_dam) %>%
    # Clezning up
    select(-yob_mod_adj) %>%
    mutate(seq_2 = seq(1, ped_length)) # would updating only new Wgs, Qs, and TQs' seq_2 save meaningful computation time? 
  
  ########################################################## BEFORE
  # Wgs of dams, New Wgs, Qs, TQs, and possibly 0sires have had their seq_2 possibly changed.
  # We will update the seq_2 in Mate of new Qs and TQs so that we can use it as seq_2 in seq_sire column for Wgroups
  
  ########################################################## NEW
  pedigree_new_offspring <- pedigree %>%
    # filter(yob %in% c(birth_year_Q_offspring)) # Wgs of dams (ordered in cluster with dam's queen offspring), offspring Qs, offspring TQs, 0sires
    filter(yob %in% birth_year_Q_offspring) # For sake of simplicity, we will update all Mates, bur for founders, which are certainly not wrong, and will preserve order of columns
  ## updating the seq_2 in Mate for Qs and TQs
  
  # For 0sires
  ########## Tables with 0sires with correct seq_2 ti put in Mates of pedigree
  if(ID_1b == '0'){
    reference_Mates <- pedigree_new_offspring %>%
      filter(Sex %in% 2, seq_dam %in% 0) %>% # we filter out 0sires
      select(seq_2, ID_mated_Q) %>%
      rename(Mate=seq_2)
    ########## Updating
    pedigree <- pedigree %>%
      left_join(reference_Mates, by = c("ID"="ID_mated_Q"), suffix = c("", ".updated")) %>%
      mutate(Mate = ifelse(!is.na(Mate.updated), Mate.updated, Mate)) %>% # if there was a Mate already, we keep it (0 or initial Mates), otherwise (if it was NA), we update it
      select(-Mate.updated)
  }
  
  # For controlled mated sires (SS only implemented so far)
  if("ID_TQ_mated" %in% colnames(pedigree_new_offspring)){
    ########## Tables with controlled mated sires with correct seq_2 to put in Mates of pedigree
    # We start by getting a table will all controlled mated Qs and the ID_TQ_mated they were mated to
    reference_Mates <- pedigree_new_offspring %>%
      filter(!(is.na(ID_TQ_mated))) %>% # 
      select(ID, ID_TQ_mated)
    # Next, we need to add to this table the seq_2 of the ID_TQ_mated in reference_Mates
    ## We first get the year of birth of the TQs
    birth_year_TQMates <- as.numeric(gsub(".*_", "", reference_Mates$ID_TQ_mated))# We match 'anything before and comprising '_' in ID_TQ_mated
    # Subtable containing among other the SS TQs and their seq_2
    seq2TQ <- pedigree %>%
      filter(yob %in% birth_year_TQMates) %>%
      select(ID, seq_2) %>%
      rename(ID_TQ_mated=ID)
    # We use the subtable seq2TQ to add the seq_2 od TQs to our reference_Mates table
    reference_Mates_seq2TQ <- reference_Mates %>%
      left_join(seq2TQ, 
                by="ID_TQ_mated") %>%
      select(-ID_TQ_mated) %>% # cleaning this column now that we have the seq_2, so as to now double that column while updating
      rename(Mate=seq_2) # For the following update
    
    ########## Updating the Mate of controlled mated queens
    pedigree <- pedigree %>%
      left_join(reference_Mates_seq2TQ, by = "ID", suffix = c("", ".updated")) %>%
      mutate(Mate = ifelse(!is.na(Mate.updated), Mate.updated, Mate)) %>% # if there was a Mate already, we keep it (0 or initial Mates), otherwise (if it was NA), we update it
      select(-Mate.updated)
  }
  
  ##########################################################
  # Seq_dam and seq_sire of Wgroup with new seq_2 
  ##########################################################
  # giving the correct seq_2 in seq_dam and seq_sire of Wgs
  ## sorting out the Wgs
  ## Because Wgroups are always added first, before other Qs offspring, they are alone with their high yob_mod and appear at the end of the pedigree
  ped_Wgs <- pedigree[(ped_length-NB_Wg_offspring+1):ped_length, ]
  ## getting the seq_2 of the dam using col ID_BQ of ped_Wgs and col ID of pedigree
  ped_Wgs <- ped_Wgs %>%
    select(-c(seq_dam, seq_sire)) %>% # cleaning out the columns with temporary seq_dam and seq_sire values
    left_join(pedigree %>%
                # we will only look in the pedigree with individuals born the year of the Wgroups' dams
                filter(yob_mod %in% birth_year_Q_offspring) %>% # saves computation time if pedigree becomes very long? Could be slices also, knowing which position to match, but less friendly to understand the code
                select(ID, seq_2, Mate) %>%
                rename(seq_dam=seq_2, seq_sire=Mate), by=c("ID_BQ"="ID"))
  ## updating the seq_dam of Wgs in pedigree based on the seq_dam in ped_Wgs
  pedigree <- pedigree[-((ped_length-NB_Wg_offspring+1):ped_length), ]
  pedigree <- pedigree %>% bind_rows(ped_Wgs) # we could have also not added Wgs to pedigree, and only add it oince their Qs have gotten a seq2, instead of adding them right away and manipulating the pedigree part about them
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              create_blocks                                          #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# creates blocks table: first col is the first seq2 of a block, second column is last seq2 of that block
# blocks are defined by the fact that the seq2 of the dam is the same
# blocks are used by the R or F90 program to produce the D matrix: a block-diagonal matrix containing mendelian sampling variances for FS
create_blocks <- function(pedigree){
  # blocks are defined by the fact that the seq2 of the dam is the same
  blocks <- pedigree %>%
    group_by(seq_dam) %>%
    filter(seq_dam != 0,# indiv from an unknown dam (therefore also sire) are not into blocks
           n() > 1) %>%  # we also filter out seq_dam that appear only once, as they when a queen has only one offspring (its Wg), that offspring does not form a block (blocks are at least of size 2)
    summarise(first_seq2 = min(seq_2), last_seq2 = max(seq_2)) %>%
    select(-seq_dam)
  
  return(blocks)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              write_AINV_input                                       #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
write_AINV_input <- function(pedigree, WD_PATH_BLUP){
  #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  #                             blocks
  # creates blocks table: first col is the first seq2 of a block, second column is last seq2 of that block
  # blocks are defined by the fact that the seq2 of the dam is the same
  # blocks are used by the R or F90 program to produce the D matrix: a block-diagonal matrix containing mendelian sampling variances for FS
  # blocks are defined by the fact that the seq2 of the dam is the same
  blocks <- pedigree %>%
    group_by(seq_dam) %>%
    filter(seq_dam != 0,# indiv from an unknown dam (therefore also sire) are not into blocks
           n() > 1) %>%  # we also filter out seq_dam that appear only once, as they when a queen has only one offspring (its Wg), that offspring does not form a block (blocks are at least of size 2)
    summarise(first_seq2 = min(seq_2), last_seq2 = max(seq_2)) %>%
    select(-seq_dam)
  
  #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  #                             singlesires
  
  # creates singlesires table: all single sires' DPQs' seq2, original ID, ND
  # IS ORDER IMPORTANT?
  singlesires <- pedigree %>%
    filter(is_SS %in% 1) %>% # only singles sires (excludes PSs), could be done with a join instead of %in%?
    select(seq_2, ID)
  
  singlesires$ND <- NB_D_FECUNDATING # ND produced by ss for each mating (not ND of how many drones mated the SS)
  
  #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  #                              ident
  
  # creates ident table: original IDs in seq1 order
  # IS ORDER CORRECT? SHOULD INCLUDE DUMMY PS? OR WILL THEY ALWAYS HAVE A identified ID_1b anyways, so that they are not dummy?
  ident <- pedigree %>%
    arrange(seq_1) %>%
    select(ID)
  
  #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  #                              steer
  # [1] default NS (in base); [2] default ND (in base); [3] check_code; [4] NTOT; [5] ncol (dummy entry); [6] singlesires present (1 or 0); [7] Nb of blocks; [8] equi (1 or 0; are founders in equilibrium);
  # [9] distribution of offspring per DPQ and per D (0: Poisson, 1: equal contrib); [10] pinv (dummy, always 0): invert A based on M, I, D decomposition);
  # [11] priF (always 0, not supported in F90 version: calculate inbreeding from possible future matings); [12] yip (always 0, not supported in F90 version: years of parents to consider for possible mating and therefore off F calculation)
  steer <- c(NS_Q_open, NB_D_FECUNDATING, 0, length(pedigree$seq_2), 0, if(MATING_SCALE == 'Q_x_DPQ'){1}else{0}, length(blocks$first_seq2), Assumed_relationship_founders,
             Assumed_dist_Poisson_EqualContrib_Doff_DPQoff, 0, 0, 0)
  
  #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  #                              writing
  # write the tables as AINV input files
  # We only keep col 1:20 of pedigree, as the rest is for internal use
  fwrite(pedigree %>% select(1:20), paste0(WD_PATH_BLUP, "/pedigree.txt"), col.names = FALSE, row.names = FALSE, sep=" ")
  fwrite(blocks, paste0(WD_PATH_BLUP, "/blocks.txt"), col.names = FALSE, row.names = FALSE, sep=" ")
  fwrite(singlesires, paste0(WD_PATH_BLUP, "/singlesires.txt"), col.names = FALSE, row.names = FALSE, sep=" ")
  fwrite(ident, paste0(WD_PATH_BLUP, "/ident.txt"), col.names = FALSE, row.names = FALSE, sep=" ")
  fwrite(as.list(steer), paste0(WD_PATH_BLUP, "/steer.txt"), col.names = FALSE, row.names = FALSE, sep=" ")
}


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Fill_ped_PS_of_Q                                       #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# NOT USED ANYMORE      ###############################################################

Fill_ped_PS_of_Q <- function(pedigree, New_virgin_Qs=WQ_n){
  
  newly_mated_Qs <- which(pedigree$ID %in% New_virgin_Qs$ID)
  
  pedigree$ID_PS[newly_mated_Qs] <- pedigree$ID_PGM[newly_mated_Qs]
  
  pedigree$ID_PS[newly_mated_Qs] <- gsub("Q", "PS", pedigree$ID_PS[newly_mated_Qs])
  
  # pedigree$ID_PS[newly_mated_Qs] <- gsub(paste0("_", (year-4)), paste0("_", (year-3)), pedigree$ID_PS[newly_mated_Qs])
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Fill_ped_F_Q                                           #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# NOT USED ANYMORE, F IS NOW IN WQ_n
#Constant-size A matrix ###############################################################
Fill_ped_F_Q <- function(pedigree, A, ped_A){
  Nb_rows_pedigree <- nrow(pedigree)
  index_indiv=(Nb_rows_pedigree-NB_W1_Q) # first line of new Qs
  NB_W1_indiv <- length(ped_A$ID)
  pedigree$F_Q[(index_indiv+1):(index_indiv+NB_W1_indiv)] <- diag(A[c(ped_A$ID), c(ped_A$ID)]) - 1
  return (pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           new_pedigree_PS                                                           #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# NOT USED ANYMORE, pedigree does not contain DUMMY INDIV anymore
# As long as PS mates are unknown: fills the pedigree file accordingly
new_pedigree_PS <- function(pedigree, BQ_n, NB_BQ, NB_W1_Q, year){
  #We create a table with PSs which we will bind_rows to the pedigree 
  pedigree_PS <- data.frame(matrix(ncol = 6, nrow = NB_BQ))
  colnames(pedigree_PS) <- c("ID", "ID_BQ", "ID_PS", "ID_PGM","ID_PGM_mated", "F_Q")
  
  pedigree_PS$ID <- BQ_n$ID #each BQ of this year will correspond to one PS
  
  #Retrieving infos about the parents, seen as animals in the pedigree
  PS_info <- pedigree %>% #we only look in lines of pedigree where the BQs could be
    # slice((max(NB_W1_Q+(NB_W1_Q+NB_BQ)*(year-3), 0) +1):(max(NB_W1_Q+(NB_W1_Q+NB_BQ)*(year-3), 0)+NB_W1_Q)) %>%
    filter(ID %in% pedigree_PS$ID)#we store infos about the BQ treated, corresponding to one PS
  
  #using this parents' info to fulfill the animals' infos in the pedigree
  pedigree_PS$ID_BQ <- PS_info$ID
  pedigree_PS$ID_PS <- PS_info$ID_PGM_mated
  pedigree_PS$ID_PGM <- PS_info$ID_PGM_mated
  
  #PS will always be mated to unknown drones
  pedigree_PS$ID_PGM_mated <- as.character(0)
  pedigree_PS$F_Q <- 0 #inbreeding of PSs will not be calculated, but the individual F of DPQs will
  
  #we give PSs an appropriated ID, which is a modification of their dam BQ IDs
  pedigree_PS$ID <- gsub("Q", "PS", pedigree_PS$ID)
  pedigree_PS$ID <- gsub(paste0("_", (year-1)), paste0("_", year), pedigree_PS$ID)
  
  pedigree_PS$ID_PS <- gsub("Q", "PS", pedigree_PS$ID_PS)
  pedigree_PS$ID_PS <- gsub(paste0("_", (year-4)), paste0("_", (year-3)), pedigree_PS$ID_PS)
  
  #We then bind_rows this table to the pedigree:
  pedigree <- bind_rows(pedigree, pedigree_PS)
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Add_unknown_PGM_mated_to_pedigree_Pim                  #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Add_unknown_PGM_mated_to_pedigree_Pim <- function(pedigree, indiv_n, ID_1b='OM'){
  
  newly_mated_TQs <- which(pedigree$ID %in% indiv_n$ID)
  #Qs of these years are fecundated by drones with unknown ascent:
  pedigree$PGM_mated_birth_year[newly_mated_TQs] <- 0
  pedigree$ID_PGM_mated[newly_mated_TQs] <- "0"
  pedigree$Q1b_birth_year[newly_mated_TQs] <- 0
  pedigree$ID_1b[newly_mated_TQs] <- ID_1b
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Add_unknown_PGM_mated_to_indiv_n                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# this function could be avoided by having the info of the pedigree mate in D_n and passing it during mating to WQ_n (function Produce_Wgroup)
Add_unknown_PGM_mated_to_indiv_n <- function(indiv_n=WQ_n, ID_1b='OM'){
  
  indiv_n$ID_PGM_mated <- ID_1b
  
  return(indiv_n)
}


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Add_unknown_PGM_mated_to_pedigree                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Add_unknown_PGM_mated_to_pedigree <- function(pedigree, indiv_n=WTQ_n, NS_Q_or_TQ_open=NS_TQ, NB_D_FECUNDATING){
  
  Opsires_index <- length(pedigree$ID) 
  NB_open_mated_indiv <- length(indiv_n$ID)
  
  # Creating the open mating sires' ped
  pedigree_0sires <- data.frame(matrix(ncol = 22, nrow = NB_open_mated_indiv))
  colnames(pedigree_0sires) <- c("ID_seq1", "Sex", "Mate",
                                 "Mating_type", "yob_mod", "seq_dam", "seq_sire", 
                                 "c8", "c9", "c10", "Test_station", "c12", "c13", "c14", "c15",
                                 "seq_1", "seq_2", "yob", "NS", "ND", "ID", "is_SS")
  
  pedigree_0sires <- pedigree_0sires %>% mutate(
    ID_seq1 = seq(Opsires_index+1, Opsires_index+NB_open_mated_indiv),
    Sex = 2,
    Mate = 0,
    Mating_type = 0, # don't know what it's supposed to code for
    yob_mod = as.numeric(gsub(".*_", "", indiv_n$ID)),
    seq_dam = 0,
    seq_sire = 0,
    
    c8 = 0,
    c9 = 0,
    c10 = 0,
    Test_station = 0,
    c12 = 0,
    c13 = 0,
    c14 = 0,
    c15 = 0,
    seq_1 = seq(Opsires_index+1, Opsires_index+NB_open_mated_indiv),
    seq_2 = seq(Opsires_index+1, Opsires_index+NB_open_mated_indiv),
    yob = as.numeric(gsub(".*_", "", indiv_n$ID)),
    NS = NS_Q_or_TQ_open,
    ND = NB_D_FECUNDATING,
    ID = gsub("(Q|TQ)", "OpSire_\\1", indiv_n$ID), # We match either Q or TQ, and append before OpSire_TQ, depending on what was matched in the firt pattern
    is_SS = NA,
    ID_mated_Q=indiv_n$ID,
    # is_Osire=1 # will be used to filter pedigree lines of 0sires in order_off_give_seq2
  )
  
  # # Updating the new sires as the mate of queens
  # ## Subtable of temporary seq_2 of mates and original ID of open mated Qs
  # ID_mates_of_Qs <- pedigree_0sires %>%
  #   select(seq_2)
  # ID_mates_of_Qs$ID_mates_Q <- indiv_n$ID
  
  
  pedigree <- bind_rows(pedigree, pedigree_0sires)
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Add_PGM_mated_to_pedigree_Pim                            #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Add_PGM_mated_to_pedigree_Pim <- function(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE){
  # Q1b_birth_year is the TQ birth year if SS mating, otherwide PGM birth year +1
  #We can now fill the pedigree acordingly for all queens of the current year:
  #we add the ID of PGM queens, mothers of DPQs constituing the Pseudo_sires to which Qs have just been attributed for mating, in column "ID_PGM_mated", useful for blup estimations.
  newly_mated_Qs <- which(pedigree$ID %in% PS_or_TQ_to_WQ_n$ID_Q)
  
  pedigree$PGM_mated_birth_year[newly_mated_Qs] <- as.numeric(gsub(".*_", "", PS_or_TQ_to_WQ_n$ID_PGM))
  pedigree$ID_PGM_mated[newly_mated_Qs] <- PS_or_TQ_to_WQ_n$ID_PGM # queens need to be orderes the same in the pedigree and the PS_or_TQ_to_WQ_n
  
  # IF single sires are used (pseudo-sires made up of only one Q), the DPQ ID is entered for ID_1b, otherwise a fabricated one from the ID_PGM and NS_Q
  if(MATING_SCALE == "Q_x_DPQ"){
    # we match everything before the "w" in the TQ ID and replace it with nothing to be left with only the birth year
    pedigree$Q1b_birth_year[newly_mated_Qs] <- gsub(".*_", "", PS_or_TQ_to_WQ_n$ID_TQ)
    pedigree$ID_1b[newly_mated_Qs] <- PS_or_TQ_to_WQ_n$ID_TQ
    
  }else if(MATING_SCALE == "Q_x_PS"){
    pedigree$Q1b_birth_year[newly_mated_Qs] <- as.numeric(gsub(".*_", "", PS_or_TQ_to_WQ_n$ID_PGM))+1
    pedigree$ID_1b[newly_mated_Qs] <- paste0("NS", NS, "_", PS_or_TQ_to_WQ_n$ID_PGM)
  }
  
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Add_PGM_mated_to_pedigree                              #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Add_PGM_mated_to_pedigree <- function(pedigree, PS_or_TQ_to_WQ_n, NS=NS_Q, MATING_SCALE){
  # NS useful for future implementation of PS mating.
  
  if(MATING_SCALE == "Q_x_DPQ"){
    # for SS mating, we change the 'Sex' in pedigree to '2' for TQs present in the PS_or_TQ_to_WQ_n PROBLEM: if used as mate, but mated Qs all die, they will not become sires (yes: they will, at least for the Wgs of the mated Q)
    Lines_pedigree_SS <- which(pedigree$ID %in% PS_or_TQ_to_WQ_n$ID_TQ)
    pedigree$Sex[Lines_pedigree_SS] <- 2
    pedigree$is_SS[Lines_pedigree_SS] <- 1 # we mark in pedigree that these sires are SS, so that we can extract them easily and efficiently when creating singlesires with create_singlesires()
    # We update the ID_mated_Q (formerly we updated Mate) of Qs in the pedigree with the ID (formerly the seq_2) of the TQ they were mated to
    ## we start by adding the seq_2 of TQs to PS_or_TQ_to_WQ_n
    ########################## BEFORE
    # PS_or_TQ_to_WQ_n <- PS_or_TQ_to_WQ_n %>%
    #   left_join(pedigree %>%
    #               filter(yob_mod %in% birth_year_mates) %>% # saves computation time?
    #               select(ID, seq_2), by=c("ID_TQ"="ID")) %>%
    #   rename(seq_2_mate=seq_2)
    # ## now that seq_2 is present in PS_or_TQ_to_WQ_n, we can update the Mate in the pedigree: their Mate becomes the seq_2 in PS_or_TQ_to_WQ_n
    # pedigree <- pedigree %>%
    #   left_join(PS_or_TQ_to_WQ_n %>% select(ID_Q, seq_2_mate), by = c("ID" = "ID_Q")) %>%
    #   mutate(Mate = ifelse(!is.na(seq_2_mate), seq_2_mate, Mate)) %>%
    #   select(-seq_2_mate)
    
    ##########################
    # NEW
    
    # Adding ID_TQ_mated to pedigree
    ## When it's the first time we ass controlled mates, there is no col ID_TQ_mate, so we put an if for that case
    if(!("ID_TQ_mated" %in% colnames(pedigree))){
      pedigree <- pedigree %>%
        left_join(PS_or_TQ_to_WQ_n %>% select(ID_Q, ID_TQ) %>%
                    rename(ID_TQ_mated = ID_TQ), by = c("ID" = "ID_Q"))
    }else{
      pedigree <- pedigree %>%
        left_join(PS_or_TQ_to_WQ_n %>% select(ID_Q, ID_TQ) %>%
                    rename(ID_TQ_mated = ID_TQ), by = c("ID" = "ID_Q"), suffix = c("", ".updated")) %>%
        mutate(ID_TQ_mated = ifelse(!is.na(ID_TQ_mated.updated), ID_TQ_mated.updated, ID_TQ_mated)) %>% # if there was a ID_TQ_mated already, we keep it (0 or initial Mates), otherwise (if it was NA), we update it
        select(-ID_TQ_mated.updated)
    }
    # pedigree$Mate[which(pedigree$ID %in% PS_or_TQ_to_WQ_n$ID_Q)] <- PS_or_TQ_to_WQ_n$seq_2 # works only if order of Qs in pedigree and PS_or_TQ_to_WQ_n is the same
    
  }else if(MATING_SCALE == "Q_x_PS"){
    # For PS mating, we first need to add pseudo-sires to the pedigree, with correct seq_dam
  }
  
  
  
  
  return(pedigree)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                                         BLUP
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Create_Perf_Bras                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Create_Perf_Bras <- function(Output_Q, Output_TQ, erase_Perf_T1_Q=FALSE, erase_Perf_T2_Q=TRUE, erase_Perf_T1_TQ=FALSE, erase_Perf_T2_TQ=FALSE,
                             year_erase_Perf_T1_Q=year, year_erase_Perf_T2_Q=year, year_erase_Perf_T1_TQ=year, year_erase_Perf_T2_TQ=year){
  # Generates a first perf file for BLUP prediction
  
  # For Qs
  perf_Bras <- Output_Q %>%
    mutate(
      ID_Wgroup = paste0("W_", ID),
      gen_mean_fix=1
    ) %>%
    select(ID_Wgroup, ID, gen_mean_fix, Beekeeper_year, Performance_year, Perf_T1, Perf_T2)
  # %>%
  #   mutate(Beekeeper_year=paste0(Beekeeper, "_", Performance_year))
  if(erase_Perf_T1_Q){
    perf_Bras <- perf_Bras %>% mutate(Perf_T1 = ifelse(Performance_year %in% eval(parse(text=year_erase_Perf_T1_Q)), -999, Perf_T1)) # potential BQs can have their perf erased to missing value to simulate partial phenotyping
  }
  if(erase_Perf_T2_Q){
    perf_Bras <- perf_Bras %>% mutate(Perf_T2 = ifelse(Performance_year %in% eval(parse(text=year_erase_Perf_T2_Q)), -999, Perf_T2)) # potential BQs can have their perf erased to missing value to simulate partial phenotyping
  }
  
  
  # For TQs
  perf_Bras_TQ <- Output_TQ %>%
    mutate(
      ID_Wgroup = paste0("W_", ID),
      gen_mean_fix=1
    ) %>%
    select(ID_Wgroup, ID, gen_mean_fix, Beekeeper_year, Performance_year, Perf_T1, Perf_T2)
  # %>%
  #   mutate(Beekeeper_year=paste0(Beekeeper, "_", Performance_year))
  if(erase_Perf_T1_TQ){
    perf_Bras_TQ <- perf_Bras_TQ %>% mutate(Perf_T1 = ifelse(Performance_year %in% eval(parse(text=year_erase_Perf_T1_TQ)), -999, Perf_T1)) # potential DPQs can have their perf erased to missing value to simulate partial phenotyping
  }
  if(erase_Perf_T2_TQ){
    perf_Bras_TQ <- perf_Bras_TQ %>% mutate(Perf_T2 = ifelse(Performance_year %in% eval(parse(text=year_erase_Perf_T2_TQ)), -999, Perf_T2)) # potential DPQs can have their perf erased to missing value to simulate partial phenotyping
  }
  
  # BQs and TQs together
  perf_Bras <- rbind(perf_Bras, perf_Bras_TQ) #homogenizing
  
  # recoding the level of the Beekeeper_year to have it in integer form
  Beekeeper_year_recode <- data.table(Beekeeper_year_ini=unique(perf_Bras$Beekeeper_year), Beekeeper_year_recode=seq(1, length(unique(perf_Bras$Beekeeper_year)) ) )
  perf_Bras <- perf_Bras %>% left_join(y=Beekeeper_year_recode, by=c("Beekeeper_year"="Beekeeper_year_ini"))
  
  # Eliminating column non-recoded and reordering
  perf_Bras <- perf_Bras %>% select(ID_Wgroup, ID, gen_mean_fix, Beekeeper_year_recode, Perf_T1, Perf_T2, Beekeeper_year)
  
  return(perf_Bras)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              perf_Bras_RenumIDsAINV                                 #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
perf_Bras_RenumIDsAINV <- function(WD_PATH_BLUP, perf_Bras){
  
  pedigree_complete <- fread(paste0(WD_PATH_BLUP, "/pedigree_complete_f.txt"), select=c("V6", "V17", "V1", "V2"))
  pedigree_complete <- pedigree_complete %>% rename (seq_dam=V6, seq_2=V17, original_ID=V1, sex=V2)
  
  # Joining columns with the queen and colony IDs as created by Brascamp programs present in pedigree_complete to Output file
  # For original queen IDs
  # Join: dam in Output (the queen in a colony with a performance) has to match the original_ID in pedigree_complete --> Adds the original ID and seq_2 of queens to Output
  perf_Bras_joinID <- perf_Bras %>% left_join(y=select(pedigree_complete, original_ID, seq_2), by=c("ID" = "original_ID")) %>% 
    rename(seq_Q=seq_2)
  
  
  # For Wgroups IDs (sex=3)
  # Join: 
  # left side: perf file, 
  # right side: pedigree_complete with only Wgroups, seq_dam and seq_2: ID of the queen of that Wgroup and of the Wgroup itsel
  # Matching: a queen in the perf file with a queen of a Wgroup in pedigree_complete --> Adds the seq_2 of Wgroups to the perf file
  # seq_2 or Wgroups are renamed seq_col
  perf_Bras_joinID <- perf_Bras_joinID %>% left_join(y= pedigree_complete %>% filter(sex %in% 3) %>% select(seq_dam, seq_2), by=c("seq_Q" = "seq_dam")) %>% 
    rename(seq_Wgroup=seq_2)
  
  
  perf_Bras_joinID <- perf_Bras_joinID %>% select(seq_Wgroup, seq_Q, 
                                                  gen_mean_fix, Beekeeper_year_recode,
                                                  Perf_T1, Perf_T2,
                                                  ID_Wgroup, ID, Beekeeper_year)# %>%
  # rename(Wgroup_OrID=ID_Wgroup,
  #        Queen_OrID=ID)
  
  
  # We recode exact 0 values with -9 (coding NA for AIREMF90) for trait 1 or 2: when it is not simulated
  # perf_Bras_joinID[is.null(perf_Bras_joinID)] <- -999 # BEWARE! Should recode 0 records for traits not phenotypes to -999
  perf_Bras_joinID$Perf_T1 <- gsub("^0$", "-999", perf_Bras_joinID$Perf_T1) # '^': means 'beginning'; '$': means 'end'
  perf_Bras_joinID$Perf_T2 <- gsub("^0$", "-999", perf_Bras_joinID$Perf_T2) # '^': means 'beginning'; '$': means 'end'
  
  return(perf_Bras_joinID)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Retrieve_BLUP_sol                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Retrieve_BLUP_sol <- function(WD_PATH_BLUP, perf_Bras_joinID, SINGLE_TRAIT, CandidateBQ_n=WQ_n, CandidateDPQ_n=WTQ_n){
  ##### Reading log file to retrieve the convergence success or not and the estimated parameters
  if(SINGLE_TRAIT){
    # reading log file
    AIREML_log <- fread(paste0(WD_PATH_BLUP, "/aireml_T1_CM_corr1.log"), fill=TRUE)
    # Convergence criterion: retaining the value in the second line, column V6 of AIRELM log file
    Conv_crit <- AIREML_log[2, 6] %>% pull(V6)
    # Convergence round: retaining the value in the second line, column V3 of AIRELM log file
    Conv_round <- AIREML_log[2, 3] %>% pull(V3)
    # estimated (CO)VARIANCES
    ## Variances
    E_VAR_A1_D <- AIREML_log[6, 1] %>% pull(V1)
    E_VAR_A1_M <- AIREML_log[7, 2] %>% pull(V2)
    ## Genetic covariances between effects
    E_COV_A1_D_M <- AIREML_log[6, 2] %>% pull(V2)
    
    ## Residual variance
    E_VAR_E1 <- AIREML_log[16, 1] %>% pull(V1)
  }else{
    # reading log file
    AIREML_log <- fread(paste0(WD_PATH_BLUP, "/aireml_T1_T2_CM_corr1.log"), fill=TRUE)
    # Convergence criterion: retaining the value in the second line, column V6 of AIRELM log file
    Conv_crit <- AIREML_log[2, 6] %>% pull(V6)
    # Convergence round: retaining the value in the second line, column V3 of AIRELM log file
    Conv_round <- AIREML_log[2, 3] %>% pull(V3)
    # estimated (CO)VARIANCES
    ## Genetic variances
    E_VAR_A1_D <- AIREML_log[6, 1] %>% pull(V1)
    E_VAR_A2_D <- AIREML_log[7, 2] %>% pull(V2)
    E_VAR_A1_M <- AIREML_log[8, 3] %>% pull(V3)
    E_VAR_A2_M <- AIREML_log[9, 4] %>% pull(V4)
    ## Genetic covariance among each of both traits
    E_COV_A1_D_M <- AIREML_log[6, 3] %>% pull(V3)
    E_COV_A2_D_M <- AIREML_log[7, 4] %>% pull(V4) #direct-maternal genetic covariance
    ## Genetic covariance between each of both traits
    E_COV_A1_A2_D <- AIREML_log[6, 2] %>% pull(V2)
    E_COV_A1_A2_M <- AIREML_log[8, 4] %>% pull(V4)
    E_COV_A1_D_A2_M <- AIREML_log[6, 4] %>% pull(V4) #direct-maternal genetic covariance between traits
    E_COV_A1_M_A2_D <- AIREML_log[7, 3] %>% pull(V3) #direct-maternal genetic covariance between traits
    
    ## Residual variances
    E_VAR_E1 <- AIREML_log[23, 1] %>% pull(V1)
    E_VAR_E2 <- AIREML_log[24, 2] %>% pull(V2)
    ## Residual covariance between traits
    E_COV_E1_E2 <- AIREML_log[23, 2] %>% pull(V2)
  }
  
  ##### Reading sol file and IDs table
  Solutions <- fread(paste0(WD_PATH_BLUP,"/solutions"), header = FALSE, select=c(1:5), col.names=c("trait", "effect", "level", "solution", "SE")) 
  # pedigree_complete <- fread(paste0(WD_PATH_BLUP, "/pedigree_complete_f.txt"), select=c("V6", "V17", "V1", "V2"))
  # pedigree_complete <- pedigree_complete %>% rename (seq_dam=V6, seq_2=V17, original_ID=V1, sex=V2)
  
  ##### Worker effect
  # Keeping only EBVs of worker effects (direct effects)
  EBVs_d <- Solutions %>%
    filter(effect %in% c(1)) %>%
    select (-effect) #we only keep animal levels, descarting solutions of non-animal effects. Effect 1 is worker effects, 2 queen effects
  
  if(SINGLE_TRAIT == FALSE){
    # Widening the table to obtain one line per animal (i.e. per 'level')
    EBVs_d_wide <- EBVs_d %>% pivot_wider(names_from = trait, values_from = c("solution", "SE")) %>%
      rename(ID_renum=level, EBV_d1=solution_1, EBV_d2=solution_2, SE_EBV_d1=SE_1, SE_EBV_d2=SE_2)
  }else{ # if there is only solution fopr trait 1
    EBVs_d_wide <- EBVs_d %>% pivot_wider(names_from = trait, values_from = c("solution", "SE")) %>%
      rename(ID_renum=level, EBV_d1=solution_1, SE_EBV_d1=SE_1)
  }
  ##### Queen effect
  # Keeping only EBVs of queen effects (maternal effects)
  EBVs_m <- Solutions %>%
    filter(effect %in% c(2)) %>%
    select (-effect) #we only keep animal levels, descarting solutions of non-animal effects. Effect 1 is worker effects, 2 queen effects
  
  if(SINGLE_TRAIT == FALSE){
    # Widening the table to obtain one line per animal (i.e. per 'level')
    EBVs_m_wide <- EBVs_m %>% pivot_wider(names_from = trait, values_from = c("solution", "SE")) %>%
      rename(EBV_m1=solution_1, EBV_m2=solution_2, SE_EBV_m1=SE_1, SE_EBV_m2=SE_2) 
  }else{ # if there is only solution fopr trait 1
    EBVs_m_wide <- EBVs_m %>% pivot_wider(names_from = trait, values_from = c("solution", "SE")) %>%
      rename(EBV_m1=solution_1, SE_EBV_m1=SE_1) 
  }
  ##### Direct an queen effect together
  EBVs_d_m <- bind_cols(EBVs_d_wide, EBVs_m_wide %>% select(-level))
  
  ##### Getting original IDs
  # WITH PEDIGREE COMPLETE: COULD NOT MAKE IT WORK, USING PERF FILE INSTEAD
  # EBVs_d_m_test <- EBVs_d_m %>% left_join(pedigree_complete %>% select(seq_2, original_ID), by=c("ID_renum"="seq_2")) %>% rename(ID_Q_original=original_ID)
  # 
  # EBVs_d_m_test_w <- EBVs_d_m_test %>% left_join(pedigree_complete %>% filter(sex %in% 3) %>% select(seq_dam, seq_2), by=c("ID_Q_original" = "seq_dam")) %>% 
  #   rename(seq_Wgroup=seq_2)
  
  # # Getting original IDs of queens 
  # EBVs_d_m_Q_test <- EBVs_d_m %>% left_join(pedigree_complete %>% select(original_ID, seq_2), by=c("ID_renum"="seq_2")) %>%
  #   filter(!(original_ID %in% NA)) %>%
  #   rename(EBV_d1_Q=EBV_d1, EBV_d2_Q=EBV_d2, SE_EBV_d1_Q=SE_EBV_d1, SE_EBV_d2_Q=SE_EBV_d2, EBV_m1_Q=EBV_m1, EBV_m2_Q=EBV_m2, SE_EBV_m1_Q=SE_EBV_m1, SE_EBV_m2_Q=SE_EBV_m2)
  # 
  # # Getting original IDs of worker groups
  # EBVs_d_m_W_test <- EBVs_d_m %>% left_join(pedigree_complete %>% select(original_ID, seq_dam), by=c("ID_renum"="seq_dam")) %>%
  #   filter(!(original_ID %in% NA)) %>%
  #   rename(EBV_d1_W=EBV_d1, EBV_d2_W=EBV_d2, SE_EBV_d1_W=SE_EBV_d1, SE_EBV_d2_W=SE_EBV_d2, EBV_m1_W=EBV_m1, EBV_m2_W=EBV_m2, SE_EBV_m1_W=SE_EBV_m1, SE_EBV_m2_W=SE_EBV_m2)
  
  
  
  if(SINGLE_TRAIT == FALSE){
    # EBVs_d_m will only contain EBVs for queens and Wgroups from colonies present in the performance file with an original (before recodeing) ID
    # Getting original IDs of queens 
    EBVs_d_m_Q <- EBVs_d_m %>% left_join(perf_Bras_joinID %>% select(ID, seq_Q), by=c("ID_renum"="seq_Q")) %>%
      filter(!(ID %in% NA)) %>%
      select(-ID_renum) %>%
      rename(EBV_d1_Q=EBV_d1, EBV_d2_Q=EBV_d2, SE_EBV_d1_Q=SE_EBV_d1, SE_EBV_d2_Q=SE_EBV_d2, EBV_m1_Q=EBV_m1, EBV_m2_Q=EBV_m2, SE_EBV_m1_Q=SE_EBV_m1, SE_EBV_m2_Q=SE_EBV_m2)
    
    # Getting original IDs of worker groups
    EBVs_d_m_W <- EBVs_d_m %>% left_join(perf_Bras_joinID %>% select(ID, seq_Wgroup), by=c("ID_renum"="seq_Wgroup")) %>%
      filter(!(ID %in% NA)) %>%
      select(-ID_renum) %>%
      rename(EBV_d1_W=EBV_d1, EBV_d2_W=EBV_d2, SE_EBV_d1_W=SE_EBV_d1, SE_EBV_d2_W=SE_EBV_d2, EBV_m1_W=EBV_m1, EBV_m2_W=EBV_m2, SE_EBV_m1_W=SE_EBV_m1, SE_EBV_m2_W=SE_EBV_m2)
    
    # Combining EBVs of Qs and Wgroups in a same table
    EBVs_d_m <- EBVs_d_m_Q %>% left_join(EBVs_d_m_W, by="ID") %>%
      select(ID, EBV_d1_Q, EBV_m1_Q, EBV_d2_Q, EBV_m2_Q, EBV_d1_W, EBV_m1_W, EBV_d2_W, EBV_m2_W,
             SE_EBV_d1_Q, SE_EBV_m1_Q, SE_EBV_d2_Q, SE_EBV_m2_Q, SE_EBV_d1_W, SE_EBV_m1_W, SE_EBV_d2_W, SE_EBV_m2_W)
  }else{ # if there is only solution fopr trait 1
    # EBVs_d_m will only contain EBVs for queens and Wgroups from colonies present in the performance file with an original (before recodeing) ID
    # Getting original IDs of queens 
    EBVs_d_m_Q <- EBVs_d_m %>% left_join(perf_Bras_joinID %>% select(ID, seq_Q), by=c("ID_renum"="seq_Q")) %>%
      filter(!(ID %in% NA)) %>%
      select(-ID_renum) %>%
      rename(EBV_d1_Q=EBV_d1, SE_EBV_d1_Q=SE_EBV_d1, EBV_m1_Q=EBV_m1, SE_EBV_m1_Q=SE_EBV_m1)
    
    # Getting original IDs of worker groups
    EBVs_d_m_W <- EBVs_d_m %>% left_join(perf_Bras_joinID %>% select(ID, seq_Wgroup), by=c("ID_renum"="seq_Wgroup")) %>%
      filter(!(ID %in% NA)) %>%
      select(-ID_renum) %>%
      rename(EBV_d1_W=EBV_d1, SE_EBV_d1_W=SE_EBV_d1, EBV_m1_W=EBV_m1, SE_EBV_m1_W=SE_EBV_m1)
    
    # Combining EBVs of Qs and Wgroups in a same table
    # WOULD IT be faster to just reorder tables based on IDs and then col bind both tables?
    EBVs_d_m <- EBVs_d_m_Q %>% left_join(EBVs_d_m_W, by="ID") %>%
      select(ID, EBV_d1_Q, EBV_m1_Q, EBV_d1_W, EBV_m1_W,
             SE_EBV_d1_Q, SE_EBV_m1_Q, SE_EBV_d1_W, SE_EBV_m1_W)
  }
  
  ##### Apiary effects
  # Getting the apiary effects
  Apiary_effects <- Solutions %>%
    filter(effect %in% c(4)) %>%
    select (-effect) #we only keep animal levels, descarting solutions of non-animal effects. Effect 1 is worker effects, 2 queen effects
  
  # pivoting the table to obtain one line per apiary (i.e. per 'level')
  if(SINGLE_TRAIT == FALSE){
    # Widening the table to obtain one line per apiary (i.e. per 'level')
    Apiary_effects_wide <- Apiary_effects %>% pivot_wider(names_from = trait, values_from = c("solution", "SE")) %>%
      rename(Beekeeper_renum=level, E_Beekeeper_x_year_effect_T1=solution_1, E_Beekeeper_x_year_effect_T2=solution_2, SE_Beekeeper_x_year_effect_T1=SE_1, SE_Beekeeper_x_year_effect_T2=SE_2) %>% 
      # joining original apiary effect coding
      left_join(perf_Bras_joinID %>% select(Beekeeper_year, Beekeeper_year_recode) %>% distinct(), by=c("Beekeeper_renum"="Beekeeper_year_recode"))
    
  }else{ # if there is only solution fopr trait 1
    # Widening the table to obtain one line per apiary (i.e. per 'level')
    Apiary_effects_wide <- Apiary_effects %>% pivot_wider(names_from = trait, values_from = c("solution", "SE")) %>%
      rename(Beekeeper_renum=level, E_Beekeeper_x_year_effect_T1=solution_1, SE_Beekeeper_x_year_effect_T1=SE_1) %>% 
      # joining original apiary effect coding
      left_join(perf_Bras_joinID %>% select(Beekeeper_year, Beekeeper_year_recode) %>% distinct(), by=c("Beekeeper_renum"="Beekeeper_year_recode"))
  }
  
  
  # CandidateBQ_n
  ## Joining these EBVs to CandidateBQ_n
  CandidateBQ_n <- CandidateBQ_n %>% left_join(EBVs_d_m, by = "ID")
  ## joining the apiary effect solutions to 
  CandidateBQ_n <- CandidateBQ_n %>% left_join(Apiary_effects_wide %>% select(-Beekeeper_renum), by = "Beekeeper_year")
  
  ## Adding the convergence criterion and Nb rounds needed to converge
  CandidateBQ_n$Conv_crit <- Conv_crit
  CandidateBQ_n$Conv_round  <- Conv_round
  ## Adding estimated (co)variances
  if(SINGLE_TRAIT){
    CandidateBQ_n$E_VAR_A1_D <- E_VAR_A1_D
    CandidateBQ_n$E_VAR_A1_M <- E_VAR_A1_M
    CandidateBQ_n$E_COV_A1_D_M <- E_COV_A1_D_M
    CandidateBQ_n$E_VAR_E1 <- E_VAR_E1
  }else{
    CandidateBQ_n$E_VAR_A1_D <- E_VAR_A1_D
    CandidateBQ_n$E_VAR_A1_M <- E_VAR_A1_M
    CandidateBQ_n$E_VAR_A2_D <- E_VAR_A2_D
    CandidateBQ_n$E_VAR_A2_M <- E_VAR_A2_M
    CandidateBQ_n$E_COV_A1_D_M <- E_COV_A1_D_M
    CandidateBQ_n$E_COV_A2_D_M <- E_COV_A2_D_M
    CandidateBQ_n$E_COV_A1_A2_D <- E_COV_A1_A2_D
    CandidateBQ_n$E_COV_A1_A2_M <- E_COV_A1_A2_M
    CandidateBQ_n$E_COV_A1_D_A2_M <- E_COV_A1_D_A2_M
    CandidateBQ_n$E_COV_A1_M_A2_D <- E_COV_A1_M_A2_D
    CandidateBQ_n$E_VAR_E1 <- E_VAR_E1
    CandidateBQ_n$E_VAR_E2 <- E_VAR_E2
    CandidateBQ_n$E_COV_E1_E2 <- E_COV_E1_E2
  }
  # CandidateDPQ_n
  ## Joining these EBVs to CandidateDPQ_n
  CandidateDPQ_n <- CandidateDPQ_n %>% left_join(EBVs_d_m, by = "ID")
  ## joining the apiary effect solutions to 
  CandidateDPQ_n <- CandidateDPQ_n %>% left_join(Apiary_effects_wide %>% select(-Beekeeper_renum), by = c("Beekeeper_year"))
  ## Adding the convergence criterion
  CandidateDPQ_n$Conv_crit <- Conv_crit
  CandidateDPQ_n$Conv_round <- Conv_round
  ## Adding estimated (co)variances
  if(SINGLE_TRAIT){
    CandidateDPQ_n$E_VAR_A1_D <- E_VAR_A1_D
    CandidateDPQ_n$E_VAR_A1_M <- E_VAR_A1_M
    CandidateDPQ_n$E_COV_A1_D_M <- E_COV_A1_D_M
    CandidateDPQ_n$E_VAR_E1 <- E_VAR_E1
  }else{
    CandidateDPQ_n$E_VAR_A1_D <- E_VAR_A1_D
    CandidateDPQ_n$E_VAR_A1_M <- E_VAR_A1_M
    CandidateDPQ_n$E_VAR_A2_D <- E_VAR_A2_D
    CandidateDPQ_n$E_VAR_A2_M <- E_VAR_A2_M
    CandidateDPQ_n$E_COV_A1_D_M <- E_COV_A1_D_M
    CandidateDPQ_n$E_COV_A2_D_M <- E_COV_A2_D_M
    CandidateDPQ_n$E_COV_A1_A2_D <- E_COV_A1_A2_D
    CandidateDPQ_n$E_COV_A1_A2_M <- E_COV_A1_A2_M
    CandidateDPQ_n$E_COV_A1_D_A2_M <- E_COV_A1_D_A2_M
    CandidateDPQ_n$E_COV_A1_M_A2_D <- E_COV_A1_M_A2_D
    CandidateDPQ_n$E_VAR_E1 <- E_VAR_E1
    CandidateDPQ_n$E_VAR_E2 <- E_VAR_E2
    CandidateDPQ_n$E_COV_E1_E2 <- E_COV_E1_E2
  }
  ## Putting both tables into a list to be returned by function
  CandidateBQ_n_DPQ_n <- list(CandidateBQ_n, CandidateDPQ_n)
  
  return(CandidateBQ_n_DPQ_n)
}
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Join_EBVs_W1_to_W2_TQ_n                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Join_EBVs_Q_n <- function(Table_to_update_n, CandidateDPQ_n, SINGLE_TRAIT){
  if(SINGLE_TRAIT == FALSE){
    Table_to_update_n <- Table_to_update_n %>% left_join(CandidateDPQ_n %>% select(ID,
                                                                                   EBV_d1_Q, EBV_m1_Q, EBV_d2_Q, EBV_m2_Q,
                                                                                   EBV_d1_W, EBV_m1_W, EBV_d2_W, EBV_m2_W,
                                                                                   SE_EBV_d1_Q, SE_EBV_m1_Q, SE_EBV_d2_Q, SE_EBV_m2_Q,
                                                                                   SE_EBV_d1_W, SE_EBV_m1_W, SE_EBV_d2_W, SE_EBV_m2_W
                                                                                   # Conv_crit, Conv_round,
                                                                                   # E_VAR_A1_D, E_VAR_A1_M, E_VAR_A2_D, E_VAR_A2_M, 
                                                                                   # E_COV_A1_D_M, E_COV_A2_D_M,
                                                                                   # E_COV_A1_A2_D, E_COV_A1_A2_M,
                                                                                   # E_COV_A1_D_A2_M, E_COV_A1_M_A2_D,
                                                                                   # E_VAR_E1, E_VAR_E2, E_COV_E1_E2
    ),
    by="ID")
  }else{ # if single trait sim
    Table_to_update_n <- Table_to_update_n %>% left_join(CandidateDPQ_n %>% select(ID,
                                                                                   EBV_d1_Q, EBV_m1_Q,
                                                                                   EBV_d1_W, EBV_m1_W,
                                                                                   SE_EBV_d1_Q, SE_EBV_m1_Q,
                                                                                   SE_EBV_d1_W, SE_EBV_m1_W
                                                                                   # Conv_crit, Conv_round,
                                                                                   # E_VAR_A1_D, E_VAR_A1_M, E_COV_A1_D_M, 
                                                                                   # E_VAR_E1
    ),
    by="ID")
  }
  return(Table_to_update_n)
}
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#                              Join_EBVs_to_Output                                     #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

Join_EBVs_to_Output <- function(Output_Q_or_TQ, Candidates_n, SINGLE_TRAIT){
  # Only works for candidates all in the same Performance_year
  
  Performance_year_candidates <- Output_Q_or_TQ[Output_Q_or_TQ$ID == Candidates_n$ID[1], "Performance_year"]
  
  if(SINGLE_TRAIT == FALSE){
    # PROBLEM: duplicates EBVs column, becauqse already present in Output_Q before joining
    # To circumvent, we work on the subtable of Output with only candidates to update EBVs on
    Output_Q_or_TQ_n <- Output_Q_or_TQ %>% filter(Performance_year %in% Performance_year_candidates) %>%
      select(-c(EBV_d1_Q, EBV_m1_Q, EBV_d2_Q, EBV_m2_Q,
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
                E_VAR_E1, E_VAR_E2, E_COV_E1_E2))
    
    
    Output_Q_or_TQ_n <- Output_Q_or_TQ_n %>% left_join(Candidates_n %>% select(ID,
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
    Output_Q_or_TQ_n <- Output_Q_or_TQ %>% filter(Performance_year %in% Performance_year_candidates) %>%
      select(-c(EBV_d1_Q, EBV_m1_Q,
                EBV_d1_W, EBV_m1_W,
                SE_EBV_d1_Q, SE_EBV_m1_Q,
                SE_EBV_d1_W, SE_EBV_m1_W,
                E_Beekeeper_x_year_effect_T1,
                SE_Beekeeper_x_year_effect_T1,
                Conv_crit, Conv_round,
                E_VAR_A1_D, E_VAR_A1_M, E_COV_A1_D_M, 
                E_VAR_E1))
    
    Output_Q_or_TQ_n <- Output_Q_or_TQ_n %>% left_join(Candidates_n %>% select(ID,
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
  # Filtering out the Qs with new EBVs from output file
  Output_Q_or_TQ <- Output_Q_or_TQ %>% filter(!(Performance_year %in% Performance_year_candidates))
  # Row_binding the subtable with new EBVs
  Output_Q_or_TQ <- Output_Q_or_TQ %>% bind_rows(
    Output_Q_or_TQ_n)
  return(Output_Q_or_TQ)
}