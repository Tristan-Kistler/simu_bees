#######################################################################################
#    These are the functions to breed of the simu_bees program                        #
#                                                                                     #
#######################################################################################

# Below two functions also need implementation to allow Qs and TQs to be on the same apiaries
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Add_colony_Performance                                                    #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Add_colony_Performance <- function(WQ_or_WTQ_n=WQ_n, indiv_type="Q", Beekeepers_fixed_effects_n, VAR_E){
  
  Nb_indiv <- dim(WQ_or_WTQ_n)[1]
  WQ_or_WTQ_n <- inner_join(WQ_or_WTQ_n, Beekeepers_fixed_effects_n, by = "Beekeeper_year")
  
  Residual_effects_T1_T2 <- rmvnorm(n=Nb_indiv, sigma=VAR_E)
  WQ_or_WTQ_n$Residual_effect_T1 <- Residual_effects_T1_T2[, 1] # res effects trait 1
  WQ_or_WTQ_n$Residual_effect_T2 <- Residual_effects_T1_T2[, 2] # res effects trait 2
  
  # Depending on if we are adding phenotypes to TQ or Q colonies, the columns of fixed effects and BVs have different names
  if (indiv_type == "Q"){
    WQ_or_WTQ_n <- WQ_or_WTQ_n %>% mutate(Perf_T1 = BV_m1_Q + BV_d1_W + Beekeeper_x_year_effect_T1 + Residual_effect_T1,
                                          Perf_T2 = BV_m2_Q + BV_d2_W + Beekeeper_x_year_effect_T2 + Residual_effect_T2)
  }else if (indiv_type == "TQ"){
    WQ_or_WTQ_n <- WQ_or_WTQ_n %>% mutate(Perf_T1 = BV_m1_TQ + BV_d1_W + Beekeeper_x_year_effect_T1 + Residual_effect_T1,
                                          Perf_T2 = BV_m2_TQ + BV_d2_W + Beekeeper_x_year_effect_T2 + Residual_effect_T2)
  }
  return(WQ_or_WTQ_n)
}


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Assign_indiv_keepers_max_connexion                                        #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Assign_indiv_keepers_max_connexion <- function(Q_or_TQ_n, Nb_keepers, Indiv_per_line, Q_AND_TQ_SAME_APIARIES, indiv_type, Performance_year=year){
  # NB_BREEDERS </= NB_BQ & NB_Q_PER_BQ must be a multiple of NB_BREEDER
  
  Indiv_per_line_per_keeper <- Indiv_per_line/Nb_keepers #we calculate how many queen of a same line will be distributed to one breeder
  
  Q_Distribution <- c()#initialising the vector
  line_distribution <- rep((1:Nb_keepers), Indiv_per_line_per_keeper) #we get a first Q distribution for the first line, repetitive in order
  
  for (i in 1:NB_BQ){ #we loop on the nb of lines to distribute the queens to
    Q_Distribution <- c(Q_Distribution, sample(line_distribution)) #we enlarge at each loop cycle the Q distribution vector with randomised line_distribution vectors
  }
  
  Q_or_TQ_n$Beekeeper_year <- Q_Distribution #we add this Q or TQ distribution in WQ_n or WTQ_n
  if(Q_AND_TQ_SAME_APIARIES == TRUE){
    Q_or_TQ_n$Beekeeper_year <- paste0("A_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "Q"){
    Q_or_TQ_n$Beekeeper_year <- paste0("B_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "TQ"){
    Q_or_TQ_n$Beekeeper_year <- paste0("T_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }
  
  return(Q_or_TQ_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Assign_indiv_keepers_min_connexion                                        #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Assign_indiv_keepers_min_connexion <- function(Q_or_TQ_n, Nb_keepers, Indiv_per_line, Q_AND_TQ_SAME_APIARIES, indiv_type, Performance_year=year){
  # NB_BREEDER == NB_BQ & NB_BQ > 1 & NB_Q_PER_BQ must be a multiple of 2
  
  Sampled_keepers <- sample(c(1:Nb_keepers)) #we get a vector of keepers in randomised order
  
  Q_first_keeper <- rep(Sampled_keepers[1], each = Indiv_per_line/2) #we create a vector which will be used to distribute the first and last line to this first randomly chosen keeper
  
  Q_Distribution <- c(Q_first_keeper, rep(Sampled_keepers[2:(Nb_keepers)], each = Indiv_per_line), Q_first_keeper) #we create the Q or TQ distribution vector
  
  Q_or_TQ_n$Beekeeper_year <- Q_Distribution #we add this Q or TQ distribution in WQ_n or WTQ_n
  
  # Depending on if Qs and TQs are on the same apiaries, we give names to the level of the beekeeper effect with different coding or not
  if(Q_AND_TQ_SAME_APIARIES == TRUE){
    Q_or_TQ_n$Beekeeper_year <- paste0("A_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "Q"){
    Q_or_TQ_n$Beekeeper_year <- paste0("B_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "TQ"){
    Q_or_TQ_n$Beekeeper_year <- paste0("T_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }
  
  return(Q_or_TQ_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Assign_indiv_keepers_any_connection                                       #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Assign_indiv_keepers_any_connection <- function(Q_or_TQ_n, indiv_type, Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_KEEPERS, Nb_sister_groups_per_keeper=NB_SISTER_GROUPS_PER_KEEPER, Nb_keepers_per_sister_group=NB_KEEPERS_PER_SISTER_GROUP, Nb_sisters_per_group_per_keeper=NB_SISTERS_PER_GROUP_PER_KEEPER, NB_BQ=NB_BQ, Performance_year=year){
  # distribute groups of sisters over keepers in order to have a certain Nb_sister_groups_per_keeper
  
  
  # Literal algorithm:
  # 1) Prepare a list of possible spots on each apiary to distribute t sisters to (so N lists of k spots)
  # 2) For each sister group i:
  # among apiaries with the most spots, draw K apiaries (if less than K apiaries with most spots available, chose them all, then chose other apaiaries among the other possibilities)
  # 3) Remove the chosen apiaries from the list of possible spots
  # 4) Repeat until all families have been distributed: it should have left no spots available.
  # This algorithm tries to prevent symmetrical attribution of sister groups, where there could be groups of apiaries never connected with each other by some families. It can still happen by random chance though. 
  
  # Initializing
  Q_or_TQ_n$Beekeeper_year <- NA
  List_possible_sis_to_apiary <- lapply(1:Nb_keepers, function(x) 1:Nb_sister_groups_per_keeper)
  itel <- 0
  
  for (i in 1:NB_BQ) { # for each sister group i
    # print(paste0("this is loop i=", i))
    # we sample from the vector of Nb_keepers vectors of Nb_sister_groups_per_keeper elements to decide to which apiary sister group of family 'i' is sent
    non_zero_positions <- which(sapply(List_possible_sis_to_apiary, function(x) length(x) > 0)) # apiaries that still have to receive a sister groups (non-zero subarrays)
    # print(paste0("non_zero_positions: "))
    # print(paste0(non_zero_positions))
    Nb_possibilities_per_apiary <- sapply(List_possible_sis_to_apiary, function(x) length(x)) # apiaries that still have to receive a sister groups (non-zero subarrays)
    # we sample max K apiaries among the ones needing to get the most families
    # We first count how many apiaries have max spots left
    keepers_needing_fam_most <- which(Nb_possibilities_per_apiary == max(Nb_possibilities_per_apiary))
    if (length(keepers_needing_fam_most) <= Nb_keepers_per_sister_group){
      chosen_keepers <- keepers_needing_fam_most # To not leave apiaries behind in family distribution, we draw the ones that need most to be drawn
    }else{
      chosen_keepers <- sample(keepers_needing_fam_most, size = Nb_keepers_per_sister_group, replace = FALSE) # If there are more than one possibilities, we randomly choose among them
    }
    # print(paste0("chosen_keepers among most needing:"))
    # print(paste0(chosen_keepers))
    if (length(chosen_keepers) < Nb_keepers_per_sister_group){ # If we still need to draw keepers
      # We eliminate the first random draw fromt he possible apiaries to draw:
      non_zero_positions <- non_zero_positions[-which(non_zero_positions %in% chosen_keepers)]
      # print(paste0("non_zero_positions without first draw"))
      # print(paste0(non_zero_positions))
      Add_keepers_to_draw <- Nb_keepers_per_sister_group-length(chosen_keepers)
      if(Add_keepers_to_draw == length(non_zero_positions)){
        chosen_keepers_add <- non_zero_positions # If keepers from which we can choose are all the ones we need to choose, we choose them
      }else{
        chosen_keepers_add <- sample(non_zero_positions, size = Add_keepers_to_draw, replace = FALSE) # If there are more than needed possibilities, we randomly choose among them
      }
      # adding chosen_keepers to chosen_keepers
      # print(paste0("chosen_keepers before adding chosen_keepers"))
      # print(paste0(chosen_keepers_add)) 
      chosen_keepers <- c(chosen_keepers, chosen_keepers_add)
    }
    # print(paste0("chosen_keepers"))
    # print(paste0(chosen_keepers))
    # we eliminate one element of the Nb_sister_groups_per_keeper-sized arrays that have been chosen
    for (j in chosen_keepers) {
      List_possible_sis_to_apiary[[j]] <- List_possible_sis_to_apiary[[j]][-1]
      # print(paste0("List_possible_sis_to_apiary"))
      # print(List_possible_sis_to_apiary)
    }
    for (keeper_x in 1:Nb_keepers_per_sister_group){
      for (col in 1:Nb_sisters_per_group_per_keeper) { # for each sister of sister-group i on apiary j
        itel <- itel + 1 # me move from one line down in the Q_or_TQ_n matrix table
        # Q_or_TQ_n[itel, 1] <- i # Family number
        Q_or_TQ_n$Beekeeper_year[itel] <- chosen_keepers[keeper_x] # we store an x in the column of the apiary where the sister of a sister group will be tested
      }
    }
  }
  
  # Depending on if Qs and TQs are on the same apiaries, we give names to the level of the beekeeper effect with different coding or not
  if(Q_AND_TQ_SAME_APIARIES == TRUE){
    Q_or_TQ_n$Beekeeper_year <- paste0("A_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "Q"){
    Q_or_TQ_n$Beekeeper_year <- paste0("B_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "TQ"){
    Q_or_TQ_n$Beekeeper_year <- paste0("T_", Q_or_TQ_n$Beekeeper_year, "_", Performance_year)
  }
  
  return(Q_or_TQ_n)
}

# NON RANDOM VERSION: CAN CREATE UNCONNECTED GROUPS OF APIARIES
# Assign_indiv_keepers_any_connection <- function(Q_or_TQ_n, indiv_type, Q_AND_TQ_SAME_APIARIES, Nb_keepers=NB_KEEPERS, Indiv_per_line=NB_Q_PER_BQ, Nb_sister_groups_per_keeper=NB_SISTER_GROUPS_PER_KEEPER, Nb_keepers_per_sister_group=NB_KEEPERS_PER_SISTER_GROUP, Apiary_size=APIARY_SIZE, Nb_sisters_per_group_per_keeper=NB_SISTERS_PER_GROUP_PER_KEEPER, NB_BQ){
#   # distribute groups of sisters over keepers in order to have a certain Nb_sister_groups_per_keeper
#   
#   # Initializing the table that will contain the Family num Id and the keepers num IDs
#   Families_to_Beekeepers <- matrix(nrow = (NB_BQ * Indiv_per_line), ncol = 2)
#   # initialization
#   keeper <- 0
#   Colony_x <- 0
#   for (i in sample(1:NB_BQ, NB_BQ, replace = F)) { # for each sister group i, taken from a order-randomized vector of families
#     for (j in 1:Nb_keepers_per_sister_group) { # for each apiary j testing sister group i
#       keeper <- keeper + 1
#       if (keeper > Nb_keepers) {
#         keeper <- 1
#       } # if we still have apiaries to distribute colonies to, if not, we start again at first apiary
#       for (k in 1:Nb_sisters_per_group_per_keeper) { # for each sister of sister-group i on apiary j
#         Colony_x <- Colony_x + 1 # me move from one line down in the Q_or_TQ_n matrix table
#         Families_to_Beekeepers[Colony_x, 1] <- i # we store the sister group number in the first column of the matrix
#         Families_to_Beekeepers[Colony_x, 2] <- keeper # we store the keeper number in the second column of the matrix
#       }
#     }
#   }
#   
#   # We order according to family nb ID
#   Families_to_Beekeepers <- Families_to_Beekeepers[order(Families_to_Beekeepers[, 1]), ]
#   
#   Q_or_TQ_n$Beekeeper_year <- Families_to_Beekeepers[, 2] #we add this Q or TQ to apiaries distribution in Q_n or TQ_n
#   
#   # Depending on if Qs and TQs are on the same apiaries, we give names to the level of the beekeeper effect with different coding or not
#   if(Q_AND_TQ_SAME_APIARIES == TRUE){
#     Q_or_TQ_n$Beekeeper_year <- paste0("A_", Q_or_TQ_n$Beekeeper_year)
#   }else if (indiv_type == "Q"){
#     Q_or_TQ_n$Beekeeper_year <- paste0("B_", Q_or_TQ_n$Beekeeper_year)
#   }else if (indiv_type == "TQ"){
#     Q_or_TQ_n$Beekeeper_year <- paste0("T_", Q_or_TQ_n$Beekeeper_year)
#   }
#   
#   return(Q_or_TQ_n)
# }

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Produce_Wgroup                                                               #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Produce_Wgroup <- function (W_indiv_n=WQ_n, Ds_n=D_n, indiv_type="Q"){
  #We first summarize breeding values of drones having mated a queen
  BV_D_moy_n <- Ds_n %>% 
    group_by(ID_mate) %>% 
    summarise(BV_d1_moy_D = mean(BV_d1_D), 
              BV_m1_moy_D = mean(BV_m1_D),
              BV_d2_moy_D = mean(BV_d2_D), 
              BV_m2_moy_D = mean(BV_m2_D)) %>%
    ungroup()
  
  if(indiv_type=="Q"){
    #We use this drones' infos and that of the queen to produce workers' BVs
    W_indiv_n <- left_join(W_indiv_n, BV_D_moy_n, by= c("ID"="ID_mate"), suffix=c(".x","")) %>% 
      mutate(BV_d1_W = 0.5*BV_d1_Q+BV_d1_moy_D, 
             BV_m1_W = 0.5*BV_m1_Q+BV_m1_moy_D,
             BV_d2_W = 0.5*BV_d2_Q+BV_d2_moy_D, 
             BV_m2_W = 0.5*BV_m2_Q+BV_m2_moy_D) #adding the data contained in BV_D_moy_n by matching of the ID column
  }else if(indiv_type=="TQ"){
    #We use this drones' infos and that of the queen to produce workers' BVs
    W_indiv_n <- left_join(W_indiv_n, BV_D_moy_n, by= c("ID"="ID_mate"), suffix=c(".x","")) %>% 
      mutate(BV_d1_W = 0.5*BV_d1_TQ+BV_d1_moy_D, 
             BV_m1_W = 0.5*BV_m1_TQ+BV_m1_moy_D,
             BV_d2_W = 0.5*BV_d2_TQ+BV_d2_moy_D, 
             BV_m2_W = 0.5*BV_m2_TQ+BV_m2_moy_D) #adding the data contained in BV_D_moy_n by matching of the ID column
    
  }
  
  return (W_indiv_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Produce_indiv_randomly                                                    #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Produce_indiv_randomly <- function (Nb_indiv, indiv_type, Q_to_mate, Var_A=VAR_A, year=year){
  # Q_to_mate is only used for males (D and poential DPQs)
  if(indiv_type=="Q"){
    indiv_n <- matrix(ncol = 6, nrow =Nb_indiv)
    colnames(indiv_n) <- c("ID", "BV_d1_Q", "BV_m1_Q", "BV_d2_Q", "BV_m2_Q", "F_ID")
    indiv_n <- as_tibble(indiv_n)
    indiv_n$F_ID <- 0
  }else if(indiv_type=="D"){
    indiv_n <- matrix(ncol = 6, nrow =Nb_indiv)
    colnames(indiv_n) <- c("ID", "BV_d1_D", "BV_m1_D", "BV_d2_D", "BV_m2_D", "ID_mate")
    indiv_n <- as_tibble(indiv_n)
    indiv_n$ID_mate <- unlist(map2(Q_to_mate$ID, NB_D_FECUNDATING, rep))
  }else if(indiv_type=="TD"){
    indiv_n <- matrix(ncol = 7, nrow =Nb_indiv)
    colnames(indiv_n) <- c("ID", "BV_d1_D", "BV_m1_D", "BV_d2_D", "BV_m2_D", "ID_mate", "F_ID")
    indiv_n <- as_tibble(indiv_n)
    indiv_n$ID_mate <- unlist(map2(Q_to_mate$ID, NB_D_FECUNDATING, rep))
    indiv_n$F_ID <- 0
  }
  
  # Naming the individuals
  indiv_n[, 1] <- paste0(indiv_type, 1:Nb_indiv, "_", year)
  
  # Drawing breeding values from a mulivariate normal distribution depending on the Var_A parameter
  BVs <- rmvnorm(n=Nb_indiv, sigma=Var_A)
  indiv_n[, 2] <- BVs[, 2] # direct (worker) effects trait 1
  indiv_n[, 3] <- BVs[, 1] # maternal (worker) effects trait 1
  indiv_n[, 4] <- BVs[, 4] # direct (worker) effects trait 2
  indiv_n[, 5] <- BVs[, 3] # maternal (worker) effects trait 2
  
  return(indiv_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#          Get_beekeepers_fixed_effects                                               #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# version for Qs and TQs on same apiaries
Get_beekeepers_fixed_effects <- function(Nb_keepers=NB_KEEPERS, indiv_type="Q", Performance_year=year, SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2, Q_AND_TQ_SAME_APIARIES){
  
  Beekeepers_fixed_effects_n <- data.frame(matrix(ncol = 3, nrow = Nb_keepers))
  colnames(Beekeepers_fixed_effects_n) <- c("Beekeeper_year",
                                            "Beekeeper_x_year_effect_T1",
                                            "Beekeeper_x_year_effect_T2")
  
  Beekeepers_fixed_effects_n$Beekeeper_year <- seq(from=1, to=Nb_keepers, by=1)
  # Depending on if Qs and TQs are on the same apiaries, we give names to the level of the beekeeper effect with different coding or not
  if(Q_AND_TQ_SAME_APIARIES == TRUE){
    Beekeepers_fixed_effects_n$Beekeeper_year <- paste0("A_", Beekeepers_fixed_effects_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "Q"){
    Beekeepers_fixed_effects_n$Beekeeper_year <- paste0("B_", Beekeepers_fixed_effects_n$Beekeeper_year, "_", Performance_year)
  }else if (indiv_type == "TQ"){
    Beekeepers_fixed_effects_n$Beekeeper_year <- paste0("T_", Beekeepers_fixed_effects_n$Beekeeper_year, "_", Performance_year)
  }
  
  Beekeepers_fixed_effects_n$Beekeeper_x_year_effect_T1 <- rnorm(Nb_keepers, mean=0, sd = SIGMA_BEEKEEPER_x_YEAR_T1)
  Beekeepers_fixed_effects_n$Beekeeper_x_year_effect_T2 <- rnorm(Nb_keepers, mean=0, sd = SIGMA_BEEKEEPER_x_YEAR_T2)
  
  return(Beekeepers_fixed_effects_n)
}
# 
# #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# #          Get_breeders_fixed_effects                                                 #
# #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# Get_breeders_fixed_effects <- function(NB_BREEDERS, SIGMA_YEAR_T1,	SIGMA_YEAR_T2, 
#                                        SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2,
#                                        Breeders_effect_T1, Breeders_effect_T2){
#   
#   Breeders_fixed_effects_n <- data.frame(matrix(ncol = 7, nrow = NB_BREEDERS))
#   colnames(Breeders_fixed_effects_n) <- c("Breeder",
#                                           "Year_effect_T1", "Breeder_effect_T1", "Breeder_x_year_effect_T1",
#                                           "Year_effect_T2", "Breeder_effect_T2", "Breeder_x_year_effect_T2")
#   
#   Breeders_fixed_effects_n$Year_effect_T1 <- rnorm(1, mean=0, sd = SIGMA_YEAR_T1)
#   Breeders_fixed_effects_n$Year_effect_T2 <- rnorm(1, mean=0, sd = SIGMA_YEAR_T2)
#   
#   Breeders_fixed_effects_n$Breeder <- seq(from=1, to=NB_BREEDERS, by=1)
#   Breeders_fixed_effects_n$Breeder_effect_T1 <- Breeders_effect_T1
#   Breeders_fixed_effects_n$Breeder_effect_T2 <- Breeders_effect_T2
#   Breeders_fixed_effects_n$Breeder_x_year_effect_T1 <- rnorm(NB_BREEDERS, mean=0, sd = SIGMA_BEEKEEPER_x_YEAR_T1)
#   Breeders_fixed_effects_n$Breeder_x_year_effect_T2 <- rnorm(NB_BREEDERS, mean=0, sd = SIGMA_BEEKEEPER_x_YEAR_T2)
#   
#   return(Breeders_fixed_effects_n)
# }
# 
# #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# #          Get_Testers_fixed_effects                                                  #
# #°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# Get_Testers_fixed_effects <- function(NB_TESTERS, SIGMA_YEAR_T1,	SIGMA_YEAR_T2,
#                                       SIGMA_BEEKEEPER_x_YEAR_T1, SIGMA_BEEKEEPER_x_YEAR_T2,
#                                       Breeders_effect_T1, Breeders_effect_T2){
#   
#   Testers_fixed_effects_n <- data.frame(matrix(ncol = 7, nrow = NB_TESTERS))
#   colnames(Testers_fixed_effects_n) <- c("Tester",
#                                          "Year_effect_T1", "Tester_effect_T1", "Tester_x_year_effect_T1",
#                                          "Year_effect_T2", "Tester_effect_T2", "Tester_x_year_effect_T2")
#   
#   Testers_fixed_effects_n$Year_effect_T1 <- rnorm(1, mean=0, sd = SIGMA_YEAR_T1)
#   Testers_fixed_effects_n$Year_effect_T2 <- rnorm(1, mean=0, sd = SIGMA_YEAR_T2)
#   
#   Testers_fixed_effects_n$Tester <- seq(from=1, to=NB_TESTERS, by=1)
#   Testers_fixed_effects_n$Tester_effect_T1 <- Testers_effect_T1 #we get the Tester's fixed effect in the previously created vector
#   Testers_fixed_effects_n$Tester_effect_T2 <- Testers_effect_T2 #we get the Tester's fixed effect in the previously created vector
#   Testers_fixed_effects_n$Tester_x_year_effect_T1 <- rnorm(NB_TESTERS, mean=0, sd = SIGMA_BEEKEEPER_x_YEAR_T1)
#   Testers_fixed_effects_n$Tester_x_year_effect_T2 <- rnorm(NB_TESTERS, mean=0, sd = SIGMA_BEEKEEPER_x_YEAR_T2)
#   
#   return(Testers_fixed_effects_n)
# }

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Winter_loss_Q                                                             #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Winter_loss_Q <- function (Q_or_TQ_n, NB_W_Q_or_TQ){
  
  NB_Q_or_TQ <- dim(Q_or_TQ_n)[1]
  
  Q_W_n <- sort.int(sample(c(1:(NB_Q_or_TQ)), NB_W_Q_or_TQ))
  
  if (Q_WINTER1_LOSS_RATE == 0){
    #if there is no loss of TQs during winter, we keep all TQs
  }else{
    Q_or_TQ_n <- Q_or_TQ_n[c(Q_W_n), ]
  }
  Q_or_TQ_and_Q_W_n <- list(Q_or_TQ_n, Q_W_n)
  return(Q_or_TQ_and_Q_W_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Create_indiv_n_blanck                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Create_indiv_n_blanck <- function(BQ_n, indiv_type, Nb_indiv, Nb_indiv_per_BQ, year){
  #Initializing the table
  
  indiv_n <- matrix(ncol = 8, nrow = Nb_indiv)
  colnames(indiv_n) <- c("ID", "ID_BQ","ID_D", 
                         paste0("BV_d1_", indiv_type), paste0("BV_m1_", indiv_type), 
                         paste0("BV_d2_", indiv_type), paste0("BV_m2_", indiv_type),
                         "ID_Line")
  
  indiv_n <- as_tibble(indiv_n)
  
  #Adding Q IDs
  indiv_n[, 1] <- paste0(indiv_type, 1:Nb_indiv, "_", year)
  
  #Adding their mother's IDs
  indiv_n$ID_BQ <- unlist(map2(BQ_n$ID, Nb_indiv_per_BQ, rep))
  indiv_n$ID_Line <- unlist(map2(BQ_n$ID_Line, Nb_indiv_per_BQ, rep))
  indiv_n$ID_Line <- factor(indiv_n$ID_Line, levels = unique(BQ_n$ID_Line)) # needed for line extinction?
  
  return(indiv_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Fill_indiv_n                                                              #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Fill_indiv_n <- function (indiv_n, indiv_type, D_n, BQ_n, VAR_A_meiose, Nb_indiv_per_BQ, NB_BQ, year){
  
  Nb_indiv <- dim(indiv_n)[1]
  
  # we create a table with only Ds having mated the BQs
  D_n_BQ <- D_n %>% filter(ID_mate %in% BQ_n$ID)#filtering preserves initial appearing order
  
  # we order lines as they are originally
  D_n_BQ$Original_order <- seq(1:(NB_D_FECUNDATING*NB_BQ))
  
  #for each BQ (for each mother), we draw Ds that mated this BQ (with replacement) to ramdomly choose father D for new Qs.
  #infos_D_n_per_Q will be disordered by the grouping, but we get the initial order back thanks to column Original_order. 
  #Nb of Ds having mated BQs could be variable.
  infos_D_n_per_Q <- D_n_BQ %>%
    group_by(ID_mate) %>% 
    sample_n(size=Nb_indiv_per_BQ, replace=T) %>%
    ungroup()
  
  # as grouping disorders lines, we reorder them in the original sequence, missing rows excluded
  infos_D_n_per_Q <- setorder(infos_D_n_per_Q, Original_order)
  
  if(identical(indiv_n$ID_BQ, infos_D_n_per_Q$ID_mate)){#We check if in fact order of BQs in indiv_n is the same as in infos_D_n_per_Q, so that we give the right infos about Ds having met Qs' BQs.
    #It should always be the case, but as this is a crucial step, a fast checking is additionally carried out.
    indiv_n$ID_D <-  infos_D_n_per_Q$ID
    indiv_n$BV_d1_D <- infos_D_n_per_Q$BV_d1_D
    indiv_n$BV_m1_D <- infos_D_n_per_Q$BV_m1_D
    indiv_n$BV_d2_D <- infos_D_n_per_Q$BV_d2_D
    indiv_n$BV_m2_D <- infos_D_n_per_Q$BV_m2_D
  }else{
    print("Appearing order of BQs in table indiv_n is different from that of infos_D_n_per_Q")
  }
  
  #Adding infos about the mothers (BQs)
  indiv_n <- indiv_n %>% 
    left_join(BQ_n %>% 
                select(ID, BV_d1_Q, BV_m1_Q, BV_d2_Q, BV_m2_Q, F_ID),  
              by = c("ID_BQ" = "ID"),  suffix=c("","_BQ") ) # %>%
  # left_join(pedigree %>% 
  #             # slice((max(NB_W1_Q+(NB_W1_Q+NB_BQ)*(year-3), 0) +1):(max(NB_W1_Q+(NB_W1_Q+NB_BQ)*(year-3), 0)+NB_W1_Q)) %>% 
  #             select(ID, F_Q), 
  #           by = c("ID_BQ" = "ID"), suffix=c("","_BQ"))
  #the slice function is used to minimize the portion of pedigree that has to be left_jointed to indiv_n, so that
  #this sliced pedigree is of constant size during all simulation years, which is computationally less demanding.
  #For this slice to give back a portion of the pedigree of 1:NB_W1_Q in year 2, a workaround was used with function max.
  #In fact, when year==2, NB_W1_Q+(NB_W1_Q+NB_BQ)*(year-3) < 0, whereas is should be 0. 
  #This is computationally faster than using a if (year==2){...} else {...}
  
  phi <- rmvnorm(n=Nb_indiv, sigma=VAR_A_meiose)#Variance of mendelian sampling terms is 1/4*VAR_A
  phi_m1 <- phi[, 1]#mendelian sampling terms of maternal effects for trait 1
  phi_d1 <- phi[, 2]#mendelian sampling terms of direct effects for trait 1
  phi_m2 <- phi[, 3]#mendelian sampling terms of maternal effects for trait 2
  phi_d2 <- phi[, 4]#mendelian sampling terms of direct effects for trait 2
  
  indiv_n <- cbind(indiv_n, phi_d1, phi_m1, phi_d2, phi_m2)
  
  if(indiv_type=="Q"){
    indiv_n <-  indiv_n %>% mutate(BV_d1_Q = 0.5*BV_d1_Q_BQ + BV_d1_D + sqrt(1 - F_ID)* phi_d1,
                                   BV_m1_Q = 0.5*BV_m1_Q_BQ + BV_m1_D + sqrt(1 - F_ID)* phi_m1,
                                   BV_d2_Q = 0.5*BV_d2_Q_BQ + BV_d2_D + sqrt(1 - F_ID)* phi_d2,
                                   BV_m2_Q = 0.5*BV_m2_Q_BQ + BV_m2_D + sqrt(1 - F_ID)* phi_m2) %>%
      rename(BV_d1_BQ=BV_d1_Q_BQ, BV_m1_BQ=BV_m1_Q_BQ, BV_d2_BQ=BV_d2_Q_BQ, BV_m2_BQ=BV_m2_Q_BQ, F_BQ=F_ID)
    
  }else if(indiv_type=="TQ"){
    indiv_n <-  indiv_n %>% mutate(BV_d1_TQ = 0.5*BV_d1_Q + BV_d1_D + sqrt(1 - F_ID)* phi_d1,
                                   BV_m1_TQ = 0.5*BV_m1_Q + BV_m1_D + sqrt(1 - F_ID)* phi_m1,
                                   BV_d2_TQ = 0.5*BV_d2_Q + BV_d2_D + sqrt(1 - F_ID)* phi_d2,
                                   BV_m2_TQ = 0.5*BV_m2_Q + BV_m2_D + sqrt(1 - F_ID)* phi_m2,) %>%
      rename(BV_d1_BQ=BV_d1_Q, BV_m1_BQ=BV_m1_Q, BV_d2_BQ=BV_d2_Q, BV_m2_BQ=BV_m2_Q, F_BQ=F_ID)
  }
  
  return(indiv_n)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Create_Q_n_blanck                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Create_Q_n_blanck                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Create_Q_n_blanck                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Create_Q_n_blanck                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Create_Q_n_blanck                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#           Create_Q_n_blanck                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°