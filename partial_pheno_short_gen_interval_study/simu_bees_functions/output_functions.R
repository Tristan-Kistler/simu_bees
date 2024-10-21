#######################################################################################
#    These are the functions which modify the output files of the simulation          #
#                                 from the simu_bees program                          #
#                                                                                     #
#######################################################################################

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Create_Output_Q                                                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Create_Output_Q <- function(WQ_n, NB_W1_Q, Rep){
  
  # gives warning 'Unknown or uninitialised column: `ID_Line`.': how to silence it?
  
  Output_Q <- WQ_n #initialization
  
  # F_Q <- pedigree %>% select(F_Q) %>% slice(1:NB_W1_Q)
  # F_Q <- rep(0, each = NB_W1_Q) # Replaced by F_ID already in WQ_n
  ID_BQ <- rep(0, each = NB_W1_Q)
  ID_D <- rep(0, each = NB_W1_Q)
  ID_PGM <- rep(0, each = NB_W1_Q)
  ID_PGM_mated <- rep(0, each = NB_W1_Q)
  Q_birth_year <- rep(1, each = NB_W1_Q) #the year indicated is the year of colony performance
  Performance_year <- rep(2, each = NB_W1_Q) #the year indicated is the year of colony performance
  Repetition <- rep(Rep, NB_W1_Q)
  
  Output_Q <- cbind(Output_Q, ID_BQ, ID_D, ID_PGM, ID_PGM_mated, Q_birth_year, Performance_year, Repetition)
  
  ID_Line <- rep(NA, each = NB_W1_Q)
  Output_Q <- cbind(Output_Q, ID_Line)
  
  Output_Q$ID_BQ <- as.character(Output_Q$ID_BQ)
  Output_Q$ID_D <- as.character(Output_Q$ID_D)
  Output_Q$ID_PGM  <- as.character(Output_Q$ID_PGM)
  Output_Q$ID_PGM_mated <- as.character(Output_Q$ID_PGM_mated)
  
  Output_Q$Q_birth_year <- as.integer(Output_Q$Q_birth_year)
  Output_Q$Performance_year <- as.integer(Output_Q$Performance_year)
  Output_Q$Repetition <- as.integer(Output_Q$Repetition)
  Output_Q$ID_type <- "Q"
  
  # Output_Q$ID_Line <- factor(Output_Q$ID_Line, levels = unique(WQ_n$ID_Line)) # Created warnings because of conflict in the factor format vs character... And not sure why I wanted it as factor
  
  return(Output_Q)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Create_Output_TQ                                                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Create_Output_TQ <- function(WTQ_n, Rep=1){
  NB_W1_TQ <- length(WTQ_n$ID)
  Output_TQ <- WTQ_n
  
  TQ_birth_year <- rep(as.integer(2), each = NB_W1_TQ)#the year indicated is the year of colony performance
  Performance_year <- rep(as.integer(3), each = NB_W1_TQ)#the year indicated is the year of colony performance
  Repetition <- rep(Rep, NB_W1_TQ)
  
  Output_TQ <- cbind(Output_TQ, TQ_birth_year, Performance_year, Repetition)
  Output_TQ$ID_type <- "TQ"
  
  # Because ID_Line is a factor in WQ_n, to not get a warning, we change it to a character vector
  Output_TQ$ID_Line <- as.character(Output_TQ$ID_Line)
  
  
  return(Output_TQ)
}
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#suppress_uninitialized_column_warning                                                #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# function "Create_Proven_BQs" creates a noisy warning.
# This function is to silence it
# suppress_uninitialized_column_warning <- function(expr) {
#   withCallingHandlers(
#     expr,
#     warning = function(w) {
#       if (grepl("uninitialised column", w$message)) {
#         # Suppress warning related to uninitialized columns
#         return(invisible())
#       }
#       # For other warnings, display normally
#       invokeRestart("muffleWarning", w)
#     }
#   )
# }

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Create_Proven_BQs                                                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Proven_BQs are BQs that survived their second overwintering and of which their TQ descendents have been phenotyped (hence the name Tester queens, testing by descent the BQs)
Create_Proven_BQs <- function(Q_or_TQ_n=BQ_n_save_one, NB_Q_or_TQ=NB_BQ, NB_W_Q_or_TQ=NB_W2_BQ, Blups_direct, Blups_maternal){
  
  # To suppress noisy warnings
  suppressWarnings({
    #Selecting wintering BQs
    Proven_BQs <- Winter_loss_Q(Q_or_TQ_n=BQ_n_save_one, NB_W_Q_or_TQ=NB_W2_BQ)[[1]]
    
    if(SEL_BLUP){
      #Updating blup values for BQs with Blups recalculated this year
      Proven_BQs <- Proven_BQs %>% select(ID, ID_cols, BV_d_Q, BV_m_Q, BV_d_moy_D, BV_m_moy_D, BV_d_W, BV_m_W, Breeder, Residual_effect, Perf, ID_Line)
      
      Proven_BQs <- Proven_BQs %>% left_join(Blups_direct, by=c("ID_cols"="ID"))
      Proven_BQs <- Proven_BQs %>% left_join(Blups_maternal, by=c("ID_cols"="ID"))
      
      Proven_BQs <- Proven_BQs %>% left_join(Blups_direct, by=c("ID"="ID"))
      Proven_BQs <- Proven_BQs %>% left_join(Blups_maternal, by=c("ID"="ID"))
      
      #We finally rename the added columns
      colnames(Proven_BQs)[(dim(Proven_BQs)[2]-7):dim(Proven_BQs)[2]] <- c("Blup_d_col", "CD_d_col", "Blup_m_col", "CD_m_col", "Blup_d_Q", "CD_d_Q", "Blup_m_Q", "CD_m_Q")
    }
    
    #adding infos about their descendent TQs
    for (i in 1:NB_W2_BQ){
      Proven_BQs$mean_desc_TQ_perf[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_Perf=mean(Perf))
      Proven_BQs$mean_desc_TQ_sd_perf[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(sd_Perf=sd(Perf))
      
      Proven_BQs$mean_desc_TQ_BV_d[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_BV_d=mean(BV_d_TQ))
      Proven_BQs$mean_desc_TQ_sd_BV_d[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(sd_BV_d=sd(BV_d_TQ))
      
      Proven_BQs$mean_desc_TQ_BV_m[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_BV_m=mean(BV_m_TQ))
      Proven_BQs$mean_desc_TQ_sd_BV_m[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(sd_BV_m=sd(BV_m_TQ))
      
      if(SEL_BLUP){#if Blups are calculated
        
        Proven_BQs$mean_desc_TQ_blup_d[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_d=mean(Blup_d_TQ))
        Proven_BQs$mean_desc_TQ_blup_d_sd[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_d_sd=sd(Blup_d_TQ))
        Proven_BQs$mean_desc_TQ_d_CD[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_CD_d_TQ=mean(CD_d_TQ))
        
        Proven_BQs$mean_desc_TQ_blup_m[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_m=mean(Blup_m_TQ))
        Proven_BQs$mean_desc_TQ_blup_m_sd[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_m_sd=sd(Blup_m_TQ))
        Proven_BQs$mean_desc_TQ_m_CD[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_CD_m_TQ=mean(CD_m_TQ))
      }
    }
    
    return(Proven_BQs)
  })
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Up_Output_proven_BQs                                                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#This function updates the output file with the new values of the proven BQs
Up_Output_proven_BQs <- function(Output_proven_BQs, Proven_BQs, Blups_direct, Blups_maternal, WTQ_n, NB_W2_BQ){
  suppressWarnings({
    if(SEL_BLUP){
      #Updating blup values for BQs with Blups recalculated this year
      Proven_BQs <- Proven_BQs %>% select(ID, ID_cols, BV_d_Q, BV_m_Q, BV_d_moy_D, BV_m_moy_D, BV_d_W, BV_m_W, Breeder, Residual_effect, Perf, ID_Line)
      
      Proven_BQs <- Proven_BQs %>% merge(x=Proven_BQs, y=Blups_direct, by.x="ID_cols", by.y="ID")
      Proven_BQs <- Proven_BQs %>% merge(x=Proven_BQs, y=Blups_maternal, by.x="ID_cols", by.y="ID")
      
      Proven_BQs <- Proven_BQs %>% merge(x=Proven_BQs, y=Blups_direct, by.x="ID", by.y="ID")
      Proven_BQs <- Proven_BQs %>% merge(x=Proven_BQs, y=Blups_maternal, by.x="ID", by.y="ID")
      
      #We finally renames the added columns correspondentely
      colnames(Proven_BQs)[(dim(Proven_BQs)[2]-7):dim(Proven_BQs)[2]] <- c("Blup_d_col", "CD_d_col", "Blup_m_col", "CD_m_col", "Blup_d_Q", "CD_d_Q", "Blup_m_Q", "CD_m_Q")
    }
    
    #adding infos about their descendent TQs
    for (i in 1:NB_W2_BQ){
      Proven_BQs$mean_desc_TQ_perf[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_Perf=mean(Perf))
      Proven_BQs$mean_desc_TQ_sd_perf[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(sd_Perf=sd(Perf))
      
      Proven_BQs$mean_desc_TQ_BV_d[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_BV_d=mean(BV_d_TQ))
      Proven_BQs$mean_desc_TQ_sd_BV_d[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(sd_BV_d=sd(BV_d_TQ))
      
      Proven_BQs$mean_desc_TQ_BV_m[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_BV_m=mean(BV_m_TQ))
      Proven_BQs$mean_desc_TQ_sd_BV_m[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(sd_BV_m=sd(BV_m_TQ))
      
      if(SEL_BLUP){#if Blups are calculated
        
        Proven_BQs$mean_desc_TQ_blup_d[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_d=mean(Blup_d_TQ))
        Proven_BQs$mean_desc_TQ_blup_d_sd[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_d_sd=sd(Blup_d_TQ))
        Proven_BQs$mean_desc_TQ_d_CD[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_CD_d_TQ=mean(CD_d_TQ))
        
        Proven_BQs$mean_desc_TQ_blup_m[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_m=mean(Blup_m_TQ))
        Proven_BQs$mean_desc_TQ_blup_m_sd[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_TQ_blup_m_sd=sd(Blup_m_TQ))
        Proven_BQs$mean_desc_TQ_m_CD[i] <- WTQ_n %>% filter(ID_BQ %in% Proven_BQs$ID[i]) %>% summarise(mean_CD_m_TQ=mean(CD_m_TQ))
      }
      
    }
    
    Output_proven_BQs <- bind_rows(Output_proven_BQs, Proven_BQs)
    
    new_Q <- (NB_W2_BQ*(year-3))#initilizing a counter
    
    Output_proven_BQs$Q_birth_year[(new_Q+1):(new_Q+NB_W2_BQ)] <- rep(year-2, NB_W2_BQ)
    Output_proven_BQs$Performance_year[(new_Q+1):(new_Q+NB_W2_BQ)] <- rep(year-1, NB_W2_BQ)
    Output_proven_BQs$Repetition[(new_Q+1):(new_Q+NB_W2_BQ)] <- rep(Rep, NB_W2_BQ)
    
    # for( i in (new_Q+1):(new_Q+NB_W2_BQ)){
    #   Output_proven_BQs$F_Q[i] <- pedigree[(NB_W1_Q+(NB_W1_Q+NB_BQ)*(year-4)+1):(NB_W1_Q+(NB_W1_Q+NB_BQ)*(year-3)), ] %>% filter(ID %in% Output_proven_BQs$ID[i]) %>% pull(F_Q)
    # }
    return(Output_proven_BQs)
  })
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Create_Output_Proven_BQs                                                              #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Create_Output_Proven_BQs <- function(Proven_BQs, NB_W2_BQ, Rep){
  Output_proven_BQs <- Proven_BQs
  
  Output_proven_BQs$ID_Line <- as.character(Output_proven_BQs$ID_Line)
  
  Repetition <- rep(Rep, NB_W2_BQ)
  Q_birth_year <- rep(1, NB_W2_BQ)
  Performance_year <- rep(2, NB_W2_BQ)
  F_Q <- rep(0, NB_W2_BQ) # useless?
  Output_proven_BQs <- cbind(Output_proven_BQs, Repetition, Q_birth_year, Performance_year, F_Q)
  
  return(Output_proven_BQs)
}
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Add_BQ_sel_dif_to_Output_Q                                                     #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# To be added directly in select BQ_n function
# Add_BQ_sel_dif_to_Output_Q <- function(BQ_n, Output_Q){
#   
#   # Which lines in Outpout_Q are those about Q in BQ_n (selected queens)?
#   which_line_in_Output_Q <- which(Output_Q$ID %in% BQ_n$ID) 
#   # Entering the selection diff for these queens in the Output_Q file
#   Output_Q$BQ_sel_diff[which_line_in_Output_Q] <- BQ_n$BQ_sel_diff[1]
#   
#   return(Output_Q)
# }

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Update_Output_Q_from_WQ_n                                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Update_Output_Q_from_WQ_n <- function(WQ_n, Output_Q, pedigree, year, Rep){
  
  NB_W1_Q <- nrow(WQ_n)
  Output_Q_n <- WQ_n %>% select(-BV_d1_BQ, -BV_m1_BQ, -phi_d1, -phi_m1,
                                -BV_d2_BQ, -BV_m2_BQ, -phi_d2, -phi_m2,
                                -F_BQ)
  
  # Sire_Mate_IDs <- pedigree %>% select(ID, ID_PGM_mated) %>%
  #   filter(ID %in% WQ_n$ID) %>%  #for queens born in this year-1, listed in the pedigree
  #   select(-ID)
    
  Q_birth_year <- rep(as.integer(year-1), each = NB_W1_Q)
  Performance_year <- rep(as.integer(year), each = NB_W1_Q)
  Repetition <- rep(as.integer(Rep), each = NB_W1_Q)
  
  Output_Q_n <- cbind(Output_Q_n, Q_birth_year, Performance_year, Repetition)
  Output_Q_n$ID_type <- "Q"
  #filling Output_Q accordingly. Couln't find a function already existing doing that, so I do it with joins and filters
  
  # Because ID_Line is a factor in WQ_n, to not get a warning, we change it to a character vector
  Output_Q_n$ID_Line <- as.character(Output_Q_n$ID_Line)
  
  # Homogenizing names: ID_PGM in WQ_n is actually the ID_PGM_mated
  # colnames(Output_Q_n)[colnames(Output_Q_n) == "ID_PGM"] <- "ID_PGM_mated" # should be avoided by giving it the right name during mating in WQ_n
  
  Output_Q <- bind_rows(Output_Q, Output_Q_n)
  
  return(Output_Q)
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Update_Output_TQ                                                      #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
Update_Output_TQ <- function(WTQ_n, Output_TQ, year){
  
  NB_W1_TQ <- nrow(WTQ_n)
  Output_TQ_n <- WTQ_n
  
  TQ_birth_year <- rep(as.integer(year-1), each = NB_W1_TQ)
  Performance_year <- rep(as.integer(year), each = NB_W1_TQ)
  Repetition <- rep(as.integer(Rep), each = NB_W1_TQ)
  
  Output_TQ_n <- cbind(TQ_birth_year, Performance_year, Repetition, Output_TQ_n)
  Output_TQ_n$ID_type <- "TQ"
  
  # Because ID_Line is a factor in WQ_n, to not get a warning, we change it to a character vector
  Output_TQ_n$ID_Line <- as.character(Output_TQ_n$ID_Line)
  
  Output_TQ <- bind_rows(Output_TQ, Output_TQ_n)
  
  return(Output_TQ)
}
