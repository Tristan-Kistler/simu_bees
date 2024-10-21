#######################################################################################
#    These are the functions which declaration                                        #
#    depends on the simulation scenario defined by the parameters                     #
#                           of the simu_bees program                                  #
#######################################################################################

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Select_BQ                                                                            #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# First selection: Select_BQ_foundress
# RANDOM VERSION
if(Q_SELECTION_CRITERION=="Random"){
  Select_BQ_foundress <- function(WQ_n, NB_BQ, Q_SELECTION_CRITERION="Random"){
    WQ_n$Original_order <- seq(1:(nrow(WQ_n)))
    BQ_n <- WQ_n %>% sample_n(NB_BQ)
    #ordering BQs as originally there were ordered in WQ_n
    BQ_n <- setorder(BQ_n, Original_order)
    WQ_n <- WQ_n %>% select(-Original_order)
    BQ_n <- BQ_n %>% select(-Original_order)
    
    BQ_n$ID_Line <- paste0("L", 1:NB_BQ)
    
    return(BQ_n)
  } # end of function declaration
  
}else{ #If Q selection is not random
  # NON-RANDOM VERSION
  Select_BQ_foundress <- function(WQ_n, NB_BQ, Q_SELECTION_CRITERION="Perf_T1/sd(Perf_T1) + Perf_T2/sd(Perf_T1)"){
    BQ_n <- WQ_n %>% top_n(NB_BQ, eval(parse(text=Q_SELECTION_CRITERION)) )
    BQ_n <- BQ_n[1:NB_BQ, ]#top_n issues the best Q's but also potentially the execos. We therefore readjust the size in this eventuality
    
    BQ_n$ID_Line <- paste0("L", 1:NB_BQ)
    
    return(BQ_n)
  }
}

# Further selections: Select_BQ
# No L selection #####################################################################
# RANDOM VERSION
  if(MAT_LINE_SELECTION == FALSE){
    if(Q_SELECTION_CRITERION=="Random"){
      Select_BQ <- function(WQ_n, NB_BQ, Q_SELECTION_CRITERION="Random"){
        WQ_n$Original_order <- seq(1:(nrow(WQ_n)))
        BQ_n <- WQ_n %>% sample_n(NB_BQ)
        #ordering BQs as originally there were ordered in WQ_n
        BQ_n <- setorder(BQ_n, Original_order)
        WQ_n <- WQ_n %>% select(-Original_order)
        BQ_n <- BQ_n %>% select(-Original_order)
        
        return(BQ_n)
      }
    }else{ #If Q selection is not random
      Select_BQ <- function(WQ_n, NB_BQ, Q_SELECTION_CRITERION="Perf_T1 + Perf_T2"){
# NON- RANDOM VERSION
        BQ_n <- WQ_n %>% top_n(NB_BQ, eval(parse(text=Q_SELECTION_CRITERION)) )
        BQ_n <- BQ_n[1:NB_BQ, ]#top_n issues the best Q's but also potentially the execos. We therefore readjust the size in this eventuality
        
        return(BQ_n)
      } # end of function declaration
    } # end of NON random selection
  } # end of NO L selection

#With L selection ##################################################################
if(MAT_LINE_SELECTION == TRUE){
  # RANDOM VERSION
  if(Q_SELECTION_CRITERION=="Random"){
    Select_BQ <- function(WQ_n, NB_BQ, Q_SELECTION_CRITERION="Random"){
      Survived_lines <- levels(factor(WQ_n$ID_Line, levels = unique(WQ_n$ID_Line)))
      # storing the original order of queens to reorder them after sampling
      WQ_n$Original_order <- seq(1:(nrow(WQ_n)))
      BQ_n <-  WQ_n %>% group_by(ID_Line) %>% sample_n(1) %>% ungroup()
      
      #for lines which had no survivor queens in WQ_n to perpetuate the line, we will divide an existing line to recreate a new line with the same ID_Line as they have before
      
      Original_lines <- levels(WQ_n$ID_Line)
      Missing_lines <- Original_lines[!(Original_lines %in% Survived_lines)]
      
      if(length(Missing_lines) == 0){
        
      }else{#if there are missing lines
        Selected_BQs <- WQ_n %>%
          # top_n(NB_BQ, eval(parse(text=Q_SELECTION_CRITERION)) ) %>% # we order WQs
          filter(!(ID %in% BQ_n$ID)) %>% # We filter out alreade chosen BQs
          sample_n(length(Missing_lines), replace = FALSE) #we select length(Missing_lines) of the top surviving Qs of WQ_n that are not used already as BQs of surving lines
        BQ_n <- bind_rows(BQ_n, Selected_BQs)#we add these queens to BQ_n
        BQ_n$ID_Line[( (length(BQ_n$ID_Line)) - (length(Selected_BQs$ID_Line)) +1) : (length(BQ_n$ID_Line))] <- Missing_lines #We change their ID_Line to match with the ID_Line she is replacing
        # If a line had to be regenerated, the order of BQs in BQ_n will not be as ordered in D_n, which will create problems when generating offspring of these BQs. We therefore reorder BQs in BQ_n as they appear in D_n
        # order BQ_n as ID_mate in D_n
        # D_n_BQ <- D_n %>% filter(ID_mate %in% BQ_n$ID) # contains repeats
        # unique_BQ_IDs_as_in_D_n <- unique(D_n_BQ$ID_mate) # original order as in D_n    (unique ID_Qs)           
        # BQ_n <- BQ_n %>% arrange(factor(ID, levels = unique_BQ_IDs_as_in_D_n))
      }
      
      #ordering BQs as originally there were ordered in WQ_n
      BQ_n <- setorder(BQ_n, Original_order)
      WQ_n <- WQ_n %>% select(-Original_order)
      BQ_n <- BQ_n %>% select(-Original_order)
      
      return(BQ_n)
    } # end of Select_BQ for random selection
  }else{ #If Q selection is not random
    #For subsequent years, BQs are selected among existing lines
    Select_BQ <- function(WQ_n, NB_BQ, Q_SELECTION_CRITERION="Perf_T1 + Perf_T2"){
      
      Survived_lines <- levels(factor(WQ_n$ID_Line, levels = unique(WQ_n$ID_Line)))
      
      # We start with the BQ chosen because they are the only survivors of a family, as top_n ignores group without possible choice
      BQ_n_sole <-  WQ_n %>% group_by(ID_Line) %>% filter(n() == 1)
      BQ_n <-  WQ_n %>% group_by(ID_Line) %>% top_n(1, eval(parse(text=Q_SELECTION_CRITERION)) ) %>% ungroup() %>%
        bind_rows(BQ_n_sole) %>% 
        arrange(factor(ID, levels = unique(WQ_n$ID) ))
      
      #for lines which had no survivor queens in WQ_n to perpetuate the line, we will divide an existing line to recreate a new line with the same ID_Line as they have before
      Original_lines <- levels(WQ_n$ID_Line)
      Missing_lines <- Original_lines[!(Original_lines %in% Survived_lines)]
      
      if(length(Missing_lines) == 0){
        
      }else{#if there are missing lines
        Selected_BQs <- WQ_n %>%
          top_n(NB_BQ, eval(parse(text=Q_SELECTION_CRITERION)) ) %>% # we order WQs
          filter(!(ID %in% BQ_n$ID)) %>% # We filter out alreade chosen BQs
          sample_n(length(Missing_lines), replace = FALSE) #we select length(Missing_lines) of the top surviving Qs of WQ_n that are not used already as BQs of surving lines
        BQ_n <- bind_rows(BQ_n, Selected_BQs)#we add these queens to BQ_n
        BQ_n$ID_Line[( (length(BQ_n$ID_Line)) - (length(Selected_BQs$ID_Line)) +1) : (length(BQ_n$ID_Line))] <- Missing_lines #We change their ID_Line to match with the ID_Line she is replacing
        # If a line had to be regenerated, the order of BQs in BQ_n will not be as ordered in D_n, which will create problems when generating offspring of these BQs. We therefore reorder BQs in BQ_n as they appear in D_n
        # order BQ_n as ID_mate in D_n
        D_n_BQ <- D_n %>% filter(ID_mate %in% BQ_n$ID) # contains repeats
        unique_BQ_IDs_as_in_D_n <- unique(D_n_BQ$ID_mate) # original order as in D_n    (unique ID_Qs)           
        BQ_n <- BQ_n %>% arrange(factor(ID, levels = unique_BQ_IDs_as_in_D_n))
      }  
      
      return(BQ_n)
    }
  } #end of If Q selection is not random
} # end of With L selection

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#Select_PS_from_WTQ                                                                   #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# Creates Q_of_PS_n
if(PAT_LINE_SELECTION==TRUE & MAT_LINE_SELECTION==TRUE){
  Select_PS_from_WTQ <- function(W1_CandidateDPQ=WTQ_n_save_two, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_BQ, Criterion){ 
    #potential DPQs should be year-2 old
    
    if(Criterion=="Random"){
      Potential_Q_of_PS_n <- W1_CandidateDPQ %>%
        group_by(ID_BQ) %>%
        summarise(avg_Perf = mean(Perf_T2 + Perf_T2)) %>% #not really useful, but to make homogeneous to non-random sel
        sample_n(NB_BQ)#We order the potential lines represented by Qs mother of PSs (paternal grand-mother queens, representing the genetic lines)
    }else{ #If PGM selection is not random
      Potential_Q_of_PS_n <- W1_CandidateDPQ %>%
        group_by(ID_BQ) %>%
        summarise(avg_Perf = mean( eval(parse(text=Criterion)) )) %>%
        arrange(desc(avg_Perf))#We order the potential lines represented by Qs mother of PSs (paternal grand-mother queens, representing the genetic lines)
    } #end of if Criterion!="Random"
    
    j=0 #we initialize a counter
    
    Q_of_PS_n <- data.frame(matrix(ncol = 2, nrow = NB_BQ))
    colnames( Q_of_PS_n) <-   c("ID_BQ", "avg_Perf")
    for(i in 1: NB_BQ){#for each potential breeding line
      if(Potential_Q_of_PS_n$ID_BQ[i] %in% W2_CandidateDPQ$ID_BQ){#we check if at least one MTQ of this lineage survived during the winter
        k=i+j
        Q_of_PS_n$ID_BQ[k] <- Potential_Q_of_PS_n$ID_BQ[i]#if it's true, we select it
        Q_of_PS_n$avg_Perf[k] <- Potential_Q_of_PS_n$avg_Perf[i]
      }else{
        j=j-1#if it is not true, we will look for the second best and will inspect one more potential Qs to still get NB_PGM Qs as we wanted
        #j-1 is helpfull so that k stays the same index of line of Q_of_PS_n. At this line, for which there was no added Q because none of the descendants of this BQ survided, we want to add a descendant from the next best BQ line.
      }
    }
    Q_of_PS_n <- Q_of_PS_n[1:NB_PGM, ]#at the end we only keep the desired number of breeding lines for paternal lineage   
    
    #We then add to Q_of_PS_n the ID_Line of retained PSs, useful to pass it to WQ_n when PSs are attributed to each WQ
    Q_of_PS_n <- Q_of_PS_n %>% left_join(select(W2_CandidateDPQ, ID_BQ, ID_Line) %>% distinct(),  by = "ID_BQ" )#getting ID_line identifier in table WQ_n
    
    return(Q_of_PS_n)
  }
}#end of if(PAT_LINE_SELECTION==TRUE & MAT_LINE_SELECTION==TRUE)
if(PAT_LINE_SELECTION==TRUE & MAT_LINE_SELECTION==FALSE){
  Select_PS_from_WTQ <- function(W1_CandidateDPQ=WTQ_n_save_two, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_BQ, Criterion){ #Q in MPQ should be aged year-2
    
    if(Criterion=="Random"){
      Potential_Q_of_PS_n <- W1_CandidateDPQ %>%
        group_by(ID_BQ) %>%
        summarise(avg_Perf = mean(Perf_T2 + Perf_T2)) %>% #not really useful, but to make homogeneous to non-random sel
        sample_n(NB_BQ)#We order the potential lines represented by Qs mother of PSs (paternal grand-mother queens, representing the genetic lines)
    }else{ #If PGM selection is not random
      Potential_Q_of_PS_n <- W1_CandidateDPQ %>%
        group_by(ID_BQ) %>%
        summarise(avg_Perf = mean( eval(parse(text=Criterion)) )) %>%
        arrange(desc(avg_Perf))#We order the potential lines represented by Qs mother of PSs (paternal grand-mother queens, representing the genetic lines)
    } #end of if Criterion!="Random"
    
    j=0 #we initialize a counter
    
    Q_of_PS_n <- data.frame(matrix(ncol = 2, nrow = NB_BQ))
    colnames( Q_of_PS_n) <-   c("ID_BQ", "avg_Perf")
    for(i in 1: NB_BQ){#for each potential breeding line
      if(Potential_Q_of_PS_n$ID_BQ[i] %in% W2_CandidateDPQ$ID_BQ){#we check if at least one MTQ of this lineage survived during the winter
        k=i+j
        Q_of_PS_n$ID_BQ[k] <- Potential_Q_of_PS_n$ID_BQ[i]#if it's true, we select it
        Q_of_PS_n$avg_Perf[k] <- Potential_Q_of_PS_n$avg_Perf[i]
      }else{
        j=j-1#if it is not true, we will look for the second best and will inspect one more potential Qs to still get NB_PGM Qs as we wanted
        #j-1 is helpfull so that k stays the same index of line of Q_of_PS_n. At this line, for which there was no added Q because none of the descendants of this BQ survided, we want to add a descendant from the next best BQ line.
      }
    }
    Q_of_PS_n <- Q_of_PS_n[1:NB_PGM, ]#at the end we only keep the desired number of breeding lines for paternal lineage   
    
    #We then add to Q_of_PS_n the ID_Line of retained PSs, useful to pass it to WQ_n when PSs are attributed to each WQ
    
    return(Q_of_PS_n)
  }
}#end of if(PAT_LINE_SELECTION==TRUE & MAT_LINE_SELECTION==FALSE)

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Create_DPQ_n_blanck                                                            #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
if(PAT_LINE_SELECTION==TRUE){
  Create_DPQ_n_blanck <- function(W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS){
    DPQ_n <- data.frame(matrix(ncol = 8, nrow =(NB_PGM*NB_DPQ_PER_PS)))
    colnames(DPQ_n) <- c("ID_TQ", "ID_PGM", "BV_d1_TQ", "BV_m1_TQ", "BV_d2_TQ", "BV_m2_TQ", "F_DPQ", "ID_Line")
    
    DPQ_n$ID_Line <- factor(DPQ_n$ID_Line, levels = unique(W2_CandidateDPQ$ID_Line))
    
    return(DPQ_n)
  }
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Fill_DPQ_n                                                                     #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# La fonction Fill_DPQ_n va rechercher les lignées retenues dans Q_of_PS_n pour choisir
# les TQs avec les meilleures performances par lignées en tant que DPQ.
# S'il y avait par hasard, dû à la mortalité hivernale des TQ dans W2_CandidateDPQ,
# moins de représentantes d'une lignées que de DPQ de cette lignées souhaitées
# (paramètre NB_DPQ_PER_PS), des lignes répéterons les reines disponibles dans le tableau DPQ_n.
# On produira alors davantage de mâles sur les DPQs restantes pour cette lignée mâle
# pour avoir un nombre de faux-bourdons total par lignée égal pour contribuer à la prochaine génération.

# Suivant que l'on soit en PAT_LINE_SELECTION avec AVOID_INLINE ou non
# 2 formulation de la fonction Fill_DPQ_n sont possible, 
# la première ajoutant une colone ID_Line à DPQ_n et l'autre non:

# First case: if we want to avoid inline pairing
if(PAT_LINE_SELECTION==TRUE & PAIRING == "AVOID_INLINE"){
  Fill_DPQ_n <- function(Q_of_PS_n, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS, DPQ_n, Criterion){
    
    #preparaing DPQ_n table with the PGM IDs:
    for (j in 1:(NB_PGM)){
      for(k in 1: NB_DPQ_PER_PS){
        DPQ_n$ID_PGM[k+NB_DPQ_PER_PS*(j-1)] <- Q_of_PS_n$ID_BQ[j]
      }
    }
    
    chosen_DPQs <- DPQ_n[0, ] #we initialize a table keeping in memory DPQs chosen
    
    for (i in seq(from = 1, to = NB_PGM*NB_DPQ_PER_PS, by = NB_DPQ_PER_PS)){
      
      if(Criterion=="Random"){ #if the DPQs constituing a PS are chosen randomly
        Chosen_DPQ <- W2_CandidateDPQ %>% 
          filter(ID_BQ %in% DPQ_n$ID_PGM[i]) 
        NB_survived_Q_per_PS <- dim(Chosen_DPQ)[1] #We count how many TQ of this PGM survived the second wintering
        if(NB_survived_Q_per_PS > NB_DPQ_PER_PS){ #If there are more surviving potential DPQ of this PGM then wanted DPQ of this lineage, we choose randomly among them 
          Chosen_DPQ <- Chosen_DPQ %>%
            sample_n(NB_DPQ_PER_PS, replace = FALSE) 
        }
      }else{ #if there is non-random selection of the DPQs inside a PS
        Chosen_DPQ <- W2_CandidateDPQ %>% 
          filter(ID_BQ %in% DPQ_n$ID_PGM[i]) %>% 
          top_n(NB_DPQ_PER_PS, eval(parse(text=paste0(Criterion))) )#we pass the Criterion as a variable and not char chain using eval(parse())
      }
      
      #if there are less than NB_DPQ_PER_PS, for example 3, representatives of this lineage, the function will still work, getting the highest number of survived representatives of this lineage (2 or even only 1)
      
      DPQ_n$ID_TQ[i] <- Chosen_DPQ$ID[1]
      DPQ_n$BV_d1_TQ[i] <- Chosen_DPQ$BV_d1_TQ[1]
      DPQ_n$BV_m1_TQ[i] <- Chosen_DPQ$BV_m1_TQ[1]
      DPQ_n$BV_d2_TQ[i] <- Chosen_DPQ$BV_d2_TQ[1]
      DPQ_n$BV_m2_TQ[i] <- Chosen_DPQ$BV_m2_TQ[1]
      DPQ_n$F_DPQ[i] <- Chosen_DPQ$F_ID[1]
      DPQ_n$ID_Line[i] <- Chosen_DPQ$ID_Line[1]
      
      chosen_DPQs <- rbind(chosen_DPQs, Chosen_DPQ)
      
      if(NB_DPQ_PER_PS > 1){ #if there could be more than 1 TQ per lineage:
        for(j in 2:NB_DPQ_PER_PS){
          if((Chosen_DPQ$ID[j]) %in% NA){#if there there are less representative of this lineage than NB_DPQ_PER_PS, than the only DPQs chosable will be among the ones that survived
            
            Survivor <- sample_n(Chosen_DPQ, 1, replace = TRUE)
            
            DPQ_n$ID_TQ[i+j-1] <- Survivor$ID[1]
            DPQ_n$BV_d1_TQ[i+j-1] <- Survivor$BV_d1_TQ[1]
            DPQ_n$BV_m1_TQ[i+j-1] <- Survivor$BV_m1_TQ[1]
            DPQ_n$BV_d2_TQ[i+j-1] <- Survivor$BV_d2_TQ[1]
            DPQ_n$BV_m2_TQ[i+j-1] <- Survivor$BV_m2_TQ[1]
            DPQ_n$F_DPQ[i+j-1] <- Survivor$F_ID[1]
            DPQ_n$ID_Line[i+j-1] <- Survivor$ID_Line[1]
          }else{#if there are other survivors of this ancestry
            DPQ_n$ID_TQ[i+j-1] <- Chosen_DPQ$ID[j]
            DPQ_n$BV_d1_TQ[i+j-1] <- Chosen_DPQ$BV_d1_TQ[j]
            DPQ_n$BV_m1_TQ[i+j-1] <- Chosen_DPQ$BV_m1_TQ[j]
            DPQ_n$BV_d2_TQ[i+j-1] <- Chosen_DPQ$BV_d2_TQ[j]
            DPQ_n$BV_m2_TQ[i+j-1] <- Chosen_DPQ$BV_m2_TQ[j]
            DPQ_n$F_DPQ[i+j-1] <- Chosen_DPQ$F_ID[j]
            DPQ_n$ID_Line[i+j-1] <- Chosen_DPQ$ID_Line[j]
            
          }
        }
      }
      
    }
    
    missing_DPQs <- which(DPQ_n$ID_TQ %in% NA) #we store eventual lines IDs for missing DPQs. Missing DPQs appear when there are less DPQ lines that survided winter than NB_PGM (high NB_PGM compared to NB_TQ_PER_BQ ans high winter loss rates). In that case, these missing DPQs will be replaced by randomly chosen DPQs among those previously selected
    if(length(missing_DPQs)==0){#if there are no missing DPQs nothing is done
      
    }else{#if there are in fact missing DPQs
      for (i in missing_DPQs){#for every missing DPQ
        replacement_DPQ <- sample_n(chosen_DPQs, 1, replace = TRUE)#we randomly chose a DPQ among the previously selected DPQs to replace a missing DPQ
        DPQ_n$ID_TQ[i] <- replacement_DPQ$ID[1]
        DPQ_n$BV_d1_TQ[i] <- replacement_DPQ$BV_d1_TQ[1]
        DPQ_n$BV_m1_TQ[i] <- replacement_DPQ$BV_m1_TQ[1]
        DPQ_n$BV_d2_TQ[i] <- replacement_DPQ$BV_d2_TQ[1]
        DPQ_n$BV_m2_TQ[i] <- replacement_DPQ$BV_m2_TQ[1]
        DPQ_n$ID_PGM[i] <- replacement_DPQ$ID_BQ[1]
        DPQ_n$F_DPQ[i] <- replacement_DPQ$F_ID[1]
        DPQ_n$ID_Line[i] <- replacement_DPQ$ID_Line[1]
      }
    }
    
    return(DPQ_n)
  }
}

# second case: without avoiding inline mating
if(PAT_LINE_SELECTION==TRUE & PAIRING != "AVOID_INLINE"){
  Fill_DPQ_n <- function(Q_of_PS_n, W2_CandidateDPQ=W2_TQ_n, NB_PGM, NB_DPQ_PER_PS, DPQ_n, Criterion){
    
    #preparaing DPQ_n table with the PGM IDs:
    for (j in 1:(NB_PGM)){
      for(k in 1: NB_DPQ_PER_PS){
        DPQ_n$ID_PGM[k+NB_DPQ_PER_PS*(j-1)] <- Q_of_PS_n$ID_BQ[j]
      }
    }
    
    chosen_DPQs <- DPQ_n[0, ] #we initialize a table keeping in memory DPQs chosen
    
    for (i in seq(from = 1, to = NB_PGM*NB_DPQ_PER_PS, by = NB_DPQ_PER_PS)){
      
      if(Criterion=="Random"){ #if the DPQs constituing a PS are chosen randomly
        Chosen_DPQ <- W2_CandidateDPQ %>% 
          filter(ID_BQ %in% DPQ_n$ID_PGM[i]) 
        NB_survived_Q_per_PS <- dim(Chosen_DPQ)[1] #We count how many TQ of this PGM survived the second wintering
        if(NB_survived_Q_per_PS > NB_DPQ_PER_PS){ #If there are more surviving potential DPQ of this PGM then wanted DPQ of this lineage, we choose randomly among them 
          Chosen_DPQ <- Chosen_DPQ %>%
            sample_n(NB_DPQ_PER_PS, replace = FALSE) 
        }
      }else{ #if there is non-random selection of the DPQs inside a PS
        Chosen_DPQ <- W2_CandidateDPQ %>% 
          filter(ID_BQ %in% DPQ_n$ID_PGM[i]) %>% 
          top_n(NB_DPQ_PER_PS, eval(parse(text=paste0(Criterion))) )#we pass the Criterion as a variable and not char chain using eval(parse())
      }
      
      #if there are less than NB_DPQ_PER_PS, for example 3, representatives of this lineage, the function will still work, getting the highest number of survived representatives of this lineage (2 or even only 1)

      DPQ_n$ID_TQ[i] <- Chosen_DPQ$ID[1]
      DPQ_n$BV_d1_TQ[i] <- Chosen_DPQ$BV_d1_TQ[1]
      DPQ_n$BV_m1_TQ[i] <- Chosen_DPQ$BV_m1_TQ[1]
      DPQ_n$BV_d2_TQ[i] <- Chosen_DPQ$BV_d2_TQ[1]
      DPQ_n$BV_m2_TQ[i] <- Chosen_DPQ$BV_m2_TQ[1]
      DPQ_n$F_DPQ[i] <- Chosen_DPQ$F_ID[1]
      DPQ_n$ID_Line[i] <- Chosen_DPQ$ID_Line[1]
      
      chosen_DPQs <- rbind(chosen_DPQs, Chosen_DPQ)
      
      if(NB_DPQ_PER_PS > 1){ #if there could be more than 1 TQ per lineage:
        for(j in 2:NB_DPQ_PER_PS){
          if((Chosen_DPQ$ID[j]) %in% NA){#if there there are less representative of this lineage than NB_DPQ_PER_PS, than the only DPQs chosable will be among the ones that survived
            
            Survivor <- sample_n(Chosen_DPQ, 1, replace = TRUE)
            
            DPQ_n$ID_TQ[i+j-1] <- Survivor$ID[1]
            DPQ_n$BV_d1_TQ[i+j-1] <- Survivor$BV_d1_TQ[1]
            DPQ_n$BV_m1_TQ[i+j-1] <- Survivor$BV_m1_TQ[1]
            DPQ_n$BV_d2_TQ[i+j-1] <- Survivor$BV_d2_TQ[1]
            DPQ_n$BV_m2_TQ[i+j-1] <- Survivor$BV_m2_TQ[1]
            DPQ_n$F_DPQ[i+j-1] <- Survivor$F_ID[1]
            DPQ_n$ID_Line[i+j-1] <- Survivor$ID_Line[1]
          }else{#if there are other survivors of this ancestry
            DPQ_n$ID_TQ[i+j-1] <- Chosen_DPQ$ID[j]
            DPQ_n$BV_d1_TQ[i+j-1] <- Chosen_DPQ$BV_d1_TQ[j]
            DPQ_n$BV_m1_TQ[i+j-1] <- Chosen_DPQ$BV_m1_TQ[j]
            DPQ_n$BV_d2_TQ[i+j-1] <- Chosen_DPQ$BV_d2_TQ[j]
            DPQ_n$BV_m2_TQ[i+j-1] <- Chosen_DPQ$BV_m2_TQ[j]
            DPQ_n$F_DPQ[i+j-1] <- Chosen_DPQ$F_ID[j]
            DPQ_n$ID_Line[i+j-1] <- Chosen_DPQ$ID_Line[j]
            
          }
        }
      }
      
    }
    
    missing_DPQs <- which(DPQ_n$ID_TQ %in% NA) #we store eventual lines IDs for missing DPQs. Missing DPQs appear when there are less DPQ lines that survided winter than NB_PGM (high NB_PGM compared to NB_TQ_PER_BQ ans high winter loss rates). In that case, these missing DPQs will be replaced by randomly chosen DPQs among those previously selected
    if(length(missing_DPQs)==0){#if there are no missing DPQs nothing is done
      
    }else{#if there are in fact missing DPQs
      for (i in missing_DPQs){#for every missing DPQ

        replacement_DPQ <- sample_n(chosen_DPQs, 1, replace = TRUE)#we randomly chose a DPQ among the previously selected DPQs to replace a missing DPQ
        DPQ_n$ID_TQ[i] <- replacement_DPQ$ID[1]
        DPQ_n$BV_d1_TQ[i] <- replacement_DPQ$BV_d1_TQ[1]
        DPQ_n$BV_m1_TQ[i] <- replacement_DPQ$BV_m1_TQ[1]
        DPQ_n$BV_d2_TQ[i] <- replacement_DPQ$BV_d2_TQ[1]
        DPQ_n$BV_m2_TQ[i] <- replacement_DPQ$BV_m2_TQ[1]
        DPQ_n$ID_PGM[i] <- replacement_DPQ$ID_BQ[1]
        DPQ_n$F_DPQ[i] <- replacement_DPQ$F_ID[1]
        DPQ_n$ID_Line[i] <- replacement_DPQ$ID_Line[1]
      }
    }
    
    return(DPQ_n)
  }
}

#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Fill_DPQ_n                                                                     #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# Selection of DPQ_n directly (no previous family selection)
if(PAT_LINE_SELECTION==FALSE){
  Select_DPQs <- function  (W1_CandidateDPQ=WTQ_n_save_two, W2_CandidateDPQ=W2_TQ_n, NB_TQ, NB_DPQ, Criterion){ 
    
    NB_W1candidateDPQ <- nrow(W1_CandidateDPQ)
    
    if(Criterion=="Random"){ #if the DPQs are chosen randomlu
      Potential_DPQs <- sample_n(W1_CandidateDPQ, NB_W1candidateDPQ, replace = FALSE)#we randomly select NB_DPQ TQs among the ones that survived the two winters
    }else{
      Potential_DPQs <- arrange(W1_CandidateDPQ, desc( eval(parse(text=paste0(Criterion))) ))#we order MTQs on their performance score
    }
    
    DPQ_n <- data.frame(matrix(ncol = 3, nrow =(NB_DPQ)))#we create a Data frame which will contain infos about selected DPQs
    colnames(DPQ_n) <- c("ID_TQ", "ID_PGM", "F_DPQ")
    
    j=0 #we initialize some counters
    k=0
    
    for(i in 1:NB_TQ){
      if( k < NB_DPQ){#if the number of DPQs selected, that survived the two winters, is not yet the desired NB_DPQ
        if(Potential_DPQs$ID[i] %in% W2_CandidateDPQ$ID){#we check if this TQ, chosen in order of Performance, survived the second wintering
          
          k=i+j#if this TQ survived, we augment a counter
          DPQ_n$ID_TQ[k] <- Potential_DPQs$ID[i]#if this TQ survived, we select it
          DPQ_n$ID_PGM[k] <- Potential_DPQs$ID_BQ[i]
          DPQ_n$BV_d1_TQ[k] <- Potential_DPQs$BV_d1_TQ[i]
          DPQ_n$BV_m1_TQ[k] <- Potential_DPQs$BV_m1_TQ[i]
          DPQ_n$BV_d2_TQ[k] <- Potential_DPQs$BV_d2_TQ[i]
          DPQ_n$BV_m2_TQ[k] <- Potential_DPQs$BV_m2_TQ[i]
          DPQ_n$F_DPQ[k] <- Potential_DPQs$F_ID[i]
          
        }else{
          j=j-1#if it is not true (the TQ did not survive), we will look for the second best and will inspect one more potential TQ to still get NB_DPQ TQs as we wanted
        }
        
      }#end of if k < NB_DPQ
    }#end of for i
    
    #if we still have not enough DPQs in DPQ_n even tough we already took all TQs that survived 2nd wintering we will complete DPQ_n by letting randomly chosen TQs appear more than 1 time in DPQ_n
    
    missing_DPQs <- which(DPQ_n$ID_TQ %in% NA) 
    if(length(missing_DPQs)!=0){#if there are no missing DPQs nothing is done
      for (i in missing_DPQs){#for every missing DPQ
        replacement_DPQ <- sample_n(Potential_DPQs, 1, replace = TRUE)#we randomly chose a DPQ among the previously selected DPQs to replace a missing DPQ
        DPQ_n$ID_TQ[i] <- replacement_DPQ$ID[1]
        DPQ_n$BV_d1_TQ[i] <- replacement_DPQ$BV_d1_TQ[1]
        DPQ_n$BV_m1_TQ[i] <- replacement_DPQ$BV_m1_TQ[1]
        DPQ_n$BV_d2_TQ[i] <- replacement_DPQ$BV_d2_TQ[1]
        DPQ_n$BV_m2_TQ[i] <- replacement_DPQ$BV_m2_TQ[1]
        DPQ_n$ID_PGM[i] <- replacement_DPQ$ID_BQ[1]
        DPQ_n$F_DPQ[i] <- replacement_DPQ$F_ID[1]
      }
    }
    
    return(DPQ_n)
  }#end of function
}
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Assign_PGM_or_TQ_to_WQ                                                         #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# WQ_n will be assigned to a TQ (here as a DPQ) or to a PS (group of sister DPQs). 
# A record is kept in the pedigree.

# Distribution is semi-random: we choose at random which TQ/PS will supply the males for which Q,
# but we distribute the TQs/PSs to the different Qs as evenly as possible.
# In TQ siblings where fewer sisters have survived than NB_DPQ_PER_PS,
# the surviving TQs are represented several times over,
# as if multiplied to maintain the same number of DPQs per pseudo-sire.
# In practice, this means that these TQs produce more Ds for fertilization,
# producing the Ds that their non-surviving sister(s) should have produced. 
# It is therefore the PSs' participation in the Ds pool that is balanced, 
# rather than that of the TQs, even if, when mortality is low, it amounts to the same thing.

# For this semi-random choice, we proceed as follows:
# step 1: randomize the order of rows in the WQ_n data frame,
# step 2: select DPQ/PS (i.e. selected TQs) for the Qs in the WQ_n data frame,
# using a random draw. When all the TQs have been drawn,
# repeat the process until enough have been drawn to fill WQ_n.


# Case 1: Semi-random QxTQ couplings (DPQ contribution balanced)
if(PAIRING == "SEMI_RANDOM" & MATING_SCALE=="Q_x_DPQ"){
  Assign_PGM_or_TQ_to_WQ <- function(DPQ_n, WQ_n, NB_DPQ){# this function semi-randomly assigns TQs to newly born Qs. It is semi random because it assigns TQs to Qs in order to even the contribution of each DPQ ancestry, so that every DPQ ancestry in DPQ_n will be assigned to a similar number of Qs.
    
    NB_W1_Q <- length(WQ_n$ID)#number of Qs to be mated
    NB_DPQ <- length(DPQ_n$ID_TQ)#number of DPQs
    
    n_samples <- floor(NB_W1_Q/(NB_DPQ))#number of samples to do
    
    DPQ_n <- sample_n(DPQ_n, NB_DPQ, replace = FALSE)#randomisation of the DPQ_n table which was before ordered by ancestry
    
    PS_or_TQ_to_WQ_n <- DPQ_n #initialisation
    
    for (i in 1:n_samples){ 
      PS_or_TQ_to_WQ_n <- rbind(PS_or_TQ_to_WQ_n, DPQ_n %>% sample_n(NB_DPQ, replace = FALSE))
    }
    PS_or_TQ_to_WQ_n <- PS_or_TQ_to_WQ_n[1:NB_W1_Q,] #we adjust to the actual size of WQ_n (for cases for which it would be needed)
    
    WQ_n$ID_TQ <- PS_or_TQ_to_WQ_n$ID_TQ #on inscnrit ensuite les ID_TQ dans WQ_n
    WQ_n$ID_PGM_mated <- PS_or_TQ_to_WQ_n$ID_PGM # we add the column with PGM mated (ID_4a) to make it accessible to Update_Output_Q_from_WQ_n
    
    PS_or_TQ_to_WQ_n$ID_Q <- WQ_n$ID #we add the ID_Q column to PS_or_TQ_to_WQ_n for use in compelting the mates in the pedigree
    PS_or_TQ_to_WQ_n_and_WQ_n <- list(PS_or_TQ_to_WQ_n, WQ_n)
    
    return(PS_or_TQ_to_WQ_n_and_WQ_n)
  }
}


# Case 2: Totally random QxTQ couplings 
# Creation of another function choosing TQs for each Q of Q_n in a totally random way. The contributions of each TQ can therefore be unbalanced:
if(PAIRING == "RANDOM" & MATING_SCALE=="Q_x_DPQ"){
  Assign_PGM_or_TQ_to_WQ <- function(DPQ_n, WQ_n){# this function semi-randmoly assigns TQs to newly born Qs. It is semi random because it assigns TQs to Qs in order to even the contribution of each DPQ ancestry, so that every DPQ ancestry in DPQ_n will be assigned to a similar number of Qs.
    
    NB_W1_Q <- length(WQ_n$ID_Q)#number of Qs in WQ_n
    NB_DPQ <- length(DPQ_n$ID_TQ)#number of DPQs
    
    PS_or_TQ_to_WQ_n <- DPQ_n %>% sample_n(NB_W1_Q, replace = TRUE)#initialisation
    
    WQ_n$ID_TQ <- PS_or_TQ_to_WQ_n$ID_TQ #on inscnrit ensuite les ID_TQ dans WQ_n
    WQ_n$ID_PGM_mated <- PS_or_TQ_to_WQ_n$ID_PGM # we add the column with PGM mated (ID_4a) to make it accessible to Update_Output_Q_from_WQ_n
    
    PS_or_TQ_to_WQ_n$ID_Q <- WQ_n$ID #we add the ID_Q column to PS_or_TQ_to_WQ_n for use in compelting the mates in the pedigree
    PS_or_TQ_to_WQ_n_and_WQ_n <- list(PS_or_TQ_to_WQ_n, WQ_n)
    
    return(PS_or_TQ_to_WQ_n_and_WQ_n)
  }
}


# Case 3: QxTQ non-intra-lineage couplings 
# Creation of another function choosing TQs for each Q of WQ_n so as to avoid TQ x Q crossings within the same line,
# with balanced participation in the total pool of Ds for each PS:
if(PAIRING == "AVOID_INLINE" & MATING_SCALE=="Q_x_DPQ"){
  Assign_PGM_or_TQ_to_WQ <- function(DPQ_n, WQ_n){# this function semi-randmoly assigns TQs to newly born Qs. It is semi random because it assigns TQs to Qs in order to even the contribution of each DPQ ancestry, so that every DPQ ancestry in DPQ_n will be assigned to a similar number of Qs. In addition, after having distributed DPQs to Qs, we check for crosses using maternally and paternally the same line and reassing other DPQs in other to avoid inline crosses.
    
    NB_W1_Q <- length(WQ_n$ID)#number of Qs in WQ_n
    NB_DPQ <- length(DPQ_n$ID_TQ)#number of DPQs
    
    n_samples <- floor(NB_W1_Q/(NB_DPQ))#number of samples to do
    
    DPQ_n <- sample_n(DPQ_n, NB_DPQ, replace = FALSE)#randomisation of the DPQ_n table which was before ordered by ancestry
    
    PS_or_TQ_to_WQ_n <- DPQ_n %>% sample_n(NB_DPQ, replace = FALSE)#initialisation
    
    for (i in 1:n_samples){ 
      PS_or_TQ_to_WQ_n <- rbind(PS_or_TQ_to_WQ_n, DPQ_n %>% sample_n(NB_DPQ, replace = FALSE))
    }
    PS_or_TQ_to_WQ_n <- PS_or_TQ_to_WQ_n[1:NB_W1_Q,] #we adjust to the actual size of WQ_n (for cases for which it would be needed)
    
    WQ_n$ID_TQ <- PS_or_TQ_to_WQ_n$ID_TQ #on inscnrit ensuite les ID_TQ dans WQ_n
    WQ_n$ID_PGM_mated <- PS_or_TQ_to_WQ_n$ID_PGM # we add the column with PGM mated (ID_4a) to make it accessible to Update_Output_Q_from_WQ_n
    
    WQ_n$ID_Line_TQ <- PS_or_TQ_to_WQ_n$ID_Line
    
    #We will now draw again TQs for pairings that involve the same line maternally and paternally
    for (i in unique(PS_or_TQ_to_WQ_n$ID_Line) ){#for all lines represented in chosen DPQs
      
      #we filter out Qs of this line to check for possible inbred crosses
      Line_i_in_WQ_n <- which(WQ_n$ID_Line %in% i) #working only on possible inbred crosses, PGM line per PGM line
      TQs_to_change <- Line_i_in_WQ_n[which(WQ_n$ID_Line_TQ[Line_i_in_WQ_n] %in% i)]#getting lines where inbred crosses have been chosen for this PGM 
      
      #We will attribute to these Qs other TQs, from anothere line as theirs
      Possible_TQs <- DPQ_n %>% filter(!(ID_Line %in% i))
      
      chosen_TQs <- Possible_TQs %>% sample_n(length(TQs_to_change), replace = TRUE) #choosing replacement TQs
      
      WQ_n$ID_TQ[TQs_to_change] <- chosen_TQs %>% pull(ID_TQ) #filling WQ_n with the newly chosen TQs
      WQ_n$ID_PGM_mated[TQs_to_change] <- chosen_TQs %>% pull(ID_PGM) #filling WQ_n with the newly chosen PGMs
      WQ_n$ID_Line_TQ[TQs_to_change] <- chosen_TQs %>% pull(ID_Line) 
    }
    
    PS_or_TQ_to_WQ_n$ID_Q <- WQ_n$ID #we add the ID_Q column to PS_or_TQ_to_WQ_n for use in compelting the mates in the pedigree
    PS_or_TQ_to_WQ_n_and_WQ_n <- list(PS_or_TQ_to_WQ_n, WQ_n)
    
    return(PS_or_TQ_to_WQ_n_and_WQ_n)
  }
}


# Case 4: Semi-random QxPS couplings (balanced PS contribution)
if(PAIRING == "SEMI_RANDOM" & MATING_SCALE=="Q_x_PS"){
  Assign_PGM_or_TQ_to_WQ <- function(DPQ_n, WQ_n){# this function semi-randomly assigns TQs to newly born Qs. It is semi random because it assigns TQs to Qs in order to even the contribution of each DPQ ancestry, so that every DPQ ancestry in DPQ_n will be assigned to a similar number of Qs.
    
    NB_W1_Q <- length(WQ_n$ID)#number of Qs in WQ_n
    NB_DPQ <- length(DPQ_n$ID_TQ)#number of DPQs
    
    n_samples <- floor(NB_W1_Q/(NB_PGM))#number of samples to do
    
    Q_of_PS_n <- sample_n(Q_of_PS_n, NB_PGM, replace = FALSE)#randomisation of the Q_of_PS_n table which was before ordered by ancestry
    
    PGM_to_WQ_n <- Q_of_PS_n #initialisation
    
    for (i in 1:n_samples){ 
      PGM_to_WQ_n <- rbind(PGM_to_WQ_n, Q_of_PS_n %>% sample_n(NB_PGM, replace = FALSE))
    }
    PGM_to_WQ_n <- PGM_to_WQ_n[1:NB_W1_Q,] #we adjust to the actual size of WQ_n (for cases for which it would be needed)
    
    colnames(PGM_to_WQ_n)[1] <- "ID_PGM" #we rename it to be homogeneous to the table PS_or_TQ_to_WQ_n generated when MATING_SCALE is set at "Q_x_DPQ". This column is used to write according mated PGM IDs to WQs in table pedigree
    WQ_n$ID_PGM_mated <- PGM_to_WQ_n$ID_PGM #we write ID_PGMs in WQ_n
    
    PGM_to_WQ_n$ID_Q <- WQ_n$ID #we add the ID_Q column to PGM_to_WQ_n for use in compelting the mates in the pedigree
    PGM_to_WQ_n_and_WQ_n <- list(PGM_to_WQ_n, WQ_n)
    
    return(PGM_to_WQ_n_and_WQ_n)
  }
}


# Case 5: Totally random QxPS couplings (contribution of PSs a priori unbalanced)
# Same as above, but this time the choice of couplings is made on a QxPS scale rather than on a QxDPQ scale.
if(PAIRING == "RANDOM" & MATING_SCALE=="Q_x_PS"){
  Assign_PGM_or_TQ_to_WQ <- function(DPQ_n, WQ_n){# this function semi-randmoly assigns TQs to newly born Qs. It is semi random because it assigns TQs to Qs in order to even the contribution of each DPQ ancestry, so that every DPQ ancestry in DPQ_n will be assigned to a similar number of Qs.
    
    NB_W1_Q <- length(WQ_n$ID)#number of Qs in WQ_n
    NB_DPQ <- length(DPQ_n$ID_TQ)#number of DPQs
    
    PGM_to_WQ_n <- Q_of_PS_n %>% sample_n(NB_W1_Q, replace = TRUE)#initialisation
    
    colnames(PGM_to_WQ_n)[1] <- "ID_PGM" #we rename it to be homogeneous to the table PS_or_TQ_to_WQ_n generated when MATING_SCALE is set at "Q_x_DPQ". This column is used to write according mated PGM IDs to WQs in table pedigree
    WQ_n$ID_PGM_mated <- PGM_to_WQ_n$ID_PGM #on inscnrit ensuite les ID_TQ dans WQ_n
    
    PGM_to_WQ_n$ID_Q <- WQ_n$ID #we add the ID_Q column to PGM_to_WQ_n for use in compelting the mates in the pedigree
    PGM_to_WQ_n_and_WQ_n <- list(PGM_to_WQ_n, WQ_n)
    
    return(PGM_to_WQ_n_and_WQ_n)
  }
}


# Case 6: Semi-random QxPS mating avoiding intra-line mating (Q and PS from the same maternal line) (balanced PS contribution)
if(PAIRING == "AVOID_INLINE" & MATING_SCALE=="Q_x_PS" & MAT_LINE_SELECTION == TRUE){
  Assign_PGM_or_TQ_to_WQ <- function(DPQ_n, WQ_n){# this function semi-randmoly assigns TQs to newly born Qs. It is semi random because it assigns TQs to Qs in order to even the contribution of each DPQ ancestry, so that every DPQ ancestry in DPQ_n will be assigned to a similar number of Qs. In addition, after having distributed DPQs to Qs, we check for crosses using maternally and paternally the same line and reassing other DPQs in other to avoid inline crosses.
    
    NB_W1_Q <- length(WQ_n$ID)#number of Qs in WQ_n
    NB_DPQ <- length(DPQ_n$ID_TQ)#number of DPQs
    
    n_samples <- floor(NB_W1_Q/(NB_PGM))#number of samples to do
    
    Q_of_PS_n <- sample_n(Q_of_PS_n, NB_PGM, replace = FALSE)#randomisation of the Q_of_PS_n table which was before ordered by ancestry
    
    PGM_to_WQ_n <- Q_of_PS_n #initialisation
    
    for (i in 1:n_samples){ 
      PGM_to_WQ_n <- rbind(PGM_to_WQ_n, Q_of_PS_n %>% sample_n(NB_PGM, replace = FALSE))
    }
    PGM_to_WQ_n <- PGM_to_WQ_n[1:NB_W1_Q,] #we adjust to the actual size of WQ_n (for cases for which it would be needed)
    
    colnames(PGM_to_WQ_n)[1] <- "ID_PGM" #we rename it to be homogeneous to the table PS_or_TQ_to_WQ_n generated when MATING_SCALE is set at "Q_x_DPQ". This column is used to write according mated PGM IDs to WQs in table pedigree
    WQ_n$ID_PGM_mated <- PGM_to_WQ_n$ID_PGM #on inscnrit ensuite les ID_TQ dans WQ_n
    
    WQ_n$ID_Line_PGM <- PGM_to_WQ_n$ID_Line
    
    #We will now draw again PSs for matings that involve the same line maternally and paternally
    for (i in unique(PGM_to_WQ_n$ID_Line) ){#for all lines represented in chosen PSs
      
      #we filter out Qs of this line to check for possible inbred crosses
      Line_i_in_WQ_n <- which(WQ_n$ID_Line %in% i) #working only on possible inbred crosses, PGM line per PGM line
      PSs_to_change <- Line_i_in_WQ_n[which(WQ_n$ID_Line_PGM[Line_i_in_WQ_n] %in% i)]#getting lines where inbred crosses have been chosen for this PGM 
      
      #We will attribute to these Qs other PSs, from another line as theirs
      Possible_PSs <- Q_of_PS_n %>% filter(!(ID_Line %in% i))
      
      chosen_PSs <- Possible_PSs %>% sample_n(length(PSs_to_change), replace = TRUE) #choosing replacement PSs
      
      WQ_n$ID_PGM_mated[PSs_to_change] <- chosen_PSs %>% pull(ID_PGM) #filling WQ_n with the newly chosen PSs
      WQ_n$ID_Line_PGM[PSs_to_change] <- chosen_PSs %>% pull(ID_Line) 
    }
    
    PGM_to_WQ_n$ID_Q <- WQ_n$ID #we add the ID_Q column to PGM_to_WQ_n for use in compelting the mates in the pedigree
    PGM_to_WQ_n_and_WQ_n <- list(PGM_to_WQ_n, WQ_n)
    
    return(PGM_to_WQ_n_and_WQ_n)
  }
}


#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
#      Produce_and_assign_D                                                               #
#°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
# Two cases are possible: either each WQ virgin queen is mated with a single DPQ (parameter MATING_SCALE == Q_x_DPQ),
# i.e. all the drones that fertilize her come from a single DPQ, 
# or the queens are mated at Pseudo-sire level and the WQs can therefore be fertilized
# by drones coming from sister NB_DPQ_PER_PS belonging to a PS.

# Case 1: MATING_SCALE == Q_x_DPQ
if(MATING_SCALE=="Q_x_DPQ"){
  Produce_and_assign_D <- function (Q_to_mate=WQ_n, DPQ_n, NB_D, NB_D_FECUNDATING, VAR_A_meiose, year){#W1_CandidateDPQ must contain queens of year n-2, WQ_n of the current year
    
    D_n <- data.frame(matrix(ncol = 6, nrow = NB_D))
    colnames(D_n) <- c("ID", "BV_d1_D", "BV_m1_D", "BV_d2_D", "BV_m2_D", "ID_mate")
    D_n <- as_tibble(D_n)
    #initialisation of D_n : inscription of the Q ID that will be fecundated and the drone's IDs that will fecundate those Qs.
    # Creating vector of W1 Qs, repeated NB_D_FECUNDATING times before coming to the next W1 Q
    D_n$ID_mate <- unlist(map2(Q_to_mate$ID, NB_D_FECUNDATING, rep))
    
    # Naming the Ds
    D_n[, 1] <- paste0("D", 1:NB_D, "_", year)
    
    #Getting BVs for the drones:
    
    D_n <- D_n %>% left_join(select(Q_to_mate, ID, ID_TQ),  by = c("ID_mate"="ID") )#getting TQ identifier in table WQ_n
    D_n <- D_n %>% left_join(select(DPQ_n, ID_TQ, BV_d1_TQ, BV_m1_TQ, BV_d2_TQ, BV_m2_TQ, F_DPQ) %>% distinct(),  by = "ID_TQ" )#Getting BVs and F coef of TQs from DPQ_n
    
    phi_D <- rmvnorm(n=NB_D, sigma=VAR_A_meiose)#Variance of mendelian sampling terms is 1/4*VAR_A
    phi_D_d1 <- phi_D[, 2]#mendelian sampling terms of direct effects trait 1
    phi_D_m1 <- phi_D[, 1]#mendelian sampling terms of maternal effects trait 1
    phi_D_d2 <- phi_D[, 4]#mendelian sampling terms of direct effects trait 2
    phi_D_m2 <- phi_D[, 3]#mendelian sampling terms of maternal effects trait 2
    
    D_n <- cbind(D_n, phi_D_d1, phi_D_m1, phi_D_d2, phi_D_m2)
    
    D_n <-  D_n %>% mutate(BV_d1_D = 0.5*BV_d1_TQ + sqrt(1 - F_DPQ ) * phi_D_d1,
                           BV_m1_D = 0.5*BV_m1_TQ + sqrt(1 - F_DPQ ) * phi_D_m1,
                           BV_d2_D = 0.5*BV_d2_TQ + sqrt(1 - F_DPQ ) * phi_D_d2,
                           BV_m2_D = 0.5*BV_m2_TQ + sqrt(1 - F_DPQ ) * phi_D_m2)
    
    D_n <- select(D_n, ID, BV_d1_D, BV_m1_D, BV_d2_D, BV_m2_D, ID_mate, ID_TQ)
    
    return(D_n)
  }
}


# Case 2: MATING_SCALE == Q_x_PS
if(MATING_SCALE=="Q_x_PS"){
  Produce_and_assign_D <- function (Q_to_mate, DPQ_n, NB_D, NB_D_FECUNDATING, VAR_A_meiose, year){ #Q_to_mate of the current year
    
    D_n <- data.frame(matrix(ncol = 7, nrow = NB_D))
    colnames(D_n) <- c("ID", "BV_d1_D", "BV_m1_D", "BV_d2_D", "BV_m2_D", "ID_mate", "ID_TQ")
    
    D_n <- as_tibble(D_n)
    #initialisation of D_n : inscription of the Q ID that will be fecundated and the drone's IDs that will fecundate those Qs.
    # Creating vector of W1 Qs, repeated NB_D_FECUNDATING times before coming to the next W1 Q
    D_n$ID_mate <- unlist(map2(Q_to_mate$ID, NB_D_FECUNDATING, rep))
    
    # Naming the Ds
    D_n[, 1] <- paste0("D", 1:NB_D, "_", year)
    
    #Getting BVs for the drones:
    D_n <- D_n %>% left_join(select(Q_to_mate, ID, ID_PGM_mated),  by = c("ID_mate"="ID") )#getting TQ identifier in table Q_to_mate
    
    #We will now draw DPQs as mothers od drones. For each ID_PGM in D_n we will draw DPQs and add these drawn IDs to D_n
    for (i in unique(D_n$ID_PGM_mated) ){#for each ID_PGM in D_n
      
      which_Q <- which(D_n$ID_PGM_mated == i)#line of D_n where we will draw DPQs of a same PS (having the same PGM)
      
      DPQ_i <- DPQ_n %>% filter(ID_PGM %in% i) %>% select(ID_TQ) %>% sample_n(length(which_Q), replace = TRUE)#we draw ramdomly IDs of DPQs for Qs mated to this same PS i
      
      D_n$ID_TQ[which_Q] <- DPQ_i %>% pull(ID_TQ)#we add these ramdomly chosen DPQs to D_n
    }
    
    D_n <- D_n %>% left_join(select(DPQ_n, ID_TQ, BV_d1_TQ, BV_m1_TQ, BV_d2_TQ, BV_m2_TQ, F_DPQ) %>% distinct(),  by = "ID_TQ" )#Getting BVs and F coef of TQs from DPQ_n
    
    phi_D <- rmvnorm(n=NB_D, sigma=VAR_A_meiose)#Variance of mendelian sampling terms is 1/4*VAR_A
    phi_D_d1 <- phi_D[, 2]#mendelian sampling terms of direct effects trait 1
    phi_D_m1 <- phi_D[, 1]#mendelian sampling terms of maternal effects trait 1
    phi_D_d2 <- phi_D[, 4]#mendelian sampling terms of direct effects trait 2
    phi_D_m2 <- phi_D[, 3]#mendelian sampling terms of maternal effects trait 2
    
    D_n <- cbind(D_n, phi_D_d1, phi_D_m1, phi_D_d2, phi_D_m2)
    
    D_n <-  D_n %>% mutate(BV_d1_D = 0.5*BV_d1_TQ + sqrt(1 - F_DPQ ) * phi_D_d1,
                           BV_m1_D = 0.5*BV_m1_TQ + sqrt(1 - F_DPQ ) * phi_D_m1,
                           BV_d2_D = 0.5*BV_d2_TQ + sqrt(1 - F_DPQ ) * phi_D_d2,
                           BV_m2_D = 0.5*BV_m2_TQ + sqrt(1 - F_DPQ ) * phi_D_m2)
    
    D_n <- select(D_n, ID, BV_d1_D, BV_m1_D, BV_d2_D, BV_m2_D, ID_mate, ID_TQ)
    
    return(D_n)
  }
}