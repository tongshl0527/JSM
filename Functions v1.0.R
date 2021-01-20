preprocess <- function(inputData) {
  # remove N/A
  df <- na.omit(inputData)
  # sort survival years
  df <- df[order(df$OVERALL_SURVIVAL_TIME_YEARS), ]
  print("Success: All N/A records are removed")
  return(df)
}

get_batchList <- function(input_df) {
  table <- count(input_df$BATCH)
  batches <- table$x
  return(batches)
}

get_best_radius_new <- function(input_data, testmethod = "chitest") {
  pvalue_list_radii <- list() # a list to store each radius's pvalue
  #create a dataframe to store all combinations and label them with an extra column 
  comb_df <- data.frame(matrix(ncol = 6))
  colnames(comb_df) <- colnames(input_data)
  colnames(comb_df)[5] <- "Stat0_radius"
  colnames(comb_df)[6] <- "Stat1_radius"
  
  actDist_0 <- get_actDist(input_data, 0)
  actDist_1 <- get_actDist(input_data, 1)
  
  radii_0 <- get_radiusRange(input_data,0) # get radii
  radii_1 <- get_radiusRange(input_data,1) # get radii
  
  #check if two radii are the same 
  if (all(radii_0 %in% radii_1) == TRUE) {
    radii <- radii_0
    count <- 0
    for (radius in radii) {
      count <- count + 1
      
      freqMatrix_0 <- get_freqMatrix(input_data, 0, radius)
      freqMatrix_1 <- get_freqMatrix(input_data, 1, radius)
      
      expVal_0 <- get_expVal(freqMatrix_0, actDist_0)
      expVal_1 <- get_expVal(freqMatrix_1, actDist_1)
      
      if (is.null(expVal_0)|is.null(expVal_1)==TRUE) {
        next
      } else {
        bestComb_all_intervals_0 <-pick_test_and_find_bestComb_all_intervals(freqMatrix_0, expVal_0, testmethod)
        bestComb_all_intervals_1 <-pick_test_and_find_bestComb_all_intervals(freqMatrix_1, expVal_1, testmethod)
        
        sampled_df_0 <- sampling_of_bestComb (input_data, bestComb_all_intervals_0, 0, radius)
        if (empty(sampled_df_0)== TRUE) {
          comb_df <- comb_df 
        } else {
          sampled_df_0$Stat0_radius <- radius
          sampled_df_0$Stat1_radius <- NA
          comb_df <- rbind(comb_df, sampled_df_0)
        }
        
        sampled_df_1 <- sampling_of_bestComb (input_data, bestComb_all_intervals_1, 1, radius)
        if (empty(sampled_df_1)== TRUE) {
          comb_df <- comb_df 
        } else {
          sampled_df_1$Stat0_radius <- NA
          sampled_df_1$Stat1_radius <- radius
          comb_df <- rbind(comb_df, sampled_df_1)
        }
      }
    }
  } else {
    #radii lengths are not the same
    #radii 0
    count <- 0
    for (radius in radii_0) {
      count <- count + 1
      freqMatrix_0 <- get_freqMatrix(input_data, 0, radius)
      expVal_0 <- get_expVal(freqMatrix_0, actDist_0)
      bestComb_all_intervals_0 <-pick_test_and_find_bestComb_all_intervals(freqMatrix_0, expVal_0, testmethod)
      sampled_df_0 <- sampling_of_bestComb (input_data, bestComb_all_intervals_0, 0, radius)
      if (empty(sampled_df_0)== TRUE) {
        comb_df <- comb_df 
      } else {
        sampled_df_0$Stat0_radius <- radius
        sampled_df_0$Stat1_radius <- NA
        comb_df <- rbind(comb_df, sampled_df_0)
      }
    }
    #radii 1
    count <- 0
    for (radius in radii_1) {
      count <- count + 1
      freqMatrix_1 <- get_freqMatrix(input_data, 1, radius)
      expVal_1 <- get_expVal(freqMatrix_1, actDist_1)
      bestComb_all_intervals_1 <-pick_test_and_find_bestComb_all_intervals(freqMatrix_1, expVal_1, testmethod)
      sampled_df_1 <- sampling_of_bestComb (input_data, bestComb_all_intervals_1, 1, radius)
      if (empty(sampled_df_1)== TRUE) {
        comb_df <- comb_df 
      } else {
        sampled_df_1$Stat0_radius <- NA
        sampled_df_1$Stat1_radius <- radius
        comb_df <- rbind(comb_df, sampled_df_1)
      }
    }
  }
  #stat0
  stat0 <- comb_df[,-6]
  stat0 <- preprocess(stat0)
  stat0_radii <- sort(unique(comb_df$Stat0_radius))
  stat0_list <- list()
  count <- 0
  for (radius in stat0_radii) {
    count <- count +1
    stat0_list[[count]] <- which(comb_df$Stat0_radius == radius)
    names(stat0_list)[count] <- radius
  }
  #stat1
  stat1 <- comb_df[,-5]
  stat1 <- preprocess(stat1)
  stat1_radii <- sort(unique(comb_df$Stat1_radius))
  stat1_list <- list()
  count <- 0
  for (radius in stat1_radii) {
    count <- count +1
    stat1_list[[count]] <- which(stat1$Stat1_radius == radius)
    names(stat1_list)[count] <- radius
  }
  #generate combinations based on the index in dataframe
  combination <- expand.grid(stat0_radii,stat1_radii)
  
  #pvalue to be stored
  p_list <- list()
  for (i in 1:nrow(combination)) {
    var1 <- combination[i,1]
    var2 <- combination[i,2]
    
    df_1 <- comb_df[which(comb_df$Stat0_radius==var1),]
    df_2 <- comb_df[which(comb_df$Stat1_radius==var2),]
    df <- rbind(df_1, df_2)
    pvalue<- get_survPlot(df)[2]
    p_list[i] <- pvalue
    names(p_list)[i] <- toString(combination[i,])
  }
  
  maxPval <- which.max(p_list)
  Pvalue <- p_list[[maxPval]]
  best_comb <- combination[maxPval,]
  best_comb[3] <- Pvalue
  names(best_comb)[1] <- "Stat0_radius"
  names(best_comb)[2] <- "Stat1_radius"
  names(best_comb)[3] <- "Pvalue" 
  #return pairs with the highest p value
  return(best_comb)
}

compareBatch_All <- function(input_df) {
  ptable <- pairwise_survdiff(Surv(OVERALL_SURVIVAL_TIME_YEARS,OVERALL_SURVIVAL_EVENT) ~ BATCH, data = input_df)
  return(ptable)
}

get_survPlot <- function(input_df) {
  surv_object <- Surv(input_df$OVERALL_SURVIVAL_TIME_YEARS, input_df$OVERALL_SURVIVAL_EVENT)
  sf.BATCH <- surv_fit(surv_object ~ input_df$BATCH, data = input_df)
  return(surv_pvalue(sf.BATCH)) # Default is "survdiff" (or "log-rank")
}

get_radiusRange <- function(input_df, #this input is specific for a event status 
                            event = 0) {
  input_df <- subset(input_df, input_df$OVERALL_SURVIVAL_EVENT == 0)
  distance_matrix <- as.matrix(dist(input_df$OVERALL_SURVIVAL_TIME_YEARS))
  sorted_disMatrix <- apply(distance_matrix, 1, sort)
  radius <- apply(sorted_disMatrix, 1, min)
  radius <- unique(radius)[2:length(radius)]
  radiusRange <- radius[!is.na(radius)]
  return(radiusRange)
}

get_freqMatrix <- function(input_df, eventStat = 0, input_radius) {
  dataframe <- subset(input_df, OVERALL_SURVIVAL_EVENT == eventStat)
  batchList <- get_batchList(input_df)
  size_segment <- round(1 + max(dataframe$OVERALL_SURVIVAL_TIME_YEARS) / input_radius)
  
  freq_Matrix <- matrix(0, ncol = size_segment, nrow = length(batchList))
  rownames(freq_Matrix) <- batchList
  colnames(freq_Matrix) <- paste("interval", c(1:size_segment), sep = "_")
  
  for (i in 1:size_segment) {
    # subset out the interval
    # based on interval, select the right range of years. Batch is a list for all batches
    segment <- subset(dataframe, OVERALL_SURVIVAL_TIME_YEARS < i * input_radius & OVERALL_SURVIVAL_TIME_YEARS >= (i - 1) * input_radius)
    freq_table <- count(segment$BATCH)
    row.names(freq_table) <- freq_table$x
    freq_Matrix[rownames(freq_table), i] <- freq_table$freq
  }
  return(freq_Matrix)
}

get_actDist <- function(input_df, eventStat = 0) {
  # count total samples
  freq_total <- count(input_df, c("BATCH", "OVERALL_SURVIVAL_EVENT"))
  # count the total samples in event 0
  freq_event <- freq_total[freq_total$OVERALL_SURVIVAL_EVENT == eventStat, ]
  sum_event <- sum(freq_event$freq)
  freq_event <- transform(freq_event, proportion = freq / sum_event)
  # standard ratio of each batch is stored in a list
  return(freq_event)
}

get_expval_at_interval <- function(input_expVal, input_at_Interval = 1) {
  valid_intervals_in_expVal <- parse_number(colnames(input_expVal))
  idx <- match(paste("interval", input_at_Interval, sep = "_"),colnames(input_expVal))
  expVal <- input_expVal[,idx,drop=FALSE]
  return(expVal)
}

get_expVal <- function(input_freqMatrix, input_actDist) {
  #remove columns that are zeros only
  input_freqMatrix <- input_freqMatrix[, colSums(input_freqMatrix != 0) > 0,drop=FALSE]
  # total number of event 0 in this segment
  n_col <- ncol(input_freqMatrix)
  expVal_matrix <- matrix(0,nrow(input_freqMatrix),ncol(input_freqMatrix))
  rownames(expVal_matrix) <- rownames(input_freqMatrix)
  colnames(expVal_matrix) <- colnames(input_freqMatrix)
  for (i in 1:n_col) {
    sumCol <- sum(input_freqMatrix[, i])
    for (j in input_actDist$BATCH) {
      expVal <- round(input_actDist[input_actDist$BATCH==j,]$proportion *sumCol)
      expVal_matrix[j, i] <- expVal
    }
  }
  #remove columns with zeros
  expVal_matrix <- expVal_matrix[, colSums(expVal_matrix != 0) > 0, drop=FALSE]
  if(dim(expVal_matrix)[2] >0 ) {
    return(expVal_matrix)
  } else {
    return(NULL)
  }
}

get_permutation <- function(input_freqMatrix, input_radiusInterval = 1) {
  #convert input_radiusInterval to colname 
  idx <- match(paste("interval", input_radiusInterval, sep = "_"),colnames(input_freqMatrix))
  column <- input_freqMatrix[,idx,drop=FALSE]
  
  a1 <- sapply(column, function(x) {c(0:x)})
  a2 <- expand.grid(a1)
  
  if(sum(a2)==0)
  {
    a2 = t(a2)
  }
  colnames(a2) <- rownames(input_freqMatrix)
  rownames(a2) <- paste("interval", input_radiusInterval,c(1:nrow(a2)), sep = "_")
  return(t(a2))
} #get permutation for one interval!!

pick_test_and_find_bestComb_all_intervals <- function(input_freqMatrix, input_expVal, testmethod) {
  if (testmethod == "chitest") {
    bestComb_all_intervals <- find_bestComb_chisq_all_intervals_new(input_freqMatrix, input_expVal)
  }
  else if (testmethod == "cortest") {
    bestComb_all_intervals <- find_bestComb_cor_all_intervals_new(input_freqMatrix, input_expVal)
  } else if (testmethod == "fishertest") {
    bestComb_all_intervals <- find_bestComb_fisher_all_intervals_new(input_freqMatrix, input_expVal)
  }
  return(bestComb_all_intervals)
}

find_bestComb_chisq_at_interval <- function(input_freqMatrix,  #remove columns with 0 
                                            input_expVal, #remove columns with 0 
                                            input_at_Interval = 1 #apply intervals that only exist in expval!
) {
  #convert input_radiusInterval to colname 
  idx <- match(paste("interval", input_at_Interval, sep = "_"),colnames(input_freqMatrix))
  #select columns match to expval
  freqMatrix <- input_freqMatrix[,idx,drop=FALSE]
  #interval_freqMatrix<- input_freqMatrix[,input_at_Interval,drop=FALSE]
  perm_at_interval <- get_permutation(freqMatrix,input_at_Interval)
  expval_at_interval <- get_expval_at_interval(input_expVal,input_at_Interval)
  pvalues <- list()
  for (i in 1:ncol(perm_at_interval)) {
    if (length(unique(perm_at_interval[, i])) == 1) {
      pvalues[[i]] <- 0
    } else if (length(unique(expval_at_interval))== 1) {
      pvalues[[i]] <- 0
    } else {
      chitest <- chisq.test(perm_at_interval[, i], expval_at_interval) #error when running event =0, check permdf and expval if there are missing values
      pvalue <- chitest$p.value
      pvalues[[i]] <- pvalue
    }
  }
  maxPval <- which.max(as.numeric(unlist(pvalues)))
  Pvalue <- pvalues[[maxPval]]
  print(Pvalue)
  if (all(perm_at_interval[, maxPval]==0)==TRUE) {
    return(NULL)
  } else {
    return(perm_at_interval[, maxPval,drop=FALSE])
  }
}

find_bestComb_chisq_all_intervals_new <- function(input_freqMatrix, input_expVal) {
  #filter the freqMatrix with only the matched one
  expVal_colnames <- colnames(input_expVal)
  freqMatrix <- input_freqMatrix[, expVal_colnames, drop = FALSE]
  
  bestComb_all_intervals <- NULL
  intervals <- parse_number(colnames(input_expVal))
  
  for (interval in intervals) {
    #at_interval <- parse_number(colnames(input_expVal[,interval,drop=FALSE]))
    # give expval
    bestComb <- find_bestComb_chisq_at_interval(freqMatrix,input_expVal,interval)
    #bestComb_all_intervals[, interval] <- bestComb
    bestComb_all_intervals <- cbind(bestComb_all_intervals, bestComb)
  }
  colnames(bestComb_all_intervals) <- colnames(input_expVal)
  return(bestComb_all_intervals)
}


find_bestComb_cor_at_interval <- function(input_freqMatrix,  #remove columns with 0 
                                            input_expVal, #remove columns with 0 
                                            input_at_Interval = 1 #apply intervals that only exist in expval!
) {
  #convert input_radiusInterval to colname 
  idx <- match(paste("interval", input_at_Interval, sep = "_"),colnames(input_freqMatrix))
  #select columns match to expval
  freqMatrix <- input_freqMatrix[,idx,drop=FALSE]
  #interval_freqMatrix<- input_freqMatrix[,input_at_Interval,drop=FALSE]
  perm_at_interval <- get_permutation(freqMatrix,input_at_Interval)
  expval_at_interval <- get_expval_at_interval(input_expVal,input_at_Interval)
  pvalues <- list()
  for (i in 1:ncol(perm_at_interval)) {
    if (length(unique(perm_at_interval[, i])) == 1) {
      pvalues[[i]] <- 0
    } else if (length(unique(expval_at_interval))== 1) {
      pvalues[[i]] <- 0
    } else {
      cortest <- cor.test(perm_at_interval[, i], expval_at_interval) #error when running event =0, check permdf and expval if there are missing values
      pvalue <- cortest$p.value
      pvalues[[i]] <- pvalue
    }
  }
  maxPval <- which.max(as.numeric(unlist(pvalues)))
  Pvalue <- pvalues[[maxPval]]
  print(Pvalue)
  if (all(perm_at_interval[, maxPval]==0)==TRUE) {
    return(NULL)
  } else {
    return(perm_at_interval[, maxPval,drop=FALSE])
  }
}

find_bestComb_cor_all_intervals_new <- function(input_freqMatrix, input_expVal) {
  #filter the freqMatrix with only the matched one
  expVal_colnames <- colnames(input_expVal)
  freqMatrix <- input_freqMatrix[, expVal_colnames, drop = FALSE]
  
  bestComb_all_intervals <- NULL
  intervals <- parse_number(colnames(input_expVal))
  
  for (interval in intervals) {
    # give expval
    bestComb <- find_bestComb_cor_at_interval(freqMatrix,input_expVal,interval)
    bestComb_all_intervals <- cbind(bestComb_all_intervals, bestComb)
  }
  colnames(bestComb_all_intervals) <- colnames(input_expVal)
  return(bestComb_all_intervals)
}

find_bestComb_fisher_at_interval <- function(input_freqMatrix,  #remove columns with 0 
                                          input_expVal, #remove columns with 0 
                                          input_at_Interval = 1 #apply intervals that only exist in expval!
) {
  #convert input_radiusInterval to colname 
  idx <- match(paste("interval", input_at_Interval, sep = "_"),colnames(input_freqMatrix))
  #select columns match to expval
  freqMatrix <- input_freqMatrix[,idx,drop=FALSE]
  perm_at_interval <- get_permutation(freqMatrix,input_at_Interval)
  expval_at_interval <- get_expval_at_interval(input_expVal,input_at_Interval)
  pvalues <- list()
  for (i in 1:ncol(perm_at_interval)) {
    if (length(unique(perm_at_interval[, i])) == 1) {
      pvalues[[i]] <- 0
    } else if (length(unique(expval_at_interval))== 1) {
      pvalues[[i]] <- 0
    } else {
      fishertest <- fisher.test(perm_at_interval[, i], expval_at_interval) #error when running event =0, check permdf and expval if there are missing values
      pvalue <- fishertest$p.value
      pvalues[[i]] <- pvalue
    }
  }
  maxPval <- which.max(as.numeric(unlist(pvalues)))
  Pvalue <- pvalues[[maxPval]]
  print(Pvalue)
  if (all(perm_at_interval[, maxPval]==0)==TRUE) {
    return(NULL)
  } else {
    return(perm_at_interval[, maxPval,drop=FALSE])
  }
}

find_bestComb_fisher_all_intervals_new <- function(input_freqMatrix, input_expVal) {
  #filter the freqMatrix with only the matched one
  expVal_colnames <- colnames(input_expVal)
  freqMatrix <- input_freqMatrix[, expVal_colnames, drop = FALSE]
  
  bestComb_all_intervals <- NULL
  intervals <- parse_number(colnames(input_expVal))
  
  for (interval in intervals) {
    # give expval
    bestComb <- find_bestComb_fisher_at_interval(freqMatrix,input_expVal,interval)
    bestComb_all_intervals <- cbind(bestComb_all_intervals, bestComb)
  }
  colnames(bestComb_all_intervals) <- colnames(input_expVal)
  return(bestComb_all_intervals)
}

sampling_event_0_new <- function(inputdata, # dataset. The naming of columns have to follow breast cancer survival analysis
                                 testmethod = "chitest", # default statistical analysis is chi square test
                                 radius = NULL # user can define radius, otherwise go for minimal radius
) {
  df <- preprocess(inputdata) # remove NA data
  if (is.null(radius) == TRUE) {
    radius <- get_best_radius(df, 0, testmethod)
  } else {
    radius <- radius # apply user defined radius
  }
  freqMatrix <- get_freqMatrix(df, 0, radius)
  actDist <- get_actDist(df,0)
  expVal <- get_expVal(freqMatrix, actDist)
  
  if (is.null(expVal)==TRUE) {
    print("Sampling is not possible at this radius")
    return(NULL) 
  } else if (is.null(freqMatrix)==TRUE) {
    print("No valid events in this dataset")
    return(NULL)
  } else {
    if (testmethod == "chitest") {
      bestComb_all_intervals <- find_bestComb_chisq_all_intervals_new(freqMatrix, expVal)
    } else if (testmethod == "cortest") {
      bestComb_all_intervals <- find_bestComb_cor_all_intervals_new(freqMatrix, expVal)
    }
    sampled_df <-  sampling_of_bestComb(df, bestComb_all_intervals, 0,radius)
    return(sampled_df)
  }
}

sampling_event_1_new <- function(inputdata, # dataset. The naming of columns have to follow breast cancer survival analysis
                                 testmethod = "chitest", # default statistical analysis is chi square test
                                 radius = NULL # user can define radius, otherwise go for minimal radius
) {
  df <- preprocess(inputdata) # remove NA data
  
  if (is.null(radius) == TRUE) {
    radius <- get_best_radius(df, 1, testmethod)
  } else {
    radius <- radius # apply user defined radius
  }
  freqMatrix <- get_freqMatrix(df, 1, radius)
  actDist <- get_actDist(df,1)
  expVal <- get_expVal(freqMatrix, actDist)
  
  if (is.null(expVal)==TRUE) {
    print("Sampling is not possible at this radius")
    return(NULL) 
  } else if (is.null(freqMatrix)==TRUE) {
    print("No valid events in this dataset")
    return(NULL)
  } else {
    if (testmethod == "chitest") {
      bestComb_all_intervals <- find_bestComb_chisq_all_intervals_new(freqMatrix, expVal)
    } else if (testmethod == "cortest") {
      bestComb_all_intervals <- find_bestComb_cor_all_intervals_new(freqMatrix, expVal)
    }
    sampled_df <-  sampling_of_bestComb(df, bestComb_all_intervals, 1,radius)
    return(sampled_df)
  }
}

# this function returns all samples in radius x and eventstat. The input_bestcomb shall include all intervals.
sampling_of_bestComb <- function(input_df, # extract segment at interval x
                                 input_bestComb, # extract bestcomb at interval x, input is a single column
                                 eventStat =0,
                                 input_radius =1) {
  if (dim(input_bestComb)[2] >0 ) {
    
    input_df <- input_df[which(input_df$OVERALL_SURVIVAL_EVENT== eventStat),]
    #check if there is any values at this radius in the best comb
    valid_intervals_in_bestcomb <- parse_number(colnames(input_bestComb))
    #subset input_df according to valid intervals
    interval_df <- data.frame()
    count <- 0
    #extrat interval
    for (i in valid_intervals_in_bestcomb) {
      # subset out the interval
      # based on interval, select the right range of years. Batch is a list for all batches
      segment <- subset(input_df, OVERALL_SURVIVAL_TIME_YEARS < i * input_radius & OVERALL_SURVIVAL_TIME_YEARS >= (i - 1) * input_radius)
      cat("interval",i,"extracted successfully")
      cat("\n")
      count <- count + 1
      
      if (empty(segment)==TRUE) {
        next
      } else {
        sampling_indicator <- input_bestComb[,count]
        valid_batches <- sampling_indicator[which(sampling_indicator>0)]
        #check non-zero batches exist in segment, select matched batch to run the following for loop
        matched_batches <- intersect(segment$BATCH, names(valid_batches))
      }
      
      if(length(matched_batches) == 0) {
        next
      } else
      {
        for (j in 1:length(matched_batches)) {
          index <- which(names(valid_batches) == matched_batches[j])
          num_of_extra <-  valid_batches[index]
          batch_name <- matched_batches[j]
          if(length(which(segment$BATCH==batch_name))==1)
          {
            sampled <- segment[which(segment$BATCH == batch_name),]
          } else {
            sampled <- segment[sample(which(segment$BATCH == batch_name), num_of_extra),] #error
          }
          interval_df <- rbind(interval_df, sampled)
        } 
      }
    }
    if (empty(interval_df) == TRUE) {
      cat("This dataset is not drawable at radius", input_radius)
      cat("\n")
      return(interval_df)
    } else {
      return(interval_df)
    }
  }
  else {
    return(NULL)
  }
}

Sampling <- function(input_unfiltered_data, #unprocessed data
                     test="chitest", #test: chitest, fishertest, cortest
                     input_radius_0=NULL, #not fill in which run based on best radius
                     input_radius_1=NULL #not fill in which run based on best radius
) {
  sampled_event_0 <- sampling_event_0_new(input_unfiltered_data, test,input_radius_0)
  sampled_event_1 <- sampling_event_1_new(input_unfiltered_data, test,input_radius_1)
  integrated_samples <- rbind(sampled_event_0, sampled_event_1)
  if (all(integrated_samples$OVERALL_SURVIVAL_EVENT==0)==TRUE) {
    print("This dataset cannot be optimized")
    return(NULL)
  } else if (all(integrated_samples$OVERALL_SURVIVAL_EVENT==1)==TRUE) {
    print("This dataset cannot be optimized")
    return(NULL)
  } else {
    return(integrated_samples)
  }
}

KM_Plot <- function(input_df) {
  surv_object <- Surv(input_df$OVERALL_SURVIVAL_TIME_YEARS, input_df$OVERALL_SURVIVAL_EVENT)
  sf.BATCH <- surv_fit(surv_object ~ input_df$BATCH, data = input_df)
  Plot <- ggsurvplot(sf.BATCH, pval = TRUE, legend.title = "")
  return(Plot) 
}
# step 2: Bar plot for all series separately and also all samples together as we have done in the GSEA R source code
# input_df: preprocessed object
# eventStat: survival event in 0 or 1

plot <- function(input_df, eventStat = 0) {
  # here you try to calculate the scale conversion between years and rank
  scale_length <- length(rank(input_df$OVERALL_SURVIVAL_TIME_YEARS, ties.method = "random"))
  surv_year_max <- max(input_df$OVERALL_SURVIVAL_TIME_YEARS)
  x_scale <- scale_length / surv_year_max
  
  # survival status: alive, OVERALL_SURVIVAL_EVENT = 0
  if (eventStat == 0) {
    event0 <- 1 - input_df$OVERALL_SURVIVAL_EVENT
    Subsurv_0 <- ggplot(input_df, aes(x = rank(input_df$OVERALL_SURVIVAL_TIME_YEARS, ties.method = "random"), y = event0, fill = BATCH)) +
      geom_bar(stat = "identity") +
      facet_grid(input_df$BATCH ~ .)
    Subsurv_0 <- Subsurv_0 + theme(
      panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", strip.text.y = element_text(angle = 0)
    ) +
      scale_x_continuous(breaks = seq(from = 0, to = scale_length, by = x_scale), labels = c(0:surv_year_max)) +
      labs(x = "Survival Time (Years)", y = "Survival Event", title = "All Survived/censored Patients")
    print(Subsurv_0)
    return(Subsurv_0)
  } else if (eventStat == 1) {
    Subsurv_1 <- ggplot(input_df, aes(x = rank(input_df$OVERALL_SURVIVAL_TIME_YEARS, ties.method = "random"), y = input_df$OVERALL_SURVIVAL_EVENT, fill = BATCH)) +
      geom_bar(stat = "identity") +
      facet_grid(input_df$BATCH ~ .)
    Subsurv_1 <- Subsurv_1 + theme(
      panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", strip.text.y = element_text(angle = 0)
    ) +
      scale_x_continuous(breaks = seq(from = 0, to = scale_length, by = x_scale), labels = c(0:surv_year_max)) +
      labs(x = "Survival Time (Years)", y = "Survival Event", title = "All deceased Patients")
    print(Subsurv_1)
    return(Subsurv_1)
  } else if (eventStat == 2) {
    Subsurv <- ggplot(input_df, aes(x = rank(input_df$OVERALL_SURVIVAL_TIME_YEARS, ties.method = "random"), y = 1, fill = BATCH)) +
      geom_bar(stat = "identity") +
      facet_grid(input_df$BATCH ~ .)
    Subsurv <- Subsurv + theme(
      panel.grid.major = element_blank(), panel.background = element_blank(), panel.grid.minor = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", strip.text.y = element_text(angle = 0)
    ) +
      scale_x_continuous(breaks = seq(from = 0, to = scale_length, by = x_scale), labels = c(0:surv_year_max)) +
      labs(x = "Survival Time (Years)", y = "Survival Event", title = "All Patients")
    print(Subsurv)
    return(Subsurv)
  }
}
