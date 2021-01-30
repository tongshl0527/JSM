library(ggplot2)
library(tibble)
library(partitions)
library(survival)
library(survminer)
library(dplyr)
library(readr)
library(plyr)
library(combinat)
library(RcppAlgos)
library(gtable)
library(gridExtra)

main_new <- function(input_unfiltered_data, # unprocessed data
                     test = "chitest", # test: chitest, fishertest, cortest
                     input_radius_0 = NULL, # not fill in which run based on best radius
                     input_radius_1 = NULL # not fill in which run based on best radius
) {
  input_data <- preprocess(input_unfiltered_data) # remove missing values
  if (is.null(input_radius_0) == TRUE & is.null(input_radius_1) == FALSE) { # no radius in parameter
    print("Error: Please fill in radius for both event status")
    return(NULL)
  } else if (is.null(input_radius_0) == FALSE & is.null(input_radius_1) == TRUE) { # no radius in parameter
    print("Error: Please fill in radius for both event status")
    return(NULL)
  } else if (is.null(input_radius_0) & is.null(input_radius_1) == TRUE) { # radii are missing
    best_pairs <- get_best_radius_new(input_data, testmethod = test) # get the best pair of radii for each status
    input_radius_0 <- best_pairs[1] # radius for status 0
    input_radius_1 <- best_pairs[2] # radius for status 1
    radius_0 <- unlist(input_radius_0)
    radius_1 <- unlist(input_radius_1)
    print(best_pairs) #show best radii for both status
  } else {
    radius_0 <- input_radius_0 # fill radius based on input radius
    radius_1 <- input_radius_1
  }
  sampled_df <- Sampling(input_unfiltered_data, test, radius_0, radius_1) # draw samples based on the given radii
  invisible(capture.output(sampled_df)) #hide printout 
  
  plot0 <- Bar_plot(sampled_df, 0)
  plot1 <- Bar_plot(sampled_df, 1)
  print(show_Two_Plots(plot0, plot1)) 
  # print(KM_Plot(sampled_df)) #sketch kaplan meier curve
  return(sampled_df) # return drawn samples
}

preprocess <- function(inputData) {
  df <- na.omit(inputData) # remove N/A
  df <- df[order(df$OVERALL_SURVIVAL_TIME_YEARS), ] # sort survival years
  print("Success: All N/A records are removed")
  return(df)
}

get_batchList <- function(input_df # input_df should be preprocessed
) { # show all batches
  table <- count(input_df$BATCH)
  batches <- table$x
  return(batches)
}

sampling_all_radii <- function(input_data, event = 0, testmethod = "chitest") {
  pvalue_list_radii <- list() # a list to store each radius's pvalue
  comb_df <- data.frame(matrix(ncol = 5)) # create a dataframe to store all combinations and label them with an extra column
  colnames(comb_df) <- colnames(input_data)
  colnames(comb_df)[5] <- "radius"

  actDist <- get_actDist(input_data, event) # Try to optimize speed: assign two cores to run status 0 and 1 separately
  radii <- get_radiusRange(input_data, event) # get all possible radii for status 0 and 1
  count <- 0
  for (radius in radii) {
    count <- count + 1

    freqMatrix <- get_freqMatrix(input_data, event, radius)
    expVal <- get_expVal(freqMatrix, actDist)

    if (is.null(expVal) == TRUE) {
      next
    } else {
      bestComb_all_intervals <- pick_test_and_find_bestComb_all_intervals(freqMatrix, expVal, testmethod) #find the best permutation for each time segment 
      sampled_df <- sampling_of_bestComb(input_data, bestComb_all_intervals, event, radius) #sampling according to the best permutations

      if (empty(sampled_df) == TRUE) {
        comb_df <- comb_df #if nothing to sample, do nothing
      } else {
        sampled_df$radius <- radius #assign its radius to each sampled dataset
        comb_df <- rbind(comb_df, sampled_df)
      }
    }
  }
  return(comb_df)
}

sampling_event <- function(inputdata, # dataset. The naming of columns have to follow breast cancer survival analysis
                           testmethod = "chitest", # default statistical analysis is chi square test
                           radius = NULL, # user can define radius, otherwise go for minimal radius
                           event = 0 #status event 0/1
) {
  df <- preprocess(inputdata) # remove NA data
  freqMatrix <- get_freqMatrix(df, event, radius) #get frequency tables for all time segments in radius X in status event 
  actDist <- get_actDist(df, event) #get actual distribution in status event 
  expVal <- get_expVal(freqMatrix, actDist) #get expected value 
  
  if (is.null(expVal) == TRUE) { #invalid expected values
    print("Sampling is not possible at this radius")
    return(NULL)
  } else if (is.null(freqMatrix) == TRUE) { #invalid frequency matrix
    print("No valid events in this dataset")
    return(NULL)
  } else {
    bestComb_all_intervals <- pick_test_and_find_bestComb_all_intervals(freqMatrix, expVal, testmethod)
    sampled_df <- sampling_of_bestComb(df, bestComb_all_intervals, event, radius) 
    return(sampled_df)
  }
}

# sample_index_all_radii <- function(input_sampled_dataset #result from sampling_all_radii()
#                                    ) {
#   stat <- preprocess(input_sampled_dataset)
#   stat_radii <- sort(unique(stat$radius))
#   stat_list <- list()
#   count <- 0
#   for (radius in stat_radii) {
#     count <- count + 1
#     stat_list[[count]] <- which(stat$radius == radius)
#     names(stat_list)[count] <- radius
#   }
#   return(stat_list)
# }

get_best_radius_new <- function(input_data, # preprocessed dataset
                                testmethod = "chitest" # select which test. chitest: chisquare test, cortest: correlation test, fishtest: fisher test
) {
  pvalue_list_radii <- list() # a list to store each radius's pvalue
  comb_df <- data.frame(matrix(ncol = 6)) # create a dataframe to store all combinations and label them with an extra column
  colnames(comb_df) <- colnames(input_data)
  colnames(comb_df)[5] <- "Stat0_radius"
  colnames(comb_df)[6] <- "Stat1_radius"

  actDist_0 <- get_actDist(input_data, 0) # Try to optimize speed: assign two cores to run status 0 and 1 separately
  actDist_1 <- get_actDist(input_data, 1) #

  radii_0 <- get_radiusRange(input_data, 0) # get all possible radii for status 0 and 1
  radii_1 <- get_radiusRange(input_data, 1)

  sampled_0 <- sampling_all_radii(input_data, 0, testmethod) #running speed: slow
  #stat0_list <- sample_index_all_radii(sampled_0)

  sampled_1 <- sampling_all_radii(input_data, 1, testmethod)
  #stat1_list <- sample_index_all_radii(sampled_1)

  stat0_radii <- sort(unique(sampled_0$radius)) #all radii for status 0
  stat1_radii <- sort(unique(sampled_1$radius)) #all radii for status 1
  combination <- expand.grid(stat0_radii, stat1_radii) # generate combinations based on the index in dataframe

  # pvalue to be stored
  p_list <- list()
  for (i in 1:nrow(combination)) { #loop through combinations
    var1 <- combination[i, 1]
    df_1 <- sampled_0[which(sampled_0$radius == var1), ] 
    
    var2 <- combination[i, 2]
    df_2 <- sampled_1[which(sampled_1$radius == var2), ]
    
    df <- rbind(df_1, df_2)
    pvalue <- get_survPlot(df)[2]
    p_list[i] <- pvalue
    names(p_list)[i] <- toString(combination[i, ])
  }

  maxPval <- which.max(p_list)
  Pvalue <- p_list[[maxPval]]
  best_comb <- combination[maxPval, ]
  best_comb[3] <- Pvalue
  names(best_comb)[1] <- "Stat0_radius"
  names(best_comb)[2] <- "Stat1_radius"
  names(best_comb)[3] <- "Pvalue"
  return(best_comb)
} #return a list of best radii for both survival status 

compareBatch_All <- function(input_df) {
  ptable <- pairwise_survdiff(Surv(OVERALL_SURVIVAL_TIME_YEARS, OVERALL_SURVIVAL_EVENT) ~ BATCH, data = input_df)
  return(ptable)
} # return a table of p-values in comparison between batches using log rank test

show_Two_Plots <- function(input_plot_1, input_plot_2) {
  return(grid.arrange(input_plot_1, input_plot_2, ncol = 2))
} # visualize any two plots on one picture

get_survPlot <- function(input_df) {
  surv_object <- Surv(input_df$OVERALL_SURVIVAL_TIME_YEARS, input_df$OVERALL_SURVIVAL_EVENT)
  sf.BATCH <- surv_fit(surv_object ~ input_df$BATCH, data = input_df)
  return(surv_pvalue(sf.BATCH))
} # Default is "survdiff" (or "log-rank") return p-value of the kaplan meiser curve

get_radiusRange <- function(input_df, # preprocesssed dataset
                            event = 0 # which status: 0 or 1
) {
  input_df <- subset(input_df, input_df$OVERALL_SURVIVAL_EVENT == event) # subset for a status
  distance_matrix <- as.matrix(dist(input_df$OVERALL_SURVIVAL_TIME_YEARS)) # returns the distance matrix computed by using the specified distance measure to compute the distances between survival time.
  sorted_disMatrix <- apply(distance_matrix, 1, sort)
  radius <- apply(sorted_disMatrix, 1, min)
  radius <- unique(radius)[2:length(radius)]
  radiusRange <- radius[!is.na(radius)] # from low radius to high radius. Lowest: the shortest time segment. Highest: the longest time segment.
  return(radiusRange)
} # get all possible radii of a status

get_freqMatrix <- function(input_df, # preprocessed dataset
                           eventStat = 0, # status event: 0/1
                           input_radius # a radius
) {
  dataframe <- subset(input_df, OVERALL_SURVIVAL_EVENT == eventStat) # subset for a status
  batchList <- get_batchList(input_df) # collect names of all batches
  size_segment <- round(1 + max(dataframe$OVERALL_SURVIVAL_TIME_YEARS) / input_radius) # get
  size_segment <- unlist(size_segment)
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
} # return a table of frequency of batches for radius x at each time segment

get_actDist <- function(input_df, eventStat = 0) {
  # count total samples
  freq_total <- count(input_df, c("BATCH", "OVERALL_SURVIVAL_EVENT"))
  # count the total samples in event 0
  freq_event <- freq_total[freq_total$OVERALL_SURVIVAL_EVENT == eventStat, ]
  sum_event <- sum(freq_event$freq)
  freq_event <- transform(freq_event, proportion = freq / sum_event)
  # standard ratio of each batch is stored in a list
  return(freq_event)
} # actual distrubution: proportion of each batch = frequency / sum of all batches' frequencies

get_expval_at_interval <- function(input_expVal, input_at_Interval = 1) {
  valid_intervals_in_expVal <- parse_number(colnames(input_expVal))
  idx <- match(paste("interval", input_at_Interval, sep = "_"), colnames(input_expVal))
  expVal <- input_expVal[, idx, drop = FALSE]
  return(expVal)
} # select a time segment to show its expected values of all batches

get_expVal <- function(input_freqMatrix, input_actDist) {
  # remove columns that are zeros only
  input_freqMatrix <- input_freqMatrix[, colSums(input_freqMatrix != 0) > 0, drop = FALSE]
  # total number of event 0 in this segment
  n_col <- ncol(input_freqMatrix)
  expVal_matrix <- matrix(0, nrow(input_freqMatrix), ncol(input_freqMatrix))
  rownames(expVal_matrix) <- rownames(input_freqMatrix)
  colnames(expVal_matrix) <- colnames(input_freqMatrix)
  for (i in 1:n_col) {
    sumCol <- sum(input_freqMatrix[, i])
    for (j in input_actDist$BATCH) {
      expVal <- round(input_actDist[input_actDist$BATCH == j, ]$proportion * sumCol)
      expVal_matrix[j, i] <- expVal
    }
  }
  # remove columns with zeros
  expVal_matrix <- expVal_matrix[, colSums(expVal_matrix != 0) > 0, drop = FALSE]
  if (dim(expVal_matrix)[2] > 0) {
    return(expVal_matrix)
  } else {
    return(NULL)
  }
} #

get_permutation <- function(input_freqMatrix, input_radiusInterval = 1) {
  # convert input_radiusInterval to colname
  idx <- match(paste("interval", input_radiusInterval, sep = "_"), colnames(input_freqMatrix))
  column <- input_freqMatrix[, idx, drop = FALSE]

  a1 <- sapply(column, function(x) {
    c(0:x)
  })
  a2 <- expand.grid(a1)

  if (sum(a2) == 0) {
    a2 <- t(a2)
  }
  colnames(a2) <- rownames(input_freqMatrix)
  rownames(a2) <- paste("interval", input_radiusInterval, c(1:nrow(a2)), sep = "_")
  return(t(a2))
} # get permutation for one interval!!

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

find_bestComb_chisq_at_interval <- function(input_freqMatrix, # remove columns with 0
                                            input_expVal, # remove columns with 0
                                            input_at_Interval = 1 # apply intervals that only exist in expval!
) {
  # convert input_radiusInterval to colname
  idx <- match(paste("interval", input_at_Interval, sep = "_"), colnames(input_freqMatrix))
  # select columns match to expval
  freqMatrix <- input_freqMatrix[, idx, drop = FALSE]
  # interval_freqMatrix<- input_freqMatrix[,input_at_Interval,drop=FALSE]
  perm_at_interval <- get_permutation(freqMatrix, input_at_Interval)
  expval_at_interval <- get_expval_at_interval(input_expVal, input_at_Interval)
  pvalues <- list()
  for (i in 1:ncol(perm_at_interval)) {
    if (length(unique(perm_at_interval[, i])) == 1) {
      pvalues[[i]] <- 0
      names(pvalues)[i] <- colnames(perm_at_interval[, i, drop = FALSE])
    } else if (length(unique(expval_at_interval)) == 1) {
      pvalues[[i]] <- 0
      names(pvalues)[i] <- colnames(perm_at_interval[, i, drop = FALSE])
    } else {
      chitest <- chisq.test(perm_at_interval[, i], expval_at_interval) # error when running event =0, check permdf and expval if there are missing values
      pvalue <- chitest$p.value
      pvalues[[i]] <- pvalue
      names(pvalues)[i] <- colnames(perm_at_interval[, i, drop = FALSE])
    }
  }
  if (any(pvalues > 0) == TRUE) {
    maxPval <- which.max(as.numeric(unlist(pvalues)))
    Pvalue <- pvalues[[maxPval]]
    Perm_interval <- names(pvalues)[maxPval]
    cat(Perm_interval, "has a P-value:", Pvalue)
    cat("\n")
    if (all(perm_at_interval[, maxPval] == 0) == TRUE) {
      return(NULL)
    } else {
      return(perm_at_interval[, maxPval, drop = FALSE])
    }
  } else {
    return(NULL)
  }
}

find_bestComb_chisq_all_intervals_new <- function(input_freqMatrix, input_expVal) {
  # filter the freqMatrix with only the matched one
  expVal_colnames <- colnames(input_expVal)
  freqMatrix <- input_freqMatrix[, expVal_colnames, drop = FALSE]

  bestComb_all_intervals <- NULL
  intervals <- parse_number(colnames(input_expVal))

  for (interval in intervals) {
    # at_interval <- parse_number(colnames(input_expVal[,interval,drop=FALSE]))
    # give expval
    bestComb <- find_bestComb_chisq_at_interval(freqMatrix, input_expVal, interval)
    # bestComb_all_intervals[, interval] <- bestComb
    bestComb_all_intervals <- cbind(bestComb_all_intervals, bestComb)
  }
  colnames(bestComb_all_intervals) <- colnames(input_expVal)
  return(bestComb_all_intervals)
}


find_bestComb_cor_at_interval <- function(input_freqMatrix, # remove columns with 0
                                          input_expVal, # remove columns with 0
                                          input_at_Interval = 1 # apply intervals that only exist in expval!
) {
  # convert input_radiusInterval to colname
  idx <- match(paste("interval", input_at_Interval, sep = "_"), colnames(input_freqMatrix))
  # select columns match to expval
  freqMatrix <- input_freqMatrix[, idx, drop = FALSE]
  # interval_freqMatrix<- input_freqMatrix[,input_at_Interval,drop=FALSE]
  perm_at_interval <- get_permutation(freqMatrix, input_at_Interval)
  expval_at_interval <- get_expval_at_interval(input_expVal, input_at_Interval)
  pvalues <- list()
  for (i in 1:ncol(perm_at_interval)) {
    if (length(unique(perm_at_interval[, i])) == 1) {
      pvalues[[i]] <- 0
    } else if (length(unique(expval_at_interval)) == 1) {
      pvalues[[i]] <- 0
    } else {
      cortest <- cor.test(perm_at_interval[, i], expval_at_interval) # error when running event =0, check permdf and expval if there are missing values
      pvalue <- cortest$p.value
      pvalues[[i]] <- pvalue
    }
  }
  maxPval <- which.max(as.numeric(unlist(pvalues)))
  Pvalue <- pvalues[[maxPval]]
  print(Pvalue)
  if (all(perm_at_interval[, maxPval] == 0) == TRUE) {
    return(NULL)
  } else {
    return(perm_at_interval[, maxPval, drop = FALSE])
  }
}

find_bestComb_cor_all_intervals_new <- function(input_freqMatrix, input_expVal) {
  # filter the freqMatrix with only the matched one
  expVal_colnames <- colnames(input_expVal)
  freqMatrix <- input_freqMatrix[, expVal_colnames, drop = FALSE]

  bestComb_all_intervals <- NULL
  intervals <- parse_number(colnames(input_expVal))

  for (interval in intervals) {
    # give expval
    bestComb <- find_bestComb_cor_at_interval(freqMatrix, input_expVal, interval)
    bestComb_all_intervals <- cbind(bestComb_all_intervals, bestComb)
  }
  colnames(bestComb_all_intervals) <- colnames(input_expVal)
  return(bestComb_all_intervals)
}

find_bestComb_fisher_at_interval <- function(input_freqMatrix, # remove columns with 0
                                             input_expVal, # remove columns with 0
                                             input_at_Interval = 1 # apply intervals that only exist in expval!
) {
  # convert input_radiusInterval to colname
  idx <- match(paste("interval", input_at_Interval, sep = "_"), colnames(input_freqMatrix))
  # select columns match to expval
  freqMatrix <- input_freqMatrix[, idx, drop = FALSE]
  perm_at_interval <- get_permutation(freqMatrix, input_at_Interval) #get all permutations of the frequency table at the time segment input_at_Interval
  expval_at_interval <- get_expval_at_interval(input_expVal, input_at_Interval) # get the expected values of the frequency table at the same time segment
  pvalues <- list()
  #compare each permutation with the expacted value
  for (i in 1:ncol(perm_at_interval)) {
    if (length(unique(perm_at_interval[, i])) == 1) { #invalid permutation, assign the pvalue 0
      pvalues[[i]] <- 0
    } else if (length(unique(expval_at_interval)) == 1) { #invalid expected values, assign the pvalue 0
      pvalues[[i]] <- 0
    } else {
      fishertest <- fisher.test(perm_at_interval[, i], expval_at_interval) # perform fisher test 
      pvalue <- fishertest$p.value 
      pvalues[[i]] <- pvalue
    }
  }
  maxPval <- which.max(as.numeric(unlist(pvalues))) #find the index of the highest pvalue
  Pvalue <- pvalues[[maxPval]] #extract that pvalue
  print(Pvalue) #show that pvalue
  if (all(perm_at_interval[, maxPval] == 0) == TRUE) { #no batches found in the permutation
    return(NULL)
  } else {
    return(perm_at_interval[, maxPval, drop = FALSE]) #extract the permutation
  }
}

find_bestComb_fisher_all_intervals_new <- function(input_freqMatrix, #result from get_freqMatrix() 
                                                   input_expVal #result from get_expVal()
                                                   ) {
  expVal_colnames <- colnames(input_expVal)
  freqMatrix <- input_freqMatrix[, expVal_colnames, drop = FALSE] # reduce dimension: filter the freqMatrix to show batches appeared in the expected value

  bestComb_all_intervals <- NULL
  intervals <- parse_number(colnames(input_expVal)) #extract time segments from the expected values

  for (interval in intervals) { #loop through all time segments in the freqMatrix
    bestComb <- find_bestComb_fisher_at_interval(freqMatrix, input_expVal, interval) #extract the permutation with the highest p-value at time segment interval
    bestComb_all_intervals <- cbind(bestComb_all_intervals, bestComb) #save it in a matrix
  }
  colnames(bestComb_all_intervals) <- colnames(input_expVal)
  return(bestComb_all_intervals) # return a matrix of permutations(highest pvalue) in each time segment 
}

# sampling_event <- function(inputdata, # dataset. The naming of columns have to follow breast cancer survival analysis
#                            testmethod = "chitest", # default statistical analysis is chi square test
#                            radius = NULL, # user can define radius, otherwise go for minimal radius
#                            event = 0 #status event 0/1
# ) {
#   df <- preprocess(inputdata) # remove NA data
#   freqMatrix <- get_freqMatrix(df, event, radius) #get frequency tables for all time segments in radius X in status event 
#   actDist <- get_actDist(df, event) #get actual distribution in status event 
#   expVal <- get_expVal(freqMatrix, actDist) #get expected value 
# 
#   if (is.null(expVal) == TRUE) { #invalid expected values
#     print("Sampling is not possible at this radius")
#     return(NULL)
#   } else if (is.null(freqMatrix) == TRUE) { #invalid frequency matrix
#     print("No valid events in this dataset")
#     return(NULL)
#   } else {
#     if (testmethod == "chitest") { #apply chi square test to find permutations for all time segments
#       bestComb_all_intervals <- find_bestComb_chisq_all_intervals_new(freqMatrix, expVal)
#     } else if (testmethod == "cortest") { #apply correlation test to find permutations for all time segments
#       bestComb_all_intervals <- find_bestComb_cor_all_intervals_new(freqMatrix, expVal)
#     }
#     sampled_df <- sampling_of_bestComb(df, bestComb_all_intervals, event, radius) 
#     return(sampled_df)
#   }
# }

# this function returns all samples in radius x and eventstat. The input_bestcomb shall include all intervals.
sampling_of_bestComb <- function(input_df, # extract segment at interval x
                                 input_bestComb, # extract bestcomb at interval x, input is a single column
                                 eventStat = 0,
                                 input_radius = 1) {
  if (dim(input_bestComb)[2] > 0) {
    input_df <- input_df[which(input_df$OVERALL_SURVIVAL_EVENT == eventStat), ]
    # check if there is any values at this radius in the best comb
    valid_intervals_in_bestcomb <- parse_number(colnames(input_bestComb))
    # subset input_df according to valid intervals
    interval_df <- data.frame()
    count <- 0
    # extrat interval
    for (i in valid_intervals_in_bestcomb) {
      segment <- subset(input_df, OVERALL_SURVIVAL_TIME_YEARS < i * input_radius & OVERALL_SURVIVAL_TIME_YEARS >= (i - 1) * input_radius) # subset out the interval based on interval, select the right range of years. Batch is a list for all batches
      cat("interval", i, "extracted successfully")
      cat("\n")
      count <- count + 1

      if (empty(segment) == TRUE) {
        next
      } else {
        sampling_indicator <- input_bestComb[, count] #
        valid_batches <- sampling_indicator[which(sampling_indicator > 0)] #extract index of best permutation that are not zero
        matched_batches <- intersect(segment$BATCH, names(valid_batches)) # check non-zero batches exist in segment, select matched batch to run the following for loop
      }

      if (length(matched_batches) == 0) { #no match
        next
      } else {
        for (j in 1:length(matched_batches)) { 
          index <- which(names(valid_batches) == matched_batches[j]) 
          num_of_extra <- valid_batches[index] 
          batch_name <- matched_batches[j] #extract valid batches
          if (length(which(segment$BATCH == batch_name)) == 1) { #one matched batch
            sampled <- segment[which(segment$BATCH == batch_name), ]
          } else {
            sampled <- segment[sample(which(segment$BATCH == batch_name), num_of_extra), ]  #more than one matched batches
          }
          interval_df <- rbind(interval_df, sampled) #extract samples
        }
      }
    }
    if (empty(interval_df) == TRUE) {
      cat("This dataset is not drawable at radius", input_radius)
      cat("\n")
      #return(interval_df)
      return(NULL)
    } else {
      return(interval_df)
    }
  }
  else {
    cat("This dataset is not drawable at radius", input_radius)
    cat("\n")
    return(NULL)
  }
}

Sampling <- function(input_unfiltered_data, # unprocessed data
                     test = "chitest", # test: chitest, fishertest, cortest
                     input_radius_0 = NULL, # not fill in which run based on best radius
                     input_radius_1 = NULL # not fill in which run based on best radius
) {
  sampled_event_0 <- sampling_event(input_unfiltered_data, test, input_radius_0, 0)
  sampled_event_1 <- sampling_event(input_unfiltered_data, test, input_radius_1, 1)
  integrated_samples <- rbind(sampled_event_0, sampled_event_1)
  if (all(integrated_samples$OVERALL_SURVIVAL_EVENT == 0) == TRUE) {
    print("This dataset cannot be optimized")
    return(NULL)
  } else if (all(integrated_samples$OVERALL_SURVIVAL_EVENT == 1) == TRUE) {
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

Bar_plot <- function(input_df, eventStat = 0) {
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
