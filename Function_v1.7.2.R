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
library(parallel)
library(tm)


# main_wo_plot <- function(input_unfiltered_data, # unprocessed data
#                          input_radius_0 = NULL, # not fill in which run based on best radius
#                          input_radius_1 = NULL # not fill in which run based on best radius
# ) {
#   input_data <- preprocess(input_unfiltered_data) # remove missing values
#   if (is.null(input_radius_0) == TRUE & is.null(input_radius_1) == FALSE) { # no radius in parameter
#     print("Error: Please fill in radius for both event status")
#     return(NULL)
#   } else if (is.null(input_radius_0) == FALSE & is.null(input_radius_1) == TRUE) { # no radius in parameter
#     print("Error: Please fill in radius for both event status")
#     return(NULL)
#   } else if (is.null(input_radius_0) & is.null(input_radius_1) == TRUE) { # radii are missing
#     best_pairs <- get_best_radius(input_data) # get the best pair of radii for each status
#     input_radius_0 <- best_pairs[1] # radius for status 0
#     input_radius_1 <- best_pairs[2] # radius for status 1
#     radius_0 <- unlist(input_radius_0)
#     radius_1 <- unlist(input_radius_1)
#     print(best_pairs) #show best radii for both status
#   } else {
#     radius_0 <- input_radius_0 # fill radius based on input radius
#     radius_1 <- input_radius_1
#   }
#   if (is.null(radius_0)|is.null(radius_1)) {
#     print("Radii = NULL. Cannot sample the dataset")
#     return(NULL)
#   } else {
#     sampled_df <- Sampling(input_unfiltered_data, radius_0, radius_1) # draw samples based on the given radii
#     invisible(capture.output(sampled_df)) #hide printout 
#     return(sampled_df) # return drawn samples
#   }
# }

# main_all_in_one <- function(input_unfiltered_data, # unprocessed data
#                          input_radius_0 = NULL, # not fill in which run based on best radius
#                          input_radius_1 = NULL # not fill in which run based on best radius
# ) {
#   all_in_one <- list()
#   input_data <- preprocess(input_unfiltered_data) # remove missing values
#   if (is.null(input_radius_0) == TRUE & is.null(input_radius_1) == FALSE) { # no radius in parameter
#     print("Error: Please fill in radius for both event status")
#     return(NULL)
#   } else if (is.null(input_radius_0) == FALSE & is.null(input_radius_1) == TRUE) { # no radius in parameter
#     print("Error: Please fill in radius for both event status")
#     return(NULL)
#   } else if (is.null(input_radius_0) & is.null(input_radius_1) == TRUE) { # radii are missing
#     best_pairs <- get_best_radius(input_data) # get the best pair of radii for each status
#     input_radius_0 <- best_pairs[1] # radius for status 0
#     input_radius_1 <- best_pairs[2] # radius for status 1
#     radius_0 <- unlist(input_radius_0)
#     radius_1 <- unlist(input_radius_1)
#     print(best_pairs) #show best radii for both status
#   } else {
#     radius_0 <- input_radius_0 # fill radius based on input radius
#     radius_1 <- input_radius_1
#   }
#   if (is.null(radius_0)|is.null(radius_1)) {
#     print("Radii = NULL. Cannot sample the dataset")
#     return(NULL)
#   } else {
#     sampled_df <- Sampling(input_unfiltered_data, radius_0, radius_1) # draw samples based on the given radii
#     invisible(capture.output(sampled_df)) #hide printout 
#     all_in_one[[1]] <- sampled_df
#     all_in_one[[2]] <- radius_0
#     all_in_one[[3]] <- radius_1
#     all_in_one[[4]] <- get_survPlot(sampled_df)[2]
#     all_in_one[[5]] <- KM_Plot(sampled_df)
#     all_in_one[[6]] <- Bar_plot(sampled_df,0)
#     all_in_one[[7]] <- Bar_plot(sampled_df,1)
#     names(all_in_one) <- c("Sampled Dataset", "Radius for event 0", "Radius for event 1", "P-value Log Rank Test", "KM Plot", "Bar plot event 0", "Bar plot event 1")
#     return(all_in_one) # return drawn samples
#   }
# }

new_main_all_in_one <- function(input_unfiltered_data, # unprocessed data
                            input_radius_0 = NULL, # not fill in which run based on best radius
                            input_radius_1 = NULL # not fill in which run based on best radius
) {
  all_in_one <- list()
  input_data <- preprocess(input_unfiltered_data) # remove missing values
  if (is.null(input_radius_0) == TRUE & is.null(input_radius_1) == FALSE) { # no radius in parameter
    print("Error: Please fill in radius for both event status")
    return(NULL)
  } else if (is.null(input_radius_0) == FALSE & is.null(input_radius_1) == TRUE) { # no radius in parameter
    print("Error: Please fill in radius for both event status")
    return(NULL)
  } else if (is.null(input_radius_0) & is.null(input_radius_1) == TRUE) { # radii are missing
    sampled_df <- get_best_radius(input_data) # get the best pair of radii for each status
    radius_0 <- unique(sampled_df$radius)[1] # fill radius based on input radius
    radius_1 <- unique(sampled_df$radius)[2]
    sampled_df <- sampled_df[,!(names(sampled_df) %in% "radius")]
      #remove radius column
  } else {
    radius_0 <- input_radius_0 # fill radius based on input radius
    radius_1 <- input_radius_1
    sampled_df <- Sampling(input_data, radius_0, radius_1)
  }
    #invisible(capture.output(sampled_df)) #hide printout 
    all_in_one[[1]] <- sampled_df
    all_in_one[[2]] <- radius_0
    all_in_one[[3]] <- radius_1
    all_in_one[[4]] <- get_survPlot(sampled_df)[2]
    all_in_one[[5]] <- KM_Plot(sampled_df)
    all_in_one[[6]] <- Bar_plot(sampled_df,0)
    all_in_one[[7]] <- Bar_plot(sampled_df,1)
    names(all_in_one) <- c("Sampled Dataset", "Radius for event 0", "Radius for event 1", "P-value Log Rank Test", "KM Plot", "Bar plot event 0", "Bar plot event 1")
    print(all_in_one[[5]])
    return(all_in_one) # return drawn samples
  }
#1 generate three datasets per size
#2 run benchmark three times on the reduced input_data

benchmark <- function(input_data) {
  # 0. sampling data
  timer <- system.time(batch_effect_free_dataset <- new_main_all_in_one(input_data))
  sampled_data <- batch_effect_free_dataset$`Sampled Dataset`
  # 1. To compare number of samples having event 0 and 1 distribution (before vs after correction)
  nrow_0_bef <- count(input_data$OVERALL_SURVIVAL_EVENT)$freq[1]
  nrow_1_bef <- count(input_data$OVERALL_SURVIVAL_EVENT)$freq[2]
  nrow_0_aft <- count(sampled_data$OVERALL_SURVIVAL_EVENT)$freq[1]
  nrow_1_aft <- count(sampled_data$OVERALL_SURVIVAL_EVENT)$freq[2]
  # 2. Total number of samples 
  total_nrow_bef <- nrow(input_data)
  total_nrow_aft <- nrow(sampled_data)
  # 3. KM curves with log rank test p values
  pval_bef <- get_survPlot(input_data)
  pval_bef <- round(pval_bef$pval, 3)
  pval_aft <- round(as.numeric(batch_effect_free_dataset$`P-value Log Rank Test`),3)
  # 4. Store radius for individual event (0/1)
  radius_0 <- batch_effect_free_dataset$`Radius for event 0`
  radius_1 <- batch_effect_free_dataset$`Radius for event 1`
  # 5. generate dataframe
  x <- data.frame("Total Sample Size before" = total_nrow_bef, "Total Sample Size after" = total_nrow_aft,"Sample Size event 0 before" = nrow_0_bef, "Sample Size event 0 after" = nrow_0_aft, "Sample Size event 1 before" = nrow_1_bef, "Sample Size event 1 after" = nrow_1_aft, "Pval before" = pval_bef, "Pval after" = pval_aft, "Event 0 radius"= radius_0, "Event 1 radius"= radius_1, "Time"= timer[1])
  rownames(x) <- NULL
  # 6. generate a list to contain the dataframe
  result <- list(x, input_data, sampled_data)
  names(result) <- c("summary","bef.corr","aft.corr")
  return(result) #list[[1]]: summary sheet; list[[2]]: before correction survival data; list[[3]]: after correction survival data
}

main_benchmark <- function(input_data, list, input_rep_times =3) {
  #input the number of regenerating datasets
  data_in_all_sizes <- replicate(input_rep_times, reduce_dataset_by_percentage(input_data, list), simplify=FALSE)
  summary<- list()
  bef.corr <- list()
  aft.corr <- list()
  count <- 0
  for (i in 1:length(list)) {
    for (j in 1:length(data_in_all_sizes)) {
      count <- count +1
      reduced_data <- data_in_all_sizes[[j]][[i]]
      benchmarked <- benchmark(reduced_data)
      sum <- as.data.frame(benchmarked[1])
      rownames(sum) <- paste(list[i], j, sep = "_")
      summary[[count]] <- sum
      bef.corr[[count]] <- as.data.frame(benchmarked[2])
      names(bef.corr)[[count]] <-paste(list[i], j, sep = "_")
      aft.corr[[count]] <- as.data.frame(benchmarked[3])
      names(aft.corr)[[count]] <-paste(list[i], j, sep = "_")
    }
  }
  benchmark_summary <- do.call(rbind,summary)
  benchmark_result <- list(benchmark_summary, bef.corr, aft.corr)
  return(benchmark_result)
}

extract_survivalData_from_main_benchmark <- function(input_main_benchmark, 
                                                     input_benchmark_ID, #ID is the rownames of summary sheet (main_benchmark[[1]])
                                                     extract_aft.corr = TRUE #False will return the uncorrected survival data
                                                     ) {
  bef_corr <- input_main_benchmark[[2]]
  aft_corr <- input_main_benchmark[[3]]
  
  if (extract_aft.corr) {
    result <- as.data.frame(aft_corr[[input_benchmark_ID]])
    colnames(result) <- removeWords(colnames(result),"aft.corr.")
  } else {
    result <- as.data.frame(bef_corr[[input_benchmark_ID]])
    colnames(result) <- removeWords(colnames(result),"bef.corr.")
  }
  return(result)
}

reduce_dataset_by_percentage <- function(input_data, input_percentage_list) {
  df <- lapply(input_percentage_list, function(x) input_data[sample(nrow(input_data), size = as.integer(nrow(input_data)*x), replace = FALSE),])
  return(df)
}

#repeat main function multiple times to generate datasets
#input_data: the reduced one per size
repeat_main <- function(input_data, input_rep_times) {
  # replicate n-times 
  a <- replicate(input_rep_times, benchmark(input_data), simplify=FALSE)
  b <- do.call("rbind",a)
  return(b)
}


#give original dataset to input_data. Define the number to run sampling per percentage. Defined the size of new dataset reduced to what percentage:Which_in_percentage_list 
# repeat_sampling_dataset_at_percentage <- function(input_data, list, input_rep_times, which_in_percentage_list) {
#   # replicate n-times 
#   a <- replicate(input_rep_times, reduce_dataset_by_percentage(input_data, list), simplify=FALSE)
#   b <- lapply(1:input_rep_times, function(x) a[[x]][[which_in_percentage_list]])
#   return(b)
# }

#convert control_data to the same format as input_data
format_inputdata <- function(input_data, control_data) {
  #preprocess data
  colnames(input_data)[3:4] <- c("OVERALL_SURVIVAL_TIME_YEARS","OVERALL_SURVIVAL_EVENT")
  colnames(input_data)[1] <- "GSM_NUMBER"
  input_data$OVERALL_SURVIVAL_TIME_YEARS <- round(input_data$OVERALL_SURVIVAL_TIME_YEARS/12,2)
  #assign batch ID to the dataset 
  input_data <- merge(input_data[,-2], control_data[,c("GSM_NUMBER","BATCH")], by = "GSM_NUMBER", sort = F, all.x = T)
  return(input_data)
  }

preprocess <- function(inputData) {
  df <- na.omit(inputData) # remove N/A
  df <- df[order(df$OVERALL_SURVIVAL_TIME_YEARS), ] # sort survival years
  print("Success: All N/A records are removed")
  return(df)
}

get_batchList <- function(input_df) { # show all batches
  table <- count(input_df$BATCH)
  batches <- table$x
  return(batches)
}

show_Two_Plots <- function(input_plot_1, input_plot_2) {
  return(grid.arrange(input_plot_1, input_plot_2, ncol = 2))
} # visualize any two plots on one picture

get_survPlot <- function(input_df) {
  surv_object <- Surv(input_df$OVERALL_SURVIVAL_TIME_YEARS, input_df$OVERALL_SURVIVAL_EVENT)
  sf.BATCH <- surv_fit(surv_object ~ input_df$BATCH, data = input_df)
  return(surv_pvalue(sf.BATCH))
} # Default is "survdiff" (or "log-rank") return p-value of the kaplan meiser curve

#used
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

#used
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
  #
  cat("The current radius is at ", input_radius)
  cat("\n")
  return(freq_Matrix)
} # return a table of frequency of batches for radius x at each time segment

#used
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

#used
get_expval_at_interval <- function(input_expVal, input_at_Interval = 1) {
  valid_intervals_in_expVal <- parse_number(colnames(input_expVal))
  idx <- match(paste("interval", input_at_Interval, sep = "_"), colnames(input_expVal))
  expVal <- input_expVal[, idx, drop = FALSE]
  return(expVal)
} # select a time segment to show its expected values of all batches

#used
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
} 

#used
get_permutation <- function(input_freqMatrix, input_expVal, input_radiusInterval = 1, threshold = 3) {
  # convert input_radiusInterval to colname
  idx <- match(paste("interval", input_radiusInterval, sep = "_"), colnames(input_freqMatrix))
  column <- input_freqMatrix[, idx, drop = FALSE]
  #subsitute 
  expVal_at_interval <- get_expval_at_interval(input_expVal,input_radiusInterval)
  #store permutation of each sample in list
  a1 <- sapply(column, function(x) {
    c(0:x)
  })
  #pick n nearest values to the expected value 
  for (i in 1:nrow(input_expVal)) {
    list <- a1[[i]]
    exp <- expVal_at_interval[i]
    #if exp is found in list
    if (is.element(exp, list)) {
      a1[[i]] <- list[which(list<= exp & list>exp-threshold)]
    } else {
      #else take the last ten values from the list
      a1[[i]] <- tail(list, n=threshold)
    }
  }
  #create permutations 
  a2 <- expand.grid(a1)
  colnames(a2) <- rownames(input_freqMatrix)
  rownames(a2) <- paste("interval", input_radiusInterval, c(1:nrow(a2)), sep = "_")
  return(t(a2))
} # get permutation for one interval!!

#used
# this function returns all samples in radius x and eventstat. The input_bestcomb shall include all intervals.
sampling_of_bestComb_new <- function(input_df, # extract segment at interval x
                                 input_bestComb, # extract bestcomb at interval x, input is a single column
                                 eventStat = 0,
                                 input_radius = 1) {
if (dim(input_bestComb)[2] > 0) {
    interval_df <- data.frame()
    # check if there is any values at this radius in the best comb
    valid_intervals_in_bestcomb <- parse_number(colnames(input_bestComb))
    #loop through bestcomb intervals
    for (i in colnames(input_bestComb)) {
      if (sum(input_bestComb[i]) == 0) {
        next
      } else {
        interval_bestComb <- setNames(input_bestComb[which(input_bestComb[,i]>0),i],rownames(input_bestComb)[which(input_bestComb[,i]>0)])
        # subset input_df according to valid intervals
        valid_intervals_in_bestcomb <- parse_number(i)
        # subset out the interval based on interval, select the right range of years. Batch is a list for all batches
        segment <- subset(input_df[which(input_df$OVERALL_SURVIVAL_EVENT == eventStat), ], OVERALL_SURVIVAL_TIME_YEARS < valid_intervals_in_bestcomb * input_radius & OVERALL_SURVIVAL_TIME_YEARS >= (valid_intervals_in_bestcomb - 1) * input_radius) 
        #draw out rows based on interval_bestComb, using apply
        for (j in 1:length(interval_bestComb)) {
          drawed_segment <- segment[sample(which(segment$BATCH == names(interval_bestComb)[j]),size=interval_bestComb[j],replace=FALSE),]
          interval_df <- rbind(interval_df,drawed_segment)
        }
      }
    }
    return(interval_df)
} else {
    return(NULL)
  }
}


find_nearest_to_expected <- function(input_freqMatrix, input_expVal, input_radiusInterval = 1) {
  # convert input_radiusInterval to colname
  idx <- match(paste("interval", input_radiusInterval, sep = "_"), colnames(input_freqMatrix))
  column <- input_freqMatrix[, idx, drop = FALSE]
  #subsitute 
  expVal_at_interval <- get_expval_at_interval(input_expVal,input_radiusInterval)
  #store permutation of each sample in list
  a1 <- lapply(column, function(x) {
    c(0:x)
  })
  #pick n nearest values to the expected value 
  for (i in 1:nrow(input_expVal)) {
    list <- a1[[i]]
    exp <- expVal_at_interval[i]
    #if exp is found in list
    if (is.element(exp, list)) {
      a1[[i]] <- exp
    } else {
      #else take the last ten values from the list
      a1[[i]] <- tail(list, n=1)
    }
  }
  a2 <- t(as.data.frame(a1))
  #create permutations 
  rownames(a2) <- rownames(expVal_at_interval)
  colnames(a2) <- colnames(expVal_at_interval)
  return(a2)
}


#used
sampling_event <- function(inputdata, # dataset. The naming of columns have to follow breast cancer survival analysis
                           #testmethod = "chitest", # default statistical analysis is chi square test
                           radius = NULL, # user can define radius, otherwise go for minimal radius
                           event = 0 #status event 0/1
) {
  df <- preprocess(inputdata) # remove NA data
  freqMatrix <- get_freqMatrix(df, eventStat = event, input_radius = radius) #get frequency tables for all time segments in radius X in status event 
  actDist <- get_actDist(df, event) #get actual distribution in status event 
  expVal <- get_expVal(freqMatrix, actDist) #get expected value 
  
  if (is.null(expVal) == TRUE) { #invalid expected values
    print("Sampling is not possible at this radius")
    return(NULL)
  } else if (is.null(freqMatrix) == TRUE) { #invalid frequency matrix
    print("No valid events in this dataset")
    return(NULL)
  } else {
    #extract expVal 
    list <- parse_number(colnames(expVal))
    matrix <- sapply(list, find_nearest_to_expected, input_freqMatrix=freqMatrix, input_expVal=expVal)
    colnames(matrix) <- colnames(expVal)
    rownames(matrix) <- rownames(expVal)
    matrix <- as.data.frame(matrix)
    sampled_df <- sampling_of_bestComb_new(df, matrix, event, radius) 
    sampled_df <- unique(sampled_df)
    return(sampled_df)
  }
}

#used
Sampling <- function(input_unfiltered_data, # unprocessed data
                     input_radius_0 = NULL, # not fill in which run based on best radius
                     input_radius_1 = NULL # not fill in which run based on best radius
) {
  sampled_event_0 <- sampling_event(input_unfiltered_data, input_radius_0, 0)
  sampled_event_1 <- sampling_event(input_unfiltered_data, input_radius_1, 1)
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

#used
KM_Plot <- function(input_df) {
  surv_object <- Surv(input_df$OVERALL_SURVIVAL_TIME_YEARS, input_df$OVERALL_SURVIVAL_EVENT)
  sf.BATCH <- surv_fit(surv_object ~ input_df$BATCH, data = input_df)
  Plot <- ggsurvplot(sf.BATCH, pval = TRUE, pval.method = TRUE, legend.title = "")
  return(Plot)
}
# step 2: Bar plot for all series separately and also all samples together as we have done in the GSEA R source code
# input_df: preprocessed object
# eventStat: survival event in 0 or 1

#used
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

#used
# get_best_radius <- function(input_data) {
#   sampled_0 <- sampling_all_radii(input_data, 0) #running speed: slow
#   sampled_1 <- sampling_all_radii(input_data, 1)
# 
#   stat0_radii <- sort(unique(sampled_0$radius)) #all radii for status 0
#   stat1_radii <- sort(unique(sampled_1$radius)) #all radii for status 1
#   combination <- expand.grid(stat0_radii, stat1_radii) # generate combinations based on the index in dataframe
#   if (nrow(combination)==0) {
#     print("No radius generated. Please try a bigger dataset.")
#     return(NULL)
#   } else {
#     # pvalue to be stored
#     p_list <- sapply(1:nrow(combination), test, input_comb=combination, input_sampled_0=sampled_0, input_sampled_1=sampled_1)
#     maxPval <- which.max(p_list)
#     Pvalue <- p_list[maxPval]
#     best_comb <- combination[maxPval, ]
#     best_comb[3] <- Pvalue
#     names(best_comb)[1] <- "Stat0_radius"
#     names(best_comb)[2] <- "Stat1_radius"
#     names(best_comb)[3] <- "Pvalue"
#     return(best_comb) #also return the samples from sampling_all_radii
#   }
# } #return a list of best radii for both survival status 

#new 
get_best_radius <- function(input_data) {
  sampled_0 <- sampling_all_radii(input_data, 0) #running speed: slow
  sampled_1 <- sampling_all_radii(input_data, 1)
  
  stat0_radii <- sort(unique(sampled_0$radius)) #all radii for status 0
  stat1_radii <- sort(unique(sampled_1$radius)) #all radii for status 1
  combination <- expand.grid(stat0_radii, stat1_radii) # generate combinations based on the index in dataframe
  if (nrow(combination)==0) {
    print("No radius generated. Please try a bigger dataset.")
    return(NULL)
  } else {
    # pvalue to be stored
    cat("Estimating pval to find the best radii. Please wait patiently...")
    cat("\n")
    p_list <- sapply(1:nrow(combination), test, input_comb=combination, input_sampled_0=sampled_0, input_sampled_1=sampled_1)
    maxPval <- which.max(p_list)
    Pvalue <- p_list[maxPval]
    best_comb <- combination[maxPval, ]
    best_comb[3] <- Pvalue
    names(best_comb)[1] <- "Stat0_radius"
    names(best_comb)[2] <- "Stat1_radius"
    names(best_comb)[3] <- "Pvalue"
    #collect the samples
    new_sample <- rbind(sampled_0[sampled_0$radius==best_comb$Stat0_radius,], sampled_1[sampled_1$radius==best_comb$Stat1_radius,])
    print(best_comb)
    return(new_sample) #also return the samples from sampling_all_radii
  }
} 

test <- function(which_row, input_comb, input_sampled_0, input_sampled_1) {
  var1 <- input_comb[which_row, 1]
  df_1 <- input_sampled_0[which(input_sampled_0$radius == var1), ] 
  
  var2 <- input_comb[which_row, 2]
  df_2 <- input_sampled_1[which(input_sampled_1$radius == var2), ]
  
  df <- rbind(df_1, df_2)
  pvalue <- get_survPlot(df)[2]
  one <- input_comb[which_row,]
  one <- paste(unlist(one), collapse=' ')
  names(pvalue) <- one
  return(pvalue)
}


#used
sampling_all_radii <- function(input_data, event = 0) {
  comb_df <- data.frame(matrix(ncol = 5)) # create a dataframe to store all combinations and label them with an extra column
  colnames(comb_df) <- colnames(input_data)
  colnames(comb_df)[5] <- "radius"
  radii <- get_radiusRange(input_data, event) # get all possible radii for status 0 and 1
  # count <- 0
  for (radius in radii) {
    sampled_df <- sampling_event(input_data,radius,event)
    if (empty(sampled_df) == TRUE) {
      next; #if nothing to sample, do nothing
    } else {
      sampled_df$radius <- radius #assign its radius to each sampled dataset
      comb_df <- rbind(comb_df, sampled_df)
    }
  }
  comb_df <- na.omit(comb_df)
  return(comb_df)
}
