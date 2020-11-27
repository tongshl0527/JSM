# Load required packages
library(survival)
library(survminer)
library(dplyr)
library(plyr)

# Import the breast cancer dataset and have a look at it
glimpse(df)

###################################### Fit survival data using the Kaplan-Meier method ###########################################
surv_object <- Surv(time = df$OVERALL_SURVIVAL_TIME_YEARS, event = df$OVERALL_SURVIVAL_EVENT)
surv_object

ptable <- pairwise_survdiff(Surv(OVERALL_SURVIVAL_TIME_YEARS,OVERALL_SURVIVAL_EVENT) ~ BATCH, data = df)
ptable
write.table(ptable$p.value, "/Users/shanlin/Documents/TTT-project Bioinformatics/pairwise comparision ptable.txt", sep="\t", row.names = TRUE)

#get a distance matrix for survival time 
distance_matrix <- as.matrix(dist(df$OVERALL_SURVIVAL_TIME_YEARS))
sorted_disMatrix <- apply(distance_matrix,1,sort)
radius <- apply(sorted_disMatrix, 1, min)
radius <- unique(radius)[2:length(radius)]
write.table(radius,"/Users/shanlin/Documents/TTT-project Bioinformatics/radius.txt", sep="\t", row.names = TRUE)

###################################Get the smallest radius################################
# this will be built as a function
# step 0, use the minimal range

get_minimal_radius <- function(radius_input) {
  n <- length(radius_input)
  for (i in 1:n) {
    if (radius_input[i] > 0) {
      minRadius <- radius_input[i]
      return(minRadius)
      break
    }
  }
}

min_rad <- get_minimal_radius(radius)

######################################Actual_distribution_samples_in_batches splitted according to event status#########################################
#calculate Actual_distribution_samples_in_batches
#count total samples 
freq_total <- count(df, c("BATCH","OVERALL_SURVIVAL_EVENT"))
#count the total samples in event 0
freq_event0 <- freq_total[freq_total$OVERALL_SURVIVAL_EVENT == "0",]
sum_event0 <- sum(freq_event0$freq)
#count the samples in event 1 for each batch
freq_event1 <- freq_total[freq_total$OVERALL_SURVIVAL_EVENT == "1",]
sum_event1 <- sum(freq_event1$freq)

#calculate ratio in event 0 for each batch 
n <- length(freq_event0$freq)
standard_ratio_0 <- list()
for (i in n) {
  freq_batch_event0 <- freq_event0$freq[i]
  ratio <- freq_batch_event0/sum_event0
  names(standard_ratio_0)[i] <- freq_event0$BATCH[i]
  standard_ratio_0[i] <- ratio
}
print(standard_ratio_0)


#calculate ratio in event 1 for each batch
m <- length(freq_event1$freq)
standard_ratio_1 <- list()
for (i in 1:m) {
  freq_batch_event1 <- freq_event1$freq[i]
  ratio <- freq_batch_event1/sum_event1
  names(standard_ratio_1)[i] <- freq_event1$BATCH[i]
  standard_ratio_1[i] <- ratio
}
print(standard_ratio_1)


#################################################Select time interval based on minimal radius############
# step 1, sum up all frequencies in that time range for event 0 and 1
df0 <- subset(df, OVERALL_SURVIVAL_EVENT==0)
df1 <- subset(df, OVERALL_SURVIVAL_EVENT==1)


#calculate the number of time intervals 
n_TI_0 <- max(df0$OVERALL_SURVIVAL_TIME_YEARS)/min_rad
n_TI_1 <- max(df1$OVERALL_SURVIVAL_TIME_YEARS)/min_rad


#scan through each time interval in df0
# create a frequency matrix
frequency_matrix_event0 <- matrix(0, ncol = n_TI_0, nrow = length(Batch))

rownames(frequency_matrix_event0) <- Batch


#create a matrix to store sample size 
for (i in seq(from = min_rad, to = n_TI_0, by= min_rad)) {
  #subset out the interval 
  #based on interval, select the right range of years. Batch is a list for all batches
  time_interval <- subset(df0, OVERALL_SURVIVAL_TIME_YEARS < i & OVERALL_SURVIVAL_TIME_YEARS >= i-min_rad) 
  #scan through the subset, count the number of batch appeared, otherwise return that batch to zero 
  for (j in length(time_interval$BATCH)) {
      #add 1
      index <- match(time_interval$BATCH[j],Batch)
      frequency_matrix_event0[index, i] <- frequency_matrix_event0[index, i] +1
  }
}
write.table(frequency_matrix_event0,"/Users/shanlin/Documents/TTT-project Bioinformatics/sample size per time interval in event 0.txt", sep="\t", row.names = TRUE, col.names = TRUE)


ratio_matrix_event0 <- matrix(0, ncol = n_TI_0, nrow = length(Batch))
rownames(ratio_matrix_event0) <- Batch
#store ratio to a new matrix
for (k in n_TI_0) {
  sumCol <- sum(frequency_matrix_event0[,k])
  for (z in length(Batch)) {
    ratio_matrix_event0[z,k] <- frequency_matrix_event0[z,k]/sumCol
  }
}
write.table(ratio_matrix_event0,"/Users/shanlin/Documents/TTT-project Bioinformatics/ratio per time interval in event 0.txt", sep="\t", row.names = TRUE, col.names = TRUE)


#test ratio of each batch against standard ratio in each time interval

check_matrix <- matrix(0, ncol(ratio_matrix_event0), nrow(ratio_matrix_event0))
for(row in 1:nrow(ratio_matrix_event0)) {
    for(col in 1:ncol(ratio_matrix_event0)) {
      if (ratio_matrix_event0[row, col] == standard_ratio_1[row]) {
        check_matrix[row, col] <- TRUE
      }
      else {
        check_matrix[row, col] <- FALSE
      }
    }
  }

################################################Follow the same as above for event 1#####################################################################

