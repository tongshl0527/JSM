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


main_new <- function(input_unfiltered_data, #unprocessed data
                     test="chitest", #test: chitest, fishertest, cortest
                     input_radius_0=NULL, #not fill in which run based on best radius
                     input_radius_1=NULL #not fill in which run based on best radius
) {
  input_data <- preprocess(input_unfiltered_data)
  if (is.null(input_radius_0)==TRUE & is.null(input_radius_1)==FALSE) {
    print("Error: Please fill in radius for both event status")
  } else if (is.null(input_radius_0)==FALSE & is.null(input_radius_1)==TRUE){
    print("Error: Please fill in radius for both event status")
  } else if (is.null(input_radius_0)&is.null(input_radius_1)==TRUE){
    best_pairs <- get_best_radius_new(input_data, testmethod = test)
    input_radius_0 <- best_pairs[1]
    input_radius_1 <- best_pairs[2]
  } else {
    input_radius_0 <- input_radius_0
    input_radius_1 <- input_radius_1
  }
  sampled_df <- Sampling(input_unfiltered_data, test, input_radius_0, input_radius_1)
  print(KM_Plot(sampled_df))
  return(sampled_df)
}