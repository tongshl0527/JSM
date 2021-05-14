
#source code
source('~/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/v1.0/Function_v1.7.2.R')

# read input
input_data <- read.delim("/Users/shanlin/Library/Mobile Documents/com~apple~CloudDocs/Pilot project/Survival_Dataset.txt")
control_data <- read.delim("/Users/shanlin/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/GEO_breast_cancer_OS_data.txt")
#preprocess
input_data <- format_inputdata(input_data, control_data)
input_data <- preprocess(input_data)
list <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
benchmarked_dataset <- main_benchmark(input_data, list, 3)
#write.csv(benchmarked_dataset,"benchmarking.csv", row.names = FALSE)
saveRDS(benchmarked_dataset, "benchmarked_dataset.rds")
