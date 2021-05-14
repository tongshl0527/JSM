#source code
source('~/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/v1.0/Function_v1.7.2.R')

# read input
input_data <- read.delim("/Users/shanlin/Library/Mobile Documents/com~apple~CloudDocs/Pilot project/Survival_Dataset.txt")
control_data <- read.delim("/Users/shanlin/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/GEO_breast_cancer_OS_data.txt")
#preprocess
input_data <- format_inputdata(input_data, control_data)
input_data <- preprocess(input_data)
# first preprocess and create a test file with 5 samples per batch group
#input_data <- input_data %>% group_by(BATCH) %>% sample_n(50,replace = TRUE)
#remove duplicated rows
input_data <- unique(input_data)
#execute function 
batch_effect_free_dataset <- new_main_all_in_one(input_data) 
#sampling inputdata randomly according to sizes in list
list <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)

saveRDS(batch_effect_free_dataset, file = "batch_effect_free_dataset_v1.7.1.rds")
# Writing mtcars data
write.table(input_data, file = "/Users/shanlintong/Library/Mobile Documents/com~apple~CloudDocs/Pilot project/Processed_Survival_Dataset.txt", sep = "\t", row.names = FALSE, col.names = TRUE)





