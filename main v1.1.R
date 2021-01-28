#source code
source('~/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/v1.0/Functions v1.0.R')
# read input
input_data <- read.delim("/Users/shanlintong/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/GEO_breast_cancer_OS_data.txt")
#preprocess
#input_data <- preprocess(input_data)
# first preprocess and create a test file with 5 samples per batch group
#input_data <- input_data %>% group_by(BATCH) %>% sample_n(50,replace = TRUE)
#remove duplicated rows
#input_data <- unique(input_data)
#execute function 
batch_effect_free_dataset <- main_new(input_data)
