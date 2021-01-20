# read input
input_data <- read.delim("/Users/shanlintong/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/GEO_breast_cancer_OS_data.txt")
#preprocess
input_data <- preprocess(input_data)
# first preprocess and create a test file with 50 samples per batch group
input_data <- input_data %>% group_by(BATCH) %>% sample_n(50,replace = TRUE)
# first preprocess and create a test file with 10 samples per batch group
input_data <- input_data %>% group_by(BATCH) %>% sample_n(10,replace = TRUE)
#remove duplicated rows
input_data <- unique(input_data)
#debugging
options(error = recover)
# first preprocess and create a test file with 5 samples per batch group
input_data <- input_data %>% group_by(BATCH) %>% sample_n(5,replace = TRUE)
