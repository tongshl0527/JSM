

# Read in txt files
GEO_breast_cancer_OS_data <- read.delim("~/Documents/TTT-project Bioinformatics/GEO_breast_cancer_OS_data.txt")
# descriptive statistics
summary(GEO_breast_cancer_OS_data) 
# Remove all missing data
dataset <- na.omit(GEO_breast_cancer_OS_data) 
#sort the order of survival years 
dataset = dataset[order(dataset[,2]),] 
dataset_event1 <- subset(dataset, )

#subset dataset to event1 and event2 dataset
dataset_event1 <- dataset[ which(dataset$OVERALL_SURVIVAL_EVENT=='1'),]
dataset_event0 <- dataset[ which(dataset$OVERALL_SURVIVAL_EVENT=='0'),]

#barplot 
barplot(dataset_event1, space = 0, border = NA ,yaxt = "n", col = "black",las=1,xaxt = "n")

#I do not understand how these part of the code works 
genelevel_data = list()
survival_time = rnorm(1000, mean = 10, sd = 3)
event = sample(c(0,1),1000, replace = TRUE)
genelevel_data$survival_time = survival_time
genelevel_data$event = event
genelevel_data = as.data.frame(genelevel_data)


sorted_vec = sort(genelevel_data[,1]
                  # , decreasing = TRUE
                  , index.return = TRUE)
sorted_vec_v1  = sorted_vec$x
sub_vec_index = which(genelevel_data$event==1)

sorted_vec_v2 = sorted_vec_v1

sorted_vec_v2[which(sorted_vec$ix%in%sub_vec_index)] = 1

sorted_vec_v2[which(!(sorted_vec$ix%in%sub_vec_index))] = 0
barplot(sorted_vec_v2
        , space = 0
        , border = NA
        # , ylim = c(-1000, 1000)
        ,yaxt = "n"
        , col = "black"
        ,las=1
        ,xaxt = "n")