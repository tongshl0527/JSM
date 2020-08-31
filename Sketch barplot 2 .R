library(ggplot2)

GEO_breast_cancer_OS_data <- read.delim("~/Documents/TTT-project Bioinformatics/GEO_breast_cancer_OS_data.txt")
GEO_breast_cancer_OS_data
# descriptive statistics
summary(GEO_breast_cancer_OS_data) 
# Remove all missing data
dataset <- na.omit(GEO_breast_cancer_OS_data) 
#sort the order of survival years 
dataset = dataset[order(dataset[,2]),] 
dataset_event1 <- subset(dataset, )

#barplot uses ggplot to show total events 
ggplot(dataset,aes(x=dataset$OVERALL_SURVIVAL_TIME_YEARS,
                                      y=dataset$OVERALL_SURVIVAL_EVENT))
                                                          +geom_bar(stat="identity") 

#subset dataset to event0 and event1 dataset
dataset_event1 <- dataset[ which(dataset$OVERALL_SURVIVAL_EVENT=='1'),]
dataset_event0 <- dataset[ which(dataset$OVERALL_SURVIVAL_EVENT=='0'),]
dataset_event0

# event 1 barplot 
ggplot_dataset_event1 <- dataset_event1
ggplot_dataset_event1$OVERALL_SURVIVAL_EVENT <- 1
ggplot(ggplot_dataset_event1,aes(x=ggplot_dataset_event1$OVERALL_SURVIVAL_TIME_YEARS,
                                 y=ggplot_dataset_event1$OVERALL_SURVIVAL_EVENT))
                                  +geom_bar(stat="identity") 
# event 0 barplot 
ggplot_dataset_event0 <- dataset_event0
ggplot_dataset_event0$OVERALL_SURVIVAL_EVENT <- 1
ggplot(ggplot_dataset_event0,aes(x=ggplot_dataset_event0$OVERALL_SURVIVAL_TIME_YEARS,
                                 y=ggplot_dataset_event0$OVERALL_SURVIVAL_EVENT))
                                  +geom_bar(stat="identity") 




############################################Ignore from here#############################################
#generate a list of indexes match to event 1
Event1_index <- which(dataset$OVERALL_SURVIVAL_EVENT==1)

#define 
sorted_vec = sort(dataset[,2], index.return = TRUE)
sorted_vec
sorted_vec_v1  = sorted_vec$x #here i am losed by Item x 
sub_vec_index = which(dataset$event==1)
sorted_vec_v2 = sorted_vec_v1
#find %IN% table.
#Find:The element(s) to look up in the vector or matrix. 
#table:The vector or matrix in which to look up the element(s).
sorted_vec_v2[which(sorted_vec$ix%in%sub_vec_index)] = 1
sorted_vec_v2[which(!(sorted_vec$ix%in%sub_vec_index))] = 0
sorted_vec_v2
########################################################################################################
barplot(dataset_event1$OVERALL_SURVIVAL_EVENT)



