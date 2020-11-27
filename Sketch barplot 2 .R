library(ggplot2)
library(tibble)

GEO_breast_cancer_OS_data <- read.delim("~/Documents/TTT-project Bioinformatics/GEO_breast_cancer_OS_data.txt")
GEO_breast_cancer_OS_data
# descriptive statistics
summary(GEO_breast_cancer_OS_data) 
# Remove all missing data
dataset <- na.omit(GEO_breast_cancer_OS_data) 
#sort the order of survival years 
dataset = dataset[order(dataset[,2]),] 
#dataset_event1 <- subset(dataset, )

# predifine color blind friendly colors in The palette with black:
color_table <- tibble(Batch = c("GPL570_GSE16446", "GPL570_GSE20685", "GPL570_GSE20711", "GPL570_GSE42568", "GPL570_GSE48390", "GPL96_GSE1456", "GPL96_GSE3494", "GPL96_GSE45255", "GPL96_GSE4922","GPL96_GSE6532", "GPL96_GSE7390"),Color = c("yellow", "darkgreen", "blue4", "lightblue", "maroon3","orange","gray","lightsalmon","red","tan","peru"))
dataset$BATCH <- factor(dataset$BATCH, levels = color_table$Batch)

# All breast cancer patients event 1 #color=factor(dataset$BATCH) is removed
All_plot <- ggplot(dataset,aes(x=rank(dataset$OVERALL_SURVIVAL_TIME_YEARS,ties.method = "random"), y=dataset$OVERALL_SURVIVAL_EVENT,fill=BATCH, alpha=0.1)) +geom_bar(stat="identity") 
#remove background and axises + add axis titles
All_plot <- All_plot + theme(panel.grid.major = element_blank(),panel.background = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank() ,axis.text.y=element_blank(),axis.ticks.y=element_blank()) + labs(x = "Survival Time", y = "Survival Event", title = "All survival Breast Cancer Patients")
All_survival <- All_plot
All_survival
# width big and short height


# All breast cancer patients event 0
# reverse event 0 to event 1
event0 <- 1- dataset$OVERALL_SURVIVAL_EVENT
All_dicease <- ggplot(dataset,aes(x=rank(dataset$OVERALL_SURVIVAL_TIME_YEARS,ties.method = "random"), y=event0, fill=BATCH, alpha=0.1)) +geom_bar(stat="identity") 
#remove background and axises + add axis titles
All_dicease <- All_dicease + theme(panel.grid.major = element_blank(), panel.background = element_blank(),panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank() ,axis.text.y=element_blank(),axis.ticks.y=element_blank()) + labs(x = "Survival Time", y = "Survival Event", title = "All Diceased Breast Cancer Patients")
All_dicease

#Export figures to pdf file 
pdf(file = "/Users/shanlin/Documents/TTT-project Bioinformatics/survival.pdf",height = 4, width = 8)
print(All_survival)
dev.off()

pdf(file = "/Users/shanlin/Documents/TTT-project Bioinformatics/diceased.pdf",height = 4, width = 8)
print(All_dicease)
dev.off()

#################################################### sketch for batch GSE16446 #########################################################################################

#Extract batch from dataset in ggplot. 

# select the first batch 
GSE16446_data <- dataset[ which(dataset$BATCH=='GPL570_GSE16446'), ]
GSE16446_data

GSE16446_plot <- ggplot(GSE16446_data,aes(x=rank(GSE16446_data$OVERALL_SURVIVAL_TIME_YEARS,ties.method = "random"), y=GSE16446_data$OVERALL_SURVIVAL_EVENT, alpha=0.5)) +geom_bar(stat="identity") 
#remove background and axises + add axis titles
GSE16446_plot <- GSE16446_plot + theme(panel.grid.major = element_blank(),panel.background = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank() ,axis.text.y=element_blank(),axis.ticks.y=element_blank()) + labs(x = "Survival Time", y = "Survival Event", title = "GSE16446 survival Breast Cancer Patients")
GSE16446_survival <- GSE16446_plot
GSE16446_survival
# width big and short height
GSE16446_event0 <- 1- GSE16446_data$OVERALL_SURVIVAL_EVENT
GSE16446_dicease <- ggplot(GSE16446_data,aes(x=rank(GSE16446_data$OVERALL_SURVIVAL_TIME_YEARS,ties.method = "random"), y=GSE16446_event0, alpha=0.5)) +geom_bar(stat="identity") 
#remove background and axises + add axis titles
GSE16446_dicease <- GSE16446_dicease + theme(panel.grid.major = element_blank(), panel.background = element_blank(),panel.grid.minor = element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank() ,axis.text.y=element_blank(),axis.ticks.y=element_blank()) + labs(x = "Survival Time", y = "Survival Event", title = "GSE16446 Diceased Breast Cancer Patients")
GSE16446_dicease

#Export figures to pdf file 
pdf(file = "/Users/shanlin/Documents/TTT-project Bioinformatics/GSE16446 survival.pdf",height = 4, width = 8)
print(GSE16446_survival)
dev.off()

pdf(file = "/Users/shanlin/Documents/TTT-project Bioinformatics/GSE16446 diceased.pdf",height = 4, width = 8)
print(GSE16446_dicease)
dev.off()
