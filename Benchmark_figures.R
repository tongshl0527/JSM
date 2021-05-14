library(reshape2)
source('~/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/v1.0/Function_v1.7.2.R')

benchmark <- read.csv("/Users/shanlin/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/Benchmarks/benchmarking.csv", header = TRUE)
benchmark_dataset <- "/Users/shanlin/Library/Mobile Documents/com~apple~CloudDocs/TTT-project Bioinformatics/Benchmarks/benchmarked_dataset.rds"
#load benchmarked dataset
file <- readRDS(benchmark_dataset, refhook = NULL)
#load summary
summary <- as.data.frame(file[1])
colnames(summary) <- removeWords(colnames(summary),"summary.")
#load 0.75_1
aft.corr <- extract_survivalData_from_main_benchmark(file, "0.75_1")
bef.corr <- extract_survivalData_from_main_benchmark(file, "0.75_1", FALSE)
#KM
KM_Plot(bef.corr)
KM_Plot(aft.corr)

df <- benchmark_preprocess(summary)
#process dataset
benchmark_preprocess <- function(benchmark) {
  df <- data.frame(benchmark$Total.Sample.Size.before, benchmark$Sample.Size.event.0.before, benchmark$Sample.Size.event.1.before, benchmark$Pval.before, benchmark$Event.0.radius, benchmark$Event.1.radius, benchmark$Time)
  colnames(df) <- sub(".before.*", "", colnames(df))
  df$benchmark.correct.status <- "before"
  df2 <- data.frame(benchmark$Total.Sample.Size.before, benchmark$Sample.Size.event.0.after, benchmark$Sample.Size.event.1.after, benchmark$Pval.after, benchmark$Event.0.radius, benchmark$Event.1.radius, benchmark$Time)
  colnames(df2) <- sub(".after.*", "", colnames(df2))
  colnames(df2) <- sub(".before.*", "", colnames(df2))
  df2$benchmark.correct.status <- "after"
  df <- rbind(df, df2)
  return(df)
}

pval_bef_vs_aft <- function(df) {
  #summary for pvalues before vs after correction 
  df.summary <- data_summary(df, varname="benchmark.Pval", groupnames = c("benchmark.Total.Sample.Size","benchmark.correct.status"))
  #dotplot
  e <- ggplot(df, aes(benchmark.Total.Sample.Size, benchmark.Pval, color = benchmark.correct.status))
  # Dotplot with summary statistics: mean +/- SD
  result <- e +   geom_jitter(position = position_jitter(0.2)) + 
    geom_line(aes(group = benchmark.correct.status),data = df.summary) +
    #geom_errorbar(aes(ymin = benchmark.Pval-sd, ymax = benchmark.Pval+sd), data = df.summary, width = 0.2)+
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    theme(legend.position = "top")  +
    scale_x_continuous(breaks = df$benchmark.Total.Sample.Size)+
    scale_y_continuous(limits = c(0, 1))+labs(title="Log Rank Test Per Sample Size Before VS After Correction", x="Total Sample Size Before Correction", y = "Pval Log Rank Test")+
    theme_classic() 
  return(result)
}


#summary for speed testing 
plot_running_time <- function(df) {
  df.summary2 <- data_summary(df, varname="benchmark.Time", groupnames = "benchmark.Total.Sample.Size")
  # timer
  p <- ggplot(df.summary2, aes(benchmark.Total.Sample.Size, benchmark.Time)) +
    geom_jitter( position = position_jitter(0.2), color = "darkgray") + 
    geom_line(aes(group = 1), data = df.summary2) +
    geom_errorbar(
      aes(ymin = benchmark.Time-sd, ymax = benchmark.Time+sd),
      data = df.summary2, width = 0.2) +
    geom_point(data = df.summary2, size = 2)+
    scale_x_continuous(breaks = df$benchmark.Total.Sample.Size)+
    labs(title="Speed test", x="Total Sample Size Before Correction", y = "User Time")+
    theme_classic() 
  return(p)
}

samp_size_bef_vs_aft <- function(benchmark) {
  dfm <- melt(benchmark[,c('Total.Sample.Size.before','Sample.Size.event.0.before','Sample.Size.event.0.after','Sample.Size.event.1.before','Sample.Size.event.1.after')],id.vars = 1)
  
  result <- ggplot(dfm,aes(x = Total.Sample.Size.before,y = value)) + 
    geom_bar(aes(fill = variable),stat = "identity",position = "dodge") +
    scale_x_continuous(breaks = df$benchmark.Total.Sample.Size)+
    labs(title="Sample size in event status before vs after sampling", x="Total Sample Size Before Correction", y = "Sample Size per event status")+
    theme_classic() 
  return(result)
}

#summary for before vs after sample size 
df.summary3 <- data_summary(benchmark, varname="Total.Sample.Size.after", groupnames = "Total.Sample.Size.before")
df.summary3$radii <- c()
# Bar plots + jittered points + error bars
ggplot(df.summary3, aes(Total.Sample.Size.before, Total.Sample.Size.after)) +
  geom_col(data = df.summary3, fill = NA, color = "black") +
  geom_jitter( position = position_jitter(0.2), color = "black") + 
  geom_errorbar( aes(ymin = Total.Sample.Size.after-sd, ymax = Total.Sample.Size.after+sd), 
                 data = df.summary3, width = 0.2) +
  scale_x_continuous(breaks = df$benchmark.Total.Sample.Size)+
  labs(title="Sample size before vs after sampling", x="Total Sample Size Before Correction", y = "Total Sample Size After Correction")+
  theme_classic() 

#+ geom_text(aes(label=df$benchmark.Event.0.radius))


#radii
#df.summary4 <- data_summary(dfm2, varname="value", groupnames = "Total.Sample.Size.before")
dfm2 <- melt(benchmark[,c('Total.Sample.Size.before','Event.0.radius','Event.1.radius')],id.vars = 1)
ggplot(dfm2,aes(x = Total.Sample.Size.before,y = value)) + 
  geom_point(aes(color = variable),stat = "identity",position = "dodge") +
  scale_x_continuous(breaks = df$benchmark.Total.Sample.Size)+
  labs(title="Optimal radii per sample size", x="Total Sample Size Before Correction", y = "Radii per sample size")+
  theme_classic() 

print(p)
# Finished line plot
p+labs(title="Log Rank Test Per Sample Size Before VS After Correction", x="Total Sample Size Before Correction", y = "Pval Log Rank Test")+
  theme_classic() 



#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  data_sum$correction.status <- varname
  return(data_sum)
}


library(survival)  
get_survPlot(batch_effect_free_dataset_1$`Sampled Dataset`)
