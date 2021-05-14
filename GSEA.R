library(dplyr)
library(reshape2)
library(pastecs)
library(ggplot2)
library(scales)
library(survival)
library(survminer)
library(parallel)
library(MASS)
library(pheatmap)
library(fgsea)
library(msigdbr)

get_gsea_all_databases <- function(species = "Homo sapiens", microarrayExp, survival_data, Gene_Map) {
  all_pathways <- get_all_pathways(species)
  geneExp_matrix <- get_geneExp_matrix(microarrayExp, Gene_Map)
  ranked_gene_list <- get_ranked_gene_list(geneExp_matrix, survival_data)
  gsea <- get_gsea(all_pathways, ranked_gene_list) #all databases are stored in a list 
  return(gsea)
}

#load genesets https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
#create a function that can be used by lapply to generate a list of all database analysis
# check for input: msigdbr_species()
# check for collection: msigdbr_collections()
get_all_pathways <- function(species = "Homo sapiens") {
  collections <- msigdbr_collections()
  l <- list()
  count <- 0
  # generate a list of subcat 
  for (i in unique(collections$gs_cat)) {
    subcategories = collections$gs_subcat[collections$gs_cat==i]
    for (j in subcategories) {
      count <- count +1
      category <- i
      subcategory <- j
      gene_sets = msigdbr(species = species, category = category, subcategory = subcategory)
      gene_sets = split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)
      l[[count]] <- gene_sets
      names(l)[[count]] <- paste(category, subcategory, sep = "")
    }
  }
  return(l)
}

#fgseaRes <- fgsea(pathways.hallmark,ranked_gene_list, nperm=1000)
get_gsea <- function(all_pathways, ranked_gene_list) {
  a <- lapply(all_pathways, fgsea, stat = ranked_gene_list, nperm=1000)
  return(a)
}

#export to csv

export_gsea_NES <- function(gsea) {
  result <- lapply(gsea, function(x) x[,c("pathway","NES")])
  for(i in names(result)){
    write.csv(result[[i]], paste0(i,".csv"))
  }
  return(result)
}

#show top pathways 
#select your favorite from all_pathways list
plot_top10_pathways <- function(ExamplePathways, ExampleRankList) {
  fgseaRes <- fgsea(pathways = ExamplePathways,
                    stats = ExampleRankList, 
                    #eps = 1e-10, #do not suggest to lower pvalue threshold: more accurate as it closes to zero
                    nperm=1000)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plot <- plotGseaTable(ExamplePathways[topPathways], ExampleRankList, fgseaRes, 
                        gseaParam=0.5)
  return(plot)
}


#
get_geneExp_matrix <- function(microarrayExp, Gene_Map) {
  #fix NA in BreastAndFat data
  microarrayExp[1,2:741] <- microarrayExp[1,1:740]
  microarrayExp[1,1] <- "GSM_NAMES"
  #microarrayExp <- microarrayExp
  rownames(microarrayExp) <- microarrayExp[,1]
  microarrayExp <- microarrayExp[,-1]
  colnames(microarrayExp) <- microarrayExp[1,]
  microarrayExp <- microarrayExp[-1,]
  
  #replace rownames from probeset to geneID 
  matched_genesID <- Gene_Map$ENTREZID[match(rownames(microarrayExp), Gene_Map$PROBESET)]
  rownames(microarrayExp) <- matched_genesID
  geneExp_matrix <- as.data.frame(sapply(microarrayExp, as.numeric))
  rownames(geneExp_matrix) <- rownames(microarrayExp)
  geneExp_matrix <- as.data.frame(t(geneExp_matrix))
  return(geneExp_matrix)
}

get_ranked_gene_list <- function(geneExp_matrix, survival_data) {
  #scale genes
  scaled_geneExp <- as.data.frame(scale(geneExp_matrix))
  #merge scaled gene data to survival 
  survival <- survival_data
  survival$geneID <- scaled_geneExp[match(rownames(scaled_geneExp), survival$GSM_IDENTIFIER),]
  #coxph function analysis on univariant genes
  covariates <- paste('survival$geneID$"',colnames(survival$geneID),'"', sep = "")
  univ_formulas <- sapply(covariates, function(x) as.formula(paste('Surv(survival$OS_time_months, survival$OS_event) ~',x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x,data=survival)}) #very slow
  #extract data
  univ_results <- lapply(univ_models, function(x){
    x <- summary(x)
    p.value<-signif(x$wald["pvalue"], digits=2)
    wald.test<-signif(x$wald["test"], digits=2)
    beta<-signif(x$coef[1], digits=2);#coeficient beta
    HR <-signif(x$coef[2], digits=2);#exp(beta)
    HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-c(beta, HR, wald.test, p.value)
    names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                  "p.value")
    return(res)
  })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  gene_survival_res <- as.data.frame(res)
  #remove non-numeric characters in a string 
  rownames(gene_survival_res) <- gsub("[^0-9.-]", "", rownames(gene_survival_res))
  res$ENTREZID <- parse_number(rownames(res))
  #give association scores
  num_matrix <- as.data.frame(sapply(gene_survival_res,as.numeric))
  rownames(num_matrix) <- rownames(gene_survival_res)
  num_matrix$ass.score <- -log10(num_matrix$p.value)*sign(num_matrix$beta)
  #rank from low to high
  ranked_genes <- num_matrix[with(num_matrix, order(num_matrix$ass.score)),]
  #add gene names
  matched_genesSym <- Gene_Map$SYMBOL[match(rownames(ranked_genes),Gene_Map$ENTREZID)]
  ranked_genes$Symbol <- matched_genesSym
  #GESA website tutorial: https://stephenturner.github.io/deseq-to-fgsea/#load_some_results_from_deseq2
  #make a gene list
  ranked_gene_list <- ranked_genes$ass.score
  names(ranked_gene_list) <- ranked_genes$Symbol
  ranked_gene_list <- sort(ranked_gene_list, decreasing = TRUE)
  ranked_gene_list <- ranked_gene_list[!duplicated(names(ranked_gene_list))]
  return(ranked_gene_list)
}


#group GSM identifiers based on batch groups
#samples_per_batch_group <- survival_data %>% group_by(GSE_number) %>% summarize(type = paste(sort(unique(GSM_IDENTIFIER)),collapse=", "))

#organize microarrayExp


#import hallmarks
# pathways.hallmark <- gmtPathways("/Users/shanlintong/Downloads/h.all.v7.4.symbols.gmt")
#run GESA
# fgseaRes <- fgsea(pathways.hallmark,ranked_gene_list, nperm=1000)
# #show in a table
# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(desc(NES))
# #plot hallmarks
# the_plot <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill=padj<0.05)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Hallmark pathways NES from GSEA") + 
#   theme_minimal()


#check top 10 hallmarks
# top10_hallmarks <- fgseaRes %>% arrange(desc(abs(NES))) %>% top_n(10, -padj)
# 
# # Enrichmentscore plot
# plotEnrichment(pathways.hallmark[["HALLMARK_MYC_TARGETS_V2"]], ranked_gene_list) + labs(title="MYC_TARGETS")
# plotEnrichment(pathways.hallmark[["HALLMARK_UV_RESPONSE_DN"]], ranked_gene_list) + labs(title="UV_RESPONSE")

# #make a ranking table, rank from top, rank from buttom
# ranking_table <- as.data.frame(rownames(ranked_genes))
# colnames(ranking_table)[1] <- "bottom_to_top"
# top_to_bottom <- ranked_genes[with(ranked_genes, order(-ranked_genes$ass.score)),]
# ranking_table$top_to_bottom <- rownames(top_to_bottom)
# 
# #export table 
# write.csv(ranking_table,"~/Library/Mobile Documents/com~apple~CloudDocs/Pilot project/ranking_table.csv", row.names = TRUE)
# 
# #draw heatmap from gene set enrichment analysis GENETICA
# rawdata_file <- "/Users/shanlintong/Library/Mobile Documents/com~apple~CloudDocs/Pilot project/GESA/Hallmark/Hallmark.csv"

draw_heatmap <- function(input_rawdata, title = "") {
  rawdata <- read.csv(input_rawdata)
  df_matrix <- acast(rawdata, gene_set~gene_name, value.var="value")
  scaled_matrix <- scale(df_matrix)
  heatmap <- pheatmap(scaled_matrix[1:20,], main = title, show_colnames =FALSE)
  return(heatmap)
}

# draw_heatmap(rawdata_file,"Transcriptome Signature of Top 100")
# 
# #draw heatmap based on z-scores 
# file <- "/Users/shanlintong/Downloads/BIOCARTA.csv"
# co_function_data <- read.csv(file)
# z_score_df <- data.frame(matrix(NA, nrow = length(unique(co_function_data$gene_set)), ncol = length(unique(co_function_data$gene_name))))
# rownames(z_score_df) <- unique(co_function_data$gene_set)
# colnames(z_score_df) <- unique(co_function_data$gene_name)
# #subset 
# for (i in 1:ncol(z_score_df)) {
#   a <- co_function_data %>% filter(gene_name == colnames(z_score_df)[i])
#   z_score_df[,i] <- a$value
# }

#descriptive statistics on samples
descrp_stat <- function(input_geneExp) {
  a <- data.frame(lapply(input_geneExp,as.numeric))
  des_stat <- as.data.frame(apply(a, 2, summary))
  return(des_stat)
}

# allOverview <- descrp_stat(microarrayExp) 
# 
# #convert the whole gene dataframe to a matrix
# geneExp_matrix <- as.matrix(sapply(geneExp, as.numeric))
# hist.obj <- hist(geneExp_matrix, freq = TRUE, main = "Histogram of Gene Expression", xlab = "Gene Expression")
# str(hist.obj)
# 
# #convert samples of the min,max,mean,median to a matrix
# plot_list = list()
# for (i in 1:6) {
#   sample_matrix <- as.matrix(sapply(allOverview[i,], as.numeric))
#   file_name <- paste("samples_histogram_", rownames(allOverview)[i], ".tiff", sep="")
#   tiff(file = paste("/Users/shanlintong/Library/Mobile Documents/com~apple~CloudDocs/Pilot project/",file_name,sep=""))
#   hist(sample_matrix , freq = TRUE, main = paste("Histogram for the", rownames(allOverview)[i], "of samples", sep=" "), xlab = "Gene Expression")
#   dev.off()
# }
# 
# #convert genes of the min,max,mean,median to a matrix
# microarrayExp <- t(microarrayExp)
# microarrayExp <- as.data.frame(microarrayExp)
# gene_overview <- descrp_stat(microarrayExp)
# for (i in 1:6) {
#   gene_matrix <- as.matrix(sapply(gene_overview[i,], as.numeric))
#   file_name <- paste("genes_histogram_", rownames(gene_overview)[i], ".tiff", sep="")
#   tiff(file = paste("/Users/shanlintong/Library/Mobile Documents/com~apple~CloudDocs/Pilot project/genes histogram",file_name,sep=""))
#   hist(gene_matrix, freq = TRUE, main = paste("Histogram for the", rownames(gene_overview)[i], "of genes", sep=" "), xlab = "Gene Expression")
#   dev.off()
# }
# 
# 
# #subset gene file based on batch group
# batch_list <- list() 
# count<-0
# for (variable in samples_per_batch_group$GSE_number) { #store each batch gene expression data to list
#   count<- count+1
#   batch_samples <- samples_per_batch_group[which(samples_per_batch_group$GSE_number==variable),]$type
#   batch_samples <- unlist(strsplit(batch_samples, ","))
#   batch_samples <- gsub(" ", "", batch_samples) 
#   batch_geneExp <- microarrayExp[batch_samples]
#   batch_list[[count]] <- batch_geneExp
#   names(batch_list)[count] <- variable
# }
# 
# #get mean for each gene 
# geneExp <- t(microarrayExp)
# geneExp <- as.data.frame(geneExp)
# des_stat_genes <- descrp_stat(geneExp) #descriptive statistics on genes
# des_stat_genes <- t(des_stat_genes)
# des_stat_genes <- as.data.frame(des_stat_genes)
# #rank genes based on mean
# ranked_genes <- des_stat_genes[order(-des_stat_genes$Mean),]
# #extract top50 genes
# top30_genes <- head(ranked_genes,30)
#   
# #histogram plot on means of top 30 genes 
# ggplot(top30_genes, aes(x= reorder(rownames(top30_genes), -Mean), y = top30_genes$Mean)) + geom_histogram(stat='identity') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(title="Mean", x ="Genes", y = "microRNA Expression")
# 
# # Make list of variable names to loop over.
# var_list = combn(names(top30_genes)[1:6], 1, simplify=FALSE)
# #make plots of top 30 genes ranked by means
# plot_list = list()
# for (i in 1:6) {
#   p <- ggplot(top30_genes, aes(x= reorder(rownames(top30_genes), -top30_genes[,i]), y = top30_genes[,i])) + geom_histogram(stat='identity') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + labs(title=var_list[[i]], x ="Genes", y = "microRNA Expression")
#   plot_list[[i]] <- p
# }
# 
# # Save plots to tiff. Makes a separate file for each plot.
# for (i in 1:6) {
#   file_name = paste("top30_genes_plot_", var_list[[i]], ".tiff", sep="")
#   tiff(file_name)
#   print(plot_list[[i]])
#   dev.off()
# }
# 
# #descriptive statistics on samples
# des_stat_samples <- descrp_stat(microarrayExp) 
# des_stat_samples <- t(des_stat_samples)
# des_stat_samples <- as.data.frame(des_stat_samples)
# #rank samples based on mean
# ranked_samples <- des_stat_samples[order(-des_stat_samples$Mean),]
# # Make list of variable names to loop over.
# var_list = combn(names(ranked_samples)[1:6], 1, simplify=FALSE)
# #make plots of samples ranked by means
# plot_list = list()
# for (i in 1:6) {
#   p <- ggplot(ranked_samples, aes(x= reorder(rownames(ranked_samples), -ranked_samples[,i]), y = ranked_samples[,i])) + geom_histogram(stat='identity') + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(title=var_list[[i]], x ="Samples", y = "microRNA Expression")
#   plot_list[[i]] <- p
# }
# 
# # Save plots to tiff. Makes a separate file for each plot.
# for (i in 1:6) {
#   file_name = paste("Samples_plot_", var_list[[i]], ".tiff", sep="")
#   tiff(file_name)
#   print(plot_list[[i]])
#   dev.off()
# }



#select batch 1: GSE16446
get_geneExp_descStat <- function(input_batch #select a batch from the batch list
) {
  batch_geneExp <- as.data.frame(t(input_batch))
  batch_geneExp_descStat <- descrp_stat(batch_geneExp)
  batch_geneExp_descStat <- as.data.frame(t(batch_geneExp_descStat))
  batch_geneExp_descStat <- batch_geneExp_descStat[order(-batch_geneExp_descStat$Mean),] #rank the mean of gene expression from high to low
  return(batch_geneExp_descStat)
}

get_samples_descStat <- function(input_batch) {
  batch_samples <- as.data.frame(input_batch)
  batch_samples_descStat <- descrp_stat(batch_samples)
  batch_samples_descStat <- as.data.frame(t(batch_samples_descStat))
  batch_samples_descStat <- batch_samples_descStat[order(-batch_samples_descStat$Mean),] #rank the mean of gene expression from high to low
  return(batch_samples_descStat)
}

get_geneExp <- function(input_batch #select a batch from the batch list
) {
  batch_geneExp <- as.data.frame(t(input_batch))
  #batch_geneExp_descStat <- descrp_stat(batch_geneExp)
  #batch_geneExp_descStat <- as.data.frame(t(batch_geneExp_descStat))
  batch_geneExp <- data.frame(lapply(batch_geneExp,as.numeric))
  return(batch_geneExp)
}

get_samples <- function(input_batch #select a batch from the batch list
) {
  batch_geneExp <- as.data.frame(input_batch)
  #batch_geneExp_descStat <- descrp_stat(batch_geneExp)
  #batch_geneExp_descStat <- as.data.frame(t(batch_geneExp_descStat))
  batch_geneExp <- data.frame(lapply(batch_geneExp,as.numeric))
  return(batch_geneExp)
}

plot_gene_distribution_of_batch <- function(input_microarrayExp, 
                                            input_samples_per_batch_group, 
                                            input_geneName) {
    #collect gene BRCA from datalist 
  #extract BRCA2
  samples_of_gene <- input_microarrayExp[input_geneName,]
  samples_of_gene_new <- data.frame(lapply(samples_of_gene,as.numeric))
  colnames(samples_of_gene_new) <- colnames(samples_of_gene)
  #make a new dataframe for each batch
  batch_samples_descStat <- descrp_stat(batch_samples)
  des_stat_df <- data.frame(nrow=6)
  batchList <- list()
  count<-0
  for (variable in input_samples_per_batch_group$GSE_number) { #store each batch gene expression data to list
    count<- count+1
    batch_samples <- input_samples_per_batch_group[which(input_samples_per_batch_group$GSE_number==variable),]$type
    batch_samples <- unlist(strsplit(batch_samples, ","))
    batch_samples <- gsub(" ", "", batch_samples) 
    batch_geneExp <- samples_of_gene_new[batch_samples]
    batch_geneExp <- t(batch_geneExp)
    batchList[count] <- batch_geneExp
    names(batchList)[count] <- variable
    #save batch gene exp to a dataframe
    
    
    #get a dis stat and save to a dataframe 
    des_stat <- as.data.frame(apply(batch_geneExp, 2, summary))
    des_stat_df[,count] <- des_stat
  }
  #plotting histogram mean vs median vs 
  return(boxplot(gene_df))
} 

