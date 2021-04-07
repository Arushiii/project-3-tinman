setwd("/projectnb/bf528/users/tinman/Project3/Analyst")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
library(limma)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_4_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
)


b<- samples[samples$chemical=='N-NITROSODIETHYLAMINE',]
control <- samples[samples$vehicle=="SALINE_100_%" & samples$chemical=="Control",]
g<- rbind(b, control)


rma_b <- rma[paste0('X',g$array_id)]
design<- model.matrix(~factor(g$chemical, levels = c("Control", "N-NITROSODIETHYLAMINE")))
colnames(design) <- c('Intercept','N-NITROSODIETHYLAMINE')

fit <- lmFit(rma_b, design)

fit <- eBayes(fit)

t <- topTable(fit, coef=2, n=nrow(rma_b), adjust='BH')

# write out the results to file
write.csv(t,'N-NITROSODIETHYLAMINE_limma_results.csv')

#sorting the data based on adjusted p-value
sorted_data<- t[order(t$adj.P.Val),]
#filtering the sorted data based on the threshold of 0.05
sorted_data1<- sorted_data[sorted_data$adj.P.Val<0.05,]

#fetching the top 10 hits.
top10<- sorted_data1[order(sorted_data1$P.Value),]
top10_hit<-top10[1:10,]
write.csv(top10_hit,'N-NITROSODIETHYLAMINE_hit.csv')
#creating the histogram
install.packages("ggplot2")


library(ggplot2)
# Basic histogram
p<-ggplot(sorted_data1, aes(x=logFC)) + 
  geom_histogram(color="black", fill="yellow")
p
#creating the scatter plot
s<-ggplot(sorted_data1, aes(x=logFC, y=P.Value)) +geom_point(shape=1,size = 1, color = "hotpink4")
s




