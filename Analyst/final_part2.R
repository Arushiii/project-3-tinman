
setwd("/projectnb/bf528/users/tinman/Project3/Analyst")

#data of rna_seq
deseq <- read.csv('/projectnb/bf528/users/tinman/Project3/programmer/deseq/dna_deseq.csv')

#data of microarray
limma <- read.csv('/projectnb/bf528/users/tinman/Project3/Analyst/N-NITROSODIETHYLAMINE_limma_results.csv')

#mapping file
base <- read.csv('/project/bf528/project_3/refseq_affy_map.csv')
#removed the na values
base<-na.omit(base)


#filtered the RNA-seq and Microarray dataset
#filtering out the genes based on adjusted p-value on deseq
DE_deseq<- deseq[deseq$padj <0.05,]

#filtering out the genes based on adjusted p-value  and logFC on limma
DE_limma<- limma[limma$adj.P.Val <0.05,]

#changing the column name
names(DE_deseq)[names(DE_deseq) == "X"] <- "REFSEQ"
names(DE_limma)[names(DE_limma)=='X'] <- 'PROBEID'

#merge the filtered DE_deseq with the base

mix <- merge(DE_deseq,base,c('REFSEQ'))
mix_f<- merge(DE_limma,mix,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
mix_f <- mix_f %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
mix_f <- mix_f %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

mix_f <- mix_f[mix_f$Regulation_deseq == mix_f$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n1<- nrow(DE_deseq)


#n2= number of gene in filtered limma
n2<- nrow(DE_limma)

# n0 represents the number of genes that are common between both methods
n0<- nrow(mix_f)

#background corrected intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2) 
#concordance calculation
c = round((2*(nx))/(n1+n2)*100) 

#xxxxx-------------------------------------------------------------------------------------xxxx


setwd("/projectnb/bf528/users/tinman/Project3/Analyst")

#data of rna_seq
deseq1 <- read.csv('/projectnb/bf528/users/tinman/Project3/programmer/deseq/er_deseq.csv')

#data of microarray
limma1 <- read.csv('/projectnb/bf528/users/tinman/Project3/Analyst/BETA-ESTRADIOL_limma_results.csv')

#mapping file
base1 <- read.csv('/project/bf528/project_3/refseq_affy_map.csv')
#removed the na values
base1<-na.omit(base1)


#filtered the RNA-seq and Microarray dataset
#filtering out the genes based on adjusted p-value on deseq
DE_deseq1<- deseq1[deseq1$padj <0.05,]

#filtering out the genes based on adjusted p-value  and logFC on limma
DE_limma1<- limma1[limma1$adj.P.Val <0.05,]

#changing the column name
names(DE_deseq1)[names(DE_deseq1) == "X"] <- "REFSEQ"
names(DE_limma1)[names(DE_limma1)=='X'] <- 'PROBEID'

#merge the filtered DE_deseq with the base

mix1 <- merge(DE_deseq1,base1,c('REFSEQ'))
mix_f1<- merge(DE_limma1,mix1,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
mix_f1 <- mix_f1 %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
mix_f1 <- mix_f1 %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

mix_f1 <- mix_f1[mix_f1$Regulation_deseq == mix_f1$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n11<- nrow(DE_deseq1)


#n2= number of gene in filtered limma
n21<- nrow(DE_limma1)

# n0 represents the number of genes that are common between both methods
n01<- nrow(mix_f1)

#background corrected intersection
nx1 <- ((n01*N) - (n11*n21))/(n01+N-n11-n21) 
#concordance calculation
c1 = round((2*(nx1))/(n11+n21)*100 )

#xxxxx-------------------------------------------------------------------------------------xxxx



setwd("/projectnb/bf528/users/tinman/Project3/Analyst")

#data of rna_seq
deseq2 <- read.csv('/projectnb/bf528/users/tinman/Project3/programmer/deseq/ppara_deseq.csv')

#data of microarray
limma2 <- read.csv('/projectnb/bf528/users/tinman/Project3/Analyst/BEZAFIBRATE_limma_results.csv')

#mapping file
base2 <- read.csv('/project/bf528/project_3/refseq_affy_map.csv')
#removed the na values
base2<-na.omit(base2)


#filtered the RNA-seq and Microarray dataset
#filtering out the genes based on adjusted p-value on deseq
DE_deseq2<- deseq2[deseq2$padj <0.05,]

#filtering out the genes based on adjusted p-value  and logFC on limma
DE_limma2<- limma2[limma2$adj.P.Val <0.05,]

#changing the column name
names(DE_deseq2)[names(DE_deseq2) == "X"] <- "REFSEQ"
names(DE_limma2)[names(DE_limma2)=='X'] <- 'PROBEID'

#merge the filtered DE_deseq with the base

mix2 <- merge(DE_deseq2,base2,c('REFSEQ'))
mix_f2<- merge(DE_limma2,mix2,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
mix_f2 <- mix_f2 %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
mix_f2 <- mix_f2 %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

mix_f2 <- mix_f2[mix_f2$Regulation_deseq == mix_f2$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n12<- nrow(DE_deseq2)


#n2= number of gene in filtered limma
n22<- nrow(DE_limma2)

# n0 represents the number of genes that are common between both methods
n02<- nrow(mix_f2)

#background corrected intersection
nx2 <- ((n02*N) - (n12*n22))/(n02+N-n12-n22) 
#concordance calculation
c2 = round((2*(nx2))/(n12+n22)*100) 

#xxxxx-------------------------------above------------------------------------------------------xxxx

#for deseq
m<- round(median(DE_deseq$baseMean,na.rm = T),digits=2)
above_d <- subset(DE_deseq, DE_deseq$baseMean >= m) #n1

#for miroarray

mm <- round(median(DE_limma$AveExpr, na.rm=T), digits=2)
above_m <- subset(DE_limma, DE_limma$AveExpr >= mm)  #n2

#merge


ab <- merge(above_d,base,c('REFSEQ'))
abo<- merge(above_m,ab,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
abo <- abo %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
abo <- abo %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

abo <- abo[abo$Regulation_deseq == abo$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n1ad<- nrow(above_d)


#n2= number of gene in filtered limma
n2am<- nrow(above_m)

# n0 represents the number of genes that are common between both methods
n0a<- nrow(abo)

#background corrected intersection
nxa <- ((n0a*N) - (n1ad*n2am))/(n0a+N-n1ad-n2am) 
#concordance calculation
c_above = round((2*(nxa))/(n1ad+n2am)*100) 

#...........below...............................................................#


m<- round(median(DE_deseq$baseMean,na.rm = T),digits=2)
below_d <- subset(DE_deseq, DE_deseq$baseMean < m) #n1

#for miroarray

mm <- round(median(DE_limma$AveExpr, na.rm=T), digits=2)
below_m <- subset(DE_limma, DE_limma$AveExpr < mm)  #n2

#merge


b <- merge(below_d,base,c('REFSEQ'))
bo<- merge(below_m,b,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
bo <- bo %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
bo <- bo %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

bo <- bo[bo$Regulation_deseq == bo$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n1bd<- nrow(below_d)


#n2= number of gene in filtered limma
n2bm<- nrow(below_m)

# n0 represents the number of genes that are common between both methods
n0b<- nrow(bo)

#background corrected intersection
nxb <- ((n0b*N) - (n1bd*n2bm))/(n0b+N-n1bd-n2bm) 
#concordance calculation
c_below = round((2*(nxb))/(n1bd+n2bm)*100) 

#==================================================================================#
#==================================================================================#


#for deseq
m1<- round(median(DE_deseq1$baseMean,na.rm = T),digits=2)
above_d1 <- subset(DE_deseq1, DE_deseq1$baseMean >= m1) #n1

#for miroarray

mm1 <- round(median(DE_limma1$AveExpr, na.rm=T), digits=2)
above_m1 <- subset(DE_limma1, DE_limma1$AveExpr >= mm1)  #n2

#merge


ab1 <- merge(above_d1,base,c('REFSEQ'))
abo1<- merge(above_m1,ab1,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
abo1 <- abo1 %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
abo1 <- abo1 %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

abo1 <- abo1[abo1$Regulation_deseq == abo1$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n1ad1<- nrow(above_d1)


#n2= number of gene in filtered limma
n2am1<- nrow(above_m1)

# n0 represents the number of genes that are common between both methods
n0a1<- nrow(abo1)

#background corrected intersection
nxa1 <- ((n0a1*N) - (n1ad1*n2am1))/(n0a1+N-n1ad1-n2am1) 
#concordance calculation
c_above1 = round((2*(nxa1))/(n1ad1+n2am1)*100) 

#...........below............................


m1<- round(median(DE_deseq1$baseMean,na.rm = T),digits=2)
below_d1 <- subset(DE_deseq1, DE_deseq1$baseMean < m1) #n1

#for miroarray

mm1 <- round(median(DE_limma1$AveExpr, na.rm=T), digits=2)
below_m1 <- subset(DE_limma1, DE_limma1$AveExpr < mm1)  #n2

#merge


b1 <- merge(below_d1,base,c('REFSEQ'))
bo1<- merge(below_m1,b1,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
bo1 <- bo1 %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
bo1 <- bo1 %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

bo1 <- bo1[bo1$Regulation_deseq == bo1$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n1bd1<- nrow(below_d1)


#n2= number of gene in filtered limma
n2bm1<- nrow(below_m1)

# n0 represents the number of genes that are common between both methods
n0b1<- nrow(bo1)

#background corrected intersection
nxb1 <- ((n0b1*N) - (n1bd1*n2bm1))/(n0b1+N-n1bd1-n2bm1) 
#concordance calculation
c_below1 = round((2*(nxb1))/(n1bd1+n2bm1)*100) 

#==================================================================================#
#==================================================================================#

#for deseq
m2<- round(median(DE_deseq2$baseMean,na.rm = T),digits=2)
above_d2 <- subset(DE_deseq2, DE_deseq2$baseMean >= m2) #n1

#for miroarray

mm2 <- round(median(DE_limma2$AveExpr, na.rm=T), digits=2)
above_m2 <- subset(DE_limma2, DE_limma2$AveExpr >= mm2)  #n2

#merge


ab2 <- merge(above_d2,base,c('REFSEQ'))
abo2<- merge(above_m2,ab2,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
abo2 <- abo2 %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
abo2 <- abo2 %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

abo2 <- abo2[abo2$Regulation_deseq == abo2$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n1ad2<- nrow(above_d2)


#n2= number of gene in filtered limma
n2am2<- nrow(above_m2)

# n0 represents the number of genes that are common between both methods
n0a2<- nrow(abo2)

#background corrected intersection
nxa2 <- ((n0a2*N) - (n1ad2*n2am2))/(n0a2+N-n1ad2-n2am2) 
#concordance calculation
c_above2 = round((2*(nxa2))/(n1ad2+n2am2)*100) 

#...........below............................


m2<- round(median(DE_deseq2$baseMean,na.rm = T),digits=2)
below_d2 <- subset(DE_deseq2, DE_deseq2$baseMean < m2) #n1

#for miroarray

mm2 <- round(median(DE_limma2$AveExpr, na.rm=T), digits=2)
below_m2 <- subset(DE_limma2, DE_limma2$AveExpr < mm2)  #n2

#merge


b2 <- merge(below_d2,base,c('REFSEQ'))
bo2<- merge(below_m2,b2,c('PROBEID'))

#we will add the direction based on the fold change for deseq and microarray
#  log2fc used for directionality in deseq
bo2 <- bo2 %>%
  mutate(Regulation_deseq = if_else(log2FoldChange < 0, '+', '-' ))

#  logFC used for directionality in limma
bo2 <- bo2 %>%
  mutate(Regulation_limma = if_else(logFC < 0, '+', '-' ))

#now filtered all the rows with different direction

bo2 <- bo2[bo2$Regulation_deseq == bo2$Regulation_limma ,]


#concordance calculation

#The sequence, revealed in Nature1, has about 25,000 Total genes in rat genome.
#thus N is the total number of gene in Rat
N =25000

#n1= number of gene in filtered deseq
n1bd2<- nrow(below_d2)


#n2= number of gene in filtered limma
n2bm2<- nrow(below_m2)

# n0 represents the number of genes that are common between both methods
n0b2<- nrow(bo2)

#background corrected intersection
nxb2 <- ((n0b2*N) - (n1bd2*n2bm2))/(n0b2+N-n1bd2-n2bm2) 
#concordance calculation
c_below2 = round((2*(nxb2))/(n1bd2+n2bm2)*100) 

#==================================================================================#
#==================================================================================#


#plots
chemical1<-c(c_above,c,c_below)
chemical2<-c(c_above1,c1,c_below1)
chemical3<-c(c_above2,c2,c_below2)

conc<-c(chemical1,chemical2,chemical3)

chem<-rep(c("N-NITROSODIMETHYLAMINE","BETA-ESTRADIOL","BEZAFIBRATE" ), each= 3)
group<-rep(c('above','overall','below'),3)

df<- data.frame(conc,chem,group)

g<- ggplot(data=df,aes(x=chem,y=conc,fill=group)) + geom_bar(stat="identity", color="hotpink",position=position_dodge())+ ylim(0,51)+scale_fill_manual("legend", values = c("above" = "#F066EA", "below" = "#BF80FF", "overall" = "#00BDD0"))

g

#======================================================================================================================================================================================================================================================
#plot for overall concordance with no.of deseq genes
conc<-c(c,c1,c2)

group<-rep(c("DNA","ER","PPARA" ))
genes<- c(n1,n11,n12)

df<- data.frame(group,conc,genes)

g<- ggplot(data=df,aes(x=genes,y=conc))+geom_point()+ geom_text(label=group)+ylim(0,100)+xlim(1000,5000)
g

#============================================================================================================================
#plot for overall concordance with no.of deseq genes
conc<-c(c,c1,c2)

group<-rep(c("DNA","ER","PPARA" ))
genes<- c(n2,n21,n22)

df<- data.frame(group,conc,genes)

g<- ggplot(data=df,aes(x=genes,y=conc))+geom_point()+ geom_text(label=group)+ylim(0,100)
g

