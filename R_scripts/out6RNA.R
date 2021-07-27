# test file changing

# load required packates
library(dplyr)
library(tidyverse)
library(edgeR)

# set working directory
setwd("~/Documents/out6RNA")

## read in all the data and add the necessary columns
AHQF1A_YB_19 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQF1A_YB_19ReadsPerGene.out.tab", sep="\t",
                         col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQF1A_YB_19 <- mutate(AHQF1A_YB_19, indv = "AHQF1A_YB_19")


AHQN_topM_6 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQN_topM_6ReadsPerGene.out.tab", sep="\t",
                        col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQN_topM_6 <- mutate(AHQN_topM_6, indv = "AHQN_topM_6")


AHQN_topM_2 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQNK_topM_2ReadsPerGene.out.tab", sep="\t",
                        col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQN_topM_2 <- mutate(AHQN_topM_2, indv = "AHQN_topM_2")


AHQN_topM_8 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQNK_topM_8ReadsPerGene.out.tab", sep="\t",
                        col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQN_topM_8 <- mutate(AHQN_topM_8, indv = "AHQN_topM_8")


AHQN_YB_9 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQNK_YB_9ReadsPerGene.out.tab", sep="\t",
                      col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQN_YB_9 <- mutate(AHQN_YB_9, indv = "AHQN_YB_9")


AHQN_YB_12 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQNK_YB_12ReadsPerGene.out.tab", sep="\t",
                       col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQN_YB_12 <- mutate(AHQN_YB_12, indv = "AHQN_YB_12")


AHQT_topM_4 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQT_topM_4ReadsPerGene.out.tab", sep="\t",
                        col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQT_topM_4 <- mutate(AHQT_topM_4, indv = "AHQT_topM_4")


AHQT_topM_5 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQT_topM_5ReadsPerGene.out.tab", sep="\t",
                        col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQT_topM_5 <- mutate(AHQT_topM_5, indv = "AHQT_topM_5")


AHQT_YB_3 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/AHQT_YB_3ReadsPerGene.out.tab", sep="\t",
                      col.names = c("gene", "count", "count_strand1", "count_strand2"))
AHQT_YB_3 <- mutate(AHQT_YB_3, indv = "AHQT_YB_3")


## collate all the individuals
all <- rbind(AHQF1A_YB_19,AHQN_topM_6)
all <- rbind(all, AHQN_topM_2)
all <- rbind(all, AHQN_topM_8)
all <- rbind(all, AHQN_YB_9)
all <- rbind(all, AHQN_YB_12)
all <- rbind(all, AHQT_topM_4)
all <- rbind(all, AHQT_topM_5)
all <- rbind(all, AHQT_YB_3)

# we'll just use the counts
all <- select(all, gene, count, indv)

#spread so columns are individuals and rows are genes
all_wide <- spread(all, indv, count)
all_wide<- filter(all_wide, gene != "N_noFeature")
all_wide<- filter(all_wide, gene != "N_multimapping")
all_wide<- filter(all_wide, gene != "N_ambiguous")
Data <- all_wide
Data <- select(Data, -gene)
# read in the annotation file
Annotation <- read.csv("Mgutt_annotations.csv",header=T)

# read in the Sample Info (treatment, etc)

Sample_Info <- read.csv("AHQ_Sample_Info.csv", header = T)

# make a DGE object
Groups <- Sample_Info$Groups
d <- DGEList(counts=Data,group=factor(Groups), genes=Annotation)

# explore & filter data to remove irrelevant sites
head(d$samples)
head(d$counts)
apply(d$counts, 2, sum)

keep <- rowSums(cpm(d) > 20) >= 2 

# only keep sites w/ at least 20 read counts per million in at least 2 indvs
d.subset <- d[keep,]
dim(d.subset) 

# reset the library sizes
d.subset$samples$lib.size <- colSums(d.subset$counts) 
head(d.subset$samples)

## Normalizing the data 

# explore the variation
barplot(d.subset$samples$lib.size*1e-6, ylab="Library size (millions)", names=rownames(d.subset$samples))  
d.subset$samples

# reset the normalization factors & normalize it
d.subset <- calcNormFactors(d.subset) 
head(d.subset$samples)

d.subset$samples

# explore the data with MDS 
plotMDS(d.subset, method="bcv", col=as.numeric(d.subset$samples$group)) 


## Explore differential expression
d.subset
#make an experimental design matrix
design <- model.matrix(~ 0 + group, data = d.subset$samples)
colnames(design) <- levels(d.subset$samples$group)
design

# estimate dispersion
d.subset <- estimateDisp(d.subset, design) ## Estimate the common and tagwise dispersions in a way that accounts for our experimental design
names(d.subset)  ## notice that the common disperion is now part of the d object

d.subset  ## Here's a look at 
d.subset$design
# exact test for differential expression between pairs

et_NvT <- exactTest(d.subset, pair=c("N_topM","T_topM"))  ## You will change "pair" and the name to the appropriate labels for your groups' questions.
et_NvF1 <- exactTest(d.subset, pair=c("F1","N"))  ## You will change "pair" and the name to the appropriate labels for your groups' questions.
et_TvF1 <- exactTest(d.subset, pair=c("T","F1"))  ## You will change "pair" and the name to the appropriate labels for your groups' questions.

# make tables of the top tags for analysis

Tbl_NvT <- topTags(et_NvT, n=nrow(et_NvT$table))$table  ## Out put the results into a table
write.table(Tbl_NvT, file="Tbl_NvT.csv", sep=",", row.names=TRUE)  ## Save the table for analyses next week!


Tbl_NvF1 <- topTags(et_NvF1, n=nrow(et_NvF1$table))$table  ## Out put the results into a table
write.table(Tbl_NvF1, file="Tbl_NvF1.csv", sep=",", row.names=TRUE)  ## Save the table for analyses next week!


Tbl_TvF1 <- topTags(et_TvF1, n=nrow(et_TvF1$table))$table  ## Out put the results into a table
write.table(Tbl_TvF1, file="Tbl_TvF1.csv", sep=",", row.names=TRUE)  ## Save the table for analyses next week!


## total number of significantly differently expressed genes for each pair: 

is.de.NvT <- decideTestsDGE(et_NvT, p = 0.05)
summary(is.de.NvT)


is.de.NvF1 <- decideTestsDGE(et_NvF1, p = 0.05)
summary(is.de.NvF1)


is.de.TvF1 <- decideTestsDGE(et_TvF1, p = 0.05)
summary(is.de.TvF1)


## plot the differential expression
plotMD(et_NvT)
abline(h=c(-1, 1), col="blue")

plotMD(et_NvF1)
abline(h=c(-1, 1), col="blue")

plotMD(et_TvF1)
abline(h=c(-1, 1), col="blue")


## make some heatmaps

library(pheatmap)
library(RColorBrewer)


# N vs T heatmap
logcpm <- cpm(d.subset, log=TRUE, prior.count = 2)  ## Calculate the counts per million
is.de.logcpm.NvT <- logcpm[is.de.NvT != 0, ]  ## Just keep the differentially expressed genes (i.e, not = 0)
pheatmap(is.de.logcpm.NvT, cluster_rows=TRUE, cutree_rows = 4)
pheatmap(is.de.logcpm.NvT, cluster_rows=TRUE, cutree_rows = 4, scale = "row")

# N vs F1 heatmap
logcpm <- cpm(d.subset, log=TRUE, prior.count = 2)  ## Calculate the counts per million
is.de.logcpm.NvF1 <- logcpm[is.de.NvF1 != 0, ]  ## Just keep the differentially expressed genes (i.e, not = 0)

pheatmap(is.de.logcpm.NvF1, cluster_rows=TRUE, cutree_rows = 4)
pheatmap(is.de.logcpm.NvF1, cluster_rows=TRUE, cutree_rows = 4, scale = "row")

# T vs F1 heatmap
logcpm <- cpm(d.subset, log=TRUE, prior.count = 2)  ## Calculate the counts per million
is.de.logcpm.TvF1 <- logcpm[is.de.TvF1 != 0, ]  ## Just keep the differentially expressed genes (i.e, not = 0)

pheatmap(is.de.logcpm.TvF1, cluster_rows=TRUE, cutree_rows = 4)
pheatmap(is.de.logcpm.TvF1, cluster_rows=TRUE, cutree_rows = 4, scale = "row")

## further analyses
Tbl_NvT 
Tbl_NvF1 
Tbl_TvF1 


# filter to only keep genes with FDR < .05
NvT.DE <-subset(Tbl_NvT, Tbl_NvT$FDR < 0.05)
NvF1.DE <-subset(Tbl_NvF1, Tbl_NvF1$FDR < 0.05)
TvF1.DE <-subset(Tbl_TvF1, Tbl_TvF1$FDR < 0.05)

dim(NvT.DE)
dim(NvF1.DE)
dim(TvF1.DE)

NvT.up <- subset(NvT.DE, NvT.DE$logFC > 0)
NvT.down <- subset(NvT.DE, NvT.DE$logFC < 0)

write.csv(NvT.DE, "NvT_DE.csv")


Tbl_NvT.DE <- mutate(Tbl_NvT, log10p= -log10(PValue))
sig <- Tbl_NvT.DE$FDR < .05
Tbl_NvT.DE <- mutate(Tbl_NvT.DE, sig = sig)

ggplot(Tbl_NvT.DE, aes(x=logFC, y=log10p, col=sig)) + geom_point() + theme_classic() + theme(legend.position = "none") + ggtitle("N vs T differential expression, topM")

write.csv(Tbl_NvT.DE, "Tbl_NvT.DE.csv")

keep <- rowSums(cpm(d) > 20) >= 2 
