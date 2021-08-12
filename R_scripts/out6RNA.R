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
#all <- rbind(AHQF1A_YB_19,AHQN_topM_6)
all <- rbind(AHQN_topM_2,AHQN_topM_6)
all <- rbind(all, AHQN_topM_8)
#all <- rbind(all, AHQN_YB_9)
#all <- rbind(all, AHQN_YB_12)
all <- rbind(all, AHQT_topM_4)
all <- rbind(all, AHQT_topM_5)
#all <- rbind(all, AHQT_YB_3)

head(all)

# we'll just use the counts
all <- select(all, gene, count, indv)
write.csv(all, "all.csv")

AHQN_topM_2

#spread so columns are individuals and rows are genes
all_wide <- spread(all, indv, count)
all_wide<- filter(all_wide, gene != "N_noFeature")
all_wide<- filter(all_wide, gene != "N_multimapping")
all_wide<- filter(all_wide, gene != "N_ambiguous")


Data <- all_wide
Data <- select(Data, -gene)

View(Data)

d$genes
all_wide

write.csv(all_wide, "all_wide.csv")
# read in the annotation file
Annotation <- read.csv("data/Mgutt_annotations.csv",header=T)
View(Annotation)
# read in the Sample Info (treatment, etc)

Sample_Info <- read.csv("data/AHQ_Sample_Info4.csv", header = T)
write.csv(Data, "Data.csv")
# make a DGE object
Groups <- Sample_Info$Groups
d <- DGEList(counts=Data,group=factor(Groups), genes=Annotation)

View(d$genes)

write.csv((cbind(d$counts, d$genes)), "test_new.csv")

# explore & filter data to remove irrelevant sites
head(d$samples)
head(d$counts)
apply(d$counts, 2, sum)

d$samples

cpm(d)

keep <- rowSums(cpm(d) > 20) >= 4
apply(d.subset$counts, 2, sum) 
# only keep sites w/ at least 20 read counts per million in at least 2 indvs
d.subset <- d[keep,]
dim(d.subset) 

# reset the library sizes
d.subset$samples$lib.size <- colSums(d.subset$counts) 
d.subset$samples

## Normalizing the data 

# explore the variation
barplot(d.subset$samples$lib.size*1e-6, ylab="Library size (millions)", names=rownames(d.subset$samples))  
d.subset$samples


# reset the normalization factors & normalize it
d.subset <- calcNormFactors(d.subset) 

# explore the data with MDS 
plotMDS(d.subset, method="bcv", col=as.numeric(d.subset$samples$group)) 


## Explore differential expression
#make an experimental design matrix
design <- model.matrix(~ 0 + group, data = d.subset$samples)
colnames(design) <- levels(d.subset$samples$group)

write.csv(d.subset$genes, "genes.csv")
# estimate dispersion
d.subset <- estimateDisp(d.subset, design)

# exact test for differential expression between pairs

et_NvT <- exactTest(d.subset, pair=c("N_topM","T_topM"))  ## You will change "pair" and the name to the appropriate labels for your groups' questions.

# make tables of the top tags for analysis
# the problem happens before this part
Tbl_NvT <- topTags(et_NvT, n=nrow(et_NvT$table))$table  ## Out put the results into a table
write.table(Tbl_NvT, file="Tbl_NvT.csv", sep=",", row.names=TRUE)  ## Save the table for analyses next week!

## total number of significantly differently expressed genes for each pair: 

is.de.NvT <- decideTestsDGE(et_NvT, p = 0.05)
summary(is.de.NvT)


## plot the differential expression
plotMD(et_NvT)
abline(h=c(-1, 1), col="blue")


# filter to only keep genes with FDR < .05
NvT.DE <-subset(Tbl_NvT, Tbl_NvT$FDR < 0.05)

dim(NvT.DE)

NvT.up <- subset(NvT.DE, NvT.DE$logFC > 0)
NvT.down <- subset(NvT.DE, NvT.DE$logFC < 0)

write.csv(NvT.DE, "NvT_DE.csv")


Tbl_NvT.DE <- mutate(Tbl_NvT, log10p= -log10(PValue))
sig <- Tbl_NvT.DE$FDR < .05
Tbl_NvT.DE <- mutate(Tbl_NvT.DE, sig = sig)

ggplot(Tbl_NvT.DE, aes(x=logFC, y=log10p, col=sig)) + geom_point() + theme_classic() + theme(legend.position = "none") + ggtitle("N vs T differential expression, topM")

write.csv(Tbl_NvT.DE, "Tbl_NvT.DE2.csv")

keep <- rowSums(cpm(d) > 20) >= 2 
