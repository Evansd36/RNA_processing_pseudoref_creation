# test file changing

# load required packates
library(dplyr)
library(tidyverse)
library(edgeR)

# set working directory
setwd("~/Documents/out6RNA")

## read in all the data and add the necessary columns
calyx_115_01_4 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/115_calyx01_4ReadsPerGene.out.tab", sep="\t",
                         col.names = c("gene", "count", "count_strand1", "count_strand2"))
calyx_115_01_4 <- mutate(calyx_115_01_4, indv = "115_calyx01_4")


calyx_627c_00_2 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/627c_calyx_00_2ReadsPerGene.out.tab", sep="\t",
                           col.names = c("gene", "count", "count_strand1", "count_strand2"))
calyx_627c_00_2 <- mutate(calyx_627c_00_2, indv = "calyx_627c_00_2")


calyx_664_11_6 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/664_calyx11_6ReadsPerGene.out.tab", sep="\t",
                           col.names = c("gene", "count", "count_strand1", "count_strand2"))
calyx_664_11_6 <- mutate(calyx_664_11_6, indv = "calyx_664_11_6")


Lf2_767_4 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/767_2ndLf_4ReadsPerGene.out.tab", sep="\t",
                           col.names = c("gene", "count", "count_strand1", "count_strand2"))
Lf2_767_4 <- mutate(Lf2_767_4, indv = "Lf2_767_4")


calyx_767_3_00 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/767_calyx_3ReadsPerGene.out.tab", sep="\t",
                      col.names = c("gene", "count", "count_strand1", "count_strand2"))
calyx_767_3_00 <- mutate(calyx_767_3_00, indv = "calyx_767_3_00")


calyx_767_T <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/767_calyx_TReadsPerGene.out.tab", sep="\t",
                        col.names = c("gene", "count", "count_strand1", "count_strand2"))
calyx_767_T <- mutate(calyx_767_T, indv = "calyx_767_T_00")


calyx_1177_3 <- read.csv("~/Dropbox/RNA2020/RNASeq0421/sample data/1177_calyx_3ReadsPerGene.out.tab", sep="\t",
                            col.names = c("gene", "count", "count_strand1", "count_strand2"))
calyx_1177_3 <- mutate(calyx_1177_3 , indv = "calyx_1177_3_00")

## collate all the individuals
all <- rbind(calyx_115_01_4,calyx_627c_00_2)
all <- rbind(all, calyx_664_11_6)
all <- rbind(all, calyx_1177_3)

# we'll just use the counts
all <- select(all, gene, count, indv)

#spread so columns are individuals and rows are genes
all_wide <- spread(all, indv, count)
all_wide<- filter(all_wide, gene != "N_noFeature")
all_wide<- filter(all_wide, gene != "N_multimapping")
all_wide<- filter(all_wide, gene != "N_ambiguous")
Data <- all_wide
Data <- select(Data, -gene)

write.csv(all_wide, "all_wide1.csv")
# read in the annotation file
Annotation <- read.csv("data/Mgutt_annotations.csv",header=T)

# read in the Sample Info (treatment, etc)
Sample_Info <- read.csv("IM_sample1.csv")

Sample_Info_1 <- filter(Sample_Info, Sample_ID != "Lf2_767_4")


# make a DGE object
Groups <- Sample_Info$trichYN
d <- DGEList(counts=Data,group=factor(Groups), genes=Annotation)

dim(Data)
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



# estimate dispersion
d.subset <- estimateDisp(d.subset, design) ## Estimate the common and tagwise dispersions in a way that accounts for our experimental design
names(d.subset)  ## notice that the common disperion is now part of the d object

d.subset  ## Here's a look at 
d.subset$design

d.subset$counts
# exact test for differential expression between pairs
d.subset$design
et_NvT <- exactTest(d.subset, pair=c("Y","N"))  ## You will change "pair" and the name to the appropriate labels for your groups' questions.

Tbl_NvT <- topTags(et_NvT, n=nrow(et_NvT$table))$table  ## Out put the results into a table
write.table(Tbl_NvT, file="Tbl_NvT.csv", sep=",", row.names=TRUE)  ## Save the table for analyses next week!

is.de.NvT <- decideTestsDGE(et_NvT, p = 0.05)
summary(is.de.NvT)

## add -log10p to make volcano plot
Tbl_NvT.DE <- mutate(Tbl_NvT, log10p= -log10(PValue))
sig <- Tbl_NvT.DE$FDR < .05
Tbl_NvT.DE <- mutate(Tbl_NvT.DE, sig = sig)


## code for making a volcano plot
ggplot(Tbl_NvT.DE, aes(x=logFC, y=log10p, col=sig)) + geom_point() + theme_classic() + theme(legend.position = "none") + ggtitle("Trich vs no trich, topM")


write.csv(Tbl_NvT.DE, "trich_noTrich.csv")
