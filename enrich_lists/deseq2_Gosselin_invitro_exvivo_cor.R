library(DESeq2)
library(edgeR)
library(data.table)

setwd("~/Documents/Genelists")
pheno <- fread("Gosselin2017_tableS1_patientinfo.txt", header=T)
GosselinCounts <- fread("Gosselin2017_tableS2_readcounts.txt", header=T)
GosselinCounts <- as.data.frame(GosselinCounts)
rownames(GosselinCounts)=make.names(GosselinCounts$V1, unique=T)
#GosselinCounts$V1 = NULL
pheno <- as.data.frame(pheno)
rownames(pheno) <- pheno$Sample
#pheno$Sample = NULL
dim(GosselinCounts)
dim(pheno)

#change colnames
colnames(GosselinCounts)
colnames(GosselinCounts) <- gsub("RNA_HMG..._s*","",colnames(GosselinCounts))
colnames(GosselinCounts) <- gsub("_Microglia", "", colnames(GosselinCounts))
colnames(GosselinCounts)

#remove monocyte columns
GosselinCounts_microglia <- GosselinCounts[, -grep("*Monocyte*", colnames(GosselinCounts))]
GosselinCounts_microglia
dim(GosselinCounts_microglia)

#remove lysate colums from count data file
GosselinCounts_microglia$S005_Lysate_ExVivo = NULL
GosselinCounts_microglia$S007_Lysate_ExVivo = NULL
GosselinCounts_microglia$S013_Lysate_ExVivo = NULL
GosselinCounts_microglia$S020_Lysate_ExVivo = NULL
GosselinCounts_microglia$S010_Lysate_ExVivo = NULL
dim(GosselinCounts_microglia)

#remove 6h/10 days and other factors
GosselinCounts_microglia$S029_InVitro_24h = NULL
GosselinCounts_microglia$S006_InVitro_10d = NULL
GosselinCounts_microglia$S006_InVitro_10d_MCSF = NULL
GosselinCounts_microglia$S006_InVitro_10d_MCSF_TGFB = NULL
GosselinCounts_microglia$S009_InVitro_10d_MCSF = NULL
GosselinCounts_microglia$S009_InVitro_10d_MCSF_TGFB = NULL
GosselinCounts_microglia$S014_InVitro_7d_MCSF = NULL
GosselinCounts_microglia$S014_InVitro_7d_MCSF_TGFB = NULL
GosselinCounts_microglia$S014_InVitro_7d_TGFB = NULL
GosselinCounts_microglia$S025_InVitro_6h_TGFBR1antagonist = NULL
GosselinCounts_microglia$S025_InVitro_6h_TGFBR = NULL
GosselinCounts_microglia$S025_InVitro_6h_TGFBR1antagonist = NULL
GosselinCounts_microglia$S031_InVitro_7d_HumanSerum = NULL
GosselinCounts_microglia$S031_InVitro_7d_StemPro = NULL
GosselinCounts_microglia$S031_InVitro_7d_StemPro_HumanSerum = NULL
GosselinCounts_microglia$S031_InVitro_7d_StemPro_HumanSerum.1 = NULL
GosselinCounts_microglia$S033_InVitro_6h = NULL
GosselinCounts_microglia$S036_InVitro_7d_ACM = NULL
GosselinCounts_microglia$S036_InVitro_7d_ACM.1 = NULL
GosselinCounts_microglia$S036_InVitro_7d_ACM_HumanSerum = NULL
dim(GosselinCounts_microglia)

#select rows in pheno file with samples present in count data file
dim(pheno)
pheno1 <- subset(pheno, rownames(pheno) %in% colnames(GosselinCounts_microglia))   #pheno should be as.data.frame
pheno1
dim(pheno1)
dim(GosselinCounts_microglia)

#Make sure that colnames in exp file match rownames in covariate file
pheno1 <- as.data.frame(pheno1)
pheno1 <- pheno1[order(pheno1$Sample),]                            #order rows pheno data
GosselinCounts_microglia$V1 = NULL
GosselinCounts_microglia <- GosselinCounts_microglia[, rownames(pheno1)]
identical(rownames(pheno1), colnames(GosselinCounts_microglia))
colnames(pheno1)[colnames(pheno1)=="ex vivo (0) / in vitro (1)"] <- "InVitro"
dds$InVitro = as.factor(dds$InVitro)                              #specify that 0 and 1 is two different groups, not continuous numerical values

#Create DESeq object 
dds <- DESeqDataSetFromMatrix(countData = GosselinCounts_microglia,
                              colData = pheno1,
                              design = ~ Patient + InVitro) #variable of interest at end of the formula

#Make sure that control group is set as the reference group
dds$InVitro <- relevel(dds$InVitro, ref="0")
summary(dds$InVitro)
head(dds)

#Pull normalized count data#
dds <- estimateSizeFactors(dds)
a <- counts(dds, normalized=TRUE)

#Run Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds)
res
sum(res$padj < 0.01, na.rm=TRUE)
resOrdered <- res[order(res$padj),]
resOrdered
write.table(resOrdered, file = "Gosselin_DGE_invitro.txt", quote = F)
Gosselin_DGE_invitro <- read.table("~/Documents/Genelists/Gosselin_DGE_invitro.txt", header=T, row.names=1)
library(data.table)
setDT(Gosselin_DGE_invitro, keep.rownames = "gene_name")
