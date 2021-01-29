library(DESeq2)
library(edgeR)
library(dplyr)
library(data.table)

setwd("C:\\Users\\Gebruiker\\Documents\\Neuroscience and cognition\\Minor Research Project\\RNA sequencing microglia\\differential expression analysis")

#read in covariate file and count file
pheno <-read.table("covariates_24samples.txt", header=T, row.names = 1)
exp <-read.table("061019_txi.rsem_genes_counts.txt", header=T, row.names = 1)

#Make sure that colnames in exp file match rownames in covariate file
colnames(exp) <- gsub("_S.*", "", colnames(exp))
exp <- exp[, rownames(pheno)]
identical(rownames(pheno), colnames(exp))

colnames(exp) %in% rownames(pheno)
dim(exp)

#select the condition to compare
pheno_LPS <- subset(pheno, stimulation == 0 | stimulation == 1)
exp_LPS <- exp[, rownames(pheno_LPS)]
identical(rownames(pheno_LPS), colnames(exp_LPS))

pheno_TGFb <- subset(pheno, stimulation ==0 | stimulation == 2)
exp_TGFb <- exp[, rownames(pheno_TGFb)]
identical(rownames(pheno_TGFb), colnames(exp_TGFb))

pheno_IFNy <- subset(pheno, stimulation ==0 | stimulation == 3)
exp_IFNy <- exp[, rownames(pheno_IFNy)]
identical(rownames(pheno_IFNy), colnames(exp_IFNy))

#without unstimulated outlier
pheno_TGFb_2 <- pheno_TGFb[c(1:8, 10:12),]
exp_TGFb_2 <- exp[, rownames(pheno_TGFb_2)]
identical(rownames(pheno_TGFb_2), colnames(exp_TGFb_2))

pheno_IFNy_2 <- pheno_IFNy[c(1:8, 10:12),]
exp_IFNy_2 <- exp[, rownames(pheno_IFNy_2)]
identical(rownames(pheno_IFNy_2), colnames(exp_IFNy_2))

#without outliers ustimulated & LPS (A3 & B3)
pheno_LPS_2 <- pheno_LPS[c(1:8, 11:12),]
exp_LPS_2 <- exp[, rownames(pheno_LPS_2)]
identical(rownames(pheno_LPS_2), colnames(exp_LPS_2))

#without outliers ustimulated & LPS (A3 & B3) and empty sample (F2)
pheno_LPS_3 <- pheno_LPS_2[c(1:7, 9:10),]
exp_LPS_3 <- exp[, rownames(pheno_LPS_3)]
identical(rownames(pheno_LPS_3), colnames(exp_LPS_3))

#round counts; deseq2 can only handle integers
exp_LPS <- round(exp_LPS, digits=0)
exp_TGFb <- round(exp_TGFb, digits=0)
exp_IFNy <- round(exp_IFNy, digits=0)

exp_LPS_2 <- round(exp_LPS_2, digits=0)
exp_TGFb_2 <- round(exp_TGFb_2, digits=0)
exp_IFNy_2 <- round(exp_IFNy_2, digits=0)

exp_LPS_3 <- round(exp_LPS_3, digits=0)

#make sure covariate variables are the right type
pheno_LPS_2$donor <- as.factor(pheno_LPS_2$donor)
pheno_LPS_2$stimulation <- as.factor(pheno_LPS_2$stimulation)
pheno_LPS_2$batch <- as.factor(pheno_LPS_2$batch)
pheno_LPS_2$gender <- as.factor(pheno_LPS_2$gender)
pheno_LPS_2$PMI <- as.integer(pheno_LPS_2$PMI)        #cannot create dds object with numeric values
pheno_LPS_2$age <- as.integer(pheno_LPS_2$age)

#Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = exp_LPS_2,
                              colData = pheno_LPS_2,
                              design = ~ donor + stimulation) #variable of interest at end of the formula

#Make sure that control group is set as the reference group
dds$stimulation <- relevel(dds$stimulation, ref="0")
summary(dds$stimulation)
head(dds)

#filter
keep.exp = rowSums(cpm(exp_LPS_3) > 1) >= 0.3*ncol(exp_LPS_3)
dds = dds[keep.exp,]

#Run Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds)
sum(res$padj < 0.05, na.rm=TRUE)
resOrdered <- res[order(res$pvalue),] 
resOrdered <- as.data.frame(resOrdered)

#add gene names
resOrdered <- read.table("LPS\\061219_genes_deseq2_1in20_donor_stimulation_LPS_all.csv", header=T)
setDT(resOrdered, keep.rownames = "GENEID")
resOrdered <- as.data.frame(resOrdered)
gene_mapping <- read.table("C:\\Users\\Gebruiker\\Documents\\Neuroscience and cognition\\Minor Research Project\\Analyse in R\\Data files\\gene_mapping_Jack.txt", header=T)
resOrdered <- left_join(resOrdered, gene_mapping, by = "GENEID")
write.table(resOrdered, file = "061219_genes_deseq2_1in30_donor_stimulation_LPS_wo_outliers_empty.csv")

#PCA
vsd <- vst(dds, blind=FALSE) 
plotPCA(vsd, intgroup = "stimulation")
dev.off()
