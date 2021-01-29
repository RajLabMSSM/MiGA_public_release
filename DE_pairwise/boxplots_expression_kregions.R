# Boxplots with expression 
# August 10, 2020 

#This command clean all variables. BE CAREFULL!!! 
rm(list = setdiff(ls(), lsf.str()))

library(tidyverse)
library(edgeR)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(ggeasy)
library(ggsci)

expression_dir = "~/ad-omics_hydra/microglia_omics/expression_tables/added_pilot_314s/expr_4brain_regions/"
work_plots = "~/pd-omics/katia/Microglia/mic_255s/figures_pdf/boxplots_expres/"
metadata_path = "~/pd-omics/katia/Microglia/mic_255s/metadata_files/"
gene_lists = "~/pd-omics/katia/scripts/GitHub_scripts/glia_omics/3rd_pass_mic_255s/brain_regions/pairwise_dream/"

## Pairwise comparison lists:
MFGxSTG = read.table(paste0(gene_lists, "MFGxSTG_dream_3rd.txt"), header = T, stringsAsFactors = F)
SVZxSTG = read.table(paste0(gene_lists, "SVZxSTG_dream_3rd.txt"), header = T, stringsAsFactors = F)
THAxSVZ = read.table(paste0(gene_lists, "THAxSVZ_dream_3rd.txt"), header = T, stringsAsFactors = F)
MFGxSVZ = read.table(paste0(gene_lists, "MFGxSVZ_dream_3rd.txt"), header = T, stringsAsFactors = F)
THAxSTG = read.table(paste0(gene_lists, "THAxSTG_dream_3rd.txt"), header = T, stringsAsFactors = F)
THAxMFG = read.table(paste0(gene_lists, "THAxMFG_dream_3rd.txt"), header = T, stringsAsFactors = F)

## Expression data
load(paste0(expression_dir, "Expression_filt_255s.Rdata")) 
dim(genes_counts_exp_3rd) 
dim(genes_counts_voom_3rd)

load(paste0(metadata_path, "metadata_255filt_eng_29jun2020.Rdata"))

## Match names
genes_voom <- as.data.frame(genes_counts_voom_3rd)

# To include cause of death and C1-C4
metadata3rd_pass$cause_of_death_categories[metadata3rd_pass$cause_of_death_categories %in% NA] <- "Other"
# table(metadata3rd_pass$cause_of_death_categories)

metadata3rd_pass$C1 = metadata3rd_pass$C1 %>% replace_na(median(metadata3rd_pass$C1, na.rm = T))
metadata3rd_pass$C2 = metadata3rd_pass$C2 %>% replace_na(median(metadata3rd_pass$C2, na.rm = T))
metadata3rd_pass$C3 = metadata3rd_pass$C3 %>% replace_na(median(metadata3rd_pass$C3, na.rm = T))
metadata3rd_pass$C4 = metadata3rd_pass$C4 %>% replace_na(median(metadata3rd_pass$C4, na.rm = T))

## Get the residuals: 
metadata4deg = metadata3rd_pass

res.pca = prcomp(t(genes_voom)) 
autoplot(res.pca, data = metadata4deg, colour = 'age') +
  scale_colour_viridis_c() +
  theme_classic()

allResiduals <- removeBatchEffect(x = genes_voom, 
                                  batch = metadata4deg$sex, 
                                  batch2 = metadata4deg$cause_of_death_categories,
                                  # design = model.matrix(~ age, data = metadata4deg), #force to not regress age. It didn't change the PCA. 
                                  covariates = as.matrix(metadata4deg[, c("picard_pct_mrna_bases", "picard_summed_median", "picard_pct_ribosomal_bases", "C1","C2","C3","C4" )]))

res.pca = prcomp(t(allResiduals)) 

autoplot(res.pca, data = metadata4deg, colour = 'age') +
  scale_colour_viridis_c() +
  theme_classic()

## Get conversion table for Gencode 30
gencode_30 = read.table("~/pd-omics/katia/ens.geneid.gencode.v30")
colnames(gencode_30) = c("ensembl","symbol")
allResiduals = as.data.frame(allResiduals)
allResiduals$ensembl = rownames(allResiduals)
allResiduals = merge(allResiduals, gencode_30, by = "ensembl")
dim(allResiduals)
head(allResiduals)

## Includes the expression in the metadata 
head(metadata4deg)
rownames(metadata4deg) = metadata4deg$donor_tissue

all(rownames(metadata4deg) == colnames(allResiduals)[2:256]) # skip the ensembl and symbol columns 

genes = c("CD83", "FCGR3B", "MRC1", "RGS1")
  
#########Boxplots with Residuals data 

for (gene in genes){
  metadata4deg$gene_1_expr = as.numeric(allResiduals[allResiduals$symbol %in% gene, as.character(metadata4deg$donor_tissue)])
  
  stat.test <- tibble::tribble(
    ~group1, ~group2,   ~p.adj,
    "MFG",     "STG", MFGxSTG$adj.P.Val[MFGxSTG$symbol == gene],
    "SVZ",     "STG", SVZxSTG$adj.P.Val[SVZxSTG$symbol == gene],
    "THA",     "SVZ", THAxSVZ$adj.P.Val[THAxSVZ$symbol == gene],
    "MFG",     "SVZ", MFGxSVZ$adj.P.Val[MFGxSVZ$symbol == gene],
    "THA",     "STG", THAxSTG$adj.P.Val[THAxSTG$symbol == gene],
    "THA",     "MFG", THAxMFG$adj.P.Val[THAxMFG$symbol == gene]
  )
  stat.test$p.adj = signif(stat.test$p.adj,2)
  #stat.test = stat.test[stat.test$p.adj<0.05,] # Only shows the significative pvalues 
  
  #pdf(file = paste0(work_plots, "pvalue_dream/boxplots_residuals/boxplot_exp_", gene ,".pdf"), width = 4, height = 4, onefile = F)  
  ggboxplot(metadata4deg, x = "tissue", y = "gene_1_expr", fill = "tissue", add = "jitter", outlier.shape = NA) +
    stat_pvalue_manual(stat.test, y.position = max(metadata4deg$gene_1_expr), step.increase = 0.1, label = "p.adj") +
    scale_fill_futurama() +
    easy_labs(x = "Region", y = bquote(italic(.(gene))~' residual expression')) +
    theme_classic() +
    easy_remove_legend()
  
  # ggsave(file = paste0(work_plots, "pvalue_dream/boxplots_residuals/boxplot_exp_", gene ,".pdf"), width = 4, height = 4)
  #dev.off()  
}




