# July 20, 2020 
# Katia Lopes

# Small heatmap 

#This command clean all variables. BE CAREFULL!!! 
rm(list = setdiff(ls(), lsf.str()))

# install.packages("wesanderson")

library(venn)
library(ggsci)
library(factoextra)
library(reshape2)
library(ggplot2)
library(Hmisc)
library(stats)
library(RColorBrewer)
library(viridis)
library(ggeasy)

gene_lists = "~/pd-omics/katia/scripts/GitHub_scripts/glia_omics/3rd_pass_mic_255s/brain_regions/pairwise_dream/"
expression_dir = "~/ad-omics_hydra/microglia_omics/expression_tables/added_pilot_314s/expr_4brain_regions/"
work_plots = "/Users/katia/Desktop/new_figures/"
metadata_path = "~/pd-omics/katia/Microglia/mic_255s/metadata_files/"
gene_lists4comp = "~/pd-omics/katia/Microglia/GeneLists_4analysis/"

## Pairwise comparison lists:
MFGxSTG = read.table(paste0(gene_lists, "MFGxSTG_dream_3rd.txt"), header = T, stringsAsFactors = F)
SVZxSTG = read.table(paste0(gene_lists, "SVZxSTG_dream_3rd.txt"), header = T, stringsAsFactors = F)
THAxSVZ = read.table(paste0(gene_lists, "THAxSVZ_dream_3rd.txt"), header = T, stringsAsFactors = F)
MFGxSVZ = read.table(paste0(gene_lists, "MFGxSVZ_dream_3rd.txt"), header = T, stringsAsFactors = F)
THAxSTG = read.table(paste0(gene_lists, "THAxSTG_dream_3rd.txt"), header = T, stringsAsFactors = F)
THAxMFG = read.table(paste0(gene_lists, "THAxMFG_dream_3rd.txt"), header = T, stringsAsFactors = F)

MFGxSTG.filt = MFGxSTG %>% dplyr::filter(adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 1)
SVZxSTG.filt = SVZxSTG %>% dplyr::filter(adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 1)
THAxSVZ.filt = THAxSVZ %>% dplyr::filter(adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 1)
MFGxSVZ.filt = MFGxSVZ %>% dplyr::filter(adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 1)
THAxSTG.filt = THAxSTG %>% dplyr::filter(adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 1)
THAxMFG.filt = THAxMFG %>% dplyr::filter(adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 1)

de_lists = list(MFGxSTG = rownames(MFGxSTG.filt),
                SVZxSTG = rownames(SVZxSTG.filt),
                THAxSVZ = rownames(THAxSVZ.filt),
                MFGxSVZ = rownames(MFGxSVZ.filt),
                THAxSTG = rownames(THAxSTG.filt),
                THAxMFG = rownames(THAxMFG.filt))

all_de_together = unique(unlist(de_lists))
length(all_de_together)

#Create the matrix
deg_matrix <- matrix(nrow = 4, ncol = 4)
colnames(deg_matrix) <- c("MFG", "STG", "SVZ", "THA")
rownames(deg_matrix) <- c("MFG", "STG", "SVZ", "THA")

#Fill the matrix 
pos_5 = dim(MFGxSTG.filt)[1]
pos_9 = dim(MFGxSVZ.filt)[1]
pos_13 = dim(THAxMFG.filt)[1]
pos_10 = dim(SVZxSTG.filt)[1]
pos_14 = dim(THAxSTG.filt)[1]
pos_15 = dim(THAxSVZ.filt)[1]

deg_matrix[5] <- pos_5
deg_matrix[9] <- pos_9
deg_matrix[10] <- pos_10
deg_matrix[13] <- pos_13
deg_matrix[14] <- pos_14
deg_matrix[15] <- pos_15
deg_matrix

deg_matrix_m = melt(deg_matrix)

ggplot(deg_matrix_m, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(size = 1, color = "white") +
  geom_text(aes(label = ifelse(Var1==Var2,"NA",value))) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  xlab("") + 
  ylab("") +
  theme(panel.grid = element_blank())+
  easy_add_legend_title("DE genes")

