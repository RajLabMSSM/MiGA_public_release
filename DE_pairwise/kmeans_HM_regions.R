# Katia Lopes
# July 09,2020
# Kmeans approach for regional differences in the microglia dataset 

#This command clean all variables. BE CAREFULL!!! 
rm(list = setdiff(ls(), lsf.str()))

# install.packages("wesanderson")

library(wesanderson)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(venn)
library(ggsci)
library(NbClust)
library(factoextra)
library(RColorBrewer)

gene_lists = "~/pd-omics/katia/scripts/GitHub_scripts/glia_omics/3rd_pass_mic_255s/brain_regions/pairwise_dream/"
expression_dir = "~/ad-omics_hydra/microglia_omics/expression_tables/added_pilot_314s/expr_4brain_regions/"
# work_plots = "~/pd-omics/katia/scripts/GitHub_scripts/glia_omics/3rd_pass_mic_255s/brain_regions/kmeans_regions/"
work_plots = "/Users/katia/OneDrive/Documentos/MountSinai/Projects/Microglia/Figures4paper/Fig_region/"
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

## Expression and metadata 
load(paste0(expression_dir, "Expression_filt_255s.Rdata")) 
dim(genes_counts_voom_3rd)

load(paste0(metadata_path, "metadata_255filt_eng_29jun2020.Rdata"))

# To include cause of death and C1-C4
metadata3rd_pass$cause_of_death_categories[metadata3rd_pass$cause_of_death_categories %in% NA] <- "Other"
# table(metadata3rd_pass$cause_of_death_categories)

metadata3rd_pass$C1 = metadata3rd_pass$C1 %>% replace_na(median(metadata3rd_pass$C1, na.rm = T))
metadata3rd_pass$C2 = metadata3rd_pass$C2 %>% replace_na(median(metadata3rd_pass$C2, na.rm = T))
metadata3rd_pass$C3 = metadata3rd_pass$C3 %>% replace_na(median(metadata3rd_pass$C3, na.rm = T))
metadata3rd_pass$C4 = metadata3rd_pass$C4 %>% replace_na(median(metadata3rd_pass$C4, na.rm = T))

## Match names
genes_voom <- as.data.frame(genes_counts_voom_3rd[rownames(genes_counts_voom_3rd) %in% all_de_together, ])
dim(genes_voom)
genes_voom_median <- list()
genes_voom_median$MFG = apply(genes_voom[, metadata3rd_pass$tissue == "MFG"], 1, median, na.rm = T)
genes_voom_median$STG = apply(genes_voom[, metadata3rd_pass$tissue == "STG"], 1, median, na.rm = T)
genes_voom_median$SVZ = apply(genes_voom[, metadata3rd_pass$tissue == "SVZ"], 1, median, na.rm = T)
genes_voom_median$THA = apply(genes_voom[, metadata3rd_pass$tissue == "THA"], 1, median, na.rm = T)

genes_voom_median <- as.data.frame(genes_voom_median)

genes_voom_median_zscore = t(scale(t(genes_voom_median))) # transpose to normalize by row  
genes_voom_median_zscore = as.data.frame(genes_voom_median_zscore)

my_bar <- data.frame(Region = colnames(genes_voom_median_zscore))
rownames(my_bar) <- my_bar$Region

my_bar_color <- list(Region = c(MFG = "#FF6F00FF",
                                STG = "#C71000FF",
                                SVZ = "#008EA0FF",
                                THA = "#8A4198FF"))

set.seed(124)
my_plot <- pheatmap::pheatmap(genes_voom_median_zscore,
                   show_rownames = T,
                   border_color = "white",
                   treeheight_col = 0,
                   annotation_col = my_bar,
                   annotation_colors = my_bar_color,
                   cluster_cols = F,
                   cluster_rows = F,
                   legend = T,
                   kmeans_k = 4)

cluster_genes <- cbind(genes_voom_median_zscore, cluster = my_plot[["kmeans"]][["cluster"]])

my_bar_color2 <- list(Region = c(MFG = "#FF6F00FF",
                                STG = "#C71000FF",
                                SVZ = "#008EA0FF",
                                THA = "#8A4198FF"),
                      cluster = c(`1` = "#DC863B",
                                  `2` = "#F8AFA8",
                                  `3` = "#FDDDA0",
                                  `4` = "#74A089"))

annot_rows = cluster_genes[,"cluster",drop=F]
rownames(annot_rows) = rownames(cluster_genes)
annot_rows$cluster = as.factor(annot_rows$cluster)

pheatmap::pheatmap(genes_voom_median_zscore[order(annot_rows$cluster, rowMeans(genes_voom_median_zscore)),],
                   annotation_row = annot_rows,
                   gaps_row = cumsum(unclass(table(annot_rows$cluster))),
                   show_rownames = F,
                   border_color = "white",
                   treeheight_col = 0,
                   annotation_col = my_bar,
                   annotation_colors = my_bar_color2,
                   cluster_cols = F,
                   cluster_rows = F,
                   legend = T)

## Get conversion table for Gencode 30
gencode_30 = read.table("~/pd-omics/katia/ens.geneid.gencode.v30")
colnames(gencode_30) = c("ensembl","symbol")
cluster_genes$ensembl = rownames(cluster_genes)
cluster_genes_symbol = merge(cluster_genes, gencode_30, by = "ensembl")

### Final plot
genes_voom_median_zscore$ensembl = rownames(genes_voom_median_zscore)
genes_voom_median_zscore_symbol = merge(genes_voom_median_zscore, gencode_30, by="ensembl")
rownames(genes_voom_median_zscore_symbol) = genes_voom_median_zscore_symbol$symbol
genes_voom_median_zscore$ensembl = NULL
genes_voom_median_zscore_symbol$ensembl = NULL
genes_voom_median_zscore_symbol$symbol = NULL 

annot_rows$ensembl = rownames(annot_rows)
annot_rows_symbol = merge(annot_rows, gencode_30, by="ensembl")
rownames(annot_rows_symbol) = annot_rows_symbol$symbol
annot_rows$ensembl = NULL
annot_rows_symbol$ensembl = NULL
annot_rows_symbol$symbol = NULL 

labels4plot = read.table("/Users/katia/OneDrive/Documentos/MountSinai/Projects/Microglia/Figures4paper/Fig_region/labels4kmeans.txt",
                         header = T, stringsAsFactors = F)

########################
# Beautiful labels
########################
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(70)

genes_voom_median_zscore_symbol_ordered = genes_voom_median_zscore_symbol[order(annot_rows_symbol$cluster, rowMeans(genes_voom_median_zscore_symbol)),]

annot_rows_symbol_df = as.data.frame(annot_rows_symbol)
annot_rows_symbol_ordered = annot_rows_symbol_df[match(rownames(genes_voom_median_zscore_symbol_ordered),rownames(annot_rows_symbol_df)),,drop=F]
all(rownames(annot_rows_symbol_ordered) == rownames(genes_voom_median_zscore_symbol_ordered))

annot_rows_symbol_ordered$cluster = as.character(annot_rows_symbol_ordered$cluster)
annot_rows_symbol_ordered$cluster = paste0("Cluster ", annot_rows_symbol_ordered$cluster)
annot_rows_symbol_ordered$symbol = rownames(annot_rows_symbol_ordered)
annot_rows_symbol_ordered = annot_rows_symbol_ordered %>% group_by(cluster) %>% mutate(clu_n = paste0(cluster,"\n(n = ",n(),")")) %>% as.data.frame()
rownames(annot_rows_symbol_ordered) = annot_rows_symbol_ordered$symbol

# This step is important to put the labels in the same order! 
labels_match = inner_join(data.frame(genes = rownames(genes_voom_median_zscore_symbol_ordered)),unique(labels4plot), by = "genes")

ha = rowAnnotation(foo = anno_mark(labels_gp = gpar(fontsize = 5), 
                                   at = which( rownames(genes_voom_median_zscore_symbol_ordered) %in% labels_match$genes ), 
                                   labels = labels_match$genes))

col_ha = columnAnnotation(Region = c("MFG", "STG", "SVZ", "THA"),
                          col = list(Region = c(MFG="#FF6F00FF", STG="#C71000FF", SVZ="#008EA0FF", THA="#8A4198FF")),
                          na_col = "white",
                          border = F)

# Distances: (“euclidean”, “maximum”, “manhattan”, “canberra”, “binary”, “minkowski”, “pearson”, “spearman”, “kendall”)

set.seed(123) 
Heatmap(as.matrix(genes_voom_median_zscore_symbol_ordered[,c("MFG", "STG", "SVZ", "THA")]),
        cluster_rows = F,
        show_row_dend = F,
        right_annotation = ha,
        col = colors,
        show_column_dend = F,
        show_row_names = F,
        name = "Z score",
        top_annotation = col_ha,
        cluster_columns = F,
        row_split = annot_rows_symbol_ordered$clu_n,
        row_title_gp = gpar(fill = "white", border="white", fontsize = 10))


