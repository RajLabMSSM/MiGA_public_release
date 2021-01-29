# Katia Lopes
# July 21, 2020
# Samples grouped by donor for all genes

#This command clean all variables. BE CAREFULL!!! 
rm(list = setdiff(ls(), lsf.str()))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library("ggeasy")
library("ggsci")
library("ggplot2")
library("factoextra")
library("grid")
library("variancePartition")
library("BiocParallel")
library("edgeR")
library("statmod")
library("gplots")
library("ggpubr")
library("readxl")
library("pheatmap")
library("reshape2")
library("kableExtra")
library("RColorBrewer")
require("gridExtra")
library("ggfortify")
library("doParallel")
library("dplyr")
library("tidyr")
library("venn")
library("viridis")

expression_dir = "~/ad-omics_hydra/microglia_omics/expression_tables/added_pilot_314s/expr_4brain_regions/"
work_plots = "~/pd-omics/katia/scripts/GitHub_scripts/glia_omics/3rd_pass_mic_255s/age_related/age_Dream_3rd/"
metadata_path = "~/pd-omics/katia/Microglia/mic_255s/metadata_files/"
gene_lists = "~/pd-omics/katia/Microglia/GeneLists_4analysis/"
plots4paper = "/Users/katia/OneDrive/Documentos/MountSinai/Projects/Microglia/Figures4paper/Fig_age/age_heatmaps/by_donor/binned_equal_groups/"

# Function to cluster the columns 
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

################# Load expression 
load(paste0(expression_dir, "Expression_filt_255s.Rdata")) 
dim(genes_counts_exp_3rd) # 19376   255

load(paste0(metadata_path, "metadata_255filt_eng_29jun2020.Rdata"))
str(metadata3rd_pass)
#################

metadata <- metadata3rd_pass
ageByDonor = unique(metadata[,c("donor_id", "age", "sex")])

ggplot(ageByDonor, aes(x=age, fill=age)) +
  geom_histogram(bins = 25, colour='black', position = "stack", fill="#5c88da99") +
  labs(x="Age", y="Donors") +
  scale_y_continuous(breaks = (1:20)) +
  scale_x_continuous(breaks=seq(20,120,10)) + 
  theme_classic()

### For the metadata 
metadata$cause_of_death_categories[metadata$cause_of_death_categories %in% NA] <- "Other"
table(metadata$cause_of_death_categories)

metadata$C1 = metadata$C1 %>% replace_na(median(metadata$C1, na.rm = T))
metadata$C2 = metadata$C2 %>% replace_na(median(metadata$C2, na.rm = T))
metadata$C3 = metadata$C3 %>% replace_na(median(metadata$C3, na.rm = T))
metadata$C4 = metadata$C4 %>% replace_na(median(metadata$C4, na.rm = T))

# Matching the names from the DE analysis to use the same code 
genes_counts4deg = genes_counts_exp_3rd
metadata4deg = metadata
all(colnames(genes_counts_exp_3rd) == metadata$donor_tissue) # Check the order of columns - TRUE 

## Load age dream analysis Rdata 
load("~/pd-omics/katia/Microglia/mic_314s/age_dream_rfiles/age_dream_objects_r.Rdata")

load("~/pd-omics/katia/scripts/GitHub_scripts/glia_omics/3rd_pass_mic_255s/age_related/age_Dream_3rd/age_dream_Residuals_3rd.Rdata")

################# Results from Github page 
age_genes_all <- read.table(paste0(work_plots, "List_age_dream_255s.txt"), header = T, stringsAsFactors = F)
age_genes_sign <- age_genes_all[age_genes_all$adj.P.Val < 0.05, ]

## Get conversion table for Gencode 30
gencode_30 = read.table("~/pd-omics/katia/ens.geneid.gencode.v30")
colnames(gencode_30) = c("ensembl","symbol")

residuals_matrix = as.data.frame(allResiduals)
residuals_matrix$ensembl = rownames(residuals_matrix)
residuals_matrix = merge(residuals_matrix, gencode_30, by="ensembl")

#Remove duplicates 
residuals_matrix <- residuals_matrix[! duplicated(residuals_matrix$symbol),]
rownames(residuals_matrix) = residuals_matrix$symbol
residuals_matrix$ensembl = NULL
residuals_matrix$symbol = NULL
dim(residuals_matrix) #19342   255

ourexpr_age_res = residuals_matrix[rownames(residuals_matrix) %in% age_genes_sign$symbol , ] # RESIDUALS TABLE 
col_age = as.data.frame(metadata4deg[c("age")] ,)
rownames(col_age) = metadata4deg$donor_tissue
summary(col_age)
ourexpr_age_zscore = t(scale(t(ourexpr_age_res))) # The scale is by column. I want the scale by row so I need to transpose the matrix. Then, I transposed again for the plot! 
ourexpr_age_zscore = as.data.frame(ourexpr_age_zscore)
ourexpr_age_zscore_sorted = ourexpr_age_zscore[, order(col_age)]

##################################################
# Binned by groups of age
# Median by donor and age binned
##################################################

dim(ourexpr_age_zscore_sorted) # Residuals 

temp_melt = melt(ourexpr_age_zscore_sorted, id.vars = 0)
temp_melt$symbol = rownames(ourexpr_age_zscore_sorted)

temp_melt_group = temp_melt %>% left_join(metadata4deg[, c("donor_tissue", "donor_id", "age")], by = c("variable" = "donor_tissue")) %>%
  #mutate(Age_bin = cut(age,breaks=c(20,30,40,50,60,70,80,90,100,110))) %>%
  mutate(Age_bin = cut(age,breaks=c(20,60,70,80,90,110))) %>%
  group_by(symbol, donor_id) %>% mutate(median_exp_per_donor = median(value)) %>%
  select(symbol, donor_id, Age_bin, median_exp_per_donor) %>%
  ungroup() %>% group_by(symbol, Age_bin) %>% summarise(median_exp = median(median_exp_per_donor))

temp_matrix = as.data.frame(temp_melt_group %>% pivot_wider(names_from = Age_bin, values_from = median_exp))
dim(temp_matrix)
rownames(temp_matrix) = temp_matrix$symbol
head(temp_matrix)
residual_matrix_age = temp_matrix %>% select(-symbol)
# ourexpr_age = residual_matrix_age[rownames(residual_matrix_age) %in% top_age_genes$symbol , ] # RESIDUALS TABLE 
ourexpr_age = residual_matrix_age[rownames(residual_matrix_age) %in% rownames(ourexpr_age_zscore_sorted) , ] # RESIDUALS TABLE 
ourexpr_age_zscore = t(scale(t(ourexpr_age))) 

col_age = data.frame(Age = colnames(residual_matrix_age))
col_age$Age = factor(col_age$Age, levels = col_age$Age)
rownames(col_age) = colnames(residual_matrix_age)
col_colors <- list(Age = magma(length(levels(col_age$Age)), direction = 1))
names(col_colors$Age) <- levels(col_age$Age)
my_pallete = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(80)
my.breaks <- seq(-2, 2, length.out = 80)
newnames <- lapply(
  rownames(ourexpr_age_zscore),
  function(x) bquote(italic(.(x))))

### Selecting the Patir gene lists
signature_genes2 = read.table(paste0(gene_lists, "Patir_core_249g.txt"), stringsAsFactors = F, check.names = F, sep = "\t", header = T)

##################################################
# Binned by groups of 20 donors 
##################################################
dim(ourexpr_age_zscore_sorted) # Residuals 

temp_melt = melt(ourexpr_age_zscore_sorted, id.vars = 0)
temp_melt$symbol = rownames(ourexpr_age_zscore_sorted)

temp_melt_group_1 = temp_melt %>% left_join(metadata4deg[, c("donor_tissue", "donor_id", "age")], by = c("variable" = "donor_tissue")) %>%
  group_by(symbol, donor_id) %>% summarise(median_exp_per_donor = median(value), age) %>% distinct()
age_g = unique(temp_melt_group_1[,c("donor_id","age")])
age_g$Age_bin = ntile(age_g$age, 5) #split in 5 quantiles
temp_melt_group_1 = temp_melt_group_1 %>% left_join(age_g[,c("donor_id","Age_bin")], by = "donor_id")
temp_melt_group = temp_melt_group_1 %>% select(symbol, donor_id, Age_bin, median_exp_per_donor) %>%
  ungroup() %>% group_by(symbol, Age_bin) %>% summarise(median_exp = median(median_exp_per_donor))

temp_matrix = as.data.frame(temp_melt_group %>% pivot_wider(names_from = Age_bin, values_from = median_exp))
dim(temp_matrix)
rownames(temp_matrix) = temp_matrix$symbol
head(temp_matrix)
residual_matrix_age = temp_matrix %>% select(-symbol)
# ourexpr_age = residual_matrix_age[rownames(residual_matrix_age) %in% top_age_genes$symbol , ] # RESIDUALS TABLE 
ourexpr_age = residual_matrix_age[rownames(residual_matrix_age) %in% rownames(ourexpr_age_zscore_sorted) , ] # RESIDUALS TABLE 
ourexpr_age_zscore = t(scale(t(ourexpr_age))) 

col_age = data.frame(Age = colnames(residual_matrix_age))
col_age$Age = factor(col_age$Age, levels = col_age$Age)
rownames(col_age) = colnames(residual_matrix_age)
col_colors <- list(Age = magma(length(levels(col_age$Age)), direction = 1))
names(col_colors$Age) <- levels(col_age$Age)
my_pallete = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(80)
my.breaks <- seq(-2, 2, length.out = 80)
newnames <- lapply(
  rownames(ourexpr_age_zscore),
  function(x) bquote(italic(.(x))))

pheatmap::pheatmap(ourexpr_age_zscore, 
                   annotation_col = col_age, 
                   annotation_colors = col_colors,
                   show_colnames = F,
                   show_rownames = F,
                   color = my_pallete,
                   breaks = my.breaks,
                   border_color = "white",
                   cluster_cols = F,
                   clustering_callback = callback,
                   clustering_method = "ward.D2",
                   treeheight_row = 0)

## Checking age 
table(age_g$Age_bin)
max(age_g$age[age_g$Age_bin == 1]) 
min(age_g$age[age_g$Age_bin == 1])

max(age_g$age[age_g$Age_bin == 2]) 
min(age_g$age[age_g$Age_bin == 2])

max(age_g$age[age_g$Age_bin == 3]) 
min(age_g$age[age_g$Age_bin == 3])

max(age_g$age[age_g$Age_bin == 4]) 
min(age_g$age[age_g$Age_bin == 4])

max(age_g$age[age_g$Age_bin == 5]) 
min(age_g$age[age_g$Age_bin == 5])

#######################################################################
# Final plot with beautiful labels
#######################################################################

labels2remove = c("PTAFR", "DHRS9", "TLR10", "SOWAHD", "SELPLG", "RIN3", "HCLS1")
signature_genes2_sel = signature_genes2[! signature_genes2$symbol %in% labels2remove ,] 

col_ha = columnAnnotation(Age = names(col_colors$Age),
                          col = list(Age = c(col_colors$Age)),
                          na_col = "white",
                          border = F)

labels_match = signature_genes2_sel[signature_genes2_sel$symbol %in% rownames(ourexpr_age_zscore) ,]
labels_match = labels_match$symbol

ha = rowAnnotation(foo = anno_mark(labels_gp = gpar(fontsize = 8), at = which(rownames(ourexpr_age_zscore) %in% signature_genes2_sel$symbol), 
                                   labels = labels_match))

colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(70)

set.seed(123) 
Heatmap(ourexpr_age_zscore,
        cluster_rows = T,
        right_annotation = ha,
        show_row_dend = F,
        show_column_dend = F,
        show_row_names = F,
        col = colors,
        show_column_names = F,
        name = "Z score",
        clustering_method_rows = "ward.D2",
        clustering_distance_rows = "pearson",
        top_annotation = col_ha) 


