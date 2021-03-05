## Omics analysis for the Microglia project

> This repository includes code and plots. Exploratory analysis and intermediate processing files are too large for this repository.

> Full nominal and permuted eQTL and sQTL summary statistics per brain region are available from Zenodo at https://doi.org/10.5281/zenodo.4118605 (eQTL) and https://doi.org/10.5281/zenodo.4118403 (sQTL). Results for eQTL and sQTL meta-analysis (mashR and METASOFT) and colocalization (COLOC) are available from Zenodo at https://doi.org/10.5281/zenodo.4118676. Allelic information for all QTL analyses is available from Zenodo at https://doi.org/10.5281/zenodo.4301005.

<p align="center">
 <img src="https://github.com/RajLabMSSM/MiGA_public_release/blob/main/Fig1.png?raw=true">
</p>

**Figure 1:** Overview of the MiGA study. 

1.Looking at covariates:
 - [Variance partition](https://rajlabmssm.github.io/MiGA_public_release/exploratory_analysis/01_VP_255s.html) plots. Includes pairwise correlation plot, VP for 255 samples, VP by region and VP for European only (N=216).

2.Exploratory Analysis:
 - [Exploratory plots after filters](https://rajlabmssm.github.io/MiGA_public_release/exploratory_analysis/02_exploratory_filtered.html). General counts for age, sex, cause of death ... with 255 samples.

 - The [PCAs are here](https://rajlabmssm.github.io/MiGA_public_release/exploratory_analysis/03_PCAs_3rd.html).

 - [Linear regression](https://rajlabmssm.github.io/MiGA_public_release/exploratory_analysis/linear_reg_pinkheatmap.html) between the first 15 PCs and the covariates.

3.[Expression of markers](https://rajlabmssm.github.io/MiGA_public_release/exploratory_analysis/04_check_markers.html) genes for 255 samples.

***************************************
Differential expression analysis - PAIRWISE
***************************************
```
Covariates: ~ sex + (1|donor_id) + age + tissue + (1|cause_of_death_categories) + C1 + C2 + C3 + C4 + picard_pct_mrna_bases + picard_summed_median + picard_pct_ribosomal_bases
```
 - [MFG x STG](https://rajlabmssm.github.io/MiGA_public_release/DE_pairwise/MFGxSTG_3rd.html)

 - [THA x SVZ](https://rajlabmssm.github.io/MiGA_public_release/DE_pairwise/THAxSVZ_3rd.html)

 - [MFG x SVZ](https://rajlabmssm.github.io/MiGA_public_release/DE_pairwise/MFGxSVZ_3rd.html)

 - [THA x STG](https://rajlabmssm.github.io/MiGA_public_release/DE_pairwise/THAxSTG_3rd.html)

 - [THA x MFG](https://rajlabmssm.github.io/MiGA_public_release/DE_pairwise/THAxMFG_3rd.html)

 - [SVZ x STG](https://rajlabmssm.github.io/MiGA_public_release/DE_pairwise/SVZxSTG_3rd.html)

Number of DE genes by comparison from Dream:
```
Filtered by: adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 1
```
| Comparison  | Genes |
| ------------- | ------ |
| MFG x STG  | 6 |
| THA x SVZ  | 387 |
| MFG x SVZ  | 609 |
| THA x STG  | 105 |
| THA x MFG  | 119 |
| SVZ x STG  | 909 |

***************************************
Sex-related analysis
***************************************

[Our sex-related analysis](https://rajlabmssm.github.io/MiGA_public_release/sex-related_analysis/sex_Dream_3rd.html) = 0 genes at FDR 5%

***************************************
Differential expression analysis for different DIAGNOSIS
***************************************
```
Covariates: ~ main_diagnosis + sex + (1|donor_id) + age + tissue + (1|cause_of_death_categories) + C1 + C2 + C3 + C4 + picard_pct_mrna_bases + picard_summed_median + picard_pct_ribosomal_bases
```
 - [Depression x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_depressionxct_dream.html)

 - [Dementia x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_dementiaxct_dream.html)

 - [Psychiatric diagnosis x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_psychiatricDiagxct_dream.html)

 - [Schizophrenia + Bipolar x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_schizophrenia_bipolarxCT_dream.html)

 - [Parkinson's x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_pdxct_dream.html)

Number of DE genes by comparison from Dream:

| Comparison  | FDR 5% |
| ------------- | ------ |
| Depression x CT  | 0 |
| Dementia x CT  | 24 |
| Psychiatric diagnosis x CT  | 0 |
| Schizophrenia + Bipolar x CT  | 0 |
| Parkinson's x CT  | 0 |

***************************************
Age-related Analysis
***************************************

[Our age-related](https://rajlabmssm.github.io/MiGA_public_release/age-related_analysis/age_dream_3ndpass.html) analysis using the tool Dream (N=255). removeBatchEffect to remove confounders from the expression table. 

[Effect size correlation](https://rajlabmssm.github.io/MiGA_public_release/age-related_analysis/effect_size_age.html) plots.

DE genes with [interaction term](https://rajlabmssm.github.io/MiGA_public_release/age-related_analysis/age_dream_3ndpass_interaction.html). Include scatter plots by region.

[Age-related for splicing](https://rajlabmssm.github.io/MiGA_public_release/age-related_analysis/age_splicing_analysis/age_related_sqtl.html) data. 

***************************************
Genotypes
***************************************

Number of samples for downstream analysis. After GSA includes European only:

| Tissue  | Before GSA filter | After GSA filter |
| ------------- | ------ | ------ |
| MFG  | 78  | 63  |
| STG  | 64  | 55  |
| THA  | 60  | 53  |
| SVZ  | 56  | 45  |
| CC  | 17  | 16  |
| CER  | 7  | 4  |
| SN  | 3  | 0  |
| Total | 285 samples | 236 samples |

Checking the numbers of samples from GSA. [Click here](https://rajlabmssm.github.io/MiGA_public_release/genotype/GSA/Microglia_Genotypes.html).

Number of eGenes by brain region, including different number of PEER factors. [Click here](https://rajlabmssm.github.io/MiGA_public_release/qtl_pages/sig_eGenes_peer_2nd.html).

You can check the lists of **eGenes**, by brain region [in this link](https://rajlabmssm.github.io/MiGA_public_release/qtl_pages/eQTLs_list_2nd.html).

Number of eGenes (includes covariate file with age, sex and 4 MDS ancestry):

| Tissue  | Sample size | eGenes at qval 5% | eGenes at qval 10% |
| ------------- | ------ | ------ | ------ |
| MFG  | 63  | 199  | 267  |
| STG  | 55  | 127  | 154  |
| THA  | 53  | 148  | 196  |
| SVZ  | 45  | 67  | 94  |

***************************************
mashR pages - eQTL
***************************************

The page with mashR results is [here](https://rajlabmssm.github.io/MiGA_public_release/mashr_pages/mashr_eqtl/pos_mashR_2nd.html). Include a table with significative genes (lfsr < 0.10), heatmap of pairwise sharing by magnitude, metaplots and posterior.beta table.

Tissue specific effects: [lfsr at 10%](https://rajlabmssm.github.io/MiGA_public_release/mashr_pages/mashr_eqtl/pos_mashR_specific_10per_2nd.html). Include lists of gene-SNP by region.

Tissue specific effects: [lfsr at 5%](https://rajlabmssm.github.io/MiGA_public_release/mashr_pages/mashr_eqtl/pos_mashR_specific_5per_2nd.html). Include lists of gene-SNP by region.

***************************************
mashR pages - sQTL
***************************************

The page with mashR results is [here](https://rajlabmssm.github.io/glia_omics/3rd_pass_mic_255s/mash_pages/pos_mashR_sQTL_3rd.html). Includes a table with significative **junction clusters** (lfsr < 0.05) and a heatmap of pairwise sharing by magnitude.

Tissue specific effects: [lfsr at 5%](https://rajlabmssm.github.io/glia_omics/3rd_pass_mic_255s/mash_pages/pos_mashR_sQTL_specific_5per3rd.html). Include lists of cluster-SNP by region.

**************************************
Colocalization and fine-mapping
***************************************

 - [Qvalue sharing](https://rajlabmssm.github.io/MiGA_public_release/qtl_analysis/qvalue_sharing.html). Includes heatmap. 

 - [Locus zoom plots](https://rajlabmssm.github.io/MiGA_public_release/qtl_analysis/locus_zoom_plots.html).

 - [Full COLOC visualization](https://rajlabmssm.github.io/MiGA_public_release/qtl_analysis/coloc_results.html).

 - Fine Mapping [COLOC PLAC-Seq](https://rajlabmssm.github.io/MiGA_public_release/qtl_analysis/coloc_fine_mapping_plac_seq.html).

 - Fine-mapping and [COLOC overlaps](https://rajlabmssm.github.io/MiGA_public_release/qtl_analysis/coloc_fine_mapping_overlaps.html).


