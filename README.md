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
Filtered by: adj.P.Val<0.05 & AveExpr > 0 & abs(logFC) > 0.5
```
| Comparison  | Genes |
| ------------- | ------ |
| MFG x STG  | 96  |
| THA x SVZ  | 1470 |
| MFG x SVZ  | 2378 |
| THA x STG  | 1109 |
| THA x MFG  | 940  |
| SVZ x STG  | 3950 |

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
 - [Depression x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_depressionxct_dream.html).

 - [Dementia x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_dementiaxct_dream.html).

 - [Psychiatric diagnosis x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_psychiatricDiagxct_dream.html).

 - [Schizophrenia + Bipolar x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_schizophrenia_bipolarxCT_dream.html).

 - [Parkinson's x CT](https://rajlabmssm.github.io/MiGA_public_release/DE_diagnosis/deg_pdxct_dream.html).

Number of DE genes by comparison from Dream:

| Comparison  | FDR 5% |
| ------------- | ------ |
| Depression x CT  | 0 |
| Dementia x CT  | 24 |
| Psychiatric diagnosis x CT  | 0 |
| Schizophrenia + Bipolar x CT  | 0 |
| Parkinson's x CT  | 0 |




