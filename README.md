# Exploring Sex Differences in Colon Adenocarcinoma: A Transcriptomic Approach

## **Background**

Colorectal cancer is the third most common cancer diagnosed and the second leading cause of cancer-related mortality worldwide [1]. The distribution of colorectal cancer differs across the worldwide population, with an evident difference in incidence and mortality by different demographic covariates, such as race, sex, and age [2].
Colon adenocarcinoma (COAD) shows clinically significant differences between sexes, with men tending to develop the disease earlier and often experiencing worse outcomes compared to women [3]. These disparities suggest that underlying biological mechanisms related to sex differences may contribute to these clinical observations [3-5].
Understanding the molecular basis for these differences is essential for improving therapeutic strategies and patient outcomes. This project aims to explore the transcriptomic differences between males and females with COAD using RNA-seq data from The Cancer Genome Atlas (TCGA), to investigate potential biomarkers for sex-specific outcomes and overall survival (OS).

## **Biological question**

The main question guiding this project is: *Are there fundamental biological differences between males and females in colon adenocarcinoma that influence disease progression and survival?*

Specifically, this study investigates whether sex-specific gene expression patterns can be linked to differences in prognosis and whether distinct gene signatures for males and females can predict OS outcomes.

## **Repository structure**

```
TCGA_COAD_sex_differences/
│-- Code/     # contains a script to reproduce the workflow and a script for the analysis
│-- Data/     # contains unprocessed and processed metadata
│-- Figures/  # contains figures produced by the script
│-- LICENSE
│-- README.md
```

## **Project dependencies**

This project is carried out using R 4.4.2 and the following R packages:
- **Bioconductor Packages**:  
  - `TCGAbiolinks`  
  - `SummarizedExperiment`  
  - `biomaRt`
- **CRAN Packages**:  
  - `gplots`
  - `ggplot2`
  - `dplyr`
  - `gtsummary`
  - `viridis`
  - `prodlim`
  - `survival`
  - `kableExtra`
  - `DiagrammeR`

## **Data sources and workflow of the analysis**

RNA-seq data of the TCGA-COAD cohort was downloaded directly into R using the Bioconductor package `TCGAbiolinks`. Metadata with survival information was downloaded from UCSC XENA and can be found in the [Data folder](https://github.com/manal-agdada/TCGA_COAD_sex_differences/tree/main/Data) alongside the processed metadata dataframe used for the analysis in this study.

Only primary tumor samples were directly downloaded to avoid biases from the presence of normal or metastatic samples.

Below, a schematic of the pipeline used in the analysis is presented. In particular, the analysis consists of:

1. Data selection, preprocessing, normalization, and filtering
2. Differential expression analysis using the pipeline edgeR
3. Functional enrichment analysis
4. Survival analysis with Cox PH regression models and Kaplan-Meier curves

![workflow](https://github.com/manal-agdada/TCGA_COAD_sex_differences/blob/main/Figures/workflow.png)


### **Selection, preprocessing, normalization, and filtering**

Only primary tumors samples were downloaded excluding normal and metastasis samples with a total of 474 samples. Also, samples from patients that were subjected to prior treatments before the collection and sequencing were excluded. 
The clinical covariantes considered in this study are the following: sex, age, race, prior malignancy, tumor site, tumor stage, OS time, and survival status. Therefore, samples with missing data were excluded, rounding up a total of 465 samples (244 male samples and 221 female samples).

Below, a demographics of the cohort under study is presented stratified by sex. All clinical data appear to be distributed homogeneously between male and female patients.

![table1](https://github.com/manal-agdada/TCGA_COAD_sex_differences/blob/main/Figures/table1_demographics.png)

The gene expression matrix under study comprises of 465 samples and 60660 genes. Normalization was applied using the `TCGAanalyze_Normalization` function from the `TCGAbiolinks` package accounting for sequencing depth and gene length. After normalization, 38691 genes were removed as lowly expressed genes, rounding up to a total of 21969 genes that will be used for subsequent analysis.

### **Differential expression analysis**

Differential expression analysis was run to identify genes that are differentially expressed between males and females in COAD. The analysis was performed using the `TCGAanalyze_DEA` function from the `TCGAbiolinks` package using `edgeR` as pipeline. A significant threshold of FDR < 0.01 and |logFC| > 1 was applied to identify 325 differentially expressed genes (DEGs), with 128 upregulated and 197 downregulated in males. Also, 276 of these DEGs are from autosomes, 19 from the X chromosome, and 26 from the Y chromosome.

Below, a barplot of the top 50 DEGs based on logFC and localization in the genome is presented.

![barplot1](https://github.com/manal-agdada/TCGA_COAD_sex_differences/blob/main/Figures/barplot_distribution_DEGs.png)

Below, the volcano plot showing the genes differentially expressed in males compared to females in terms of statistical significance (-log10FDR) and logFC is presented, highlighting the most significantly changing genes.

![volcanoplot](https://github.com/manal-agdada/TCGA_COAD_sex_differences/blob/main/Figures/volcanoplot_DEGs.png)

Also, the heatmap of the top 100 deregulated genes is shown below.

![heatmap](https://github.com/manal-agdada/TCGA_COAD_sex_differences/blob/main/Figures/heatmap_top100DEGs.png)

### **Functional enrichment analysis**

### **Survival analysis**

## **Discussion and Conclusion**

### **Limitations of the study**

## **References**

1. Wu Z, Huang Y, Zhang R, Zheng C, You F, Wang M, Xiao C, Li X. **Sex differences in colorectal cancer: with a focus on sex hormone-gut microbiome axis.** Cell Commun Signal. 2024 Mar 7;22(1):167. doi: 10.1186/s12964-024-01549-2. PMID: 38454453; PMCID: PMC10921775.

2. Nwaokorie A, Kolch W, Fey D. **A Systems Biology Approach to Understand the Racial Disparities in Colorectal Cancer.** Cancer Res Commun. 2024 Jan 12;4(1):103-117. doi: 10.1158/2767-9764.CRC-22-0464. PMID: 38051091; PMCID: PMC10785768.

3. Lopes-Ramos CM, Kuijjer ML, Ogino S, Fuchs CS, DeMeo DL, Glass K, Quackenbush J. **Gene Regulatory Network Analysis Identifies Sex-Linked Differences in Colon Cancer Drug Metabolism.** Cancer Res. 2018 Oct 1;78(19):5538-5547. doi: 10.1158/0008-5472.CAN-18-0454. Erratum in: Cancer Res. 2019 Apr 15;79(8):2084. doi: 10.1158/0008-5472.CAN-19-0678. PMID: 30275053; PMCID: PMC6169995.

4. Abancens M, Bustos V, Harvey H, McBryan J, Harvey BJ. **Sexual Dimorphism in Colon Cancer.** Front Oncol. 2020 Dec 9;10:607909. doi: 10.3389/fonc.2020.607909. PMID: 33363037; PMCID: PMC7759153.

5. Baraibar I, Ros J, Saoudi N, Salvà F, García A, Castells MR, Tabernero J, Élez E. **Sex and gender perspectives in colorectal cancer.** ESMO Open. 2023 Apr;8(2):101204. doi: 10.1016/j.esmoop.2023.101204. Epub 2023 Apr 3. PMID: 37018873; PMCID: PMC10163160.

