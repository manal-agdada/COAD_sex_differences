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
│-- Code/     # script for the analysis
│-- Figures/  # figures produced by the script
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

## **Data sources and workflow of the analysis**

[workflow.pdf](https://github.com/manal-agdada/TCGA_COAD_sex_differences/blob/main/Figures/workflow.pdf)


### **Preprocessing, normalization, and filtering**

### **Differential expression analysis**

### **Functional enrichment analysis**

### **Survival analysis**

## **Discussion and Conclusion**

## **References**

1. Wu Z, Huang Y, Zhang R, Zheng C, You F, Wang M, Xiao C, Li X. **Sex differences in colorectal cancer: with a focus on sex hormone-gut microbiome axis.** Cell Commun Signal. 2024 Mar 7;22(1):167. doi: 10.1186/s12964-024-01549-2. PMID: 38454453; PMCID: PMC10921775.

2. Nwaokorie A, Kolch W, Fey D. **A Systems Biology Approach to Understand the Racial Disparities in Colorectal Cancer.** Cancer Res Commun. 2024 Jan 12;4(1):103-117. doi: 10.1158/2767-9764.CRC-22-0464. PMID: 38051091; PMCID: PMC10785768.

3. Lopes-Ramos CM, Kuijjer ML, Ogino S, Fuchs CS, DeMeo DL, Glass K, Quackenbush J. **Gene Regulatory Network Analysis Identifies Sex-Linked Differences in Colon Cancer Drug Metabolism.** Cancer Res. 2018 Oct 1;78(19):5538-5547. doi: 10.1158/0008-5472.CAN-18-0454. Erratum in: Cancer Res. 2019 Apr 15;79(8):2084. doi: 10.1158/0008-5472.CAN-19-0678. PMID: 30275053; PMCID: PMC6169995.

4. Abancens M, Bustos V, Harvey H, McBryan J, Harvey BJ. **Sexual Dimorphism in Colon Cancer.** Front Oncol. 2020 Dec 9;10:607909. doi: 10.3389/fonc.2020.607909. PMID: 33363037; PMCID: PMC7759153.

5. Baraibar I, Ros J, Saoudi N, Salvà F, García A, Castells MR, Tabernero J, Élez E. **Sex and gender perspectives in colorectal cancer.** ESMO Open. 2023 Apr;8(2):101204. doi: 10.1016/j.esmoop.2023.101204. Epub 2023 Apr 3. PMID: 37018873; PMCID: PMC10163160.

