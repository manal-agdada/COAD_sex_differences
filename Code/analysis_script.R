######################################################################################
#### Exploring Sex Differences in Colon Adenocarcinoma: A Transcriptomic Approach #### 
######################################################################################


# set working directory 
setwd("C:/Users/agdad/Desktop/New folder/HackBio/project_COAD")


# load libraries ####
library("TCGAbiolinks")
library(SummarizedExperiment)
library(dplyr)
library(gtsummary)
library(biomaRt)
library(ggplot2)
library(gplots)
library(viridis)
library(prodlim)
library(survival)


# download the data ####
tcga_coad <- GDCquery(project = 'TCGA-COAD',
                      data.category = 'Transcriptome Profiling',
                      experimental.strategy = 'RNA-Seq',
                      data.type = 'Gene Expression Quantification',
                      sample.type = "Primary Tumor")
GDCdownload(tcga_coad) 
coad_data <- GDCprepare(tcga_coad) 

# get the raw counts
coad_rawdata <- assays(coad_data) 
dim(coad_rawdata$unstranded) # 481 samples, 60660 genes
coad_rawdata <- coad_rawdata$unstranded 
table(is.na(coad_rawdata)) # no NAs

# load the metadata file downloaded from UCSC XENA
metadata <- read.delim(gzfile("TCGA-COAD.clinical.tsv.gz"), sep="\t", header=TRUE)
View(metadata) # 562 samples

# filter only primary tumors
metadata <- metadata %>% 
  filter(metadata$sample_type.samples == "Primary Tumor")
dim(metadata) # 474 samples

metadata <- data.frame("barcode" = metadata$sample,
                       "Race" = metadata$race.demographic,
                       "Gender" = metadata$gender.demographic,
                       "days_to_last_follow_up" = metadata$days_to_last_follow_up.diagnoses,
                       "Age (Years)" = metadata$age_at_index.demographic,
                       "days_to_death" = metadata$days_to_death.demographic,
                       "Overall Survival Status" = metadata$vital_status.demographic,
                       "Prior Malignancy" = metadata$prior_malignancy.diagnoses,
                       "prior_treatment" = metadata$prior_treatment.diagnoses,
                       "tumor_stage" = metadata$ajcc_pathologic_stage.diagnoses,
                       "Primary Site" = metadata$site_of_resection_or_biopsy.diagnoses)

# categorization of tumor_stage: I, IA, II, IIA, IIB, IIC, III, IIIA, IIIB, IIIC, IV, IVA, IVB in I, II, III, IV
metadata$tumor_stage <- ifelse(metadata$tumor_stage %in% c("Stage I", "Stage IA"), "Stage I",
                               ifelse(metadata$tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"), "Stage II",
                                      ifelse(metadata$tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"), "Stage III",
                                             ifelse(metadata$tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB"), "Stage IV",
                                                    "Not Reported"))))
metadata <- metadata %>%
  rename(`Tumor Stage` = tumor_stage)

# create overall survival variable
metadata <- metadata %>% 
  mutate(
    Overall.Survival..Months = ifelse(
      Overall.Survival.Status == "Dead",
      days_to_death / 30.44,
      days_to_last_follow_up / 30.44
    )
  )

# convert the missing values with NA
metadata[metadata == ""] <- NA

# data quality check ####
table(is.na(metadata$Gender)) # 2 NAs
table(is.na(metadata$Overall.Survival..Months)) # 3 NAs
table(is.na(metadata$Overall.Survival.Status)) # 2 NAs
table(is.na(metadata$prior_treatment)) # 2 NAs + 3 yes prior treatment
table(is.na(metadata$Race)) # 2 NAs
table(is.na(metadata$Age..Years.)) # 2 NAs
table(is.na(metadata$Prior.Malignancy)) # 2 NAs
table(is.na(metadata$Primary.Site)) # 2 NAs

# filter the samples not annotated for sex and the samples that received prior treatment
metadata <- metadata %>% 
  filter(!is.na(metadata$Gender) &
           metadata$prior_treatment != "Yes" &
           !is.na(metadata$prior_treatment) &
           !is.na(metadata$Overall.Survival..Months) &
           !is.na(metadata$Overall.Survival.Status) &
           !is.na(metadata$Race) &
           !is.na(metadata$Age..Years.) &
           !is.na(metadata$Prior.Malignancy) &
           !is.na(metadata$Primary.Site))
dim(metadata) # 468 samples
table(metadata$Gender) # 221 female, 247 male

# change names of the covariates to prepare for the patients table
metadata <- metadata %>%
  rename(`Age (Years)` = Age..Years.)
metadata <- metadata %>%
  rename(`Overall Survival Status` = Overall.Survival.Status)
metadata <- metadata %>%
  rename(`Prior Malignancy` = Prior.Malignancy)
metadata <- metadata %>%
  rename(`Primary Site` = Primary.Site)
metadata <- metadata %>%
  rename(`Overall Survival (Months)` = Overall.Survival..Months)
metadata <- metadata %>%
  mutate(Race = case_when(
    Race == "american indian or alaska native" ~ "American Indian or Alaska Native",
    Race == "asian" ~ "Asian",
    Race == "black or african american" ~ "Black or African American",
    Race == "white" ~ "White",
    Race == "not reported" ~ "Not Reported"  
  ))
metadata <- metadata %>%
  mutate(Gender = case_when(
    Gender == "female" ~ "Female",
    Gender == "male" ~ "Male",
  ))

# match samples in coad_rawdata and metadata 
metadata$full_barcode <- colnames(coad_rawdata)[match(metadata$barcode, substr(colnames(coad_rawdata), 1, 16))]
dim(metadata) # 468
table(is.na(metadata$full_barcode)) # 3 NAs

metadata <- metadata %>% 
  filter(!is.na(metadata$full_barcode))
dim(metadata) # 465
table(metadata$Gender) # 221 female, 244 male

coad_rawdata <- coad_rawdata[, colnames(coad_rawdata) %in% metadata$full_barcode]
dim(coad_rawdata) # 465 samples, 60660 genes

setdiff(colnames(coad_rawdata), metadata$full_barcode)  
setdiff(metadata$full_barcode, colnames(coad_rawdata))  # all samples in the matrix and metadata coincide

metadata$barcode <- metadata$full_barcode

# patient table clinical characteristics
metadata %>%
  dplyr::select(Gender, `Age (Years)`, `Overall Survival Status`, `Overall Survival (Months)`, `Prior Malignancy`, `Tumor Stage`, Race, `Primary Site`) %>%
  tbl_summary(by = Gender) %>%
  add_overall() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Gender**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  modify_caption("**Table 1. Patient Characteristics**") %>%
  bold_labels()

# save metadata file
write.csv(metadata, "TCGA_COAD_processed_metadata.csv")

# normalize and filter the raw counts
norm_data <- TCGAanalyze_Normalization(tabDF = coad_rawdata, geneInfo = geneInfoHT, method = "geneLength")

filt_data <- TCGAanalyze_Filtering(tabDF = norm_data,
                                   method = "quantile",
                                   qnt.cut = 0.25)
coad <- filt_data
dim(coad) # filtered out 38691 genes 


# differential expression analysis ####
female <- metadata$barcode[metadata$Gender == "Female"]
male <- metadata$barcode[metadata$Gender == "Male"]

DEA <- TCGAanalyze_DEA(mat1 = coad[ , female], # control
                       mat2 = coad[ , male], # case
                       Cond1type = "Female",
                       Cond2type = "Male",
                       pipeline = "edgeR")
DEA.Level <- 
  TCGAanalyze_LevelTab(DEA, "Female", "Male",
                       coad[ , female],
                       coad[ , male])

# annotation
gene_names <- rownames(DEA.Level)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description', 'chromosome_name'),
  filters = 'ensembl_gene_id',
  values = gene_names,
  mart = mart)

# merge 
DEA.Level$ensembl_gene_id <- rownames(DEA.Level)
DEA.Level <-  DEA.Level %>%
  left_join(annot, by = 'ensembl_gene_id')

# filtering significant genes
DEGs <- DEA.Level
DEGs$status <- "NO"
DEGs$status[DEGs$logFC > 1 & DEGs$FDR < 0.01] <- "UP"
DEGs$status[DEGs$logFC < (-1) & DEGs$FDR < 0.01] <- "DOWN"
table(DEGs$status) # down 128, up 197
DEGs_df <- DEGs[DEGs$status == "UP"| DEGs$status == "DOWN", ] # 325 DEGs 

# classify the DEGs based on the position
DEGs_df$chromosome_name <- ifelse(DEGs_df$chromosome_name == "X", "Sex Chromosome X",
                           ifelse(DEGs_df$chromosome_name == "Y", "Sex Chromosome Y", "Autosome"))
table(DEGs_df$chromosome_name) # 276 autosome, 45 sex chromosome (26 Y, 19 X)

# visualization autosome-sex chromosome
DEGs_df <- DEGs_df %>%
  distinct(hgnc_symbol, .keep_all = TRUE)
DEGs_df$Gene <- factor(DEGs_df$hgnc_symbol, levels = DEGs_df$hgnc_symbol[order(DEGs_df$logFC)])
top50_DEGs <- DEGs_df %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 50)
ggplot(top50_DEGs, aes(x = logFC, y = Gene)) +
  geom_bar(stat = "identity", aes(fill = chromosome_name)) + 
  scale_fill_manual(values = c("Sex Chromosome X" = "red", "Sex Chromosome Y" = "blue", "Autosome" = "yellow")) +  
  theme_minimal() +
  labs(title = "Top 50 DEGs\nMale vs Female",
       x = "Log2 Fold Change (logFC)",
       y = "Genes",
       fill = "Chromosome Type") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.53, face = "bold"))

# volcano plot
Up10 <- DEA.Level %>%
  slice_max(order_by = logFC, n = 10)
Down10 <- DEA.Level %>%
  slice_min(order_by = logFC, n = 10)
Top20 <- bind_rows(Up10, Down10)

ggplot(DEA.Level, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(FDR < 0.01 & abs(logFC) > 1,
                                ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not significant")), size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR", title = "Volcano Plot of DEGs\nMale vs Female", color = "Gene Regulation") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top") +
  ggrepel::geom_label_repel(data = Top20, aes(label = hgnc_symbol), max.overlaps = Inf)

# heatmap
# select the top 100 deregulated genes
Up50 <- DEA %>%
  slice_max(order_by = logFC, n = 50)
Down50 <- DEA %>%
  slice_min(order_by = logFC, n = 50)
Top100 <- bind_rows(Up50, Down50)
heat.DEGs <- coad[rownames(Top100), ]
heat.DEGs <- log2(heat.DEGs + 1)
ccodes <- ifelse(colnames(heat.DEGs) %in% female, "red", "blue") # female red, male blue
heatmap.2(x = as.matrix(heat.DEGs),
          col = viridis(10),
          Rowv = F, Colv = T,
          sepcolor = 'black',
          trace = "none",
          scale = 'row',
          key = TRUE,
          key.title = "Expression",
          key.xlab = "Z-score",
          dendrogram = "column",
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap of top 100 DEGs\nMale vs Female\n(N=465)",
          na.color = 'black',
          ColSideColors = ccodes,
          labRow = NA,
          labCol = NA,
          margins = c(11,10))
legend("left",                       
       legend = c("Female", "Male"),     
       col = c("red", "blue"),           
       lty = 1,                          
       lwd = 4,                          
       cex = 0.55,
       xpd = TRUE,
       title = "Gender",
       inset = c(-0.5, 0)) 

# visualization of distribution of genes (only on autosome and X except Xist) 
DEGs_df_filtered <- DEGs_df %>% 
  filter(!(chromosome_name == "Sex Chromosome Y" | Gene == "XIST"))
View(DEGs_df_filtered)
DEGs_df_filtered$Gene <- factor(DEGs_df_filtered$hgnc_symbol, levels = DEGs_df_filtered$hgnc_symbol[order(DEGs_df_filtered$logFC)])
top50_DEGs_filtered <- DEGs_df_filtered %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 50)
ggplot(top50_DEGs_filtered, aes(x = logFC, y = Gene)) +
  geom_bar(stat = "identity", aes(fill = chromosome_name)) + 
  scale_fill_manual(values = c("Sex Chromosome X" = "red", "Autosome" = "yellow")) +  
  theme_minimal() +
  labs(title = "Top 50 DEGs\nMale vs Female",
       x = "Log2 Fold Change (logFC)",
       y = "Genes",
       fill = "Chromosome Type") +
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.53, face = "bold"))

# volcano plot
DEA.Level_filtered <- DEA.Level %>% 
  filter(!(chromosome_name == "Y" | hgnc_symbol == "XIST"))
dim(DEA.Level_filtered)
Up10_filtered <- DEA.Level_filtered %>%
  slice_max(order_by = logFC, n = 10)
Down10_filtered <- DEA.Level_filtered %>%
  slice_min(order_by = logFC, n = 10)
Top20_filtered <- bind_rows(Up10_filtered, Down10_filtered)

ggplot(DEA.Level_filtered, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(FDR < 0.01 & abs(logFC) > 1,
                                ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not significant")), size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR", title = "Volcano Plot of DEGs\nMale vs Female", color = "Gene Regulation") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top") +
  ggrepel::geom_label_repel(data = Top20_filtered, aes(label = hgnc_symbol), max.overlaps = Inf)

# heatmap
# select the top 100 deregulated genes filtered
DEA_filtered <- DEA
gene_names <- rownames(DEA_filtered)
annot <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name'),
  filters = 'ensembl_gene_id',
  values = gene_names,
  mart = mart)
DEA_filtered$ensembl_gene_id <- rownames(DEA_filtered)
DEA_filtered <-  DEA_filtered %>%
  left_join(annot, by = 'ensembl_gene_id')
DEA_filtered <- DEA_filtered %>%
  distinct(hgnc_symbol, .keep_all = TRUE)
DEA_filtered <- DEA_filtered %>% 
  filter(!(chromosome_name == "Y" | hgnc_symbol == "XIST"))
Up50_filtered <- DEA_filtered %>%
  slice_max(order_by = logFC, n = 50)
Down50_filtered <- DEA_filtered %>%
  slice_min(order_by = logFC, n = 50)
Top100_filtered <- bind_rows(Up50_filtered, Down50_filtered)
heat.DEGs_filtered <- coad[rownames(coad) %in% Top100_filtered$ensembl_gene_id, ]
heat.DEGs_filtered <- log2(heat.DEGs_filtered + 1)
ccodes <- ifelse(colnames(heat.DEGs_filtered) %in% female, "red", "blue") # female red, male blue
heatmap.2(x = as.matrix(heat.DEGs_filtered),
          col = viridis(10),
          Rowv = F, Colv = T,
          sepcolor = 'black',
          trace = "none",
          scale = 'row',
          key = TRUE,
          key.title = "Expression",
          key.xlab = "Z-score",
          dendrogram = "column",
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap of top 100 DEGs\nMale vs Female\n(N=465)",
          na.color = 'black',
          ColSideColors = ccodes,
          labRow = NA,
          labCol = NA,
          margins = c(11,10))
legend("left",                       
       legend = c("Female", "Male"),     
       col = c("red", "blue"),           
       lty = 1,                          
       lwd = 4,                          
       cex = 0.55,
       xpd = TRUE,
       title = "Gender",
       inset = c(-0.5, 0)) 

# functional enrichment analysis ####
# selection of up- and down-regulated genes from the DEA (1.5 times of difference)
upreg.genes <- rownames(subset(DEA, logFC > 1 & FDR < 0.01)) # 205
dnreg.genes <- rownames(subset(DEA, logFC < -1 & FDR < 0.01)) # 122

# convert ensemble IDs to gene IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
upreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = upreg.genes,
                     mart = mart)$hgnc_symbol
dnreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = dnreg.genes,
                     mart = mart)$hgnc_symbol

# EA for up- and down-regulated genes 
up.EA <- TCGAanalyze_EAcomplete(TFname = "Upregulated", upreg.genes) 
dn.EA <- TCGAanalyze_EAcomplete(TFname = "Downregulated", dnreg.genes)

# barplot for enriched pathways in up-regulated genes 
TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP), 
                        GOBPTab = up.EA$ResBP, 
                        GOCCTab = up.EA$ResCC,
                        GOMFTab = up.EA$ResMF,
                        PathTab = up.EA$ResPat, 
                        nRGTab = upreg.genes, 
                        nBar = 10, 
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

# barplot for enriched pathways in down-regulated genes
TCGAvisualize_EAbarplot(tf = rownames(dn.EA$ResBP), 
                        GOBPTab = dn.EA$ResBP, 
                        GOCCTab = dn.EA$ResCC,
                        GOMFTab = dn.EA$ResMF,
                        PathTab = dn.EA$ResPat, 
                        nRGTab = dnreg.genes, 
                        nBar = 10, 
                        text.size = 2, 
                        fig.width = 30,
                        fig.height = 15)


# survival analysis ####
metadata_surv <- metadata
metadata_surv$event <- ifelse(metadata_surv$`Overall Survival Status` == "Dead", 1, 0) # dead = 1, alive = 0
metadata_surv$time <- metadata_surv$`Overall Survival (Months)`
metadata_surv <- metadata_surv %>% 
  dplyr::select(-days_to_last_follow_up, -days_to_death, -prior_treatment, -`Overall Survival Status`, -`Overall Survival (Months)`)
metadata_surv <- metadata_surv %>% 
  rename(prior_malignancy = `Prior Malignancy`, tumor_stage = `Tumor Stage`, primary_site = `Primary Site`, age = `Age (Years)`, race = Race, gender = Gender)

# code the clinical covariate for easier interpretation 
# gender M-F
metadata_surv$gender <- as.factor(metadata_surv$gender) # ref = female

# race white-black/AA-other
metadata_surv$race <- ifelse(metadata_surv$race == "White", "White",
                    ifelse(metadata_surv$race == "Black or African American", "Black/AA", "Other"))
metadata_surv$race <- as.factor(metadata_surv$race)
metadata_surv$race <- relevel(metadata_surv$race, ref = "White") # ref = white

# age (<=50, >50)
metadata_surv$age <- ifelse(metadata_surv$age <= 50, 0, 1)
metadata_surv$age <- as.factor(metadata_surv$age) # ref = 0 so age <= 50

# prior malignancy
metadata_surv$prior_malignancy <- as.factor(metadata_surv$prior_malignancy) # ref = no

# tumor stage
# meta$tumor_stage <- as.character(meta$tumor_stage)
metadata_surv$tumor_stage <- as.factor(metadata_surv$tumor_stage)
metadata_surv$tumor_stage <- relevel(metadata_surv$tumor_stage, ref = "Stage I") # ref = stage I

# primary site
metadata_surv$primary_site <- ifelse(metadata_surv$primary_site %in% c("Ascending colon", "Cecum", "Hepatic flexure of colon"), "Right Colon (Proximal)",
                            ifelse(metadata_surv$primary_site == "Transverse colon", "Transverse colon", "Left Colon (Distal)"))
metadata_surv$primary_site <- as.factor(metadata_surv$primary_site)
metadata_surv$primary_site <- relevel(metadata_surv$primary_site, ref = "Right Colon (Proximal)") # ref = right colon (proximal)

# association between gender and outcome 
fit.pl <- prodlim(Hist(time, event) ~ gender, data = metadata_surv)

model <- coxph(Surv(time, event) ~ gender, data = metadata_surv)
summary(model) # association with gender: in male, the hazard of outcome is 1.11 times higher compared to female (p=0.6)

plot(fit.pl, 
     lwd = 3,             
     legend.x = 0.8, legend.y = 0.3, 
     legend.cex = 1.2, 
     legend.title = "",
     xlab = "Time (days)",  
     ylab = "Survival Probability")
title("Kaplan-Meier Survival Curve: Male vs Female")

# add HR and p-value to the plot
hr <- round(exp(coef(model)), 2)  
ci <- round(exp(confint(model)), 2)
pval <- summary(model)$sctest[3]

text(x = max(metadata_surv$time) * 1, y = 0.9, 
     labels = paste0("HR = ", hr, " (", ci[1], "-", ci[2], ")\n",
                     "p = ", signif(pval, 3)),
     cex = 1.2, pos = 2, col = "black")


# association between DEGs and outcome 
# univariate analysis with all 259 DEGs identified
# subset genes expect XIST and genes on the Y chromosome
dim(DEGs_df_filtered) # 259 filtered DEGs
View(metadata_surv)

# subset the coad matrix with DEGs_df_filtered
coad_surv <- coad[DEGs_df_filtered$mRNA, ]
dim(coad_surv) # 259 genes, 465 samples 
coad_surv <- as.data.frame(coad_surv)

# use the gene names and not ensembl names
coad_surv_names <- rownames(coad_surv)
annot <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = coad_surv_names,
  mart = mart)
coad_surv$ensembl_gene_id <- rownames(coad_surv)
coad_surv <-  coad_surv %>%
  left_join(annot, by = 'ensembl_gene_id')
rownames(coad_surv) <- coad_surv$hgnc_symbol
coad_surv <- coad_surv %>% 
  dplyr::select(-ensembl_gene_id) 
coad_surv <- coad_surv %>% 
  dplyr::select(-hgnc_symbol) 

# classification for each gene into high and low 
coad_surv <- coad_surv %>%
  apply(1, function(x) ifelse(x > median(x, na.rm = TRUE), "high", "low")) %>%
  t() %>%
  as.data.frame()

# merge coad_surv to metadata_surv
coad_surv <- as.data.frame(t(coad_surv))
coad_surv$barcode <- rownames(coad_surv)
combined_metadata_surv <- merge(metadata_surv, coad_surv, by = "barcode")
View(combined_metadata_surv)


# univariate analysis for all 259 DEGs
# female
female_cohort <- subset(combined_metadata_surv, gender == "Female")
female_cohort <- female_cohort %>% 
  dplyr::select(-gender)
female_cohort <- female_cohort %>% 
  dplyr::select(-full_barcode)
uva_results_female <- data.frame(gene = character(),
                          HR = numeric(),
                          CI_lower = numeric(),
                          CI_upper = numeric(),
                          p_value = numeric(),
                          stringsAsFactors = FALSE)

uva_columns <- setdiff(colnames(female_cohort), c("barcode", "time", "event", "race", "age", "prior_malignancy", "tumor_stage", "primary_site")) 

for (gene in uva_columns) { 
  cox_model <- coxph(Surv(time, event) ~ female_cohort[[gene]], data = female_cohort)
  summary_model <- summary(cox_model)
  
  HR <- summary_model$coefficients[1, "exp(coef)"]
  CI_lower <- summary_model$conf.int[1, "lower .95"]
  CI_upper <- summary_model$conf.int[1, "upper .95"]
  p_value <- summary_model$coefficients[1, "Pr(>|z|)"]
  
  uva_results_female <- rbind(uva_results_female, data.frame(gene = gene, 
                                                             HR = HR, 
                                                             CI_lower = CI_lower, 
                                                             CI_upper = CI_upper, 
                                                             p_value = p_value))
}

# use BH correction to get FDR
uva_results_female$FDR <- p.adjust(uva_results_female$p_value, method = "BH")

# filter out genes with a p-value < 0.05
uva_results_female_filtered <- uva_results_female %>% 
  filter(p_value < 0.05)

dim(uva_results_female_filtered) # 9 genes
significant_genes_female <- uva_results_female$gene[uva_results_female$p_value < 0.05]


# multivariate model
# subset female_cohort with the genes of interest
female_cohort_mva <- female_cohort[, c("time", "event", "age", "race", "prior_malignancy", 
                                       "tumor_stage", "primary_site", significant_genes_female)]

multi_model_female <- coxph(Surv(time, event) ~ KLK5 + SPOCK3 + PWRN1 + LHX3 + ADAD2 + SCN1A + KRT38 + FGF23 + PNMA5 + age + race + prior_malignancy + tumor_stage + primary_site, data = female_cohort_mva)
summary(multi_model_female)

# backward selection
multi_model_female_1 <- coxph(Surv(time, event) ~ KLK5 + PWRN1 + LHX3 + ADAD2 + KRT38 + age + race + prior_malignancy + tumor_stage + primary_site, data = female_cohort_mva)
summary(multi_model_female_1)

# check PH assumptions for categorical variables
ph_test_female <- cox.zph(multi_model_female_1)
ph_test_female
plot(ph_test_female) # no PH violations

# forest plot
CI <- exp(confint(multi_model_female_1))
plot_multi <- sjPlot::plot_model(multi_model_female_1,
                                 axis.lim = c(min(CI), max(CI)),
                                 auto.label = F,show.values = T,show.p = T,vline.color = "black",title="Hazard Ratio")
plot_multi + theme(plot.title = element_text(hjust = 0.5))

# KM of the significant genes: KLK5, PWRN1, LHX3, ADAD2, KRT38

par(mar = c(4, 4, 2, 1), mfrow = c(1, 3))
# KLK5
fit.pl <- prodlim(Hist(time, event) ~ KLK5, data = female_cohort)
model <- coxph(Surv(time, event) ~ KLK5, data = female_cohort)
summary(model) 

plot(fit.pl, 
     lwd = 3,             
     legend.x = 0, legend.y = 0.15, 
     legend.cex = 1.1, 
     legend.title = "",
     xlab = "Time (days)",  
     ylab = "Survival Probability")
title("Kaplan-Meier Survival Curve: KLK5")

# add HR and p-value to the plot
hr <- round(exp(coef(model)), 2)  
ci <- round(exp(confint(model)), 2)
pval <- summary(model)$sctest[3]

text(x = max(metadata_surv$time) * 1, y = 0.02, 
     labels = paste0("HR = ", hr, " (", ci[1], "-", ci[2], ")\n",
                     "p = ", signif(pval, 3)),
     cex = 1.1, pos = 2, col = "black")

# PWRN1
fit.pl <- prodlim(Hist(time, event) ~ PWRN1, data = female_cohort)

model <- coxph(Surv(time, event) ~ PWRN1, data = female_cohort)
summary(model) 

plot(fit.pl, 
     lwd = 3,             
     legend.x = 0.8, legend.y = 0.15, 
     legend.cex = 1.1, 
     legend.title = "",
     xlab = "Time (days)",  
     ylab = "Survival Probability")
title("Kaplan-Meier Survival Curve: PWRN1")

# add HR and p-value to the plot
hr <- round(exp(coef(model)), 2)  
ci <- round(exp(confint(model)), 2)
pval <- summary(model)$sctest[3]

text(x = max(metadata_surv$time) * 1, y = 0.02, 
     labels = paste0("HR = ", hr, " (", ci[1], "-", ci[2], ")\n",
                     "p = ", signif(pval, 3)),
     cex = 1.1, pos = 2, col = "black")

# LHX3
fit.pl <- prodlim(Hist(time, event) ~ LHX3, data = female_cohort)

model <- coxph(Surv(time, event) ~ LHX3, data = female_cohort)
summary(model) 

plot(fit.pl, 
     lwd = 3,             
     legend.x = 0.8, legend.y = 0.15, 
     legend.cex = 1.1, 
     legend.title = "",
     xlab = "Time (days)",  
     ylab = "Survival Probability")
title("Kaplan-Meier Survival Curve: LHX3")

# add HR and p-value to the plot
hr <- round(exp(coef(model)), 2)  
ci <- round(exp(confint(model)), 2)
pval <- summary(model)$sctest[3]

text(x = max(metadata_surv$time) * 1, y = 0.02, 
     labels = paste0("HR = ", hr, " (", ci[1], "-", ci[2], ")\n",
                     "p = ", signif(pval, 3)),
     cex = 1.1, pos = 2, col = "black")

par(mar = c(4, 4, 2, 1), mfrow = c(1, 2))
# ADAD2
fit.pl <- prodlim(Hist(time, event) ~ ADAD2, data = female_cohort)

model <- coxph(Surv(time, event) ~ ADAD2, data = female_cohort)
summary(model) 

plot(fit.pl, 
     lwd = 3,             
     legend.x = 0.8, legend.y = 0.2, 
     legend.cex = 1.1, 
     legend.title = "",
     xlab = "Time (days)",  
     ylab = "Survival Probability")
title("Kaplan-Meier Survival Curve: ADAD2")

# add HR and p-value to the plot
hr <- round(exp(coef(model)), 2)  
ci <- round(exp(confint(model)), 2)
pval <- summary(model)$sctest[3]

text(x = max(metadata_surv$time) * 1, y = 0.02, 
     labels = paste0("HR = ", hr, " (", ci[1], "-", ci[2], ")\n",
                     "p = ", signif(pval, 3)),
     cex = 1.1, pos = 2, col = "black")

# KRT38
fit.pl <- prodlim(Hist(time, event) ~ KRT38, data = female_cohort)

model <- coxph(Surv(time, event) ~ KRT38, data = female_cohort)
summary(model) 

plot(fit.pl, 
     lwd = 3,             
     legend.x = 0.8, legend.y = 0.2, 
     legend.cex = 1.1, 
     legend.title = "",
     xlab = "Time (days)",  
     ylab = "Survival Probability")
title("Kaplan-Meier Survival Curve: KRT38")

# add HR and p-value to the plot
hr <- round(exp(coef(model)), 2)  
ci <- round(exp(confint(model)), 2)
pval <- summary(model)$sctest[3]

text(x = max(metadata_surv$time) * 1, y = 0.02, 
     labels = paste0("HR = ", hr, " (", ci[1], "-", ci[2], ")\n",
                     "p = ", signif(pval, 3)),
     cex = 1.1, pos = 2, col = "black")


# male
male_cohort <- subset(combined_metadata_surv, gender == "Male")
male_cohort <- male_cohort %>% 
  dplyr::select(-gender)
male_cohort <- male_cohort %>% 
  dplyr::select(-full_barcode)
uva_results_male <- data.frame(gene = character(),
                               HR = numeric(),
                               CI_lower = numeric(),
                               CI_upper = numeric(),
                               p_value = numeric(),
                               stringsAsFactors = FALSE)

uva_columns <- setdiff(colnames(male_cohort), c("barcode", "time", "event", "race", "age", "prior_malignancy", "tumor_stage", "primary_site")) 

for (gene in uva_columns) { 
  cox_model <- coxph(Surv(time, event) ~ male_cohort[[gene]], data = male_cohort)
  summary_model <- summary(cox_model)
  
  HR <- summary_model$coefficients[1, "exp(coef)"]
  CI_lower <- summary_model$conf.int[1, "lower .95"]
  CI_upper <- summary_model$conf.int[1, "upper .95"]
  p_value <- summary_model$coefficients[1, "Pr(>|z|)"]
  
  uva_results_male <- rbind(uva_results_male, data.frame(gene = gene, 
                                                         HR = HR, 
                                                         CI_lower = CI_lower, 
                                                         CI_upper = CI_upper, 
                                                         p_value = p_value))
}

# use BH correction to get FDR
uva_results_male$FDR <- p.adjust(uva_results_male$p_value, method = "BH")

# filter out genes with a p-value < 0.05
uva_results_male_filtered <- uva_results_male %>% 
  filter(p_value < 0.05)

dim(uva_results_male_filtered) # 14 genes
significant_genes_male <- uva_results_male$gene[uva_results_male$p_value < 0.05]


# multivariate model
# subset female_cohort with the genes of interest
male_cohort_mva <- male_cohort[, c("time", "event", "age", "race", "prior_malignancy", 
                                    "tumor_stage", "primary_site", significant_genes_male)]

multi_model_male <- coxph(Surv(time, event) ~ NNAT + KRT16 + TCHH + CDK5R2 + KRT6C + ADAD2 + CHRNA4 + LINC00645 + DSCAM + AIRE + LRRTM1 + SBSN + PTPRN + PAEP + age + race + prior_malignancy + tumor_stage + primary_site, data = male_cohort_mva)
summary(multi_model_male)

# backward selection
multi_model_male_1 <- coxph(Surv(time, event) ~ PAEP + age + race + prior_malignancy + tumor_stage + primary_site, data = male_cohort_mva)
summary(multi_model_male_1)

# check PH assumptions for categorical variables
ph_test_male <- cox.zph(multi_model_male_1)
ph_test_male
plot(ph_test_male) # no PH violations

# forest plot
CI <- exp(confint(multi_model_male_1))
plot_multi <- sjPlot::plot_model(multi_model_male_1,
                                 axis.lim = c(min(CI), max(CI)),
                                 auto.label = F,show.values = T,show.p = T,vline.color = "black",title="Hazard Ratio")
plot_multi + theme(plot.title = element_text(hjust = 0.5))

# KM of the significant genes: PAEP

par(mar = c(4, 4, 2, 1), mfrow = c(1, 1))
# PAEP
fit.pl <- prodlim(Hist(time, event) ~ PAEP, data = male_cohort)
model <- coxph(Surv(time, event) ~ PAEP, data = male_cohort)
summary(model) 

plot(fit.pl, 
     lwd = 3,             
     legend.x = 0.8, legend.y = 0.3, 
     legend.cex = 1.2, 
     legend.title = "",
     xlab = "Time (days)",  
     ylab = "Survival Probability")
title("Kaplan-Meier Survival Curve: PAEP")

# add HR and p-value to the plot
hr <- round(exp(coef(model)), 2)  
ci <- round(exp(confint(model)), 2)
pval <- summary(model)$sctest[3]

text(x = max(metadata_surv$time) * 0.9, y = 0.9, 
     labels = paste0("HR = ", hr, " (", ci[1], "-", ci[2], ")\n",
                     "p = ", signif(pval, 3)),
     cex = 1.2, pos = 2, col = "black")