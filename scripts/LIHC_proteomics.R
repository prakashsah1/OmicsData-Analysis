##===============================================================================
## Recreate Figure 3 and related supplementary figures from Gao et al., Cell 2019
##===============================================================================

# Load packages
library(data.table)
library(readxl)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(clusterProfiler)
library(ggplot2)
library(survival)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ===================
# Recreate Figure 3A 
# ===================

# 1. Read input files (files downloaded from NCI proteomics data commons)
tsv_path  <- "/Users/prakashsah/Github/LIHC_Proteomics/Data/Zhou_Liver_Cancer_Proteome.tmt11.tsv"
xlsx_path <- "/Users/prakashsah/Github/LIHC_Proteomics/Gao et al 2019/mmc3.xlsx"

# Load proteome matrix (log2 sample/reference values)
prot <- fread(tsv_path, sep="\t", header=TRUE)

# Remove summary rows
prot <- prot %>% filter(!(Gene %in% c("Mean","Median","StdDev")))
# Remove last 6 annotation columns
prot <- prot[, -( (ncol(prot)-5):ncol(prot) )]

# 2. Extract DE gene list (1,274 proteins from Sheet 2 of mmc3.xlsx)
de_genes <- read_excel(xlsx_path, sheet = 2)   # sheet 2 = DE protein list
gene_col <- grep("Gene|Symbol", names(de_genes), value=TRUE)[1]
de_list  <- unique(de_genes[[gene_col]])

# 3. Subset to tumor samples only
# Tumor sample columns start with "Txxx Log Ratio" or "Txxx Unshared Log Ratio"
tumor_cols <- grep("^T\\d+", names(prot), value=TRUE)
expr <- prot %>%
  filter(Gene %in% de_list) %>%
  select(Gene, all_of(tumor_cols))

# Subset to tumor samples (Unshared Log Ratios only)
tumor_cols <- grep("^T\\d+ Unshared Log Ratio$", names(prot), value=TRUE)
expr <- prot %>%
  filter(Gene %in% de_list) %>%
  select(Gene, all_of(tumor_cols))

# Subset to tumor samples (Log Ratios only, not Unshared)
tumor_cols <- grep("^T\\d+ Log Ratio$", names(prot), value = TRUE)

expr <- prot %>%
  filter(Gene %in% de_list) %>%
  select(Gene, all_of(tumor_cols))
# Note: The Txxx Log Ratio (shared) was used for this analysis. However, codes for extracting both together and separately are given above.  

# Set rownames to gene
expr_mat <- as.data.frame(expr)
rownames(expr_mat) <- expr_mat$Gene
expr_mat <- expr_mat[,-1]
expr_mat <- as.matrix(expr_mat)

#Impute any missing values
anyNA(expr_mat)
imputed_matrix <- impute.knn(expr_mat)$data
anyNA(imputed_matrix)

# 4. Row-wise Z-score (per protein)
z_expr <- t(scale(t(imputed_matrix), center=TRUE, scale=TRUE))

# 5. Consensus clustering of tumor samples
set.seed(123)
cc_res <- ConsensusClusterPlus(
  z_expr,
  maxK=6,
  reps=1000,
  pItem=0.8,
  pFeature=0.8,
  clusterAlg="km",
  distance="euclidean",
  title="Consensus_Proteome",
  plot="pdf"
)

# Choose k=3 as per paper
subgroups <- cc_res[[3]]$consensusClass
subgroups <- factor(subgroups, labels=c("S-Mb","S-Me","S-Pf"))
# order matrix and annotation by proteome subgroup
sample_order <- order(subgroups)
z_expr_ordered <- z_expr[, sample_order]

ann_col <- data.frame(Subgroup=subgroups[sample_order])
rownames(ann_col) <- colnames(z_expr_ordered)


# 6. Heatmap (Figure 3A)
palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)

pheatmap(
  z_expr_ordered,
  cluster_rows=TRUE,
  cluster_cols=FALSE,       # don’t re-cluster samples
  annotation_col=ann_col,
  show_colnames=FALSE,
  show_rownames=FALSE,
  color=palette
)


# ==================
# Recreate Figure 3B
# ==================

# 1. Read clinical data (downloaded from NCI proteomics data commons)
clinical_data1 = read_excel("/Users/prakashsah/Github/LIHC_Proteomics/Data/S049_HCC_Clinical_Information_and_TMT11_Sample_Mapping_Gao2019.r1.xlsx", sheet = "Table 1. Clinical Information")
clinical_data2 = read_excel("/Users/prakashsah/Github/LIHC_Proteomics/Data/S049_HCC_Clinical_Information_and_TMT11_Sample_Mapping_Gao2019.r1.xlsx", sheet = "Table 2. Follow-up -PHI Removed")

# 2. create a mapping and annotation data frame
map_df <- clinical_data1 %>%
select(`Tumor ID in Proteome experiment`, `Case ID`) %>%
rename(TumorID = `Tumor ID in Proteome experiment`,
CaseID = `Case ID`)

ann_col1 = data.frame(Subgroup=subgroups)
ann_col1$TumorID <- sub("^T([0-9]+).*", "\\1", rownames(ann_col1))

# clinical data only 160 entries and the notes say 5 patients had their data removed. 
# T724 is mistakenly named. it should T723 (the paired normal sample is 724)
setdiff(ann_col1$TumorID, map_df$TumorID)
ann_col1 <- ann_col1[!ann_col1$TumorID %in% c("213", "351", "447", "535", "697"), ]
ann_col1$TumorID[ann_col1$TumorID == "724"] <- "723"
rownames(ann_col1)[rownames(ann_col1) == "T724 Log Ratio"] <- "T723 Log Ratio"
setdiff(ann_col1$TumorID, map_df$TumorID) # check again

# Ensure TumorID is character in both data frames
map_df <- map_df %>% mutate(TumorID = as.character(TumorID))
ann_col1 <- ann_col1 %>% mutate(TumorID = as.character(TumorID))

# 3. Merge clinical follow-up with mapping (CaseID) and subgroup assignment
survival_df <- clinical_data2 %>%
rename(CaseID = `Case ID`,
OS_time = `Survival time (month)`,
OS_event = `Survival Test`,
DFS_time = `The time to relapse (month)`,
DFS_event = `DFS test (disease free survival)`) %>%
left_join(map_df, by="CaseID") %>%
left_join(ann_col1, by="TumorID")
# Final check
head(survival_df)


# 4. survival analysis and plot
survival_df$Subgroup <- factor(survival_df$Subgroup, levels = c("S-Mb", "S-Pf", "S-Me")) #Reorder Subgroup levels
levels(survival_df$Subgroup)
fit_OS <- survfit(Surv(OS_time, OS_event) ~ Subgroup, data=survival_df)
ggsurvplot(fit_OS, data=survival_df, pval=TRUE, risk.table=TRUE,
palette=c("red", "blue", "green"))


# =======================================
# Differential analysis for Figure 3E & F
# =======================================

# 1. Select tumor samples only for all genes
expr_tumor <- prot %>% select(Gene, all_of(tumor_cols))
colnames(expr_tumor) <- sub(" Log Ratio$", "", colnames(expr_tumor))
colnames(expr_tumor)[colnames(expr_tumor) == "T724"] <- "T723" # The T723 is mislabelled as T724 

# 2. Get thrombus status
thrombus <- clinical_data1 %>%
  select(`Tumor ID in Proteome experiment`, `Tumor thrombus`) %>%
  rename(TumorID = `Tumor ID in Proteome experiment`,
         Thrombus = `Tumor thrombus`) %>%
  mutate(TumorID = paste0("T", TumorID))   # add "T" prefix

# 3. Match samples between expr_tumor and thrombus
common <- intersect(colnames(expr_tumor)[-1], thrombus$TumorID)

expr_tumor <- expr_tumor %>% select(Gene, all_of(common))
thrombus <- thrombus %>% filter(TumorID %in% common)

grp <- thrombus$Thrombus    # 0 = no thrombus, 1 = thrombus
names(grp) <- thrombus$TumorID

cat("Number of tumor samples in proteome:", ncol(expr_tumor)-1, "\n")
cat("Number of tumor samples with thrombus info:", length(grp), "\n")

# 4. Build matrix for tumor samples
expr_tumor_mat <- as.data.frame(expr_tumor)
rownames(expr_tumor_mat) <- expr_tumor_mat$Gene
expr_tumor_mat <- expr_tumor_mat[,-1]
expr_tumor_mat <- as.matrix(expr_tumor_mat)

expr_filt <- expr_tumor_mat[rowMeans(is.na(expr_tumor_mat)) <= 0.5, ] # Filter proteins with >50% missing values
cat("Proteins after filtering (>50% missing removed):", nrow(expr_filt), "\n")

expr_imputed <- impute.knn(expr_filt)$data #k-NN imputation for missing values
grp <- grp[colnames(expr_imputed)] # Align thrombus labels with columns

# 5. Differential test (two-sample t-test per protein)
t_tests <- apply(expr_imputed, 1, function(x) {
  yes <- x[grp == 1]
  no  <- x[grp == 0]
  test <- t.test(yes, no)   # Welch’s t-test
  c(mean_yes = mean(yes),
    mean_no  = mean(no),
    logFC    = mean(yes) - mean(no),
    pval     = test$p.value)
})

res_df <- as.data.frame(t(t_tests))
res_df$qval <- p.adjust(res_df$pval, method="BH")
res_df$Gene <- rownames(res_df)

# 5. Get significant proteins
sig_df <- res_df %>% filter(qval < 0.1) 

DE_thrombus_mat <- expr_imputed[sig_df$Gene, , drop = FALSE] # subset matrix for significant genes

# 6. create annotation data frame with both thrombus status and proteome subgroups
ann_col3 <- merge(
  ann_col2[, c("TumorID", "Subgroup")],
  thrombus[, c("TumorID", "Thrombus")],
  by = "TumorID"
) #ann_col2 was created for figure S7A first in this analysis(see below).
ann_col3 <- ann_col3[match(colnames(DE_thrombus_mat), ann_col3$TumorID)] # match with matrix order
rownames(ann_col3) = ann_col3$TumorID
identical(colnames(DE_thrombus_mat), rownames(ann_col3))# should return TRUE

#order matrix by thrombus status
ord1 <- order(ann_col3$Thrombus)
DE_thrombus_mat <- DE_thrombus_mat[, ord1]
ann_col3 <- ann_col3[ord1, ]

ann_col3$Thrombus <- factor(ann_col3$Thrombus, levels = c("0","1")) # Make sure thrombus is categorical

# set annotation colors manually to match other figures
ann_colors <- list(
  Thrombus = c("0" = "grey70", "1" = "black"),
  Subgroup = c("S-Mb" = "green", 
               "S-Me" = "salmon", 
               "S-Pf" = "skyblue")
)

# 7. plot heatmap (Figure 3E)
pheatmap(
  t(scale(t(DE_thrombus_mat))),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  show_colnames = FALSE,
  show_rownames = FALSE,
  annotation_col = ann_col3[, c("Thrombus", "Subgroup")],
  annotation_colors = ann_colors,
  color = palette
  )

# 8. Create list of downregulated genes in tumors with thrombus
down_genes_thrombus <- sig_df$Gene[sig_df$logFC < 0]

#clean up down regulated genes
down_genes_clean <- down_genes_thrombus
down_genes_clean <- setdiff(down_genes_clean, c("ATP5MF-PTCD1", "ABHD14A-ACY1")) # Remove the read-through genes
down_genes_clean <- c(down_genes_clean, "ATP5MF", "PTCD1", "ABHD14A", "ACY1") # Add their parent genes
down_genes_clean <- unique(down_genes_clean) # Remove any accidental duplicates

# 9. perform GO enrichment analysis using clusterprofiler
thrombus_BP <- enrichGO(
gene          = down_genes_clean,
OrgDb         = org.Hs.eg.db,
keyType       = "SYMBOL",
ont           = "BP", 
pAdjustMethod = "BH",
pvalueCutoff  = 0.05
)
# 10. plot top 30 enriched GO terms (Figure 3F)
barplot(
thrombus_BP,
showCategory = 30,
title        = "GO Biological Processes (Down-regulated)",
font.size    = 10
)

#save the enriched GO terms as csv
write.csv(as.data.frame(thrombus_BP), "GO_BP_results.csv", row.names = FALSE)


# ====================
# Recreate Figure S7A
# ====================

# 1. Create matrix for HCC relevant genes
HCC_relevant_genes <- scan(pipe("pbpaste"), what = character()) # make sure to copy the list of 9 HCC relevant genes in plain text format. 
HCC_rel_mat <- expr_imputed[HCC_relevant_genes, , drop=FALSE]
HCC_rel_mat

# 2. create annotation data frame 
ann_col2 <- data.frame(Subgroup = subgroups)
ann_col2$TumorID <- sub("^(T[0-9]+).*", "\\1", rownames(ann_col2)) # Keep TumorID as "T123"
ann_col2 <- ann_col2[!ann_col2$TumorID %in% c("T213", "T351", "T447", "T535", "T697"), ]
ann_col2$TumorID[ann_col2$TumorID == "T724"] <- "T723"
rownames(ann_col2)[rownames(ann_col2) == "T724 Log Ratio"] <- "T723 Log Ratio"

common <- intersect(colnames(HCC_rel_mat), ann_col2$TumorID)

HCC_rel_mat <- HCC_rel_mat[, common]
ann_col2 <- ann_col2[match(common, ann_col2$TumorID), , drop=FALSE]

# 3. order matrix and annotation by proteome subgroup
ord2 <- order(ann_col2$Subgroup)
HCC_rel_mat <- HCC_rel_mat[, ord2]
ann_col2 <- ann_col2[ord2, , drop=FALSE]

# 4. plot heatmap for HCC relevant genes (Figure S7A)
pheatmap(
  (t(scale(t(HCC_rel_mat)))),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = data.frame(Subgroup = ann_col2$Subgroup,
                              row.names = ann_col2$TumorID),
  show_colnames = FALSE,
  show_rownames = TRUE,
  color = palette
)


# ===================
# Recreate Figure S7B
# ===================

# 1. Create list of clinically relevant HCC biomarkers
HCC_biomarkers <- scan(pipe("pbpaste"), what = character()) # make sure to copy the list of HCC relevant biomarkers in plain text format. 
HCC_biomarkers

# 2. create matrix for heatmap (Keep only "Log Ratio" columns for each sample)
prot_log_only <- prot %>%
  select(Gene, contains("Log Ratio")) %>%
  select(-contains("Unshared"))

prot_log_only <- prot_log_only %>% column_to_rownames("Gene") #Set Gene as rownames
colnames(prot_log_only) <- gsub(" Log Ratio", "", colnames(prot_log_only)) # Trim sample IDs

HCC_bm_mat <- as.matrix(prot_log_only[rownames(prot_log_only) %in% HCC_biomarkers, ])
sum(is.na(HCC_bm_mat)) # check any missing values
HCC_bm_mat <- impute.knn(HCC_bm_mat)$data # impute missing values

# 3. create annotation data frame
SampleID <- gsub(" Log Ratio", "", colnames(prot_log_only))
Sample_type <- ifelse(grepl("^T", SampleID), "Tumor", "Normal")
subgroup_info <- ann_col
rownames(subgroup_info) <- gsub(" Log Ratio$", "", rownames(subgroup_info))
subgroup_info$TumorID <- rownames(subgroup_info)
subgroup <- subgroup_info$Subgroup[match(SampleID, subgroup_info$TumorID)]
ann_col4 <- data.frame(
  SampleType = Sample_type,
  Subgroup = subgroup
)
rownames(ann_col4) <- SampleID

# 4. Order samples by Normal followed by tumors (by proteome subgroups)
Sample_order <- order(ann_col4$SampleType, ann_col4$Subgroup)

ann_col4 <- ann_col4[Sample_order, , drop = FALSE]
HCC_bm_mat  <- HCC_bm_mat[, Sample_order]
identical(colnames(HCC_bm_mat), rownames(ann_col4)) #check sample order

# 5. plot heatmap (Figure S7B)
summary(as.vector(t(scale(t(HCC_bm_mat)))))
# since heatmap color scale is being pulled by few extreme values, the z scores are capped.
HCC_bm_mat_z <- t(scale(t(HCC_bm_mat)))
HCC_bm_mat_z[HCC_bm_mat_z > 3] <- 3
HCC_bm_mat_z[HCC_bm_mat_z < -3] <- -3

breaks <- seq(-3, 3, length.out = length(palette) + 1) # Define breaks that match the palette length
pheatmap(
  HCC_bm_mat_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col4,
  show_colnames = FALSE,
  annotation_colors = ann_colors,
  color = palette, 
  breaks = breaks
)

# 6. Barplot of logFC for HCC relevant biomarkers
tumor_samples <- grep("^T", colnames(HCC_bm_mat))
normal_samples <- grep("^P", colnames(HCC_bm_mat))

tumor_means  <- rowMeans(HCC_bm_mat[, tumor_samples], na.rm = TRUE)
normal_means <- rowMeans(HCC_bm_mat[, normal_samples], na.rm = TRUE)

log2FC <- tumor_means - normal_means   # difference of log means
df_log2FC <- data.frame(Gene = rownames(HCC_bm_mat), log2FC = log2FC)

ggplot(df_log2FC, aes(x = Gene, y = log2FC)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, color = "black") + # axis line through 0
  coord_cartesian(ylim = c(-1.5, 1.5)) +        # cap range
  scale_x_discrete(limits = rev(sort(unique(df_log2FC$Gene)))) + 
  labs(
    x = NULL,
    y = "log2(Tumor / Normal)",
    title = "HCC relevant biomarkers"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

# ===================
# Recreate Figure S7C
# ===================

# For a list of potential drug targets, the authors did the following:
  # 1. list of proteins upregulated in tumors based on proteomics, 2. drug targets file from Drugbank database, 3. intersect the list of upregulated proteins and drug targets to get list of potential drug targets.
  # An account on the Drug database is needed to access drug target data. 
  # To perform what authors did, following code (commented out) can be used.
    #de_up <- de_genes %>% filter(`Log2 FC` > 0)
    #drug_targets <- read_csv("filname of the drug targets from Drungbank") #use appropriate function based on file type
    #drug_targets_genes, extract gene name/symbols for the approved/under clinical trial targets
    #HCC_drug_targets <- intersect(de_up$`Gene symbol`, drug_targets_genes)

# 1. Create the list of drug targets (list of drug targets from the plot in the paper was used)
HCC_drug_targets <- c(
  "PLA2G2A", "IMPDH2", "ATIC", "TYMS", "PCSK9", "SOAT1", "G6PD", "NPC1L1", "CPD",
  "PYCR1", "PYCR2", "BCAT1", "ACSL3", "ATF1", "TOP1", "TOP2A", "TOP2B", "POLE3",
  "RRM2", "POLE4", "HDAC1", "HDAC2", "SFRP4", "RRM1", "CDK1", "CDK2", "MAPK13",
  "LIG3", "HDAC3", "PRKDC", "PARP1", "CHEK2", "APEX1", "TK1", "P2RX4", "MMP14",
  "STK24", "VCAN", "STK39", "COQ8B", "ATOX1", "DDX5", "P4HA2", "STIP1")
HCC_drug_targets

# 2. create matrix for drug targets
HCC_dt_mat <- as.matrix(prot_log_only[rownames(prot_log_only) %in% HCC_drug_targets, ])
sum(is.na(HCC_dt_mat)) # check any missing values
HCC_dt_mat <- impute.knn(HCC_dt_mat)$data # impute missing values

HCC_dt_mat <- HCC_dt_mat[, Sample_order] # order samples in same way as for figure S7 B.
identical(colnames(HCC_dt_mat), rownames(ann_col4))

# 3. plot heatmap (Figure S7C)
summary(as.vector(t(scale(t(HCC_bm_mat)))))
# since heatmap color scale is being pulled by few extreme values, the z scores are capped.
HCC_dt_mat_z <- t(scale(t(HCC_dt_mat)))
HCC_dt_mat_z[HCC_bm_mat_z > 3] <- 3
HCC_dt_mat_z[HCC_bm_mat_z < -3] <- -3

breaks <- seq(-3, 3, length.out = length(palette) + 1) # Define breaks that match the palette length
pheatmap(
  HCC_dt_mat_z,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = ann_col4,
  show_colnames = FALSE,
  annotation_colors = ann_colors,
  color = palette, 
  breaks = breaks
)

# 4. plot relative intensity of proteins that are potential drug targets
# make a long format table
HCC_dt_long <- HCC_dt_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>% 
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

HCC_dt_long$Subgroup <- ann_col4[HCC_dt_long$Sample, "Subgroup"] # add subgroup info
HCC_dt_long <- HCC_dt_long %>% filter(!is.na(Subgroup)) # drop normal samples

# Convert to linear intensity
HCC_dt_long <- HCC_dt_long %>%
  mutate(Intensity = 2^Expression)

# calculate mean expression of genes across the proteome subgroups
HCC_dt_meanInt <- HCC_dt_long %>%
  group_by(Gene, Subgroup) %>%
  summarise(MeanInt = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Convert to relative %
HCC_dt_meanInt_relative <- HCC_dt_meanInt %>%
  group_by(Gene) %>%
  mutate(RelativeInt = 100 * MeanInt / sum(MeanInt, na.rm = TRUE)) %>%
  ungroup()
HCC_dt_meanInt_relative$Subgroup <- factor(HCC_dt_meanInt_relative$Subgroup, levels = c("S-Mb", "S-Me", "S-Pf"))

# Plot barplot
ggplot(HCC_dt_meanInt_relative, aes(x = Gene, y = RelativeInt, fill = Subgroup)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  coord_flip() + 
  scale_x_discrete(limits = rev(sort(unique(HCC_dt_meanInt_relative$Gene)))) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("S-Mb" = "green", "S-Me" = "red", "S-Pf" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Relative intensity (%)") +
  xlab("")
