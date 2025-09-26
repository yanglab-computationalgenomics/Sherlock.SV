library(dplyr)
library(data.table)
library(readxl)
library(tidyverse)
library(corrplot)
library(car)

### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
sample_info_path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
sv_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_all_sv_with_signature.tsv")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
kras.G12C.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12C_samples.txt")
kras.G12V.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12V_samples.txt")
kras.G12D.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12D_samples.txt")
kras.Q61H.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pQ61H_samples.txt")
kras.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_other_mutations.txt")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
fusion.detail.path  <- file.path(sherlock.path,"scratch/fusion/end_annotation/tier1_driver_fusions.txt")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
cluster.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_CGR_signature.tsv")
noncluster.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_nonclustered_cmplx_sv.tsv")
simple.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_simple_sv_with_signature.tsv")
cluster_path <- file.path(sherlock.path,paste0("scratch/07-cluster/1217_samples_hg38/Sherlock_SV49_simple/","8","signature_myannotation/","hc_pearson","/sig_cluster","5",".tsv"))
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.suppl.txt")
supple.table1.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/Paper/yangyang_sherlock/Sherlock.SV.draft.20250702/supplementary.tables/table.S1.sample.info.1209.xlsx"
supple.table2.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/Paper/yangyang_sherlock/Sherlock.SV.draft.20250702/supplementary.tables/table.S2.driver.fusions.1209.xlsx"
plot.path <- file.path(scratch.path,"29-vif","plot")
if (dir.exists(plot.path)==FALSE){
  dir.create(plot.path,showWarnings = F,recursive = T)
}

driver_set <- c("EGFR","TP53","RBM10","KRAS",
                "PIK3CA","ERBB2","SETD2","CTNNB1","ARID1A",
                "HUWE1","MET","SMAD4","GNAS","RB1",
                "PTEN","CDKN2A","ARID2","KEAP1","NKX2-1",
                "NF1","SMARCA4","ATM","BRAF","ATRX",
                "STK11","UBA1","ARID1B","PDGFRA","PIK3R1",
                "ALK","TERT","CHD2","NFE2L2","RET",
                "MEN1","AKT1","EIF1AX","MDM2","ROS1",
                "MAP2K1","CCNE1","CDK4","NRAS","RIT1",
                "CUL3","U2AF1","HRAS","CCND1")
# driver_set <- c("EGFR","TP53","RBM10","KRAS",
#                 "PIK3CA","SETD2","ARID1A","HUWE1","GNAS",
#                 "RB1","PTEN","CDKN2A","ARID2","KEAP1",
#                 "NF1","SMARCA4","ATM","STK11","ARID1B",
#                 "PDGFRA","ALK")

cnv_set <- c("CDKN2A","TERT","MDM2","EGFR","PTPRD",
             "MYC","CTNND2","NKX2-1","PTK6","MCL1",
             "FGFR1",
             # "GRIN2A",
             "CCNE1","STK11","MET","BRAF")
loss_set <- c("CDKN2A","PTPRD","STK11")
gain_set <- cnv_set[cnv_set %in% loss_set==FALSE]

fusion_set <- c("ALK","RET","ROS1","MET","NRG1",
                "FGFR1","FGFR2","FGFR3","FGFR4",
                "NTRK1","NTRK2","NTRK3")

cnv_values <- paste0("cnv.",cnv_set)
gistic_values <- paste0("gistic.",cnv_set)
loss_values <- paste0(loss_set,".loss")
gain_values <- paste0(gain_set,".gain")
loss_gistic <- paste0("gistic.",loss_set)
gain_gistic <- paste0("gistic.",gain_set)



### read
supple.table1 <- readxl::read_excel(supple.table1.path)
supple.table2 <- readxl::read_excel(supple.table2.path)
all_sv  <- read.delim(sv_path,stringsAsFactors = FALSE)
sbs4 <- read.delim2(sbs4.path)
fusion.detail <- read.delim2(fusion.detail.path)
fusion.detail <- fusion.detail %>%
  dplyr::mutate(
    Fusion_Gene = case_when(
      GENE1 %in% fusion_set ~ GENE1,
      GENE2 %in% fusion_set ~ GENE2,
      TRUE ~ NA_character_
    )
  ) %>% filter(Known_Fusion=="Yes")
fusion_samples <- fusion.detail$SAMPLE %>% unique()
tp53_mutated_samples=read.delim(tp53.mut.path, header=F, check.names=F, stringsAsFactors=F)
kras_mutated_samples=read.delim(kras.mut.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12c_samples=read.delim(kras.G12C.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12v_samples=read.delim(kras.G12V.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12d_samples=read.delim(kras.G12D.path, header=F, check.names=F, stringsAsFactors=F)
kras_q61h_samples=read.delim(kras.Q61H.path, header=F, check.names=F, stringsAsFactors=F)
kras_others_samples=read.delim(kras.others.path, header=F, check.names=F, stringsAsFactors=F)
all_sample <- read.delim(sample_info_path, header=T, sep="\t", stringsAsFactors = F) 
all_sample <- merge(all_sample,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
all_sample <- all_sample %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
all_sample$TP53 <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(tp53_mutated_samples)),"Mut","WT")
all_sample$KRAS <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(kras_g12c_samples)),"G12C","WT")
all_sample$KRAS <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(kras_g12v_samples)),"G12V",all_sample$KRAS)
all_sample$KRAS <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(kras_g12d_samples)),"G12D",all_sample$KRAS)
all_sample$KRAS <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(kras_q61h_samples)),"Q61H",all_sample$KRAS)
all_sample$KRAS <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(kras_others_samples)),"Others",all_sample$KRAS)
all_sample$EGFR <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(egfr_mutated_E479_A483del_samples)),"Exon 19 del","WT")
all_sample$EGFR <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(egfr_mutated_L591R_samples)),"L858R",all_sample$EGFR)
all_sample$EGFR <- ifelse(all_sample$Tumor_Barcode %in% as.vector(unlist(egfr_mutated_others_samples)),"Others",all_sample$EGFR)
all_sample$Fusion_status <- ifelse(all_sample$Tumor_Barcode %in% fusion_samples,"+","-")
all_sample <- all_sample %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
histology <- read.csv(histology.path)
all_sample <- merge(all_sample,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")


#### consistent smoking
all_sample_consistent_smoking <- all_sample %>% 
  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")) %>%
  filter(Assigned_Population %in% c("EAS","EUR"))
all_sample_consistent_smoking %>% select(Assigned_Population,Smoking_SBS4,purity,wgd,TP53,EGFR,KRAS,Fusion_status)


## ---- Input (edit if your data frame name differs) ----
# Using your example object:
df_raw <- all_sample_consistent_smoking %>%
  dplyr::select(Assigned_Population, Smoking_SBS4, purity, wgd, TP53, EGFR, KRAS, Fusion_status)

## ---- Recode to numeric/binary variables ----
# Smoking: Smoker-SBS4 = 1, Non-Smoker-noSBS4 = 0
# Ancestry (EUR vs other): EUR = 1, other = 0
# WGD: WGD = 1, nWGD = 0
# Gene mutation presence: anything not "WT" (and not NA) = 1, else 0
# Fusion status: "+" = 1, "-" (or others) = 0
df_num <- df_raw %>%
  mutate(
    Smoking = case_when(
      Smoking_SBS4 == "Smoker-SBS4" ~ 1,
      Smoking_SBS4 == "Non-Smoker-noSBS4" ~ 0,
      TRUE ~ NA_real_
    ),
    Ancestry = if_else(Assigned_Population == "EUR", 1, 0, missing = NA_real_),
    WGD = case_when(
      wgd == "WGD" ~ 1,
      wgd == "nWGD" ~ 0,
      TRUE ~ NA_real_
    ),
    TP53mut = if_else(!is.na(TP53) & TP53 != "WT", 1, 0, missing = NA_real_),
    EGFRmut = if_else(!is.na(EGFR) & EGFR != "WT", 1, 0, missing = NA_real_),
    KRASmut = if_else(!is.na(KRAS) & KRAS != "WT", 1, 0, missing = NA_real_),
    Fusion = if_else(Fusion_status == "+", 1, 0, missing = NA_real_)
  ) %>%
  dplyr::select(purity, Smoking, Ancestry, WGD, TP53mut, EGFRmut, KRASmut, Fusion)

## ---- Correlation matrix (pairwise complete obs) ----
corr_mat <- cor(df_num, use = "pairwise.complete.obs", method = "pearson")

# Identify the largest absolute off-diagonal correlation
abs_corr <- abs(corr_mat)
diag(abs_corr) <- NA
max_abs_corr_value <- max(abs_corr, na.rm = TRUE)

# Retrieve the variable pair achieving the maximum absolute correlation
which_max <- which(abs_corr == max_abs_corr_value, arr.ind = TRUE)[1, ]
max_pair <- c(rownames(abs_corr)[which_max[1]], colnames(abs_corr)[which_max[2]])

cat("Largest absolute correlation:", round(max_abs_corr_value, 3),
    "between", max_pair[1], "and", max_pair[2], "\n")

## ---- Correlation plot ----
# Positive correlations in blue, negative in red (corrplot default maps this)
pdf(file.path(plot.path,"correlations.pdf"),height = 3,width = 3)
corrplot(corr_mat,
         method = "circle",
         type = "upper",
         tl.col = "black",    # text color
         tl.cex = 0.6,        # text size
         tl.srt = 90,
         addCoef.col = NA,    # hide numbers
         cl.cex = 0.6,        # color legend size
         number.cex = 0.6,
         col = colorRampPalette(c("blue", "white", "red"))(200))

## Manually add a legend for circle size
# Choose a few representative correlation values
legend_vals <- c(0.2, 0.5, 0.8)

# Compute circle radii used by corrplot (sqrt of abs(correlation))
circle_sizes <- sqrt(legend_vals)

legend("bottomleft",
       legend = paste0("r = ??", legend_vals),
       pt.cex = circle_sizes * 2,  # scale factor for points
       pch = 21,
       col = "black",
       pt.bg = "grey",
       bty = "n",
       title = "Circle size\n= |correlation|",
       cex = 0.6, text.col = "black")
dev.off()

## ---- VIF analysis ----
# Use purity as the numeric response (only continuous variable available here).
# This mirrors the ???compute VIFs for variables in the multivariable linear model??? step.
fit <- lm(purity ~ Smoking + Ancestry + WGD + TP53mut + EGFRmut + KRASmut + Fusion,
          data = df_num)

vifs <- car::vif(fit)
print(vifs)

# Bar plot of VIFs
vif_df <- enframe(vifs, name = "variable", value = "VIF")

ggplot(vif_df, aes(x = reorder(variable, VIF), y = VIF)) +
  geom_col(fill = "black") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", size = 0.3) +
  coord_flip() +
  labs(x = "Variable", y = "Variance Inflation Factor (VIF)") +
  theme_minimal(base_size = 6, base_family = "sans") +  # set base size
  theme(
    axis.text = element_text(size = 6, colour = "black", face = "plain"),
    axis.title = element_text(size = 6, colour = "black", face = "plain"),
    plot.title = element_text(size = 6, colour = "black", face = "plain"),
    legend.text = element_text(size = 6, colour = "black", face = "plain"),
    legend.title = element_text(size = 6, colour = "black", face = "plain")
  )
ggsave(file.path(plot.path,"vif.pdf"),width = 2,height = 3)
