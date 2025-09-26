library(ggplot2)
library(grid)
library(dplyr)
library(data.table)
library(ggpubr)
library(tidyr)
library(purrr)
library(stats)

## path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
cnv_path <- file.path(sherlock.path,"Data/CNV","cnv_starfish.tsv")
cnv_gistic_path <- file.path(scratch.path,"22-comprehensive-oncoplot/arm-level-and-gene-level-cn/revised-gistic-cn.txt")
snv_path <- file.path(scratch.path,"06-SNV/sample_snv.txt")
sv.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
outpath <- file.path(scratch.path,"22-comprehensive-oncoplot","plot")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

### read
other_sample_infomation=read.delim(sampleinfo.path, header=T, check.names=F, stringsAsFactors=F, na.strings="no_na")
histology <- read.csv(histology.path)
other_sample_infomation <- merge(other_sample_infomation,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")
samplebsid <- other_sample_infomation %>% pull(Tumor_Barcode)
snv  <-  read.delim(snv_path,row.names = 1,check.names = F)
cnv_gistic <- read.delim2(cnv_gistic_path,check.names = FALSE)
sbs4 <- read.delim2(sbs4.path)

### driver genes
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
driver_set <- c("EGFR","TP53","RBM10","KRAS",
                "PIK3CA","SETD2","ARID1A","HUWE1","GNAS",
                "RB1","PTEN","CDKN2A","ARID2","KEAP1",
                "NF1","SMARCA4","ATM","STK11","ARID1B",
                "PDGFRA","ALK")

cnv_set <- c("CDKN2A","TERT","MDM2","EGFR","PTPRD",
             "MYC","CTNND2","NKX2-1","PTK6","MCL1",
             "FGFR1","CCNE1","STK11","MET","BRAF")
loss_set <- c("CDKN2A","PTPRD","STK11")
gain_set <- cnv_set[cnv_set %in% loss_set==FALSE]
gistic_values <- paste0("gistic.",cnv_set)
loss_values <- paste0("gistic.",loss_set)
gain_values <- paste0("gistic.",gain_set)

### driver genes
driver_df <- snv[,driver_set] %>% mutate(across(everything(), ~ na_if(., "None")))
### gistic genes
gistic_df <- cnv_gistic[,c("SampleID",cnv_set)]
# merge
all_df <- other_sample_infomation
all_df <- merge(all_df,driver_df %>% mutate(sample=rownames(.)),by.x="Tumor_Barcode",by.y="sample",all.x=T) # add driver
all_df <- merge(all_df,gistic_df %>% setnames(.,cnv_set,paste0("gistic.",cnv_set)),by.x="Tumor_Barcode",by.y="SampleID",all.x=T) # add cnv
all_df <- merge(all_df,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
all_df <- all_df %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
all_df <- all_df %>%
  mutate(across(all_of(gain_values), ~ ifelse(.x ==3,"gain",NA))) %>%
  mutate(across(all_of(loss_values), ~ ifelse(.x ==-3,"loss",NA))) %>%
  arrange_at(c(driver_set,loss_values,gain_values))
# simplize all_df
all_df_simple <- all_df %>%
  mutate(across(all_of(driver_set), ~ ifelse(is.na(.x),"WT","Mut"))) %>%
  mutate(across(all_of(gistic_values), ~ ifelse(is.na(.x),"WT","Mut")))


# Prepare data
filter_df_simple <- all_df_simple %>%
  filter(Smoking_SBS4 == "Non-Smoker-noSBS4") %>%
  filter(Histology == "Adenocarcinoma") 

EGFR.total.mut <- table(filter_df_simple$EGFR)["Mut"] %>% as.numeric()
EGFR.total.wt <- table(filter_df_simple$EGFR)["WT"] %>% as.numeric()

plot_df <- filter_df_simple %>%
  dplyr::select(Tumor_Barcode, EGFR, driver_set,gistic_values) %>%
  unnest(c(driver_set,gistic_values)) %>%
  pivot_longer(-c(Tumor_Barcode, EGFR), names_to = "Gene", values_to = "Status") %>%
  mutate(EGFR_status = ifelse(EGFR == "Mut", "EGFR_mut", "EGFR_wt")) %>%
  group_by(Gene, EGFR_status) %>%
  summarise(n_altered = sum(Status == "Mut"), .groups = "drop") %>%
  pivot_wider(names_from = EGFR_status, values_from = n_altered, values_fill = 0) %>%
  mutate(EGFR_wt = -EGFR_wt/EGFR.total.wt) %>%  # flip direction for wildtype
  mutate(EGFR_mut = EGFR_mut/EGFR.total.mut) %>%  # flip direction for wildtype
  mutate(type = ifelse(Gene %in% driver_set,1,
                        ifelse(Gene %in% loss_values,2,3))) %>%
  arrange(type,desc(EGFR_mut))
gene_order <- plot_df$Gene
gene_order <- c("TP53","RBM10","KRAS","PIK3CA","SETD2",
                "ARID1A","HUWE1","GNAS","RB1","ARID2",
                "PTEN","ATM","NF1","CDKN2A","SMARCA4",
                "ARID1B","PDGFRA","ALK","KEAP1","STK11",
                "gistic.CDKN2A","gistic.PTPRD","gistic.STK11",
                "gistic.TERT","gistic.MDM2","gistic.MYC","gistic.NKX2-1","gistic.PTK6",
                "gistic.EGFR","gistic.CTNND2","gistic.MCL1","gistic.CCNE1","gistic.BRAF",
                "gistic.MET","gistic.FGFR1")

# Get sample-level long format
fisher_input_df <- filter_df_simple %>%
  dplyr::select(Tumor_Barcode, EGFR, driver_set, gistic_values) %>%
  unnest(c(driver_set, gistic_values)) %>%
  pivot_longer(-c(Tumor_Barcode, EGFR), names_to = "Gene", values_to = "Status") %>%
  filter(Gene != "EGFR") %>%
  mutate(EGFR_status = ifelse(EGFR == "Mut", "EGFR_mut", "EGFR_wt"),
         Altered = ifelse(Status == "Mut", 1, 0))
# Fisher's test per gene
fisher_results <- fisher_input_df %>%
  group_by(Gene) %>%
  summarise(
    mut_alt = sum(EGFR_status == "EGFR_mut" & Altered == 1),
    mut_wt  = sum(EGFR_status == "EGFR_mut" & Altered == 0),
    wt_alt  = sum(EGFR_status == "EGFR_wt" & Altered == 1),
    wt_wt   = sum(EGFR_status == "EGFR_wt" & Altered == 0),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    pval = fisher.test(matrix(c(mut_alt, mut_wt, wt_alt, wt_wt), nrow = 2))$p.value,
    or = fisher.test(matrix(c(mut_alt, mut_wt, wt_alt, wt_wt), nrow = 2))$estimate
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(pval, method = "fdr"),
         label = case_when(
           fdr < 0.001 ~ "***",
           fdr < 0.01 ~ "**",
           fdr < 0.05 ~ "*",
           TRUE ~ ""
         ))

plot_df_annotated <- plot_df %>%
  left_join(fisher_results %>% dplyr::select(Gene, pval, or,fdr, label), by = "Gene")
plot_df_annotated <- plot_df_annotated %>% 
  mutate(Gene=factor(Gene,levels=gene_order)) %>% 
  mutate(y_label_pos = ifelse(or > 1, EGFR_mut + 0.02, EGFR_wt - 0.02))
# Plot
ggplot(plot_df_annotated, aes(x = Gene)) +
  geom_bar(aes(y = EGFR_mut), stat = "identity", fill = "black",color = "black",width = 0.7) +
  geom_bar(aes(y = EGFR_wt), stat = "identity", fill = "white", color = "black",width = 0.7) +
  geom_text(aes(y = y_label_pos, label = label), size = 2, vjust = 1, angle=90) +
  theme_bw() +
  labs(x = "Gene", y = "Percentage") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=6, colour = "black", angle = 90, hjust = 0, vjust = 0.5),
    axis.text.y = element_text(size=6, colour = "black"),
    axis.title.y = element_text(size=6, colour = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(size=0.1, colour = "black"),
    axis.line.x.top = element_line(size=0.1, colour = "black"),
    axis.line.x.bottom = element_line(size=0.1, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.234),
    axis.ticks.length.y = unit(0.8, "mm"),
    axis.ticks.y = element_line(linewidth = 0.234),
    plot.title = element_text(size=6, colour = "black"),
    panel.background = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.title.x = element_text(size=6, angle = 0, colour = "black", hjust = 0.5, vjust = 0.5),
    strip.title.y = element_text(size=6, angle = 0, colour = "black", hjust = 1, vjust = 0),
    strip.text.x = element_text(size=6, angle = 0, colour = "black", hjust = 0.5, vjust = 0.5),
    strip.text.y = element_text(size=6, angle = 0, colour = "black", hjust = 1, vjust = 0),
    strip.background = element_blank(),
    strip.placement = "top",
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.title = element_text(size=6,color = "black"),
    legend.text = element_text(size=6,color = "black"),
    legend.position = "bottom",
    plot.margin = grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(outpath,"luad.consistent.never.smoker.egfr.cooccurence.pdf"),width = 5,height = 3)
