library(dplyr)
library(ggplot2)
library(ggpattern)
library(ComplexHeatmap)
library(purrr)
library(broom)
library(Barnard)
library(DescTools)

mytheme <- theme_bw() +
  theme(
    axis.text.x = element_text(size=6, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
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

### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.suppl.txt")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.main.txt")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
sv.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")

outpath <- file.path(scratch.path,"26-sv-and-cnv","plot")
if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

### read
cnv_set <- c("CDKN2A","TERT","MDM2","EGFR","PTPRD",
             "MYC","CTNND2","NKX2-1","PTK6","MCL1",
             "FGFR1",
             # "GRIN2A",
             "CCNE1","STK11","MET","BRAF")
loss_set <- c("CDKN2A","PTPRD","STK11")
gain_set <- cnv_set[cnv_set %in% loss_set==FALSE]
loss_gistic <- paste0("gistic.",loss_set)
gain_gistic <- paste0("gistic.",gain_set)

sv <- read.delim(sv.path)
sbs4 <- read.delim2(sbs4.path)
histology <- read.csv(histology.path)
oncoprint_table <- read.delim(oncoprint_table_path,check.names = F) %>%
  mutate(across(all_of(gain_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==3,"gain",NA)))) %>%
  mutate(across(all_of(loss_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==-3,"loss",NA))))
oncoprint_table <- merge(oncoprint_table,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
oncoprint_table <- oncoprint_table %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
oncoprint_table <- merge(oncoprint_table,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

samplebsid <- oncoprint_table %>% pull(Tumor_Barcode)
egfr_mutated_E479_A483del_samples <- read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_L591R_samples <- read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_others_samples <- read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)

### 
mygene <- "EGFR" #NKX2-1 EGFR
mysvregion=data.frame(7,55000000,56000000,"EGFR")
# mysvregion=data.frame(14,36000000,37000000,"NKX2-1")

if (mygene %in% names(oncoprint_table) == FALSE){
  mutgene = "EGFR"
} else {
  mutgene = mygene
}
mytable <- oncoprint_table %>% dplyr::select(Tumor_Barcode,Smoking_SBS4,mutgene,paste0("gistic.",mygene))

# sample with SVs around
sv.positive.sample <- sv %>% filter((CHROM1==mysvregion[1,1] & POS1 >= mysvregion[1,2] & POS1 <= mysvregion[1,3])|
                                      (CHROM2==mysvregion[1,1] & POS2 >= mysvregion[1,2] & POS2 <= mysvregion[1,3])) %>%
  pull(SAMPLE) %>% unique()
sv.negative.sample <- samplebsid[!samplebsid %in% sv.positive.sample]

# Step 1: Create the group label
plot_data <- mytable %>%
  mutate(
    mut = ifelse(is.na(.data[[mutgene]]), "mut-", "mut+"),
    svmut = ifelse(Tumor_Barcode %in% sv.positive.sample,"sv+","sv-"),
    cnv = ifelse(.data[[paste0("gistic.", mygene)]] %in% c("gain","loss"), "cnv+", "cnv-")
  ) %>%
  mutate(mut_gain_group = paste0(mut,svmut,cnv))

# Step 2: Summarize counts and percentages
summary_data <- plot_data %>%
  group_by(Smoking_SBS4, mut_gain_group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Smoking_SBS4) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(
    smk_group = ifelse(grepl("Non-Smoker", Smoking_SBS4), "Non-Smoker", "Smoker"),
    mut_group = ifelse(grepl("mut\\+", mut_gain_group), "Mutation+", "Mutation-"),
    sv_group = ifelse(grepl("sv\\+", mut_gain_group), "SV+", "SV-"),
    cn_group = ifelse(grepl("cnv\\+", mut_gain_group), "CNV+", "CNV-")
  ) %>%
  arrange(smk_group,,cn_group,sv_group,mut_group)


# plot
ggplot(summary_data, aes(x = Smoking_SBS4, y = percentage, 
                         fill = smk_group, 
                         alpha = cn_group,
                         pattern = mut_group),color="red") +
  geom_bar_pattern(
    stat = "identity", 
    position = "stack", 
    color = NA,
    pattern_fill = "black", 
    pattern_angle = 45, 
    pattern_density = 0.1,
    pattern_spacing = 0.1,
    pattern_key_scale_factor = 0.1
  ) +
  facet_grid(.~Smoking_SBS4, scales = "free_x") +
  ylab("Percentage of samples") +
  xlab("") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_fill_manual(values = c(
    "Non-Smoker" = "#00a800",
    "Smoker" = "#FF6EC7"
  )) +
  scale_alpha_manual(values = c(
    "CNV+" = 1,
    "CNV-" = 0.2
  )) +
  scale_pattern_manual(values = c(
    "Mutation+" = "stripe",  # Stripe fill for mutation+
    "Mutation-" = "none",     # No pattern for mutation-
    "SV+" = "stripe",  # Stripe fill for mutation+
    "SV-" = "none"     # No pattern for mutation-
  )) + mytheme
ggsave(file.path(outpath,"EGFR_gain_by_EGFR_SV.pdf"),width = 5,height = 5)
ggsave(file.path(outpath,"EGFR_gain_by_EGFR_SNV.pdf"),width = 5,height = 5)


out = summary_data %>% 
  group_by(Smoking_SBS4,cn_group,sv_group) %>%
  mutate(freq = sum(count)) %>%
  group_by(Smoking_SBS4,cn_group) %>%
  mutate(pct=freq / sum(count) * 100) %>% 
  dplyr::select(Smoking_SBS4,cn_group,sv_group,pct) %>% unique()

summary_data %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")) 


########### oncoprint #############
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*1, gp = gpar(fill = "#f7f9f9", col = NA))  # Background (light gray)
  },
  L858R = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#1B9E77", col = NA))
  },
  Exon19del = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#D95F02", col = NA))
  },
  Others = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#7570B3", col = NA))
  },
  WT = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "white", col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#33a02c", col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#1f78b4", col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#6a3d9a", col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#b15928", col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#df65b0", col = NA))
  },
  
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#a50f15", col = NA))
  },
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#8dd3c7", col = NA))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "green", col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#ff7f00", col = NA))
  },
  none = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "white", col = NA))
  },
  `cnv+` = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#e41a1c", col = NA))
  },
  `cnv-` = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#253494", col = NA))
  },
  `Non-Smoker-noSBS4` = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#00a800", col = NA))
  },
  `Smoker-SBS4` = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*1, gp = gpar(fill = "#FF6EC7", col = NA))
  }
)


# Prepare the annotations
plot_data <- plot_data %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")) 
plot_data <- plot_data %>% mutate(Annotation=ifelse(Tumor_Barcode %in% egfr_mutated_E479_A483del_samples,"Exon19del",
                                            ifelse(Tumor_Barcode %in% egfr_mutated_L591R_samples,"L858R",
                                                   ifelse(Tumor_Barcode %in% egfr_mutated_others_samples,"Others",
                                                          "WT")))) %>% 
  mutate(Annotation=factor(Annotation,levels = c("L858R","Exon19del","Others","WT")))

oncoprint_matrix <- plot_data %>% dplyr::select(Tumor_Barcode,Smoking_SBS4,cnv,EGFR,Annotation) %>% 
  mutate(EGFR=ifelse(is.na(EGFR),"",EGFR)) %>%
  as.matrix()
oncoprint_matrix[oncoprint_matrix == ""] <- "none"
rownames(oncoprint_matrix) <- oncoprint_matrix[,"Tumor_Barcode"]
# Create an ordering dataframe
sample_order_df <- data.frame(
  Sample = oncoprint_matrix[, "Tumor_Barcode"],
  Smoking_SBS4 = oncoprint_matrix[, "Smoking_SBS4"],
  cnv = oncoprint_matrix[, "cnv"],
  Annotation = oncoprint_matrix[,"Annotation"],
  EGFR = oncoprint_matrix[, "EGFR"],
  stringsAsFactors = FALSE
) %>%
  mutate(
    Smoking_SBS4 = factor(Smoking_SBS4, levels = c("Non-Smoker-noSBS4", "Smoker-SBS4")),
    cnv = factor(cnv, levels = c("cnv+", "cnv-")),
    Annotation=factor(Annotation,levels = c("L858R","Exon19del","Others","WT"))) %>%
  arrange(Smoking_SBS4, cnv,Annotation, EGFR)


# Now get the sample order
sample_order <- sample_order_df$Sample

pdf(file.path(outpath,"oncoprint.smk.cnv.egfr.pdf"),height = 7,width = 4)
oncoPrint(
  oncoprint_matrix[sample_order, c("Smoking_SBS4", "cnv", "Annotation")],
  alter_fun = alter_fun,
  col = c(
    "L858R"="#1B9E77",
    "Exon19del"="#D95F02",
    "Others"="#7570B3",
    "WT"="white",
    "Frame_Shift_Del" = "#1f78b4",
    "Frame_Shift_Ins" = "#6a3d9a",
    "In_Frame_Del" = "#b15928",
    "In_Frame_Ins" = "#df65b0",
    "Missense_Mutation" = "#33a02c",
    "Nonsense_Mutation" = "#a50f15",
    "Translation_Start_Site" = "#8dd3c7",
    "Nonstop_Mutation" = "green",
    "Splice_Site" = "#ff7f00",
    "cnv+" = "#e41a1c",
    "cnv-" = "#253494",
    "Non-Smoker-noSBS4" = "#00a800",
    "Smoker-SBS4" = "#FF6EC7",
    "none" = "white"
  ),
  row_order = sample_order,
  show_row_names = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6)
)
dev.off()

############################# fisher's test co-occurence and mutually exclusive ###################
# Step 1: Prepare your data
plot_data2 <- plot_data %>%
  mutate(
    EGFRmut = ifelse(Annotation %in% c("L858R", "Exon19del", "Others"), "Mut", "WT"),
    CNVstatus = case_when(
      cnv == "cnv+" ~ "CNV+",
      cnv == "cnv-" ~ "CNV-",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CNVstatus))  # remove samples without CNV info

# Fisher's test at each Smoking_SBS4 group
fisher_results <- plot_data2 %>%
  group_by(Smoking_SBS4) %>%
  group_split() %>%
  map_dfr(~ {
    temp <- .x
    tbl <- table(temp$EGFRmut, temp$CNVstatus)
    print(tbl)  # Print the table to check
    
    if (all(dim(tbl) == c(2,2))) {  # Only if table is 2x2
      # Correct order: (row = Mut or WT, column = CNV+ or CNV-)
      fisher_test <- fisher.test(matrix(c(
        tbl["WT", "CNV-"],
        tbl["WT", "CNV+"],
        tbl["Mut", "CNV-"],
        tbl["Mut", "CNV+"]
        ),
        nrow = 2, byrow = TRUE  # must specify byrow = TRUE
      ))
      print(fisher_test)  # Print the Fisher test result
      
      tibble(
        Smoking_SBS4 = unique(temp$Smoking_SBS4),
        p_value = fisher_test$p.value,
        odds_ratio = fisher_test$estimate,
        log2.odds_ratio = log2(fisher_test$estimate),
        CNVplus_EGFRmut = tbl["Mut", "CNV+"],
        CNVplus_EGFRwt = tbl["WT", "CNV+"],
        CNVminus_EGFRmut = tbl["Mut", "CNV-"],
        CNVminus_EGFRwt = tbl["WT", "CNV-"]
      )
    } else {
      tibble(
        Smoking_SBS4 = unique(temp$Smoking_SBS4),
        p_value = NA,
        odds_ratio = NA,
        CNVplus_EGFRmut = NA,
        CNVplus_EGFRwt = NA,
        CNVminus_EGFRmut = NA,
        CNVminus_EGFRwt = NA
      )
    }
  })


############################### Breslow-Day test ######################
# Smoking_SBS4 == "Non-Smoker-noSBS4"
non_smoker_tbl <- matrix(c(
  plot_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4") %>% filter(cnv=="cnv+"&mut=="mut+") %>% nrow(),
  plot_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4") %>% filter(cnv=="cnv+"&mut=="mut-") %>% nrow(),
  plot_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4") %>% filter(cnv=="cnv-"&mut=="mut+") %>% nrow(),
  plot_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4") %>% filter(cnv=="cnv-"&mut=="mut-") %>% nrow()
), nrow = 2, byrow = TRUE)

# Smoking_SBS4 == "Smoker-SBS4"
smoker_tbl <- matrix(c(
  plot_data %>% filter(Smoking_SBS4 == "Smoker-SBS4") %>% filter(cnv=="cnv+"&mut=="mut+") %>% nrow(),
  plot_data %>% filter(Smoking_SBS4 == "Smoker-SBS4") %>% filter(cnv=="cnv+"&mut=="mut-") %>% nrow(),
  plot_data %>% filter(Smoking_SBS4 == "Smoker-SBS4") %>% filter(cnv=="cnv-"&mut=="mut+") %>% nrow(),
  plot_data %>% filter(Smoking_SBS4 == "Smoker-SBS4") %>% filter(cnv=="cnv-"&mut=="mut-") %>% nrow()
), nrow = 2, byrow = TRUE)

# Combine into an array
tbl_array <- array(c(non_smoker_tbl, smoker_tbl), 
                   dim = c(2, 2, 2),
                   dimnames = list(
                     CNV = c("CNV+", "CNV-"),
                     EGFRmut = c("Mut", "WT"),
                     SmokingStatus = c("Non-Smoker", "Smoker")
                   ))


BreslowDayTest(tbl_array)
mantelhaen.test(tbl_array)


