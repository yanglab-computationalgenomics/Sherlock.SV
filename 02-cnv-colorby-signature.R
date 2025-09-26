library(ggplot2)
library(dplyr)
library(data.table)
library(regioneR)
library(tidyr)
library(purrr)
library(GenomicRanges)
library(ggrepel)
library(ggpattern)
### function
plot_svsig_enrichment_bar <- function(
    plot_df,
    fisher_df=fisher_plot_df_capped,
    group_var = "Smoking",         # "Smoking" or "EGFR"
    feature_var = "SVtype",        # "SVtype" or "signature_closest"
    group_levels = c("Non-Smoker", "Smoker"),  # or c("WT", "Mut")
    siglevel = siglevel        # named color vector
) {
  # 1. Prepare percent table
  plot_pct <- plot_df %>%
    filter(.data[[group_var]] %in% group_levels) %>%
    group_by(gene, feature_value = .data[[feature_var]], group = .data[[group_var]]) %>%
    count() %>%
    group_by(gene, group) %>%
    mutate(
      percent = n / sum(n) * 100,
      percent_signed = percent
    ) %>%
    ungroup() %>%
    mutate(
      gene = factor(gene, levels = unique(plot_df$gene)),
      feature_value = factor(feature_value, levels = siglevel),
      group = factor(group, levels = group_levels)
    )
  
  # 2. Join fisher results and label enrichment
  plot_pct <- plot_pct %>%
    left_join(
      fisher_df %>%
        filter(group_by == group_var) %>%
        transmute(
          gene,
          feature_value = variable,
          odds_ratio = as.numeric(odds_ratio),
          FDR,
          enrichment_raw = case_when(
            FDR < 0.05 & odds_ratio > 1 ~ "enriched",
            FDR < 0.05 & odds_ratio < 1 ~ "depleted",
            TRUE ~ "none"
          )
        ),
      by = c("gene", "feature_value")
    ) %>%
    mutate(
      enrichment_type = case_when(
        enrichment_raw == "enriched" & group == group_levels[2] ~ "enriched",
        enrichment_raw == "depleted" & group == group_levels[1] ~ "enriched",
        TRUE ~ "none"
      )
    ) %>%
    mutate(
      gene = factor(gene, levels = c("CDKN2A", "PTPRD", "STK11",
                                     "TERT", "MDM2", "MYC", "NKX2-1",
                                     "PTK6", "EGFR", "CTNND2", "MCL1",
                                     "CCNE1", "BRAF", "MET", "FGFR1")),
      feature_value = factor(feature_value, levels = siglevel)) %>% 
    arrange(group,gene,feature_value)
  
  # 3. Label counts per gene and group
  gene_counts <- plot_pct %>%
    group_by(gene, group) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(
      y_pos = ifelse(group == group_levels[2], -2, 2),
      gene = factor(gene, levels = levels(plot_pct$gene)),
      group = factor(group, levels = group_levels)
    )
  
  # 4. Plot
  p <- ggplot(plot_pct, 
              aes(x = gene, y = percent_signed, fill = feature_value, color = enrichment_type,linewidth = enrichment_type)) +
    geom_bar(stat = "identity", width = 0.8,position = position_stack(reverse = FALSE)) +
    # geom_text(
    #   data = gene_counts,
    #   aes(x = gene, y = -10, label = n),
    #   inherit.aes = FALSE,
    #   size = 2,angle=90
    # ) +
    scale_fill_manual(values = sig_colors, drop = FALSE) +
    scale_color_manual(
      breaks = c("none", "enriched"),
      values = c("white", "black")
    ) +
    scale_linewidth_manual(
      breaks = c("none", "enriched"),
      values = c(0, 0.2)
    ) +
    scale_y_continuous(labels = function(x) paste0(x, "%"))+
    facet_wrap( .~group,nrow = 2,scales = "free") +
    labs(
      title = paste0(feature_var, " Supporting CNV Events by Gene"),
      x = "Gene",
      y = paste("Percentage (", group_levels[1], "up / ", group_levels[2], "down)", sep = ""),
      fill = feature_var
    ) +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
  p
  return(p)
}

run_fisher_by_gene <- function(df, feature_col, group_col) {
  feature_vals <- unique(df[[feature_col]])
  
  fisher_results <- map_dfr(feature_vals, function(val) {
    tbl <- table(df[[group_col]], df[[feature_col]] == val)
    
    if (all(dim(tbl) == c(2, 2))) {
      test <- fisher.test(tbl)
      
      tibble(
        gene = unique(df$gene),  # assumes each df is 1 gene
        variable = val,
        group_by = group_col,
        p_value = test$p.value,
        odds_ratio = test$estimate,
        n_total = nrow(df),
        table_11 = tbl[1,1],
        table_12 = tbl[1,2],
        table_21 = tbl[2,1],
        table_22 = tbl[2,2]
      )
    } else {
      tibble(
        gene = unique(df$gene),
        variable = val,
        group_by = group_col,
        p_value = NA_real_,
        odds_ratio = NA_real_,
        n_total = nrow(df),
        table_11 = NA_integer_,
        table_12 = NA_integer_,
        table_21 = NA_integer_,
        table_22 = NA_integer_
      )
    }
  }) %>%
    group_by(gene) %>%
    mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
    ungroup()
  
  return(fisher_results)
}


find_closest_sv_for_pos <- function(sample_id, chr, pos, sv_table) {
  sv_sample <- sv_table %>% filter(SAMPLE == sample_id)
  
  # Combine POS1 and POS2 breakpoints from CHROM1 and CHROM2
  sv_all <- bind_rows(
    sv_sample %>%
      transmute(
        SVID = as.character(SVID),
        Signature = as.character(Signature),
        SVTYPE = as.character(SVTYPE),
        pos_sv = as.integer(POS1),
        chrom = as.character(CHROM1)
      ),
    sv_sample %>%
      transmute(
        SVID = as.character(SVID),
        Signature = as.character(Signature),
        SVTYPE = as.character(SVTYPE),
        pos_sv = as.integer(POS2),
        chrom = as.character(CHROM2)
      )
  ) %>%
    filter(chrom == chr) %>%
    mutate(dist = abs(pos_sv - pos))
  
  if (nrow(sv_all) == 0) {
    return(tibble(
      svid = NA_character_,
      signature = NA_character_,
      dist = NA_integer_,
      rawsvtype = NA_character_
    ))
  }
  
  closest_idx <- which.min(sv_all$dist)
  
  return(tibble(
    svid = sv_all$SVID[closest_idx],
    signature = sv_all$Signature[closest_idx],
    dist = sv_all$dist[closest_idx],
    rawsvtype = sv_all$SVTYPE[closest_idx]
  ))
}

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
get_sample_gene_cn <- function(sample_name, cnv_df, gene_gr) {
  sample_cnv <- cnv_df %>% filter(sample == sample_name)
  
  if (nrow(sample_cnv) == 0) {
    return(tibble(sample = sample_name, gene = gene_positions$hgnc_symbol,
                  cn = NA_integer_, cn_start = NA_integer_, cn_end = NA_integer_))
  }
  
  cnv_gr <- GRanges(
    seqnames = sample_cnv$chromosome,
    ranges = IRanges(start = sample_cnv$start, end = sample_cnv$end),
    total_cn = sample_cnv$total_cn
  )
  
  overlaps <- findOverlaps(gene_gr, cnv_gr)
  
  overlap_tbl <- tibble(
    gene = as.character(mcols(gene_gr)$gene[queryHits(overlaps)]),
    cn = mcols(cnv_gr)$total_cn[subjectHits(overlaps)],
    cn_start = start(cnv_gr)[subjectHits(overlaps)],
    cn_end = end(cnv_gr)[subjectHits(overlaps)]
  ) %>% 
    group_by(gene) %>%
    arrange(cn) %>%
    dplyr::slice(1) %>%  # take the row with minimum CN
    ungroup()
  
  result <- tibble(gene = gene_positions$hgnc_symbol,
                   chr = gene_positions$chromosome_name) %>%
    left_join(overlap_tbl, by = "gene") %>%
    mutate(sample = sample_name) %>%
    dplyr::select(sample, gene, cn,chr, cn_start, cn_end)
  
  return(result)
}

### colors
sig_colors <- c(
  "1 ecDNA/double minutes" = "#e78ac3",
  "2 BFB cycles/chromatin bridge" = "#fc8d62",
  "3 Large loss" = "#a6d854",
  "4 Micronuclei" = "#66c2a5",
  "5 Large gain" = "#ffd92f",
  "6 Hourglass" = "#8da0cb",
  "chromoplexy" = "#5DADE2",
  "cycle_templated_ins" = "#F8C471",
  "complex_unclear" = "#C89EC7",
  "Del1" = "#D4E9F8",
  "Del2" = "#85C0E8",
  "Del3" = "#5EABDE",
  "TD1" = "#F8DAD7",
  "TD2" = "#F09389",
  "Fb inv" = "#3DBFC3",
  "Large intra" = "#EFC319",
  "Tra" = "#A06AAC",
  "ambiguous" = "grey90",
  "Unknown" = "grey90",
  "cluster cplx"="#F4924A",
  "noncluster cplx"="#91E83D",
  "simple"="#6E9EEA"
)


### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
cnv_path <- file.path(sherlock.path,"Data/CNV","cnv_starfish.tsv")
cnv_gistic_path <- file.path(sherlock.path,"Data/Significant_Focal_SCNAs/all_thresholded.by_genes.txt")
chrlength.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/Data/REF/hg38/genome_length/genome_length_hg38.tsv"
centromeres.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/Data/REF/hg38/centromeres/hg38.centromeres.bed"
table.path <- file.path(sherlock.path,"scratch","22-comprehensive-oncoplot","arm-level-and-gene-level-cn")
sv.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.main.txt")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
outpath <- file.path(scratch.path,"22-comprehensive-oncoplot","plot-cnv-color-by-signature")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

# choose my gene
gain.gene <- c("TERT","MDM2","EGFR","MYC","CTNND2","NKX2-1","PTK6","MCL1","FGFR1",
               # "GRIN2A",
               "CCNE1","MET","BRAF")
loss.gene <- c("CDKN2A","PTPRD","STK11")
mygene <- c(gain.gene,loss.gene)
loss_gistic <- paste0("gistic.",loss.gene)
gain_gistic <- paste0("gistic.",gain.gene)
loss_values <- paste0(loss.gene,".loss")
gain_values <- paste0(gain.gene,".gain")

### read
load("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/R/ensembl_hg38.RData")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
other_sample_infomation=read.delim(sampleinfo.path, header=T, check.names=F, stringsAsFactors=F, na.strings="no_na")
samplebsid <- read.delim(sampleinfo.path,header = TRUE) %>% pull(Tumor_Barcode)
cnv <- read.delim(cnv_path,header = F) %>% dplyr::select(V1,V2,V3,V4,V5) %>% 
  setnames(.,old=c("V1","V2","V3","V4","V5"),
           new=c("chromosome","start","end","total_cn","sample")) %>% filter(sample %in% samplebsid)
armcp <- read.delim2(file.path(table.path,"arm-level-major-cn.txt"),check.names = F)
chrcp <- read.delim2(file.path(table.path,"chr-level-major-cn.txt"),check.names = F)
samplecp <- read.delim2(file.path(table.path,"sample-level-major-cn.txt"),check.names = F)
cnv_gistic <- read.delim2(cnv_gistic_path,check.names = FALSE)
cnv_selected <- cnv_gistic %>% dplyr::select(-c(`Locus ID`, Cytoband))
cnv_transposed <- as.data.frame(t(cnv_selected[,-1]))
colnames(cnv_transposed) <- cnv_selected$`Gene Symbol`
cnv_transposed <- cbind(SampleID = rownames(cnv_transposed), cnv_transposed)
sv <- read.delim(sv.path)
sbs4 <- read.delim2(sbs4.path)
oncoprint_table <- read.delim(oncoprint_table_path,check.names = F) %>%
  mutate(across(all_of(gain_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==3,"gain",NA)))) %>%
  mutate(across(all_of(loss_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==-3,"loss",NA))))
oncoprint_table <- oncoprint_table %>% setnames(.,old=c(loss_gistic,gain_gistic),new=c(loss_values,gain_values))
oncoprint_table <- merge(oncoprint_table,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
oncoprint_table <- oncoprint_table %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
histology <- read.csv(histology.path)

### prepare
gene_positions <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position","hgnc_symbol"),
  filters = "hgnc_symbol",
  values = mygene,
  mart = ensembl
) %>% dplyr::select("chromosome_name", "start_position", "end_position","hgnc_symbol") %>%
  mutate(chromosome_name=as.character(chromosome_name))


###################### gene-level cn ###################
# gene GRanges
gene_gr <- GRanges(
  seqnames = gene_positions$chromosome_name,
  ranges = IRanges(chr = gene_positions$chromosome_name,start = gene_positions$start_position, end = gene_positions$end_position),
  gene = as.character(gene_positions$hgnc_symbol)
)

# gene cp and CNV position
genecp_with_pos <- unique(cnv$sample) %>%
  map_dfr(get_sample_gene_cn, cnv_df = cnv, gene_gr = gene_gr)

# merge all cns
# Step 1: gistic_long
gistic_long <- cnv_transposed %>%
  pivot_longer(-SampleID, names_to = "gene", values_to = "cnv_val") %>%
  setnames(.,old="SampleID",new="sample")

# Step 2: merge chromosome-level CN（chrcp is wide format）
chrcp_long <- chrcp %>%
  pivot_longer(-sample, names_to = "gene", values_to = "chr_cp")

# Step 3: merge arm-level CN（armcp is wide format）
armcp_long <- armcp %>%
  pivot_longer(-sample, names_to = "gene", values_to = "arm_cp")

# Step 4: merge sample-level modal CN
samplecp <- samplecp %>%
  dplyr::select(sample, chr_cn) %>% setnames(.,old="chr_cn",new="baseline_cp")

# Step 5: merge all
genecp_merged <- genecp_with_pos %>%
  left_join(gistic_long, by = c("sample", "gene")) %>%
  left_join(chrcp_long, by = c("sample", "gene")) %>%
  left_join(armcp_long, by = c("sample", "gene")) %>%
  left_join(samplecp, by = "sample")

# determin gain and loss
genecp_annotated <- genecp_merged %>%
  mutate(
    status = case_when(
      # Gain criteria
      gene %in% gain.gene &
        cnv_val %in% c(0, 1, 2) &
        !is.na(cn) & !is.na(arm_cp) & !is.na(chr_cp) & !is.na(baseline_cp) &
        cn > arm_cp & cn >= 3 & cn > chr_cp & cn > baseline_cp ~ "gain",
      
      # Loss criteria
      gene %in% loss.gene &
        cnv_val %in% c(0, -1, -2) &
        !is.na(cn) & !is.na(arm_cp) & !is.na(chr_cp) & !is.na(baseline_cp) &
        cn < 2 & cn < arm_cp & cn < chr_cp & cn < baseline_cp ~ "loss",
      
      # Otherwise
      TRUE ~ NA_character_
    )
  )

# add closest SV breakpoint
genecp_sv <- genecp_annotated %>% filter(status %in% c("gain", "loss")) 
genecp_sv_annotated <- genecp_sv %>%
  rowwise() %>%
  mutate(
    start = list(find_closest_sv_for_pos(sample, chr, cn_start, sv)),
    end   = list(find_closest_sv_for_pos(sample, chr, cn_end, sv))
  ) %>%
  ungroup() %>%
  unnest_wider(start, names_sep = "_") %>%
  unnest_wider(end, names_sep = "_") %>%
  mutate(
    svid_closest = ifelse(start_dist <= end_dist, start_svid, end_svid),
    signature_closest = ifelse(start_dist <= end_dist, start_signature, end_signature),
    rawsvtype_closest = ifelse(start_dist <= end_dist, start_rawsvtype, end_rawsvtype),
    dist_closest = pmin(start_dist, end_dist, na.rm = TRUE)
  )


# count  the percentage of start_signature == end_signature
genecp_sv_annotated %>%
  filter(!is.na(start_signature), !is.na(end_signature)) %>%
  mutate(signature_match = start_signature == end_signature) %>%
  summarise(
    total = n(),
    match = sum(signature_match),
    ratio = match / total
  )

# count # of pairs
signature_pairs <- genecp_sv_annotated %>%
  filter(!is.na(start_signature), !is.na(end_signature)) %>%
  mutate(
    sig_pair = map2_chr(start_signature, end_signature, ~ {
      sorted <- sort(c(.x, .y))
      paste(sorted, collapse = ":")
    })
  )
signature_pair_count <- signature_pairs %>%
  count(sig_pair, sort = TRUE)

### prepare plot table and test table
other_sample_infomation <- merge(other_sample_infomation,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")


plot_df <- merge(genecp_sv_annotated,
                 other_sample_infomation[,c("Tumor_Barcode","Smoking","Inclusion")],
                 by.x="sample",by.y="Tumor_Barcode",all.x=T) %>% filter(Inclusion=="In")
plot_df <- merge(plot_df,sbs4,by.x = "sample",by.y = "Tumor_Barcode",all.x = T)
plot_df <- plot_df %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4"))) # add SBS
plot_df <- plot_df %>% mutate(signature_closest = ifelse(is.na(signature_closest),"Unknown",signature_closest)) # add unknown
plot_df$EGFR <- ifelse(plot_df$sample %in% c(egfr_mutated_E479_A483del_samples,egfr_mutated_L591R_samples,egfr.others.path),"Mut","WT") # add egfr
plot_df <- plot_df %>% mutate(SVtype = ifelse(signature_closest %in% c("1 ecDNA/double minutes",
                                                                       "2 BFB cycles/chromatin bridge",
                                                                       "3 Large loss",
                                                                       "4 Micronuclei",
                                                                       "5 Large gain",
                                                                       "6 Hourglass"),"cluster cplx",
                                              ifelse(signature_closest %in% c("chromoplexy",
                                                                              "cycle_templated_ins",
                                                                              "complex_unclear"),"noncluster cplx",
                                                     ifelse(signature_closest %in% c("Del1",
                                                                                     "Del2",
                                                                                     "Del3",
                                                                                     "TD1",
                                                                                     "TD2",
                                                                                     "Fb inv",
                                                                                     "Large intra",
                                                                                     "Tra"),"simple","Unknown"))))
plot_df <- plot_df %>% mutate(SVtype=factor(SVtype,levels=c("cluster cplx","noncluster cplx","simple","Unknown")))
test_df <- plot_df

####################### test #########################
res_svtype_smoking <- test_df %>% 
  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")) %>% 
  filter(SVtype != "Unknown") %>%
  group_split(gene) %>% 
  map_dfr(~ run_fisher_by_gene(.x, feature_col = "SVtype", group_col = "Smoking"))
                                         
res_svtype_egfr    <- test_df %>% 
  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4")) %>%
  filter(SVtype != "Unknown") %>%
  group_split(gene) %>% 
  map_dfr(~ run_fisher_by_gene(.x, feature_col = "SVtype", group_col = "EGFR"))

res_sig_smoking    <- test_df %>% 
  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")) %>% 
  filter(SVtype != "Unknown") %>%
  group_split(gene) %>% 
  map_dfr(~ run_fisher_by_gene(.x, feature_col = "signature_closest", group_col = "Smoking"))

res_sig_egfr       <- test_df %>% 
  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4")) %>%
  filter(SVtype != "Unknown") %>%
  group_split(gene) %>% 
  map_dfr(~ run_fisher_by_gene(.x, feature_col = "signature_closest", group_col = "EGFR"))

# add labels
res_svtype_smoking$test <- "SVtype x Smoking"
res_svtype_egfr$test    <- "SVtype x EGFR"
res_sig_smoking$test    <- "Signature x Smoking"
res_sig_egfr$test       <- "Signature x EGFR"

# merge 
fisher_plot_df <- bind_rows(res_svtype_smoking, res_svtype_egfr, res_sig_smoking, res_sig_egfr) 

fisher_plot_df_capped <- fisher_plot_df %>%
  mutate(
    gene = factor(gene),
    raw_log2OR = log2(odds_ratio),
    log2OR = case_when(
      is.infinite(log2(odds_ratio)) & log2(odds_ratio) > 0 ~ 6,  # right cap for +Inf
      is.infinite(log2(odds_ratio)) & log2(odds_ratio) < 0 ~ -6, # left cap for -Inf
      TRUE ~ log2(odds_ratio)
    ),
    neglog10FDR = -log10(FDR + 1e-300),
    stroke_val = ifelse(FDR < 0.05, 0.5, 0),
    label_inf = case_when(
      is.infinite(raw_log2OR) ~ paste0(variable, " (Inf)"),
      TRUE ~ variable
    )
  )
# plot
ggplot(
  fisher_plot_df_capped,
  aes(x = log2OR, y = gene)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70",linewidth=0.1189) +
  geom_point(
    aes(fill = variable, size = neglog10FDR, stroke = stroke_val,
        color = ifelse(FDR < 0.05, "black", "white")),
    shape = 21
  ) +
  geom_text_repel(
    data = filter(fisher_plot_df_capped, FDR < 0.05),
    aes(label = label_inf),
    size = 1.5,
    max.overlaps = Inf,
    box.padding = 0.05,
    point.padding = 0.05
  ) +
  scale_color_manual(breaks = c("black","white"),values = c("black","white"))+
  scale_fill_manual(values = sig_colors, name = "SVtype / Signature") +
  scale_size_continuous(range = c(1, 2), name = "-log10(FDR)") +
  # set limits wider to leave space for capped ±Inf
  scale_x_continuous(limits = c(-6, 6)) +
  facet_wrap(~ test, scales = "free", ncol = 2) +
  labs(
    title = "Per-Gene Fisher Test: OR with ±Inf Capped",
    x = "log2(Odds Ratio)",
    y = "Gene"
  ) +
  mytheme +
  theme(axis.title.x = element_text(size=6,color="black"))
ggsave(file.path(outpath,"plot.fisher.svtype.and.svsig.with.smk.and.egfr.pdf"),width = 7,height = 5)




################### plot #################
# For SVtype x Smoking
plot_svsig_enrichment_bar(plot_df %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")), fisher_df=fisher_plot_df_capped, group_var = "Smoking", feature_var = "SVtype",siglevel = c("Unknown","cluster cplx","noncluster cplx","simple"))
ggsave(file.path(outpath,"pct.svtype.support.cnv.by.gene.consistentsmk.smoking.pdf"),width = 1.8,height = 3)

# For SVtype x EGFR
plot_svsig_enrichment_bar(plot_df %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4")), fisher_plot_df_capped, group_var = "EGFR", feature_var = "SVtype", group_levels = c("Mut","WT"),siglevel = c("Unknown","cluster cplx","noncluster cplx","simple"))
ggsave(file.path(outpath,"pct.svtype.support.cnv.by.gene.alltumor.egfr.pdf"),width = 1.8,height = 3)

# For SV signature x Smoking
plot_svsig_enrichment_bar(plot_df %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")), fisher_plot_df_capped, group_var = "Smoking", feature_var = "signature_closest",siglevel = c("Unknown","1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass","chromoplexy","cycle_templated_ins","complex_unclear","Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra"))
ggsave(file.path(outpath,"pct.svsig.support.cnv.by.gene.consistentsmk.smoking.pdf"),width = 1.8,height = 3) #3.464

# For SV signature x EGFR
plot_svsig_enrichment_bar(plot_df %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4")), fisher_plot_df_capped, group_var = "EGFR", feature_var = "signature_closest", group_levels = c("WT", "Mut"),siglevel = c("Unknown","1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass","chromoplexy","cycle_templated_ins","complex_unclear","Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra"))
ggsave(file.path(outpath,"pct.svsig.support.cnv.by.gene.alltumor.egfr.pdf"),width = 1.8,height = 3) #3.464

##############
out=plot_df %>%  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")) %>%
  filter(gene=="EGFR") %>% 
  group_by(Smoking_SBS4,SVtype,signature_closest,EGFR) %>% count() %>%
  mutate(signature_closest=factor(signature_closest,levels=names(sig_colors)))

ggplot(data=out) +
  geom_bar(aes(x=signature_closest,y=n,fill = EGFR),stat = "identity")+
  facet_grid(Smoking_SBS4~.)+
  scale_fill_manual(breaks = c("Mut","WT"),values = c("red","grey"))+
  ylab("Number of sample")+
  mytheme
ggsave(file.path(outpath,"egfr.gain.by.sig.egfrmut.annotation.pdf"),width = 5,height = 5)
