### starfish generate cmplx and unclustered  sv

### install starfish
Packages <-
  c(
    "ShatterSeeky",
    "GenomeInfoDb",
    "plyr",
    "data.table",
    "GenomicRanges",
    "IRanges",
    "MASS",
    "ggplot2",
    "grid",
    "gridExtra",
    "dplyr",
    "ConsensusClusterPlus",
    "factoextra",
    "gplots",
    "ggpubr",
    "reshape2",
    "cowplot",
    "scales",
    "patchwork",
    "Cairo",
    "ggforce"
  )
lapply(Packages, library, character.only = TRUE)



###test starfish
library(Starfish)
library(data.table)
library(dplyr)
library(neuralnet)

library(biomaRt)


### INPUT
scratch.path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch"
all_sv_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
cnv_path <- file.path("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/Data/CNV","cnv_starfish.tsv")
sampleinfo.path <- file.path(scratch.path,"01-SV-files/1217_samples_hg38","sherlock_lung_final_wgs_1217_groups.txt")
cgr.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38","sherlock_connected_CGR_event.tsv")

outpath <- file.path(scratch.path,"12-individual-plot")
if (dir.exists(outpath)==FALSE){
  dir.create(outpath,showWarnings = F,recursive = T)
}
setwd(outpath)

### read files
all_sv <- read.delim(all_sv_path) %>% 
  dplyr::select(CHROM1,POS1,CHROM2,POS2,STRAND1,STRAND2,SVTYPE,SAMPLE) %>% 
  data.table::setnames(.,old=c("CHROM1","POS1","CHROM2","POS2","STRAND1","STRAND2","SVTYPE","SAMPLE"),
                       new=c("chrom1","pos1","chrom2","pos2","strand1","strand2","svtype","sample"))
cnv <- read.delim(cnv_path,header = F) %>% 
  setnames(.,old=c("V1","V2","V3","V4","V5"),
           new=c("chromosome","start","end","total_cn","sample"))
cgr <- read.delim(cgr.path)




### run modified starfish plot individual 
# source("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/OpenPBTA/analyses/03-starfish/04-starfish-plot-modified-individual.R")
# starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
#                                   mysample="BS_6T87HSPX",mycluster="BS_6T87HSPX_10_11",mychr=c("chr18","chr10","chr11"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
#                                   pdfname,genome_v = "hg38")



### run modified starfish plot individual  zoom in
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") #hg38
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="asia.ensembl.org")

# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org") #hg19



regions <- c("7:55000000:56000000", "14:37000000:38000000","14:36000000:37000000")
genes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
               filters = "chromosomal_region",
               values = regions,
               mart = mart) %>% as.data.frame() %>% pull(hgnc_symbol)
genes <- genes[genes!=""]
# gene_info <- getBM(mart, attributes=c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "ensembl_transcript_id", "strand","hgnc_symbol"), 
#                    filters = "hgnc_symbol", values = genes, uniqueRows=TRUE)
# exon_info <- getBM(mart, attributes=c("chromosome_name","exon_chrom_start", "exon_chrom_end", "ensembl_transcript_id", "ensembl_exon_id",
#                                       # "transcript_count",
#                                       "strand"), 
#                    filters = "ensembl_transcript_id", values = c("ENST00000275493","ENST00000250448"), uniqueRows=TRUE)
gene_info <- getBM(attributes = c("hgnc_symbol", 
                                  "ensembl_transcript_id", 
                                  "chromosome_name",
                                  "transcript_start", 
                                  "transcript_end", 
                                  "transcript_biotype"), 
                   filters = "hgnc_symbol", 
                   values = genes, 
                   mart = mart)
classic_transcripts <- gene_info[gene_info$transcript_biotype == "protein_coding", ]
classic_transcripts$length <- abs(classic_transcripts$transcript_start-classic_transcripts$transcript_end)
classic_transcripts <- classic_transcripts[order(classic_transcripts$hgnc_symbol, 
                                                 classic_transcripts$length, 
                                                 decreasing = TRUE), ]
classic_transcripts <- classic_transcripts[!duplicated(classic_transcripts$hgnc_symbol), ]




# # add exon id
# for (i in 1:nrow(exon_info)){
#   ti <- exon_info$ensembl_transcript_id[i]
#   starti <- exon_info$exon_chrom_start[i]
#   strandi <- exon_info$strand[i]
#   if(strandi=="-1"){
#     start.sort <- sort(exon_info[exon_info$ensembl_transcript_id==ti,"exon_chrom_start"],decreasing = T)
#   }else{
#     start.sort <- sort(exon_info[exon_info$ensembl_transcript_id==ti,"exon_chrom_start"])
#   }
#   
#   rank <- which(start.sort == starti)
#   exon_info$order[i] <- rank
# }
# plot
source("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/analysis/05-shatterseek/starfish-plot-modified-individual-zoomin.R")
# EGFR
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0358-T01",
  mycluster = "NSLC-0358-T01_7_10",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0358-T01_ecDNA_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "IGC-02-1105-T01",
  mycluster = "IGC-02-1105-T01_7_14_18",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/IGC-02-1105-T01_BFB_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0054-T01",
  mycluster = "NSLC-0054-T01_7_9_12_16_19",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0054-T01_BFB_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0067-T01",
  mycluster = "NSLC-0067-T01_7_12",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0067-T01_BFB_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0295-T01",
  mycluster = "NSLC-0295-T01_7_9_13",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0295-T01_BFB_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0336-T01",
  mycluster = "NSLC-0336-T01_1_5_7_9_11_20_22",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0336-T01_BFB_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0017-T01",
  mycluster = "NSLC-0017-T01_7",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0017-T01_largegain_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0678-T01",
  mycluster = NA,
  mychr = c("chr7"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0678-T01_cti_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "LU-FF1_tumor",
  mycluster = NA,
  mychr = c("chr7"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/LU-FF1_tumor_cplxunclear_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0685-T01",
  mycluster = NA,
  mychr = c("chr1","chr5","chr7","chr15","chr22"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0685-T01_cplxunclear_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0201-T01",
  mycluster = NA,
  mychr = c("chr11","chr12","chr14","chr20","chr7","chrX"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0201-T01_cplxunclear_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0607-T01",
  mycluster = NA,
  mychr = c("chr7"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0607-T01_fbinv_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0859-T01",
  mycluster = NA,
  mychr = c("chr7"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr7", 55000000,56000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 7),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 7) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0859-T01_fbinv_EGFR_zoomin.pdf"),
  genome_v = "hg38"
)

# FOXA1
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0225-T01",
  mycluster = "NSLC-0225-T01_3_5_14",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0225-T01_ecDNA_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "E00934102",
  mycluster = "E00934102_3_5_12_14_17_21",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/E00934102_ecDNA_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "E00934102",
  mycluster = "E00934102_3_5_12_14_17_21",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/E00934102_ecDNA_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)

starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-1049-T01",
  mycluster = "NSLC-1049-T01_12_14_16",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-1049-T01_bfb_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0539-T01",
  mycluster = "NSLC-0539-T01_14_15",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0539-T01_largegain_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-1186-T01",
  mycluster = "NSLC-1186-T01_2_6_14_X",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-1186-T01_largegain_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0011-T01",
  mycluster = "NSLC-0011-T01_14_19_20",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0011-T01_largegain_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)

starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "LU-FF62_tumor",
  mycluster = NA,
  mychr = "chr14",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 37000000,38000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/LU-FF62_tumor_fbinv_FOXA1_zoomin.pdf"),
  genome_v = "hg38"
)

# NKX2-1
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0225-T01",
  mycluster = "NSLC-0225-T01_3_5_14",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 36000000,37000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0225-T01_ecDNA_NKX2-1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0088-T01",
  mycluster = "NSLC-0088-T01_2_6_14",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 36000000,37000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0088-T01_BFB_NKX2-1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0192-T01",
  mycluster = "NSLC-0192-T01_12_14",
  mychr = NA,# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 36000000,37000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0192-T01_largegain_NKX2-1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0750-T01",
  mycluster = NA,
  mychr = paste0("chr","14"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 36000000,37000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0750-T01_cti_NKX2-1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "EBUS030713B",
  mycluster = NA,
  mychr = paste0("chr","14"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 36000000,37000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/EBUS030713B_TD2_NKX2-1_zoomin.pdf"),
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0121-T01",
  mycluster = NA,
  mychr = paste0("chr","14"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr14", 36000000,37000000),
  transcript_df = classic_transcripts %>% filter(chromosome_name == 14),
  mytranscriptid=classic_transcripts %>% filter(chromosome_name == 14) %>% pull(ensembl_transcript_id),
  pdfname = file.path(outpath,"plot/NSLC-0121-T01_fbinv_NKX2-1_zoomin.pdf"),
  genome_v = "hg38"
)

