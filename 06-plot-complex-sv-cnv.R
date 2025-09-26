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

library("Starfish")
library(data.table)
library(dplyr)
library(neuralnet)
library(biomaRt)
source("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/OpenPBTA/analyses/03-starfish/04-starfish-plot-modified.R")
source("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/OpenPBTA/analyses/03-starfish/04-starfish-plot-modified-individual.R")
source("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/OpenPBTA/analyses/03-starfish/04-starfish-plot-modified-individual-zoomin.R")

### path
scratch.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/scratch"
all_sv_path <- file.path(scratch.path,"01-SV-files","1217_samples_hg38","1217_manta_meerkat_union_window50_highquality.txt")
all_sv_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
cnv_path <- file.path("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/Data/CNV","cnv_starfish.tsv")
sampleinfo.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/Data/SAMPLE_INFO/wgs_1217_info.20250128.txt"
cgr.sig.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_pcawg_6signatures_class.tsv")
cgr.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_connected_CGR_event.tsv")

### read files
samplebsid <- read.delim(sampleinfo.path,header = TRUE) %>% pull(Tumor_Barcode)
all_sv <- read.delim(all_sv_path) %>% 
  dplyr::select(CHROM1,POS1,CHROM2,POS2,STRAND1,STRAND2,SVTYPE,SAMPLE,Signature) %>% 
  data.table::setnames(.,old=c("CHROM1","POS1","CHROM2","POS2","STRAND1","STRAND2","SVTYPE","SAMPLE"),
                       new=c("chrom1","pos1","chrom2","pos2","strand1","strand2","svtype","sample"))
cnv <- read.delim(cnv_path,header = F) %>% dplyr::select(V1,V2,V3,V4,V5) %>% 
  setnames(.,old=c("V1","V2","V3","V4","V5"),
           new=c("chromosome","start","end","total_cn","sample")) %>% filter(sample %in% samplebsid)
cgr.sig <-  read.delim(cgr.sig.path)
cgr <- read.delim(cgr.path)
cgr <- merge(cgr,cgr.sig[,c("cluster_id","CGR_signature")],all.x=T)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host="www.ensembl.org") #hg38
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org") #hg19

#######################################
### run modified starfish plot ########
#######################################
# 1. ecDNA
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","1.ecDNA")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified(sv_file=all_sv,cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="1 ecDNA/double minutes"),genome_v = "hg38")
starfish_plot_modified(sv_file=all_sv %>% filter(sample=="NSLC-0012-T01"),cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="1 ecDNA/double minutes"),genome_v = "hg38")
# 2. BFB
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","2.BFB")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified(sv_file=all_sv,cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="2 BFB cycles/chromatin bridge"),genome_v = "hg38")
starfish_plot_modified(sv_file=all_sv %>% filter(sample=="FH5PJAZF"),cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="2 BFB cycles/chromatin bridge"),genome_v = "hg38")
# 3. Large loss
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","3.Largeloss")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified(sv_file=all_sv,cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="3 Large loss"),genome_v = "hg38")
starfish_plot_modified(sv_file=all_sv %>% filter(sample=="NSLC-0066-T01"),cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="3 Large loss"),genome_v = "hg38")
# 4. Micronuclei
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","4.Micronuclei")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified(sv_file=all_sv,cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="4 Micronuclei"),genome_v = "hg38")
starfish_plot_modified(sv_file=all_sv %>% filter(sample=="LU-FF21_tumor"),cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="4 Micronuclei"),genome_v = "hg38")
# 5. Large gain
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","5.Largegain")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified(sv_file=all_sv,cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="5 Large gain"),genome_v = "hg38")
starfish_plot_modified(sv_file=all_sv %>% filter(sample=="NSLC-0823-T01"),cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="5 Large gain"),genome_v = "hg38")
# 6. Hourglass
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","6.Hourglass")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified(sv_file=all_sv,cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="6 Hourglass"),genome_v = "hg38")
starfish_plot_modified(sv_file=all_sv %>% filter(sample=="NSLC-0174-T01"),cnv_file=cnv,cgr=cgr %>% filter(CGR_signature=="6 Hourglass"),genome_v = "hg38")

# 7. chromoplexy
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","chromoplexy")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0098-T02",mycluster=NA,
                                  mychr=c(all_sv %>% filter(sample=="NSLC-0098-T02" &Signature=="chromoplexy") %>% pull(chrom1) %>% as.vector(),
                                          all_sv %>% filter(sample=="NSLC-0098-T02" &Signature=="chromoplexy") %>% pull(chrom2) %>% as.vector())%>% unique(),
                                  # mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0098-T02.chromoplexy.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv %>% filter(Signature=="chromoplexy"),cnv_file=cnv,cgr,
                                  mysample="NSLC-0098-T02",mycluster=NA,
                                  mychr=c(all_sv %>% filter(sample=="NSLC-0098-T02" &Signature=="chromoplexy") %>% pull(chrom1) %>% as.vector(),
                                          all_sv %>% filter(sample=="NSLC-0098-T02" &Signature=="chromoplexy") %>% pull(chrom2) %>% as.vector())%>% unique(),
                                  # mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0098-T02.chromoplexy.only.pdf",genome_v = "hg38")

# 8. Complex unclear
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","Complex_unclear")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0035-T01",mycluster=NA,
                                  mychr=c(all_sv %>% filter(sample=="NSLC-0035-T01" &Signature=="complex_unclear") %>% pull(chrom1) %>% as.vector(),
                                          all_sv %>% filter(sample=="NSLC-0035-T01" &Signature=="complex_unclear") %>% pull(chrom2) %>% as.vector())%>% unique(),
                                  # mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0035-T01.complex_unclear.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv %>% filter(Signature=="complex_unclear"),cnv_file=cnv,cgr,
                                  mysample="NSLC-0035-T01",mycluster=NA,
                                  mychr=c(all_sv %>% filter(sample=="NSLC-0035-T01" &Signature=="complex_unclear") %>% pull(chrom1) %>% as.vector(),
                                          all_sv %>% filter(sample=="NSLC-0035-T01" &Signature=="complex_unclear") %>% pull(chrom2) %>% as.vector())%>% unique(),
                                  # mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0035-T01.complex_unclear.only.pdf",genome_v = "hg38")
# 9. Complex unclear
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","cti")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0684-T01",mycluster=NA,
                                  mychr=c(all_sv %>% filter(sample=="NSLC-0684-T01" &Signature=="cycle_templated_ins") %>% pull(chrom1) %>% as.vector(),
                                          all_sv %>% filter(sample=="NSLC-0684-T01" &Signature=="cycle_templated_ins") %>% pull(chrom2) %>% as.vector())%>% unique(),
                                  # mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0684-T01.cti.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv %>% filter(Signature=="cycle_templated_ins"),cnv_file=cnv,cgr,
                                  mysample="NSLC-0684-T01",mycluster=NA,
                                  mychr=c(all_sv %>% filter(sample=="NSLC-0684-T01" &Signature=="cycle_templated_ins") %>% pull(chrom1) %>% as.vector(),
                                          all_sv %>% filter(sample=="NSLC-0684-T01" &Signature=="cycle_templated_ins") %>% pull(chrom2) %>% as.vector())%>% unique(),
                                  # mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0684-T01.cti.only.pdf",genome_v = "hg38")
###############################################
### run modified individual starfish plot #####
###############################################
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","individual")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="IGC-09-1042-T02",mycluster=NA,mychr=c("chr16","chr5","chr9","chr22"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="IGC-09-1042-T02.chr5.9.16.22.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0265-T01",mycluster=NA,mychr=c("chr1","chr22"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0265-T01.chr1.22.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0440-T01",mycluster=NA,mychr=c("chr20","chr22"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0440-T01.chr20.22.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0090-T01-R1",mycluster=NA,mychr=c("chr20","chr22"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0090-T01-R1.chr20.22.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0421-T01",mycluster=NA,mychr=c("chr20","chr22"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0421-T01.chr20.22.pdf",genome_v = "hg38")
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0461-T01",mycluster=NA,mychr=c("chr20","chr22"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0461-T01.chr20.22.pdf",genome_v = "hg38")

# simple SVs
# del
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-0417-T01",mycluster=NA,mychr=c("chr9"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-0417-T01.chr9.pdf",genome_v = "hg38")
# dup
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="S01356_20100601_TPD_01.01",mycluster=NA,mychr=c("chr15"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="S01356_20100601_TPD_01.01.chr15.pdf",genome_v = "hg38")
# inv
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="S01356_20100601_TPD_01.01",mycluster=NA,mychr=c("chr8"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="S01356_20100601_TPD_01.01.chr8.pdf",genome_v = "hg38")
# tra
starfish_plot_modified_individual(sv_file=all_sv,cnv_file=cnv,cgr,
                                  mysample="NSLC-1133-T01",mycluster=NA,mychr=c("chr4","chr11"),# mycluster has higher priority than mychr, mychr will use when mycluster is NA
                                  pdfname="NSLC-1133-T01.chr4.11.pdf",genome_v = "hg38")


######################################################
### run modified individual zoon in starfish plot ####
######################################################
genes <- c("TTC28","CHEK2")
enst <- c("ENST00000397906","ENST00000404276")
gene_info <- getBM(mart, attributes=c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "ensembl_transcript_id", "strand","hgnc_symbol"), 
                   filters = "hgnc_symbol", values = genes, uniqueRows=TRUE)
exon_info <- getBM(mart, attributes=c("chromosome_name","exon_chrom_start", "exon_chrom_end", "ensembl_transcript_id", "ensembl_exon_id",
                                      "transcript_count","strand"), 
                   filters = "ensembl_transcript_id", values = enst, uniqueRows=TRUE)
# add exon id
for (i in 1:nrow(exon_info)){
  ti <- exon_info$ensembl_transcript_id[i]
  starti <- exon_info$exon_chrom_start[i]
  strandi <- exon_info$strand[i]
  if(strandi=="-1"){
    start.sort <- sort(exon_info[exon_info$ensembl_transcript_id==ti,"exon_chrom_start"],decreasing = T)
  }else{
    start.sort <- sort(exon_info[exon_info$ensembl_transcript_id==ti,"exon_chrom_start"])
  }
  
  rank <- which(start.sort == starti)
  exon_info$order[i] <- rank
}

outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sv.cnv.plot","individual.zoomin")
if (dir.exists(outpath)==FALSE){dir.create(outpath,showWarnings = F,recursive = T)}
setwd(outpath)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "IGC-09-1042-T02",
  mycluster = NA,
  mychr = "chr22",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr22", 27863278,28885264),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst),
  mytranscriptid=enst,
  pdfname = "IGC-09-1042-T02.chr22.27863278.28885264.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0265-T01",
  mycluster = NA,
  mychr = "chr22",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr22", 27863278,28885264),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst),
  mytranscriptid=enst,
  pdfname = "NSLC-0265-T01.chr22.27863278.28885264.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "IGC-10-1185-T01",
  mycluster = NA,
  mychr = "chr22",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr22", 27863278,28885264),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst),
  mytranscriptid=enst,
  pdfname = "IGC-10-1185-T01.chr22.27863278.28885264.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "LU-51_tumor",
  mycluster = NA,
  mychr = "chr22",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr22", 27863278,28885264),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst),
  mytranscriptid=enst,
  pdfname = "LU-51_tumor.chr22.27863278.28885264.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0018-T01",
  mycluster = NA,
  mychr = "chr22",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr22", 27863278,28885264),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst),
  mytranscriptid=enst,
  pdfname = "NSLC-0018-T01.chr22.27863278.28885264.pdf",
  genome_v = "hg38"
)






genes <- c("MSL3","FRMPD4")
enst <- c("ENST00000361672","ENST00000656302")
gene_info <- getBM(mart, attributes=c("chromosome_name", "start_position", "end_position", "ensembl_gene_id", "ensembl_transcript_id", "strand","hgnc_symbol"), 
                   filters = "hgnc_symbol", values = genes, uniqueRows=TRUE)
exon_info <- getBM(mart, attributes=c("chromosome_name","exon_chrom_start", "exon_chrom_end", "ensembl_transcript_id", "ensembl_exon_id",
                                      "transcript_count","strand"), 
                   filters = "ensembl_transcript_id", values = enst, uniqueRows=TRUE)
# add exon id
for (i in 1:nrow(exon_info)){
  ti <- exon_info$ensembl_transcript_id[i]
  starti <- exon_info$exon_chrom_start[i]
  strandi <- exon_info$strand[i]
  if(strandi=="-1"){
    start.sort <- sort(exon_info[exon_info$ensembl_transcript_id==ti,"exon_chrom_start"],decreasing = T)
  }else{
    start.sort <- sort(exon_info[exon_info$ensembl_transcript_id==ti,"exon_chrom_start"])
  }
  
  rank <- which(start.sort == starti)
  exon_info$order[i] <- rank
}

starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0002-T01",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "NSLC-0002-T01.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0027-T01",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "NSLC-0027-T01.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0148-T01",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "NSLC-0148-T01.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0180-T01",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "NSLC-0180-T01.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)

starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0238-T01",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "NSLC-0238-T01.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "NSLC-0310-T01",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "NSLC-0310-T01.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "TCGA-49-4507-01A",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "TCGA-49-4507-01A.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)
starfish_plot_modified_individual_zoomin(
  sv_file = all_sv,
  cnv_file = cnv,
  cgr,
  mysample = "TCGA-MP-A5C7-01A",
  mycluster = NA,
  mychr = "chrX",# mycluster has higher priority than mychr, mychr will use when mycluster is NA
  zoomin_region = c("chr23", 11700000,12738198),
  transcript_df = exon_info %>% filter(ensembl_transcript_id %in% enst) %>% mutate(chromosome_name=23),
  mytranscriptid=enst,
  pdfname = "TCGA-MP-A5C7-01A.chrX.11700000.12738198.pdf",
  genome_v = "hg38"
)
