# modified from Xiaoming's script
# plot data (data_info)
# index   Sample_id       Status              group
# 1       2GUGKAEN        EUR                 population
# 310     LU-FF62_tumor   EAS                 population
# 8       IGC-02-1097-T0  Smoker_SBS4         smoking
# 311     NSLC-0466-T02   NonSmoker_NonSBS4   smoking


library(ggplot2)
library(grid)
library(gtable)
library(egg)
library(dplyr)
library(data.table)
library(tidyr)
library(readxl)
library(gridExtra)

### function
universal_theme = theme(axis.title.x=element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.line.x = element_blank(),
                        axis.text.y = element_text(size = 6),
                        axis.title.y = element_text(angle = 0, size = 7, vjust = 0.5),
                        axis.line.y = element_line(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        legend.text = element_text(size = 6),
                        legend.title = element_blank(),
                        legend.key.size = (unit(0.3, "cm")),
                        plot.margin = unit(c(0,0.5,0,0.5), "cm")
)

sort_by_count <- function(input_data_frame, select_column){
  select_column=enquo(select_column)
  temp_count <- input_data_frame %>% group_by(!!select_column) %>% count(sort=T)
  temp_count=as.data.frame(temp_count)
  colnames(temp_count)=c("select_column", "n")
  return(temp_count$select_column)
}


### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
fusion.detail.path  <- file.path(sherlock.path,"scratch/fusion/end_annotation/tier1_driver_fusions.txt")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
hq.path <- file.path(sherlock.path,"/Data/SAMPLE_INFO/HQ_samples.csv")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
kras.G12C.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12C_samples.txt")
kras.G12V.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12V_samples.txt")
kras.G12D.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12D_samples.txt")
kras.Q61H.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pQ61H_samples.txt")
kras.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_other_mutations.txt")
signature.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38","sherlock_CGR_signature.tsv")
sv_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_all_sv_with_signature.tsv")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","plot")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

### read
histology <- read.csv(histology.path)
fusion_set <- c("ALK","RET","ROS1","MET","NRG1",
                "FGFR1","FGFR2","FGFR3","FGFR4",
                "NTRK1","NTRK2","NTRK3")
all_sv  <- read.delim(sv_path,stringsAsFactors = FALSE)
fusion.detail <- read.delim2(fusion.detail.path)
fusion.detail <- fusion.detail %>%
  dplyr::mutate(
    Fusion_Gene = case_when(
      GENE1 %in% fusion_set ~ GENE1,
      GENE2 %in% fusion_set ~ GENE2,
      TRUE ~ NA_character_
    )
  ) %>% filter(Known_Fusion=="Yes")
other_sample_infomation=read.delim(sampleinfo.path, header=T, check.names=F, stringsAsFactors=F, na.strings="no_na")
tp53_mutated_samples=read.delim(tp53.mut.path, header=F, check.names=F, stringsAsFactors=F)
kras_mutated_samples=read.delim(kras.mut.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12c_samples=read.delim(kras.G12C.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12v_samples=read.delim(kras.G12V.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12d_samples=read.delim(kras.G12D.path, header=F, check.names=F, stringsAsFactors=F)
kras_q61h_samples=read.delim(kras.Q61H.path, header=F, check.names=F, stringsAsFactors=F)
kras_others_samples=read.delim(kras.others.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F)
cluster=read.delim(signature.path, header=T, sep="\t", stringsAsFactors = F)
colnames(cluster)=c("cgr_chr", "chr","start", "end", "cgr_status", "link_chr", "mechanism", "sample")
cluster <- cluster %>% 
  mutate(mechanism=ifelse(mechanism=="1 ecDNA/double minutes","1 ecDNA",cluster$mechanism)) %>%
  mutate(mechanism=ifelse(mechanism=="2 BFB cycles/chromatin bridge","2 BFB",mechanism)) %>% 
  mutate(mechanism=ifelse(mechanism=="","Unassigned",mechanism))
sbs4 <- read.delim2(sbs4.path)
hqsample <- read.csv(hq.path) %>% pull(Tumor_Barcode)

### add everything
other_sample_infomation <- merge(other_sample_infomation,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")
merge_sample_info <- other_sample_infomation
merge_sample_info <- merge_sample_info %>% mutate(fusion_status=ifelse(Tumor_Barcode %in% fusion.detail$SAMPLE,NA,"wildtype"))
merge_sample_info$Population <- merge_sample_info$Assigned_Population
merge_sample_info$Population <- ifelse(merge_sample_info$Population=="EUR","a.European",merge_sample_info$Population)
merge_sample_info$Population <- ifelse(merge_sample_info$Population=="EAS","b.Asian",merge_sample_info$Population)
merge_sample_info$Population <- ifelse(merge_sample_info$Population %in% c("AFR","AMR or Mixed"),"c.Others",merge_sample_info$Population)
merge_sample_info$Smoking <- ifelse(merge_sample_info$Smoking=="Non-Smoker","a.Never smoker",merge_sample_info$Smoking)
merge_sample_info$Smoking <- ifelse(merge_sample_info$Smoking=="Smoker","b.Smoker",merge_sample_info$Smoking)
merge_sample_info$Smoking <- ifelse(merge_sample_info$Smoking=="Unknown","c.Unknown",merge_sample_info$Smoking)
merge_sample_info$Histology <- ifelse(merge_sample_info$Histology %in% c("Adenocarcinoma"),"a.Adenocarcinoma",merge_sample_info$Histology)
merge_sample_info$Histology <- ifelse(merge_sample_info$Histology %in% c("Squamous cell carcinoma"),"b.Squamous cell carcinoma",merge_sample_info$Histology)
merge_sample_info$Histology <- ifelse(merge_sample_info$Histology %in% c("Carcinoid tumor"),"c.Carcinoid tumor",merge_sample_info$Histology)
merge_sample_info$Histology <- ifelse(merge_sample_info$Histology %in% c("Adenosquamous carcinoma"),"d.Adenosquamous carcinoma",merge_sample_info$Histology)
merge_sample_info$Histology <- ifelse(merge_sample_info$Histology %in% c("Others"),"e.Others",merge_sample_info$Histology)
merge_sample_info$TP53 <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(tp53_mutated_samples)),"Mut","WT")
merge_sample_info$KRAS <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(kras_g12c_samples)),"G12C","WT")
merge_sample_info$KRAS <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(kras_g12v_samples)),"G12V",merge_sample_info$KRAS)
merge_sample_info$KRAS <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(kras_g12d_samples)),"G12D",merge_sample_info$KRAS)
merge_sample_info$KRAS <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(kras_q61h_samples)),"Q61H",merge_sample_info$KRAS)
merge_sample_info$KRAS <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(kras_others_samples)),"Other_mut",merge_sample_info$KRAS)
merge_sample_info$EGFR <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(egfr_mutated_E479_A483del_samples)),"E479_A483del","WT")
merge_sample_info$EGFR <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(egfr_mutated_L591R_samples)),"L591R",merge_sample_info$EGFR)
merge_sample_info$EGFR <- ifelse(merge_sample_info$Tumor_Barcode %in% as.vector(unlist(egfr_mutated_others_samples)),"Other_mut",merge_sample_info$EGFR)
merge_sample_info <- merge(merge_sample_info,sbs4)

# add CGR signature
sample_cluster=unique(cluster[c("sample","mechanism")]) %>% data.table::setnames(.,old="mechanism",new="clustercplx")
sample_cluster_count=as.data.frame(sample_cluster %>% group_by(sample) %>% count())
sample_cluster <- sample_cluster %>% 
  mutate(clustercplx=ifelse(sample %in% sample_cluster_count[sample_cluster_count$n>1,"sample"],"Multiple",clustercplx)) %>% 
  unique()
merge_sample_info=merge(merge_sample_info, sample_cluster, by.x="Tumor_Barcode", by.y="sample", all.x=TRUE)
merge_sample_info <- merge_sample_info %>% mutate(clustercplx=ifelse(is.na(clustercplx),"no_CGR",clustercplx))
# add non-clustered cplx 
sample_noncluster <- all_sv %>% 
  dplyr::filter(Signature %in% c("chromoplexy","cycle_templated_ins","complex_unclear")) %>% dplyr::select(SAMPLE,Signature) %>% unique()
sample_noncluster_count <- all_sv %>% filter(Signature %in% c("chromoplexy","cycle_templated_ins","complex_unclear")) %>% 
  group_by(SAMPLE,Signature) %>% count() %>% group_by(SAMPLE) %>% count()
sample_noncluster <- sample_noncluster %>% 
  mutate(nonclustercplx=ifelse(SAMPLE %in% (sample_noncluster_count %>% filter(n>1) %>% pull(SAMPLE)),
                          "Multiple",Signature)) %>% 
  dplyr::select(SAMPLE,nonclustercplx) %>% unique()
merge_sample_info=merge(merge_sample_info, sample_noncluster, by.x="Tumor_Barcode", by.y="SAMPLE", all.x=TRUE)
merge_sample_info <- merge_sample_info %>% mutate(nonclustercplx=ifelse(is.na(nonclustercplx),"no_CGR",nonclustercplx))
### add TTC28 intro 1 transpose
ttc28_sample <- all_sv %>% filter((CHROM1==22 & POS1 >= 28600000 & POS1 <= 28700000 & Signature =="Tra")|(CHROM2==22 & POS2 >= 28600000 & POS2 <= 28700000  & Signature =="Tra")) %>% pull(SAMPLE)
merge_sample_info$TTC28 <- ifelse(merge_sample_info$Tumor_Barcode %in% ttc28_sample,"L1","No")
### add fusion detail
fusion_df <- data.frame(Sample = merge_sample_info$Tumor_Barcode, stringsAsFactors = FALSE)
fusion_df[fusion_set] <- NA
for (genei in fusion_set) {
  matched_samples <- fusion.detail %>%
    filter(Fusion_Gene == genei) %>%
    pull(SAMPLE)
  fusion_df[fusion_df$Sample %in% matched_samples, genei] <- "fusion"
}
fusion_df_long <- gather(fusion_df,fusion,gene,fusion_set) %>%
  filter(gene=="fusion") %>% dplyr::select(Sample,fusion)
merge_sample_info <- merge(merge_sample_info,fusion_df_long %>% setnames(.,old="fusion",new="fusion_detail"),
                           by.x = "Tumor_Barcode",by.y = "Sample",all.x = T)
merge_sample_info$fusion_detail <- ifelse(is.na(merge_sample_info$fusion_detail),"WT",merge_sample_info$fusion_detail)
### combine smoking and sbs4
merge_sample_info$Smoking_SBS4 <- paste0(merge_sample_info$Smoking,"-",ifelse(merge_sample_info$SBS4==0,"noSBS4","SBS4"))
# a.Never smoker-noSBS4   a.Never smoker-SBS4       b.Smoker-noSBS4         b.Smoker-SBS4      c.Unknown-noSBS4 
# 815                    56                    40                   305                     1
### filter 
#merge_sample_info <- merge_sample_info %>% filter(smoking != "Smoking_others")

#########################
## PLOT IN ALL SAMPLE ###
#########################
### order
merge_sample_info_all <- merge_sample_info %>% 
  #filter(Tumor_Barcode %in% hqsample) %>% # filter high quality
  mutate(Smoking=factor(Smoking,levels=c("b.Smoker","a.Never smoker","c.Unknown")),
         Population=factor(Population,levels=c("a.European","b.Asian","c.Others")),
         Smoking_SBS4=factor(Smoking_SBS4,levels=c("a.Never smoker-noSBS4","a.Never smoker-SBS4","b.Smoker-noSBS4","b.Smoker-SBS4","c.Unknown-noSBS4"))) %>%
  arrange(#Smoking_SBS4,
          Population,
          Smoking,
          clustercplx,
          nonclustercplx,
          TP53,EGFR,KRAS,fusion_status,Gender,Histology) 

merge_sample_info_all$index=seq(1,nrow(merge_sample_info_all))
saveRDS(merge_sample_info_all, file.path(outpath,"merge_sample_info_all.rds"))

### prepare plot info
data_info <- rbind(
  merge_sample_info_all[,c("index", "Tumor_Barcode","Population")] %>% setnames(.,old="Population",new="status") %>% mutate(group="Population"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","Smoking")] %>% setnames(.,old="Smoking",new="status") %>% mutate(group="Smoking"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","Smoking_SBS4")] %>% setnames(.,old="Smoking_SBS4",new="status") %>% mutate(group="Smoking_SBS4"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","EGFR")] %>% setnames(.,old="EGFR",new="status") %>% mutate(group="EGFR"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","TP53")] %>% setnames(.,old="TP53",new="status") %>% mutate(group="TP53"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","KRAS")] %>% setnames(.,old="KRAS",new="status") %>% mutate(group="KRAS"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","fusion_detail")] %>% setnames(.,old="fusion_detail",new="status") %>% mutate(group="Fusion"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","Gender")] %>% setnames(.,old="Gender",new="status") %>% mutate(group="Sex"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","Histology")] %>% setnames(.,old="Histology",new="status") %>% mutate(group="Tumor type"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","clustercplx")] %>% setnames(.,old="clustercplx",new="status") %>% mutate(group="Clustered complex SV signature"),
  merge_sample_info_all[,c("index", "Tumor_Barcode","nonclustercplx")] %>% setnames(.,old="nonclustercplx",new="status") %>% mutate(group="Nonclustered complex SV")
  #merge_sample_info_all[,c("index", "Tumor_Barcode","TTC28")] %>% setnames(.,old="TTC28",new="status") %>% mutate(group="TTC28")
)
data_info$group = factor(data_info$group, 
                         levels=c("Tumor type", 
                                  "Sex", 
                                  "Fusion", 
                                  "KRAS", 
                                  "EGFR",
                                  # "TTC28", 
                                  "TP53", 
                                  "Nonclustered complex SV" ,
                                  "Clustered complex SV signature", 
                                  "Smoking",
                                  "Population",
                                  #"Fusion",
                                  "Smoking_SBS4"
                                  ))
data_info$status = factor(data_info$status, levels=unique(data_info$status))

### plot
output_plot<-ggplot(data_info, aes(x=index, y=group)) + 
  geom_tile(aes(fill = factor(status),width=1,height=0.7)) + 
  scale_fill_manual("",values=c("a.Never smoker-noSBS4"="#00a800","a.Never smoker-SBS4"="#00a800","b.Smoker-noSBS4"="#FF6EC7","b.Smoker-SBS4"="#FF6EC7","c.Unknown-noSBS4"="#7C7C7C",
                                "a.European"="#FFD92F", "b.Asian"="#66C2A5","c.Others"="#7C7C7C",
                                "a.Never smoker"="#00a800","b.Smoker"="#FF6EC7","c.Unknown"="#7C7C7C",
                                "L591R"="#1B9E77","E479_A483del"="#D95F02","Other_mut"="#7570B3",
                                "G12C"="#DC0073","G12V"="#4A6FA5","G12D"="#23F0C7","Q61H"="#F9B3D1",
                                "Mut"="black", "WT"="white",
                                "fusion"="black", "wildtype"="white",
                                "ALK"="#ff7043","ROS1"="#cddc39","RET"="#26c6da","MET"="#80cbc4",
                                "FGFR1"="#6C5B7B","FGFR2"="#6C5B7B","FGFR3"="#6C5B7B","FGFR4"="#6C5B7B",
                                "NTRK1"="#CAD7B2","NTRK2"="#CAD7B2","NTRK3"="#CAD7B2",
                                "Male"="#8DA0CB","Female" = "#FC8D62",
                                "a.Adenocarcinoma" = "#4CAF50",  
                                "b.Squamous cell carcinoma" = "#2196F3",  
                                "c.Carcinoid tumor" = "#FF9800",  
                                "d.Adenosquamous carcinoma" = "#9C27B0",  
                                "e.Others" = "#F0F0F0",
                                "1 ecDNA"="#e78ac3",
                                "2 BFB"="#fc8d62", 
                                "3 Large loss"="#a6d854",
                                "4 Micronuclei"="#66c2a5",
                                "5 Large gain"="#ffd92f",
                                "6 Hourglass"="#8da0cb", 
                                "Multiple"="black", 
                                "no_CGR"="white",
                                "Unassigned"="white",
                                "chromoplexy"="#60ABDE",
                                "cycle_templated_ins"="#F7C472",
                                "complex_unclear"="#C99DC5"
                                # "L1"="black", "No"="white"
                                )) + 
  guides(fill=guide_legend(ncol=3)) +
  universal_theme + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = (unit(0.2, "cm")), legend.position="bottom"
  ) + 
  ylab("")

#pdf(file.path(outpath,"smoking_CGR_fusion_sample_onprint_draf.pdf"), useDingbats=FALSE, width=5, height=3)
#pdf(file.path(outpath,"smoking_CGR_fusion_sample_onprint_draf_hq.pdf"), useDingbats=FALSE, width=5, height=3)
pdf(file.path(outpath,"smoking_CGR_fusion_sample_onprint_draf_sbs4.pdf"), useDingbats=FALSE, width=5, height=5)
output_plot
dev.off()


####################################################################
data_info %>% count(group, status)

p.info.function <- function(myinfo,mylevel){
  merge_sample_info_all$Test=merge_sample_info_all[,myinfo]
  print(merge_sample_info_all)
  out=merge_sample_info_all %>% 
    filter(Smoking!="c.Unknown") %>%
    count(Smoking_SBS4,Test) %>% 
    group_by(Smoking_SBS4) %>% 
    mutate(sum=sum(n)) %>% 
    mutate(pct=n/sum) %>%
    mutate(Group=paste0(Smoking_SBS4,"(",sum,")")) %>%
    mutate(Test=as.character(Test))%>%
    mutate(Test = factor(Test, levels = mylevel)) %>%
    ggplot(data=.) + 
    geom_bar(aes(x=Group,y=pct,fill=Test),size=0.234,stat = "identity")+
    scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
    ggtitle(myinfo)+
    labs(y="Proportion of tumor")+
    scale_fill_manual("",values=c("a.Never smoker-noSBS4"="grey90","a.Never smoker-SBS4"="grey50","b.Smoker-noSBS4"="grey30","b.Smoker-SBS4"="grey10","c.Unknown-noSBS4"="white",
                                  "a.European"="#FFD92F", "b.Asian"="#66C2A5","c.Others"="#7C7C7C",
                                  "a.Never smoker"="#006400","b.Smoker"="#FF6EC7","c.Unknown"="#F0F0F0",
                                  "L591R"="#1B9E77","E479_A483del"="#D95F02","Other_mut"="#7570B3",
                                  "G12C"="#DC0073","G12V"="#4A6FA5","G12D"="#23F0C7","Q61H"="#F9B3D1",
                                  "Mut"="black", "WT"="#F0F0F0",
                                  "fusion"="black", "wildtype"="white",
                                  "ALK"="#ff7043","ROS1"="#cddc39","RET"="#26c6da","MET"="#80cbc4",
                                  "FGFR1"="#6C5B7B","FGFR2"="#6C5B7B","FGFR3"="#6C5B7B","FGFR4"="#6C5B7B",
                                  "NTRK1"="#CAD7B2","NTRK2"="#CAD7B2","NTRK3"="#CAD7B2",
                                  "Male"="#8DA0CB","Female" = "#FC8D62",
                                  "a.Adenocarcinoma" = "#4CAF50",  
                                  "b.Squamous cell carcinoma" = "#2196F3",  
                                  "c.Carcinoid tumor" = "#FF9800",  
                                  "d.Adenosquamous carcinoma" = "#9C27B0",  
                                  "e.Others" = "#F0F0F0",
                                  "1 ecDNA"="#e78ac3",
                                  "2 BFB"="#fc8d62", 
                                  "3 Large loss"="#a6d854",
                                  "4 Micronuclei"="#66c2a5",
                                  "5 Large gain"="#ffd92f",
                                  "6 Hourglass"="#8da0cb", 
                                  "Multiple"="black", 
                                  "no_CGR"="#F0F0F0",
                                  "Unassigned"="#F0F0F0",
                                  "chromoplexy"="#60ABDE",
                                  "cycle_templated_ins"="#F7C472",
                                  "complex_unclear"="#C99DC5"
    ),guide = guide_legend(nrow = 2, byrow = TRUE)) + 
    theme_bw()+
    theme(axis.text.x = element_text(size=6,colour = "black",angle = 90),
          axis.text.y = element_text(size=6,colour = "black"),
          axis.title.y = element_text(size=6,colour = "black"),
          axis.title.x = element_blank(),
          axis.line.y = element_line(size=0.1,colour = "black"), 
          axis.line.x.top = element_line(size=0.1,colour = "black"), 
          axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
          axis.ticks.x =  element_line(linewidth = 0.234),
          axis.ticks.length.y =unit(0.8, "mm"),
          axis.ticks.y =  element_line(linewidth = 0.234),
          panel.background = element_blank(),
          panel.spacing.x = unit(1, "mm"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0.5),
          strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0.5),
          strip.text.y.left = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0),
          strip.background =  element_blank(),
          strip.placement = "top",
          legend.key.height = unit(2,"mm"),
          legend.key.width = unit(2,"mm"),
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.position = "bottom",
          plot.margin=grid::unit(c(1,2,0,2), "mm")
    )
  return(out)
}
p.fusion=p.info.function(myinfo="fusion_detail",mylevel=c("WT","FGFR1","FGFR2","FGFR3","NTRK1","NTRK2","NTRK3",
                                                          "MET","RET","ROS1","ALK"))
p.gender=p.info.function(myinfo="Gender",mylevel = c("Male","Female"))
p.pop=p.info.function(myinfo="Population",mylevel = c("c.Others","a.European","b.Asian"))
p.tp53=p.info.function(myinfo="TP53",mylevel = c("WT","Mut"))
p.kras=p.info.function(myinfo="KRAS",mylevel = c("WT","Other_mut","Q61H","G12D","G12V","G12C"))
p.egfr=p.info.function(myinfo="EGFR",mylevel = c("WT","Other_mut","E479_A483del","L591R"))
p.clustercplx=p.info.function(myinfo="clustercplx",mylevel = c("no_CGR","Unassigned","Multiple",
                                                               "6 Hourglass","5 Large gain","4 Micronuclei",
                                                               "3 Large loss","2 BFB","1 ecDNA"))
p.nonclustercplx=p.info.function(myinfo="nonclustercplx",mylevel = c("no_CGR","Multiple",
                                                                     "cycle_templated_ins","complex_unclear","chromoplexy"))
p=grid.arrange(p.gender,p.pop,p.clustercplx,p.nonclustercplx,
          p.egfr,p.kras,p.tp53,p.fusion,#nullGrob(),
          nrow = 2)
ggsave(plot=p,file.path("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/scratch/24-smoking-sbs4/plot",
                        "info_diff_smokingsbs.pdf"),width = 7,height = 6)


############## violin plot of SBS4 ############
mytheme <- theme(axis.text.x = element_text(size=6,colour = "black",angle = 0),
                 axis.text.y = element_text(size=6,colour = "black"),
                 axis.title.y = element_text(size=6,colour = "black"),
                 axis.title.x = element_blank(),
                 axis.line.y = element_line(size=0.1,colour = "black"), 
                 axis.line.x.top = element_line(size=0.1,colour = "black"), 
                 axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
                 axis.ticks.x =  element_line(linewidth = 0.234),
                 axis.ticks.length.y =unit(0.8, "mm"),
                 axis.ticks.y =  element_line(linewidth = 0.234),
                 panel.background = element_blank(),
                 panel.spacing.x = unit(1, "mm"),
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 plot.title = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0.5),
                 strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0.5),
                 strip.text.y.left = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0),
                 strip.background =  element_blank(),
                 strip.placement = "top",
                 legend.key.height = unit(2,"mm"),
                 legend.key.width = unit(2,"mm"),
                 legend.title = element_blank(),
                 legend.text = element_text(size=6),
                 legend.position = "bottom",
                 plot.margin=grid::unit(c(1,2,0,2), "mm")
)

ggplot(data=merge_sample_info_all %>% 
         mutate(SBS4_mdf=ifelse(SBS4==0,0.1,SBS4)) %>%
         filter(Smoking != "c.Unknown") %>%
         mutate(Smoking=factor(Smoking,levels=c("a.Never smoker","b.Smoker"))))+
  geom_violin(aes(x=Smoking,y=SBS4_mdf,fill=Smoking),color=NA)+
  geom_boxplot(aes(x=Smoking,y=SBS4_mdf),outlier.shape = NA,width=0.2)+
  scale_y_log10(breaks=c(0.1,1,100,10000,1000000),
                labels=c(0.1,1,100,10000,1000000))+
  scale_fill_manual(breaks = c("a.Never smoker","b.Smoker"),
                    values = c("#00a800","#FF6EC7"),
                    labels = c("Never smoker","Smoker"))+
  ylab("SBS4")+
  theme_bw()+
  mytheme
ggsave(file.path("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/scratch/24-smoking-sbs4/plot",
                        "sbs4_value_smoking_violine.pdf"),width = 2,height = 2)


# ggplot(
#   data = merge_sample_info_all %>%
#     mutate(SBS4_mdf = ifelse(SBS4 == 0, 0.1, SBS4)) %>%
#     filter(Smoking != "c.Unknown") %>%
#     mutate(Smoking = factor(Smoking, levels = c("a.Never smoker", "b.Smoker")))
# ) +
#   geom_density(aes(x = SBS4_mdf, fill = Smoking), alpha = 0.4) +
#   scale_x_log10(breaks = c(0.1, 1, 100, 10000, 1000000),
#                 labels = c(0.1, 1, 100, 10000, 1000000)) +
#   scale_fill_manual(
#     breaks = c("a.Never smoker", "b.Smoker"),
#     values = c("#00a800", "#FF6EC7"),
#     labels = c("Never smoker", "Smoker")
#   ) +
#   xlab("SBS4") + ylab("Density") +
#   theme_bw()+
#   mytheme
# ggsave(file.path("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/scratch/24-smoking-sbs4/plot",
#                  "sbs4_value_smoking_density.pdf"),width = 2,height = 2)
# 
# 
# merge_sample_info_all %>%
#   mutate(SBS4_mdf = ifelse(SBS4 == 0, 0.1, SBS4)) %>%
#   filter(Smoking != "c.Unknown") %>%
#   mutate(Smoking = factor(Smoking, levels = c("a.Never smoker", "b.Smoker"))) %>%
#   ggplot(aes(x = SBS4_mdf, fill = Smoking)) +
#   geom_histogram(
#     data = . %>% filter(Smoking == "a.Never smoker"),
#     aes(y = after_stat(count)),
#     bins = 30,
#     position = "identity"
#   ) +
#   geom_histogram(
#     data = . %>% filter(Smoking == "b.Smoker"),
#     aes(y = -after_stat(count)),
#     bins = 30,
#     position = "identity"
#   ) +
#   scale_x_log10(
#     breaks = c(0.1, 1, 100, 10000, 1000000),
#     labels = c(0.1, 1, 100, 10000, 1000000)
#   ) +
#   scale_y_continuous(
#     breaks = c(-100, 0, 100,400,800),
#     labels = abs(c(-100, 0, 100,400,800))
#   ) +
#   scale_fill_manual(
#     values = c("a.Never smoker" = "#00a800", "b.Smoker" = "#FF6EC7"),
#     labels = c("Never smoker", "Smoker")
#   ) +
#   xlab("SBS4") +
#   ylab("Count") +
#   geom_hline(yintercept = 0, color = "black") +
#   theme_bw() +
#   mytheme+
#   theme(axis.title.x = element_text(size=6,colour = "black"))
# ggsave(file.path("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/scratch/24-smoking-sbs4/plot",
#                  "sbs4_value_smoking_hist.pdf"),width = 2,height = 2)
