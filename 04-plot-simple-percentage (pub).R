#### 2023-11-15
# adapted from xiaoming
# input data (rely on mechanism, Gender, mechanism_status to plot):
# sample        mechanism     mechanism_status    somking_status  population  Gender
# IGC-02-1001   3 Large loss  Y                   smoker          EUR         Female
# do not generate multiple signature


### library
library(dplyr)
library(data.table)
library("ggplot2")
library(tidyr)
library(ggpubr)
library(purrr)


###function
p.adjust_function <- function(df,col=3,p.bfrn,p.fdr,name.bfrn,name.fdr){
  df[,p.bfrn] <- p.adjust(df[,col],method = "bonferroni",n=nrow(df))
  df[,p.fdr] <- p.adjust(df[,col],method = "BH")
  names(df)[c(p.bfrn,p.fdr)] <- c(name.bfrn,name.fdr)
  return(df)
}

generate_data <- function(all_samples=all_sample,
                          cluster=cluster,
                          category_name="All Samples"){
  
  
  all_sample_count=nrow(all_samples)
  
  cluster=cluster[cluster$sample %in% all_samples$Tumor_Barcode,]
  sample_with_signature_count=length(unique(cluster$sample))
  
  category_count=data.frame(category=category_name, total_sample=all_sample_count, signature_sample=sample_with_signature_count)
  category_count$percentage=category_count$signature_sample/category_count$total_sample
  category_count$category_label=paste0(category_count$category," ","(",category_count$total_sample,")")
  
  unique_sample_with_signature=cluster %>% count(sample, mechanism)
  # unique_sample_with_signature=unique(cluster[,c("sample", "mechanism")])
  all_sample_id=as.data.frame(all_samples[,c("Tumor_Barcode")])
  colnames(all_sample_id)=c("sample")
  all_sample_signature=merge(unique_sample_with_signature, all_sample_id,by.x=c("sample"),by.y=c("sample"),all.y=T)
  all_sample_signature$mechanism=ifelse(is.na(all_sample_signature$mechanism),"Non-CGR",all_sample_signature$mechanism)
  all_sample_signature$category=category_name
  all_sample_signature=merge(all_sample_signature, unique(category_count[c("category","category_label","percentage")]),by.x="category",by.y="category")
  
  each_sample_signature_count=as.data.frame(all_sample_signature %>% group_by(sample) %>% count())
  each_sample_signature_count$percent=1/each_sample_signature_count$n
  
  all_sample_signature=merge(all_sample_signature, each_sample_signature_count[c("sample","percent")],by.x="sample",by.y="sample")
  all_sample_signature$mech_all=""
  for (i in 1:nrow(all_sample_signature)){
    all_sample_signature$mech_all[i]=paste(sort(unique(all_sample_signature[all_sample_signature$sample==all_sample_signature$sample[i],]$mechanism)),collapse=";")
  }
  all_sample_signature <- all_sample_signature %>% group_by(sample) %>% mutate(proportion=n/sum(n))
  all_sample_signature <- all_sample_signature %>% group_by(sample) %>% mutate(sum=sum(n))
  # add cluster
  all_sample_signature <- merge(all_sample_signature,mycluster[,c("sample","consensusClass")],by.x = "sample",by.y = "sample",all.x = T)
  return(all_sample_signature)
}

plot.function <- function(all_sample_signature,errorbar.df,variable.name){
  category.label <- unique(all_sample_signature$category_label)
  names(category.label) <- unique(all_sample_signature$category)
  all_sample_signature$n <- ifelse(all_sample_signature$mechanism=="Non-CGR",1,all_sample_signature$n)
  p <- ggplot(all_sample_signature,aes(x=n,y=sample,fill=mechanism))+
    geom_bar(stat = "identity",position = "stack",width = 1.2)+
    facet_wrap(category~.,
               labeller = labeller(category=category.label),
               scale="free", #free
               drop=TRUE,strip.position = "top",nrow=1)+
    scale_fill_manual(values = c("Del1"="#D6EAF8",
                                 "Del2"="#85C1E9", 
                                 "Del3"="#5DADE2", 
                                 "TD1"="#FADBD8",
                                 "TD2"="#F1948A",
                                 "Fb inv"="#14CCCC",
                                 "Large intra"="#F1C40F",
                                 "Tra"="#A569BD",
                                 "Non-CGR"="lightgrey",
                                 "Unassigned"="black"
    ))+
    theme(panel.spacing = unit(0, "lines"), 
          axis.title=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.ticks.y=element_blank(), 
          axis.text.x =element_text(size = 5,color = "black",angle = 90),
          axis.text.y =element_blank(), 
          legend.position="bottom", 
          legend.text = element_text(size = 6,color = "black"),
          legend.title = element_text(size = 6,color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border=element_blank(), 
          axis.line.x=element_blank(), 
          axis.line.y=element_blank(), 
          strip.text.x = element_text(angle=90,size = 6,color = "black"), 
          strip.background = element_blank())
  p2 <-  ggplot(all_sample_signature %>% 
                  select(sample,category,category_label,consensusClass,percentage,mech_all) %>% 
                  unique() %>%
                  mutate(consensusClass=as.factor(consensusClass))
                )+
    geom_bar(aes(x=1,y=sample,fill=consensusClass),stat = "identity",position = "fill",width = 1.2)+
    facet_wrap(category~.,
               labeller = labeller(category=category.label),
               scale="free",drop=TRUE,strip.position = "bottom",nrow=1)+
    scale_fill_manual(values = c("1"="#F1948A",
                                 "2"="#85C1E9",
                                 "3"="#A569BD",
                                 "4"="#FADBD8",
                                 "0"="#EDEDED"
    ))+
    theme(panel.spacing = unit(0, "lines"), 
          axis.title=element_blank(), 
          axis.ticks=element_blank(), 
          axis.text=element_blank(), 
          legend.position="bottom", 
          legend.text = element_text(size = 6,color = "black"),
          legend.title = element_text(size = 6,color = "black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border=element_blank(), 
          axis.line.x=element_blank(), 
          axis.line.y=element_blank(), 
          strip.text.x = element_blank(),
          strip.background = element_blank())
  if (is.null(errorbar.df)==FALSE){
    code.nchar <- nchar(variable.name[1])
    errorbar.df$Code1.xposition <- apply(errorbar.df, 1,function(x) which(variable.name == x["Code1"]))
    errorbar.df$Code2.xposition <- apply(errorbar.df, 1,function(x) which(variable.name == x["Code2"]))
    errorbar.df$yposition <- apply(errorbar.df, 1, function(x) code.nchar+1-grep(TRUE,strsplit(x[1], "")[[1]] != strsplit(x[2], "")[[1]]))
    errorbar.df$color=ifelse(errorbar.df$value<0.1,"red","black")
    p.error <- ggplot(errorbar.df) +
      geom_errorbar(aes(xmin=Code1.xposition+0.1,xmax=Code2.xposition-0.1,y=yposition),size=0.1189)+
      geom_text(aes(x=(Code1.xposition+Code2.xposition)/2,y=yposition,label=signif(value,2),color=color),size=2,vjust=-1)+
      scale_x_continuous(labels=variable.name,breaks=seq(1,length(variable.name),1),limits = c(1,length(variable.name)))+
      scale_color_manual(breaks = c("black","red"),values = c("black","red"))+
      facet_grid(FDR~.)+
      # facet_grid(Bonferroni~.)+
      theme(panel.spacing = unit(0, "lines"), 
            axis.title=element_blank(), 
            axis.ticks=element_blank(), 
            axis.text.x =element_text(size = 6,color = "black",angle = 90),
            axis.text.y =element_blank(), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            panel.border=element_blank(), 
            axis.line.x=element_blank(), 
            axis.line.y=element_blank(), 
            legend.position = "none")
    p <- ggpubr::ggarrange(p.error,p,p2,nrow = 3,heights = c(5,1,1))
  }
  return(p)
}
### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
hq.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/HQ_samples.csv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
simple.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_simple_sv_with_signature.tsv")
outpath <- file.path(scratch.path,"05-sigprofiler","1217_samples_hg38","plot")
cluster_path <- file.path(scratch.path,"07-cluster/1217_samples_hg38/Sherlock_SV49_simple/8signature_myannotation/hc_pearson/sig_cluster4.tsv")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

### read
# all sample
all_sample <- read.delim(sampleinfo.path, header=T, sep="\t", stringsAsFactors = F) %>% 
  dplyr::select(Tumor_Barcode, Assigned_Population, Gender, Smoking,high_quality,Histology) 
histology <- read.csv(histology.path)
all_sample <- merge(all_sample,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, ) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")
# simple SV signature
simple=read.delim(simple.path, header=T, sep="\t", stringsAsFactors = F) %>%
  setnames(.,old=c("SAMPLE","Signature"),new=c("sample","mechanism"))
# sbs4
sbs4 <- read.delim2(sbs4.path)
# mutation
mycluster <- read.delim(cluster_path,row.names = 1) %>% mutate(sample=row.names(.))
tp53_mutated_samples=read.delim(tp53.mut.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
kras_mutated_samples=read.delim(kras.mut.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F)  %>% pull(V1)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
# add information
all_sample <- all_sample %>% mutate(tp53=ifelse(Tumor_Barcode %in% tp53_mutated_samples,"mut","wt"))
all_sample <- all_sample %>% mutate(kras=ifelse(Tumor_Barcode %in% kras_mutated_samples,"mut","wt"))
all_sample <- all_sample %>% mutate(egfr=ifelse(Tumor_Barcode %in% c(egfr_mutated_E479_A483del_samples,egfr_mutated_L591R_samples,egfr_mutated_others_samples),"mut","wt"))
all_sample <- all_sample %>% mutate(egfr.detail=ifelse(Tumor_Barcode %in% c(egfr_mutated_L591R_samples),"L591R",
                                                       ifelse(Tumor_Barcode %in% c(egfr_mutated_E479_A483del_samples),"A483del",
                                                              ifelse(Tumor_Barcode %in% c(egfr_mutated_others_samples),"Others","wt"))))
all_sample <- all_sample %>% mutate(egfr.kras=ifelse(egfr=="mut"&kras=="mut","mut-mut",
                                                     ifelse(egfr=="mut"&kras=="wt","mut-wt",
                                                            ifelse(egfr=="wt"&kras=="mut","wt-mut","wt-wt"))
))
all_sample <- all_sample %>% mutate(SP_Group_new=ifelse(Assigned_Population=="EUR" & Smoking=="Smoker","S_U",
                                              ifelse(Assigned_Population=="EUR" & Smoking=="Non-Smoker", "N_U",
                                                     ifelse(Assigned_Population=="EAS" & Smoking=="Non-Smoker", "N_A","Others"))
))
all_sample <- merge(all_sample,sbs4)
all_sample <- all_sample %>% mutate(SP_Group_SBS4=ifelse(Assigned_Population=="EUR" & Smoking=="Smoker" & SBS4>0,"S_U",
                                                        ifelse(Assigned_Population=="EUR" & Smoking=="Non-Smoker" & SBS4==0, "N_U",
                                                               ifelse(Assigned_Population=="EAS" & Smoking=="Non-Smoker" & SBS4==0, "N_A","Others"))
))
all_sample <- all_sample %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
######################################
#### simple SV ####
######################################
#### choose high quality
all_sample_input=all_sample %>% filter(high_quality==1)
all_sample_input=all_sample
all_sample_input=all_sample %>% filter(Histology %in% c("Adenocarcinoma"))
all_sample_input=all_sample %>% filter(Histology %in% c("Adenocarcinoma")) %>% filter(high_quality==1)

# 1. All samples
all_sample_data=generate_data(all_samples=all_sample_input,cluster=simple,category_name="All Samples")
# 2. EUR smoker
eur_smoking_data=generate_data(all_samples=all_sample_input %>% filter(SP_Group_new=="S_U"),cluster=simple,category_name="EUR Smoker")
# 3. EUR non-smoker
eur_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(SP_Group_new=="N_U"),cluster=simple,category_name="EUR Never Smoker")
# 4. EAS non-smoker
eas_nonsmokint_data=generate_data(all_samples=all_sample_input %>% filter(SP_Group_new=="N_A"),cluster=simple,category_name="EAS Never Smoker")
# 5. Female
female_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Female"),cluster=simple,category_name="Females")
# 6. Males
male_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Male"),cluster=simple,category_name="Males")
# 7.1. TP53
tp53wt_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "mut"),cluster=simple,category_name="TP53 WT")
tp53mut_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "wt"),cluster=simple,category_name="TP53 Mut")
# 7.2. TP53 smoking and nonsmoking
tp53wt_smoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "wt" & Smoking == "Smoker"),cluster=simple,category_name="TP53 WT Smoker")
tp53mut_smoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "mut" & Smoking == "Smoker"),cluster=simple,category_name="TP53 Mut Smoker")
tp53wt_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "wt" & Smoking == "Non-Smoker"),cluster=simple,category_name="TP53 WT Smoker")
tp53mut_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "mut" & Smoking == "Non-Smoker"),cluster=simple,category_name="TP53 Mut Smoker")
# 8.1. EGFR
egfr.L591R_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="L591R"),cluster=simple,category_name="EGFR L591R")
egfr.E479_A483del_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="A483del"),cluster=simple,category_name="EGFR Exon 19 del")
egfr.Others_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="Others"),cluster=simple,category_name="EGFR Others")
egfr.wt_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="wt"),cluster=simple,category_name="EGFR WT")
# 8.2. EGFR smoking and nonsmoking
egfr.L591R_smoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="L591R" & Smoking == "Smoker"),cluster=simple,category_name="EGFR L591R Smoker")
egfr.E479_smoking_A483del_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="A483del" & Smoking == "Smoker"),cluster=simple,category_name="EGFR Exon 19 del Smoker")
egfr.Others_smoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="Others" & Smoking == "Smoker"),cluster=simple,category_name="EGFR Others Smoker")
egfr.wt_smoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="wt" & Smoking == "Smoker"),cluster=simple,category_name="EGFR WT Smoker")
egfr.L591R_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="L591R" & Smoking == "Non-Smoker"),cluster=simple,category_name="EGFR L591R Never Smoker")
egfr.E479_nonsmoking_A483del_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="A483del" & Smoking == "Non-Smoker"),cluster=simple,category_name="EGFR Exon 19 del Never Smoker")
egfr.Others_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="Others" & Smoking == "Non-Smoker"),cluster=simple,category_name="EGFR Others Never Smoker")
egfr.wt_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="wt" & Smoking == "Non-Smoker"),cluster=simple,category_name="EGFR WT Never Smoker")
# 9.1. KRAS
kraswt_data=generate_data(all_samples=all_sample_input %>% filter(kras=="wt"),cluster=simple,category_name="KRAS WT")
krasmut_data=generate_data(all_samples=all_sample_input %>% filter(kras=="mut"),cluster=simple,category_name="KRAS Mut")
# 9.2. KRAS smoking and nonsmoking
kraswt_smoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="wt" & Smoking == "Smoker"),cluster=simple,category_name="KRAS WT Smoker")
krasmut_smoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="mut" & Smoking == "Smoker"),cluster=simple,category_name="KRAS Mut Smoker")
kraswt_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="wt" & Smoking == "Non-Smoker"),cluster=simple,category_name="KRAS WT Never Smoker")
krasmut_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="mut" & Smoking == "Non-Smoker"),cluster=simple,category_name="KRAS Mut Never Smoker")
# 10. EGFR-KRAS by SP_Group_new, Gender, tp53
ek1111 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"S_U-F-TP53wt-EGFRmutKRASwt") #0
ek1112 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"S_U-F-TP53wt-EGFRwtKRASwt") #22
ek1113 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"S_U-F-TP53wt-EGFRwtKRASmut") #11
ek1121 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"S_U-F-TP53mut-EGFRmutKRASwt") #1
ek1122 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"S_U-F-TP53mut-EGFRwtKRASwt") #47
ek1123 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"S_U-F-TP53mut-EGFRwtKRASmut") #18
ek1211 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"S_U-M-TP53wt-EGFRmutKRASwt") #4
ek1212 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"S_U-M-TP53wt-EGFRwtKRASwt") #52
ek1213 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"S_U-M-TP53wt-EGFRwtKRASmut") #51
ek1221 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"S_U-M-TP53mut-EGFRmutKRASwt") #13
ek1222 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"S_U-M-TP53mut-EGFRwtKRASwt") #138
ek1223 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"S_U-M-TP53mut-EGFRwtKRASmut") #31

ek2111 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"N_U-F-TP53wt-EGFRmutKRASwt") #133
ek2112 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"N_U-F-TP53wt-EGFRwtKRASwt") #202
ek2113 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"N_U-F-TP53wt-EGFRwtKRASmut") #28
ek2121 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"N_U-F-TP53mut-EGFRmutKRASwt") #64
ek2122 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"N_U-F-TP53mut-EGFRwtKRASwt") #49
ek2123 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"N_U-F-TP53mut-EGFRwtKRASmut") #6
ek2211 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"N_U-M-TP53wt-EGFRmutKRASwt") #23
ek2212 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"N_U-M-TP53wt-EGFRwtKRASwt") #51
ek2213 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"N_U-M-TP53wt-EGFRwtKRASmut") #5
ek2221 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"N_U-M-TP53mut-EGFRmutKRASwt") #12
ek2222 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"N_U-M-TP53mut-EGFRwtKRASwt") #12
ek2223 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"N_U-M-TP53mut-EGFRwtKRASmut") #1

ek3111 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"N_A-F-TP53wt-EGFRmutKRASwt") #167
ek3112 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"N_A-F-TP53wt-EGFRwtKRASwt") #82
ek3113 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"N_A-F-TP53wt-EGFRwtKRASmut") #6
ek3121 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"N_A-F-TP53mut-EGFRmutKRASwt") #119
ek3122 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"N_A-F-TP53mut-EGFRwtKRASwt") #30
ek3123 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"N_A-F-TP53mut-EGFRwtKRASmut") #1
ek3211 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"N_A-M-TP53wt-EGFRmutKRASwt") #27
ek3212 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"N_A-M-TP53wt-EGFRwtKRASwt") #22
ek3213 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"N_A-M-TP53wt-EGFRwtKRASmut") #1
ek3221 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"N_A-M-TP53mut-EGFRmutKRASwt") #18
ek3222 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"N_A-M-TP53mut-EGFRwtKRASwt") #6
ek3223 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"N_A-M-TP53mut-EGFRwtKRASmut") #0

# 11. only EGFR/KRAS/TP53 AND smoking
ns.tw.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"NeverSmoker-TP53wt-EGFRwtKRASwt") #370
ns.tw.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"NeverSmoker-TP53wt-EGFRmutKRASwt") #354
ns.tw.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"NeverSmoker-TP53wt-EGFRwtKRASmut") #44
ns.tm.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"NeverSmoker-TP53mut-EGFRwtKRASwt") #105
ns.tm.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"NeverSmoker-TP53mut-EGFRmutKRASwt") #217
ns.tm.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"NeverSmoker-TP53mut-EGFRwtKRASmut") #10
ss.tw.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "wt" & egfr.kras=="wt-wt"),simple,"Smoker-TP53wt-EGFRwtKRASwt") #77
ss.tw.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "wt" & egfr.kras=="mut-wt"),simple,"Smoker-TP53wt-EGFRmutKRASwt") #4
ss.tw.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "wt" & egfr.kras=="wt-mut"),simple,"Smoker-TP53wt-EGFRwtKRASmut") #64
ss.tm.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "mut" & egfr.kras=="wt-wt"),simple,"Smoker-TP53mut-EGFRwtKRASwt") #186
ss.tm.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "mut" & egfr.kras=="mut-wt"),simple,"Smoker-TP53mut-EGFRmutKRASwt") #16
ss.tm.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "mut" & egfr.kras=="wt-mut"),simple,"Smoker-TP53mut-EGFRwtKRASmut") #49

nrow(ns.tw.ew.kw)+nrow(ns.tw.em.kw)+nrow(ns.tw.ew.km)+
  nrow(ns.tm.ew.kw)+nrow(ns.tm.em.kw)+nrow(ns.tm.ew.km)+
  nrow(ss.tw.ew.kw)+nrow(ss.tw.em.kw)+nrow(ss.tw.ew.km)+
  nrow(ss.tm.ew.kw)+nrow(ss.tm.em.kw)+nrow(ss.tm.ew.km)

##################### Plot :freq #################################################
variable.name <- c("ns.tw.ew.kw","ns.tw.em.kw","ns.tw.ew.km","ns.tm.ew.kw","ns.tm.em.kw","ns.tm.ew.km",
                   "ss.tw.ew.kw",
                   "ss.tw.em.kw",
                   "ss.tw.ew.km","ss.tm.ew.kw","ss.tm.em.kw","ss.tm.ew.km")
full.name <- c("NeverSmoker-TP53wt-EGFRwtKRASwt","NeverSmoker-TP53wt-EGFRmutKRASwt","NeverSmoker-TP53wt-EGFRwtKRASmut",
               "NeverSmoker-TP53mut-EGFRwtKRASwt","NeverSmoker-TP53mut-EGFRmutKRASwt","NeverSmoker-TP53mut-EGFRwtKRASmut",
               "Smoker-TP53wt-EGFRwtKRASwt",
               "Smoker-TP53wt-EGFRmutKRASwt",
               "Smoker-TP53wt-EGFRwtKRASmut",
               "Smoker-TP53mut-EGFRwtKRASwt","Smoker-TP53mut-EGFRmutKRASwt","Smoker-TP53mut-EGFRwtKRASmut")
names(full.name) <- variable.name
all_sample_signature=rbind(ns.tw.ew.kw,ns.tw.em.kw,ns.tw.ew.km,ns.tm.ew.kw,ns.tm.em.kw,ns.tm.ew.km,
                           ss.tw.ew.kw,
                           ss.tw.em.kw,
                           ss.tw.ew.km,ss.tm.ew.kw,ss.tm.em.kw,ss.tm.ew.km)

plot.groups <- full.name
mysig=c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra")

all_sample_signature=rbind(ns.tw.ew.kw %>% mutate(category.short="ns.tw.ew.kw"),
                           ns.tw.em.kw %>% mutate(category.short="ns.tw.em.kw"),
                           ns.tw.ew.km %>% mutate(category.short="ns.tw.ew.km"),
                           ns.tm.ew.kw %>% mutate(category.short="ns.tm.ew.kw"),
                           ns.tm.em.kw %>% mutate(category.short="ns.tm.em.kw"),
                           ns.tm.ew.km %>% mutate(category.short="ns.tm.ew.km"),
                           ss.tw.ew.kw %>% mutate(category.short="ss.tw.ew.kw"),
                           ss.tw.em.kw %>% mutate(category.short="ss.tw.em.kw"),
                           ss.tw.ew.km %>% mutate(category.short="ss.tw.ew.km"),
                           ss.tm.ew.kw %>% mutate(category.short="ss.tm.ew.kw"),
                           ss.tm.em.kw %>% mutate(category.short="ss.tm.em.kw"),
                           ss.tm.ew.km %>% mutate(category.short="ss.tm.ew.km"))
freq.df <- merge(all_sample_input,
                 all_sample_signature %>% 
                   select("sample","category","category.short","mechanism","n") %>% 
                   tidyr::spread(.,mechanism,n),
                 all.x = T,
                 by.x = "Tumor_Barcode",by.y = "sample",
) %>% filter(is.na(category)==FALSE) %>%
  dplyr::select("Tumor_Barcode","category","category.short","Smoking_SBS4",mysig) %>%
  replace(is.na(.), 0) %>%
  tidyr::gather(.,Sig,Num,Del1:Tra) %>%
  mutate(category=factor(category,levels = plot.groups))

# test
sample_sizes <- freq.df %>%
  group_by(category.short, Sig) %>%
  summarise(n = n(), .groups = "drop") %>% 
  mutate(category.short=factor(category.short,levels=variable.name)) %>%
  arrange(category.short)
filtered_groups <- sample_sizes %>% filter(n>=20) %>% pull(category.short) %>% unique() %>% as.character()
all_pairs <- combn(filtered_groups, 
                   2, simplify = TRUE) # Use combn to generate all possible pairs
test_pairs <- all_pairs[, apply(all_pairs, 2, function(pair) sum(strsplit(pair[1], "")[[1]] != strsplit(pair[2], "")[[1]]) == 1)] # Filter pairs with one character difference
test_pairs <- matrix(test_pairs, ncol = 2, byrow = TRUE) 

# Initialize results list
test_results <- list()
# Unique signatures
signatures <- unique(freq.df$Sig)
# Iterate over each signature
for (sig in signatures) {
  # Subset data by signature
  sig_df <- freq.df %>% filter(Sig == sig)
  # Iterate over each pair
  for (i in 1:nrow(test_pairs)) {
    group1 <- test_pairs[i, 1]
    group2 <- test_pairs[i, 2]
    # Get data for both groups
    df1 <- sig_df %>% filter(category.short == group1)
    df2 <- sig_df %>% filter(category.short == group2)
    # Check if both groups have at least 20 samples
    if (nrow(df1) >= 20 & nrow(df2) >= 20) {
      # Perform t-test
      test <- wilcox.test(df1$Num, df2$Num)
      # Store result
      test_results[[length(test_results) + 1]] <- data.frame(
        Signature = sig,
        Group1 = group1,
        Group2 = group2,
        n1 = nrow(df1),
        n2 = nrow(df2),
        p.value = test$p.value,
        mean1 = mean(df1$Num),
        mean2 = mean(df2$Num),
        median1 = median(df1$Num),
        median2 = median(df2$Num)
      )
    }
  }
}
# Combine all results into one data frame
test_df <- do.call(rbind, test_results)
# FDR adjust p-values within each Signature
test_df <- test_df %>%
  group_by(Signature) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  ungroup()
# add n
freq.df <- freq.df %>%
  left_join(sample_sizes, by = c("category.short", "Sig")) %>%
  mutate(category.label = paste0(category.short, "\n(n=", n, ")")) %>% 
  mutate(category.short=factor(category.short,levels=variable.name[variable.name %in% filtered_groups])) %>% 
  arrange(category.short)
# Ensure category.label order matches plot
label_map <- freq.df %>% distinct(category.short, category.label)
# Add pair_id and structure for ggpubr
test_pairs_df <- as_tibble(test_pairs) %>%
  setnames(.,old=c("V1","V2"),new=c("Group1","Group2")) %>%
  mutate(pair_id = row_number())

p_label_df <- test_df %>%
  filter(p.adj < 0.05) %>%
  inner_join(test_pairs_df, by = c("Group1", "Group2")) %>%
  left_join(label_map, by = c("Group1" = "category.short")) %>%
  setnames(.,old=c("category.label"),new=c("group1")) %>%
  left_join(label_map, by = c("Group2" = "category.short")) %>%
  setnames(.,old=c("category.label"),new=c("group2")) %>%
  mutate(
    Sig = factor(Signature, levels = mysig),
    label = signif(p.adj, 2),
    base_y = map2_dbl(group1, group2, ~max(freq.df %>% 
                                             filter(category.label %in% c(.x, .y) & Sig == Signature) %>% 
                                             pull(Num), na.rm = TRUE)),
    y.position =  pair_id) %>%
  dplyr::select(Sig, group1, group2, y.position, label)
p_label_df$Sig <- factor(p_label_df$Sig, levels = mysig)
p_label_df <- p_label_df %>% as.data.frame()


# Prepare plot data with zero replaced for log10
plot_df <- freq.df %>%
  filter(n >= 20) %>%
  mutate(Num = ifelse(Num == 0, 0.5, Num),
         Sig = factor(Sig,levels = mysig),
         category.short = factor(category.short,filtered_groups)) %>% 
  arrange(category.short)
plot_df <- plot_df %>% mutate(category.label=factor(category.label,levels=unique(plot_df$category.label)))

ggplot() + 
  # Violin plot
  geom_violin(data=plot_df,aes(x = category.label, y = Num, fill = Smoking_SBS4),
              trim = TRUE, color = NA) +
  scale_fill_manual(values = c("Non-Smoker-noSBS4" = "#00a800",
                               "Smoker-SBS4" = "#FF6EC7")) +
  
  # Boxplot overlay
  geom_boxplot(data=plot_df,aes(x = category.label, y = Num),
               width = 0.2, size = 0.1, outlier.shape = NA) +
  
  # Add FDR annotations via stat_pvalue_manual
  stat_pvalue_manual(
    data = p_label_df,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.1,
    size = 1.8,
    bracket.size = 0.2
  ) +
  # Facet by Signature
  facet_wrap(Sig ~ ., scales = "free", nrow = 8) +
  # Axes and styling
  scale_x_discrete(drop = FALSE) +
  scale_y_log10(breaks = c(0.5, 1, 10, 100)) +
  ylab("Number of signature") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.line.y = element_line(size = 0.1, colour = "black"), 
    axis.line.x.top = element_line(size = 0.1, colour = "black"), 
    axis.line.x.bottom = element_line(size = 0.1, colour = "black"), 
    axis.ticks.length.x = unit(0.8, "mm"),
    axis.ticks.length.y = unit(0.8, "mm"),
    panel.background = element_blank(),
    panel.spacing.x = unit(0.2, "cm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text.x = element_text(size = 6, angle = 0, colour = "black", hjust = 1, vjust = 0.5),
    strip.background = element_blank(),
    strip.placement = "top",
    plot.margin = grid::unit(c(0, 2, 0, 2), "mm"),
    legend.key = element_rect(fill = "transparent"),
    legend.key.size = unit(0.15, "cm"),
    legend.text = element_text(size = 6, colour = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -20),
    legend.box.margin = margin(0, 0, 0, 0)
  )
ggsave(file.path(outpath,"simple.egfr.kras.freq.9groups.pdf"),width = 2,height = 25)
ggsave(file.path(outpath,"simple.egfr.kras.freq.9groups.luad.pdf"),width = 2,height = 25)
ggsave(file.path(outpath,"simple.egfr.kras.freq.9groups.luad.hq.pdf"),width = 2,height = 25)
ggsave(file.path(outpath,"simple.egfr.kras.freq.9groups.hq.pdf"),width = 2,height = 25)
##################### Plot1:others #################################################
all_sample_signature=rbind(all_sample_data, 
                           eur_smoking_data,eur_nonsmoking_data,eas_nonsmokint_data,
                           female_data,male_data,
                           tp53wt_data,tp53mut_data,
                           egfr.L591R_data,egfr.E479_A483del_data,egfr.Others_data,egfr.wt_data,
                           kraswt_data,krasmut_data,
                           tp53wt_smoking_data,tp53mut_smoking_data,tp53wt_nonsmoking_data,tp53mut_nonsmoking_data,
                           egfr.wt_smoking_data,
                           # egfr.L591R_smoking_data,
                           egfr.E479_smoking_A483del_data,egfr.Others_smoking_data,
                           egfr.L591R_nonsmoking_data,egfr.E479_nonsmoking_A483del_data,egfr.Others_nonsmoking_data,egfr.wt_nonsmoking_data,
                           kraswt_smoking_data,krasmut_smoking_data,kraswt_nonsmoking_data,krasmut_nonsmoking_data
)
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = c("All Samples",
                                             "EAS Never Smoker","EUR Never Smoker","EUR Smoker",
                                             "Females","Males",
                                             "TP53 WT","TP53 Mut",
                                             "EGFR WT","EGFR L591R","EGFR Exon 19 del","EGFR Others",
                                             "KRAS WT","KRAS Mut",
                                             "TP53 WT Smoker","TP53 Mut Smoker","TP53 WT Never Smoker","TP53 Mut Never Smoker",
                                             "EGFR WT Smoker",
                                             # "EGFR L591R Smoker",
                                             "EGFR Exon 19 del Smoker","EGFR Others Smoker",
                                             "EGFR WT Never Smoker","EGFR L591R Never Smoker","EGFR Exon 19 del Never Smoker","EGFR Others Never Smoker",
                                             "KRAS WT Smoker","KRAS Mut Smoker","KRAS WT Never Smoker","KRAS Mut Never Smoker"
  ))) %>% 
  mutate(consensusClass=factor(consensusClass,levels = c(3,2,1,4,NA)),)
all_sample_signature=all_sample_signature %>% 
  mutate(mech_all=factor(mech_all,levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR")),
         mechanism=factor(mechanism,levels = c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra","Non-CGR"))) %>% 
  arrange(
    consensusClass,
    sum)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
plot.function(all_sample_signature,NULL,NULL)
ggsave(file.path(outpath,"simple.percentage.pdf"),width = 6,height = 4)
#############################  plot2  egfr-kras ###################################
variable.name <- c("ek1111","ek1112","ek1113","ek1121","ek1122","ek1123","ek1211","ek1212","ek1213","ek1221","ek1222","ek1223",
                   "ek2111","ek2112","ek2113","ek2121","ek2122","ek2123","ek2211","ek2212","ek2213","ek2221","ek2222","ek2223",
                   "ek3111","ek3112","ek3113","ek3121","ek3122","ek3123","ek3211","ek3212","ek3213","ek3221","ek3222","ek3223")
full.name <- c("S_U-F-TP53wt-EGFRmutKRASwt","S_U-F-TP53wt-EGFRwtKRASwt","S_U-F-TP53wt-EGFRwtKRASmut",
               "S_U-F-TP53mut-EGFRmutKRASwt","S_U-F-TP53mut-EGFRwtKRASwt","S_U-F-TP53mut-EGFRwtKRASmut",
               "S_U-M-TP53wt-EGFRmutKRASwt","S_U-M-TP53wt-EGFRwtKRASwt","S_U-M-TP53wt-EGFRwtKRASmut",
               "S_U-M-TP53mut-EGFRmutKRASwt","S_U-M-TP53mut-EGFRwtKRASwt","S_U-M-TP53mut-EGFRwtKRASmut",
               "N_U-F-TP53wt-EGFRmutKRASwt","N_U-F-TP53wt-EGFRwtKRASwt","N_U-F-TP53wt-EGFRwtKRASmut",
               "N_U-F-TP53mut-EGFRmutKRASwt","N_U-F-TP53mut-EGFRwtKRASwt","N_U-F-TP53mut-EGFRwtKRASmut",
               "N_U-M-TP53wt-EGFRmutKRASwt","N_U-M-TP53wt-EGFRwtKRASwt","N_U-M-TP53wt-EGFRwtKRASmut",
               "N_U-M-TP53mut-EGFRmutKRASwt","N_U-M-TP53mut-EGFRwtKRASwt","N_U-M-TP53mut-EGFRwtKRASmut",
               "N_A-F-TP53wt-EGFRmutKRASwt","N_A-F-TP53wt-EGFRwtKRASwt","N_A-F-TP53wt-EGFRwtKRASmut",
               "N_A-F-TP53mut-EGFRmutKRASwt","N_A-F-TP53mut-EGFRwtKRASwt","N_A-F-TP53mut-EGFRwtKRASmut",
               "N_A-M-TP53wt-EGFRmutKRASwt","N_A-M-TP53wt-EGFRwtKRASwt","N_A-M-TP53wt-EGFRwtKRASmut",
               "N_A-M-TP53mut-EGFRmutKRASwt","N_A-M-TP53mut-EGFRwtKRASwt","N_A-M-TP53mut-EGFRwtKRASmut")
names(full.name) <- variable.name
all_sample_signature=rbind(ek1111,ek1112,ek1113,ek1121,ek1122,ek1123,ek1211,ek1212,ek1213,ek1221,ek1222,ek1223,
                           ek2111,ek2112,ek2113,ek2121,ek2122,ek2123,ek2211,ek2212,ek2213,ek2221,ek2222,ek2223,
                           ek3111,ek3112,ek3113,ek3121,ek3122,ek3123,ek3211,ek3212,ek3213,ek3221,ek3222,NULL)
# order sample
all_sample_signature$consensusClass <- ifelse(is.na(all_sample_signature$consensusClass),0,all_sample_signature$consensusClass)
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = full.name)) %>% 
  dplyr::arrange(.,category,desc(percentage),mech_all) 
all_sample_signature=all_sample_signature %>% 
  mutate(mech_all=factor(mech_all,
                         levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR"))) %>% 
  mutate(consensusClass=factor(consensusClass,levels = c(4,1,2,3,0))) %>%
  arrange(desc(sum),consensusClass,desc(percentage),mech_all)
  # arrange(consensusClass,desc(percentage),mech_all)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
# fisher's test: cluster1/notcluster1 cluster2/notcluster2 ...
all_pairs <- combn(variable.name, 2, simplify = TRUE) # Use combn to generate all possible pairs
test_pairs <- all_pairs[, apply(all_pairs, 2, function(pair) sum(strsplit(pair[1], "")[[1]] != strsplit(pair[2], "")[[1]]) == 1)] # Filter pairs with one character difference
test_pairs <- matrix(test_pairs, ncol = 2, byrow = TRUE) #t
fisher_result <- t(apply(test_pairs, 1, function(x){
  v1 <- full.name[x[1]]
  v2 <- full.name[x[2]]
  
  v1.cluster1 <- all_sample_signature %>% filter(category==v1 & consensusClass ==1) %>% select(sample) %>% unique() %>% nrow()
  v1.cluster0234 <- all_sample_signature %>% filter(category==v1 & consensusClass !=1) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster1 <- all_sample_signature %>% filter(category==v2 & consensusClass ==1) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster0234 <- all_sample_signature %>% filter(category==v2 & consensusClass !=1) %>% select(sample) %>% unique() %>% nrow()
  
  v1.cluster2 <- all_sample_signature %>% filter(category==v1 & consensusClass ==2) %>% select(sample) %>% unique() %>% nrow()
  v1.cluster0134 <- all_sample_signature %>% filter(category==v1 & consensusClass !=2) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster2 <- all_sample_signature %>% filter(category==v2 & consensusClass ==2) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster0134 <- all_sample_signature %>% filter(category==v2 & consensusClass !=2) %>% select(sample) %>% unique() %>% nrow()
  
  v1.cluster3 <- all_sample_signature %>% filter(category==v1 & consensusClass ==3) %>% select(sample) %>% unique() %>% nrow()
  v1.cluster0124 <- all_sample_signature %>% filter(category==v1 & consensusClass !=3) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster3 <- all_sample_signature %>% filter(category==v2 & consensusClass ==3) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster0124 <- all_sample_signature %>% filter(category==v2 & consensusClass !=3) %>% select(sample) %>% unique() %>% nrow()
  
  v1.cluster4 <- all_sample_signature %>% filter(category==v1 & consensusClass ==4) %>% select(sample) %>% unique() %>% nrow()
  v1.cluster0123 <- all_sample_signature %>% filter(category==v1 & consensusClass !=4) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster4 <- all_sample_signature %>% filter(category==v2 & consensusClass ==4) %>% select(sample) %>% unique() %>% nrow()
  v2.cluster0123 <- all_sample_signature %>% filter(category==v2 & consensusClass !=4) %>% select(sample) %>% unique() %>% nrow()
  
  v1.total <- v1.cluster1 + v1.cluster0234
  v2.total <- v2.cluster1 + v2.cluster0234
  
  v1.del1 <- c(all_sample_signature %>% filter(category==v1 & mechanism =="Del1") %>% pull(n),rep(0,2000))[1:v1.total] #rep(0,2000) used to add 0 to samples without specific signature
  v2.del1 <- c(all_sample_signature %>% filter(category==v2 & mechanism =="Del1") %>% pull(n),rep(0,2000))[1:v2.total]
  
  v1.del2 <- c(all_sample_signature %>% filter(category==v1 & mechanism =="Del2") %>% pull(n),rep(0,2000))[1:v1.total]
  v2.del2 <- c(all_sample_signature %>% filter(category==v2 & mechanism =="Del2") %>% pull(n),rep(0,2000))[1:v2.total]
  
  v1.del3 <- c(all_sample_signature %>% filter(category==v1 & mechanism =="Del3") %>% pull(n),rep(0,2000))[1:v1.total]
  v2.del3 <- c(all_sample_signature %>% filter(category==v2 & mechanism =="Del3") %>% pull(n),rep(0,2000))[1:v2.total]
  
  v1.td1 <- c(all_sample_signature %>% filter(category==v1 & mechanism =="TD1") %>% pull(n),rep(0,2000))[1:v1.total]
  v2.td1 <- c(all_sample_signature %>% filter(category==v2 & mechanism =="TD1") %>% pull(n),rep(0,2000))[1:v2.total]
  
  v1.td2 <- c(all_sample_signature %>% filter(category==v1 & mechanism =="TD2") %>% pull(n),rep(0,2000))[1:v1.total]
  v2.td2 <- c(all_sample_signature %>% filter(category==v2 & mechanism =="TD2") %>% pull(n),rep(0,2000))[1:v2.total]
  
  v1.fbinv <- c(all_sample_signature %>% filter(category==v1 & mechanism =="Fb inv") %>% pull(n),rep(0,2000))[1:v1.total]
  v2.fbinv <- c(all_sample_signature %>% filter(category==v2 & mechanism =="Fb inv") %>% pull(n),rep(0,2000))[1:v2.total]
  
  v1.largeintra <- c(all_sample_signature %>% filter(category==v1 & mechanism =="Large intra") %>% pull(n),rep(0,2000))[1:v1.total]
  v2.largeintra <- c(all_sample_signature %>% filter(category==v2 & mechanism =="Large intra") %>% pull(n),rep(0,2000))[1:v2.total]
  
  v1.tra <- c(all_sample_signature %>% filter(category==v1 & mechanism =="Tra") %>% pull(n),rep(0,2000))[1:v1.total]
  v2.tra <- c(all_sample_signature %>% filter(category==v2 & mechanism =="Tra") %>% pull(n),rep(0,2000))[1:v2.total]
  
  if (v1.total >= 10 & v2.total >= 10){
    fisher1vs0234 <- fisher.test(matrix(c(v1.cluster1,v1.cluster0234,v2.cluster1,v2.cluster0234),nrow = 2,ncol = 2))
    fisher2vs0134 <- fisher.test(matrix(c(v1.cluster2,v1.cluster0134,v2.cluster2,v2.cluster0134),nrow = 2,ncol = 2))
    fisher3vs0124 <- fisher.test(matrix(c(v1.cluster3,v1.cluster0124,v2.cluster3,v2.cluster0124),nrow = 2,ncol = 2))
    fisher4vs0123 <- fisher.test(matrix(c(v1.cluster4,v1.cluster0123,v2.cluster4,v2.cluster0123),nrow = 2,ncol = 2))

    t.del1 <- t.test(v1.del1, v2.del1)
    t.del2 <- t.test(v1.del2, v2.del2)
    t.del3 <- t.test(v1.del3, v2.del3)
    t.td1 <- t.test(v1.td1, v2.td1)
    t.td2 <- t.test(v1.td2, v2.td2)
    t.fbinv <- t.test(v1.fbinv, v2.fbinv)
    t.largeintra <- t.test(v1.largeintra, v2.largeintra)
    t.tra <- t.test(v1.tra, v2.tra)
    
    # t.del1 <- wilcox.test(v1.del1,v2.del1,paired=FALSE)
    # t.del2 <- wilcox.test(v1.del2,v2.del2,paired=FALSE)
    # t.del3 <- wilcox.test(v1.del3,v2.del3,paired=FALSE)
    # t.td1 <- wilcox.test(v1.td1,v2.td1,paired=FALSE)
    # t.td2 <- wilcox.test(v1.td2,v2.td2,paired=FALSE)
    # t.fbinv <- wilcox.test(v1.fbinv,v2.fbinv,paired=FALSE)
    # t.largeintra <- wilcox.test(v1.largeintra,v2.largeintra,paired=FALSE)
    # t.tra <- wilcox.test(v1.tra,v2.tra,paired=FALSE)
    
    out <- c(x[1],x[2],v1,v2,v1.total,v2.total,
             v1.cluster1,v1.cluster0234,v2.cluster1,v2.cluster0234,fisher1vs0234$p.value,fisher1vs0234$estimate,
             v1.cluster2,v1.cluster0134,v2.cluster2,v2.cluster0134,fisher2vs0134$p.value,fisher2vs0134$estimate,
             v1.cluster3,v1.cluster0124,v2.cluster3,v2.cluster0124,fisher3vs0124$p.value,fisher3vs0124$estimate,
             v1.cluster4,v1.cluster0123,v2.cluster4,v2.cluster0123,fisher4vs0123$p.value,fisher4vs0123$estimate,
             t.del1$p.value,t.del1$statistic,
             t.del2$p.value,t.del2$statistic,
             t.del3$p.value,t.del3$statistic,
             t.td1$p.value,t.td1$statistic,
             t.td2$p.value,t.td2$statistic,
             t.fbinv$p.value,t.fbinv$statistic,
             t.largeintra$p.value,t.largeintra$statistic,
             t.tra$p.value,t.tra$statistic
             )
  } else {
    out <- c(x[1],x[2],v1,v2,v1.total,v2.total,
             v1.cluster1,v1.cluster0234,v2.cluster1,v2.cluster0234,"Not tested","Not tested",
             v1.cluster2,v1.cluster0134,v2.cluster2,v2.cluster0134,"Not tested","Not tested",
             v1.cluster3,v1.cluster0124,v2.cluster3,v2.cluster0124,"Not tested","Not tested",
             v1.cluster4,v1.cluster0123,v2.cluster4,v2.cluster0123,"Not tested","Not tested",
             "Not tested","Not tested",
             "Not tested","Not tested",
             "Not tested","Not tested",
             "Not tested","Not tested",
             "Not tested","Not tested",
             "Not tested","Not tested",
             "Not tested","Not tested",
             "Not tested","Not tested"
             )
  }
  return(out)
})) %>% as.data.frame()#fisher's test
names(fisher_result) <- c("Code1","Code2","Group1","Group2","Group1.num","Group2.num",
                          "Group1.cluster1.num","Group1.cluster0234.num","Group2.cluster1.num","Group2.cluster0234.num","P.1vs0234","OR.1vs0234",
                          "Group1.cluster2.num","Group1.cluster0134.num","Group2.cluster2.num","Group2.cluster0134.num","P.2vs0134","OR.2vs0134",
                          "Group1.cluster3.num","Group1.cluster0124.num","Group2.cluster3.num","Group2.cluster0124.num","P.3vs0124","OR.3vs0124",
                          "Group1.cluster4.num","Group1.cluster0123.num","Group2.cluster4.num","Group2.cluster0123.num","P.4vs0123","OR.4vs0123",
                          "P.Del1","t.Del1",
                          "P.Del2","t.Del2",
                          "P.Del3","t.Del3",
                          "P.TD1","t.TD1",
                          "P.TD2","t.TD2",
                          "P.Fbinv","t.Fbinv",
                          "P.Largeintra","t.Largeintra",
                          "P.Tra","t.Tra"
                          )
# filter fisher
fisher_result_tested <- fisher_result %>% filter(P.1vs0234 !="Not tested") %>% 
  filter((grepl("EGFRmutKRASwt",Group1)&grepl("EGFRwtKRASmut",Group2))==FALSE) %>% #remove pairs EGFRmutKRASwt and EGFRwtKRASmut
  filter((grepl("S_U",Group1)&grepl("N_A",Group2))==FALSE) %>% #remove pairs S_U and N_A
  mutate(P.1vs0234=as.numeric(P.1vs0234),OR.1vs0234=as.numeric(OR.1vs0234),
         P.2vs0134=as.numeric(P.2vs0134),OR.2vs0134=as.numeric(OR.2vs0134),
         P.3vs0124=as.numeric(P.3vs0124),OR.3vs0124=as.numeric(OR.3vs0124),
         P.4vs0123=as.numeric(P.4vs0123),OR.4vs0123=as.numeric(OR.4vs0123),
         P.Del1=as.numeric(P.Del1),t.Del1=as.numeric(t.Del1),
         P.Del2=as.numeric(P.Del2),t.Del2=as.numeric(t.Del2),
         P.Del3=as.numeric(P.Del3),t.Del3=as.numeric(t.Del3),
         P.TD1=as.numeric(P.TD1),t.TD1=as.numeric(t.TD1),
         P.TD2=as.numeric(P.TD2),t.TD2=as.numeric(t.TD2),
         P.Fbinv=as.numeric(P.Fbinv),t.Fbinv=as.numeric(t.Fbinv),
         P.Largeintra=as.numeric(P.Largeintra),t.Largeintra=as.numeric(t.Largeintra),
         P.Tra=as.numeric(P.Tra),t.Tra=as.numeric(t.Tra)
         )
# adjust p value
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=11,47,48,"Bonferroni.1vs0234","FDR.1vs0234")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=17,49,50,"Bonferroni.2vs0134","FDR.2vs0134")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=23,51,52,"Bonferroni.3vs0124","FDR.3vs0124")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=29,53,54,"Bonferroni.4vs0123","FDR.4vs0123")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=31,55,56,"Bonferroni.Del1","FDR.Del1")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=33,57,58,"Bonferroni.Del2","FDR.Del2")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=35,59,60,"Bonferroni.Del3","FDR.Del3")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=37,61,62,"Bonferroni.TD1","FDR.TD1")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=39,63,64,"Bonferroni.TD2","FDR.TD2")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=41,65,66,"Bonferroni.Fbinv","FDR.Fbinv")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=43,67,68,"Bonferroni.Largeintra","FDR.Largeintra")
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=45,69,70,"Bonferroni.Tra","FDR.Tra")
# bar plot
# pay attention to line 143-144 FDR and Bonferroni
errorbar.df <- fisher_result_tested %>% select("Code1","Code2","Group1","Group2","Group1.num","Group2.num",
                                               "FDR.1vs0234","FDR.2vs0134","FDR.3vs0124","FDR.4vs0123",
                                               "FDR.Del1","FDR.Del2","FDR.Del3","FDR.TD1","FDR.TD2",
                                               "FDR.Fbinv","FDR.Largeintra","FDR.Tra") %>%
  tidyr::gather(.,FDR,value,FDR.1vs0234:FDR.Tra,factor_key = T)
plot.function(all_sample_signature,errorbar.df,variable.name)
ggsave(file.path(outpath,"simple.egfr.kras.percentage.ttest.FDR.pdf"),width = 6,height = 24) #pay attention to line 429-445, choose ttest or ranktest
ggsave(file.path(outpath,"simple.egfr.kras.percentage.ranktest.FDR.pdf"),width = 6,height = 24)  #pay attention to line 429-445, choose ttest or ranktest

# pay attention to line 143-144 FDR and Bonferroni
errorbar.df <- fisher_result_tested %>% dplyr::select("Code1","Code2","Group1","Group2",
                                               "Bonferroni.1vs0234","Bonferroni.2vs0134","Bonferroni.3vs0124","Bonferroni.4vs0123",
                                               "Bonferroni.Del1","Bonferroni.Del2","Bonferroni.Del3","Bonferroni.TD1","Bonferroni.TD2",
                                               "Bonferroni.Fbinv","Bonferroni.Largeintra","Bonferroni.Tra") %>%
  tidyr::gather(.,Bonferroni,value,Bonferroni.1vs0234:Bonferroni.Tra,factor_key = T)
plot.function(all_sample_signature,errorbar.df,variable.name)
ggsave(file.path(outpath,"simple.egfr.kras.percentage.ttest.Bonferroni.pdf"),width = 6,height = 24)  #pay attention to line 429-445, choose ttest or ranktest
ggsave(file.path(outpath,"simple.egfr.kras.percentage.ranktest.Bonferroni.pdf"),width = 6,height = 24)  #pay attention to line 429-445, choose ttest or ranktest
##################### Plot3:triangle #################################################
triangle.df.fdr <- fisher_result_tested %>% select("Code1","Code2","Group1","Group2",
                                               "FDR.Del1","FDR.Del2","FDR.Del3","FDR.TD1","FDR.TD2",
                                               "FDR.Fbinv","FDR.Largeintra","FDR.Tra") %>%
  tidyr::gather(.,compare,FDR,FDR.Del1:FDR.Tra,factor_key = T) %>% 
  mutate(FDR.simplify=ifelse(FDR<0.001,3,
                       ifelse(FDR<0.01,2,
                              ifelse(FDR<0.1,1,0)))) %>%
  mutate(FDR.simplify=factor(FDR.simplify)) 
triangle.df.t <-fisher_result_tested %>% select("Code1","Code2","Group1","Group2",
                                                    "t.Del1","t.Del2","t.Del3","t.TD1","t.TD2",
                                                    "t.Fbinv","t.Largeintra","t.Tra") %>%
  tidyr::gather(.,t,tvalue,t.Del1:t.Tra,factor_key = T)
triangle.df <- cbind(triangle.df.fdr,triangle.df.t[,c("t","tvalue")])

  
triangle.df$diff <- apply(triangle.df, 1, function(pair) grep(TRUE,strsplit(pair[1], "")[[1]] != strsplit(pair[2], "")[[1]]))
triangle.df$comparision <- ifelse(triangle.df$diff==3 & (grepl("S_U",triangle.df$Group1)|grepl("S_U",triangle.df$Group2)),"Smoking",
                                  ifelse(triangle.df$diff==3 & (grepl("N_A",triangle.df$Group1)|grepl("N_A",triangle.df$Group2)),"Population",
                                         ifelse(triangle.df$diff==4,"Gender",
                                                ifelse(triangle.df$diff==5,"TP53",
                                                       ifelse(triangle.df$diff==6 & (grepl("EGFRmut",triangle.df$Group1)|grepl("EGFRmut",triangle.df$Group2)),"EGFR","KRAS")))
                                         )
                                  )
triangle.df <- triangle.df %>% mutate(Code1=factor(Code1,levels = variable.name),
                                      Code2=factor(Code2,levels = rev(variable.name)))


ggplot(data=triangle.df)+
  geom_point(aes(x=Code1,y=Code2,color=FDR.simplify,shape=comparision),size=2) +
  facet_grid(compare~.,scales ="fixed")+
  scale_color_manual(breaks = c(0,1,2,3),
                    values = c("grey80","#F8C471","#DC7633","red"),
                    labels = c("NS","<0.1","<0.001","<0.000001"),drop=FALSE)+
  scale_shape_manual(values = c(15,16,17,7,8,13),drop=FALSE)+
  scale_y_discrete(
                   breaks = variable.name,
                   labels = full.name,
                   drop=FALSE
                   ) +
  scale_x_discrete(
                   breaks = variable.name,
                   labels = full.name,
                   position = "top",
                   drop=FALSE
                   ) +
  coord_cartesian(clip = "off") +
  theme(axis.text.x.top = element_text(size=6,colour = "black",angle=90,hjust=0.5,vjust=0.3),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x =  element_blank(),
        axis.ticks.y =  element_blank(),
        axis.ticks.length.y =unit(0.8, "mm"),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        panel.grid.major = element_line(size=0.1189,color="grey90"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size=6,colour = "black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        # strip.position="top",
        legend.key = element_blank(),
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(4,"mm"),
        legend.title = element_text(size=6,colour = "black"),
        legend.text = element_text(size=6),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        plot.margin=grid::unit(c(2,2,0,2), "mm"))
ggsave(file.path(outpath,"simple.egfr.kras.ttest.pdf"),width = 4,height = 18)


##################### Plot4 :freq #################################################
variable.name <- c("ek2111","ek2112","ek1223","ek1222",
                   "ek1212","ek2212","ek1122","ek3111")
full.name <- c("N_U-F-TP53wt-EGFRmutKRASwt","N_U-F-TP53wt-EGFRwtKRASwt","S_U-M-TP53mut-EGFRwtKRASmut","S_U-M-TP53mut-EGFRwtKRASwt",
               "S_U-M-TP53wt-EGFRwtKRASwt","N_U-M-TP53wt-EGFRwtKRASwt","S_U-F-TP53mut-EGFRwtKRASwt","N_A-F-TP53wt-EGFRmutKRASwt")
all_sample_signature=rbind(ek2111,ek2112,ek1223,ek1222,
                           ek1212,ek2212,ek1122,ek3111)

variable.name <- c("ek1111","ek1112","ek1113","ek1121","ek1122","ek1123","ek1211","ek1212","ek1213","ek1221","ek1222","ek1223",
                   "ek2111","ek2112","ek2113","ek2121","ek2122","ek2123","ek2211","ek2212","ek2213","ek2221","ek2222","ek2223",
                   "ek3111","ek3112","ek3113","ek3121","ek3122","ek3123","ek3211","ek3212","ek3213","ek3221","ek3222","ek3223")
full.name <- c("S_U-F-TP53wt-EGFRmutKRASwt","S_U-F-TP53wt-EGFRwtKRASwt","S_U-F-TP53wt-EGFRwtKRASmut",
               "S_U-F-TP53mut-EGFRmutKRASwt","S_U-F-TP53mut-EGFRwtKRASwt","S_U-F-TP53mut-EGFRwtKRASmut",
               "S_U-M-TP53wt-EGFRmutKRASwt","S_U-M-TP53wt-EGFRwtKRASwt","S_U-M-TP53wt-EGFRwtKRASmut",
               "S_U-M-TP53mut-EGFRmutKRASwt","S_U-M-TP53mut-EGFRwtKRASwt","S_U-M-TP53mut-EGFRwtKRASmut",
               
               "N_U-F-TP53wt-EGFRmutKRASwt","N_U-F-TP53wt-EGFRwtKRASwt","N_U-F-TP53wt-EGFRwtKRASmut",
               "N_U-F-TP53mut-EGFRmutKRASwt","N_U-F-TP53mut-EGFRwtKRASwt","N_U-F-TP53mut-EGFRwtKRASmut",
               "N_U-M-TP53wt-EGFRmutKRASwt","N_U-M-TP53wt-EGFRwtKRASwt","N_U-M-TP53wt-EGFRwtKRASmut",
               "N_U-M-TP53mut-EGFRmutKRASwt","N_U-M-TP53mut-EGFRwtKRASwt","N_U-M-TP53mut-EGFRwtKRASmut",
               
               "N_A-F-TP53wt-EGFRmutKRASwt","N_A-F-TP53wt-EGFRwtKRASwt","N_A-F-TP53wt-EGFRwtKRASmut",
               "N_A-F-TP53mut-EGFRmutKRASwt","N_A-F-TP53mut-EGFRwtKRASwt","N_A-F-TP53mut-EGFRwtKRASmut",
               "N_A-M-TP53wt-EGFRmutKRASwt","N_A-M-TP53wt-EGFRwtKRASwt","N_A-M-TP53wt-EGFRwtKRASmut",
               "N_A-M-TP53mut-EGFRmutKRASwt","N_A-M-TP53mut-EGFRwtKRASwt","N_A-M-TP53mut-EGFRwtKRASmut")
plot.groups <- c("N_U-F-TP53wt-EGFRmutKRASwt","N_U-F-TP53wt-EGFRwtKRASwt",
                 "S_U-M-TP53mut-EGFRwtKRASmut","S_U-M-TP53mut-EGFRwtKRASwt",
                 "S_U-M-TP53wt-EGFRwtKRASwt","N_U-M-TP53wt-EGFRwtKRASwt",
                 "S_U-F-TP53mut-EGFRwtKRASwt","N_A-F-TP53wt-EGFRmutKRASwt") # for main figure

mysig=c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra")
plot.groups <- full.name # for supple figure
names(full.name) <- variable.name
all_sample_signature=rbind(ek1111,ek1112,ek1113,ek1121,ek1122,ek1123,ek1211,ek1212,ek1213,ek1221,ek1222,ek1223,
                           ek2111,ek2112,ek2113,ek2121,ek2122,ek2123,ek2211,ek2212,ek2213,ek2221,ek2222,ek2223,
                           ek3111,ek3112,ek3113,ek3121,ek3122,ek3123,ek3211,ek3212,ek3213,ek3221,ek3222,NULL)
freq.df <- merge(all_sample_input,
                 all_sample_signature %>% 
                   select("sample","category","mechanism","n") %>% 
                   tidyr::spread(.,mechanism,n),
                 all.x = T,
                 by.x = "Tumor_Barcode",by.y = "sample",
                 ) %>% filter(SP_Group_new %in% c("N_A","S_U","N_U")) %>%
  mutate(category=paste0(SP_Group_new,"-",substring(Gender,1,1),"-","TP53",tp53,"-","EGFR",egfr,"KRAS",kras)) %>%
  filter(category %in% plot.groups) %>%
  dplyr::select("Tumor_Barcode","category",mysig) %>%
  replace(is.na(.), 0) %>%
  tidyr::gather(.,Sig,Num,Del1:Tra) %>%
  mutate(category=factor(category,levels = plot.groups))
ggplot(freq.df) + 
  geom_boxplot(aes(x=category,y=Num),size=0.1,outlier.size = 0.1
               # outlier.shape = NA
               ) +
  facet_wrap(factor(Sig,levels=mysig)~.,scales = "free",nrow=3)+
  scale_x_discrete(drop=FALSE)+
  scale_y_log10()+
  ylab("Number of signature") + 
  theme_bw()+
  theme(axis.text.x = element_text(size=6,colour = "black",angle = 90,vjust=0.5),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.1,colour = "black"), 
        axis.line.x.top = element_line(size=0.1,colour = "black"), 
        axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
        axis.ticks.length.x =unit(0.8, "mm"),
        axis.ticks.length.y =unit(0.8, "mm"),
        panel.background = element_blank(),
        # panel.margin.x=unit(0.25, "cm"),
        panel.spacing.x = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size=6, angle = 90, colour = "black",hjust = 1,vjust = 0.5),
        strip.background =  element_blank(),
        strip.placement = "top",
        plot.margin=grid::unit(c(0,2,0,2), "mm"),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.15, "cm"),
        legend.text = element_text(size=6,colour = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.margin=margin(t=-20),
        legend.box.margin=margin(0,0,0,0))
ggsave(file.path(outpath,"simple.egfr.kras.freq.8groups.pdf"),width = 4,height = 8)
ggsave(file.path(outpath,"simple.egfr.kras.freq.36groups.pdf"),width = 7,height = 10)
