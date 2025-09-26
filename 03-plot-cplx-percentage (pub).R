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
library(ggpubr)


###function
p.adjust_function <- function(df,col=3,p.bfrn,p.fdr){
  df[,p.bfrn] <- p.adjust(df[,col],method = "bonferroni",n=nrow(df))
  df[,p.fdr] <- p.adjust(df[,col],method = "BH")
  names(df)[c(p.bfrn,p.fdr)] <- c("Bonferroni","FDR")
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
  
  unique_sample_with_signature=unique(cluster[,c("sample", "mechanism")])
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
  all_sample_signature$num.cgr=sample_with_signature_count
  all_sample_signature$num.nocgr=all_sample_count-sample_with_signature_count
  return(all_sample_signature)
}

plot.function <- function(all_sample_signature,fisher_result_tested,variable.name){
  category.label <- unique(all_sample_signature$category_label)
  names(category.label) <- unique(all_sample_signature$category)
  p <- ggplot()+
    geom_bar(data=all_sample_signature,
             aes(x=percent,y=sample,fill=mechanism),
             stat = "identity",position = "stack",width = 1.2)+
    geom_text(data=all_sample_signature %>% dplyr::select(category,num.cgr,num.nocgr) %>% unique(),
              aes(x=0.5,y=Inf,label=paste0(num.cgr,"/",num.nocgr)),size=2)+
    facet_wrap(category~.,
               labeller = labeller(category=category.label),
               scale="free",drop=TRUE,strip.position = "bottom",nrow=1)+
    scale_fill_manual(values = c("1 ecDNA/double minutes"="#e78ac3",
                                 "2 BFB cycles/chromatin bridge"="#fc8d62", 
                                 "3 Large loss"="#a6d854", 
                                 "4 Micronuclei"="#66c2a5",
                                 "5 Large gain"="#ffd92f",
                                 "6 Hourglass"="#8da0cb",
                                 "chromoplexy"="#5DADE2",
                                 "cycle_templated_ins"="#F8C471",
                                 "complex_unclear"="#C89EC7",
                                 "Non-CGR"="lightgrey",
                                 "Unassigned"="black"))+
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
          strip.text.x = element_text(angle=90,size = 6,color = "black"), 
          strip.background = element_blank())
  if (is.null(fisher_result_tested)==FALSE){
    code.nchar <- nchar(variable.name[1])
    errorbar.df <- fisher_result_tested
    errorbar.df$Code1.xposition <- apply(errorbar.df, 1,function(x) which(variable.name == x["Code1"]))
    errorbar.df$Code2.xposition <- apply(errorbar.df, 1,function(x) which(variable.name == x["Code2"]))
    errorbar.df$yposition <- apply(errorbar.df, 1, function(x) code.nchar+1-grep(TRUE,strsplit(x[1], "")[[1]] != strsplit(x[2], "")[[1]]))
    errorbar.df$color=ifelse(errorbar.df$FDR<0.05,"red","black")
    p.error <- ggplot(errorbar.df) +
      geom_errorbar(aes(xmin=Code1.xposition+0.1,xmax=Code2.xposition-0.1,y=yposition),size=0.1189)+
      geom_text(aes(x=(Code1.xposition+Code2.xposition)/2,y=yposition,label=signif(FDR,2),color=color),size=2,vjust=-1)+
      scale_x_continuous(labels=variable.name,breaks=seq(1,length(variable.name),1),limits = c(1,length(variable.name)))+
      scale_color_manual(breaks = c("black","red"),values = c("black","red"))+
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
    p <- ggpubr::ggarrange(p.error,p,nrow = 2)

  }
  return(p)
}
### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
hq.path <- file.path(sherlock.path,"/Data/SAMPLE_INFO/HQ_samples.csv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
cluster.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_CGR_signature.tsv")
noncluster.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_nonclustered_cmplx_sv.tsv")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","plot")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

### read
tp53_mutated_samples=read.delim(tp53.mut.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
kras_mutated_samples=read.delim(kras.mut.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F)  %>% pull(V1)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
sbs4 <- read.delim2(sbs4.path)
histology <- read.csv(histology.path)
# all sample
all_sample <- read.delim(sampleinfo.path, header=T, sep="\t", stringsAsFactors = F) %>% 
  dplyr::select(Tumor_Barcode, Assigned_Population, Gender, SP_Group,Smoking,Histology,high_quality)
all_sample <- merge(all_sample,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

# cluster complex SV
cluster=read.delim(cluster.path, header=T, sep="\t", stringsAsFactors = F)
colnames(cluster)=c("cgr_chr", "chr","start", "end", "cgr_status", "link_chr", "mechanism", "sample")
cluster$mechanism <- ifelse(cluster$mechanism=="","Unassigned",cluster$mechanism)
# noncluster complex sv
noncluster <- read.delim(noncluster.path) %>% dplyr::select(SAMPLE,FULL_SV_TYPE) %>% 
  setnames(.,old=c("SAMPLE","FULL_SV_TYPE"),
           new=c("sample","mechanism")) %>% 
  mutate(mechanism=ifelse(mechanism %in% c("chromoplexy_del","chromoplexy_bal"),
                          "chromoplexy",mechanism)) 
# mutation
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
# add sbs4
all_sample <- merge(all_sample,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
all_sample <- all_sample %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
######################################
#### clustered cplx sv signatures ####
######################################
### choose cluster
mycluster <- cluster
mycluster <- noncluster
#### choose high quality
all_sample_input=all_sample 
all_sample_input=all_sample %>% filter(high_quality ==1)
all_sample_input=all_sample %>% filter(Histology %in% c("Adenocarcinoma"))
all_sample_input=all_sample %>% filter(Histology %in% c("Adenocarcinoma")) %>% filter(high_quality ==1)


# 1. All samples
all_sample_data=generate_data(all_samples=all_sample_input,cluster=mycluster,category_name="All Samples")
# 2. EUR smoker
eur_smoking_data=generate_data(all_samples=all_sample_input %>% filter(SP_Group_new=="S_U"),cluster=mycluster,category_name="EUR Smoker")
# 3. EUR non-smoker
eur_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(SP_Group_new=="N_U"),cluster=mycluster,category_name="EUR Never Smoker")
# 4. EAS non-smoker
eas_nonsmokint_data=generate_data(all_samples=all_sample_input %>% filter(SP_Group_new=="N_A"),cluster=mycluster,category_name="EAS Never Smoker")
# 5. Female
female_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Female"),cluster=mycluster,category_name="Females")
female_luad_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Female" & Histology=="Adenocarcinoma"),cluster=mycluster,category_name="Female_luad")
female_luad_nonsmoker_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Female" & Histology=="Adenocarcinoma"& Smoking == "Non-Smoker"),cluster=mycluster,category_name="Female_luad_Never Smoker")
female_nonsmoker_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Female" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="Female_Never Smoker")
# 6. Males
male_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Male"),cluster=mycluster,category_name="Males")
male_luad_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Male" & Histology=="Adenocarcinoma"),cluster=mycluster,category_name="Male_luad")
male_luad_nonsmoker_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Male" & Histology=="Adenocarcinoma"& Smoking == "Non-Smoker"),cluster=mycluster,category_name="Male_luad_Never Smoker")
male_nonsmoker_data=generate_data(all_samples=all_sample_input %>% filter(Gender=="Male" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="Male_Never Smoker")
# 7.1. TP53
tp53wt_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "mut"),cluster=mycluster,category_name="TP53 WT")
tp53mut_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "wt"),cluster=mycluster,category_name="TP53 Mut")
# 7.2. TP53 smoking and nonsmoking
tp53wt_smoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "wt" & Smoking == "Smoker"),cluster=mycluster,category_name="TP53 WT Smoker")
tp53mut_smoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "mut" & Smoking == "Smoker"),cluster=mycluster,category_name="TP53 Mut Smoker")
tp53wt_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "wt" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="TP53 WT Never Smoker")
tp53mut_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(tp53 == "mut" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="TP53 Mut Never Smoker")
# 8.1. EGFR
egfr.L591R_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="L591R"),cluster=mycluster,category_name="EGFR L591R")
egfr.E479_A483del_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="A483del"),cluster=mycluster,category_name="EGFR Exon 19 del")
egfr.Others_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="Others"),cluster=mycluster,category_name="EGFR Others")
egfr.wt_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="wt"),cluster=mycluster,category_name="EGFR WT")
# 8.2. EGFR smoking and nonsmoking
egfr.L591R_smoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="L591R" & Smoking == "Smoker"),cluster=mycluster,category_name="EGFR L591R Smoker")
egfr.E479_smoking_A483del_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="A483del" & Smoking == "Smoker"),cluster=mycluster,category_name="EGFR Exon 19 del Smoker")
egfr.Others_smoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="Others" & Smoking == "Smoker"),cluster=mycluster,category_name="EGFR Others Smoker")
egfr.wt_smoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="wt" & Smoking == "Smoker"),cluster=mycluster,category_name="EGFR WT Smoker")
egfr.L591R_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="L591R" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="EGFR L591R Never Smoker")
egfr.E479_nonsmoking_A483del_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="A483del" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="EGFR Exon 19 del Never Smoker")
egfr.Others_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="Others" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="EGFR Others Never Smoker")
egfr.wt_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(egfr.detail=="wt" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="EGFR WT Never Smoker")
# 9.1. KRAS
kraswt_data=generate_data(all_samples=all_sample_input %>% filter(kras=="wt"),cluster=mycluster,category_name="KRAS WT")
krasmut_data=generate_data(all_samples=all_sample_input %>% filter(kras=="mut"),cluster=mycluster,category_name="KRAS Mut")
# 9.2. KRAS smoking and nonsmoking
kraswt_smoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="wt" & Smoking == "Smoker"),cluster=mycluster,category_name="KRAS WT Smoker")
krasmut_smoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="mut" & Smoking == "Smoker"),cluster=mycluster,category_name="KRAS Mut Smoker")
kraswt_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="wt" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="KRAS WT Never Smoker")
krasmut_nonsmoking_data=generate_data(all_samples=all_sample_input %>% filter(kras=="mut" & Smoking == "Non-Smoker"),cluster=mycluster,category_name="KRAS Mut Never Smoker")
# 10. EGFR by SP_Group_new, gender, tp53, egfr
e1 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr=="wt"),mycluster,"S_U-F-TP53wt-EGFRwt") #33
e2 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr=="wt"),mycluster,"S_U-M-TP53wt-EGFRwt") #103
e3 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr=="wt"),mycluster,"S_U-F-TP53mut-EGFRwt") #65
e4 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr=="wt"),mycluster,"S_U-M-TP53mut-EGFRwt") #165
e5 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr=="mut"),mycluster,"S_U-F-TP53wt-EGFRmut") #2
e6 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr=="mut"),mycluster,"S_U-M-TP53wt-EGFRmut") #5
e7 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr=="mut"),mycluster,"S_U-F-TP53mut-EGFRmut") #2
e8 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr=="mut"),mycluster,"S_U-M-TP53mut-EGFRmut") #13
e9 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr=="wt"),mycluster,"N_U-F-TP53wt-EGFRwt") #229
e10 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr=="wt"),mycluster,"N_U-M-TP53wt-EGFRwt") #55
e11 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr=="wt"),mycluster,"N_U-F-TP53mut-EGFRwt") #53
e12 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr=="wt"),mycluster,"N_U-M-TP53mut-EGFRwt") #13
e13 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr=="mut"),mycluster,"N_U-F-TP53wt-EGFRmut") #133
e14 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr=="mut"),mycluster,"N_U-M-TP53wt-EGFRmut") #23
e15 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr=="mut"),mycluster,"N_U-F-TP53mut-EGFRmut") #63
e16 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr=="mut"),mycluster,"N_U-M-TP53mut-EGFRmut") #12
e17 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr=="wt"),mycluster,"N_A-F-TP53wt-EGFRwt") #88
e18 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr=="wt"),mycluster,"N_A-M-TP53wt-EGFRwt") #23
e19 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr=="wt"),mycluster,"N_A-F-TP53mut-EGFRwt") #32
e20 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr=="wt"),mycluster,"N_A-M-TP53mut-EGFRwt") #6
e21 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr=="mut"),mycluster,"N_A-F-TP53wt-EGFRmut") #167
e22 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr=="mut"),mycluster,"N_A-M-TP53wt-EGFRmut") #27
e23 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr=="mut"),mycluster,"N_A-F-TP53mut-EGFRmut") #120
e24 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr=="mut"),mycluster,"N_A-M-TP53mut-EGFRmut") #18
# 11. KRAS by SP_Group_new, Gender, tp53, egfr
k1 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & kras=="wt"),mycluster,"S_U-F-TP53wt-KRASwt") #22
k2 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & kras=="wt"),mycluster,"S_U-M-TP53wt-KRASwt") #56
k3 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & kras=="wt"),mycluster,"S_U-F-TP53mut-KRASwt") #48
k4 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & kras=="wt"),mycluster,"S_U-M-TP53mut-KRASwt") #147
k5 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & kras=="mut"),mycluster,"S_U-F-TP53wt-KRASmut") #13
k6 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & kras=="mut"),mycluster,"S_U-M-TP53wt-KRASmut") #52
k7 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & kras=="mut"),mycluster,"S_U-F-TP53mut-KRASmut") #19
k8 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & kras=="mut"),mycluster,"S_U-M-TP53mut-KRASmut") #31
k9 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & kras=="wt"),mycluster,"N_U-F-TP53wt-KRASwt") #334
k10 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & kras=="wt"),mycluster,"N_U-M-TP53wt-KRASwt") #73
k11 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & kras=="wt"),mycluster,"N_U-F-TP53mut-KRASwt") #110
k12 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & kras=="wt"),mycluster,"N_U-M-TP53mut-KRASwt") #24
k13 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & kras=="mut"),mycluster,"N_U-F-TP53wt-KRASmut") #28
k14 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & kras=="mut"),mycluster,"N_U-M-TP53wt-KRASmut") #5
k15 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & kras=="mut"),mycluster,"N_U-F-TP53mut-KRASmut") #6
k16 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & kras=="mut"),mycluster,"N_U-M-TP53mut-KRASmut") #1
k17 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & kras=="wt"),mycluster,"N_A-F-TP53wt-KRASwt") #249
k18 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & kras=="wt"),mycluster,"N_A-M-TP53wt-KRASwt") #49
k19 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & kras=="wt"),mycluster,"N_A-F-TP53mut-KRASwt") #24
k20 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & kras=="wt"),mycluster,"N_A-M-TP53mut-KRASwt") #6
k21 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & kras=="mut"),mycluster,"N_A-F-TP53wt-KRASmut") #6
k22 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & kras=="mut"),mycluster,"N_A-M-TP53wt-KRASmut") #6
k23 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & kras=="mut"),mycluster,"N_A-F-TP53mut-KRASmut") #6
k24 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & kras=="mut"),mycluster,"N_A-M-TP53mut-KRASmut") #
# 12. EGFR-KRAS by SP_Group_new, Gender, tp53
ek1111 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"S_U-F-TP53wt-EGFRmutKRASwt") #0
ek1112 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"S_U-F-TP53wt-EGFRwtKRASwt") #22
ek1113 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"S_U-F-TP53wt-EGFRwtKRASmut") #11
ek1121 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"S_U-F-TP53mut-EGFRmutKRASwt") #1
ek1122 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"S_U-F-TP53mut-EGFRwtKRASwt") #47
ek1123 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"S_U-F-TP53mut-EGFRwtKRASmut") #18
ek1211 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"S_U-M-TP53wt-EGFRmutKRASwt") #4
ek1212 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"S_U-M-TP53wt-EGFRwtKRASwt") #52
ek1213 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"S_U-M-TP53wt-EGFRwtKRASmut") #51
ek1221 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"S_U-M-TP53mut-EGFRmutKRASwt") #13
ek1222 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"S_U-M-TP53mut-EGFRwtKRASwt") #134
ek1223 <- generate_data(all_sample_input %>% filter(SP_Group_new=="S_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"S_U-M-TP53mut-EGFRwtKRASmut") #31

ek2111 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"N_U-F-TP53wt-EGFRmutKRASwt") #133
ek2112 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"N_U-F-TP53wt-EGFRwtKRASwt") #201
ek2113 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"N_U-F-TP53wt-EGFRwtKRASmut") #28
ek2121 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"N_U-F-TP53mut-EGFRmutKRASwt") #63
ek2122 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"N_U-F-TP53mut-EGFRwtKRASwt") #47
ek2123 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"N_U-F-TP53mut-EGFRwtKRASmut") #6
ek2211 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"N_U-M-TP53wt-EGFRmutKRASwt") #23
ek2212 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"N_U-M-TP53wt-EGFRwtKRASwt") #50
ek2213 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"N_U-M-TP53wt-EGFRwtKRASmut") #5
ek2221 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"N_U-M-TP53mut-EGFRmutKRASwt") #12
ek2222 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"N_U-M-TP53mut-EGFRwtKRASwt") #12
ek2223 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_U" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"N_U-M-TP53mut-EGFRwtKRASmut") #1

ek3111 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"N_A-F-TP53wt-EGFRmutKRASwt") #167
ek3112 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"N_A-F-TP53wt-EGFRwtKRASwt") #82
ek3113 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"N_A-F-TP53wt-EGFRwtKRASmut") #6
ek3121 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"N_A-F-TP53mut-EGFRmutKRASwt") #120
ek3122 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"N_A-F-TP53mut-EGFRwtKRASwt") #31
ek3123 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Female" & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"N_A-F-TP53mut-EGFRwtKRASmut") #1
ek3211 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"N_A-M-TP53wt-EGFRmutKRASwt") #27
ek3212 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"N_A-M-TP53wt-EGFRwtKRASwt") #22
ek3213 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"N_A-M-TP53wt-EGFRwtKRASmut") #1
ek3221 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"N_A-M-TP53mut-EGFRmutKRASwt") #18
ek3222 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"N_A-M-TP53mut-EGFRwtKRASwt") #6
ek3223 <- generate_data(all_sample_input %>% filter(SP_Group_new=="N_A" & Gender=="Male" & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"N_A-M-TP53mut-EGFRwtKRASmut") #0

# 13. only EGFR/KRAS/TP53 AND smoking
ns.tw.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"NeverSmoker-TP53wt-EGFRwtKRASwt") #370
ns.tw.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"NeverSmoker-TP53wt-EGFRmutKRASwt") #354
ns.tw.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"NeverSmoker-TP53wt-EGFRwtKRASmut") #44
ns.tm.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"NeverSmoker-TP53mut-EGFRwtKRASwt") #105
ns.tm.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"NeverSmoker-TP53mut-EGFRmutKRASwt") #217
ns.tm.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4") & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"NeverSmoker-TP53mut-EGFRwtKRASmut") #10
ss.tw.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "wt" & egfr.kras=="wt-wt"),mycluster,"Smoker-TP53wt-EGFRwtKRASwt") #77
ss.tw.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "wt" & egfr.kras=="mut-wt"),mycluster,"Smoker-TP53wt-EGFRmutKRASwt") #4
ss.tw.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "wt" & egfr.kras=="wt-mut"),mycluster,"Smoker-TP53wt-EGFRwtKRASmut") #64
ss.tm.ew.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "mut" & egfr.kras=="wt-wt"),mycluster,"Smoker-TP53mut-EGFRwtKRASwt") #186
ss.tm.em.kw <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "mut" & egfr.kras=="mut-wt"),mycluster,"Smoker-TP53mut-EGFRmutKRASwt") #16
ss.tm.ew.km <- generate_data(all_sample_input %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4") & tp53 == "mut" & egfr.kras=="wt-mut"),mycluster,"Smoker-TP53mut-EGFRwtKRASmut") #49

nrow(ns.tw.ew.kw)+nrow(ns.tw.em.kw)+nrow(ns.tw.ew.km)+
  nrow(ns.tm.ew.kw)+nrow(ns.tm.em.kw)+nrow(ns.tm.ew.km)+
  nrow(ss.tw.ew.kw)+nrow(ss.tw.em.kw)+nrow(ss.tw.ew.km)+
  nrow(ss.tm.ew.kw)+nrow(ss.tm.em.kw)+nrow(ss.tm.ew.km)
#############################  plot0 only EGFR/KRAS/TP53 AND smoking  ###################################
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
# order sample
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = full.name)) %>% 
  dplyr::arrange(.,category,desc(percentage),mech_all) 
all_sample_signature=all_sample_signature %>% 
  mutate(mech_all=factor(mech_all,
                         levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR"))) %>% 
  arrange(category,desc(percentage),mech_all)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
# fisher test
all_pairs <- combn(variable.name, 2, simplify = TRUE) # Use combn to generate all possible pairs
test_pairs <- all_pairs[, apply(all_pairs, 2, function(pair) sum(strsplit(pair[1], "")[[1]] != strsplit(pair[2], "")[[1]]) == 1)] # Filter pairs with one character difference
test_pairs <- matrix(test_pairs, ncol = 2, byrow = TRUE) #t
fisher_result <- t(apply(test_pairs, 1, function(x){
  v1 <- full.name[x[1]]
  v2 <- full.name[x[2]]
  v1.cgr <- all_sample_signature %>% filter(category==v1) %>% dplyr::select(num.cgr) %>% unique() %>% as.numeric()
  v1.nocgr <- all_sample_signature %>% filter(category==v1) %>% dplyr::select(num.nocgr) %>% unique() %>% as.numeric()
  v2.cgr <- all_sample_signature %>% filter(category==v2) %>% dplyr::select(num.cgr) %>% unique() %>% as.numeric()
  v2.nocgr <- all_sample_signature %>% filter(category==v2) %>% dplyr::select(num.nocgr) %>% unique() %>% as.numeric()
  v1.total <- sum(c(v1.cgr,v1.nocgr),na.rm = T)
  v2.total <- sum(c(v2.cgr,v2.nocgr),na.rm = T)
  if (v1.total >= 20 & v2.total >= 20){
    fisher <- fisher.test(matrix(c(v1.cgr,v1.nocgr,v2.cgr,v2.nocgr),nrow = 2,ncol = 2))
    out <- c(x[1],x[2],v1,v2,v1.total,v2.total,v1.cgr,v1.nocgr,v2.cgr,v2.nocgr,fisher$p.value,fisher$estimate)
  } else {
    out <- c(x[1],x[2],v1,v2,v1.total,v2.total,v1.cgr,v1.nocgr,v2.cgr,v2.nocgr,"Not tested","Not tested")
  }
  return(out)
})) %>% as.data.frame()#fisher's test
names(fisher_result) <- c("Code1","Code2","Group1","Group2","Group1.num","Group2.num","Group1.cgr.num","Group1.nocgr.num","Group2.cgr.num","Group2.nocgr.num","P","OR")
# filter fisher, only keep tested
fisher_result_tested <- fisher_result %>% filter(P !="Not tested") %>% mutate(P=as.numeric(P),OR=as.numeric(OR))
# filter fisher, remove pairs EGFRmutKRASwt and EGFRwtKRASmut
fisher_result_tested <- fisher_result_tested %>% 
  filter((grepl("EGFRmutKRASwt",Group1)&grepl("EGFRwtKRASmut",Group2))==FALSE)
# adjust p value
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=11,13,14)

# plot
plot.function(all_sample_signature,fisher_result_tested,variable.name)
ggsave(file.path(outpath,"clustered.smk.tp53.egfr.kras.percentage.luad.pdf"),width = 6,height = 7)
ggsave(file.path(outpath,"nonclustered.smk.tp53.egfr.kras.percentage.luad.pdf"),width = 6,height = 7)
ggsave(file.path(outpath,"nonclustered.smk.tp53.egfr.kras.percentage.luad.hq.pdf"),width = 6,height = 7)


#############################  plot1 others  ###################################
all_sample_signature=rbind(all_sample_data, 
                           eur_smoking_data,eur_nonsmoking_data,eas_nonsmokint_data,
                           female_data,male_data,
                           tp53wt_data,tp53mut_data,
                           egfr.L591R_data,egfr.E479_A483del_data,egfr.Others_data,egfr.wt_data,
                           kraswt_data,krasmut_data,
                           tp53wt_smoking_data,tp53mut_smoking_data,tp53wt_nonsmoking_data,tp53mut_nonsmoking_data,
                           egfr.wt_smoking_data,
                           #egfr.L591R_smoking_data,
                           egfr.E479_smoking_A483del_data,egfr.Others_smoking_data,
                           egfr.L591R_nonsmoking_data,egfr.E479_nonsmoking_A483del_data,egfr.Others_nonsmoking_data,egfr.wt_nonsmoking_data,
                           kraswt_smoking_data,krasmut_smoking_data,kraswt_nonsmoking_data,krasmut_nonsmoking_data
                           )
# order sample
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = c("All Samples",
                                             "EAS Never Smoker","EUR Never Smoker","EUR Smoker",
                                             "Females","Males",
                                             "TP53 WT","TP53 Mut",
                                             "EGFR WT","EGFR L591R","EGFR Exon 19 del","EGFR Others",
                                             "KRAS WT","KRAS Mut",
                                             "TP53 WT Smoker","TP53 Mut Smoker","TP53 WT Never Smoker","TP53 Mut Never Smoker",
                                             "EGFR WT Smoker",
                                             #"EGFR L591R Smoker",
                                             "EGFR Exon 19 del Smoker","EGFR Others Smoker",
                                             "EGFR WT Never Smoker","EGFR L591R Never Smoker","EGFR Exon 19 del Never Smoker","EGFR Others Never Smoker",
                                             "KRAS WT Smoker","KRAS Mut Smoker","KRAS WT Never Smoker","KRAS Mut Never Smoker"
                                             ))) %>% 
  arrange(category,desc(percentage),mech_all) 
all_sample_signature=all_sample_signature %>% mutate(mech_all=factor(mech_all,levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR"))) %>% arrange(category,desc(percentage),mech_all)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
# plot
plot.function(all_sample_signature,NULL,NULL)
ggsave(file.path(outpath,"nonclustered.cplx.percentage.pdf"),width = 8,height = 4)
ggsave(file.path(outpath,"clustered.cplx.percentage.pdf"),width = 8,height = 4)
#ggsave(file.path(outpath,"clustered.cplx.percentage.hq.luad.pdf"),width = 8,height = 4)
#ggsave(file.path(outpath,"clustered.cplx.percentage.luad.pdf"),width = 8,height = 4)
#ggsave(file.path(outpath,"nonclustered.cplx.percentage.hq.luad.pdf"),width = 8,height = 4)
#ggsave(file.path(outpath,"nonclustered.cplx.percentage.luad.pdf"),width = 8,height = 4)
####################################################################

#############################  plot2  egfr ###################################
all_sample_signature=rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19,e20,e21,e22,e23,e24)
# order sample
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = c("S_U-F-TP53wt-EGFRwt","S_U-F-TP53wt-EGFRmut","S_U-F-TP53mut-EGFRwt","S_U-F-TP53mut-EGFRmut",
                                             "S_U-M-TP53wt-EGFRwt","S_U-M-TP53wt-EGFRmut","S_U-M-TP53mut-EGFRwt","S_U-M-TP53mut-EGFRmut",
                                             "N_U-F-TP53wt-EGFRwt","N_U-F-TP53wt-EGFRmut","N_U-F-TP53mut-EGFRwt","N_U-F-TP53mut-EGFRmut",
                                             "N_U-M-TP53wt-EGFRwt","N_U-M-TP53wt-EGFRmut","N_U-M-TP53mut-EGFRwt","N_U-M-TP53mut-EGFRmut",
                                             "N_A-F-TP53wt-EGFRwt","N_A-F-TP53wt-EGFRmut","N_A-F-TP53mut-EGFRwt","N_A-F-TP53mut-EGFRmut",
                                             "N_A-M-TP53wt-EGFRwt","N_A-M-TP53wt-EGFRmut","N_A-M-TP53mut-EGFRwt","N_A-M-TP53mut-EGFRmut"
                                             ))) %>% 
  arrange(category,desc(percentage),mech_all) 
all_sample_signature=all_sample_signature %>% mutate(mech_all=factor(mech_all,levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR"))) %>% arrange(category,desc(percentage),mech_all)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
# plot
plot.function(all_sample_signature,NULL,NULL)
ggsave(file.path(outpath,"clustered.egfr.percentage.pdf"),width = 5,height = 4)
ggsave(file.path(outpath,"nonclustered.egfr.percentage.pdf"),width = 5,height = 4)
####################################################################

#############################  plot3  kras ###################################
all_sample_signature=rbind(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20,k21,k22,k23)
# order sample
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = c("S_U-F-TP53wt-KRASwt","S_U-F-TP53wt-KRASmut","S_U-F-TP53mut-KRASwt","S_U-F-TP53mut-KRASmut",
                                             "S_U-M-TP53wt-KRASwt","S_U-M-TP53wt-KRASmut","S_U-M-TP53mut-KRASwt","S_U-M-TP53mut-KRASmut",
                                             "N_U-F-TP53wt-KRASwt","N_U-F-TP53wt-KRASmut","N_U-F-TP53mut-KRASwt","N_U-F-TP53mut-KRASmut",
                                             "N_U-M-TP53wt-KRASwt","N_U-M-TP53wt-KRASmut","N_U-M-TP53mut-KRASwt","N_U-M-TP53mut-KRASmut",
                                             "N_A-F-TP53wt-KRASwt","N_A-F-TP53wt-KRASmut","N_A-F-TP53mut-KRASwt","N_A-F-TP53mut-KRASmut",
                                             "N_A-M-TP53wt-KRASwt","N_A-M-TP53wt-KRASmut","N_A-M-TP53mut-KRASwt","N_A-M-TP53mut-KRASmut"
  ))) %>% 
  arrange(category,desc(percentage),mech_all) 
all_sample_signature=all_sample_signature %>% mutate(mech_all=factor(mech_all,levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR"))) %>% arrange(category,desc(percentage),mech_all)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
# plot
plot.function(all_sample_signature,NULL,NULL)
ggsave(file.path(outpath,"clustered.kras.percentage.pdf"),width = 5,height = 4)
ggsave(file.path(outpath,"nonclustered.kras.percentage.pdf"),width = 5,height = 4)
####################################################################

#############################  plot3  egfr-kras ###################################
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
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = full.name)) %>% 
  dplyr::arrange(.,category,desc(percentage),mech_all) 
all_sample_signature=all_sample_signature %>% 
  mutate(mech_all=factor(mech_all,
                         levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR"))) %>% 
  arrange(category,desc(percentage),mech_all)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
# fisher test
all_pairs <- combn(variable.name, 2, simplify = TRUE) # Use combn to generate all possible pairs
test_pairs <- all_pairs[, apply(all_pairs, 2, function(pair) sum(strsplit(pair[1], "")[[1]] != strsplit(pair[2], "")[[1]]) == 1)] # Filter pairs with one character difference
test_pairs <- matrix(test_pairs, ncol = 2, byrow = TRUE) #t
fisher_result <- t(apply(test_pairs, 1, function(x){
  v1 <- full.name[x[1]]
  v2 <- full.name[x[2]]
  v1.cgr <- all_sample_signature %>% filter(category==v1) %>% dplyr::select(num.cgr) %>% unique() %>% as.numeric()
  v1.nocgr <- all_sample_signature %>% filter(category==v1) %>% dplyr::select(num.nocgr) %>% unique() %>% as.numeric()
  v2.cgr <- all_sample_signature %>% filter(category==v2) %>% dplyr::select(num.cgr) %>% unique() %>% as.numeric()
  v2.nocgr <- all_sample_signature %>% filter(category==v2) %>% dplyr::select(num.nocgr) %>% unique() %>% as.numeric()
  v1.total <- sum(c(v1.cgr,v1.nocgr),na.rm = T)
  v2.total <- sum(c(v2.cgr,v2.nocgr),na.rm = T)
  if (v1.total >= 10 & v2.total >= 10){
    fisher <- fisher.test(matrix(c(v1.cgr,v1.nocgr,v2.cgr,v2.nocgr),nrow = 2,ncol = 2))
    out <- c(x[1],x[2],v1,v2,v1.total,v2.total,v1.cgr,v1.nocgr,v2.cgr,v2.nocgr,fisher$p.value,fisher$estimate)
  } else {
    out <- c(x[1],x[2],v1,v2,v1.total,v2.total,v1.cgr,v1.nocgr,v2.cgr,v2.nocgr,"Not tested","Not tested")
  }
  return(out)
})) %>% as.data.frame()#fisher's test
names(fisher_result) <- c("Code1","Code2","Group1","Group2","Group1.num","Group2.num","Group1.cgr.num","Group1.nocgr.num","Group2.cgr.num","Group2.nocgr.num","P","OR")
# filter fisher, only keep tested
fisher_result_tested <- fisher_result %>% filter(P !="Not tested") %>% mutate(P=as.numeric(P),OR=as.numeric(OR))
# filter fisher, remove pairs EGFRmutKRASwt and EGFRwtKRASmut
fisher_result_tested <- fisher_result_tested %>% 
  filter((grepl("EGFRmutKRASwt",Group1)&grepl("EGFRwtKRASmut",Group2))==FALSE) %>%
  filter((grepl("S_U",Group1)&grepl("N_A",Group2))==FALSE)
# adjust p value
fisher_result_tested <- p.adjust_function(fisher_result_tested,col=11,13,14)

# plot
plot.function(all_sample_signature,fisher_result_tested,variable.name)
ggsave(file.path(outpath,"clustered.egfr.kras.percentage.pdf"),width = 6,height = 7)
ggsave(file.path(outpath,"nonclustered.egfr.kras.percentage.pdf"),width = 6,height = 7)
####################################################################

#############################  plot4 sex in nonsmoker/luad  ###################################
all_sample_signature=rbind(female_data,male_data,
                           female_nonsmoker_data,male_nonsmoker_data,
                           female_luad_data,male_luad_data,
                           female_luad_nonsmoker_data,male_luad_nonsmoker_data)
# order sample
all_sample_signature=all_sample_signature %>% 
  mutate(category=factor(category,levels = c("Females","Males",
                                             "Female_Never Smoker","Male_Never Smoker",
                                             "Female_luad","Male_luad",
                                             "Female_luad_Never Smoker","Male_luad_Never Smoker"
  ))) %>% 
  arrange(category,desc(percentage),mech_all) 
all_sample_signature=all_sample_signature %>% mutate(mech_all=factor(mech_all,levels = c(all_sample_signature$mech_all[all_sample_signature$mech_all != "Non-CGR"] %>% unique(),"Non-CGR"))) %>% arrange(category,desc(percentage),mech_all)
all_sample_signature$sample=factor(all_sample_signature$sample,levels=unique(all_sample_signature$sample))
# plot
plot.function(all_sample_signature,NULL,NULL)
ggsave(file.path(outpath,"clustered.cplx.percentage.sex.pdf"),width = 3,height = 4)
ggsave(file.path(outpath,"clustered.cplx.percentage.sex.hq.pdf"),width = 3,height = 4)
ggsave(file.path(outpath,"nonclustered.cplx.percentage.sex.hq.pdf"),width = 3,height = 4)
ggsave(file.path(outpath,"nonclustered.cplx.percentage.sex.pdf"),width = 3,height = 4)

fisher.test(matrix(c(514,156,109,39),nrow = 2))
fisher.test(matrix(c(551,130,271,70),nrow = 2))
fisher.test(matrix(c(480,110,95,20),nrow = 2))
fisher.test(matrix(c(563,107,117,31),nrow = 2))
fisher.test(matrix(c(603,78,310,31),nrow = 2))
fisher.test(matrix(c(520,70,99,16),nrow = 2))

fisher.test(matrix(c(285,98,63,30),nrow = 2))
fisher.test(matrix(c(294,68,147,32),nrow = 2))
fisher.test(matrix(c(259,60,52,13),nrow = 2))
fisher.test(matrix(c(318,65,71,22),nrow = 2))
fisher.test(matrix(c(329,33,167,12),nrow = 2))
fisher.test(matrix(c(288,31,56,9),nrow = 2))

fisher.test(matrix(c(384,37,312,44),nrow = 2))
fisher.test(matrix(c(312,44,229,16),nrow = 2))
fisher.test(matrix(c(603,78,310,31),nrow = 2))

####################################################################

#############################  plot 5 Fig 2a bars ###################################
merge_sample_info_all <- readRDS(file.path(scratch.path,"04-shatterseek","1217_samples_hg38","plot","merge_sample_info_all.rds"))
sample.order <- merge_sample_info_all$Tumor_Barcode

## 1. choose clustered
mycluster <- cluster
# generate data
all_sample_data=generate_data(all_samples=all_sample,
                              cluster=mycluster,
                              category_name="Fig2a.clustered")
# order
all_sample_data=all_sample_data %>% 
  mutate(mech_all=factor(mech_all,levels = c(all_sample_data$mech_all[all_sample_data$mech_all != "Non-CGR"] %>% 
                                               unique(),"Non-CGR"))) %>% 
  mutate(sample=factor(sample,levels=sample.order)) %>%
  arrange(sample,desc(percentage),mech_all)
all_sample_data$sample=factor(all_sample_data$sample,levels=unique(all_sample_data$sample))
# plot
clustered.plot <- plot.function(all_sample_data,NULL,NULL)

## 2. choose non clustered
mycluster <- noncluster
# generate data
all_sample_data=generate_data(all_samples=all_sample,
                              cluster=mycluster,
                              category_name="Fig2a.nonclustered")
# order
all_sample_data=all_sample_data %>% 
  mutate(mech_all=factor(mech_all,levels = c(all_sample_data$mech_all[all_sample_data$mech_all != "Non-CGR"] %>% 
                                               unique(),"Non-CGR"))) %>% 
  mutate(sample=factor(sample,levels=sample.order)) %>%
  arrange(sample,desc(percentage),mech_all)
all_sample_data$sample=factor(all_sample_data$sample,levels=unique(all_sample_data$sample))
# plot
nonclustered.plot <- plot.function(all_sample_data,NULL,NULL)

p <- ggarrange(clustered.plot,nonclustered.plot,nrow=1)
ggsave(plot=p,
       file.path(outpath,"clustered.nonclustered.bars.pdf"),width =3,height = 4)
