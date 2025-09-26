### modifeied from Xiaoming's script
library(dplyr)
library(sjPlot)
library(ggplot2)
library(MASS)
library(jtools)
library(data.table)

### function
myfit_df_function <- function(fitresult){
  fit=fitresult
  star_p <- c(0.05,0.01,0.001)
  myfit_df <- as.data.frame(exp(cbind(OR = coef(fit), confint(fit))))
  myfit_df$variable <- rownames(myfit_df)
  p_df  <-  data.frame(coef(summary(fit)))
  p_df$variable <- rownames(p_df)
  names(p_df) <- c("estimate","Stderror","tvalue","pvalue","variable")
  myfit_df <- merge(myfit_df,p_df[,c("variable","pvalue")],all.x = T)
  myfit_df$star <- apply(myfit_df,1,function(x) paste0(rep("*",sum(as.numeric(x["pvalue"]) < star_p,na.rm = T)),collapse = ""))
  # sort
  myfit_df$variable <- factor(myfit_df$variable,
                              levels = rev(c("Assigned_PopulationEUR",
                                             "SmokingSmoker",
                                             "SmokingUnknown",
                                             "GenderMale",
                                             "HistologyAdenosquamous carcinoma",
                                             "HistologyCarcinoid tumor",
                                             "HistologySquamous cell carcinoma",
                                             "HistologyOthers",
                                             "purity",
                                             "wgdWGD",
                                             "scna_groupMezzo-forte",
                                             "scna_groupForte",
                                             "StageII","StageIII","StageIV",
                                             c(paste0("factor(",mygene,")","Mut")),
                                             "fusionMutated",
                                             "driver.freq",
                                             "(Intercept)")))
  return(myfit_df)
}

analysis_data_function <- function(){
  analysis_data <- merge(other_sample_infomation,snv[,c("sample",mygene)],by.x = "Tumor_Barcode",by.y="sample",all.x = T) %>% 
    dplyr::select(Tumor_Barcode, Assigned_Population, Smoking, Gender, Histology, purity, wgd, scna_group, Stage, mut_load,mygene) %>%
    mutate(
      ecDNA=ifelse(Tumor_Barcode %in% ecDNA_samples,1,0),
      bfb=ifelse(Tumor_Barcode %in% bfb_samples,1,0),
      large_loss=ifelse(Tumor_Barcode %in% large_loss_samples,1,0),
      micronuclei=ifelse(Tumor_Barcode %in% micronuclei_samples,1,0),
      large_gain=ifelse(Tumor_Barcode %in% large_gain_samples,1,0),
      hourglass=ifelse(Tumor_Barcode %in% hourglass_samples,1,0),
      any=ifelse(Tumor_Barcode %in% no_cgr_samples,0,1),
      
      chromoplexy=ifelse(Tumor_Barcode %in% chromoplexy_samples,1,0),
      cycle_templated_ins=ifelse(Tumor_Barcode %in% cycle_templated_ins_samples,1,0),
      complex_unclear=ifelse(Tumor_Barcode %in% complex_unclear_samples,1,0),
      any_nonclustered=ifelse(Tumor_Barcode %in% no_nonclustered_samples,0,1),
      
      fusion=ifelse(Tumor_Barcode %in% fusion_samples,"Mutated","Wildtype"),
      Stage=ifelse(Stage %in% c("I","IA","IA1","IA2","IA3","IB"),"I",
                   ifelse(Stage %in% c("II","IIA","IIB"),"II",
                          ifelse(Stage %in% c("III","IIIA","IIIB"),"III",
                                 ifelse(Stage %in% c("IV","IVA"),"IV",NA))))
    ) %>% 
    filter(Assigned_Population %in% c("EUR","EAS")) %>% # only in EUR and EAS
    mutate(across(all_of(mygene), ~ ifelse(. == "None", "Wild", "Mut"))) %>%
    mutate(
      across(all_of(mygene), ~ factor(.,levels = c("Wild","Mut"))),
      fusion=factor(fusion, levels=c("Wildtype", "Mutated")),
      scna_group=factor(scna_group, levels = c("Piano", "Mezzo-forte", "Forte")),
      histology=factor(Histology, levels = c("Adenocarcinoma","Squamous cell carcinoma","Carcinoid tumor","Adenosquamous carcinoma","Others")),
      population=factor(Assigned_Population,levels = c("EAS","EUR")),
      Stage=factor(Stage,levels=c("I","II","III","IV"))
    )
  # add driver SV freq
  oncoprint_table$driver.freq <- rowSums(
    !is.na(oncoprint_table[, c(loss_gistic, gain_gistic)]) &
      oncoprint_table[, c(loss_gistic, gain_gistic)] != "NA"
  )
  analysis_data <- merge(analysis_data,oncoprint_table[,c("Tumor_Barcode","driver.freq")],all.x=T)
  # add Smoking_SBS4
  analysis_data <- merge(analysis_data,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
  analysis_data <- analysis_data %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
  return(analysis_data)
}

clustered_function <- function(){
  myfit_df_all <- rbind(myfit_df_function(any_regression_result) %>% mutate(sig="Any Complex SVs"),
                        myfit_df_function(ecDNA_regression_result) %>% mutate(sig="1 ecDNA"),
                        myfit_df_function(bfb_regression_result) %>% mutate(sig="2 Chr bridge"),
                        myfit_df_function(large_loss_regression_result) %>% mutate(sig="3 Large loss"),
                        myfit_df_function(micronuclei_regression_result) %>% mutate(sig="4 Micronulei"),
                        myfit_df_function(large_gain_regression_result) %>% mutate(sig="5 Large gain"),
                        myfit_df_function(hourglass_regression_result) %>% mutate(sig="6 Hourglass"))
  
  myfit_df_all$sig <- gsub("\\."," ",myfit_df_all$sig)
  myfit_df_all$sig <- factor(myfit_df_all$sig,levels = c("Any Complex SVs","1 ecDNA","2 Chr bridge","3 Large loss","4 Micronulei","5 Large gain","6 Hourglass"))
  myfit_df_all$ORstar <- ifelse(grepl("\\*",myfit_df_all$star),paste0(signif(myfit_df_all$OR,2)," ",myfit_df_all$star) ,NA)  
  myfit_df_all$OR_color <- ifelse(myfit_df_all$OR>1,"red","blue")
  myfit_df_all$OR_color <- ifelse(grepl("\\*",myfit_df_all$star),myfit_df_all$OR_color ,"black")
  myfit_df_all <- myfit_df_all[myfit_df_all$variable %in% "(Intercept)" ==FALSE,]
  myvariable <- unique(myfit_df_all$variable)
  variable_fullname <- c("European VS Asian",
                         "Male VS Female",
                         "Smoker VS Never smoker",
                         "Unknown VS Never smoker",
                         "Adenosquamous carcinoma",
                         "Carcinoid tumor",
                         "Squamous cell carcinoma",
                         "Others",
                         "Tumor purity",
                         "Whole-genome duplication",
                         "Mezzo-forte",
                         "Forte",
                         "StageII","StageIII","StageIV",
                         paste0(mygene," mutant VS Wildtype"),
                         "Fusion VS Wildtype",
                         "Driver SV num.")
  names(variable_fullname) <- c("Assigned_PopulationEUR",
                                "GenderMale",
                                "SmokingSmoker",
                                "SmokingUnknown",
                                "HistologyAdenosquamous carcinoma",
                                "HistologyCarcinoid tumor",
                                "HistologySquamous cell carcinoma",
                                "HistologyOthers",
                                "purity",
                                "wgdWGD",
                                "scna_groupMezzo-forte",
                                "scna_groupForte",
                                "StageII","StageIII","StageIV",
                                paste0("factor(",mygene,")Mut"),
                                "fusionMutated",
                                "driver.freq")
  myfit_df_all$variablefullname <- variable_fullname[as.character(myfit_df_all$variable)]
  xmax <- ifelse(max(myfit_df_all[,4],na.rm=T)>5,5,ceiling(max(myfit_df_all[,4],na.rm=T)))
  print(myfit_df_all)
  
  myfit_df_all_plot <- myfit_df_all %>% filter(variable != "SmokingUnknown") %>% mutate(
    variable=factor(variable,
                    rev(c("Assigned_PopulationEUR",
                          "SmokingSmoker",
                          "GenderMale",
                          "HistologySquamous cell carcinoma",
                          "HistologyCarcinoid tumor",
                          "HistologyAdenosquamous carcinoma",
                          "HistologyOthers",
                          "purity",
                          "wgdWGD",
                          "scna_groupMezzo-forte",
                          "scna_groupForte",
                          "StageII","StageIII","StageIV",
                          paste0("factor(",mygene,")Mut"),
                          "fusionMutated",
                          "driver.freq"))
    )
  )
  p <- ggplot(data=myfit_df_all_plot, 
              aes(y=variable, x=OR, xmin=`2.5 %`, xmax=`97.5 %`)) +
    geom_point(size=0.00001) +
    expand_limits(x = 5) + 
    scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                  labels=c(0.01,0.01,0.1,1,10,100),
                  limits=c(0.01,100),
                  oob = scales::oob_keep
    ) +
    # scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000),
    #               labels=c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000),
    #               limits=c(0.0001,10000)
    # ) +
    # scale_x_log10()+
    scale_y_discrete(breaks=unique(myfit_df_all_plot$variable),
                     labels=unique(myfit_df_all_plot$variablefullname))+
    facet_grid(.~sig,scales = "free") +
    scale_size_area(max_size = 1)+
    coord_cartesian(clip = "on") +
    geom_vline(xintercept = 1, linetype="solid",color="black",lwd=0.2) + 
    geom_errorbarh(height=0,size=0.1) + 
    geom_text(data=myfit_df_all_plot, aes(color=OR_color,x=ifelse(OR>1,`2.5 %`,`97.5 %`),y=variable,label=ORstar,hjust = ifelse(OR>1,0,1)),vjust = -0.6,size=1.8)+
    scale_color_manual(values = c("red","blue","black"),
                       breaks = c("red","blue","black")) +
    xlab("Odds Ratio") + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=6,colour = "black"),
          axis.text.x = element_text(size=6,colour = "black",angle=90,hjust = 1),
          axis.text.y = element_text(size=6,colour = "black"),
          strip.text.x = element_text(size=6,colour = "black"),
          axis.ticks = element_line(color = "black",size=0.1),
          axis.ticks.length=unit(0.02,"inch"),
          legend.position = "none",
          plot.title =element_text(size=6, face='bold',hjust = 0,vjust = 0),
          panel.background =element_blank(),
          panel.border=element_blank(), 
          axis.line=element_line(color = "black",size=0.1))
  return(list(p,myfit_df_all_plot))
}

nonclustered_funciton <- function(){
  myfit_df_all <- rbind(myfit_df_function(any_nonclustered_regression_result) %>% mutate(sig="Any non-clustered cplx SVs"),
                        myfit_df_function(chromoplexy_regression_result) %>% mutate(sig="Chromoplexy"),
                        myfit_df_function(cycle_templated_ins_regression_result) %>% mutate(sig="Cycle of templated insertion"),
                        myfit_df_function(complex_unclear_regression_result) %>% mutate(sig="Complex unclear"))
  
  myfit_df_all$sig <- gsub("\\."," ",myfit_df_all$sig)
  myfit_df_all$sig <- factor(myfit_df_all$sig,levels = c("Any non-clustered cplx SVs","Chromoplexy","Cycle of templated insertion","Complex unclear"))
  myfit_df_all$ORstar <- ifelse(grepl("\\*",myfit_df_all$star),paste0(signif(myfit_df_all$OR,2)," ",myfit_df_all$star) ,NA)  
  myfit_df_all$OR_color <- ifelse(myfit_df_all$OR>1,"red","blue")
  myfit_df_all$OR_color <- ifelse(grepl("\\*",myfit_df_all$star),myfit_df_all$OR_color ,"black")
  myfit_df_all <- myfit_df_all[myfit_df_all$variable %in% "(Intercept)" ==FALSE,]
  myvariable <- unique(myfit_df_all$variable)
  variable_fullname <- c("European VS Asian",
                         "Male VS Female",
                         "Smoker VS Never smoker",
                         "Unknown VS Never smoker",
                         "Adenosquamous carcinoma",
                         "Carcinoid tumor",
                         "Squamous cell carcinoma",
                         "Others",
                         "Tumor purity",
                         "Whole-genome duplication",
                         "Mezzo-forte",
                         "Forte",
                         "StageII","StageIII","StageIV",
                         paste0(mygene," mutant VS Wildtype"),
                         "Fusion VS Wildtype",
                         "Driver SV num.")
  names(variable_fullname) <- c("Assigned_PopulationEUR",
                                "GenderMale",
                                "SmokingSmoker",
                                "SmokingUnknown",
                                "HistologyAdenosquamous carcinoma",
                                "HistologyCarcinoid tumor",
                                "HistologySquamous cell carcinoma",
                                "HistologyOthers",
                                "purity",
                                "wgdWGD",
                                "scna_groupMezzo-forte",
                                "scna_groupForte",
                                "StageII","StageIII","StageIV",
                                paste0("factor(",mygene,")Mut"),
                                "fusionMutated",
                                "driver.freq")
  myfit_df_all$variablefullname <- variable_fullname[as.character(myfit_df_all$variable)]
  xmax <- ifelse(max(myfit_df_all[,4],na.rm=T)>5,5,ceiling(max(myfit_df_all[,4],na.rm=T)))
  print(myfit_df_all)
  
  myfit_df_all_plot <- myfit_df_all %>% filter(variable != "SmokingUnknown") %>% mutate(
    variable=factor(variable,
                    rev(c("Assigned_PopulationEUR",
                          "SmokingSmoker",
                          "GenderMale",
                          "HistologySquamous cell carcinoma",
                          "HistologyCarcinoid tumor",
                          "HistologyAdenosquamous carcinoma",
                          "HistologyOthers",
                          "purity",
                          "wgdWGD",
                          "scna_groupMezzo-forte",
                          "scna_groupForte",
                          "StageII","StageIII","StageIV",
                          paste0("factor(",mygene,")Mut"),
                          "fusionMutated",
                          "driver.freq"))
    )
  )
  p <- ggplot(data=myfit_df_all_plot, 
              aes(y=variable, x=OR, xmin=`2.5 %`, xmax=`97.5 %`)) +
    geom_point(size=0.00001) +
    expand_limits(x = 5) + 
    scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                  labels=c(0.01,0.01,0.1,1,10,100),
                  limits=c(0.01,100),
                  oob = scales::oob_keep) +
    
    # scale_x_log10(breaks=c(1E-4,1E-2,1,1E2,1E4,1E6),
    #               labels=c(1E-4,1E-2,1,1E2,1E4,1E6),
    #               limits=c(1E-4,1e6),
    #               oob = scales::oob_keep) +

    # scale_x_log10(breaks=c(1E-6, 1E-4,1E-2,1,1E2,1E4,1E6,1E8),
    #               labels=c(1E-6, 1E-4,1E-2,1,1E2,1E4,1E6,1E8),
    #               limits=c(1E-7, 1E9),
    #               oob = scales::oob_keep) +
    
    # scale_x_log10(breaks=c(1E-2, 1,1E2,1E4),
    #               labels=c(1E-2, 1,1E2,1E4),
    #               limits=c(1E-2, 1000),
    # oob = scales::oob_keep) +
    scale_y_discrete(breaks=unique(myfit_df_all_plot$variable),
                     labels=unique(myfit_df_all_plot$variablefullname))+
    facet_grid(.~sig,scales = "free") +
    scale_size_area(max_size = 1)+
    coord_cartesian(clip = "on") +
    geom_vline(xintercept = 1, linetype="solid",color="black",lwd=0.2) + 
    geom_errorbarh(height=0,size=0.1) + 
    geom_text(data=myfit_df_all_plot, aes(color=OR_color,x=ifelse(OR>1,`2.5 %`,`97.5 %`),y=variable,label=ORstar,hjust = ifelse(OR>1,0,1)),vjust = -0.6,size=1.8)+
    scale_color_manual(values = c("red","blue","black"),
                       breaks = c("red","blue","black")) +
    xlab("Odds Ratio") + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=6,colour = "black"),
          axis.text.x = element_text(size=6,colour = "black",angle=90,hjust = 1),
          axis.text.y = element_text(size=6,colour = "black"),
          strip.text.x = element_text(size=6,colour = "black"),
          axis.ticks = element_line(color = "black",size=0.1),
          axis.ticks.length=unit(0.02,"inch"),
          legend.position = "none",
          plot.title =element_text(size=6, face='bold',hjust = 0,vjust = 0),
          panel.background =element_blank(),
          panel.border=element_blank(), 
          axis.line=element_line(color = "black",size=0.1))
  return(list(p,myfit_df_all_plot))
}

### path 
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
signature.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38","sherlock_CGR_signature.tsv")
hq.path <- file.path(sherlock.path,"/Data/SAMPLE_INFO/HQ_samples.csv")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
egfr.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/regression/EGFR_mutated_samples.tsv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
fusion_path  <- file.path(sherlock.path,"scratch/fusion/end_annotation/tier1_driver_fusions.txt")
noncluster.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_nonclustered_cmplx_sv.tsv")
snv_path <- file.path(sherlock.path,"scratch/06-SNV/sample_snv.txt")
sv.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.suppl.txt")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
outpath <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","plot")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

cnv_set <- c("CDKN2A","TERT","MDM2","EGFR","PTPRD",
             "MYC","CTNND2","NKX2-1","PTK6","MCL1",
             "FGFR1",
             # "GRIN2A",
             "CCNE1","STK11","MET","BRAF")
loss_set <- c("CDKN2A","PTPRD","STK11")
gain_set <- cnv_set[cnv_set %in% loss_set==FALSE]

cnv_values <- paste0("cnv.",cnv_set)
gistic_values <- paste0("gistic.",cnv_set)
loss_values <- paste0(loss_set,".loss")
gain_values <- paste0(gain_set,".gain")
loss_gistic <- paste0("gistic.",loss_set)
gain_gistic <- paste0("gistic.",gain_set)
### read
histology <- read.csv(histology.path)
sv <- read.delim(sv.path)
snv  <-  read.delim(snv_path,row.names = 1)
snv_copy <- snv
cluster=read.delim(signature.path, header=T, sep="\t", stringsAsFactors = F)
colnames(cluster)=c("cgr_chr", "chr","start", "end", "cgr_status", "link_chr", "mechanism", "sample")
other_sample_infomation=read.delim(sampleinfo.path, header=T, check.names=F, stringsAsFactors=F, na.strings="no_na")
fusion <- read.delim2(fusion_path)
fusion_samples=fusion %>% filter(Known_Fusion=="Yes") %>% pull(SAMPLE) %>% unique()
noncluster <- read.delim(noncluster.path) %>% dplyr::select(SAMPLE,FULL_SV_TYPE) %>% 
  setnames(.,old=c("SAMPLE","FULL_SV_TYPE"),
           new=c("sample","mechanism")) %>% 
  mutate(mechanism=ifelse(mechanism %in% c("chromoplexy_del","chromoplexy_bal"),
                          "chromoplexy",mechanism)) 
hqsample <- read.csv(hq.path) %>% pull(Tumor_Barcode)
snv$sample <- row.names(snv)
oncoprint_table <- read.delim(oncoprint_table_path,check.names = F) %>%
  mutate(across(all_of(cnv_values), ~ ifelse((.x >= 5 & baseline_cn <= 2.7)|(.x >= 9 & baseline_cn > 2.7), "gain",
                                             ifelse((.x == 0 & baseline_cn <= 2.7)|(.x < (baseline_cn-2.7) & baseline_cn > 2.7),
                                                    "loss",
                                                    NA)))) %>%
  mutate(across(all_of(gain_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==3,"gain",NA)))) %>%
  mutate(across(all_of(loss_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==-3,"loss",NA))))
sbs4 <- read.delim2(sbs4.path)
other_sample_infomation <- merge(other_sample_infomation,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

### add information
ecDNA_samples=cluster %>% filter(mechanism=="1 ecDNA/double minutes") %>% pull(sample)  %>% unique()
bfb_samples=cluster %>% filter(mechanism=="2 BFB cycles/chromatin bridge") %>% pull(sample) %>% unique()
large_loss_samples=cluster %>% filter(mechanism=="3 Large loss") %>% pull(sample)  %>% unique()
micronuclei_samples=cluster %>% filter(mechanism=="4 Micronuclei") %>% pull(sample) %>% unique()
large_gain_samples=cluster %>% filter(mechanism=="5 Large gain") %>% pull(sample) %>% unique()
hourglass_samples=cluster %>% filter(mechanism=="6 Hourglass") %>% pull(sample) %>% unique()
no_cgr_samples=other_sample_infomation %>% filter(Tumor_Barcode %in% cluster$sample==FALSE) %>% pull(Tumor_Barcode)
chromoplexy_samples=noncluster %>% filter(mechanism=="chromoplexy") %>% pull(sample) %>% unique()
cycle_templated_ins_samples=noncluster %>% filter(mechanism=="cycle_templated_ins") %>% pull(sample) %>% unique()
complex_unclear_samples=noncluster %>% filter(mechanism=="complex_unclear") %>% pull(sample) %>% unique()
no_nonclustered_samples=other_sample_infomation %>% filter(Tumor_Barcode %in% noncluster$sample==FALSE) %>% pull(Tumor_Barcode)

### calculate mutation rate
# Create a copy of the original data (optional)
# Replace "None" with 0 (no mutation) and any other value with 1 (mutation)
snv_binary <- as.data.frame(apply(snv_copy, 2, function(x) ifelse(x == "None", 0, 1)))
# Convert the result to numeric for proper summation
snv_binary <- as.data.frame(lapply(snv_binary, as.numeric))
# Calculate the mutation rate for each gene
mutation_rate <- colSums(snv_binary) / nrow(snv_binary)
# Convert to data.frame to view results clearly
mutation_rate_df <- data.frame(Gene = names(mutation_rate), MutationRate = mutation_rate)
# View the mutation rate for each gene
print(mutation_rate_df)
### genes that mut rate is higher than 5% 10%
priority_genes <- c("TP53", "EGFR", "KRAS")
# genes.5pct <- mutation_rate_df %>% filter(MutationRate >= 0.05) %>% pull(Gene) #144
# genes.5pct <- c(priority_genes[priority_genes %in% genes.5pct], genes.5pct[!genes.5pct %in% priority_genes])
# genes.10pct <- mutation_rate_df %>% filter(MutationRate >= 0.1) %>% pull(Gene) #17
# genes.10pct <- c(priority_genes[priority_genes %in% genes.10pct], genes.10pct[!genes.10pct %in% priority_genes])
mygene <- c("TP53","EGFR","KRAS")
# mygene <- genes.5pct
# mygene <- genes.10pct



################################################################
## logistic regression for clustered complex SV ###
analysis_data <- analysis_data_function() %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4"))
analysis_data <- analysis_data %>% filter(Tumor_Barcode %in% hqsample) #for high quality
ecDNA_data=analysis_data[analysis_data$Tumor_Barcode %in% c(ecDNA_samples, no_cgr_samples),]
myformula <- as.formula(paste0("ecDNA","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
ecDNA_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=ecDNA_data)
bfb_data=analysis_data[analysis_data$Tumor_Barcode %in% c(bfb_samples, no_cgr_samples),]
myformula <- as.formula(paste0("bfb","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
bfb_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=bfb_data)
large_loss_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_loss_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_loss","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_loss_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_loss_data)
micronuclei_data=analysis_data[analysis_data$Tumor_Barcode %in% c(micronuclei_samples, no_cgr_samples),]
myformula <- as.formula(paste0("micronuclei","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
micronuclei_regression_result <- glm(formula=myformula,  family=binomial(link="logit"), data=micronuclei_data)
large_gain_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_gain_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_gain","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_gain_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_gain_data)
hourglass_data=analysis_data[analysis_data$Tumor_Barcode %in% c(hourglass_samples, no_cgr_samples),]
myformula <- as.formula(paste0("hourglass","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
hourglass_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=hourglass_data)
myformula <- as.formula(paste0("any","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=analysis_data)
clustered_function()[1]
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression_hq.pdf"),height = 2.7,width = 5)
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression.pdf"),height = 2.7,width = 5)
write.csv(clustered_function()[2],file.path(outpath,"pcawg_signature_multivariate_logit_regression_hq.csv"))
write.csv(clustered_function()[2],file.path(outpath,"pcawg_signature_multivariate_logit_regression.csv"))
################################################################
### logistic regression for non-clustered complex SV ##
analysis_data <- analysis_data_function() %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4"))
analysis_data <- analysis_data %>% filter(Tumor_Barcode %in% hqsample) #for high quality
chromoplexy_data=analysis_data[analysis_data$Tumor_Barcode %in% c(chromoplexy_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("chromoplexy","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
chromoplexy_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=chromoplexy_data)
cycle_templated_ins_data=analysis_data[analysis_data$Tumor_Barcode %in% c(cycle_templated_ins_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("cycle_templated_ins","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
cycle_templated_ins_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=cycle_templated_ins_data)
complex_unclear_data=analysis_data[analysis_data$Tumor_Barcode %in% c(complex_unclear_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("complex_unclear","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
complex_unclear_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=complex_unclear_data)
any_nonclustered_data <- analysis_data
myformula <- as.formula(paste0("any_nonclustered","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_nonclustered_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=any_nonclustered_data)
nonclustered_funciton()[1]
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression_hq.pdf"),height = 2.7,width = 3.5)
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression.pdf"),height = 2.7,width = 3.5)
write.csv(nonclustered_funciton()[2],file.path(outpath,"nonclustered_cplx_multivariate_logit_regression_hq.csv"))
write.csv(nonclustered_funciton()[2],file.path(outpath,"nonclustered_cplx_multivariate_logit_regression.csv"))
##################################################################################################


##################################################################################################
### breakdown into smoking ###
analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Smoker-SBS4") #for smoking only
ecDNA_data=analysis_data[analysis_data$Tumor_Barcode %in% c(ecDNA_samples, no_cgr_samples),]
myformula <- as.formula(paste0("ecDNA","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
ecDNA_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=ecDNA_data)
bfb_data=analysis_data[analysis_data$Tumor_Barcode %in% c(bfb_samples, no_cgr_samples),]
myformula <- as.formula(paste0("bfb","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
bfb_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=bfb_data)
large_loss_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_loss_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_loss","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
large_loss_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_loss_data)
micronuclei_data=analysis_data[analysis_data$Tumor_Barcode %in% c(micronuclei_samples, no_cgr_samples),]
myformula <- as.formula(paste0("micronuclei","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
micronuclei_regression_result <- glm(formula=myformula,  family=binomial(link="logit"), data=micronuclei_data)
large_gain_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_gain_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_gain","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
large_gain_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_gain_data)
hourglass_data=analysis_data[analysis_data$Tumor_Barcode %in% c(hourglass_samples, no_cgr_samples),]
myformula <- as.formula(paste0("hourglass","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
hourglass_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=hourglass_data)
myformula <- as.formula(paste0("any","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
any_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=analysis_data)
clustered_function()[[1]]+
  scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                labels=c(0.01,0.01,0.1,1,10,100),
                limits=c(0.01,100),
                oob = scales::oob_keep
  )
  
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression.Smoker.pdf"),height = 2,width = 5)

analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Smoker-SBS4") #for smoking only
chromoplexy_data=analysis_data[analysis_data$Tumor_Barcode %in% c(chromoplexy_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("chromoplexy","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
chromoplexy_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=chromoplexy_data)
cycle_templated_ins_data=analysis_data[analysis_data$Tumor_Barcode %in% c(cycle_templated_ins_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("cycle_templated_ins","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
cycle_templated_ins_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=cycle_templated_ins_data)
complex_unclear_data=analysis_data[analysis_data$Tumor_Barcode %in% c(complex_unclear_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("complex_unclear","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
complex_unclear_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=complex_unclear_data)
any_nonclustered_data <- analysis_data
myformula <- as.formula(paste0("any_nonclustered","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
any_nonclustered_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=any_nonclustered_data)
nonclustered_funciton()[[1]]+
  scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100,1000,10000),
                labels=c(0.01,0.01,0.1,1,10,100,100,100),
                limits=c(0.01,10000),
                oob = scales::oob_keep
  )
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression.Smoker.pdf"),height = 2,width = 3.5)


##################################################################################################
### breakdown into non-smoking ###
analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4")
ecDNA_data=analysis_data[analysis_data$Tumor_Barcode %in% c(ecDNA_samples, no_cgr_samples),]
myformula <- as.formula(paste0("ecDNA","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
ecDNA_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=ecDNA_data)
bfb_data=analysis_data[analysis_data$Tumor_Barcode %in% c(bfb_samples, no_cgr_samples),]
myformula <- as.formula(paste0("bfb","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
bfb_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=bfb_data)
large_loss_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_loss_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_loss","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_loss_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_loss_data)
micronuclei_data=analysis_data[analysis_data$Tumor_Barcode %in% c(micronuclei_samples, no_cgr_samples),]
myformula <- as.formula(paste0("micronuclei","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
micronuclei_regression_result <- glm(formula=myformula,  family=binomial(link="logit"), data=micronuclei_data)
large_gain_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_gain_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_gain","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_gain_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_gain_data)
hourglass_data=analysis_data[analysis_data$Tumor_Barcode %in% c(hourglass_samples, no_cgr_samples),]
myformula <- as.formula(paste0("hourglass","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
hourglass_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=hourglass_data)
myformula <- as.formula(paste0("any","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=analysis_data)
clustered_function()[[1]]+
  scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                labels=c(0.01,0.01,0.1,1,10,100),
                limits=c(0.01,100),
                oob = scales::oob_keep
  )
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression.NeverSmoker.pdf"),height = 2.3,width = 5)

analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4") #for smoking only
chromoplexy_data=analysis_data[analysis_data$Tumor_Barcode %in% c(chromoplexy_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("chromoplexy","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
chromoplexy_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=chromoplexy_data)
cycle_templated_ins_data=analysis_data[analysis_data$Tumor_Barcode %in% c(cycle_templated_ins_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("cycle_templated_ins","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
cycle_templated_ins_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=cycle_templated_ins_data)
complex_unclear_data=analysis_data[analysis_data$Tumor_Barcode %in% c(complex_unclear_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("complex_unclear","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
complex_unclear_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=complex_unclear_data)
any_nonclustered_data <- analysis_data
myformula <- as.formula(paste0("any_nonclustered","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_nonclustered_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=any_nonclustered_data)
nonclustered_funciton()[[1]]+
  scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                labels=c(0.01,0.01,0.1,1,10,100),
                limits=c(0.01,100),
                oob = scales::oob_keep
  )
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression.NeverSmoker.pdf"),height = 2.3,width = 3.5)
##################################################################################################


##################################################################################################
### breakdown into high quality smoking ###
analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Smoker-SBS4") %>% filter(Tumor_Barcode %in% hqsample) #for smoking only
ecDNA_data=analysis_data[analysis_data$Tumor_Barcode %in% c(ecDNA_samples, no_cgr_samples),]
myformula <- as.formula(paste0("ecDNA","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
ecDNA_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=ecDNA_data)
bfb_data=analysis_data[analysis_data$Tumor_Barcode %in% c(bfb_samples, no_cgr_samples),]
myformula <- as.formula(paste0("bfb","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
bfb_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=bfb_data)
large_loss_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_loss_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_loss","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
large_loss_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_loss_data)
micronuclei_data=analysis_data[analysis_data$Tumor_Barcode %in% c(micronuclei_samples, no_cgr_samples),]
myformula <- as.formula(paste0("micronuclei","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
micronuclei_regression_result <- glm(formula=myformula,  family=binomial(link="logit"), data=micronuclei_data)
large_gain_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_gain_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_gain","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
large_gain_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_gain_data)
hourglass_data=analysis_data[analysis_data$Tumor_Barcode %in% c(hourglass_samples, no_cgr_samples),]
myformula <- as.formula(paste0("hourglass","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
hourglass_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=hourglass_data)
myformula <- as.formula(paste0("any","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
any_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=analysis_data)
clustered_function()[[1]]+
  scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                labels=c(0.01,0.01,0.1,1,10,100),
                limits=c(0.01,100),
                oob = scales::oob_keep
  )
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression.Smoker.hq.pdf"),height = 1.8,width = 5)

analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Smoker-SBS4")   %>% filter(Tumor_Barcode %in% hqsample) #for smoking only
chromoplexy_data=analysis_data[analysis_data$Tumor_Barcode %in% c(chromoplexy_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("chromoplexy","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
chromoplexy_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=chromoplexy_data)
cycle_templated_ins_data=analysis_data[analysis_data$Tumor_Barcode %in% c(cycle_templated_ins_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("cycle_templated_ins","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
cycle_templated_ins_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=cycle_templated_ins_data)
complex_unclear_data=analysis_data[analysis_data$Tumor_Barcode %in% c(complex_unclear_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("complex_unclear","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
complex_unclear_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=complex_unclear_data)
any_nonclustered_data <- analysis_data
myformula <- as.formula(paste0("any_nonclustered","~Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+")))
any_nonclustered_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=any_nonclustered_data)
nonclustered_funciton()[[1]]+
  scale_x_log10(breaks=c(0.000001,0.0001,0.01,1,100,10000,1000000),
                labels=c(0.01,0.01,0.01,1,100,100,100),
                limits=c(0.0000001,10000000),
                oob = scales::oob_keep
  )
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression.Smoker.hq.pdf"),height = 1.8,width = 3.5)
##################################################################################################


##################################################################################################
### breakdown into high-quality non-smoking ###
analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4")  %>% filter(Tumor_Barcode %in% hqsample)
ecDNA_data=analysis_data[analysis_data$Tumor_Barcode %in% c(ecDNA_samples, no_cgr_samples),]
myformula <- as.formula(paste0("ecDNA","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
ecDNA_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=ecDNA_data)
bfb_data=analysis_data[analysis_data$Tumor_Barcode %in% c(bfb_samples, no_cgr_samples),]
myformula <- as.formula(paste0("bfb","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
bfb_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=bfb_data)
large_loss_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_loss_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_loss","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_loss_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_loss_data)
micronuclei_data=analysis_data[analysis_data$Tumor_Barcode %in% c(micronuclei_samples, no_cgr_samples),]
myformula <- as.formula(paste0("micronuclei","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
micronuclei_regression_result <- glm(formula=myformula,  family=binomial(link="logit"), data=micronuclei_data)
large_gain_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_gain_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_gain","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_gain_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_gain_data)
hourglass_data=analysis_data[analysis_data$Tumor_Barcode %in% c(hourglass_samples, no_cgr_samples),]
myformula <- as.formula(paste0("hourglass","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
hourglass_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=hourglass_data)
myformula <- as.formula(paste0("any","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=analysis_data)
clustered_function()[[1]]+
  scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                labels=c(0.01,0.01,0.1,1,10,100),
                limits=c(0.01,100),
                oob = scales::oob_keep
  )
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression.Neversmoker.hq.pdf"),height = 2.3,width = 5)

analysis_data <- analysis_data_function()
analysis_data <- analysis_data %>% filter(Smoking_SBS4 == "Non-Smoker-noSBS4")   %>% filter(Tumor_Barcode %in% hqsample) #for smoking only
chromoplexy_data=analysis_data[analysis_data$Tumor_Barcode %in% c(chromoplexy_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("chromoplexy","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
chromoplexy_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=chromoplexy_data)
cycle_templated_ins_data=analysis_data[analysis_data$Tumor_Barcode %in% c(cycle_templated_ins_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("cycle_templated_ins","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
cycle_templated_ins_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=cycle_templated_ins_data)
complex_unclear_data=analysis_data[analysis_data$Tumor_Barcode %in% c(complex_unclear_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("complex_unclear","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
complex_unclear_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=complex_unclear_data)
any_nonclustered_data <- analysis_data
myformula <- as.formula(paste0("any_nonclustered","~Assigned_Population+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_nonclustered_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=any_nonclustered_data)
nonclustered_funciton()[[1]]+
  scale_x_log10(breaks=c(0.01,0.01,0.1,1,10,100),
                labels=c(0.01,0.01,0.1,1,10,100),
                limits=c(0.01,100),
                oob = scales::oob_keep
  )
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression.Neversmoker.hq.pdf"),height = 2.3,width = 3.5)
##################################################################################################



############################### +Stage #################################
analysis_data <- analysis_data_function() %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4"))
analysis_data <- analysis_data %>% filter(Tumor_Barcode %in% hqsample) #for high quality
ecDNA_data=analysis_data[analysis_data$Tumor_Barcode %in% c(ecDNA_samples, no_cgr_samples),]
myformula <- as.formula(paste0("ecDNA","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
ecDNA_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=ecDNA_data)
bfb_data=analysis_data[analysis_data$Tumor_Barcode %in% c(bfb_samples, no_cgr_samples),]
myformula <- as.formula(paste0("bfb","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
bfb_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=bfb_data)
large_loss_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_loss_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_loss","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_loss_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_loss_data)
micronuclei_data=analysis_data[analysis_data$Tumor_Barcode %in% c(micronuclei_samples, no_cgr_samples),]
myformula <- as.formula(paste0("micronuclei","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
micronuclei_regression_result <- glm(formula=myformula,  family=binomial(link="logit"), data=micronuclei_data)
large_gain_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_gain_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_gain","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
large_gain_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_gain_data)
hourglass_data=analysis_data[analysis_data$Tumor_Barcode %in% c(hourglass_samples, no_cgr_samples),]
myformula <- as.formula(paste0("hourglass","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
hourglass_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=hourglass_data)
myformula <- as.formula(paste0("any","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=analysis_data)
clustered_function()[1]
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression_addstage_hq.pdf"),height = 2.7,width = 5)

analysis_data <- analysis_data_function() %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4"))
analysis_data <- analysis_data %>% filter(Tumor_Barcode %in% hqsample) #for high quality
chromoplexy_data=analysis_data[analysis_data$Tumor_Barcode %in% c(chromoplexy_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("chromoplexy","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
chromoplexy_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=chromoplexy_data)
cycle_templated_ins_data=analysis_data[analysis_data$Tumor_Barcode %in% c(cycle_templated_ins_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("cycle_templated_ins","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
cycle_templated_ins_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=cycle_templated_ins_data)
complex_unclear_data=analysis_data[analysis_data$Tumor_Barcode %in% c(complex_unclear_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("complex_unclear","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
complex_unclear_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=complex_unclear_data)
any_nonclustered_data <- analysis_data
myformula <- as.formula(paste0("any_nonclustered","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion"))
any_nonclustered_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=any_nonclustered_data)
nonclustered_funciton()[1]
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression_addstage_hq.pdf"),height = 2.7,width = 3.5)



################################# +driver.freq #################################
analysis_data <- analysis_data_function() %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4"))
analysis_data <- analysis_data %>% filter(Tumor_Barcode %in% hqsample) #for high quality
ecDNA_data=analysis_data[analysis_data$Tumor_Barcode %in% c(ecDNA_samples, no_cgr_samples),]
myformula <- as.formula(paste0("ecDNA","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
ecDNA_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=ecDNA_data)
bfb_data=analysis_data[analysis_data$Tumor_Barcode %in% c(bfb_samples, no_cgr_samples),]
myformula <- as.formula(paste0("bfb","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
bfb_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=bfb_data)
large_loss_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_loss_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_loss","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
large_loss_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_loss_data)
micronuclei_data=analysis_data[analysis_data$Tumor_Barcode %in% c(micronuclei_samples, no_cgr_samples),]
myformula <- as.formula(paste0("micronuclei","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
micronuclei_regression_result <- glm(formula=myformula,  family=binomial(link="logit"), data=micronuclei_data)
large_gain_data=analysis_data[analysis_data$Tumor_Barcode %in% c(large_gain_samples, no_cgr_samples),]
myformula <- as.formula(paste0("large_gain","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
large_gain_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=large_gain_data)
hourglass_data=analysis_data[analysis_data$Tumor_Barcode %in% c(hourglass_samples, no_cgr_samples),]
myformula <- as.formula(paste0("hourglass","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
hourglass_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=hourglass_data)
myformula <- as.formula(paste0("any","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
any_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=analysis_data)
clustered_function()[1]
ggsave(file.path(outpath,"pcawg_signature_multivariate_logit_regression_driversv_hq.pdf"),height = 2.7,width = 5)

analysis_data <- analysis_data_function() %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4"))
analysis_data <- analysis_data %>% filter(Tumor_Barcode %in% hqsample) #for high quality
chromoplexy_data=analysis_data[analysis_data$Tumor_Barcode %in% c(chromoplexy_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("chromoplexy","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
chromoplexy_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=chromoplexy_data)
cycle_templated_ins_data=analysis_data[analysis_data$Tumor_Barcode %in% c(cycle_templated_ins_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("cycle_templated_ins","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
cycle_templated_ins_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=cycle_templated_ins_data)
complex_unclear_data=analysis_data[analysis_data$Tumor_Barcode %in% c(complex_unclear_samples, no_nonclustered_samples),]
myformula <- as.formula(paste0("complex_unclear","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
complex_unclear_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=complex_unclear_data)
any_nonclustered_data <- analysis_data
myformula <- as.formula(paste0("any_nonclustered","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+",paste0(paste0("factor(",mygene,")"), collapse = "+"),"+fusion","+driver.freq"))
any_nonclustered_regression_result <- glm(formula=myformula, family=binomial(link="logit"), data=any_nonclustered_data)
nonclustered_funciton()[1]
ggsave(file.path(outpath,"nonclustered_cplx_multivariate_logit_regression_driversv_hq.pdf"),height = 2.7,width = 3.5)





