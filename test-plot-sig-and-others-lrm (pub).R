library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape)
library(car)
library(tidyr)
library(data.table)
library(readxl)

######################################## FUNCTION ######################################
logistic_regression_plot_function <- function(
  signature.num =  "8",
  cutoff  =  "10",
  wgdscna = TRUE,
  smk = TRUE,
  mutationload  = TRUE,
  cosmic = FALSE,
  smkfilter=NA,
  stage=FALSE,
  sample_signature,
  mygenei=c("TP53","EGFR","KRAS"),
  outpath = file.path(
    sherlock_path,"scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature"),
  outputname  =  "all_sig_cutoff_10_with_WGDSCNA",
  hq=FALSE,
  plotheight=3.3,
  myscalex=scale_x_continuous(breaks = c(0, 1, 2),labels = c(0, 1, 2),limits= c(0, 2))
){
  cutoff <- ifelse(cutoff!="linear",as.numeric(cutoff),cutoff)
  if (!dir.exists(outpath)){
    dir.create(outpath,showWarnings = F,recursive = T)
  }
  
 
  sample_cgr <- complex_sv_with_signature %>% group_by(sample_id,CGR_signature,cluster_id) %>% count() %>% group_by(sample_id,CGR_signature) %>% count()
  sample_cgr_number_wide <- cast(sample_cgr,sample_id ~ CGR_signature)
  sample_cgr_number_wide[is.na(sample_cgr_number_wide)] <- 0
  
  unclustered_cmplx$FULL_SV_TYPE <- gsub("chromoplexy_bal|chromoplexy_del","chromoplexy",unclustered_cmplx$FULL_SV_TYPE)
  unclustered_cmplx <- unclustered_cmplx %>% group_by(SAMPLE,FULL_SV_TYPE,CLUSTER_ID) %>% count() %>% group_by(SAMPLE,FULL_SV_TYPE) %>% count()
  unclustered_cmplx_wide <- cast(unclustered_cmplx,SAMPLE ~ FULL_SV_TYPE)
  unclustered_cmplx_wide[is.na(unclustered_cmplx_wide)] <- 0
  
  ### merge files
  sample_signature <- sample_signature %>% count(SAMPLE,Signature) %>% spread(.,Signature,n) %>% replace(is.na(.), 0)
  rownames(sample_signature) <- sample_signature$SAMPLE
  sample_signature <- sample_signature[,2:(as.numeric(signature.num)+1)]
  plot_df <- sample_signature
  #  add no signature samples
  sample_with0_sv <- sample_info$Tumor_Barcode[sample_info$Tumor_Barcode %in% rownames(plot_df) == FALSE]
  plot_df_add <- data.frame(matrix(data=0,ncol = ncol(plot_df),nrow = length(sample_with0_sv)))
  colnames(plot_df_add) <- colnames(plot_df)
  rownames(plot_df_add) <- sample_with0_sv
  plot_df <- rbind(plot_df,
                   plot_df_add
  )
  signames  <-  names(plot_df)
  plot_df$sum <- apply(plot_df[,1:as.numeric(signature.num)],1,sum)
  # add sample name
  plot_df$sample <- rownames(plot_df)
  # add sample info
  plot_df <- merge(plot_df,sample_info,by.x = "sample",by.y="Tumor_Barcode",all.y = T)
  # add mygenes
  #mygene <- c("TP53","EGFR","KRAS")
  mygene <- mygenei
  plot_df <- merge(plot_df,snv[,c("sample",mygene)],by.x = "sample",by.y="sample",all.x = T)
  # add fusion
  plot_df <- plot_df %>% mutate(fusion=ifelse(sample %in% fusion.detail$SAMPLE,"Mut","WT"))
  # add cgr
  plot_df <- merge(plot_df,sample_cgr_number_wide,by.x = "sample",by.y = "sample_id",all.x = T,sort = F)
  # add unclustered complex
  plot_df <- merge(plot_df,unclustered_cmplx_wide,by.x = "sample",by.y = "SAMPLE",all.x = T,sort = F)
  # add sbs signature
  sig_cosmic_names_all <- names(sig_cosmic)
  sig_cosmic_clear <- sig_cosmic
  sig_cosmic_clear$Tumor_Barcode <- rownames(sig_cosmic_clear)
  sig_cosmic_names <- names(sig_cosmic_clear)[names(sig_cosmic_clear) %in% "Tumor_Barcode" ==  FALSE]
  plot_df <- merge(plot_df,sig_cosmic_clear,by.x = "sample",by.y = "Tumor_Barcode",all.x = T,sort = F)
  # add number of driver sv
  oncoprint_table$driver.freq <- rowSums(
    !is.na(oncoprint_table[, c(loss_gistic, gain_gistic)]) &
      oncoprint_table[, c(loss_gistic, gain_gistic)] != "NA"
  )
  plot_df <- merge(plot_df,oncoprint_table[,c("Tumor_Barcode","driver.freq")],all.x=T,by.x = "sample",by.y = "Tumor_Barcode")
  # add Smoking_SBS4
  plot_df <- merge(plot_df,sbs4 ,by.x = "sample",by.y = "Tumor_Barcode",all.x = T)
  plot_df <- plot_df %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4_new>0,"SBS4","noSBS4")))
  
  # filter sample
  plot_df <- plot_df %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4"))
  if (hq==TRUE){
    plot_df <- plot_df %>% filter(sample %in% hqsample)
  }
  if (smkfilter %in% c("Smoker")){
    plot_df <- plot_df %>% filter(Smoking=="Smoker") 
  }
  if (smkfilter %in% c("Non-Smoker")){
    plot_df <- plot_df %>% filter(Smoking=="Non-Smoker")
  }
  plot_df <- plot_df %>%
    filter(Assigned_Population %in% c("EUR","EAS"))
  print(paste0(nrow(plot_df)," samples in EUR/EAS with consistent smoking"))
  print(table(plot_df$Assigned_Population))
  # print(table(plot_df$fusion))
  
  sig_cosmic_names_keep <- sig_cosmic_names_all[apply(plot_df[,sig_cosmic_names_all], 2, function(x) sum(x>0)>10)]

  head(plot_df)
  # test
  regression_df <- plot_df %>% dplyr::select(sample,signames,sig_cosmic_names,
                                      Assigned_Population,Gender,Smoking,wgd,scna_group,mut_load,indel_load,all_of(mygene),
                                      fusion,Histology,purity,driver.freq,Stage
                                      )
    
  regression_df[,sig_cosmic_names] <- log10(regression_df[,sig_cosmic_names]+1)
  regression_df$mut_load <- log10(regression_df$mut_load)
  regression_df <- regression_df %>% 
    mutate(
           across(all_of(mygene), ~ ifelse(. == "None", "Wild", "Mut")),
           fusion=ifelse(fusion=="Mut","Mut","Wild"),
           Stage=ifelse(Stage %in% c("I","IA","IA1","IA2","IA3","IB"),"I",
                        ifelse(Stage %in% c("II","IIA","IIB"),"II",
                               ifelse(Stage %in% c("III","IIIA","IIIB"),"III",
                                      ifelse(Stage %in% c("IV","IVA"),"IV",NA))))
           ) %>%
    mutate(Assigned_Population=factor(Assigned_Population,levels = c("EAS","EUR")),
           Gender=factor(Gender,levels = c("Female","Male")),
           Smoking=factor(Smoking,levels = c("Non-Smoker","Smoker","Unknown")),
           wgd=factor(wgd,levels = c("nWGD","WGD")),
           scna_group=factor(scna_group,levels = c("Piano","Mezzo-forte","Forte")),
           across(all_of(mygene), ~ factor(.,levels = c("Wild","Mut"))),
           fusion=factor(fusion,levels = c("Wild","Mut")),
           Histology=factor(Histology,levels = c("Adenocarcinoma",
                                                 "Squamous cell carcinoma",
                                                 "Carcinoid tumor",
                                                 "Adenosquamous carcinoma",
                                                 "Others")),
           Stage=factor(Stage,levels=c("I","II","III","IV"))
           )

 
  if (cutoff=="linear"){
    regression_df[,signames] <- log10(regression_df[,signames]+1)
  }
  if (is.numeric(cutoff)){
    regression_df[,signames] <- ifelse(plot_df[,signames]>cutoff,
                                       1,
                                       0# ifelse(plot_df[,signames]==0,0,NA)
                                       )
  }
  
  star_p <- c(0.05,0.01,0.001)
  sig_name_list <- list(
    "8"= paste0("SBS49",c("B","A","D","H","F","C","E","G")),
    "9"= paste0("SBS49",c("B","A","C","H","I","F","E","D","G")))
  sig_fullname_list <- list(
    "8"= c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra"),
    "9"= c("Del1","Del2","Del3","Del4 + recip inv","TD1","TD2","Fb inv","Large mixed","Tra"))
  sig_name <-  unlist(sig_name_list[signature.num])
  sig_full <- unlist(sig_fullname_list[signature.num])
  names(sig_full) <- sig_name
  
  # logistic
  p_list <- list()
  count <- 0
  myfit_df_all <- data.frame()
  names(regression_df) <- gsub(" ",".",names(regression_df))

  for (i in signames){
    i=gsub(" ",".",i)
    count <- count + 1
    wgdscna_formula="+wgd+scna_group"
    smk_formula="+Smoking"
    mutationload_formula="+mut_load"
    cosmic_formula = paste0(sig_cosmic_names_keep, collapse = "+")
    stage_formula="+Stage"
    
    myformula <-
      as.formula(
        paste0(
          i,
          ifelse (smkfilter %in% c("Smoker") & hq==TRUE,"~Gender","~Assigned_Population+Gender"),
          ifelse(smk == FALSE, "", smk_formula),
          ifelse(wgdscna == TRUE, wgdscna_formula, ""),
          ifelse(mutationload == TRUE, mutationload_formula,""),
          ifelse(stage == TRUE, stage_formula,""),
          "+",paste0(mygene, collapse = "+"),
          "+",ifelse(cosmic == TRUE, cosmic_formula,""),
          "purity",
          ifelse (smkfilter %in% c("Smoker") & hq==TRUE,"","+fusion"),
          "+Histology"
          # "+driver.freq"
        )
      )
    
    # sample size
    sig_tell <- regression_df[,i] == 1
    num_sig_tell <- sum(sig_tell)
    num_sig_tell_rest <- nrow(regression_df)-num_sig_tell
    
    # logistiC
    if (is.numeric(cutoff)){
      fit <- glm(myformula, data = regression_df, family = "binomial")
    }
    # linear
    if ((cutoff=="linear")){
      fit <- glm(myformula, data = regression_df)
    }
    myfit_df <- as.data.frame(exp(cbind(OR = coef(fit), confint(fit))))
    myfit_df$variable <- rownames(myfit_df)
    p_df  <-  data.frame(coef(summary(fit)))
    p_df$variable <- rownames(p_df)
    names(p_df) <- c("estimate","Stderror","tvalue","pvalue","variable")
    myfit_df <- merge(myfit_df,p_df[,c("variable","pvalue")],all.x = T)
    myfit_df$star <- apply(myfit_df,1,function(x) paste0(rep("*",sum(as.numeric(x["pvalue"]) < star_p,na.rm = T)),collapse = ""))
    
    # add  tested sample num
    if (is.numeric(cutoff)){
      myfit_df$sample.num <- paste0(num_sig_tell,":",num_sig_tell_rest)
    } else {
      myfit_df$sample.num <- nrow(regression_df)
    }
    
    # multicollinearity
    print(myformula)
    # print(vif(fit))

    # sort
    myfit_df$variable <- factor(myfit_df$variable,
                                levels = rev(c("Assigned_PopulationEUR",
                                               "SmokingSmoker",
                                               "SmokingUnknown",
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
                                               "mut_load",
                                               "indel_load",
                                               # "TP53Mut",
                                               # "EGFRMut",
                                               # "KRASMut",
                                               c(paste0(mygene,"Mut")),
                                               "fusionMut",
                                               "driver.freq",
                                               sig_cosmic_names,
                                               "(Intercept)")))
    
    # rbind
    myfit_df$sig <-  i
    myfit_df_all <- rbind(myfit_df_all,myfit_df)
    
  }
  
  myfit_df_all$sig <- gsub("\\."," ",myfit_df_all$sig)
  myfit_df_all$sig <- factor(myfit_df_all$sig,levels = sig_full)
  myfit_df_all$OR <- round(myfit_df_all$OR,2)
  myfit_df_all$ORstar <- ifelse(grepl("\\*",myfit_df_all$star),paste0(myfit_df_all$OR," ",myfit_df_all$star) ,NA)  
  myfit_df_all$OR_color <- ifelse(myfit_df_all$OR>1,"red","blue")
  myfit_df_all$OR_color <- ifelse(grepl("\\*",myfit_df_all$star),myfit_df_all$OR_color ,"black")
  myfit_df_all <- myfit_df_all[myfit_df_all$variable %in% "(Intercept)" ==FALSE,]
  myvariable <- unique(myfit_df_all$variable)
  variable_fullname <- c("European VS Asian",
                         "Smoker VS Never smoker",
                         "Others VS Never smoker",
                         "Male VS Female",
                         "Adenosquamous carcinoma",
                         "Carcinoid tumor",
                         "Squamous cell carcinoma",
                         "Others",
                         "Tumor purity",
                         "Whole-genome duplication",
                         "Mezzo-forte",
                         "Forte",
                         "StageII","StageIII","StageIV",
                         "Mutation load",
                         "Indel load",
                         paste0(mygene," mutant VS Wildtype"),
                         "Fusion VS Wildtype",
                         "Driver SV num.",
                         sig_cosmic_names,
                         "Intercept")
  names(variable_fullname) <- c("Assigned_PopulationEUR",
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
                                "mut_load",
                                "indel_load",
                                paste0(mygene,"Mut"),
                                "fusionMut",
                                "driver.freq",
                                sig_cosmic_names,
                                "(Intercept)")
  myfit_df_all$variablefullname <- variable_fullname[as.character(myfit_df_all$variable)]
  xmax <- ifelse(max(myfit_df_all[,4],na.rm=T)>5,5,ceiling(max(myfit_df_all[,4],na.rm=T)))
  
  print(myfit_df_all)
  write.csv(myfit_df_all,file.path(outpath,paste0(outputname,".csv")),row.names = F,quote = F)
  
  # plot
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
                      "mut_load",
                      "indel_load",
                      paste0(mygene,"Mut"),
                      "fusionMut",
                      "driver.freq",
                      sig_cosmic_names))
                    )
  )
  ggplot(data=myfit_df_all_plot, 
         aes(y=variable, x=OR, xmin=`2.5 %`, xmax=`97.5 %`)) +
    geom_point(size=0.00001) +
    expand_limits(x = 5) + 
    myscalex+
    scale_y_discrete(breaks=unique(myfit_df_all_plot$variable),
                     labels=unique(myfit_df_all_plot$variablefullname))+
    facet_grid(.~sig,scales = "free") +
    scale_size_area(max_size = 1)+
    coord_cartesian(clip = 'on') +
    geom_vline(xintercept = 1, linetype="solid",color="black",lwd=0.2) + 
    geom_errorbarh(height=0,size=0.1) + 
    geom_text(data=myfit_df_all_plot, aes(color=OR_color,x=ifelse(OR>1,`2.5 %`,`97.5 %`),y=variable,label=ORstar,hjust = ifelse(OR>1,0,1)),vjust = -0.6,size=1.8)+
    scale_color_manual(values = c("red","blue","black"),
                       breaks = c("red","blue","black")) +
    xlab("Odds Ratio") + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=6,colour = "black"),
          axis.text.x = element_text(size=6,colour = "black"),
          axis.text.y = element_text(size=6,colour = "black"),
          strip.text.x = element_text(size=6,colour = "black"),
          axis.ticks = element_line(color = "black",size=0.1),
          axis.ticks.length=unit(0.02,"inch"),
          legend.position = "none",
          plot.title =element_text(size=6, face='bold',hjust = 0,vjust = 0),
          panel.background =element_blank(),
          panel.border=element_blank(), 
          axis.line=element_line(color = "black",size=0.1))
  
  ggsave(file.path(outpath,paste0(outputname,".pdf")),height = plotheight,width = 6)
}


####################################################
### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
sig8_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_simple_sv_with_signature.tsv")
hq_sample_path  <- file.path(sherlock.path,"Data/SAMPLE_INFO/HQ_samples.csv")
snv_path <- file.path(scratch.path,"06-SNV/sample_snv.txt")
sample_info_path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
fusion.detail.path  <- file.path(sherlock.path,"scratch/fusion/end_annotation/tier1_driver_fusions.txt")
complex_sv_with_signature_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/complex_sv_with_signature.tsv")
unclustered_cmplx_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_nonclustered_cmplx_sv.tsv")
sig_cosmic_path <-  file.path(sherlock.path,"Data/SNV/sherlock_sig_cosmic.txt")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.suppl.txt")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")

fusion_set <- c("ALK","RET","ROS1","MET","NRG1",
                "FGFR1","FGFR2","FGFR3","FGFR4",
                "NTRK1","NTRK2","NTRK3")

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
# read
hqsample <- read.csv(hq_sample_path) %>% pull(Tumor_Barcode)
snv  <-  read.delim(snv_path,row.names = 1)
sample_info <-  read.delim(sample_info_path)
histology <- read.csv(histology.path)
sample_info <- merge(sample_info,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")
complex_sv_with_signature <- read.delim(complex_sv_with_signature_path)
unclustered_cmplx <- read.delim(unclustered_cmplx_path)
fusion.detail <- read.delim2(fusion.detail.path)
fusion.detail <- fusion.detail %>%
  dplyr::mutate(
    Fusion_Gene = case_when(
      GENE1 %in% fusion_set ~ GENE1,
      GENE2 %in% fusion_set ~ GENE2,
      TRUE ~ NA_character_
    )
  ) %>% filter(Known_Fusion=="Yes")
snv_copy <- snv
snv$sample <- row.names(snv)
sig_cosmic <- read.delim(sig_cosmic_path,row.names = 1)
sample_signature <- read.delim(sig8_path)
oncoprint_table <- read.delim(oncoprint_table_path,check.names = F) %>%
  mutate(across(all_of(cnv_values), ~ ifelse((.x >= 5 & baseline_cn <= 2.7)|(.x >= 9 & baseline_cn > 2.7), "gain",
                                             ifelse((.x == 0 & baseline_cn <= 2.7)|(.x < (baseline_cn-2.7) & baseline_cn > 2.7),
                                                    "loss",
                                                    NA)))) %>%
  mutate(across(all_of(gain_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==3,"gain",NA)))) %>%
  mutate(across(all_of(loss_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==-3,"loss",NA))))
sbs4 <- read.delim2(sbs4.path) 
sbs4 <- sbs4 %>% setnames(.,old="SBS4",new="SBS4_new")

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
genes.5pct <- mutation_rate_df %>% filter(MutationRate >= 0.05) %>% pull(Gene) #144
genes.5pct <- c(priority_genes[priority_genes %in% genes.5pct], genes.5pct[!genes.5pct %in% priority_genes])
genes.10pct <- mutation_rate_df %>% filter(MutationRate >= 0.1) %>% pull(Gene) #17
genes.10pct <- c(priority_genes[priority_genes %in% genes.10pct], genes.10pct[!genes.10pct %in% priority_genes])

################################## RUN ################
# #modify the function manually first, about line 122
# for (cutoffi  in  c(5,10,20,"linear")[4]){
#   for (wgdscnai in c(FALSE,TRUE)[2]){
#     for (mutationloadi in c(FALSE,TRUE)[2]) {
#       for (smki in  c(FALSE,TRUE)[1]){
#         for (cosmici in  c(FALSE,TRUE)[1])
#         
#         hqi=FALSE
#         logistic_regression_plot_function(
#           signature.num =  "8",
#           cutoff  =  cutoffi,
#           wgdscna = wgdscnai,
#           smk = smki,
#           mutationload=mutationloadi,
#           smkfilter=NA,
#           sample_signature = sample_signature,
#           mygenei = c("TP53","EGFR","KRAS"),
#           cosmic=cosmici,
#           hq=hqi,
#           plotheight=3.3,
#           outpath = file.path(sherlock_path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
#           outputname  =  paste0("all_sig_cutoff_",cutoffi,
#                                 "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
#                                 "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
#                                 "_smk",grep(smki,c(FALSE,TRUE))-1,
#                                 "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
#                                 "_hq",grep(hqi,c(FALSE,TRUE))-1)
#         )
#         
#       }
#     }
#   }
# }
###################################################################################
# 1. priority_genes only in high quality samples
cutoffi="linear"
wgdscnai=TRUE
smki=TRUE
mutationloadi=FALSE
cosmici=FALSE
hqi=TRUE
stagei=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter=NA,
  sample_signature = sample_signature,
  mygenei = priority_genes,
  cosmic=cosmici,
  hq=hqi,
  stage = stagei,
  plotheight=2,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_stage",grep(stagei,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1
  ),
  myscalex = scale_x_continuous(breaks = c(0, 1, 2),labels = c(0, 1, 2),limits= c(0, 2))
)

# 2. priority_genes only in all samples
cutoffi="linear"
wgdscnai=TRUE
smki=TRUE
mutationloadi=FALSE
cosmici=FALSE
hqi=FALSE
stagei=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter=NA,
  sample_signature = sample_signature,
  mygenei = priority_genes,
  cosmic=cosmici,
  hq=hqi,
  plotheight=2.5,
  stage = stagei,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_stage",grep(stagei,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1
  ),
  myscalex = scale_x_continuous(breaks = c(0, 1, 2),labels = c(0, 1, 2),limits= c(0, 2))
)


# 3. priority_genes only in smokers
# 0.5,1,5,10
cutoffi="linear"
wgdscnai=TRUE
smki=FALSE
mutationloadi=FALSE
cosmici=FALSE
hqi=FALSE
stagei=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter="Smoker",
  sample_signature = sample_signature,
  mygenei = priority_genes,
  cosmic=cosmici,
  hq=hqi,
  stage = stagei,
  plotheight=1.7,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_stage",grep(stagei,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1,
                        "_smokers"
  ),
  myscalex = scale_x_log10(breaks=c(0.5,1,5,10),labels=c(0.5,1,5,10),limits=c(0.5,10))
)

# 4. priority_genes only in nonsmokers
# 0.5,1,2
cutoffi="linear"
wgdscnai=TRUE
smki=FALSE
mutationloadi=FALSE
cosmici=FALSE
hqi=FALSE
stagei=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter="Non-Smoker",
  sample_signature = sample_signature,
  mygenei = priority_genes,
  cosmic=cosmici,
  hq=hqi,
  stage = stagei,
  plotheight=1.9,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_stage",grep(stagei,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1,
                        "_nonsmokers"
  ),
  myscalex = scale_x_log10(breaks=c(0.5,1,2),labels=c(0.5,1,2),limits=c(0.5,2))
)

# 5. priority_genes only in high quality smokers
# 0.2,1,5
cutoffi="linear"
wgdscnai=TRUE
smki=FALSE
mutationloadi=FALSE
cosmici=FALSE
hqi=TRUE
stagei=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter="Smoker",
  sample_signature = sample_signature,
  mygenei = priority_genes,
  cosmic=cosmici,
  hq=hqi,
  stage = stagei,
  plotheight=1.4,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_stage",grep(stagei,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1,
                        "_smokers" 
  ),
  myscalex = scale_x_log10(breaks=c(0.2,1,5),labels=c(0.2,1,5),limits=c(0.2,5))
)


# 6. priority_genes only in high quality nonsmokers
# 0.5,1,3
cutoffi="linear"
wgdscnai=TRUE
smki=FALSE
mutationloadi=FALSE
cosmici=FALSE
hqi=TRUE
stagei=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter="Non-Smoker",
  sample_signature = sample_signature,
  mygenei = priority_genes,
  cosmic=cosmici,
  hq=hqi,
  stage = stagei,
  plotheight=1.9,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_stage",grep(stagei,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1,
                        "_nonsmokers"
  ),
  myscalex = scale_x_log10(breaks=c(0.5,1,3),labels=c(0.5,1,3),limits=c(0.5,3))
)

###################################################################################
# genes10pct
cutoffi="linear"
wgdscnai=TRUE
smki=TRUE
mutationloadi=FALSE
cosmici=FALSE
hqi=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter=NA,
  sample_signature = sample_signature,
  mygenei = genes.10pct,
  cosmic=cosmici,
  hq=hqi,
  plotheight=5,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1,
                        "_genes10pct"
                        )
)

# genes5pct
cutoffi="linear"
wgdscnai=TRUE
smki=TRUE
mutationloadi=FALSE
cosmici=FALSE
hqi=FALSE
logistic_regression_plot_function(
  signature.num =  "8",
  cutoff  =  cutoffi,
  wgdscna = wgdscnai,
  smk = smki,
  mutationload=mutationloadi,
  smkfilter=NA,
  sample_signature = sample_signature,
  mygenei = genes.5pct,
  cosmic=cosmici,
  hq=hqi,
  plotheight=20,
  outpath = file.path(sherlock.path,paste0("scratch/07-signature/1217_samples_hg38/Sherlock_SV49_simple_regressionplot/8signature_myannotation/cutoff",cutoffi)),
  outputname  =  paste0("all_sig_cutoff_",cutoffi,
                        "_wgdcsnc",grep(wgdscnai,c(FALSE,TRUE))-1,
                        "_mutload",grep(mutationloadi,c(FALSE,TRUE))-1,
                        "_smk",grep(smki,c(FALSE,TRUE))-1,
                        "_cosmic",grep(cosmici,c(FALSE,TRUE))-1,
                        "_hq",grep(hqi,c(FALSE,TRUE))-1,
                        "_genes5pct"
  )
)

