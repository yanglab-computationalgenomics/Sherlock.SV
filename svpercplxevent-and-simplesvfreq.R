#### library
library(dplyr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(gridExtra)
library(sigminer)
library(tidyr)
library(grid)

p.adjust_function <- function(df,col=3,p.bfrn,p.fdr){
  df[,p.bfrn] <- p.adjust(df[,col],method = "bonferroni",n=nrow(df))
  df[,p.fdr] <- p.adjust(df[,col],method = "BH")
  return(df)
}
### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
simplesv.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_simple_sv_with_signature.tsv")
cluster.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_clustered_cmplx_sv.tsv")
noncluster.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_nonclustered_cmplx_sv.tsv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
plot.path <- file.path(scratch.path,"24-smoking-sbs4","plot")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")

if (!dir.exists(plot.path)){
  dir.create(plot.path,showWarnings = F,recursive = T)
}
### read
simplesv <- read.delim(simplesv.path)
cluster <- read.delim(cluster.path)
noncluster <- read.delim(noncluster.path)
sampleinfo <-  read.delim(sampleinfo.path)
sbs4 <- read.delim2(sbs4.path)
histology <- read.csv(histology.path)
sampleinfo <- merge(sampleinfo,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")


### add sbs4 to sampleinfo
sampleinfo <- merge(sampleinfo,sbs4)
sampleinfo$Smoking_SBS4 <- paste0(sampleinfo$Smoking,"-",ifelse(sampleinfo$SBS4==0,"noSBS4","SBS4"))
sampleinfo$Smoking_SBS4 <- ifelse(sampleinfo$Smoking_SBS4=="Non-Smoker-noSBS4","1.Non-Smoker-noSBS4",
                                  ifelse(sampleinfo$Smoking_SBS4=="Non-Smoker-SBS4","2.Non-Smoker-SBS4",
                                         ifelse(sampleinfo$Smoking_SBS4=="Smoker-noSBS4","3.Smoker-noSBS4",
                                                ifelse(sampleinfo$Smoking_SBS4=="Smoker-SBS4","4.Smoker-SBS4",sampleinfo$Smoking_SBS4)))
                                  )
### merge
allsv <- rbind(simplesv %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                          "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","Signature") %>% 
                 setnames(.,old="Signature",new="Signature") %>% mutate(CLUSTER_ID=paste0(SAMPLE,"_",Signature)),
               cluster %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                         "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","CGR_signature","CLUSTER_ID") %>% 
                 setnames(.,old="CGR_signature",new="Signature") ,
               noncluster %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                            "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","FULL_SV_TYPE","CLUSTER_ID")%>% 
                 setnames(.,old="FULL_SV_TYPE",new="Signature") %>%
                 mutate(Signature=ifelse(Signature %in% c("chromoplexy_del","chromoplexy_bal"),"chromoplexy",Signature),
                        CLUSTER_ID=paste0(SAMPLE,"_",CLUSTER_ID))) 
allsv <- merge(allsv,sampleinfo[,c("Tumor_Barcode","Assigned_Population","Gender","Smoking","SP_Group",
                                   "Histology","purity","wgd","scna_group","mut_load","indel_load","Smoking_SBS4","Inclusion")],
               by.x = "SAMPLE",by.y = "Tumor_Barcode",all.x = T) %>% filter(Inclusion=="In")
#write.table(allsv,file.path(scratch.path,"24-smoking-sbs4","allsv_signature_cluster_annotation.txt"),quote = F,row.names = F,col.names = T,sep = "\t")
### SV to breakends
breakpoint <- rbind(allsv %>% dplyr::select(SAMPLE,CHROM1,POS1,SVTYPE,STRAND1,Signature,CLUSTER_ID,Assigned_Population,Gender,Smoking,
                                              Histology,purity,wgd,scna_group,mut_load,indel_load,Smoking_SBS4) %>%
                      setnames(.,old = c("CHROM1","POS1","STRAND1"),new = c("Chromosome","POS","STRAND")),
                    allsv %>% dplyr::select(SAMPLE,CHROM2,POS2,SVTYPE,STRAND2,Signature,CLUSTER_ID,Assigned_Population,Gender,Smoking,
                                              Histology,purity,wgd,scna_group,mut_load,indel_load,Smoking_SBS4) %>%
                      setnames(.,old = c("CHROM2","POS2","STRAND2"),new = c("Chromosome","POS","STRAND")))
############################################# plot ############################################# 
breakpoint <- breakpoint %>% filter(Smoking_SBS4 %in% c("Unknown-noSBS4")==FALSE)
plot_brk <- breakpoint %>% 
  group_by(Smoking_SBS4,Signature,CLUSTER_ID) %>% 
  count() %>% 
  mutate(SP_Group_new=factor(Smoking_SBS4,levels = c("Non-Smoker-noSBS4","Non-Smoker-SBS4","Smoker-noSBS4","Smoker-SBS4")))
plot_chr <- breakpoint %>% 
  group_by(Smoking_SBS4,Signature,CLUSTER_ID,Chromosome) %>% count() %>% 
  group_by(Smoking_SBS4,Signature,CLUSTER_ID) %>% count() %>% 
  mutate(Smoking_SBS4=factor(Smoking_SBS4,levels = c("Non-Smoker-noSBS4","Non-Smoker-SBS4","Smoker-noSBS4","Smoker-SBS4")))

brk_plot_function <- function(myplot){
  out <- ggplot(data=myplot)+
    geom_violin(aes(x = Smoking_SBS4, y = n,fill=Smoking_SBS4))+
    geom_boxplot(aes(x = Smoking_SBS4, y = n),size=0.1,outlier.shape = NA,width=0.4) +
    ylab("Number of breakpoints")+
    scale_fill_manual(breaks = c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4"),
                      values = c("#00a800","#00a800","#FF6EC7","#FF6EC7"))+
    theme_bw()+
    theme(axis.text.x = element_text(size=6,colour = "black"),
          axis.text.y = element_text(size=6,colour = "black"),
          axis.title.y = element_text(size=6,colour = "black"),
          axis.title.x = element_blank(),
          axis.line.y = element_line(size=0.1,colour = "black"), 
          axis.line.x.top = element_line(size=0.1,colour = "black"), 
          axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
          axis.ticks.x =  element_line(linewidth = 0.234),
          axis.ticks.length.y =unit(0.8, "mm"),
          axis.ticks.y =  element_line(linewidth = 0.234),
          plot.title = element_text(size=6,colour = "black"),
          panel.background = element_blank(),
          panel.spacing.x = unit(1, "mm"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0.5),
          strip.text.y.left = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0),
          strip.background =  element_blank(),
          strip.placement = "top",
          legend.key.height = unit(2,"mm"),
          legend.key.width = unit(2,"mm"),
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.position = "bottom",
          plot.margin=grid::unit(c(1,2,0,2), "mm"))
  return(out)
}

comparisons_list <- list(
  c("1.Non-Smoker-noSBS4", "2.Non-Smoker-SBS4"),
  c("3.Smoker-noSBS4", "4.Smoker-SBS4"),
  c("1.Non-Smoker-noSBS4", "3.Smoker-noSBS4"),
  c("2.Non-Smoker-SBS4", "4.Smoker-SBS4"),
  c("1.Non-Smoker-noSBS4", "4.Smoker-SBS4")
)
############################ 1. clustered cplx sigs breakpoints #######################
clustered.cplx.sig <- c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass")
p.clusteredcplx=brk_plot_function(myplot=plot_brk %>% 
                                    filter(Signature %in% clustered.cplx.sig)) +
  stat_pvalue_manual(get_adj_p(plot_brk %>% filter(Signature %in% clustered.cplx.sig),
                               comparisons=comparisons_list,
                               method = "wilcox.test",
                               p.adjust.method = "fdr",
                               .col = "n", .grp = "Smoking_SBS4") %>% as.data.frame() %>%
                       p.adjust_function(.,col=4,p.bfrn="p.bfrn",p.fdr="p.fdr") %>%
                       mutate(p.fdr = ifelse(as.numeric(p.fdr) < 0.001,
                                             formatC(as.numeric(p.fdr), format = "e", digits = 1),
                                             formatC(as.numeric(p.fdr), format = "g", digits = 2))
                       ) %>% 
                       mutate(y.position = log10(y.position)+0.2),
                     size=2,
                     label = "p.fdr")+
  ggtitle("Clutered cplx SVs")+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 90))
plot_brk %>% filter(Signature %in% clustered.cplx.sig) %>% 
  dplyr::group_by(Smoking_SBS4) %>% count()
plot_brk %>% filter(Signature %in% clustered.cplx.sig) %>% 
  dplyr::group_by(Smoking_SBS4) %>% mutate(median=median(n)) %>% dplyr::select(Smoking_SBS4,median) %>% unique()

############################ 2. non-clustered cplx sig breakpoints ############################ 
noncluster.sig <- c("cycle_templated_ins","chromoplexy","complex_unclear")
p.nonclusteredcplx=brk_plot_function(myplot=plot_brk %>% 
                                       filter(Signature %in% noncluster.sig)) +
  stat_pvalue_manual(get_adj_p(plot_brk %>% filter(Signature %in% noncluster.sig),
                               comparisons=comparisons_list,
                               method = "wilcox.test",
                               p.adjust.method = "fdr",
                               .col = "n", .grp = "Smoking_SBS4") %>% as.data.frame() %>%
                       p.adjust_function(.,col=4,p.bfrn="p.bfrn",p.fdr="p.fdr") %>%
                       mutate(p.fdr = ifelse(as.numeric(p.fdr) < 0.001,
                                             formatC(as.numeric(p.fdr), format = "e", digits = 1),
                                             formatC(as.numeric(p.fdr), format = "g", digits = 2))
                       )%>% 
                       mutate(y.position = log10(y.position)+0.2),
                     size=2,
                     label = "p.fdr")+
  ggtitle("Non-clutered cplx SVs")+
  scale_y_log10(breaks=c(6,10,30,100,300))+
  theme(axis.text.x = element_text(angle = 90))
plot_brk %>% filter(Signature %in% noncluster.sig) %>% 
  dplyr::group_by(Smoking_SBS4) %>% count()
plot_brk %>% filter(Signature %in% noncluster.sig) %>% 
  dplyr::group_by(Smoking_SBS4) %>% mutate(median=median(n)) %>% dplyr::select(Smoking_SBS4,median) %>% unique()

############################ 3. simple sv signatures ############################ 
simple.sig <- c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra")
simple.df <- merge(sampleinfo,
                   allsv %>% filter(Signature %in% simple.sig) %>%
                     count(SAMPLE,Signature) %>% tidyr::spread(.,Signature,n),
                   by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across(simple.sig, ~replace(., is.na(.), 0.6))) %>%
  tidyr::gather(.,Signature,n,simple.sig) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4")),
         Signature=factor(Signature,levels=simple.sig))
p.simple=brk_plot_function(myplot= simple.df) +
  coord_cartesian()+
  facet_grid(.~Signature)+
  stat_pvalue_manual(simple.df %>%
                       group_by(Signature) %>%
                       group_modify( ~ {
                         # Apply the test on the subset for each Signature
                         res <- get_adj_p(
                           .x,
                           comparisons = comparisons_list,
                           method = "wilcox.test",
                           p.adjust.method = "fdr",
                           .col = "n",
                           .grp = "Smoking_SBS4"
                         ) %>%
                           as.data.frame() %>%
                           p.adjust_function(col = 4,
                                             p.bfrn = "p.bfrn",
                                             p.fdr = "p.fdr")
                         res
                       }) %>% mutate(p.fdr = ifelse(as.numeric(p.fdr) < 0.001,
                                                    formatC(as.numeric(p.fdr), format = "e", digits = 1),
                                                    formatC(as.numeric(p.fdr), format = "g", digits = 2))
                       ) %>% as.data.frame()%>% 
                       mutate(y.position = log10(y.position)+0.2),
                     size=2,
                     label = "p.fdr")+
  ggtitle("Simple SV Signatures")+
  labs(y="Number of SVs")+
  scale_y_log10(breaks=c(0.6,1,10,100,1000))+
  theme(axis.text.x = element_text(angle = 90))

x=merge(sampleinfo,
      allsv %>% filter(Signature %in% simple.sig) %>%
        count(SAMPLE,Signature) %>% spread(.,Signature,n),
      by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across(simple.sig, ~replace(., is.na(.), 0))) %>%
  gather(.,Signature,n,simple.sig) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4"))) %>%
  filter(Signature %in% simple.sig) %>% 
  dplyr::group_by(Smoking_SBS4,Signature) %>% 
  mutate(median=median(n)) %>% 
  dplyr::select(Smoking_SBS4,median) %>% unique() %>%
  arrange(Signature,Smoking_SBS4)

############################ 4. all sv  ############################ 
allsv.df <- merge(sampleinfo,
                   allsv %>% count(SAMPLE),
                   by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across("n", ~replace(., is.na(.), 0.5))) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4")))
p.allsv=brk_plot_function(myplot= allsv.df) +
  coord_cartesian()+
  stat_pvalue_manual(get_adj_p(allsv.df,
                               comparisons=comparisons_list,
                               method = "wilcox.test",
                               p.adjust.method = "fdr",
                               .col = "n", .grp = "Smoking_SBS4") %>% as.data.frame() %>%
                       p.adjust_function(.,col=4,p.bfrn="p.bfrn",p.fdr="p.fdr") %>%
                       mutate(p.fdr = ifelse(as.numeric(p.fdr) < 0.001,
                                             formatC(as.numeric(p.fdr), format = "e", digits = 1),
                                             formatC(as.numeric(p.fdr), format = "g", digits = 2))
                       )%>% 
                       mutate(y.position = log10(y.position)+0.2),
                     size=2,
                     label = "p.fdr")+
  ggtitle("All SVs")+
  labs(y="Number of SVs")+
  scale_y_log10(breaks=c(0.5,1,10,100,1000))+
  theme(axis.text.x = element_text(angle = 90))

merge(sampleinfo,
      allsv %>% count(SAMPLE),
      by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across("n", ~replace(., is.na(.), 0))) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4"))) %>%
  group_by(Smoking_SBS4) %>%
  mutate(median=median(n)) %>%
  dplyr::select(Smoking_SBS4,median) %>% unique()
############################ 5. all simple SVs  ############################ 
simple.sig <- c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra")
allsimple.df <- merge(sampleinfo,
                  allsv %>% filter(Signature %in% simple.sig) %>% count(SAMPLE),
                  by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across("n", ~replace(., is.na(.), 0.5))) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4")))
p.allsimple=brk_plot_function(myplot= allsimple.df) +
  coord_cartesian()+
  stat_pvalue_manual(get_adj_p(allsimple.df,
                               comparisons=comparisons_list,
                               method = "wilcox.test",
                               p.adjust.method = "fdr",
                               .col = "n", .grp = "Smoking_SBS4") %>% as.data.frame() %>%
                       p.adjust_function(.,col=4,p.bfrn="p.bfrn",p.fdr="p.fdr") %>%
                       mutate(p.fdr = ifelse(as.numeric(p.fdr) < 0.001,
                                             formatC(as.numeric(p.fdr), format = "e", digits = 1),
                                             formatC(as.numeric(p.fdr), format = "g", digits = 2))
                       )%>% 
                       mutate(y.position = log10(y.position)+0.2),
                     size=2,
                     label = "p.fdr")+
  ggtitle("All Simple SVs")+
  labs(y="Number of SVs")+
  scale_y_log10(breaks=c(0.5,1,10,100,1000))+
  theme(axis.text.x = element_text(angle = 90))

merge(sampleinfo,
      allsv %>% filter(Signature %in% simple.sig) %>% count(SAMPLE),
      by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across("n", ~replace(., is.na(.), 0))) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4"))) %>%
  group_by(Smoking_SBS4) %>%
  mutate(median=median(n)) %>%
  dplyr::select(Smoking_SBS4,median) %>% unique()
############################ 6. all cplx SVs  ############################ 
cplx.sig <- c(clustered.cplx.sig,noncluster.sig)
allcplx.df <- merge(sampleinfo,
                      allsv %>% filter(Signature %in% cplx.sig) %>% count(SAMPLE),
                      by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across("n", ~replace(., is.na(.), 0.5))) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4")))
p.allcplx=brk_plot_function(myplot= allcplx.df) +
  coord_cartesian()+
  stat_pvalue_manual(get_adj_p(allcplx.df,
                               comparisons=comparisons_list,
                               method = "wilcox.test",
                               p.adjust.method = "fdr",
                               .col = "n", .grp = "Smoking_SBS4") %>% as.data.frame() %>%
                       p.adjust_function(.,col=4,p.bfrn="p.bfrn",p.fdr="p.fdr") %>%
                       mutate(p.fdr = ifelse(as.numeric(p.fdr) < 0.001,
                                             formatC(as.numeric(p.fdr), format = "e", digits = 1),
                                             formatC(as.numeric(p.fdr), format = "g", digits = 2))
                       )%>% 
                       mutate(y.position = log10(y.position)+0.2),
                     size=2,
                     label = "p.fdr")+
  ggtitle("All Complex SVs")+
  labs(y="Number of SVs")+
  scale_y_log10(breaks=c(0.5,1,10,100,1000))+
  theme(axis.text.x = element_text(angle = 90))
merge(sampleinfo,
      allsv %>% filter(Signature %in% cplx.sig) %>% count(SAMPLE),
      by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across("n", ~replace(., is.na(.), 0))) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4"))) %>%
  group_by(Smoking_SBS4) %>%
  dplyr::mutate(median=median(n)) %>%
  dplyr::select(Smoking_SBS4,median) %>% unique()
############################ 7. complex sv signatures ############################ 
cplx.sig <- c(clustered.cplx.sig,noncluster.sig)
cplx.df <- merge(sampleinfo,
                   allsv %>% filter(Signature %in% cplx.sig) %>%
                     count(SAMPLE,Signature) %>% spread(.,Signature,n),
                   by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across(cplx.sig, ~replace(., is.na(.), 0.5))) %>%
  gather(.,Signature,n,cplx.sig) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,
                             levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4")),
         Signature=factor(Signature,levels=cplx.sig))
p.cplx=brk_plot_function(myplot= cplx.df %>% mutate(n=ifelse(n==0,n+0.1,n))) +
  facet_grid(.~Signature)+
  stat_pvalue_manual(cplx.df %>%
                       group_by(Signature) %>%
                       group_modify( ~ {
                         # Apply the test on the subset for each Signature
                         res <- get_adj_p(
                           .x,
                           comparisons = comparisons_list,
                           method = "wilcox.test",
                           p.adjust.method = "fdr",
                           .col = "n",
                           .grp = "Smoking_SBS4"
                         ) %>%
                           as.data.frame() %>%
                           p.adjust_function(col = 4,
                                             p.bfrn = "p.bfrn",
                                             p.fdr = "p.fdr")
                         res
                       }) %>% 
                       mutate(p.fdr = ifelse(as.numeric(p.fdr) < 0.001,
                                             formatC(as.numeric(p.fdr), format = "e", digits = 1),
                                             formatC(as.numeric(p.fdr), format = "g", digits = 2))
                       ) %>% as.data.frame()%>% 
                       mutate(y.position = log10(y.position)+0.2),
                     size=2,
                     label = "p.fdr")+
  ggtitle("Complex SV Signature")+
  labs(y="Number of SVs")+
  scale_y_log10(breaks=c(0.5,1,10,100,1000))+
  theme(axis.text.x = element_text(angle = 90))

x=merge(sampleinfo,
        allsv %>% filter(Signature %in% cplx.sig) %>%
          count(SAMPLE,Signature) %>% spread(.,Signature,n),
        by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T) %>%
  mutate(across(cplx.sig, ~replace(., is.na(.), 0))) %>%
  gather(.,Signature,n,cplx.sig) %>%
  filter(Smoking_SBS4!="Unknown-noSBS4") %>%
  mutate(Smoking_SBS4=factor(Smoking_SBS4,levels=c("1.Non-Smoker-noSBS4","2.Non-Smoker-SBS4","3.Smoker-noSBS4","4.Smoker-SBS4"))) %>%
  filter(Signature %in% cplx.sig) %>% 
  dplyr::group_by(Smoking_SBS4,Signature) %>% 
  mutate(median=median(n)) %>% 
  dplyr::select(Smoking_SBS4,median) %>% unique() %>%
  arrange(Signature,Smoking_SBS4)

### ggarrage
p <- grid.arrange(
  grid.arrange(p.allsv,p.allsimple,p.allcplx,nullGrob(),nrow=1),
  p.simple,
  nrow = 2,
  heights = c(2, 3)
)
ggsave(plot=p,file.path(plot.path,paste0("sv.smoking_sbs4.400",".pdf")),width = 170,height = 400,units = c("mm")) 
ggsave(plot=p,file.path(plot.path,paste0("sv.smoking_sbs4.300",".pdf")),width = 170,height = 300,units = c("mm")) 


### ggarrage
p <- grid.arrange(
  p.simple,
  grid.arrange(p.clusteredcplx, p.nonclusteredcplx, p.clusteredcplx, p.nonclusteredcplx, nrow = 1),
  nrow = 2
)
ggsave(plot=p,file.path(plot.path,paste0("breakpoint.per.event.smoking_sbs4",".pdf")),width = 160,height = 120,units = c("mm")) 
ggsave(plot=p,file.path(plot.path,paste0("breakpoint.per.event.smoking_sbs4.height300",".pdf")),width = 160,height = 300,units = c("mm")) 
