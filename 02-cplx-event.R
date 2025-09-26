#### library
library(dplyr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(ggpubr)
library(rstatix)
library(broom)

### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
simplesv.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_simple_sv_with_signature.tsv")
cluster.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_clustered_cmplx_sv.tsv")
noncluster.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_nonclustered_cmplx_sv.tsv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
plot.path <- file.path(scratch.path,"13-hotspot","plot")

### read
simplesv <- read.delim(simplesv.path)
cluster <- read.delim(cluster.path)
noncluster <- read.delim(noncluster.path)
sampleinfo <-  read.delim(sampleinfo.path)
sbs4 <- read.delim2(sbs4.path)
histology <- read.csv(histology.path)
sampleinfo <- merge(sampleinfo,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

### merge
allsv <- rbind(simplesv %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                   "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","Signature") %>% 
                 setnames(.,old="Signature",new="Signature") %>% mutate(CLUSTER_ID=row.names(.)),
               cluster %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                  "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","CGR_signature","CLUSTER_ID") %>% 
                 setnames(.,old="CGR_signature",new="Signature") ,
               noncluster %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                     "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","FULL_SV_TYPE","CLUSTER_ID")%>% 
                 setnames(.,old="FULL_SV_TYPE",new="Signature") %>%
                 mutate(Signature=ifelse(Signature %in% c("chromoplexy_del","chromoplexy_bal"),"chromoplexy",Signature),
                        CLUSTER_ID=paste0(SAMPLE,"_",CLUSTER_ID))) 
allsv <- allsv %>% filter(SAMPLE %in% sampleinfo$Tumor_Barcode)
allsv <- merge(allsv,sampleinfo[,c("Tumor_Barcode","Assigned_Population","Gender","Smoking","SP_Group",
                                     "Histology","purity","wgd","scna_group","mut_load","indel_load")],
               by.x = "SAMPLE",by.y = "Tumor_Barcode",all.x = T)
allsv <- merge(allsv,sbs4,by.x = "SAMPLE",by.y = "Tumor_Barcode",all.x = T)
allsv <- allsv %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))

allsv <- allsv %>% mutate(SP_Group_new=ifelse(Assigned_Population=="EUR" & Smoking_SBS4=="Smoker-SBS4","S_U",
                                              ifelse(Assigned_Population=="EUR" & Smoking_SBS4=="Non-Smoker-noSBS4", "N_U",
                                                     ifelse(Assigned_Population=="EAS" & Smoking_SBS4=="Non-Smoker-noSBS4", "N_A","Others"))
                                              ))


sampleinfo <- merge(sampleinfo,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
sampleinfo <- sampleinfo %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
sampleinfo %>% mutate(SP_Group_new=ifelse(Assigned_Population=="EUR" & Smoking_SBS4=="Smoker-SBS4","S_U",
                                          ifelse(Assigned_Population=="EUR" & Smoking_SBS4=="Non-Smoker-noSBS4", "N_U",
                                                 ifelse(Assigned_Population=="EAS" & Smoking_SBS4=="Non-Smoker-noSBS4", "N_A","Others")))) %>%
  count(SP_Group_new)

### complex event
cplx_event <- c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss",
                "4 Micronuclei","5 Large gain","6 Hourglass",
                "cycle_templated_ins","chromoplexy","complex_unclear")
cplx_sv <- allsv %>% filter(Signature %in% cplx_event )

### SV to breakends
breakpoint <- rbind(cplx_sv %>% dplyr::select(SAMPLE,CHROM1,POS1,SVTYPE,STRAND1,Signature,CLUSTER_ID,Assigned_Population,Gender,Smoking,SP_Group_new,
                                              Histology,purity,wgd,scna_group,mut_load,indel_load) %>%
                      setnames(.,old = c("CHROM1","POS1","STRAND1"),new = c("Chromosome","POS","STRAND")),
                    cplx_sv %>% dplyr::select(SAMPLE,CHROM2,POS2,SVTYPE,STRAND2,Signature,CLUSTER_ID,Assigned_Population,Gender,Smoking,SP_Group_new,
                                              Histology,purity,wgd,scna_group,mut_load,indel_load) %>%
                      setnames(.,old = c("CHROM2","POS2","STRAND2"),new = c("Chromosome","POS","STRAND")))


### plot_df
plot_brk <- breakpoint %>% 
  filter(SP_Group_new %in% c("N_A","N_U","S_U")) %>% 
  group_by(SP_Group_new,Signature,CLUSTER_ID) %>% 
  count() %>% 
  mutate(SP_Group_new=factor(SP_Group_new,levels = c("N_A","N_U","S_U")))
plot_chr <- breakpoint %>% 
  filter(SP_Group_new %in% c("N_A","N_U","S_U")) %>% 
  group_by(SP_Group_new,Signature,CLUSTER_ID,Chromosome) %>% 
  count() %>% 
  group_by(SP_Group_new,Signature,CLUSTER_ID) %>% count() %>% 
  mutate(SP_Group_new=factor(SP_Group_new,levels = c("N_A","N_U","S_U")))


### plot
brk_plot_function <- function(myplot, stat.test){
  ymix=min(myplot$n)
  out <- ggplot(data=myplot, aes(x = SP_Group_new, y = n))+
    geom_violin(aes(fill = SP_Group_new),trim = TRUE, color = NA) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    stat_pvalue_manual(
      stat.test,
      label = "p.adj.label", # FDR-corrected p-value
      group.by = "Signature",
      tip.length = 0.1,
      size = 2,
      color = "label.color" 
    ) +
    ylab("Number of breakpoints")+
    coord_trans(y = "log10")+
    scale_y_continuous(
      breaks = c(1, 10, 25,50,100,1000),
      labels = c(1, 10, 25,50,100,1000)
      # limits = c(1,1000)
    ) +
    scale_fill_manual(values = c("N_A"="#00a800","N_U"="#00a800","S_U"="#FF6EC7","N_A+N_U"="#00a800"))+
    scale_color_manual(values = c("black"="black","red"="red"))+
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
stat.test_function <- function(myfilter){
  testdata <- bind_rows(
    # N_A vs N_U
    myfilter %>%
      filter(SP_Group_new %in% c("N_A", "N_U")) %>%
      mutate(compare_group = "N_A vs N_U",
             test_group = SP_Group_new),
    # N_U vs S_U
    myfilter %>%
      filter(SP_Group_new %in% c("N_U", "S_U")) %>%
      mutate(compare_group = "N_U vs S_U",
             test_group = SP_Group_new),
    # (N_A + N_U) vs S_U
    myfilter %>%
      filter(SP_Group_new %in% c("N_A", "N_U", "S_U")) %>%
      mutate(compare_group = "(N_A+N_U) vs S_U",
             test_group = ifelse(SP_Group_new %in% c("N_A", "N_U"), "N_total", "S_U"))
  )
  stat.test <- testdata %>%
    group_by(compare_group) %>%
    group_modify(~{
      test <- wilcox.test(n ~ test_group, data = .x)
      tibble(
        group1 = unique(.x$test_group)[1],
        group2 = unique(.x$test_group)[2],
        p = test$p.value
      )
    }) %>%
    ungroup() %>%
    mutate(
      p.adj = p.adjust(p, method = "fdr"),
      p.adj.label = ifelse(
        p.adj < 1e-4,
        format(p.adj, scientific = TRUE, digits = 2),
        format(round(p.adj, 3), nsmall = 2)
      ),
      y.position = seq(from = max(myfilter$n, na.rm = TRUE) * 1.1,
                       by = max(myfilter$n, na.rm = TRUE) * 0.1,
                       length.out = n())
    ) %>%
    add_significance(p.col = "p.adj") %>% 
    mutate(label.color = ifelse(p.adj < 0.05, "red", "black"))
  return(stat.test)
}

plot_brk_filtered <- plot_brk %>% filter(Signature %in% c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass"))
stat.test <- stat.test_function(plot_brk_filtered)
p1 <- brk_plot_function(myplot = plot_brk_filtered, stat.test = stat.test) +
  ggtitle("Clustered cplx SVs") 
  

plot_brk_filtered <- plot_brk %>% filter(Signature %in% c("cycle_templated_ins","chromoplexy","complex_unclear"))
stat.test <- stat.test_function(plot_brk_filtered)
p2 <- brk_plot_function(myplot = plot_brk_filtered, stat.test = stat.test) +
  ggtitle("Non-clutered cplx SVs") 

pout=ggarrange(p1,p2,nrow = 1)
ggsave( plot = pout,file.path(plot.path,paste0("breakpoint.num.comparision.pop_smk",".pdf")),width = 80,height = 90,units = c("mm")) 

plot_brk %>% 
  filter(Signature %in% c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass")) %>% 
  dplyr::group_by(SP_Group_new) %>% count()
plot_brk %>% filter(Signature %in% c("cycle_templated_ins","chromoplexy","complex_unclear")) %>% 
  dplyr::group_by(SP_Group_new) %>% count()


plot_brk %>% filter(Signature %in% c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("cycle_templated_ins","chromoplexy","complex_unclear")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()

breakpoint %>% 
  #filter(SP_Group_new %in% c("N_A","N_U","S_U")) %>% 
  group_by(SP_Group_new,Signature,CLUSTER_ID,SAMPLE) %>% 
  count() %>%
  filter(Signature %in% c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge",
                          "3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass")) %>% as.data.frame() %>%
  dplyr::select(SP_Group_new,SAMPLE) %>% unique() %>% count(SP_Group_new)

breakpoint %>% 
  #filter(SP_Group_new %in% c("N_A","N_U","S_U")) %>% 
  group_by(SP_Group_new,Signature,CLUSTER_ID,SAMPLE) %>% 
  count() %>%
  filter(Signature %in% c("cycle_templated_ins","chromoplexy","complex_unclear")) %>% as.data.frame() %>%
  dplyr::select(SP_Group_new,SAMPLE) %>% unique() %>% count(SP_Group_new)


####################### signature by signature ######################
stat.test.signature_function <- function(myfilter){
  stat.test <- bind_rows(
    # N_A vs N_U
    myfilter %>%
      filter(SP_Group_new %in% c("N_A", "N_U")) %>%
      group_by(Signature) %>%
      mutate(compare_group = "N_A vs N_U", test_group = SP_Group_new),
    # N_U vs S_U
    myfilter %>%
      filter(SP_Group_new %in% c("N_U", "S_U")) %>%
      group_by(Signature) %>%
      mutate(compare_group = "N_U vs S_U", test_group = SP_Group_new),
    # (N_A + N_U) vs S_U
    myfilter %>%
      filter(SP_Group_new %in% c("N_A", "N_U", "S_U")) %>%
      group_by(Signature) %>%
      mutate(compare_group = "(N_A+N_U) vs S_U",
             test_group = ifelse(SP_Group_new %in% c("N_A", "N_U"), "N_total", "S_U"))) %>%
    group_by(Signature, compare_group) %>%
    group_modify(~ {
      test <- wilcox.test(n ~ test_group, data = .x)
      tibble(
        group1 = unique(.x$test_group)[1],
        group2 = unique(.x$test_group)[2],
        p = test$p.value
      )}) %>%
    group_by(Signature) %>%
    mutate(
      p.adj = p.adjust(p, method = "fdr")
    ) %>%
    ungroup() %>%
    mutate(
      p.adj.label = ifelse(
        p.adj < 1e-4,
        format(p.adj, scientific = TRUE, digits = 2),
        format(round(p.adj, 3), nsmall = 2)
      ),
      label.color = ifelse(p.adj < 0.05, "red", "black")
    ) %>%
    group_by(Signature) %>%
    mutate(y.position = seq(from = max(myfilter$n, na.rm = TRUE) * 1.1,
                            by = max(myfilter$n, na.rm = TRUE) * 0.1,
                            length.out = n()))
  return(stat.test)
}

plot_brk_filtered <- plot_brk %>% filter(Signature %in% c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass"))
stat.test <- stat.test.signature_function(plot_brk_filtered)
p1 <- brk_plot_function(myplot = plot_brk_filtered, stat.test = stat.test %>% as.data.frame()) +
  facet_grid(.~Signature)+
  ggtitle("Clustered cplx SVs") +
  scale_y_continuous(breaks = c(1,10,25,50,100,1000), labels = c(1,10,25,50,100,1000))

plot_brk_filtered <- plot_brk %>% filter(Signature %in% c("cycle_templated_ins","chromoplexy","complex_unclear"))
stat.test <- stat.test.signature_function(plot_brk_filtered)
p2 <- brk_plot_function(myplot = plot_brk_filtered, stat.test = stat.test %>% as.data.frame()) +
  facet_grid(.~Signature)+
  ggtitle("Non-clutered cplx SVs") +
  scale_y_continuous(breaks = c(1,10,25,50,100,1000), labels = c(1,10,25,50,100,1000))

ggarrange(p1,p2,widths = c(1,0.5),nrow = 1)
ggsave(file.path(plot.path,
                 paste0("breakpoint.num.comparision.pop_smk.signatures",".pdf")),width = 180,height = 50,units = c("mm")) 
ggsave(file.path(plot.path,
                 paste0("breakpoint.num.comparision.pop_smk.signatures.height300",".pdf")),width = 180,height = 300,units = c("mm"))


plot_brk %>% 
  filter(Signature %in% c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass")) %>% 
  dplyr::group_by(SP_Group_new,Signature) %>% count()
plot_brk %>% filter(Signature %in% c("cycle_templated_ins","chromoplexy","complex_unclear")) %>% 
  dplyr::group_by(SP_Group_new,Signature) %>% count()

plot_brk %>% filter(Signature %in% c("chromoplexy")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("cycle_templated_ins")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("1 ecDNA/double minutes")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("2 BFB cycles/chromatin bridge")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("3 Large loss")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("4 Micronuclei")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("5 Large gain")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("complex_unclear")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()
plot_brk %>% filter(Signature %in% c("6 Hourglass")) %>% 
  group_by(SP_Group_new) %>% mutate(median=median(n)) %>% dplyr::select(SP_Group_new,median) %>% unique()

# myplot=plot_chr
# ggboxplot(data=myplot, x = "pop_smk_group", y = "n",size=0.1,outlier.shape = 20,width=0.5,outlier.size=0.01)+
#   #geom_jitter(data=myplot, aes(x = pop_smk_group, y = n),fill=0.1,size=0.1,width = 0.2,alpha=0.1)+
#   # coord_trans(y = "log10")+
#   # scale_y_continuous(breaks = c(0,1,2,4,10,20,100))+
#   ylab("Number of Chromosomes")+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=6,colour = "black"),
#         axis.text.y = element_text(size=6,colour = "black"),
#         axis.title.y = element_text(size=6,colour = "black"),
#         axis.title.x = element_blank(),
#         axis.line.y = element_line(size=0.1,colour = "black"), 
#         axis.line.x.top = element_line(size=0.1,colour = "black"), 
#         axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
#         axis.ticks.x =  element_line(linewidth = 0.234),
#         axis.ticks.length.y =unit(0.8, "mm"),
#         axis.ticks.y =  element_line(linewidth = 0.234),
#         panel.background = element_blank(),
#         panel.spacing.x = unit(1, "mm"),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0.5),
#         strip.text.y.left = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0),
#         strip.background =  element_blank(),
#         strip.placement = "top",
#         legend.key.height = unit(2,"mm"),
#         legend.key.width = unit(2,"mm"),
#         legend.title = element_blank(),
#         legend.text = element_text(size=6),
#         legend.position = "bottom",
#         plot.margin=grid::unit(c(1,2,0,2), "mm"))+
#   stat_compare_means(method = "wilcox.test", label = "p.format", size=2,
#                      comparisons = list(
#                        c("N_A","N_U"),
#                        c("N_U","S_U"),
#                        c("N_A","S_U"))) # Add test result with p-value
# ggsave(file.path(plot.path,
#                  paste0("chr.num.comparision.pop_smk",".pdf")),width = 50,height = 30,units = c("mm")) 
# 
# myplot=plot_chr
# ggboxplot(data=myplot, x = "pop_smk_group", y = "n",size=0.1,outlier.shape = 20,width=0.5,outlier.size=0.01)+
#   #geom_jitter(data=myplot, aes(x = pop_smk_group, y = n),fill=0.1,size=0.1,width = 0.2,alpha=0.1)+
#   facet_grid(.~Signature)+
#   # coord_trans(y = "log10")+
#   # scale_y_continuous(breaks = c(0,1,2,4,10,20,100))+
#   ylab("Number of Chromosomes")+
#   theme_bw()+
#   theme(axis.text.x = element_text(size=6,colour = "black"),
#         axis.text.y = element_text(size=6,colour = "black"),
#         axis.title.y = element_text(size=6,colour = "black"),
#         axis.title.x = element_blank(),
#         axis.line.y = element_line(size=0.1,colour = "black"), 
#         axis.line.x.top = element_line(size=0.1,colour = "black"), 
#         axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
#         axis.ticks.x =  element_line(linewidth = 0.234),
#         axis.ticks.length.y =unit(0.8, "mm"),
#         axis.ticks.y =  element_line(linewidth = 0.234),
#         panel.background = element_blank(),
#         panel.spacing.x = unit(1, "mm"),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0.5),
#         strip.text.y.left = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0),
#         strip.background =  element_blank(),
#         strip.placement = "top",
#         legend.key.height = unit(2,"mm"),
#         legend.key.width = unit(2,"mm"),
#         legend.title = element_blank(),
#         legend.text = element_text(size=6),
#         legend.position = "bottom",
#         plot.margin=grid::unit(c(1,2,0,2), "mm"))+
#   stat_compare_means(method = "wilcox.test", label = "p.format", size=2,
#                      comparisons = list(
#                        c("N_A","N_U"),
#                        c("N_U","S_U"),
#                        c("N_A","S_U"))) # Add test result with p-value
# ggsave(file.path(plot.path,
#                  paste0("chr.num.comparision.pop_smk.signatures",".pdf")),width = 180,height = 30,units = c("mm")) 
# 
