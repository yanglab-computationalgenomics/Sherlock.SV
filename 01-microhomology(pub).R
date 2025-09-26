### microhomology
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(ggh4x)

### function
homology.plot.function <- function(
  allsv_plot,
  category="smoking",
  category.name="Smoking",
  category.order=c("Smoker_SBS4","NonSmoker_NonSBS4","Others"),
  category.color=c("red","blue","grey50"),
  category.fill=c("red","blue","grey50"),
  myposition="stack"
){
  mysig <- c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass",
             "chromoplexy","cycle_templated_ins","complex_unclear",
             "Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra","un_assigned")
  mysig.full <- c("1 ecDNA","2 chr bridge","3 Large loss","4 Micronuclei","5 Large gain","6 Hourglass",
                  "Chromoplexy","Cycle of\ntemplated insertion","Complex unclear",
                  "Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra","un_assigned")
  names(mysig.full) <- mysig
  
  summaryplot <- allsv_plot %>% group_by_at(c("Signature","HOMINS",category)) %>% count() %>%
    mutate(Signature=factor(Signature,levels = mysig)) %>%
    setnames(.,old=c(category),new = c("category")) %>% as.data.frame() %>%
    mutate(category=factor(category,levels = category.order))
  
  p <- ggplot(data=summaryplot %>% filter(Signature!="un_assigned")) + 
    geom_bar(aes(x=HOMINS,y=n,fill=category,color=category),stat = "identity",position=myposition,width = 0.8,size=0.234/2)+
    facet_grid2(Signature~.,
                scales = "free",
                space = "free_x",
                switch = "both",
                labeller=as_labeller(mysig.full),
                strip = strip_vanilla(clip = "off"))+
    scale_x_continuous(expand = c(0.01, 0.01))+
    scale_fill_manual(name=category.name,
                      breaks = category.order,
                      values = category.fill,
                      drop=FALSE)+
    scale_color_manual(breaks = category.order,
                      values = category.color,
                      drop=FALSE)+
    guides(fill = guide_legend(nrow = 5))+
    ylab("Number of SVs")+
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
          panel.background = element_blank(),
          panel.spacing.x = unit(1, "mm"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size=6, angle = 45, colour = "black",hjust = 1,vjust = 0.5),
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
  if (myposition=="stack"){
    p <- p + scale_y_continuous(
      limits = ~ c(0, ceiling(max(.x)/100)*100),
      breaks = c(0,500,1000,2000,3000,4000,5000,6000,7000,8000),
      expand = c(0, 0)
    )
  }
  return(p)
}

### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
simplesv.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_simple_sv_with_signature.tsv")
cluster.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_clustered_cmplx_sv.tsv")
noncluster.path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_nonclustered_cmplx_sv.tsv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
outpath <- file.path(scratch.path,"11-microhomology","plot")
hq.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/HQ_samples.csv")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")

if (dir.exists(outpath)==FALSE){
  dir.create(outpath,showWarnings = F,recursive = T)
}
### read
simplesv <- read.delim(simplesv.path)
cluster <- read.delim(cluster.path)
noncluster <- read.delim(noncluster.path)
sampleinfo <-  read.delim(sampleinfo.path)
histology <- read.csv(histology.path)
sampleinfo <- merge(sampleinfo,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")
hqsample <- read.csv(hq.path) %>% pull(Tumor_Barcode)
tp53_mutated_samples=read.delim(tp53.mut.path, header=F, check.names=F, stringsAsFactors=F)
kras_mutated_samples=read.delim(kras.mut.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F)
sbs4 <- read.delim2(sbs4.path)
### merge
allsv <- rbind(simplesv %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                   "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","Signature") %>% 
                setnames(.,old="Signature",new="Signature") ,
               cluster %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                   "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","CGR_signature") %>% 
                 setnames(.,old="CGR_signature",new="Signature") ,
               noncluster %>% dplyr::select("SAMPLE","CHROM1","POS1","STRAND1","CHROM2","POS2","STRAND2","SVTYPE",
                                  "manta_HOMLEN","manta_SVINSLEN","meerkat_HOMLENG","meerkat_SVINSLENG","ALGORITHMS","FULL_SV_TYPE")%>% 
                 setnames(.,old="FULL_SV_TYPE",new="Signature") %>%
                 mutate(Signature=ifelse(Signature %in% c("chromoplexy_del","chromoplexy_bal"),"chromoplexy",Signature))
               ) %>% 
  mutate(HOMLEN=ifelse(ALGORITHMS %in% c("Meerkat-Manta","Meerkat"),meerkat_HOMLENG,manta_HOMLEN),
         SVINSLEN=ifelse(ALGORITHMS %in% c("Meerkat-Manta","Meerkat"),meerkat_SVINSLENG,manta_SVINSLEN)) %>% #use meerkat first
  mutate(HOMLEN=ifelse((ALGORITHMS %in% c("Manta")) & is.na(HOMLEN) & is.na(SVINSLEN),0,HOMLEN)) %>%   #manta na na > 0 na
  mutate(HOMINS=ifelse(is.na(HOMLEN),-SVINSLEN,HOMLEN)) %>%
  mutate(HOMINS=ifelse(HOMINS > 10,10,HOMINS)) %>%
  mutate(HOMINS=ifelse(HOMINS < -10,-10,HOMINS))

allsv_out <- allsv %>% mutate(SVID=paste0("SVID",rownames(.)))
write.table(allsv_out,
            file.path(scratch,"04-shatterseek/1217_samples_hg38","sherlock_all_sv_with_signature.tsv"),row.names = F,
            quote = FALSE,sep = "\t")

### add other info
allsv_plot <- merge(allsv,sampleinfo,by.x = "SAMPLE",by.y = "Tumor_Barcode",all.x = T) %>% filter(Inclusion=="In")
allsv_plot$TP53 <- ifelse(allsv_plot$SAMPLE %in% as.vector(unlist(tp53_mutated_samples)),"mut","WT")
allsv_plot$KRAS <- ifelse(allsv_plot$SAMPLE %in% as.vector(unlist(kras_mutated_samples)),"mut","WT")
allsv_plot$EGFR <- ifelse(allsv_plot$SAMPLE %in% as.vector(unlist(egfr_mutated_E479_A483del_samples)),"E479_A483del","WT")
allsv_plot$EGFR <- ifelse(allsv_plot$SAMPLE %in% as.vector(unlist(egfr_mutated_L591R_samples)),"L591R",allsv_plot$EGFR)
allsv_plot$EGFR <- ifelse(allsv_plot$SAMPLE %in% as.vector(unlist(egfr_mutated_others_samples)),"Other_mut",allsv_plot$EGFR)
allsv_plot <- merge(allsv_plot,sbs4,by.x = "SAMPLE",by.y = "Tumor_Barcode",all.x = T)
allsv_plot <- allsv_plot %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))


### plot
# Smoking_SBS4
homology.plot.function(allsv_plot %>% filter(Smoking_SBS4 %in%  c("Non-Smoker-noSBS4","Smoker-SBS4")),
                       category="Smoking_SBS4",
                       category.name="Smoking_SBS4",
                       category.order=c("Non-Smoker-noSBS4","Smoker-SBS4"),
                       category.color=c("#00a800","#FF6EC7"),
                       category.fill=c("#00a800","#FF6EC7"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.consistent.smoking.pdf"),width = 6.5,height = 6.5)


# 1. smoking
homology.plot.function(allsv_plot,
                       category="smoking",
                       category.name="Smoking",
                       category.order=c("Others","NonSmoker_NonSBS4","Smoker_SBS4"),
                       category.color=c("grey80","black","black"),
                       category.fill=c("grey80","white","black"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.smoking.pdf"),width = 3,height = 4)
# 2. population
homology.plot.function(allsv_plot,
                       category="population",
                       category.name="Population",
                       category.order=c("EUR","EAS","AMR or Mixed","AFR"),
                       category.color=c("#e96293","#7a69a2","#66bdbb","#f9ed5a"),
                       category.fill=c("#e96293","#7a69a2","#66bdbb","#f9ed5a"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.population.pdf"),width = 3,height = 4)
# 3. pop_smk_group
homology.plot.function(allsv_plot,
                       category="pop_smk_group",
                       category.name="pop_smk_group",
                       category.order=c("Others","N_A","N_U","S_U"),
                       category.color=c("grey80","#006400","#00c2af","#ff6ec7"),
                       category.fill=c("grey80","#006400","#00c2af","#ff6ec7"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.popsmk.pdf"),width = 3,height = 4)
# 4. histology
homology.plot.function(allsv_plot,
                       category="histology",
                       category.name="histology",
                       category.order=c("Others","Squamous carcinoma","Carcinoid tumor","Adenosquamous carcinoma","Adenocarcinoma"),
                       category.color=c("grey80","#c0d849","#bb56a0","#ff9540","#3cb44b"),
                       category.fill=c("grey80","#c0d849","#bb56a0","#ff9540","#3cb44b"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.histology.pdf"),width = 3,height = 4)
# 5. gender
homology.plot.function(allsv_plot,
                       category="gender",
                       category.name="Gender",
                       category.order=c("Female","Male"),
                       category.color=c("#b73377","#2da9d9"),
                       category.fill=c("#b73377","#2da9d9"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.gender.pdf"),width = 3,height = 4)
# 6. tp53
homology.plot.function(allsv_plot,
                       category="TP53",
                       category.name="TP53",
                       category.order=c("WT","mut"),
                       category.color=c("black","black"),
                       category.fill=c("white","black"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.tp53.pdf"),width = 3,height = 4)
homology.plot.function(allsv_plot %>% filter(TP53.smoking %in% c("mut-NonSmoker_NonSBS4","mut-Smoker_SBS4","WT-NonSmoker_NonSBS4","WT-Smoker_SBS4")),
                       category="TP53.smoking",
                       category.name="TP53",
                       category.order=c("mut-NonSmoker_NonSBS4","mut-Smoker_SBS4","WT-NonSmoker_NonSBS4","WT-Smoker_SBS4"),
                       category.color=c("pink","red","lightblue","blue"),
                       category.fill=c("pink","red","lightblue","blue"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.tp53.smoking.pdf"),width = 3,height = 4)
# 7. EGFR
homology.plot.function(allsv_plot %>% filter(EGFR.smoking %in% c("Other_mut-NonSmoker_NonSBS4","L591R-NonSmoker_NonSBS4","E479_A483del-NonSmoker_NonSBS4","WT-NonSmoker_NonSBS4",
                                                                 "Other_mut-Smoker_SBS4","L591R-Smoker_SBS4","E479_A483del-Smoker_SBS4","WT-Smoker_SBS4")),
                       category="EGFR.smoking",
                       category.name="EGFR",
                       category.order=c("WT-NonSmoker_NonSBS4",
                                        "WT-Smoker_SBS4","Other_mut-NonSmoker_NonSBS4","L591R-NonSmoker_NonSBS4","E479_A483del-NonSmoker_NonSBS4",
                                        "Other_mut-Smoker_SBS4","L591R-Smoker_SBS4","E479_A483del-Smoker_SBS4"
                       ),
                       category.color=c("grey80","grey80","#F13F19","#F13F19","#F13F19",
                                        "#F13F19","#F13F19","#F13F19"),
                       category.fill=c("grey80","grey80","#F13F19","#F13F19","#F13F19",
                                       "#F13F19","#F13F19","#F13F19"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.egfr.smoking.pdf"),width = 3,height = 4)
homology.plot.function(allsv_plot %>% filter(EGFR.smoking %in% c("Other_mut-NonSmoker_NonSBS4","L591R-NonSmoker_NonSBS4","E479_A483del-NonSmoker_NonSBS4","WT-NonSmoker_NonSBS4",
                                                                 "Other_mut-Smoker_SBS4","L591R-Smoker_SBS4","E479_A483del-Smoker_SBS4","WT-Smoker_SBS4")),
                       category="EGFR.smoking",
                       category.name="EGFR",
                       category.order=c("WT-NonSmoker_NonSBS4",
                                        "WT-Smoker_SBS4","Other_mut-NonSmoker_NonSBS4","L591R-NonSmoker_NonSBS4","E479_A483del-NonSmoker_NonSBS4",
                                        "Other_mut-Smoker_SBS4","L591R-Smoker_SBS4","E479_A483del-Smoker_SBS4"
                                        ),
                       category.color=c("grey80","grey80","#F13F19","#F13F19","#F13F19",
                                        "#F13F19","#F13F19","#F13F19"),
                       category.fill=c("grey80","grey80","#F13F19","#F13F19","#F13F19",
                                       "#F13F19","#F13F19","#F13F19"),
                       myposition="fill")
ggsave(file.path(outpath,"microhomology.egfr.smoking.fill.pdf"),width = 3,height = 4)
# 8. kras
homology.plot.function(allsv_plot,
                       category="KRAS",
                       category.name="KRAS",
                       category.order=c("WT","mut"),
                       category.color=c("black","black"),
                       category.fill=c("white","black"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.kras.pdf"),width = 3,height = 4)
# 9. smoking tp53 kras egfr
homology.plot.function(allsv_plot,
                       category="smoking.tp53.egfr.kras",
                       category.name="smoking.tp53.egfr.kras",
                       category.order=c("Nonsmoker-TP53mut-EGFRwt-KRASmut","Nonsmoker-TP53wt-EGFRwt-KRASmut","smoker-TP53mut-EGFRmut-KRASmut","smoker-TP53mut-EGFRmut-KRASwt","smoker-TP53mut-EGFRwt-KRASmut",
                                        "smoker-TP53wt-EGFRmut-KRASmut","smoker-TP53wt-EGFRmut-KRASwt",
                                        "smoker-TP53wt-EGFRwt-KRASmut","smoker-TP53wt-EGFRwt-KRASwt",
                                        "Nonsmoker-TP53mut-EGFRmut-KRASwt",
                                        "Nonsmoker-TP53mut-EGFRwt-KRASwt","Nonsmoker-TP53wt-EGFRmut-KRASwt",
                                        "Nonsmoker-TP53wt-EGFRwt-KRASwt",
                                        "smoker-TP53mut-EGFRwt-KRASwt"
                                        ),
                       category.color=c("grey90","grey90","grey90", "grey90", "grey90", "grey90", "grey90",
                                        "grey90", "grey90",
                                        "#1f78b4", 
                                        "yellow", "#ff7f00",
                                        "#A8EEF0", 
                                        "#cab2d6"
                                        ),
                       category.fill=c("grey90","grey90","grey90", "grey90", "grey90", "grey90", "grey90",
                                       "grey90", "grey90",
                                       "#1f78b4", 
                                       "yellow", "#ff7f00",
                                       "#A8EEF0", 
                                       "#cab2d6"),
                       # category.color=c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a", 
                       #                  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99", "#b15928", "#8dd3c7", "#fb8072"),
                       # category.fill=c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a", 
                       #                 "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99", "#b15928", "#8dd3c7", "#fb8072"),
                       myposition="fill")
ggsave(file.path(outpath,"microhomology.smoking.tp53.egfr.kras.fill.pdf"),width = 3,height = 4)
homology.plot.function(allsv_plot,
                       category="smoking.tp53.egfr.kras",
                       category.name="smoking.tp53.egfr.kras",
                       category.order=c("Nonsmoker-TP53mut-EGFRwt-KRASmut","Nonsmoker-TP53wt-EGFRwt-KRASmut","smoker-TP53mut-EGFRmut-KRASmut","smoker-TP53mut-EGFRmut-KRASwt","smoker-TP53mut-EGFRwt-KRASmut",
                                        "smoker-TP53wt-EGFRmut-KRASmut","smoker-TP53wt-EGFRmut-KRASwt",
                                        "smoker-TP53wt-EGFRwt-KRASmut","smoker-TP53wt-EGFRwt-KRASwt",
                                        "Nonsmoker-TP53mut-EGFRmut-KRASwt",
                                        "Nonsmoker-TP53mut-EGFRwt-KRASwt","Nonsmoker-TP53wt-EGFRmut-KRASwt",
                                        "Nonsmoker-TP53wt-EGFRwt-KRASwt",
                                        "smoker-TP53mut-EGFRwt-KRASwt"
                       ),
                       category.color=c("grey90","grey90","grey90", "grey90", "grey90", "grey90", "grey90",
                                        "grey90", "grey90",
                                        "#1f78b4", 
                                        "yellow", "#ff7f00",
                                        "#A8EEF0", 
                                        "#cab2d6"
                       ),
                       category.fill=c("grey90","grey90","grey90", "grey90", "grey90", "grey90", "grey90",
                                       "grey90", "grey90",
                                       "#1f78b4", 
                                       "yellow", "#ff7f00",
                                       "#A8EEF0", 
                                       "#cab2d6"),
                       # category.color=c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a", 
                       #                  "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99", "#b15928", "#8dd3c7", "#fb8072"),
                       # category.fill=c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a", 
                       #                 "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99", "#b15928", "#8dd3c7", "#fb8072"),
                       myposition="stack")
ggsave(file.path(outpath,"microhomology.smoking.tp53.egfr.kras.stack.pdf"),width = 3,height = 4)
