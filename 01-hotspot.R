#### library
library(dplyr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(ggpubr)

#### function
fragile.site.function <- function(x){
  chr1 <- as.character(x["CHROM1"])
  pos1 <- as.numeric(x["POS1"])
  chr2 <- as.character(x["CHROM2"])
  pos2 <- as.numeric(x["POS2"])
  t1 <- sum(fragile_site$chrom==chr1 & fragile_site$start<= pos1 & fragile_site$end >= pos1)>=1
  t2 <- sum(fragile_site$chrom==chr2 & fragile_site$start<= pos2 & fragile_site$end >= pos2)>=1
  if (sum(t1&t2)>=1){
    return("Fragile")
  } else {
    return("Not fragile")
  }
}


#### path
mychr <- paste0("chr",c(seq(1,22,1),"X","Y"))
options(scipen = 10)
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
hq.path <- file.path(sherlock.path,"/Data/SAMPLE_INFO/HQ_samples.csv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
sv.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
centromeres.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/Data/REF/hg38/centromeres/hg38.centromeres.bed"
chrlength.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/Data/REF/hg38/genome_length/genome_length_hg38.tsv"
fragile.path <- file.path("/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/OpenPBTA/scratch","fragile_site_hg38.csv")
plot.path <- file.path(scratch.path,"13-hotspot","plot")
out.path <- file.path(scratch.path,"13-hotspot")
cosmic.cancer.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/Data/COSMIC/cancergene/Census_allTue Sep  6 22_27_43 2022.csv"

if (dir.exists(plot.path)==FALSE){
  dir.create(plot.path,showWarnings = F,recursive = T)
}


#### read
bin = 1000000
sv <- read.delim(sv.path)
sampleinfo <-  read.delim(sampleinfo.path)
histology <- read.csv(histology.path)
sampleinfo <- merge(sampleinfo,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")
hqsample <- read.csv(hq.path) %>% pull(Tumor_Barcode)
tp53_mutated_samples=read.delim(tp53.mut.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
kras_mutated_samples=read.delim(kras.mut.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F) %>% pull(V1)
centromeres <- read.delim(centromeres.path,header = F)
names(centromeres) <- c("Chromosome","Pos1","Pos2")
centromeres$Chromosome <- factor(centromeres$Chromosome,levels = mychr)
chrlength <- read.delim(chrlength.path)
names(chrlength) <- c("Chromosome","Length")
chrlength$Chromosome <- paste0("chr",chrlength$Chromosome)
chrlength$Chromosome <- factor(chrlength$Chromosome,levels = mychr)
fragile_site <- read.csv(fragile.path,fileEncoding = "UTF-8-BOM")
cosmic.cancer <- read.csv(cosmic.cancer.path)
sbs4 <- read.delim2(sbs4.path)

####### add information
svmerge <- sv
# add fragile
svmerge$fragile <- apply(svmerge, 1, fragile.site.function)
# add histology etc.
sampleinfo <- merge(sampleinfo,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
svmerge <- merge(svmerge,sampleinfo[,c("Tumor_Barcode","Assigned_Population","Gender","Smoking","SP_Group",
                                       "Histology","purity","wgd","scna_group","mut_load","indel_load","SBS4")],
                 by.x = "SAMPLE",by.y = "Tumor_Barcode",all.x = T) %>% 
  filter(SAMPLE %in% sampleinfo$Tumor_Barcode)
# # add sv size
# svmerge$SIZE <- ifelse(svmerge$SVTYPE=="TRA",0,abs(svmerge$POS1-svmerge$POS2))
# cancer gene
cosmic.cancer$chr <- sapply(strsplit(cosmic.cancer$Genome.Location, ":"), function(x) paste0("chr",x[1]))
cosmic.cancer$pos1 <- regmatches(cosmic.cancer$Genome.Location, regexpr(":(.*?)-", cosmic.cancer$Genome.Location))
cosmic.cancer$pos1 <- substr(cosmic.cancer$pos1, 2, nchar(cosmic.cancer$pos1) - 1)
cosmic.cancer$pos2 <- sapply(regmatches(cosmic.cancer$Genome.Location, regexpr(":(.*?)-", cosmic.cancer$Genome.Location), invert = TRUE),function(x) x[2])
cosmic.cancer <- cosmic.cancer %>% filter(pos1 != "") %>% 
  mutate(pos1=as.numeric(pos1),pos2=as.numeric(pos2)) %>%
  mutate(pos1.bin=ceiling(pos1/bin),pos2.bin=ceiling(pos2/bin))

####### SV to breakpoint
# breakpoint level
breakpoint <- rbind(svmerge %>% dplyr::select(SAMPLE,CHROM1,POS1,SVTYPE,STRAND1,Signature,"Assigned_Population","Gender","Smoking","SP_Group",
                                              "Histology","SBS4",purity,wgd,scna_group,mut_load,indel_load,fragile) %>%
                      setnames(.,old = c("CHROM1","POS1","STRAND1"),new = c("Chromosome","POS","STRAND")),
                    svmerge %>% dplyr::select(SAMPLE,CHROM2,POS2,SVTYPE,STRAND2,Signature,"Assigned_Population","Gender","Smoking","SP_Group",
                                              "Histology","SBS4",purity,wgd,scna_group,mut_load,indel_load,fragile) %>%
                      setnames(.,old = c("CHROM2","POS2","STRAND2"),new = c("Chromosome","POS","STRAND")))
# add bin rank
breakpoint$bin <- ceiling(breakpoint$POS / bin)
# add EGFR mut
breakpoint$EGFR <- ifelse(breakpoint$SAMPLE %in% c(egfr_mutated_E479_A483del_samples,
                                                   egfr_mutated_L591R_samples,
                                                   egfr.others.path),"Mut","WT")
# add Smoking_SBS4
breakpoint <- breakpoint %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
sampleinfo <- sampleinfo %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))

####### filter
breakpoint <- breakpoint %>% filter(Smoking_SBS4 %in% c("Smoker-SBS4","Non-Smoker-noSBS4"))
####### plot
# my.category <- "EGFR"
# my.order <- c("Mut","WT")
# my.color <- c("red","grey60")

my.category <- "Smoking"
my.order <- c("Non-Smoker","Smoker","Unknown")
my.color <- c("#b3b3b3", "black", "#F0F0F0")
my.color <- c("black", "black", "black")

breakpoint$my.category <- breakpoint[,my.category]

# prepare chromosome
chr.df <- data.frame("Chromosome"=mychr)
chr.df$Chromosome <- factor(chr.df$Chromosome,levels = mychr)
# prepare plot df
plot.breakpoint.df <- breakpoint %>%
  group_by(my.category,SAMPLE,Chromosome,Signature,bin) %>% 
  count() %>% 
  mutate(Chromosome=factor(paste0("chr",Chromosome),levels = mychr),
         Signature=factor(Signature,levels = c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss",
                                               "4 Micronuclei","5 Large gain","6 Hourglass","un_assigned",
                                               "chromoplexy","cycle_templated_ins","complex_unclear",
                                               "Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra"
         ))) %>%
  group_by(my.category,Chromosome,Signature,bin) %>% mutate(binsum.at.category=sum(n),n2=1) %>%
  group_by(Chromosome,Signature,bin) %>% mutate(binsum.at.alltumor=sum(n)) %>%
  # group_by(my.category,Chromosome,Signature,bin) %>% mutate(binsample.at.tumor.type=sum(n2)) %>% 
  group_by(Chromosome,Signature,bin) %>% mutate(freq=sum(n2)) %>% 
  # arrange(desc(binsum.at.category),desc(bin))
  arrange(Signature,Chromosome,bin,desc(bin),desc(my.category))

plot.breakpoint.df_smoker <- breakpoint %>%
  filter(Smoking_SBS4 %in% c("Smoker-SBS4")) %>%
  group_by(my.category,SAMPLE,Chromosome,Signature,bin) %>% 
  count() %>% 
  mutate(Chromosome=factor(paste0("chr",Chromosome),levels = mychr),
         Signature=factor(Signature,levels = c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss",
                                               "4 Micronuclei","5 Large gain","6 Hourglass","un_assigned",
                                               "chromoplexy","cycle_templated_ins","complex_unclear",
                                               "Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra"
         ))) %>%
  group_by(my.category,Chromosome,Signature,bin) %>% mutate(binsum.at.category=sum(n),n2=1) %>%
  group_by(Chromosome,Signature,bin) %>% mutate(binsum.at.alltumor=sum(n)) %>%
  # group_by(my.category,Chromosome,Signature,bin) %>% mutate(binsample.at.tumor.type=sum(n2)) %>% 
  group_by(Chromosome,Signature,bin) %>% mutate(freq=sum(n2)) %>% 
  # arrange(desc(binsum.at.category),desc(bin))
  arrange(Signature,Chromosome,bin,desc(bin),desc(my.category))
plot.breakpoint.df_nonsmoker <- breakpoint %>%
  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4")) %>%
  group_by(my.category,SAMPLE,Chromosome,Signature,bin) %>% 
  count() %>% 
  mutate(Chromosome=factor(paste0("chr",Chromosome),levels = mychr),
         Signature=factor(Signature,levels = c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss",
                                               "4 Micronuclei","5 Large gain","6 Hourglass","un_assigned",
                                               "chromoplexy","cycle_templated_ins","complex_unclear",
                                               "Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra"
         ))) %>%
  group_by(my.category,Chromosome,Signature,bin) %>% mutate(binsum.at.category=sum(n),n2=1) %>%
  group_by(Chromosome,Signature,bin) %>% mutate(binsum.at.alltumor=sum(n)) %>%
  # group_by(my.category,Chromosome,Signature,bin) %>% mutate(binsample.at.tumor.type=sum(n2)) %>% 
  group_by(Chromosome,Signature,bin) %>% mutate(freq=sum(n2)) %>% 
  # arrange(desc(binsum.at.category),desc(bin))
  arrange(Signature,Chromosome,bin,desc(bin),desc(my.category))
smoker.n <- sampleinfo %>% 
  filter(Smoking_SBS4 %in% c("Smoker-SBS4")) %>%
  count() %>% as.numeric()
nonsmoker.n <- sampleinfo %>% 
  filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4")) %>%
  count() %>% as.numeric()
plot.breakpoint.df <- rbind(plot.breakpoint.df_nonsmoker %>% mutate(my.category="nonsmoker"),
                            plot.breakpoint.df_smoker %>% mutate(n2=-n2,my.category="smoker"))

plot.breakpoint.df.pct <- plot.breakpoint.df %>% dplyr::select(my.category,Chromosome,bin,freq) %>%
  unique() %>%
  mutate(pct=ifelse(my.category=="nonsmoker",freq/nonsmoker.n,-freq/smoker.n))


my.order <- c("nonsmoker","smoker")
my.color <- c("#00a800","#FF6EC7")

# add cancer gene
plot.breakpoint.df$cancer.gene <- apply(plot.breakpoint.df,1,function(x){
  chri=x["Chromosome"]
  bini=x["bin"] %>% as.numeric()
  tell=grep(TRUE,cosmic.cancer$chr==chri & cosmic.cancer$pos1.bin <= bini & cosmic.cancer$pos2.bin >= bini)
  gene=paste(cosmic.cancer$Gene.Symbol[tell],collapse = ";")
  return(gene)
})
plot.breakpoint.df.unique <- plot.breakpoint.df %>% 
  select(Signature,Chromosome,bin,freq,cancer.gene) %>% 
  unique() %>% filter(cancer.gene != "")
  
# plot
hotspot.bar <- ggplot() + 
  geom_rect(data=chrlength,
            aes(fill = Chromosome,
                               xmin = 1,xmax = ceiling(Length/bin),
                               ymin = -Inf,ymax = Inf),
            alpha = 0.3) +
  geom_bar(data=plot.breakpoint.df %>% filter(Signature != "un_assigned"),
           aes(x=bin,y=n2,fill=my.category,color=my.category),
           stat = "identity",size=0.01)+
  scale_color_manual(name=my.category,
                     breaks = my.order,
                     values = my.color,
                     drop=FALSE
  )+
  scale_fill_manual(name="",
                    breaks = c(mychr,my.order),
                    values = c(rep(c("white","grey90"),12),my.color),
                    drop=FALSE
  )+
  #geom_hline(data=centromeres,aes(yintercept = 3),alpha=0.5,linetype="dashed",color="grey80",size=0.1)+
  geom_vline(data=centromeres,aes(xintercept = (ceiling(Pos1/bin) + ceiling(Pos2/bin))/2),alpha=0.5,linetype="dashed",color="grey80",size=0.1)+
  geom_vline(data=chrlength,aes(xintercept = 0),alpha=0.5,linetype="solid",color="white",size=0.1)+
  geom_vline(data=chrlength,aes(xintercept = ceiling(Length/bin)),alpha=0.5,linetype="solid",color="white",size=0.1)+
  scale_x_continuous(expand = c(0,0))+
  facet_grid2(Signature~Chromosome,
              scales = "free",
              space = "free_x",
              switch = "both",
              strip = strip_vanilla(clip = "off"))+
  guides(fill = guide_legend(nrow = 4),
         colour = "none")+
  theme_bw()+
  ylab("Number of sample")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.1,colour = "black"), 
        axis.line.x.top = element_line(size=0.1,colour = "black"), 
        axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
        axis.ticks.x =  element_blank(),
        axis.ticks.length.y =unit(0.8, "mm"),
        axis.ticks.y =  element_line(linewidth = 0.234),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 0.5,vjust = 0.5),
        strip.text.y.left = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0),
        strip.background =  element_blank(),
        strip.placement = "top",
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(2,"mm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        legend.position = "bottom",
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )

pct.bar <- ggplot() + 
  geom_rect(data=chrlength,
            aes(fill = Chromosome,
                xmin = 1,xmax = ceiling(Length/bin),
                ymin = -Inf,ymax = Inf),
            alpha = 0.3) +
  geom_bar(data=plot.breakpoint.df.pct %>% filter(Signature != "un_assigned"),
           aes(x=bin,y=pct*100,fill=my.category,color=my.category),
           stat = "identity",size=0.01)+
  # scale_y_continuous(labels = scales::percent) +
  scale_color_manual(name=my.category,
                     breaks = my.order,
                     values = my.color,
                     drop=FALSE
  )+
  scale_fill_manual(name="",
                    breaks = c(mychr,my.order),
                    values = c(rep(c("white","grey90"),12),my.color),
                    drop=FALSE
  )+
  #geom_hline(data=centromeres,aes(yintercept = 3),alpha=0.5,linetype="dashed",color="grey80",size=0.1)+
  geom_vline(data=centromeres,aes(xintercept = (ceiling(Pos1/bin) + ceiling(Pos2/bin))/2),alpha=0.5,linetype="dashed",color="grey80",size=0.1)+
  geom_vline(data=chrlength,aes(xintercept = 0),alpha=0.5,linetype="solid",color="white",size=0.1)+
  geom_vline(data=chrlength,aes(xintercept = ceiling(Length/bin)),alpha=0.5,linetype="solid",color="white",size=0.1)+
  scale_x_continuous(expand = c(0,0))+
  facet_grid2(Signature~Chromosome,
              scales = "free",
              space = "free_x",
              switch = "both",
              strip = strip_vanilla(clip = "off"))+
  guides(fill = guide_legend(nrow = 4),
         colour = "none")+
  theme_bw()+
  ylab("Number of sample")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.1,colour = "black"), 
        axis.line.x.top = element_line(size=0.1,colour = "black"), 
        axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
        axis.ticks.x =  element_blank(),
        axis.ticks.length.y =unit(0.8, "mm"),
        axis.ticks.y =  element_line(linewidth = 0.234),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 0.5,vjust = 0.5),
        strip.text.y.left = element_text(size=6, angle = 0, colour = "black",hjust = 1,vjust = 0),
        strip.background =  element_blank(),
        strip.placement = "top",
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(2,"mm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        legend.position = "bottom",
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )

chr.bar <- ggplot() + geom_rect(data=chrlength,aes(xmin = 1,xmax = ceiling(Length/bin),ymin = 0,ymax = 1),
                                alpha = 0.3,fill="grey20") +
  geom_rect(data=centromeres,aes(xmin = ceiling(Pos1/bin),xmax = ceiling(Pos2/bin)+1,ymin = 0,ymax = 1),
            alpha = 1,fill="red") +
  facet_grid2(.~Chromosome,
              scales = "free",
              space = "free_x",
              switch = "both",
              strip = strip_vanilla(clip = "off"))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.ticks.x =  element_blank(),
        axis.ticks.length.y =unit(0, "mm"),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y.left = element_blank(),
        strip.background =  element_blank(),
        strip.placement = "top",
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(2,"mm"),
        legend.title = element_text(size=6),
        legend.text = element_text(size=6),
        legend.position = "bottom",
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )

# ggarrange(hotspot.bar,chr.bar,ncol = 1,heights = c(5,0.1))
# ggsave(file.path(plot.path,
#                  paste0("hotspot.across.chr.allsample.popsmoking",bin,".pdf")),width = 170,height = 200,units = c("mm")) 
# ggsave(file.path(plot.path,
#                  paste0("hotspot.across.chr.allsample.popsmoking.smoker",bin,".pdf")),width = 170,height = 200,units = c("mm")) 
# ggsave(file.path(plot.path,
#                  paste0("hotspot.across.chr.allsample.popsmoking.nonsmoker",bin,".pdf")),width = 170,height = 200,units = c("mm")) 
# ggsave(file.path(plot.path,
#                  paste0("hotspot.across.chr.allsample.egfr",bin,".pdf")),width = 170,height = 200,units = c("mm")) 

ggarrange(pct.bar,chr.bar,ncol = 1,heights = c(5,0.1))
ggsave(file.path(plot.path,
                 paste0("hotspot.across.chr.allsample.popsmoking.advanced.sbs4consistent",bin,".pdf")),width = 170,height = 200,units = c("mm")) 
ggsave(file.path(plot.path,
                 paste0("hotspot.across.chr.allsample.popsmoking.advanced.sbs4inconsistent",bin,".pdf")),width = 170,height = 200,units = c("mm")) 

# "1 ecDNA/double minutes"
# "2 BFB cycles/chromatin bridge"
# "3 Large loss"
# "4 Micronuclei"
# "5 Large gain"
# "6 Hourglass"
# "cycle_templated_ins"
# "chromoplexy"
# "complex_unclear"
# "Del1"
# "Del2"
# "Del3"
# "TD1"
# "TD2"
# "Large intra"
# "Fb inv"
# "Tra"

plot.breakpoint.df.pct.caner <- merge(plot.breakpoint.df.pct,cosmic.cancer[,c("chr","pos1.bin","Gene.Symbol")],
                                      by.x = c("Chromosome","bin"),by.y = c("chr","pos1.bin"),all.x = T)
findsig=c("Tra")
findchr="X"
x0=plot.breakpoint.df.pct.caner %>% filter(Signature %in% findsig) %>% arrange(desc(pct))
x=plot.breakpoint.df.pct.caner %>% filter(Signature %in% findsig & Chromosome==paste0("chr",findchr)) %>% arrange(desc(pct))
y=plot.breakpoint.df %>% filter(Signature %in% findsig & Chromosome==paste0("chr",findchr) & bin %in% c(60))
z=breakpoint %>% filter(Signature %in% findsig & Chromosome==findchr & bin %in% c(8)) 

x=breakpoint %>% filter(Signature %in% findsig & Chromosome==findchr) %>% group_by(SAMPLE,bin) %>% count() %>% group_by(bin) %>% count() %>% arrange(desc(n))
y=plot.breakpoint.df %>% filter(Signature %in% findsig & Chromosome==paste0("chr",findchr) & bin %in% c(147))
z=breakpoint %>% filter(Signature %in% findsig & Chromosome==findchr & bin %in% c(147)) 
z1=z %>% select(SAMPLE,Signature,POS) %>% group_by(SAMPLE,Signature) %>% count() %>% group_by(SAMPLE) %>% count()
z1=sv %>% filter((CHROM2==findchr|CHROM2==findchr) & Signature %in% findsig)
unique(c(sv %>% filter(SAMPLE == "NSLC-0201-T01" & Signature %in% findsig) %>% pull(CHROM1),
sv %>% filter(SAMPLE == "NSLC-0201-T01" & Signature %in% findsig) %>% pull(CHROM2)))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr9" &Signature=="Del3")
y=plot.breakpoint.df %>% filter(Chromosome=="chr9" &Signature=="Del3" & bin %in% c(10))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr9" &Signature=="3 Large loss")
y=plot.breakpoint.df %>% filter(Chromosome=="chr9" &Signature=="3 Large loss" & bin %in% c(22))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr20" &Signature=="6 Hourglass")
y=plot.breakpoint.df %>% filter(Chromosome=="chr20" &Signature=="6 Hourglass" & bin %in% c(39))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr7" &Signature=="Del1")
y=plot.breakpoint.df %>% filter(Chromosome=="chr7" &Signature=="Del1" & bin %in% c(2))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chrX" &Signature=="Del1")
y=plot.breakpoint.df %>% filter(Chromosome=="chrX" &Signature=="Del1" & bin %in% c(2,90,91,92))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr22" &Signature=="Tra")
y=plot.breakpoint.df %>% filter(Chromosome=="chr22" &Signature=="Tra" & bin %in% c(29,13))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr20" &Signature=="Del3")
x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr2" &Signature=="TD1")

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr11" &Signature=="TD2")
y=plot.breakpoint.df %>% filter(Chromosome=="chr11" &Signature=="TD2" & bin %in% c(66))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr1" &Signature=="3 Large loss")
y=plot.breakpoint.df %>% filter(Chromosome=="chr1" &Signature=="3 Large loss" & bin %in% c(25,30))

x=plot.breakpoint.df.pct %>% filter(Chromosome=="chrX" &Signature=="Tra")

z=sv %>% filter(((CHROM1=="22" & POS1 >= 28000000 & POS1 <= 29000000)| (CHROM2=="22" & POS2 >= 28000000 & POS2 <= 29000000))& Signature=="Tra")
z=sv %>% filter(((CHROM1=="22" & POS1 >= 12000000 & POS1 <= 13000000)| (CHROM2=="22" & POS2 >= 12000000 & POS2 <= 13000000))& Signature=="Tra")
z=sv %>% filter(((CHROM1=="X" & POS1 >= 11000000 & POS1 <= 12000000)| (CHROM2=="X" & POS2 >= 11000000 & POS2 <= 12000000))& Signature=="Tra")
x=plot.breakpoint.df.pct %>% filter(Chromosome=="chr22" &Signature=="Tra")
y=plot.breakpoint.df %>% filter(Chromosome=="chr21" &Signature=="Tra" & bin %in% c(11))

z1=sv %>% filter(((CHROM1=="22" & POS1 >= 12000000 & POS1 <= 13000000)| (CHROM2=="22" & POS2 >= 12000000 & POS2 <= 13000000))& Signature=="Tra")
z2=sv %>% filter(((CHROM1=="20" & POS1 >= 29000000 & POS1 <= 30000000)| (CHROM2=="20" & POS2 >= 29000000 & POS2 <= 30000000))& Signature=="Tra")
z3=sv %>% filter(((CHROM1=="21" & POS1 >= 10000000 & POS1 <= 11000000)| (CHROM2=="21" & POS2 >= 10000000 & POS2 <= 11000000))& Signature=="Tra")
z1.s=z1 %>% pull(SAMPLE) %>% unique()
z2.s=z2 %>% pull(SAMPLE) %>% unique()
z3.s=z3 %>% pull(SAMPLE) %>% unique()

common_elements <- intersect(z1.s, z2.s)
percentage_overlap <- (length(common_elements) / min(length(z1.s), length(z2.s))) * 100
percentage_overlap

common_elements <- intersect(z1.s, z3.s)
percentage_overlap <- (length(common_elements) / min(length(z1.s), length(z3.s))) * 100
percentage_overlap

common_elements <- intersect(z2.s, z3.s)
percentage_overlap <- (length(common_elements) / min(length(z2.s), length(z3.s))) * 100
percentage_overlap
