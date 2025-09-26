setwd("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock")

########################### library ##########################
library(ggplot2)
library(reshape)
library(scales)
library(dplyr)

########################### save path #########################
out_path <- file.path("D:","OneDrive","OneDrive - The University of Chicago","office","GitHub","Sherlock","scratch","PCAWG_LUAD","Battenberg_and_sv_50bp")
if (!dir.exists(out_path)){
  dir.create(out_path)
}

########################### file path #########################
all_sv_path <-
  file.path(
    "D:",
    "OneDrive",
    "OneDrive - The University of Chicago",
    "office",
    "GitHub",
    "Sherlock",
    "scratch",
    "PCAWG_LUAD",
    "compare_manta.meerkat.tnscope.50bp",
    "tables",
    "all",
    "PCAWG_SV_by_orginal_manta_meerkat_tnscope.csv"
  )

bp_pcawg_bb_path <-  
  file.path(
    "D:",
    "OneDrive",
    "OneDrive - The University of Chicago",
    "office",
    "GitHub",
    "Sherlock",
    "scratch",
    "PCAWG_LUAD",
    "Battenberg",
    "BP_PCAWG_BBsegments.txt"
  )

bp_sherlock_bb_path <-  
  file.path(
    "D:",
    "OneDrive",
    "OneDrive - The University of Chicago",
    "office",
    "GitHub",
    "Sherlock",
    "scratch",
    "PCAWG_LUAD",
    "Battenberg",
    "BP_Sherlock_BBsegments.txt"
  )


########################### read file #########################
all_sv <- read.csv(all_sv_path)
bp_pcawg_bb <- read.delim(bp_pcawg_bb_path)
bp_sherlock_bb <- read.delim(bp_sherlock_bb_path)
pcawg_sv <- all_sv[all_sv$pcawg=="1",]
sherlock_sv <- all_sv[all_sv$manta=="1" | all_sv$meerkat=="1",]
tnscope_sv <- all_sv[all_sv$tnscope=="1",]


########################## calculate distance ######################
# sv <- pcawg_sv
# bp <- bp_pcawg_bb

distance_function <- function(sv,bb){
  for (i in 1:nrow(sv)) {
    chr1 <- as.character(sv$CHROM[i])
    pos1 <- sv$POS[i]
    chr2 <- as.character(sv$CHROM2[i])
    pos2 <- sv$POS2[i]
    samplei <- as.character(sv$sample[i])
    distance1 <- min(abs(as.numeric(unlist(bb[bb$Tumor_Barcode==samplei & bb$chr==chr1,c("startpos","endpos")])) - pos1))
    distance2 <- min(abs(as.numeric(unlist(bb[bb$Tumor_Barcode==samplei & bb$chr==chr2,c("startpos","endpos")])) - pos2))
    sv$distance1[i] <- distance1
    sv$distance2[i] <- distance2
  }
  return(sv)
}

pcawg_sv <- distance_function(pcawg_sv,bp_pcawg_bb)
sherlock_sv <- distance_function(sherlock_sv,bp_sherlock_bb)
tnscope_sv <- distance_function(tnscope_sv,bp_sherlock_bb)
# add size
pcawg_sv$size <- ifelse(pcawg_sv$SVTYPE != "TRA",abs(pcawg_sv$POS-pcawg_sv$POS2),0)
sherlock_sv$size <- ifelse(sherlock_sv$SVTYPE != "TRA",abs(sherlock_sv$POS-sherlock_sv$POS2),0)
tnscope_sv$size <- ifelse(tnscope_sv$SVTYPE != "TRA",abs(tnscope_sv$POS-tnscope_sv$POS2),0)
# add size group
sizecutoff <- 10000
pcawg_sv$sizegroup <- ifelse(pcawg_sv$size>0 & pcawg_sv$size<sizecutoff,"small","large")
sherlock_sv$sizegroup <- ifelse(sherlock_sv$size>0 & sherlock_sv$size<sizecutoff,"small","large")
tnscope_sv$sizegroup <- ifelse(tnscope_sv$size>0 & tnscope_sv$size<sizecutoff,"small","large")
# turn inf to max
inf_value <- 1000000000
pcawg_sv$distance1 <- ifelse(is.infinite(pcawg_sv$distance1),inf_value,pcawg_sv$distance1)
pcawg_sv$distance2 <- ifelse(is.infinite(pcawg_sv$distance2),inf_value,pcawg_sv$distance2)
sherlock_sv$distance1 <- ifelse(is.infinite(sherlock_sv$distance1),inf_value,sherlock_sv$distance1)
sherlock_sv$distance2 <- ifelse(is.infinite(sherlock_sv$distance2),inf_value,sherlock_sv$distance2)
tnscope_sv$distance1 <- ifelse(is.infinite(tnscope_sv$distance1),inf_value,tnscope_sv$distance1)
tnscope_sv$distance2 <- ifelse(is.infinite(tnscope_sv$distance2),inf_value,tnscope_sv$distance2)

####################### which is supported (bar plot) #######################
# both distances less than 5000
cutoff <- 10000

# add support
sherlock_sv$support <- ifelse(sherlock_sv$distance1<cutoff & sherlock_sv$distance2<cutoff,"support","no")
pcawg_sv$support <- ifelse(pcawg_sv$distance1<cutoff & pcawg_sv$distance2<cutoff,"support","no")
tnscope_sv$support <- ifelse(tnscope_sv$distance1<cutoff & tnscope_sv$distance2<cutoff,"support","no")

# PCAWG
pcawg_sv$manta_meerkat_union <- pcawg_sv$manta | pcawg_sv$meerkat
pcawg_sv_plot <- rbind(
  "support_overlap"=data.frame("support"="support","calling"="overlap","value"=sum(pcawg_sv$support=="support" & pcawg_sv$manta_meerkat_union=="TRUE")),
  "support_specific"=data.frame("support"="support","calling"="specific","value"=sum(pcawg_sv$support=="support" & pcawg_sv$manta_meerkat_union=="FALSE")),
  "unsupport_overlap"=data.frame("support"="unsupport","calling"="overlap","value"=sum(pcawg_sv$support=="no" & pcawg_sv$manta_meerkat_union=="TRUE")),
  "unsupport_specific"=data.frame("support"="unsupport","calling"="specific","value"=sum(pcawg_sv$support=="no" & pcawg_sv$manta_meerkat_union=="FALSE"))
)
pcawg_sv_plot$support <- factor(pcawg_sv_plot$support,levels=c("unsupport","support"))
pcawg_sv_plot$calling <- factor(pcawg_sv_plot$calling,levels = c("specific","overlap"))
specitic_unsupport_ratio <- percent(pcawg_sv_plot$value[4]/(pcawg_sv_plot$value[4]+pcawg_sv_plot$value[2]))
specitic_support_ratio <- percent(pcawg_sv_plot$value[2]/(pcawg_sv_plot$value[2]+pcawg_sv_plot$value[4]))
overlap_unsupport_ratio <- percent(pcawg_sv_plot$value[3]/(pcawg_sv_plot$value[3]+pcawg_sv_plot$value[1]))
overlap_support_ratio <- percent(pcawg_sv_plot$value[1]/(pcawg_sv_plot$value[1]+pcawg_sv_plot$value[3]))

ggplot(data=pcawg_sv_plot,aes(x=calling,
                              y=value,
                              alpha=support,
                              fill=calling))+ 
  geom_col(position="stack",width = 0.8)+
  annotate(geom="text", x=2, y=pcawg_sv_plot$value[1]+200, label=overlap_unsupport_ratio,color="blue",size=2)+
  annotate(geom="text", x=1, y=pcawg_sv_plot$value[2]+200, label=specitic_unsupport_ratio,color="blue",size=2)+
  annotate(geom="text", x=2, y=pcawg_sv_plot$value[1]-200, label=overlap_support_ratio,color="blue",size=2)+
  annotate(geom="text", x=1, y=pcawg_sv_plot$value[2]-200, label=specitic_support_ratio,color="blue",size=2)+
  scale_fill_manual(values=c('#CB4335', "#F2AE10"), 
                    breaks = c("overlap","specific"),
                    labels = c("Overlap","PCAWG-specific"))+
  scale_x_discrete(breaks = c("overlap","specific"),
                   labels = c("Overlap","PCAWG-specific"))+
  scale_alpha_manual(values=c(0.5,1), 
                     breaks = c("unsupport","support"))+
  ylab("SV Number")+
  ggtitle(paste0("PCAWG\n","cutoff is ",cutoff))+
  ylim(0,4500)+
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.text.y = element_text(angle=0,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.1,colour = "black"), 
        axis.line.x.top = element_line(size=0.1,colour = "black"), 
        axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
        axis.ticks.length.x =unit(0.8, "mm"),
        axis.ticks.x =  element_line(linewidth = 0.234),
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
        title =  element_text(size=6,colour = "black"),
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )

ggsave(file.path(out_path,"pcawg.sv.supported.by.cnv.1kb.pdf"),width = 1.5,height = 2.25)
ggsave(file.path(out_path,"pcawg.sv.supported.by.cnv.5kb.pdf"),width = 1.5,height = 2.25)
ggsave(file.path(out_path,"pcawg.sv.supported.by.cnv.10kb.pdf"),width = 1.5,height = 2.25)


ggplot(data=pcawg_sv_plot %>% group_by(calling) %>% mutate(pct=value/sum(value)) %>% filter(support=="support"))+ 
  geom_col(aes(x=calling,y=pct,alpha=support,fill=calling))+
  geom_text(aes(x=calling,y=pct,label=paste0(round(pct,2)*100,"%")),hjust=0.5,vjust=-1,size=2)+
  scale_fill_manual(values=c('#CB4335', "#F2AE10"), 
                    breaks = c("overlap","specific"),
                    labels = c("Overlap","PCAWG-specific"))+
  scale_x_discrete(breaks = c("overlap","specific"),
                   labels = c("Overlap","PCAWG-specific"))+
  scale_alpha_manual(values=c(0.5,1), 
                     breaks = c("unsupport","support"))+
  scale_y_continuous(limits = c(0,1))+
  ylab("Percentage of SVs supported by CNV")+
  ggtitle(paste0("PCAWG\n","cutoff is ",cutoff))+
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.text.y = element_text(angle=0,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.1,colour = "black"), 
        axis.line.x.top = element_line(size=0.1,colour = "black"), 
        axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
        axis.ticks.length.x =unit(0.8, "mm"),
        axis.ticks.x =  element_line(linewidth = 0.234),
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
        title =  element_text(size=6,colour = "black"),
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(out_path,"pcawg.sv.supported.by.cnv.pct.1kb.pdf"),width = 1.5,height = 2.25)
ggsave(file.path(out_path,"pcawg.sv.supported.by.cnv.pct.5kb.pdf"),width = 1.5,height = 2.25)
ggsave(file.path(out_path,"pcawg.sv.supported.by.cnv.pct.10kb.pdf"),width = 1.5,height = 2.25)

# SHERLOCK
sherlock_sv_plot <- rbind(
  "support_overlap"=data.frame("support"="support","calling"="overlap","value"=sum(sherlock_sv$support=="support" & sherlock_sv$pcawg=="1")),
  "support_specific"=data.frame("support"="support","calling"="specific","value"=sum(sherlock_sv$support=="support" & sherlock_sv$pcawg=="0")),
  "unsupport_overlap"=data.frame("support"="unsupport","calling"="overlap","value"=sum(sherlock_sv$support=="no" & sherlock_sv$pcawg=="1")),
  "unsupport_specific"=data.frame("support"="unsupport","calling"="specific","value"=sum(sherlock_sv$support=="no" & sherlock_sv$pcawg=="0"))
)
sherlock_sv_plot$support <- factor(sherlock_sv_plot$support,levels=c("unsupport","support"))
sherlock_sv_plot$calling <- factor(sherlock_sv_plot$calling,levels = c("specific","overlap"))
overlap_unsupport_ratio <- percent(sherlock_sv_plot$value[3]/(sherlock_sv_plot$value[3]+sherlock_sv_plot$value[1]))
overlap_support_ratio <- percent(sherlock_sv_plot$value[1]/(sherlock_sv_plot$value[3]+sherlock_sv_plot$value[1]))
specitic_unsupport_ratio <- percent(sherlock_sv_plot$value[4]/(sherlock_sv_plot$value[4]+sherlock_sv_plot$value[2]))
specitic_support_ratio <- percent(sherlock_sv_plot$value[2]/(sherlock_sv_plot$value[4]+sherlock_sv_plot$value[2]))

ggplot(data=sherlock_sv_plot,aes(x=calling,
                                 y=value,
                                 alpha=support,
                                 fill=calling))+ 
  geom_col(position="stack",width = 0.8)+
  annotate(geom="text", x=2, y=sherlock_sv_plot$value[1]+200, label=overlap_unsupport_ratio,color="blue",size=2)+
  annotate(geom="text", x=1, y=sherlock_sv_plot$value[2]+200, label=specitic_unsupport_ratio,color="blue",size=2)+
  annotate(geom="text", x=2, y=sherlock_sv_plot$value[1]-200, label=overlap_support_ratio,color="blue",size=2)+
  annotate(geom="text", x=1, y=sherlock_sv_plot$value[2]-200, label=specitic_support_ratio,color="blue",size=2)+
  scale_fill_manual(values=c('#CB4335', "#3498DB"), 
                    breaks = c("overlap","specific"),
                    labels = c("Overlap","specific"))+
  scale_alpha_manual(values=c(0.5,1), 
                     breaks = c("unsupport","support"),
                     labels=c("no supported CNV","CNV support")
  )+
  scale_x_discrete(breaks = c("overlap","specific"),
                   labels = c("Overlap","Manta/Meerkat-specific"))+
  ylab("SV Number") +
  ggtitle(paste0("Manta+Meerkat Union\n","cutoff is ",cutoff))+
  ylim(0,5500)+
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.text.y = element_text(angle=0,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.1,colour = "black"), 
        axis.line.x.top = element_line(size=0.1,colour = "black"), 
        axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
        axis.ticks.length.x =unit(0.8, "mm"),
        axis.ticks.x =  element_line(linewidth = 0.234),
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
        title =  element_text(size=6,colour = "black"),
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(out_path,"sherlock.sv.supported.by.cnv.1kb.pdf"),width = 1.5,height = 2.5)
ggsave(file.path(out_path,"sherlock.sv.supported.by.cnv.5kb.pdf"),width = 1.5,height = 2.5)
ggsave(file.path(out_path,"sherlock.sv.supported.by.cnv.10kb.pdf"),width = 1.5,height = 2.5)

ggplot(data=sherlock_sv_plot %>% group_by(calling) %>% mutate(pct=value/sum(value)) %>% filter(support=="support"))+ 
  geom_col(aes(x=calling,y=pct,alpha=support,fill=calling))+
  geom_text(aes(x=calling,y=pct,label=paste0(round(pct,2)*100,"%")),hjust=0.5,vjust=-1,size=2)+
  scale_fill_manual(values=c('#CB4335', "#3498DB"), 
                    breaks = c("overlap","specific"),
                    labels = c("Overlap","specific"))+
  scale_alpha_manual(values=c(0.5,1), 
                     breaks = c("unsupport","support"),
                     labels=c("no supported CNV","CNV support")
  )+
  scale_x_discrete(breaks = c("overlap","specific"),
                   labels = c("Overlap","Manta/Meerkat-specific"))+
  scale_y_continuous(limits = c(0,1))+
  ylab("Percentage of SVs supported by CNV")+
  ggtitle(paste0("Manta+Meerkat Union\n","cutoff is ",cutoff))+
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.text.y = element_text(angle=0,hjust = 0,vjust = 0,size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_line(size=0.1,colour = "black"), 
        axis.line.x.top = element_line(size=0.1,colour = "black"), 
        axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
        axis.ticks.length.x =unit(0.8, "mm"),
        axis.ticks.x =  element_line(linewidth = 0.234),
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
        title =  element_text(size=6,colour = "black"),
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(out_path,"sherlock.sv.supported.by.cnv.pct.1kb.pdf"),width = 1.5,height = 2.5)
ggsave(file.path(out_path,"sherlock.sv.supported.by.cnv.pct.5kb.pdf"),width = 1.5,height = 2.5)
ggsave(file.path(out_path,"sherlock.sv.supported.by.cnv.pct.10kb.pdf"),width = 1.5,height = 2.5)
