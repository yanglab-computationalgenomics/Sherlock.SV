library(dplyr)
library(ggplot2)

sv_path <-
  file.path(
    "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock",
    "scratch",
    "01-SV-files",
    "1217_samples_hg38",
    "1217_manta_meerkat_union_window50.txt"
  )
sigprofile13_path <-
  file.path(
    "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV53_clusterasone/CH53/All_Solutions/SBS53_13_Signatures/Activities/SBS53_S13_NMF_Activities.txt"
  )
sample_info_path <-
  file.path(
    "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/Data/SAMPLE_INFO/sherlock_sampleinfo_merged.txt"
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
    "01-SV-files",
    "1217_samples_hg38",
    "BP_Sherlock_BBsegments.txt",
    "BP_Sherlock_BBsegments_nogaps.txt"
  )
rmsk_path <- file.path("D:/Data/REF/hg38/hg38_2020_rmsk.bed")
outpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/99-qc"


# read files
all_sv  <- read.delim(sv_path,stringsAsFactors = FALSE)
all_sv$size <- ifelse(all_sv$SVTYPE != "TRA",abs(all_sv$POS1-all_sv$POS2),0)
sample_info <- read.delim(sample_info_path,stringsAsFactors = FALSE)
sigprofile13 <- read.delim(sigprofile13_path)
bp_sherlock_bb <- read.delim(bp_sherlock_bb_path,stringsAsFactors = FALSE)
rmsk <- read.delim(rmsk_path,header = F,stringsAsFactors = FALSE)

#combine small  and  large
mytype <- c("DEL")
test <- all_sv[all_sv$SVTYPE %in% mytype,]
test$size <- abs(test$POS1 - test$POS2)
test$sizegroup <- ifelse(test$size>1000,"large","small")
length(unique(test[test$sizegroup=="small","SAMPLE"]))
nrow(test[test$sizegroup=="small",])

output_sample  <-  "TCGA-38-4630-01A"
output <- all_sv[all_sv$SAMPLE==output_sample,]
write.table(output,file.path("D:/OneDrive/OneDrive - The University of Chicago/AskForHelp/Sherlock",paste0(output_sample,"_allsv",".tsv")),row.names = F,quote = F,sep = "\t")
output_sample_bed <- output[output$SVTYPE!="TRA",]
output_sample_bed <- cbind(paste0("chr",output_sample_bed$CHROM1),
                           output_sample_bed$POS1-3000,
                           output_sample_bed$POS2+3000,
                           output_sample_bed$SVTYPE )
write.table(output_sample_bed,file.path("D:/OneDrive/OneDrive - The University of Chicago/AskForHelp/Sherlock",paste0(output_sample,"_region",".bed")),
            row.names = F,quote = F,sep = "\t",col.names = F)

############sample level
sample_small <- data.frame(test %>% group_by(SAMPLE,sizegroup) %>% count() %>% filter(sizegroup=="small"))
nrow(sample_small[sample_small$n<5,])

output_sample  <-  "NSLC-1219-T01"
output <- test[test$SAMPLE==output_sample,]
write.table(output,file.path("D:/OneDrive/OneDrive - The University of Chicago/AskForHelp/Sherlock",paste0(output_sample,"_small_",mytype,".tsv")),row.names = F,quote = F,sep = "\t")


############### make a cnv  bp table ##########
cnv_bp <- bind_rows(bp_sherlock_bb[,c("Tumor_Barcode","chr","startpos")],
                    bp_sherlock_bb[,c("Tumor_Barcode","chr","endpos")])
cnv_bp$pos <- ifelse(is.na(cnv_bp$startpos),cnv_bp$endpos,cnv_bp$startpos)
for (i in 1:nrow(all_sv)){
  samplei <- as.character(all_sv$SAMPLE[i])
  chr1 <- all_sv$CHROM1[i]
  chr2 <- all_sv$CHROM2[i]
  pos1 <- all_sv$POS1[i]
  pos2 <- all_sv$POS2[i]
  
  distance1 <-  min(abs(cnv_bp[cnv_bp$Tumor_Barcode == samplei &
                                 cnv_bp$chr == chr1, "pos"] - pos1))
  distance2 <-  min(abs(cnv_bp[cnv_bp$Tumor_Barcode == samplei &
                                 cnv_bp$chr == chr2, "pos"] - pos2))
  
  all_sv[i,"distance1"] <- distance1
  all_sv[i,"distance2"] <- distance2
  
}



ggplot(data = test,aes(x=ALGORITHMS,fill=sizegroup)) + 
  geom_bar() +
  ylab(paste(mytype,collapse = "-"))

test2 <- test[test$ALGORITHMS  %in% c("Manta","Meerkat-Manta"),]
ggplot(data = test2,aes(x=sizegroup,y=manta_PR)) + 
  geom_violin(size=0.5)+
  scale_y_log10() +
  xlab(paste(mytype,collapse = "-"))+
  theme(axis.text.x = element_text(size=6,colour = "black"),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_text(size=6,colour = "black"),
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
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(outpath,paste0("manta_pr_",paste(mytype,collapse = "-"),".pdf")),width = 1.3,height = 1.3)
ggplot(data = test2,aes(x=sizegroup,y=manta_SR)) + 
  geom_violin(size=0.5)+
  scale_y_log10() +
  xlab(paste(mytype,collapse = "-"))+
  theme(axis.text.x = element_text(size=6,colour = "black"),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_text(size=6,colour = "black"),
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
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(outpath,paste0("manta_sr_",paste(mytype,collapse = "-"),".pdf")),width = 1.3,height = 1.3)
ggplot(data = test2,aes(x=manta_PR,y=manta_SR,color=sizegroup,alpha=0.05)) + 
  geom_point()+
  scale_y_log10() +
  scale_x_log10()

test3 <- test[test$ALGORITHMS  %in% c("Meerkat","Meerkat-Manta"),]
ggplot(data = test3,aes(x=sizegroup,y=meerkat_DISC_PAIR)) + 
  geom_violin(size=0.5)+
  scale_y_log10() +
  xlab(paste(mytype,collapse = "-"))+
  theme(axis.text.x = element_text(size=6,colour = "black"),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_text(size=6,colour = "black"),
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
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(outpath,paste0("meerkat_pr_",paste(mytype,collapse = "-"),".pdf")),width = 1.3,height = 1.3)
ggplot(data = test3,aes(x=sizegroup,y=meerkat_SPLIT_READ)) + 
  geom_violin(size=0.5)+
  scale_y_log10() +
  xlab(paste(mytype,collapse = "-"))+
  theme(axis.text.x = element_text(size=6,colour = "black"),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_text(size=6,colour = "black"),
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
        plot.margin=grid::unit(c(1,2,0,2), "mm")
  )
ggsave(file.path(outpath,paste0("meerkat_sr_",paste(mytype,collapse = "-"),".pdf")),width = 1.3,height = 1.3)

ggplot(data = test,aes(x=manta_PR,y=manta_SR,color=sizegroup,alpha=0.05)) + 
  geom_point()+
  scale_y_log10() +
  scale_x_log10()

# with other information
all_samples_mysvtype <- test %>% group_by(SAMPLE,sizegroup) %>% count() %>% filter(sizegroup=="small")
all_samples_mysvtype <- merge(all_samples_mysvtype,sample_info,by.x = "SAMPLE",by.y = "tumor_barcode",all.x = T)
ggplot(data = all_samples_mysvtype,aes(x=histology,y=n)) + 
  geom_violin()+
  geom_jitter()+
  scale_y_log10() +
  ylab(paste0("Small ",mytype))


# sample  level
test_small <- test[test$sizegroup=="small",]

ggplot(data = test_small,aes(x=SAMPLE,fill=ALGORITHMS))+
  geom_bar() +
  ylab(paste0("Num of small ",mytype))+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank())


####### plot distance
plot_point <- function(df,xcol,ycol){
  ggplot(df,
         aes(
           x = df[,xcol],
           y = df[,ycol],
           color = sizegroup,
           alpha = 0.0001
         )) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    xlab(xcol) +
    ylab(ycol) +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))
}
plot_point(test,"distance1","distance2")

#######  rmsk
rmsk_annotation_function <- function(test){
  rmsk_test <- test
  rmsk_test$CHROM1 <- paste0("chr",rmsk_test$CHROM1)
  rmsk_test$CHROM2 <- paste0("chr",rmsk_test$CHROM2)
  for (i in 1:nrow(rmsk_test)){ #nrow(rmsk_test)
    # print(i)
    # print(rmsk_test[i,])
    chr1 <- rmsk_test$CHROM1[i]
    chr2 <- rmsk_test$CHROM2[i]
    pos1 <- rmsk_test$POS1[i]
    pos2 <- rmsk_test$POS2[i]
    tell1 <- rmsk$V1 == chr1 & rmsk$V2<= pos1 & rmsk$V3 >= pos1
    tell2 <- rmsk$V1 == chr2 & rmsk$V2<= pos2 & rmsk$V3 >= pos2
    if (sum(tell1)>0){
      rmsk_test[i,"rmsk1"] <- rmsk$V4[tell1][1]
    }
    if (sum(tell2)>0){
      rmsk_test[i,"rmsk2"] <- rmsk$V4[tell2][1]
    }
  }
  write.table(rmsk_test,file.path("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/99-qc",
                                  paste0(mytype,"rmsk_annotation.tsv")),
              sep = "\t",row.names = F,quote = F)
  return(rmsk_test)
}

mytype <- "DEL"
test <- all_sv[all_sv$SVTYPE==mytype,]
test$sizegroup <- ifelse(test$size>1000,"large","small")
del_rmsk  <-  rmsk_annotation_function(test)

mytype <- "DUP"
test <- all_sv[all_sv$SVTYPE==mytype,]
test$sizegroup <- ifelse(test$size>1000,"large","small")
dup_rmsk  <-  rmsk_annotation_function(test)

#########################################
library(ggplot2)

mytype <- "DUP"

mytype_annotation <-
  read.delim(
    file.path(
      "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/99-qc",
      paste0(mytype, "rmsk_annotation.tsv")
    )
  )

mytype_annotation$rmsk <- ifelse(is.na(mytype_annotation$rmsk1),0,1) + ifelse(is.na(mytype_annotation$rmsk2),0,1)
mytype_annotation$rmsk <- as.factor(mytype_annotation$rmsk)

mytype_annotation$alu <- ifelse(grepl("^Alu",mytype_annotation$rmsk1),1,0) +  ifelse(grepl("^Alu",mytype_annotation$rmsk2),1,0)
mytype_annotation$alu <- as.factor(mytype_annotation$alu)

rmsk_plot_function <- function(mytype_annotation,x){
  print(x)
  ggplot(mytype_annotation,aes(x=sizegroup,fill=mytype_annotation[,x])) +  
    geom_bar(position="fill")  +
    xlab(mytype) + 
    scale_fill_manual(breaks = c(0,1,2),labels = c("None", "One-end", "Both"),values = c("grey","#F1948A","#CB4335"))+
    guides(fill=guide_legend(title=x))
}

rmsk_plot_function(mytype_annotation,"rmsk")
rmsk_plot_function(mytype_annotation,"alu")
rmsk_plot_function(mytype_annotation[mytype_annotation$SAMPLE=="NSLC-1219-T01",],"rmsk")
rmsk_plot_function(mytype_annotation[mytype_annotation$SAMPLE=="NSLC-1219-T01",],"alu")


##### check samples
sizelist1 <- c(0,1000,10000,100000,1000000,10000000,100000000)
sizelist2 <- c(1000,10000,100000,1000000,10000000,100000000,1000000000)

mysample <- "NSLC-1219-T01"
mysample_sv <- data.frame(all_sv %>% filter(SAMPLE == mysample))
for (i in 1:nrow(mysample_sv)){
  sizei <-  mysample_sv$size[i]
  sizegroup <- paste0("sizegroup",grep(TRUE,sizei >=  sizelist1 & sizei < sizelist2))
  mysample_sv[i,"sizegroup"] <-  sizegroup
}

# plot manta PR SR
ggplot(mysample_sv[mysample_sv$ALGORITHMS %in% c("Manta","Meerkat-Manta"),],aes(x=manta_PR,y=manta_SR,color=SVTYPE,size=sizegroup,alpha=0.05)) + 
  geom_jitter(width=0.1,height = 0.1) +
  ggtitle(mysample)+
  scale_color_manual(values = c("#459FF4","#F2353A","#F5E115","#94F515","#F578F7"))+
  scale_x_log10()+
  scale_y_log10()
ggplot(mysample_sv[mysample_sv$ALGORITHMS %in% c("Meerkat","Meerkat-Manta"),],aes(x=meerkat_DISC_PAIR,y=meerkat_SPLIT_READ,color=SVTYPE,size=sizegroup,alpha=0.05)) + 
  geom_jitter(width=0.1,height = 0.1) +
  ggtitle(mysample)+
  scale_color_manual(values = c("#459FF4","#F2353A","#F5E115","#94F515","#F578F7"))+
  scale_x_log10()+
  scale_y_log10()



#### check hg39 and hg19  deletion samples
hg19_path <- file.path("D:/Data/Sherlock/PCAWG_LUAD/SV/meerkat/reshape")
small_del_samples <-  unique(sample_small$SAMPLE)
chrom_df <- data.frame()
for (samplei in small_del_samples){
  samplei_hg19_path <- file.path(hg19_path,paste0(samplei,"_meerkat.tsv"))
  if (file.exists(samplei_hg19_path)){
    print(samplei)
    print(sample_small[sample_small$SAMPLE==samplei,"n"])
    samplei_hg19  <- read.delim(samplei_hg19_path)
    
    small_del_hg38  <- all_sv[all_sv$SAMPLE==samplei & all_sv$SVTYPE=="DEL" & all_sv$size <1000,]
    small_del_hg19  <- samplei_hg19[samplei_hg19$SVTYPE=="DEL" & samplei_hg19$SVLEN <1000,]
    for (chri in c(seq(1,22,1),"X","Y")){
      chrom_df <- rbind(chrom_df,data.frame(
        "sample"=samplei,
        "chr"=chri,
        "hg19"=nrow(small_del_hg19[small_del_hg19$CHROM==chri,]),
        "hg38"=nrow(small_del_hg38[small_del_hg38$CHROM1==chri,])
      )) 
    }
  }
  
}

sum((chrom_df$hg38-chrom_df$hg19)[chrom_df$hg38-chrom_df$hg19 > 0])
sum(chrom_df$hg38)

ggplot(data=chrom_df,aes(x=hg38,y=hg19)) +  
  geom_jitter(width = 0.1,height = 0.1)
