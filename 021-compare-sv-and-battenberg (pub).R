# Yang Yang
# This is a script for compare bateernburg breakpoints with SV breakpoints of Sherlock


setwd("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock")

########################### library ##########################
library(ggplot2)
library(reshape)
library(scales)
library(dplyr)

########################### save path #########################
out_path <- file.path("D:","OneDrive","OneDrive - The University of Chicago","office","GitHub","Sherlock","scratch","012-Battenberg-and-sv")
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
    "01-SV-files",
    "1217_samples_hg38",
    "1217_manta_meerkat_union_window50.txt"
  )

# 
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



#################  read file ################
all_sv  <- read.delim(all_sv_path,stringsAsFactors = FALSE)
bp_sherlock_bb <- read.delim(bp_sherlock_bb_path,stringsAsFactors = FALSE)



############### make a cnv  bp table ##########
cnv_bp <- bind_rows(bp_sherlock_bb[,c("Tumor_Barcode","chr","startpos")],
                    bp_sherlock_bb[,c("Tumor_Barcode","chr","endpos")])
cnv_bp$pos <- ifelse(is.na(cnv_bp$startpos),cnv_bp$endpos,cnv_bp$startpos)



############## measure distances #############
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

# write.table(
#   all_sv,
#   file.path(
#     out_path,
#     "1217_manta_meerkat_union_window50(battenberg_annotation).txt"
#   ),
#   quote = FALSE,
#   row.names = FALSE,
#   sep = "\t"
# )


######### read distance #############
all_sv <- read.delim(file.path(
  out_path,
  "1217_manta_meerkat_union_window50(battenberg_annotation).txt"),stringsAsFactors = FALSE)

######## count ratio of meerkat support and cnv  support #############
range <- seq(0,50,1)
distance_cutoff1000 <- seq(0,1000,1)
distance_cutoff10000 <- seq(0,10000,1)
count_df <-  data.frame()
for (PRi  in  c(NA, range)) {
  for (SRi in c(NA, range)) {
    meerkat_support_count  <-
      as.numeric(
        all_sv  %>% filter(
          manta_PR %in% PRi &
            manta_SR  %in% SRi &
            ALGORITHMS %in% c("Meerkat-Manta")
        ) %>% count(.drop = FALSE)
      )
    distance_cutoff1000_support_count  <-
      as.numeric(
        all_sv  %>% filter(
          manta_PR %in% PRi &
            manta_SR  %in% SRi &
            ALGORITHMS %in% c("Meerkat-Manta", "Manta") &
            distance1 %in% distance_cutoff1000 &
            distance2 %in% distance_cutoff1000
        ) %>% count(.drop = FALSE)
      )
    distance_cutoff10000_support_count  <-
      as.numeric(
        all_sv  %>% filter(
          manta_PR %in% PRi &
            manta_SR  %in% SRi &
            ALGORITHMS %in% c("Meerkat-Manta", "Manta") &
            distance1 %in% distance_cutoff10000 &
            distance2 %in% distance_cutoff10000
        ) %>% count(.drop = FALSE)
      )
    all_count <- as.numeric(
      all_sv  %>% filter(
        manta_PR %in% PRi &
          manta_SR  %in% SRi &
          ALGORITHMS %in% c("Meerkat-Manta", "Manta")
      ) %>% count(.drop = FALSE)
    )
    count_df  <-  rbind(
      count_df,
      data.frame(
        "manta_PR" = PRi,
        "manta_SR" = SRi,
        "manta_count" = all_count,
        "meerkat_support_count" = meerkat_support_count,
        "distance_cutoff1000_support_count" = distance_cutoff1000_support_count,
        "distance_cutoff10000_support_count" = distance_cutoff10000_support_count
      )
    )
  }
}
count_df[is.na(count_df)]<- "NA"
count_df$manta_PR <- factor(count_df$manta_PR,levels = c("NA", range))
count_df$manta_SR <- factor(count_df$manta_SR,levels = c("NA", range))

#write.csv(count_df,file.path(out_path,"manta_SV_supported_by_meerkat_or_batternberg.csv"),quote = FALSE,row.names = FALSE)

######## count ratio of manta support and cnv  support #############
range <- seq(0,50,1)
distance_cutoff1000 <- seq(0,1000,1)
distance_cutoff10000 <- seq(0,10000,1)
count_df_meerkat <-  data.frame()
for (PRi  in  c(NA, range)) {
  for (SRi in c(NA, range)) {
    manta_support_count  <-
      as.numeric(
        all_sv  %>% filter(
          meerkat_DISC_PAIR %in% PRi &
            meerkat_SPLIT_READ  %in% SRi &
            ALGORITHMS %in% c("Meerkat-Manta")
        ) %>% count(.drop = FALSE)
      )
    distance_cutoff1000_support_count  <-
      as.numeric(
        all_sv  %>% filter(
          meerkat_DISC_PAIR %in% PRi &
            meerkat_SPLIT_READ  %in% SRi &
            ALGORITHMS %in% c("Meerkat-Manta", "Meerkat") &
            distance1 %in% distance_cutoff1000 &
            distance2 %in% distance_cutoff1000
        ) %>% count(.drop = FALSE)
      )
    distance_cutoff10000_support_count  <-
      as.numeric(
        all_sv  %>% filter(
          meerkat_DISC_PAIR %in% PRi &
            meerkat_SPLIT_READ  %in% SRi &
            ALGORITHMS %in% c("Meerkat-Manta", "Meerkat") &
            distance1 %in% distance_cutoff10000 &
            distance2 %in% distance_cutoff10000
        ) %>% count(.drop = FALSE)
      )
    all_count <- as.numeric(
      all_sv  %>% filter(
        meerkat_DISC_PAIR %in% PRi &
          meerkat_SPLIT_READ  %in% SRi &
          ALGORITHMS %in% c("Meerkat-Manta", "Meerkat")
      ) %>% count(.drop = FALSE)
    )
    count_df_meerkat  <-  rbind(
      count_df_meerkat,
      data.frame(
        "meerkat_PR" = PRi,
        "meerkat_SR" = SRi,
        "meerkat_count" = all_count,
        "manta_support_count" = manta_support_count,
        "distance_cutoff1000_support_count" = distance_cutoff1000_support_count,
        "distance_cutoff10000_support_count" = distance_cutoff10000_support_count
      )
    )
  }
}
count_df_meerkat[is.na(count_df_meerkat)]<- "NA"
count_df_meerkat$meerkat_PR <- factor(count_df_meerkat$meerkat_PR,levels = c("NA", range))
count_df_meerkat$meerkat_SR <- factor(count_df_meerkat$meerkat_SR,levels = c("NA", range))

#write.csv(count_df_meerkat,file.path(out_path,"meerkat_SV_supported_by_manta_or_batternberg.csv"),quote = FALSE,row.names = FALSE)
####### plot ###################
plot_point <- function(df,xcol,ycol){
  ggplot(df,
         aes(
           x = df[,xcol],
           y = df[,ycol],
           color = ALGORITHMS,
           alpha = 0.01
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
# all_sv$manta_PR <- ifelse(is.na(all_sv$manta_PR),"NA",all_sv$manta_PR)
# all_sv$manta_SR <- ifelse(is.na(all_sv$manta_PR),"NA",all_sv$manta_SR)
plot_point(all_sv,"distance1","distance2")
plot_point(all_sv,"manta_PR","manta_SR")

range <- seq(0,50,1)
count_df <- read.csv(file.path(out_path,"manta_SV_supported_by_meerkat_or_batternberg.csv"))
count_df[is.na(count_df)]<- "NA"
count_df$manta_PR <- factor(count_df$manta_PR,levels = c("NA", range))
count_df$manta_SR <- factor(count_df$manta_SR,levels = c("NA", range))
count_df_meerkat <- read.csv(file.path(out_path,"meerkat_SV_supported_by_manta_or_batternberg.csv"))
count_df_meerkat[is.na(count_df_meerkat)]<- "NA"
count_df_meerkat$meerkat_PR <- factor(count_df_meerkat$meerkat_PR,levels = c("NA", range))
count_df_meerkat$meerkat_SR <- factor(count_df_meerkat$meerkat_SR,levels = c("NA", range))

plot_count <-
  function(count_df,
           xcol,
           ycol,
           sizecol,
           colorcol1,
           colorcol2,
           minvalue,
           maxvalue,
           # midpointvalue,
           cololegendname,
           sizelegendname,
           xtitle,
           ytitle) {
    ggplot(
      count_df,
      aes(
        x = count_df[, xcol],
        y = count_df[, ycol],
        #size = 4,
        size = count_df[, sizecol],
        fill = count_df[, colorcol1] / count_df[, colorcol2]
      )
    ) + geom_point(shape = 21) +
      xlab(xcol) +
      ylab(ycol) +
      scale_size_continuous(range=c(0,4),breaks = c(1,10,100,1000),trans = "sqrt")+
      scale_fill_gradient2(
        low = "white",
        mid = "yellow",
        high = "black",
        midpoint = 0.1,
        limits = c(minvalue, maxvalue),
        # space = "Lab",
        na.value = "white",
        oob=squish,
        # guide = "colourbar",
        # aesthetics = "colour",
        name = cololegendname) + 
      labs(x=xtitle,y=ytitle) +
      theme(axis.text.x = element_text(size=6,colour = "black",angle = 90),
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
      )+
      guides(size = guide_legend(title = sizelegendname))
  }

#  for manta
plot_count(
  count_df,
  "manta_PR",
  "manta_SR",
  "manta_count",
  "meerkat_support_count",
  "manta_count",
  0,
  1,
  # 0.4,
  "Ratio of Meerkat-supported SV",
  "Number of Manta-called SV",
  "Paired-end Reads (Manta)",
  "Split Reads (Manta)"
)

plot_count(
  count_df %>% filter(distance_cutoff1000_support_count>0),
  "manta_PR",
  "manta_SR",
  "manta_count",
  "distance_cutoff1000_support_count",
  "manta_count",
  0,
  0.2,
  # 0.1,
  "Ratio of batternberg-CNV-supported SV (cutoff=1K)",
  "Number of Manta-called SV",
  "Paired-end Reads (Manta)",
  "Split Reads (Manta)"
)
ggsave(file.path(out_path,"mantasv.supporttingreads.with.1k.cnvsupport.pdf"),width = 3.8,height = 4)
plot_count(
  count_df %>% filter(distance_cutoff10000_support_count>0),
  "manta_PR",
  "manta_SR",
  "manta_count",
  "distance_cutoff10000_support_count",
  "manta_count",
  0,
  0.6,
  # 0.3,
  "Ratio of batternberg-CNV-supported SV (cutoff=10K)",
  "Number of Manta-called SV",
  "Paired-end Reads (Manta)",
  "Split Reads (Manta)"
)
ggsave(file.path(out_path,"mantasv.supporttingreads.with.10k.cnvsupport.pdf"),width = 3.8,height = 4)

# for meerkat
plot_count(
  count_df_meerkat,
  "meerkat_PR",
  "meerkat_SR",
  "meerkat_count",
  "manta_support_count",
  "meerkat_count",
  0,
  1,
  # 0.5,
  "Ratio of manta-supported",
  "Number of Manta-called SV",
  "Paired-end Reads (Meerkat)",
  "Split Reads (Meerkat)"
)

plot_count(
  count_df_meerkat,
  "meerkat_PR",
  "meerkat_SR",
  "meerkat_count",
  "distance_cutoff1000_support_count",
  "meerkat_count",
  0,
  0.3,
  # 0.15,
  "Ratio of batternberg-CNV-supported SV (cutoff=1K)",
  "Number of Manta-called SV",
  "Paired-end Reads (Meerkat)",
  "Split Reads (Meerkat)"
)

plot_count(
  count_df_meerkat,
  "meerkat_PR",
  "meerkat_SR",
  "meerkat_count",
  "distance_cutoff10000_support_count",
  "meerkat_count",
  0,
  0.7,
  # 0.3,
  "Ratio of batternberg-CNV-supported SV (cutoff=10K)",
  "Number of Manta-called SV",
  "Paired-end Reads (Meerkat)",
  "Split Reads (Meerkat)"
)




# count
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta")) %>% count() #280508
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0)) %>% count()  #14058
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(1,2,3)) %>% count()  #26367
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0) & manta_SR  %in% c(0)) %>% count() #66
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(2) & manta_SR  %in% c(2)) %>% count() #3329
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & is.na(manta_PR)) %>% count()  #0

all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & is.na(manta_SR)) %>% count()  #80670

all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0)) %>% group_by(SVTYPE) %>% count()
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(2)  & manta_SR %in% c(2)) %>% group_by(SVTYPE) %>% count()
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(1,2,3)) %>% group_by(SVTYPE) %>% count()
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & is.na(manta_SR)) %>% group_by(SVTYPE) %>% count()
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0) & manta_SR  %in% c(0)) %>% group_by(SVTYPE) %>% count()


a=all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0) & SVTYPE %in% c("DEL"))
a=all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0) & SVTYPE %in% c("DUP"))
a=all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0) & SVTYPE %in% c("TRA"))
all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(2) & manta_SR  %in% c(2) & SVTYPE %in% c("TRA")) %>% count() #3329
a=all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & manta_PR %in% c(0,1,2,3))

all_sv  %>% filter(ALGORITHMS %in% c("Manta","Meerkat-Manta") & is.na(manta_HOMLEN) & is.na(manta_SVINSLEN)) %>% count() #120270


