# library
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(ggpubr)
library(data.table)
library(rstatix)

#### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock/"
scratch.path <- file.path(sherlock.path,"scratch")
sv_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_all_sv_with_signature.tsv")
sample_info_path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
hq.path <- file.path(sherlock.path,"/Data/SAMPLE_INFO/HQ_samples.csv")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")

sample_info <- read.delim(sample_info_path,stringsAsFactors = FALSE)
output_dir <- file.path(scratch.path,"02-sv-discription","1217_samples_hg38")
if (!dir.exists(output_dir)){
  dir.create(output_dir,recursive = T,showWarnings = FALSE)
}

### read
histology <- read.csv(histology.path)
sample_info <- merge(sample_info,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

all_sv  <- read.delim(sv_path,stringsAsFactors = FALSE)
allsamples <-  unique(sample_info$Tumor_Barcode)
hqsample <- read.csv(hq.path) %>% pull(Tumor_Barcode)

### prepare plot df
plot.df <- all_sv %>% count(SAMPLE)
# add zero sample
plot.df <- merge(plot.df,sample_info[,c("Tumor_Barcode","Histology")],all.y = T,by.x = "SAMPLE",by.y = "Tumor_Barcode")
plot.df <- plot.df %>%
  mutate(n=ifelse(is.na(n),0,n)) %>%
  mutate(Tumortype=Histology)
# filter hg quality
#plot.df <- plot.df %>% filter(SAMPLE %in% hqsample)
# add histology number
plot.df <- merge(plot.df,
                 plot.df %>% group_by(Tumortype) %>% dplyr::count() %>% data.table::setnames(.,new="Tumortype.num",old="n"),
                 all.x = T)
plot.df$Tumortype.num.full <- paste0(plot.df$Tumortype," (",plot.df$Tumortype.num,")")  
plot.df <- plot.df %>% dplyr::group_by(Tumortype.num.full) %>% dplyr::mutate(median.tumor=median(n)) %>% as.data.frame()
# order
plot.df <- plot.df %>% arrange(desc(median.tumor),desc(n))
tumor.order <- unique(plot.df$Tumortype.num.full)
tumor.order <- c(tumor.order[grepl("Others|Undetermined ",tumor.order)==FALSE],
                 tumor.order[grepl("Others|Undetermined ",tumor.order)==TRUE])
sample.order <- unique(plot.df$SAMPLE)
plot.df <- plot.df %>% mutate(SAMPLE=factor(SAMPLE,levels = sample.order),
                              Tumortype.num.full=factor(Tumortype.num.full,levels = tumor.order))
# add simple/complex/uncluster sv number
plot.df <- merge(plot.df,
                 all_sv %>% filter(Signature %in% c("1 ecDNA/double minutes","2 BFB cycles/chromatin bridge","3 Large loss",
                                                    "4 Micronuclei","5 Large gain","6 Hourglass","un_assigned")) %>% dplyr::group_by(SAMPLE) %>% dplyr::count() %>% setnames(.,old="n",new="cluster cplx") %>% as.data.frame(),
                 all.x=T) # complex sv number
plot.df <- merge(plot.df,
                     all_sv %>% filter(Signature %in% c("chromoplexy","complex_unclear","cycle_templated_ins"))  %>% dplyr::group_by(SAMPLE) %>% dplyr::count() %>% setnames(.,old="n",new="uncluster cplx"),
                     all.x=T) # uncluster sv number
plot.df <- merge(plot.df,
                 all_sv %>% filter(Signature %in% c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra"))  %>% dplyr::group_by(SAMPLE) %>% dplyr::count() %>% setnames(.,old="n",new="simple"),
                     all.x=T) # simple sv number
plot.df <- plot.df %>% mutate_at(c("cluster cplx","uncluster cplx","simple"),~replace_na(.,0)) %>% as.data.frame()

### 1.plot freq
plotfile.median <- unique(plot.df[,c("Tumortype.num.full","Tumortype.num","median.tumor")]) %>%setnames(old="median.tumor",new="median")
plotfile.iqr <- plot.df %>%
  dplyr::group_by(Tumortype.num.full, Tumortype.num) %>%
  dplyr::summarise(
    Q1 = quantile(n, 0.25, na.rm = TRUE),
    Q3 = quantile(n, 0.75, na.rm = TRUE),
    median = median(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(IQR = Q3 - Q1)
p.freq<-ggplot(data=plot.df,aes(x=SAMPLE,y=n))  + 
  geom_point(size=0.5/2.83)+
  facet_grid2(.~Tumortype.num.full,
              scales = "free_x",
              space = "free_x",
              switch = "both",
              strip = strip_vanilla(clip = "off")) +
  geom_hline(data = plotfile.median, aes(yintercept = median),color="red",size=0.5/2.83)+
  geom_hline(data = plotfile.iqr, aes(yintercept = Q1), color = "blue", linetype = "dashed", linewidth = 0.3) +
  geom_hline(data = plotfile.iqr, aes(yintercept = Q3), color = "blue", linetype = "dashed", linewidth = 0.3)+
  geom_text(data = plotfile.median,aes(x=Tumortype.num/2,y=median,label = median),color="red",size=2,vjust=-1)+
  scale_y_log10() +
  scale_x_discrete(expand = c(0, 0))+
  ylab("Number of somatic SVs") + 
  xlab("Sample") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.234,colour = "black"),
        axis.line.x.top = element_line(linewidth = 0.234,colour = "black"), 
        axis.line.x.bottom = element_line(linewidth = 0.234,colour = "black"), 
        axis.ticks.x =  element_blank(),
        axis.ticks.y =  element_line(linewidth = 0.234),
        axis.ticks.length.y =unit(0.8, "mm"),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        # strip.position="top",
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(2,"mm"),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = "right",
        plot.margin=grid::unit(c(2,2,0,2), "mm")
  ) + theme(strip.text.x = element_text(size=6, angle = 90, colour = "black",hjust = 0.5,vjust=0.5))

# 2. plot.percentage
p.percentage <- ggplot(data=plot.df %>% tidyr::gather(.,Signaturetype,svnum,"cluster cplx":"simple") %>% 
         mutate(Signaturetype=factor(Signaturetype,levels = c("cluster cplx","uncluster cplx","simple")),
                SAMPLE=factor(SAMPLE,levels = sample.order)),
       aes(x=SAMPLE,y=svnum,fill=Signaturetype))  + geom_bar(stat = "identity",position = "fill") +
  facet_grid2(.~Tumortype.num.full,
              scales = "free_x",
              space = "free_x",
              switch = "both",
              strip = strip_vanilla(clip = "off")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0),position = "top")+
  xlab("Sample") + 
  ylab("SV percentage") +
  scale_fill_manual(breaks = c("cluster cplx","uncluster cplx","simple"),
                    values = c("#F4924A","#91E83D","#6E9EEA"),
                    labels = c("Clustered complex SVs","Non-clustered complex SVs","Simple SVs"),
                    name=fill) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size=6,colour = "black"),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.234),
        axis.ticks.x =  element_blank(),
        axis.ticks.y =  element_line(linewidth = 0.234),
        axis.ticks.length.y =unit(0.8, "mm"),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.key.height = unit(2,"mm"),
        legend.key.width = unit(2,"mm"),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        plot.margin=grid::unit(c(0,2,2,2), "mm")
  )
# merge plots
p=ggarrange(p.freq,p.percentage,align = c("hv"),ncol = 1,nrow = 2,heights = c(0.5,0.5))
ggsave(file.path(output_dir,"complexandsimplesv.freq.across.tumortypes.pdf"),width = 85,height = 130,units = c("mm"))
#ggsave(file.path(output_dir,"complexandsimplesv.freq.across.tumortypes_hq.pdf"),width = 110,height = 200,units = c("mm"))


### a Wilcoxon nonparametric test to look at whether these medians (n) between different tumor types are different
# Step 1: Filter out "Others"
plot.df.filtered <- plot.df %>%
  filter(Tumortype != "Others") %>%
  mutate(Tumortype = factor(Tumortype))

# Step 2: Get all pairwise comparisons
tumor_pairs <- combn(levels(plot.df.filtered$Tumortype), 2, simplify = FALSE)

# Step 3: Perform Wilcoxon tests with FDR adjustment
wilcox_results <- plot.df.filtered %>%
  wilcox_test(n ~ Tumortype, comparisons = tumor_pairs) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

# Step 4: Add median n values for each group
medians <- plot.df.filtered %>%
  group_by(Tumortype) %>%
  summarise(median_n = median(n))

# Step 5: Join medians for each group1 and group2
wilcox_results <- wilcox_results %>%
  left_join(medians, by = c("group1" = "Tumortype")) %>%
  dplyr::rename(median_group1 = median_n) %>%
  left_join(medians, by = c("group2" = "Tumortype")) %>%
  dplyr::rename(median_group2 = median_n)

# Step 6: View the final table
print(wilcox_results)
