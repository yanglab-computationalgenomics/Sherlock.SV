#### library
library(dplyr)
library(data.table)
library(ggplot2)
library(ggh4x)
library(ggpubr)


### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
sv.path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_all_sv_with_signature.tsv")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
outpath <- file.path(scratch.path,"24-smoking-sbs4","plot")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

### read
other_sample_infomation <- read.delim(sampleinfo.path, header=T, check.names=F, stringsAsFactors=F, na.strings="no_na")
sv <- read.delim(sv.path)
histology <- read.csv(histology.path)
other_sample_infomation <- merge(other_sample_infomation,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

### prepare
plot_df <- merge(other_sample_infomation,
                 sv %>% count(SAMPLE),
                 by.x = "Tumor_Barcode",by.y = "SAMPLE",all.x = T)
plot_df$n <- ifelse(is.na(plot_df$n),0,plot_df$n )
plot_df <- plot_df %>% filter(Smoking %in% c("Non-Smoker","Smoker"))

### plot
ggplot(plot_df, aes(x = Smoking, y = n, fill = Smoking)) +
  geom_violin(data=plot_df, aes(x = Smoking, y = n,fill = Smoking),trim = TRUE, color = "black") +
  geom_boxplot(data=plot_df, aes(x = Smoking, y = n),width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("Non-Smoker" = "#00a800", "Smoker" = "#FF6EC7"))+
  stat_compare_means(method = "wilcox.test", label.y = max(plot_df$n) * 1.05,size=2) + 
  scale_y_continuous(breaks = c(0,500,1000))+
  labs(y = "Number of total SVs") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 6, colour = "black", angle = 90),
    axis.text.y = element_text(size = 6, colour = "black"),
    axis.title.y = element_text(size = 6, colour = "black"),
    axis.title.x = element_blank(),
    axis.line.y = element_line(size = 0.1, colour = "black"),
    axis.line.x.top = element_line(size = 0.1, colour = "black"),
    axis.line.x.bottom = element_line(size = 0.1, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.234),
    axis.ticks.length.y = unit(0.8, "mm"),
    axis.ticks.y = element_line(linewidth = 0.234),
    plot.title = element_text(size = 6, colour = "black"),
    panel.background = element_blank(),
    panel.spacing.x = unit(1, "mm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text.x = element_text(size = 6, angle = 0, colour = "black", hjust = 1, vjust = 0.5),
    strip.text.y.left = element_text(size = 6, angle = 0, colour = "black", hjust = 1, vjust = 0),
    strip.background = element_blank(),
    strip.placement = "top",
    legend.key.height = unit(2, "mm"),
    legend.key.width = unit(2, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.margin = grid::unit(c(1, 2, 0, 2), "mm")
  )
ggsave(filename = file.path(outpath,"boxplot.ttest.totalsv.by.smoking.pdf"),width = 2,height = 2)





