### alluvival plot
install.packages("devtools")
devtools::install_github("davidsjoberg/ggsankey")
### library
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(ggsankey)

sankey.function <- function(plot.df,myorder=c("Histology","Population","Sex","Smoking")){
  plot.df.sankey <- plot.df %>% 
    select(Histology,Population,Sex,Smoking) %>% 
    mutate(
      population=factor(Population,levels = c("a.European","b.Asian","c.Others")),
      smoking=factor(Smoking,levels = c("a.Never smoker","b.Smoker","c.Unknown")),
      gender=factor(Sex,levels = c("Male","Female")),
      histology=factor(Histology,levels = c("a.Adenocarcinoma",
                                            "b.Squamous cell carcinoma",
                                            "c.Carcinoid tumor",
                                            "d.Adenosquamous carcinoma",
                                            "e.Others"))
    ) %>% arrange(.,Histology,Population,Sex,Smoking) %>%
    ggsankey::make_long(myorder)
  
  ggplot(plot.df.sankey, aes(x = x, next_x = next_x, 
                             node = node, next_node = next_node, 
                             fill = factor(node), label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "grey30",width = 0.2) +
    geom_sankey_label(size = 2, 
                      fill = NA) +
    scale_fill_manual("",values=c("a.European"="#FFD92F", 
                                  "b.Asian"="#66C2A5", 
                                  "c.Others"="#7C7C7C", 
                                  "a.Never smoker"="#00a800", 
                                  "b.Smoker"="#FF6EC7", 
                                  "c.Unknown"="#F0F0F0",
                                  "a.Adenocarcinoma" = "#4CAF50",  
                                  "b.Squamous cell carcinoma" = "#2196F3",  
                                  "c.Carcinoid tumor" = "#FF9800",  
                                  "d.Adenosquamous carcinoma" = "#9C27B0",  
                                  "e.Others" = "#F0F0F0", 
                                  "Male"="#8DA0CB", 
                                  "Female" = "#FC8D62")) + 
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5))
}

### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
hq.path <- file.path(sherlock.path,"/Data/SAMPLE_INFO/HQ_samples.csv")
sampleinfo.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
outpath <- file.path(sherlock.path,"scratch","17-alluvival")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}
### read
sampleinfo <-  read.delim(sampleinfo.path)
hqsample <- read.csv(hq.path) %>% pull(Tumor_Barcode)
histology <- read.csv(histology.path)
sampleinfo <- merge(sampleinfo,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

### prepare
# all sample
plot.df <- sampleinfo %>% 
  select(Assigned_Population,Smoking,Gender,Histology) %>% 
  mutate(Assigned_Population=ifelse(Assigned_Population %in% c("AFR","AMR or Mixed"),"Others",Assigned_Population))
# # high quality sample
# plot.df <- sampleinfo %>% 
#   filter(Tumor_Barcode %in% hqsample) %>%
#   select(Assigned_Population,Smoking,Gender,Histology)  %>% 
#   mutate(Assigned_Population=ifelse(Assigned_Population %in% c("AFR","AMR or Mixed"),"Others",Assigned_Population))



# 1. population
plot.df$Population <- plot.df$Assigned_Population
plot.df$Population <- ifelse(plot.df$Population=="EUR","a.European",plot.df$Population)
plot.df$Population <- ifelse(plot.df$Population=="EAS","b.Asian",plot.df$Population)
plot.df$Population <- ifelse(plot.df$Population=="Others","c.Others",plot.df$Population)
# 2. smoking
plot.df$Smoking <- ifelse(plot.df$Smoking=="Non-Smoker","a.Never smoker",plot.df$Smoking)
plot.df$Smoking <- ifelse(plot.df$Smoking=="Smoker","b.Smoker",plot.df$Smoking)
plot.df$Smoking <- ifelse(plot.df$Smoking=="Unknown","c.Unknown",plot.df$Smoking)
# 3. Histology
plot.df$Histology <- ifelse(plot.df$Histology %in% c("Adenocarcinoma"),"a.Adenocarcinoma",plot.df$Histology)
plot.df$Histology <- ifelse(plot.df$Histology %in% c("Squamous cell carcinoma"),"b.Squamous cell carcinoma",plot.df$Histology)
plot.df$Histology <- ifelse(plot.df$Histology %in% c("Carcinoid tumor"),"c.Carcinoid tumor",plot.df$Histology)
plot.df$Histology <- ifelse(plot.df$Histology %in% c("Adenosquamous carcinoma"),"d.Adenosquamous carcinoma",plot.df$Histology)
plot.df$Histology <- ifelse(plot.df$Histology %in% c("Others"),"e.Others",plot.df$Histology)
# 4. gender
plot.df$Sex <- plot.df$Gender

# plot
# 1. all sample
sankey.function(plot.df,myorder=c("Histology","Population","Sex","Smoking"))
ggsave(file.path(outpath,"alluvia.plot.sample.info_his.pop.sex.smk.pdf"),width = 120,height = 70,units="mm",dpi = 300)
sankey.function(plot.df,myorder=c("Population","Sex","Smoking","Histology"))
ggsave(file.path(outpath,"alluvia.plot.sample.info_pop.sex.smk.his.pdf"),width = 80,height = 70,units="mm",dpi = 300)
# # 2. high quality sample
# sankey.function(plot.df,myorder=c("Histology","Population","Sex","Smoking"))
# ggsave(file.path(outpath,"alluvia.plot.sample.info_his.pop.sex.smk_hq.pdf"),width = 120,height = 70,units="mm",dpi = 300)
# sankey.function(plot.df,myorder=c("Population","Sex","Smoking","Histology"))
# ggsave(file.path(outpath,"alluvia.plot.sample.info_pop.sex.smk.his_hq.pdf"),width = 80,height = 70,units="mm",dpi = 300)


plot.df %>% count(Histology)
plot.df %>% count(Assigned_Population)
plot.df %>% count(Gender)
plot.df %>% count(Smoking)













# 
# 
# 
# is_alluvia_form(as.data.frame(plot.df), axes = 1:4, silent = TRUE)
# 
# ### plot
# ggplot(as.data.frame(plot.df),
#        aes(y = n,
#            axis1 = histology, axis3 = population)) +
#   geom_alluvium(aes(fill = histology),width = 1/4, knot.pos = 0.2, reverse = FALSE) +
#   scale_fill_manual(values = c("#3cb44b", "#c0d849", "#ff9540", "#bb56a0", "#ededed", "#b3b3b3"),
#                     breaks = c("Adenocarcinoma","Squamous carcinoma","Adenosquamous carcinoma","Carcinoid tumor","NA","Others")
#                     ) +
#   guides(fill = "none") +
#   geom_stratum(alpha = 0.5, width = 1/4, reverse = FALSE,fill = "white", color = "grey30")+
#   geom_text(stat = "stratum",aes(label = after_stat(stratum)),reverse = FALSE,size=2,colour="black") +
#   scale_x_continuous(breaks = 1:4, labels = c("Tumor type","Sex","Smoking","Population"),expand = c(0,0)) +
#   # coord_flip() +
#   theme_minimal()+
#   theme(
#     axis.ticks = element_blank(),
#     axis.title = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(size=6,color="black"),
#     panel.background = element_blank(), # remove the gray background
#     panel.grid.major = element_blank(), # remove the major grid lines
#     panel.grid.minor = element_blank()# remove the minor grid lines
#   ) +
#   theme(axis.line=element_blank())+
#   theme(panel.background=element_blank())+
#   theme(panel.border=element_blank())
# ggsave(file.path(outpath,"alluvia.plot.sample.info.pdf"),width = 200,height = 80,units="mm",dpi = 300)
