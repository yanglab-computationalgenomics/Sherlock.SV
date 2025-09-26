setwd("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock")

library(ggplot2)
library(reshape2)
library(customLayout)
library(dplyr)
library(tidyr)

root_dir <- getwd()

out_path <- file.path(root_dir,"scratch","PCAWG_LUAD","compare")


######################################################FUNCTION##############################
alg_test <- function(gold, screen) {
  a <- sum(gold & screen)
  b <- sum(gold==FALSE & screen)
  c <- sum(gold & screen==FALSE)
  # a2 <- sum(gold & screen & brass_snowman_dranger_union_find)
  # b2 <- sum(gold==FALSE & screen & brass_snowman_dranger_union_find)
  # c2 <- sum(gold & screen==FALSE & brass_snowman_dranger_union_find)
  t_p <- sum(gold & screen) / sum(gold)
  t_n <- sum(gold==FALSE & screen==FALSE) / sum(gold==FALSE)
  #print(sum(gold==FALSE & screen==FALSE))
  accuracy <-
    (sum(gold * screen) + sum(gold == FALSE &
                                screen == FALSE)) / (
                                  sum(gold * screen) + sum(gold == FALSE &
                                                             screen == FALSE) + sum(gold == FALSE &
                                                                                      screen) + sum(gold & screen == FALSE)
                                )
  out <- c(t_p, t_n,accuracy,
           a,b,c
           # a2,b2,c2
           )
  #print(c("sensitivity", "specificity","accuracy"))
  return(out)
}
list2df <- function(x){
  out = data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE))
  names(out) <- c("Sensitivity","Specificity","Accuracy","TP","FP","FN"
                  # "TTP","FFP","TFN"
                  )
  return(out)
}
##############################################################################################

compare_list <- list.files(file.path(root_dir,"scratch","PCAWG_LUAD","compare","tables"))
compare_list <- compare_list[grepl("_compare500.csv",compare_list)]
compare_list <- compare_list[grepl("_compare10.csv",compare_list)]
compare_list <- compare_list[grepl("_compare50.csv",compare_list)]

all_compare <- data.frame()

manta_list <- list()
meerkat_list <- list()
tnscope_list <- list()

manta_meerkat_intersection_list <- list()
manta_tnscope_intersection_list <- list()
meerkat_tnscope_intersection_list <- list()
manta_meerkat_tnscope_intersection_list <- list()

manta_meerkat_union_list <- list()
manta_tnscope_union_list <- list()
meerkat_tnscope_union_list <- list()
manta_meerkat_tnscope_union_list <- list()

goldmanta_screenmeerkat_list <- list()
goldmanta_screentnscopre_list <- list()
goldmanta_screenmanta_list <- list()
goldmeerkat_screenmanta_list <- list()
goldmeerkat_screentnscope_list <- list()
goldmeerkat_screenmeerkat_list <- list()
goldtnscope_screenmanta_list <- list()
goldtnscope_screenmeerkat_list <- list()
goldtnscope_screentnscope_list <- list()

any2_list <- list()
for (i in compare_list){
  sample <- substr(i,1,16)
  print(sample)
  compare <- read.csv(file.path(root_dir,"scratch","PCAWG_LUAD","compare","tables",i))
  compare$sample <- i
  all_compare <- rbind(all_compare,compare)
  
  pcawg_find <- compare$pcawg != "0"
  manta_find <- compare$manta != "0"
  meerkat_find <- compare$meerkat != "0"
  tnscope_find <- compare$tnscope != "0"
  # brass_snowman_dranger_union_find <- compare$brass != "0" | compare$snowman != "0" | compare$dranger != "0"
  manta_list[[sample]] <- alg_test(pcawg_find,manta_find)
  meerkat_list[[sample]] <- alg_test(pcawg_find,meerkat_find)
  tnscope_list[[sample]] <- alg_test(pcawg_find,tnscope_find)
  manta_meerkat_intersection_list[[sample]] <- alg_test(pcawg_find,manta_find*meerkat_find)
  manta_tnscope_intersection_list[[sample]] <- alg_test(pcawg_find,manta_find*tnscope_find)
  meerkat_tnscope_intersection_list[[sample]] <- alg_test(pcawg_find,meerkat_find*tnscope_find)
  manta_meerkat_tnscope_intersection_list[[sample]] <- alg_test(pcawg_find,manta_find*meerkat_find*tnscope_find)
  manta_meerkat_union_list[[sample]] <- alg_test(pcawg_find,manta_find|meerkat_find)
  manta_tnscope_union_list[[sample]] <- alg_test(pcawg_find,manta_find|tnscope_find)
  meerkat_tnscope_union_list[[sample]] <- alg_test(pcawg_find,meerkat_find|tnscope_find)
  manta_meerkat_tnscope_union_list[[sample]] <- alg_test(pcawg_find,manta_find|meerkat_find|tnscope_find)
  any2_list[[sample]] <- alg_test(pcawg_find,manta_find*meerkat_find|manta_find*tnscope_find|meerkat_find*tnscope_find)
  goldmanta_screenmeerkat_list[[sample]] <- alg_test(manta_find,meerkat_find)
  goldmanta_screentnscopre_list[[sample]] <- alg_test(manta_find,tnscope_find)
  goldmanta_screenmanta_list[[sample]] <- alg_test(manta_find,manta_find)
  goldmeerkat_screenmanta_list[[sample]] <- alg_test(meerkat_find,manta_find)
  goldmeerkat_screentnscope_list[[sample]] <- alg_test(meerkat_find,tnscope_find)
  goldmeerkat_screenmeerkat_list[[sample]] <- alg_test(meerkat_find,meerkat_find)
  goldtnscope_screenmanta_list[[sample]] <- alg_test(tnscope_find,manta_find)
  goldtnscope_screenmeerkat_list[[sample]]<- alg_test(tnscope_find,meerkat_find)
  goldtnscope_screentnscope_list[[sample]]<- alg_test(tnscope_find,tnscope_find)
}


manta_df <- list2df(manta_list)
meerkat_df <- list2df(meerkat_list)
tnscope_df <- list2df(tnscope_list)
manta_meerkat_intersection_df <- list2df(manta_meerkat_intersection_list)
manta_tnscope_intersection_df <- list2df(manta_tnscope_intersection_list)
meerkat_tnscope_intersection_df <- list2df(meerkat_tnscope_intersection_list)
manta_meerkat_tnscope_intersection_df <- list2df(manta_meerkat_tnscope_intersection_list)
manta_meerkat_union_df <- list2df(manta_meerkat_union_list)
manta_tnscope_union_df <- list2df(manta_tnscope_union_list)
meerkat_tnscope_union_df <- list2df(meerkat_tnscope_union_list)
manta_meerkat_tnscope_union_df <- list2df(manta_meerkat_tnscope_union_list)
any2_df <- list2df(any2_list)

goldmanta_screenmeerkat_df <- list2df(goldmanta_screenmeerkat_list)
goldmanta_screentnscope_df <- list2df(goldmanta_screentnscopre_list)
goldmanta_screenmanta_df <- list2df(goldmanta_screenmanta_list)
goldmeerkat_screenmanta_df <- list2df(goldmeerkat_screenmanta_list)
goldmeerkat_screentnscope_df <- list2df(goldmeerkat_screentnscope_list)
goldmeerkat_screenmeerkat_df <- list2df(goldmeerkat_screenmeerkat_list)
goldtnscope_screenmanta_df <- list2df(goldtnscope_screenmanta_list)
goldtnscope_screenmeerkat_df<- list2df(goldtnscope_screenmeerkat_list)
goldtnscope_screentnscope_df<- list2df(goldtnscope_screentnscope_list)

df=list(manta_df,meerkat_df,tnscope_df,any2_df,
        manta_meerkat_intersection_df,manta_tnscope_intersection_df,meerkat_tnscope_intersection_df,manta_meerkat_tnscope_intersection_df,
        manta_meerkat_union_df,manta_tnscope_union_df,meerkat_tnscope_union_df,manta_meerkat_tnscope_union_df)
plot_name=c("01.Manta","02.Meerkat","03.TNscope","04.Any 2 of_Manta/Meerkat/TNscope",
            "05.Manta+Meerkat_intersection","06.Manta+TNscope_intersection","07.Meerkat+TNscope_intersection","08.Manta+Meerkat+TNscope_intersection",
            "09.Manta+Meerkat_union","10.Manta+TNscope_union","11.Meerkat+TNscope_union","12.Manta+Meerkat+TNscope_union")
plot_df <- as.data.frame(t(data.frame(matrix(unlist(lapply(df, function(x) apply(x[,4:6],2,sum))),nrow = 3)))) 
rownames(plot_df) <- plot_name
colnames(plot_df) <- c("TP","FP","FN")
plot_df$combination <- rownames(plot_df)
ggplot(data=plot_df %>% tidyr::gather(.,compare,number,TP:FN) %>% mutate(compare=factor(compare,levels = c("FP","TP","FN"))))+
  geom_bar(aes(x=compare,y=number,fill=compare),stat = "identity")+
  geom_text(aes(x=compare,y=number,label=number),size=2,hjust=0.5,vjust=-0.5)+
  facet_wrap(.~combination,nrow=3)+
  scale_fill_manual(values=c('#3498DB', "#CB4335","#F2AE10"))+
  ylab("SV Number")+
  theme_bw() + 
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.234),
        axis.ticks.x =  element_blank(),
        axis.ticks.y =  element_line(linewidth = 0.234),
        axis.ticks.length.y =unit(0.8, "mm"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6,angle = 0,colour = "black"),
        axis.text.y = element_text(size=6,colour = "black"),
        axis.title.y = element_text(size = 6),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "mm"),
        panel.grid = element_blank(),
        #panel.border = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(color="black", size=10, face="bold.italic",hjust = 0.5),
        strip.text.x = element_text(size = 6, colour = "black"),
        legend.position = "none"
  )

ggsave(file.path(out_path,paste0("alg_compare_500bp.pdf")),width = 6,height = 4)
ggsave(file.path(out_path,paste0("alg_compare_50bp.pdf")),width = 6,height = 4)
# pcawg_find <- all_compare$pcawg != "0"
# manta_find <- all_compare$manta != "0"
# meerkat_find <- all_compare$meerkat != "0"
# tnscope_find <- all_compare$tnscope != "0"
# 
# alg_test(pcawg_find,manta_find)
# alg_test(pcawg_find,meerkat_find)
# alg_test(pcawg_find,tnscope_find)
# alg_test(pcawg_find,manta_find*meerkat_find)
# alg_test(pcawg_find,manta_find*tnscope_find)
# alg_test(pcawg_find,meerkat_find*tnscope_find)
# alg_test(pcawg_find,manta_find*meerkat_find*tnscope_find)
# alg_test(pcawg_find,manta_find|meerkat_find)
# alg_test(pcawg_find,manta_find|tnscope_find)
# alg_test(pcawg_find,meerkat_find|tnscope_find)
# alg_test(pcawg_find,manta_find|meerkat_find|tnscope_find)
# alg_test(pcawg_find,manta_find*meerkat_find|manta_find*tnscope_find|meerkat_find*tnscope_find)


#############################################plot1##################################
# plot_function <- function(df,plot_name){
#   plot_df <- df[,c("Sensitivity","Specificity","Accuracy")]
#   plot_df <- melt(plot_df)
#   p <- ggplot(data=plot_df,aes(x=variable,y=value,fill=variable))+
#     geom_violin()+
#     scale_fill_manual(values=c('#CB4335', "#3498DB","#0BC417"))+
#     geom_jitter(width=0.2,stat = "identity")+
#     ylab("Rate")+
#     ggtitle(plot_name)+
#     theme_bw() + 
#     theme(panel.background = element_blank(),
#           axis.title.x = element_blank(),
#           axis.text.x = element_text(size=10,angle = 00),
#           axis.text.y = element_text(size=13),
#           axis.title.y = element_text(size = 15),
#           plot.title = element_text(color="black", size=10, face="bold.italic",hjust = 0.5),
#           legend.position = "none"
#     )
#   return(p)
# }
# 
# p1 <- plot_function(manta_df,"1.Manta")
# p2 <- plot_function(meerkat_df,"2.Meerkat")
# p3 <- plot_function(tnscope_df,"3.TNscope")
# p4 <- plot_function(manta_meerkat_intersection_df,"4.Manta+Meerkat_intersection")
# p5 <- plot_function(manta_tnscope_intersection_df,"5.Manta+TNscope_intersection")
# p6 <- plot_function(meerkat_tnscope_intersection_df,"6.Meerkat+TNscope_intersection")
# p7 <- plot_function(manta_meerkat_tnscope_intersection_df,"7.Manta+Meerkat+TNscope_intersection")
# p8 <- plot_function(manta_meerkat_union_df,"8.Manta+Meerkat_union")
# p9 <- plot_function(manta_tnscope_union_df,"9.Manta+TNscope_union")
# p10 <- plot_function(meerkat_tnscope_union_df,"10.Meerkat+TNscope_union")
# p11 <- plot_function(manta_meerkat_tnscope_union_df,"11.Manta+Meerkat+TNscope_union")
# p12 <- plot_function(any2_df,"12.Any 2 of Manta/Meerkat/TNscope")
# 
# # merge 7 sig
# lay1 <- lay_new(matrix(1:12, ncol = 4))
# cl <- lay1
# combine <- lay_grid(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12), cl)  
# ggsave(file.path(out_path,paste0("4alg_compare.png")),plot=combine,width = 11.5,height = 8,units = "in")
# ggsave(file.path(out_path,paste0("4alg_compare.pdf")),plot=combine,width = 11.5,height =8,units = "in")


#############################################plot2##################################
plot_function2 <- function(df,plot_name,facetcol){
  plot_df <- as.data.frame(t(data.frame(matrix(unlist(lapply(df, function(x) apply(x[,4:9],2,sum))),nrow = 6))))
  rownames(plot_df) <- plot_name
  colnames(plot_df) <- c("TP","FP","FN","TTP","FFP","TFN")
  plot_df$name <- plot_name
  plot_df_new <- data.frame()
  new_matrix <- data.frame(plot_df$TP-plot_df$TTP,
                      plot_df$TP-plot_df$TP+plot_df$TTP,
                      plot_df$FP-plot_df$FFP,
                      plot_df$FP-plot_df$FP+plot_df$FFP,
                      plot_df$FN-plot_df$TFN,
                      plot_df$FN-plot_df$FN+plot_df$TFN
  )
  names(new_matrix) <- c("TP-notDetect","TP-Detect","FP-notDetect","FP-Detect","FN-notDetect","FN-Detect")
  new_matrix <- c(t(new_matrix))
  new_matrix <- as.data.frame(new_matrix)
  new_matrix$V1 <- rep(c("TP","TP","FP","FP","FN","FN"),nrow(new_matrix)/6)
  new_matrix$V2 <- rep(c("notDetect","Detect"),nrow(new_matrix)/2)
  new_matrix$method <- gsub("_","\n",rep(plot_name,each=6))
  new_matrix$V1 <- factor(new_matrix$V1,levels = c("FP","TP","FN"))
  new_matrix$V2 <- factor(new_matrix$V2,levels = c("notDetect","Detect"))
  
  
  # plot_df <- melt(plot_df)
  # plot_df$variable <- factor(plot_df$variable,levels = c("FP","TP","FN","FFP","TTP","TFN"))
  # print(plot_df)
  p <- ggplot(data=new_matrix,aes(x=V1,y=new_matrix,fill=V1
                                  # alpha=V2
                                  ))+
    geom_bar(stat = "identity")+
    facet_wrap(~method,ncol = facetcol
               )+
    scale_fill_manual(values=c('#3498DB', "#CB4335","#F2AE10"))+
    #scale_alpha_discrete(range=c(0.5,1))+
    ylab("SV Number")+
    #ggtitle(plot_name)+
    theme_bw() + 
    theme(panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=13,angle = 0,colour = "black"),
          axis.text.y = element_text(size=13,colour = "black"),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(color="black", size=10, face="bold.italic",hjust = 0.5),
          strip.text.x = element_text(size = 14, colour = "black"),
          legend.position = "none"
    )
  return(p)
}

df=list(manta_df,meerkat_df,tnscope_df,any2_df,
        manta_meerkat_intersection_df,manta_tnscope_intersection_df,meerkat_tnscope_intersection_df,manta_meerkat_tnscope_intersection_df,
        manta_meerkat_union_df,manta_tnscope_union_df,meerkat_tnscope_union_df,manta_meerkat_tnscope_union_df)
plot_name=c("01.Manta","02.Meerkat","03.TNscope","04.Any 2 of_Manta/Meerkat/TNscope",
            "05.Manta+Meerkat_intersection","06.Manta+TNscope_intersection","07.Meerkat+TNscope_intersection","08.Manta+Meerkat+TNscope_intersection",
            "09.Manta+Meerkat_union","10.Manta+TNscope_union","11.Meerkat+TNscope_union","12.Manta+Meerkat+TNscope_union")



raw=plot_function3(df,plot_name,4)


ggsave(file.path(out_path,paste0("4alg_compare(num)2.png")),plot=raw,width = 11.5,height = 8,units = "in")
ggsave(file.path(out_path,paste0("4alg_compare(num)2.pdf")),plot=raw,width = 11.5,height =8,units = "in")

################################################################################
# plot
plot_function3 <- function(df,plot_name){
  plot_df <- as.data.frame(t(data.frame(matrix(unlist(lapply(df, function(x) apply(x[,4:6],2,sum))),nrow = 3))))
  rownames(plot_df) <- plot_name
  colnames(plot_df) <- c("TP","FP","FN"
                         # "TTP","FFP","TFN"
                         )
  plot_df$name <- plot_name
  plot_df_new <- data.frame()
  new_matrix <- data.frame(plot_df$TP-plot_df$TTP,
                           plot_df$TP-plot_df$TP+plot_df$TTP,
                           plot_df$FP-plot_df$FFP,
                           plot_df$FP-plot_df$FP+plot_df$FFP,
                           plot_df$FN-plot_df$TFN,
                           plot_df$FN-plot_df$FN+plot_df$TFN
  )
  names(new_matrix) <- c("TP-notDetect","TP-Detect","FP-notDetect","FP-Detect","FN-notDetect","FN-Detect")
  new_matrix <- c(t(new_matrix))
  new_matrix <- as.data.frame(new_matrix)
  new_matrix$V1 <- rep(c("TP","TP","FP","FP","FN","FN"),nrow(new_matrix)/6)
  new_matrix$V2 <- rep(c("notDetect","Detect"),nrow(new_matrix)/2)
  new_matrix$method <- gsub("_","\n",rep(plot_name,each=6))
  new_matrix$V1 <- factor(new_matrix$V1,levels = c("FP","TP","FN"))
  new_matrix$V2 <- factor(new_matrix$V2,levels = c("notDetect","Detect"))
  new_matrix$gold <- substr(matrix(unlist(strsplit(new_matrix$method,split="\n")),ncol=2,byrow = T)[,1],
                            8,
                            nchar(matrix(unlist(strsplit(new_matrix$method,split="\n")),ncol=2,byrow = T)[,1]))
  new_matrix$screen <- substr(matrix(unlist(strsplit(new_matrix$method,split="\n")),ncol=2,byrow = T)[,2],
                              7,
                              nchar(matrix(unlist(strsplit(new_matrix$method,split="\n")),ncol=2,byrow = T)[,2]))
  
  
  # plot_df <- melt(plot_df)
  # plot_df$variable <- factor(plot_df$variable,levels = c("FP","TP","FN","FFP","TTP","TFN"))
  # print(plot_df)
  p <- ggplot(data=new_matrix,aes(x=V1,y=new_matrix,fill=V1
                                  # alpha=V2
  ))+
    geom_bar(stat = "identity")+
    facet_grid(gold~screen,
    )+
    scale_fill_manual(values=c('#3498DB', "#CB4335","#F2AE10"))+
    #scale_alpha_discrete(range=c(0.5,1))+
    ylab("SV Number")+
    #ggtitle(plot_name)+
    theme_bw() + 
    theme(panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=13,angle = 0,colour = "black"),
          axis.text.y = element_text(size=13,colour = "black"),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(color="black", size=10, face="bold.italic",hjust = 0.5),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.background.x = element_rect(fill="#3498DB"),
          strip.background.y = element_rect(fill="#F2AE10"),
          strip.text.y = element_text(size = 14, colour = "black"),
          strip.placement = "outside",
          
          legend.position = "none"
    )
  return(p)
}




# plot
df = list(
  goldmanta_screenmanta_df,
  goldmanta_screenmeerkat_df,
  goldmanta_screentnscope_df,
  goldmeerkat_screenmanta_df,
  goldmeerkat_screenmeerkat_df,
  goldmeerkat_screentnscope_df,
  goldtnscope_screenmanta_df,
  goldtnscope_screenmeerkat_df,
  goldtnscope_screentnscope_df
)
plot_name=c("01.goldmanta_screenmanta",
            "02.goldmanta_screenmeerkat",
            "03.goldmanta_screentnscope",
            "04.goldmeerkat_screenmanta",
            "05.goldmeerkat_screenmeerkat",
            "06.goldmeerkat_screentnscope",
            "07.goldtnscope_screenmanta",
            "08.goldtnscope_screenmeerkat",
            "09.goldtnscope_screentnscope")
plot_function3(df,plot_name)
raw2=plot_function3(df,plot_name)
ggsave(file.path(out_path,paste0("3alg_compare(num).png")),plot=raw2,width = 11.5,height = 8,units = "in")
ggsave(file.path(out_path,paste0("3alg_compare(num).pdf")),plot=raw2,width = 11.5,height =8,units = "in")


df = list(
  # goldmanta_screenmanta_df,
  goldmanta_screenmeerkat_df,
  goldmanta_screentnscope_df,
  # goldmeerkat_screenmanta_df,
  # goldmeerkat_screenmeerkat_df,
  goldmeerkat_screentnscope_df
  # goldtnscope_screenmanta_df,
  # goldtnscope_screenmeerkat_df,
  # goldtnscope_screentnscope_df
)
plot_name=c(
  # "01.goldmanta_screenmanta",
            "02.goldmanta_screenmeerkat",
            "03.goldmanta_screentnscope",
            # "04.goldmeerkat_screenmanta",
            # "05.goldmeerkat_screenmeerkat",
            "06.goldmeerkat_screentnscope"
            # "07.goldtnscope_screenmanta",
            # "08.goldtnscope_screenmeerkat",
            # "09.goldtnscope_screentnscope"
            )
plot_function3(df,plot_name)
ggsave(file.path(out_path,paste0("3alg_compare_simple(num).png")),width = 6,height = 4,units = "in")
ggsave(file.path(out_path,paste0("3alg_compare_simple(num).pdf")),width = 6,height =4,units = "in")
