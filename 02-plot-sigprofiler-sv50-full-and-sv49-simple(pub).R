# library
library(reshape)
library(magrittr)
library(ggplot2)
library(dplyr)
library(ggh4x)


plot_sigprofiler_function <- function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_8_Signatures/Signatures/SBS50_S8_Signatures.txt",
                                      output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                                      class_num="50",
                                      color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                                      figure_name="SV50_full_8_Signatures",
                                      signaturenum="8",
                                      signame_list=NA,
                                      order=NA){
  if (!dir.exists(output_signature_dir)) {
    dir.create(output_signature_dir, recursive = TRUE)
  }
  
  ####################### annotate row ##############################
  if  (class_num=="49"){
    yinter <- c(2,7,10,11,12,30,31)+0.5
  }
  if  (class_num=="50"){
    yinter <- c(4,7,10,11,12,30,31)+0.5
  }
  ############################ READ FILES ##########################
  phat <- read.delim(phat_path,row.names = 1)
  colnames(phat) <- gsub(paste0("SBS",class_num),"",colnames(phat))
  
  
  ############################# plot signatures ########################
  phat$group  <- row.names(phat)
  # phat$group_fullname <- ifelse(phat$group %in% names(del_dup_rename),
  #                               del_dup_rename[phat$group],
  #                               phat$group)
  phat <- melt(phat, value.name="percent", variable.name="sig", na.rm=TRUE)
  phat_matrix_color<-read.csv(color_path,header=TRUE, fileEncoding="UTF-8-BOM")
  phat_color <- merge(phat,phat_matrix_color,by.x="group",by.y="sig")
  signew=dplyr::group_by(phat_color, group) %>% dplyr::mutate(percentage = value/sum(value))
  signew_data=as.data.frame(signew)
  
  signew_data_order <- signew_data[order(signew_data$order,signew_data$group), ]
  signew_data_order$svtype = factor(signew_data_order$svtype, levels = c("DEL","DUP","INV","TRA","C_A","C_G","C_T","T_A","T_C","T_G"))
  
  
  signew_data_order$group_fullname = factor(
    signew_data_order$group_fullname,
    levels = rev(unique(phat_matrix_color$group_fullname))
  )
  
  ################### plot signatures ###################
  if  (!is.na(signame_list)){
    names(signame_list) <- order
    signew_data_order$signame <- signame_list[as.character(signew_data_order$variable)]
    signew_data_order$signame = factor(signew_data_order$signame,levels = rev(signame_list))
  } else{
    signew_data_order$signame <- signew_data_order$variable
  }
  
  
  
  signature_plot <- ggplot(data = signew_data_order, aes(y=group_fullname,x=percentage,fill=group,
                                                         # color=group,
                                                         alpha=1)) + 
    geom_col(width = 0.9) +
    geom_hline(yintercept = yinter,
               color="grey")+
    scale_fill_manual(values=as.vector(signew_data_order$sigcolor),
                      breaks = as.vector(signew_data_order$group)) +
    # scale_color_manual(values=as.vector(signew_data_order$sigcolor),
    #                    breaks = as.vector(signew_data_order$group)) +
    scale_x_continuous(breaks = c(0,0.5,1),
                       labels = c(0,0.5,1)) +
    facet_wrap2(~signame, 
                nrow = 1,
                scales = "free_x",
                strip = strip_vanilla(clip = "off")
                # labeller=labeller(variable = ss_labels)
    ) +
    ggtitle(figure_name)+
    theme(axis.text.x = element_text(size = 5,color = "black"),
          axis.text.y = element_text(size = 6,color = "black",vjust = 0.5,hjust = 0),
          axis.title = element_blank(),
          axis.ticks =  element_line(linewidth = 0.234),
          panel.background = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(size=6, angle = 45, colour = "black",hjust = 0,vjust = 0),
          strip.background = element_blank(),
          strip.placement = "outside",
          plot.title = element_text(size=6,color = "black"),
          plot.margin=grid::unit(c(2,5,2,2), "mm")
    )
  
  signature_plot
  ggsave(file.path(output_signature_dir,paste0(figure_name,".pdf")), height = 150, width =30+6*as.numeric(signaturenum),units = c("mm"))
}

##################### RUN ######################
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_8_Signatures/Signatures/SBS50_S8_Signatures.txt",
                                      output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                                      class_num="50",
                                      color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                                      figure_name="SV50_full_8_Signatures",
                                      signaturenum="8",
                                      signame_list=c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large mixed","Tra"),
                                      order=c("B","A","D","H","F","E","C","G"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_9_Signatures/Signatures/SBS50_S9_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_9_Signatures",
                          signaturenum="9",
                          signame_list=c("Del1","Del2","Del3","Del4","TD1","TD2","Fb inv","Large mixed","Tra"),
                          order=c("B","A","E","C","I","G","F","D","H"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_10_Signatures/Signatures/SBS50_S10_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_10_Signatures",
                          signaturenum="10",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","TD1","TD2","Fb inv","Large mixed","Tra"),
                          order=c("C","A","B","E","I","J","H","D","F","G"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_11_Signatures/Signatures/SBS50_S11_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_11_Signatures",
                          signaturenum="11",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","TD1","TD2","TD3","Fb inv","Large mixed","Tra"),
                          order=c("D","A","B","E","I","J","K","H","C","G","F"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_12_Signatures/Signatures/SBS50_S12_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_12_Signatures",
                          signaturenum="12",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","TD1","TD2","TD3","Fb inv","Unbal inv","Tra"),
                          order=c("E","A","C","F","I","L","K","J","G","B","H","D"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_13_Signatures/Signatures/SBS50_S13_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_13_Signatures",
                          signaturenum="13",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","TD1","TD2","TD3","TD4","Fb inv","Unbal inv","Tra"),
                          order=c("E","C","D","G","J","M","L","K","F","I","A","H","B"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_14_Signatures/Signatures/SBS50_S14_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_14_Signatures",
                          signaturenum="14",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","TD1","TD2","TD3","TD4","Fragile site","Fb inv","Unbal inv","Tra"),
                          order=c("D","C","F","H","M","N","L","K","G","J","E","A","I","B"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_15_Signatures/Signatures/SBS50_S15_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_15_Signatures",
                          signaturenum="15",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","Del7","TD1","TD2","TD3","TD4","Fragile site","Fb inv","Unbal inv","Tra"),
                          order=c("C","E","G","D","H","M","O","K","L","F","N","J","A","I","B"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV50_full/CH50/All_Solutions/SBS50_16_Signatures/Signatures/SBS50_S16_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV50_full",
                          class_num="50",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/50sv_full_matrix_color_svtype.csv",
                          figure_name="SV50_full_16_Signatures",
                          signaturenum="16",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","Del7","TD1","TD2","TD3","TD4","Fragile site","Fb inv","Unbal inv1","Unbal inv2","Tra"),
                          order=c("C","D","F","E","H","N","O","M","K","G","P","J","A","L","I","B"))


plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_7_Signatures/Signatures/SBS49_S7_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_7_Signatures",
                          signaturenum="7",
                          signame_list=c("Del1","Del2","Del3","TD1","TD2","Fb inv","Tra"),
                          order=c("B","A","C","G","E","D","F"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_8_Signatures/Signatures/SBS49_S8_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_8_Signatures",
                          signaturenum="8",
                          signame_list=c("Del1","Del2","Del3","TD1","TD2","Fb inv","Large mixed","Tra"),
                          order=c("B","A","D","H","F","E","C","G"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_9_Signatures/Signatures/SBS49_S9_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_9_Signatures",
                          signaturenum="9",
                          signame_list=c("Del1","Del2","Del3","Del4","TD1","TD2","Fb inv","Large mixed","Tra"),
                          order=c("B","A","C","H","I","F","E","D","G"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_10_Signatures/Signatures/SBS49_S10_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_10_Signatures",
                          signaturenum="10",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","TD1","TD2","Fb inv","Large mixed","Tra"),
                          order=c("C","B","A","D","I","J","G","E","F","H"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_11_Signatures/Signatures/SBS49_S11_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_11_Signatures",
                          signaturenum="11",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","TD1","TD2","TD3","Fb inv","Large mixed","Tra"),
                          order=c("D","A","B","G","K","J","I","E","C","H","F"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_12_Signatures/Signatures/SBS49_S12_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_12_Signatures",
                          signaturenum="12",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","TD1","TD2","TD3","Fb inv","Large mixed","Tra"),
                          order=c("E","A","B","G","J","L","K","I","F","C","H","D"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_13_Signatures/Signatures/SBS49_S13_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_13_Signatures",
                          signaturenum="13",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","TD1","TD2","TD3","TD4","Fb inv","Large mixed","Tra"),
                          order=c("E","A","C","G","J","M","K","I","F","L","B","H","D"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_14_Signatures/Signatures/SBS49_S14_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_14_Signatures",
                          signaturenum="14",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","TD1","TD2","TD3","TD4","Fragile site","Fb inv","Large mixed","Tra"),
                          order=c("E","C","D","H","M","N","K","J","F","L","G","A","I","B"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_15_Signatures/Signatures/SBS49_S15_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_15_Signatures",
                          signaturenum="15",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","Del7","TD1","TD2","TD3","TD4","Fragile site","Fb inv","Large mixed","Tra"),
                          order=c("C","D","F","E","H","M","N","K","L","G","O","J","A","I","B"))
plot_sigprofiler_function(phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_16_Signatures/Signatures/SBS49_S16_Signatures.txt",
                          output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler-plot/1217_samples_hg38/Sherlock_SV49_simple",
                          class_num="49",
                          color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/49sv_simple_matrix_color_svtype.csv",
                          figure_name="SV49_simple_16_Signatures",
                          signaturenum="16",
                          signame_list=c("Del1","Del2","Del3","Del4","Del5","Del6","Del7","Del8","TD1","TD2","TD3","TD4","Fragile site","Fb inv","Large mixed","Tra"),
                          order=c("C","F","D","G","J","H","N","O","K","L","E","M","A","I","P","B"))
