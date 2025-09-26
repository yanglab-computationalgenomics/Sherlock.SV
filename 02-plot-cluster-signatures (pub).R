library(ComplexHeatmap)
library(dplyr)
library(reshape)
library(readxl)
library(data.table)
############################### functions ##################################
plot_function <- function(signature.num,cluster.num,method,mysort,hq=FALSE){
  cluster_path=file.path(sherlock.path,
                         paste0("scratch/07-cluster/1217_samples_hg38/Sherlock_SV49_simple/",signature.num,"signature_myannotation/",method,"/sig_cluster",cluster.num,".tsv")
  )
  plotpath=file.path(sherlock.path,
                     paste0("scratch/07-cluster/1217_samples_hg38/Sherlock_SV49_simple_plot/",signature.num,"signature_myannotation/",method)
  )
  
  ### make dir
  if (!dir.exists(plotpath)){
    dir.create(plotpath,showWarnings = F,recursive = T)
  }
  
  ### read files
  cluster <- read.delim(cluster_path,row.names = 1)
  
  ### add samples that not in the cluster
  cluster.tmp <- merge(cluster %>% mutate(sample=rownames(.)),
                   sample_info %>% dplyr::select(Tumor_Barcode),
                   by.x = "sample",by.y = "Tumor_Barcode",
                   all.y = T
                   )
  cluster.tmp[cluster.tmp$sample %in% rownames(cluster)==FALSE,
              colnames(cluster)[grepl("consensusClass|order|col",colnames(cluster))==FALSE]] <- 0
  cluster.tmp[cluster.tmp$sample %in% rownames(cluster)==FALSE,
              colnames(cluster)[grepl("consensusClass",colnames(cluster))]] <- 0
  cluster.tmp[cluster.tmp$sample %in% rownames(cluster)==FALSE,
              colnames(cluster)[grepl("col",colnames(cluster))]] <- "#e5e5e5"
  cluster.tmp[cluster.tmp$sample %in% rownames(cluster)==FALSE,
              colnames(cluster)[grepl("order",colnames(cluster))]] <- seq(max(cluster$order)+1,nrow(sample_info),1)
  rownames(cluster.tmp) <- cluster.tmp$sample
  cluster.tmp <- cluster.tmp %>% dplyr::select(colnames(cluster))
  cluster <- cluster.tmp
  
  sample_cgr <- complex_sv_with_signature %>% group_by(sample_id,CGR_signature,cluster_id) %>% count() %>% group_by(sample_id,CGR_signature) %>% count()
  sample_cgr_number_wide <- cast(sample_cgr,sample_id ~ CGR_signature)
  sample_cgr_number_wide[is.na(sample_cgr_number_wide)] <- 0
  
  unclustered_cmplx$FULL_SV_TYPE <- gsub("chromoplexy_bal|chromoplexy_del","chromoplexy",unclustered_cmplx$FULL_SV_TYPE)
  unclustered_cmplx <- unclustered_cmplx %>% group_by(SAMPLE,FULL_SV_TYPE,CLUSTER_ID) %>% count() %>% group_by(SAMPLE,FULL_SV_TYPE) %>% count()
  unclustered_cmplx_wide <- cast(unclustered_cmplx,SAMPLE ~ FULL_SV_TYPE)
  unclustered_cmplx_wide[is.na(unclustered_cmplx_wide)] <- 0
  
  ############################### plot file ##################################
  ### merge files
  plot_df <- cluster
  plot_df$sum <- apply(plot_df[,1:as.numeric(signature.num)],1,sum)
  # add sample name
  plot_df$sample <- rownames(plot_df)
  if (hq==TRUE){
    plot_df <- plot_df %>% filter(sample %in% hqsample)
  }
  
  # add sample info
  plot_df <- merge(plot_df,sample_info,by.x = "sample",by.y="Tumor_Barcode",all.x = T)
  # add mygenes
  plot_df$TP53 <- ifelse(plot_df$sample %in% as.vector(unlist(tp53_mutated_samples)),"Mut","WT")
  plot_df$KRAS <- ifelse(plot_df$sample %in% as.vector(unlist(kras_g12c_samples)),"G12C","WT")
  plot_df$KRAS <- ifelse(plot_df$sample %in% as.vector(unlist(kras_g12v_samples)),"G12V",plot_df$KRAS)
  plot_df$KRAS <- ifelse(plot_df$sample %in% as.vector(unlist(kras_g12d_samples)),"G12D",plot_df$KRAS)
  plot_df$KRAS <- ifelse(plot_df$sample %in% as.vector(unlist(kras_q61h_samples)),"Q61H",plot_df$KRAS)
  plot_df$KRAS <- ifelse(plot_df$sample %in% as.vector(unlist(kras_others_samples)),"Other_mut",plot_df$KRAS)
  plot_df$EGFR <- ifelse(plot_df$sample %in% as.vector(unlist(egfr_mutated_E479_A483del_samples)),"E479_A483del","WT")
  plot_df$EGFR <- ifelse(plot_df$sample %in% as.vector(unlist(egfr_mutated_L591R_samples)),"L591R",plot_df$EGFR)
  plot_df$EGFR <- ifelse(plot_df$sample %in% as.vector(unlist(egfr_mutated_others_samples)),"Other_mut",plot_df$EGFR)
  # add fusion
  plot_df <- merge(plot_df,fusion_df_long %>% setnames(.,old="fusion",new="fusion_detail"),
                             by.x = "sample",by.y = "Sample",all.x = T)
  plot_df$fusion_detail <- ifelse(is.na(plot_df$fusion_detail),"WT",plot_df$fusion_detail)
  # add cgr
  plot_df <- merge(plot_df,
                   sample_cgr_number_wide,
                   by.x = "sample",
                   by.y = "sample_id",
                   all.x = T,
                   sort = F)
  # add unclustered complex
  plot_df <- merge(plot_df,
                   unclustered_cmplx_wide,
                   by.x = "sample",
                   by.y = "SAMPLE",
                   all.x = T,
                   sort = F)
  #  sort
  if (mysort=="cluster"){
    plot_df <- plot_df %>% mutate(consensusClass=factor(consensusClass,levels = c(4,1,2,3,seq(5,as.numeric(signature.num)),0)))  %>%arrange(consensusClass,desc(sum))
    #plot_df <- plot_df %>% mutate(consensusClass=factor(consensusClass,levels = c(3,1,2,0)))  %>%arrange(consensusClass,desc(sum))
    #plot_df <- plot_df %>% mutate(consensusClass=factor(consensusClass,levels = c(5,1,2,4,3,0)))  %>%arrange(consensusClass,desc(sum))
    #plot_df <- plot_df %>% mutate(consensusClass=factor(consensusClass,levels = c(6,1,2,4,5,3,0)))  %>%arrange(consensusClass,desc(sum))
  }
  if (mysort=="smoking"){
    plot_df <- plot_df[order(plot_df$Smoking,plot_df$consensusClass,plot_df$sum,decreasing =T),]
  }
  # Population
  plot_df <- plot_df %>% mutate(Assigned_Population=ifelse(Assigned_Population %in% c("AFR","AMR or Mixed"),"Others",Assigned_Population))
  # SP_Group_new
  plot_df <- plot_df %>% mutate(SP_Group_new=ifelse(Assigned_Population=="EUR" & Smoking=="Smoker","S_U",
                                        ifelse(Assigned_Population=="EUR" & Smoking=="Non-Smoker", "N_U",
                                               ifelse(Assigned_Population=="EAS" & Smoking=="Non-Smoker", "N_A","Others"))
  ))
  
  ############################### colors ##################################
  # mutation color
  mutation_color <- c("L591R"="#1B9E77","E479_A483del"="#D95F02","Other_mut"="#7570B3",
                      "G12C"="#DC0073","G12V"="#4A6FA5","G12D"="#23F0C7","Q61H"="#F9B3D1",
                      "Mut"="black", "WT"="white")
  # group color
  population_color <- c("Others"="#7C7C7C", "EAS"="#66C2A5", "EUR"="#FFD92F")
  # smoking color
  smoking_color <- c("Non-Smoker"="#00a800","Smoker"="#FF6EC7","Unknown"="#7C7C7C")
  # gender color
  gender_color <-  c("Male"="#8DA0CB", "Female" = "#FC8D62")
  # fusion color
  fusion_color <- c("ALK"="#ff7043","ROS1"="#cddc39","RET"="#26c6da","MET"="#80cbc4",
                    "FGFR1"="#6C5B7B","FGFR2"="#6C5B7B","FGFR3"="#6C5B7B","FGFR4"="#6C5B7B",
                    "NTRK1"="#CAD7B2","NTRK2"="#CAD7B2","NTRK3"="#CAD7B2","WT"="white")
  # histology color
  histology_color <- c("Adenocarcinoma" = "#4CAF50", 
                       "Squamous cell carcinoma" = "#2196F3" ,
                       "Carcinoid tumor" = "#FF9800" , 
                       "Adenosquamous carcinoma" = "#9C27B0"  ,
                       "Others" = "#F0F0F0" )
  # signature col
  sig_name_list <- list(
    "8"= c("B","A","D","C","H","F","E","G"),
    "9"= c("B","A","C","H","I","F","E","D","G")
  )
  sig_col_list <- list(
    "8" = c("#D6EAF8","#85C1E9","#5DADE2","#F1C40F","#FADBD8","#F1948A","#14CCCC","#A569BD"),
    "9" = c("#D6EAF8","#85C1E9","#5DADE2", "#3498DB","#FADBD8","#F1948A","#14CCCC","#F1C40F","#A569BD")
  )
  sig_name <-  unlist(sig_name_list[signature.num])
  sig_col <- unlist(sig_col_list[signature.num])
  names(sig_col) <- sig_name
  # sig_name_full <- paste0("SBS49",sig_name)
  sig_name_full <- c("Del1","Del2","Del3","Large.intra","TD1","TD2","Fb.inv","Tra")
  
  # signature col  legend
  sig_name_list_legend <- list(
    "8"= c("Del1","Del2","Del3","Large intra","TD1","TD2","Fb inv","Tra"),
    "9"= c("Del1","Del2","Del3","Del4+TD3+recip inv","TD1","TD2","Fb inv","Del4+TD3+unbal inv","Tra")
  )
  sig_col_list_legend <- list(
    "8" = c("#D6EAF8","#85C1E9","#5DADE2","#F1C40F","#FADBD8","#F1948A","#14CCCC","#A569BD"),
    "9" = c("#D6EAF8","#85C1E9","#5DADE2", "#3498DB","#FADBD8","#F1948A","#14CCCC","#F1C40F","#A569BD")
  )
  sig_name_legend <-  unlist(sig_name_list_legend[signature.num])
  sig_col_legend <- unlist(sig_col_list_legend[signature.num])
  names(sig_col_legend) <- sig_name_legend

  
  
  
  # cluster col
  class_col_df <- cluster %>% group_by(consensusClass,col) %>% dplyr::count()
  consensusClass_col <- class_col_df$col
  # consensusClass_col <- c("#F1C40F","#85C1E9","#F1948A","#FADBD8")
  names(consensusClass_col) <- class_col_df$consensusClass

  
  
  complexheatmap_function <- function(plot_df,savepath){
    #plot_df <- plot_df %>% filter(smoking %in% c("Smoker_SBS4","NonSmoker_NonSBS4"))
    column_ha <- HeatmapAnnotation(
     
      "Cluster"=anno_simple(plot_df[,"consensusClass"] %>% as.character(),
                            simple_anno_size = unit(4, "mm"),col = consensusClass_col),
      "Simple SV freq." = anno_barplot(
        t(apply(plot_df[, sig_name_full], 1, as.numeric)),
        height = unit(8, "mm"),
        width = 1,
        gp = gpar(col = sig_col,
                  boder = NA,
                  fill = sig_col,
                  just =  "left"),
        #show_annotation_name = TRUE,
        axis_param = list(gp=gpar(fontsize=6))
      ),
      
      "Simple SV signature pct." = anno_barplot(
        t(apply(plot_df[, sig_name_full], 1, function(x) as.numeric(x)/sum(x))) %>% {replace(., is.na(.), 0)},
        height = unit(8, "mm"),
        width = 1,
        gp = gpar(col = sig_col,
                  boder = sig_col,
                  fill = sig_col),
        #show_annotation_name = TRUE,
        axis_param = list(gp=gpar(fontsize=6))
      ),
      
      "Del1" = anno_barplot(plot_df[, "Del1"],height = unit(4, "mm"),gp = gpar(col = sig_col["B"],boder = NA,fill = sig_col["B"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      "Del2" = anno_barplot(plot_df[, "Del2"],height = unit(4, "mm"),gp = gpar(col = sig_col["A"],boder = NA,fill = sig_col["A"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      "Del3" = anno_barplot(plot_df[, "Del3"],height = unit(4, "mm"),gp = gpar(col = sig_col["D"],boder = NA,fill = sig_col["D"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      "TD1" = anno_barplot(plot_df[, "TD1"],height = unit(4, "mm"),gp = gpar(col = sig_col["H"],boder = NA,fill = sig_col["H"]),axis_param = list(gp=gpar(fontsize=6),at = c(150,300),labels = c(150,300))),
      "TD2" = anno_barplot(plot_df[, "TD2"],height = unit(4, "mm"),gp = gpar(col = sig_col["F"],boder = NA,fill = sig_col["F"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      "Fb inv" = anno_barplot(plot_df[, "Fb.inv"],height = unit(4, "mm"),gp = gpar(col = sig_col["E"],boder = NA,fill = sig_col["E"]),axis_param = list(gp=gpar(fontsize=6),at = c(50,100),labels = c(50,100))),
      "Large intra" = anno_barplot(plot_df[, "Large.intra"],height = unit(4, "mm"),gp = gpar(col = sig_col["C"],boder = NA,fill = sig_col["C"]),axis_param = list(gp=gpar(fontsize=6),at = c(50),labels = c(50))),
      "Tra" = anno_barplot(plot_df[, "Tra"],height = unit(4, "mm"),gp = gpar(col = sig_col["G"],boder = NA,fill = sig_col["G"]),axis_param = list(gp=gpar(fontsize=6),at = c(40),labels = c(40))),
      
      
      # "Del1" = anno_barplot(plot_df[, "SBS49B"],height = unit(4, "mm"),gp = gpar(col = sig_col["B"],boder = NA,fill = sig_col["B"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      # "Del2" = anno_barplot(plot_df[, "SBS49A"],height = unit(4, "mm"),gp = gpar(col = sig_col["A"],boder = NA,fill = sig_col["A"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      # "Del3" = anno_barplot(plot_df[, "SBS49D"],height = unit(4, "mm"),gp = gpar(col = sig_col["D"],boder = NA,fill = sig_col["D"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      # "TD1" = anno_barplot(plot_df[, "SBS49H"],height = unit(4, "mm"),gp = gpar(col = sig_col["H"],boder = NA,fill = sig_col["H"]),axis_param = list(gp=gpar(fontsize=6),at = c(150,300),labels = c(150,300))),
      # "TD2" = anno_barplot(plot_df[, "SBS49F"],height = unit(4, "mm"),gp = gpar(col = sig_col["F"],boder = NA,fill = sig_col["F"]),axis_param = list(gp=gpar(fontsize=6),at = c(100,200),labels = c(100,200))),
      # "Fb inv" = anno_barplot(plot_df[, "SBS49E"],height = unit(4, "mm"),gp = gpar(col = sig_col["E"],boder = NA,fill = sig_col["E"]),axis_param = list(gp=gpar(fontsize=6),at = c(50,100),labels = c(50,100))),
      # "Large intra" = anno_barplot(plot_df[, "SBS49C"],height = unit(4, "mm"),gp = gpar(col = sig_col["C"],boder = NA,fill = sig_col["C"]),axis_param = list(gp=gpar(fontsize=6),at = c(50),labels = c(50))),
      # "Tra" = anno_barplot(plot_df[, "SBS49G"],height = unit(4, "mm"),gp = gpar(col = sig_col["G"],boder = NA,fill = sig_col["G"]),axis_param = list(gp=gpar(fontsize=6),at = c(40),labels = c(40))),

      
      # "Del1" = anno_barplot(plot_df[, "SBS49B"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["B"],boder = NA,fill = sig_col["B"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Del2" = anno_barplot(plot_df[, "SBS49A"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["A"],boder = NA,fill = sig_col["A"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Del3" = anno_barplot(plot_df[, "SBS49C"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["C"],boder = NA,fill = sig_col["C"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Del4+TD3+recip inv" = anno_barplot(plot_df[, "SBS49H"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["H"],boder = NA,fill = sig_col["H"]),axis_param = list(gp=gpar(fontsize=15))),
      # "TD1" = anno_barplot(plot_df[, "SBS49I"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["I"],boder = NA,fill = sig_col["I"]),axis_param = list(gp=gpar(fontsize=15))),
      # "TD2" = anno_barplot(plot_df[, "SBS49F"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["F"],boder = NA,fill = sig_col["F"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Fb inv" = anno_barplot(plot_df[, "SBS49E"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["E"],boder = NA,fill = sig_col["E"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Del4+TD3+unbal inv" = anno_barplot(plot_df[, "SBS49D"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["D"],boder = NA,fill = sig_col["D"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Tra" = anno_barplot(plot_df[, "SBS49G"],height = unit(1.5, "cm"),gp = gpar(col = sig_col["G"],boder = NA,fill = sig_col["G"]),axis_param = list(gp=gpar(fontsize=15))),
      # 
      # "1 ecDNA" = anno_barplot(plot_df[, "1 ecDNA/double minutes"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["CGR_ecDNA"],boder = NA,fill = fusion_cgr_col["CGR_ecDNA"]),axis_param = list(gp=gpar(fontsize=15))),
      # "2 Chr bridge" = anno_barplot(plot_df[, "2 BFB cycles/chromatin bridge"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["CGR_chr_bridge"],boder = NA,fill = fusion_cgr_col["CGR_chr_bridge"]),axis_param = list(gp=gpar(fontsize=15))),
      # "3 Large loss" = anno_barplot(plot_df[, "3 Large loss"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["CGR_large_loss"],boder = NA,fill = fusion_cgr_col["CGR_large_loss"]),axis_param = list(gp=gpar(fontsize=15))),
      # "4 Micronuclei" = anno_barplot(plot_df[, "4 Micronuclei"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["CGR_micronuclei"],boder = NA,fill = fusion_cgr_col["CGR_micronuclei"]),axis_param = list(gp=gpar(fontsize=15))),
      # "5 Large gain" = anno_barplot(plot_df[, "5 Large gain"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["CGR_large_gain"],boder = NA,fill = fusion_cgr_col["CGR_large_gain"]),axis_param = list(gp=gpar(fontsize=15))),
      # "6 Hourglass" = anno_barplot(plot_df[, "6 Hourglass"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["CGR_hourglass"],boder = NA,fill = fusion_cgr_col["CGR_hourglass"]),axis_param = list(gp=gpar(fontsize=15))),
      
      # "Chromoplexy" = anno_barplot(plot_df[, "chromoplexy"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["chromoplexy"],boder = NA,fill = fusion_cgr_col["chromoplexy"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Cycle templated ins" = anno_barplot(plot_df[, "cycle_templated_ins"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["cycle_templated_ins"],boder = NA,fill = fusion_cgr_col["cycle_templated_ins"]),axis_param = list(gp=gpar(fontsize=15))),
      # "Complex unclear" = anno_barplot(plot_df[, "complex_unclear"],height = unit(1.5, "cm"),gp = gpar(col = fusion_cgr_col["complex_unclear"],boder = NA,fill = fusion_cgr_col["complex_unclear"]),axis_param = list(gp=gpar(fontsize=15))),
      
      "Population" = anno_simple(plot_df[,"Assigned_Population"],simple_anno_size = unit(3, "mm"),col = population_color),
      "Smoking" = anno_simple(plot_df[,"Smoking"],simple_anno_size = unit(3, "mm"),col = smoking_color),

      "TP53" = anno_simple(plot_df[,"TP53"],simple_anno_size = unit(3, "mm"),col = mutation_color),
      "EGFR" = anno_simple(plot_df[,"EGFR"],simple_anno_size = unit(3, "mm"),col = mutation_color),
      "KRAS" = anno_simple(plot_df[,"KRAS"],simple_anno_size = unit(3, "mm"),col = mutation_color),

      # "Fusion complex SV" = anno_simple(plot_df[,"fusion_CGR"],simple_anno_size = unit(1.5, "cm"),col = fusion_cgr_col),
      # "Fusion simple SV" = anno_simple(plot_df[,"fusion_simplesv"],simple_anno_size = unit(1.5, "cm"),col = fusion_simplesv_col),
      "Fusion" = anno_simple(plot_df[,"fusion_detail"],simple_anno_size = unit(3, "mm"),col = fusion_color),
      
      # "SCNA" = anno_simple(plot_df[,"scna_group"],simple_anno_size = unit(1.5, "cm"),col = scna_color),
      # "WGD" = anno_simple(plot_df[,"wgd"],simple_anno_size = unit(1.5, "cm"),col = wgd_color),
      "Sex" = anno_simple(plot_df[,"Gender"],simple_anno_size = unit(3, "mm"),col = gender_color),
      "Tumor type" = anno_simple(plot_df[,"Histology"],simple_anno_size = unit(3, "mm"),col = histology_color),

      annotation_name_gp = gpar(fontsize = 6),
      # column_names_gp = gpar(fontsize = 6),
      # row_names_gp = gpar(fontsize = 6),
      annotation_name_rot = 0,
      annotation_name_side = "left",
      annotation_name_align = F,
      gap = unit(c(1.5,
                   # 1.5,
                   1.5,
                   1.5,
                   rep(0.45,as.numeric(signature.num)-1),1.5,
                   # rep(4.5,8),15, # complex SV
                   rep(0.45,7)
      ),
      "mm")
    )
    
    lgd_list <- packLegend(
      Legend(labels = names(sig_col_legend),
             title = "Simple SV signatures",
             #type = "lines",
             # pch = 50,
             #background = "white",
             legend_gp = gpar(col=sig_col_legend,
                              fill=sig_col_legend,
                              pch=50
                              # lty = 2,
                              # lwd = 6
             ),
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6,fontface="plain"),
             grid_height = unit(3,"mm"),
             grid_width = unit(3,"mm"),
             title_position = "leftcenter",
             nrow = 1
      ),
      Legend(labels = names(smoking_color),
             title = "Smoking",
             #type = "lines",
             # pch = 50,
             #background = "white",
             legend_gp = gpar(col=smoking_color,
                              fill=smoking_color,
                              pch=50
                              # lty = 2,
                              # lwd = 6
             ),
             border = smoking_color,
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6,fontface="plain"),
             grid_height = unit(3,"mm"),
             grid_width = unit(3,"mm"),
             title_position = "leftcenter",
             nrow = 1
             
      ),
      Legend(labels = names(population_color),
             title = "Population",
             #type = "lines",
             # pch = 50,
             #background = "white",
             legend_gp = gpar(col=population_color,
                              fill=population_color,
                              pch=50
                              # lty = 2,
                              # lwd = 6
             ),
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6,fontface="plain"),
             grid_height = unit(3,"mm"),
             grid_width = unit(3,"mm"),
             title_position = "leftcenter",
             nrow = 1

      ),
      Legend(labels = names(mutation_color),
             title = "TP53/KRAS/EGFR",
             #type = "lines",
             # pch = 50,
             #background = "white",
             legend_gp = gpar(color=mutation_color,
                              fill=mutation_color,
                              pch=50
                              # lty = 2,
                              # lwd = 6
             ),
             # border =mutation_color_legend_border[c("L591R","E479_A483del","Other Mut","WT")],
             border =c("black"),
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6,fontface="plain"),
             grid_height = unit(3,"mm"),
             grid_width = unit(3,"mm"),
             title_position = "leftcenter",
             nrow =1
             
      ),
      Legend(labels = names(fusion_color),
             title = "Fusion",
             #type = "lines",
             # pch = 50,
             #background = "white",
             legend_gp = gpar(col=fusion_color,
                              fill=fusion_color,
                              pch=50
                              # lty = 2,
                              # lwd = 6
             ),
             border = fusion_color,
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6,fontface="plain"),
             grid_height = unit(3,"mm"),
             grid_width = unit(3,"mm"),
             title_position = "leftcenter",
             nrow = 1

      ),
      Legend(labels = names(gender_color),
             title = "Gender",
             #type = "lines",
             # pch = 50,
             #background = "white",
             legend_gp = gpar(col=gender_color,
                              fill=gender_color,
                              pch=50
                              # lty = 2,
                              # lwd = 6
             ),
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6,fontface="plain"),
             grid_height = unit(3,"mm"),
             grid_width = unit(3,"mm"),
             title_position = "leftcenter",
             nrow = 1

      ),
      Legend(labels = names(histology_color),
             title = "Histology",
             #type = "lines",
             # pch = 50,
             #background = "white",
             legend_gp = gpar(col=histology_color,
                              fill=histology_color,
                              pch=50
                              # lty = 2,
                              # lwd = 6
             ),
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6,fontface="plain"),
             grid_height = unit(3,"mm"),
             grid_width = unit(3,"mm"),
             title_position = "leftcenter",
             nrow = 1
      ),
      column_gap = unit(1, "mm"), row_gap = unit(1, "mm")
     
    )
    
    zero_row_mat <- matrix(nrow = 0, ncol = nrow(plot_df))
    
    pdf(file=savepath,width = 1.080*10,height = 0.540*10)
    
    ht <- Heatmap(zero_row_mat,
                  top_annotation = column_ha,
                  row_gap = unit(5, "mm")
    )
    
    draw(ht,
         annotation_legend_side = "bottom",
         annotation_legend_list = lgd_list
    )
    
    dev.off()
  }

  if (hq==TRUE){
    addname <- ".hq"
  } else {
    addname <- ""
  }
  
  complexheatmap_function(plot_df = plot_df,file.path(plotpath,paste0(cluster.num,"clusters_","sortby_",mysort,"_all",addname,".pdf"))) #all
  complexheatmap_function(plot_df = plot_df[plot_df$SP_Group_new %in% c("S_U"),],file.path(plotpath,paste0(cluster.num,"clusters_","sortby_",mysort,"_smoking",addname,".pdf"))) # smoking population
  complexheatmap_function(plot_df = plot_df[plot_df$SP_Group_new %in% c("N_U","N_A"),],file.path(plotpath,paste0(cluster.num,"clusters_","sortby_",mysort,"_nonsmoking",addname,".pdf"))) # non_smoking population

}



#################################################### prepare #################
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
snv_path <- file.path(sherlock.path,"scratch/06-SNV/sample_snv.txt")
sample_info_path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")
fusion_path  <- file.path(sherlock.path,"scratch/fusion/sherlock1217.driver.fusion.summary.csv")
complex_sv_with_signature_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/complex_sv_with_signature.tsv")
unclustered_cmplx_path <- file.path(scratch.path,"04-shatterseek/1217_samples_hg38/sherlock_nonclustered_cmplx_sv.tsv")
gene_annotation_path <- file.path(sherlock.path,"Data/SNV/gene_annotation.txt")
hq.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/HQ_samples.csv")
fusion.detail.path  <- file.path(sherlock.path,"scratch/fusion/end_annotation/tier1_driver_fusions.txt")
tp53.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/TP53_mutated_samples.tsv")
kras.mut.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_mutated_samples.tsv")
egfr.E479_A483del.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_E479_A483del_samples.tsv")
egfr.L591R.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_L591R_samples.tsv")
egfr.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/EGFR_mutated_others_samples.tsv")
kras.G12C.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12C_samples.txt")
kras.G12V.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12V_samples.txt")
kras.G12D.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pG12D_samples.txt")
kras.Q61H.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_pQ61H_samples.txt")
kras.others.path <- file.path(sherlock.path,"xiaomingsherlock/fusion/oncoplot/KRAS_other_mutations.txt")

fusion_set <- c("ALK","RET","ROS1","MET","NRG1",
                "FGFR1","FGFR2","FGFR3","FGFR4",
                "NTRK1","NTRK2","NTRK3")
snv  <-  read.delim(snv_path,row.names = 1)
sample_info <-  read.delim(sample_info_path)
histology <- read.csv(histology.path)
sample_info <- merge(sample_info,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  dplyr::select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

complex_sv_with_signature <- read.delim(complex_sv_with_signature_path)
unclustered_cmplx <- read.delim(unclustered_cmplx_path)
fusion.detail <- read.delim2(fusion.detail.path)
fusion.detail <- fusion.detail %>%
  dplyr::mutate(
    Fusion_Gene = case_when(
      GENE1 %in% fusion_set ~ GENE1,
      GENE2 %in% fusion_set ~ GENE2,
      TRUE ~ NA_character_
    )
  ) %>% filter(Known_Fusion=="Yes")
snv$sample <- row.names(snv)
gene_annotation <- read.delim(gene_annotation_path)
hqsample <- read.csv(hq.path) %>% pull(Tumor_Barcode)
tp53_mutated_samples=read.delim(tp53.mut.path, header=F, check.names=F, stringsAsFactors=F)
kras_mutated_samples=read.delim(kras.mut.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12c_samples=read.delim(kras.G12C.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12v_samples=read.delim(kras.G12V.path, header=F, check.names=F, stringsAsFactors=F)
kras_g12d_samples=read.delim(kras.G12D.path, header=F, check.names=F, stringsAsFactors=F)
kras_q61h_samples=read.delim(kras.Q61H.path, header=F, check.names=F, stringsAsFactors=F)
kras_others_samples=read.delim(kras.others.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_E479_A483del_samples=read.delim(egfr.E479_A483del.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_L591R_samples=read.delim(egfr.L591R.path, header=F, check.names=F, stringsAsFactors=F)
egfr_mutated_others_samples=read.delim(egfr.others.path, header=F, check.names=F, stringsAsFactors=F)
fusion_df <- data.frame(Sample = sample_info$Tumor_Barcode, stringsAsFactors = FALSE)
fusion_df[fusion_set] <- NA
for (genei in fusion_set) {
  matched_samples <- fusion.detail %>%
    filter(Fusion_Gene == genei) %>%
    pull(SAMPLE)
  fusion_df[fusion_df$Sample %in% matched_samples, genei] <- "fusion"
}
fusion_df_long <- tidyr::gather(fusion_df,fusion,gene,fusion_set) %>%
  filter(gene=="fusion") %>% dplyr::select(Sample,fusion)


# check  function before run, the signature part
# for  (mycluster  in 2:15){
#   for  (mymethod in c("hc_pearson","hc_euclidean","km","pam_euclidean","pam_pearson")[1]){
#     plot_function(signature.num="9",
#                   cluster.num=as.character(mycluster),
#                   method=mymethod
#     )
#   }
# }

# check  function before run, the signature part
for  (mycluster  in 2:15){
  for  (mymethod in c("hc_pearson","hc_euclidean","km","pam_euclidean","pam_pearson")){
    plot_function(signature.num="8",
                  cluster.num=as.character(mycluster),
                  method=mymethod,
                  mysort="cluster",hq=FALSE
    )
  }
}
for  (mycluster  in 2:15){
  for  (mymethod in c("hc_pearson","hc_euclidean","km","pam_euclidean","pam_pearson")){
    plot_function(signature.num="8",
                  cluster.num=as.character(mycluster),
                  method=mymethod,
                  mysort="smoking",hq=FALSE
    )
  }
}

# plot_function(signature.num="8",cluster.num=as.character(3),method="hc_pearson",mysort="cluster",hq=FALSE)
# plot_function(signature.num="8",cluster.num=as.character(4),method="hc_pearson",mysort="cluster",hq=FALSE)
# plot_function(signature.num="8",cluster.num=as.character(5),method="hc_pearson",mysort="cluster",hq=FALSE)
# plot_function(signature.num="8",cluster.num=as.character(6),method="hc_pearson",mysort="cluster",hq=FALSE)

plot_function(signature.num="8",cluster.num=as.character(4),method="hc_pearson",mysort="cluster",hq=FALSE)
plot_function(signature.num="8",cluster.num=as.character(5),method="hc_pearson",mysort="cluster",hq=FALSE)
plot_function(signature.num="8",cluster.num=as.character(4),method="hc_pearson",mysort="cluster",hq=TRUE)
