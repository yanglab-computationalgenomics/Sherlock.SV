library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggpubr)
library(broom)
library(purrr)
library(data.table)
library(rstatix)
library(car)

# Function: if mutated, 0 if wildtype, NA if unknown
mutation_score <- function(x) ifelse(is.na(x) | x == "NA", 0, 1)

forest.plot <- function(model=model.snv,mytitle){
  # Extract and tidy coefficients
  myfit_df_all <- tidy(model, conf.int = TRUE, exponentiate = TRUE)
  
  # Drop intercept
  myfit_df_all <- myfit_df_all %>%
    filter(term != "(Intercept)") %>%
    rename(variable = term, OR = estimate, `2.5 %` = conf.low, `97.5 %` = conf.high)
  
  # Add stars based on p-value
  myfit_df_all <- myfit_df_all %>%
    mutate(
      star = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      ORstar = ifelse(grepl("\\*", star), paste0(signif(OR, 2), " ", star), NA),
      OR_color = case_when(
        OR > 1 & star != "" ~ "red",
        OR < 1 & star != "" ~ "blue",
        TRUE ~ "black"
      )
    )
  variable_fullname <- c("European VS Asian",
                         "Male VS Female",
                         "Smoker VS Never smoker",
                         "Unknown VS Never smoker",
                         "Adenosquamous carcinoma",
                         "Carcinoid tumor",
                         "Squamous cell carcinoma",
                         "Others",
                         "Tumor purity",
                         "Whole-genome duplication",
                         "Mezzo-forte",
                         "Forte",
                         "StageII","StageIII","StageIV",
                         "TP53 mutant VS Wildtype","EGFR mutant VS Wildtype","KRAS mutant VS Wildtype",
                         "Fusion VS Wildtype",
                         "Mutation load",
                         "Indel load",
                         "SV burden"
  )
  names(variable_fullname) <- c("Assigned_PopulationEUR",
                                "GenderMale",
                                "SmokingSmoker",
                                "SmokingUnknown",
                                "HistologyAdenosquamous carcinoma",
                                "HistologyCarcinoid tumor",
                                "HistologySquamous cell carcinoma",
                                "HistologyOthers",
                                "purity",
                                "wgdWGD",
                                "scna_groupMezzo-forte",
                                "scna_groupForte",
                                "StageII","StageIII","StageIV",
                                "TP53_statusMut","EGFR_statusMut","KRAS_statusMut",
                                "Fusion_statusMut",
                                "log10(mut_load + 1)",
                                "log10(indel_load + 1)",
                                "svburden"
  )
  myfit_df_all$variablefullname <- variable_fullname[as.character(myfit_df_all$variable)]
  # Order factors
  myfit_df_all <- myfit_df_all %>%
    mutate(variable = factor(variable, levels = rev(unique(variable))))
  
  # Plot
  p <- ggplot(data = myfit_df_all,
              aes(y = variable, x = OR, xmin = `2.5 %`, xmax = `97.5 %`)) +
    geom_point(size = 0.5) +
    geom_errorbarh(height = 0.15, size = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", lwd = 0.3) +
    geom_text(aes(color = OR_color,
                  x = OR, #ifelse(OR > 1, `2.5 %`, `97.5 %`),
                  y = variable,
                  label = ORstar,
                  hjust = ifelse(OR > 1, 0, 1)),
              vjust = -0.5,
              size = 2) +
    scale_x_log10(
      breaks = c(0.01, 0.1, 0.5,1,2,4,8,16),
      #limits = c(0.01, xmax),
      oob = scales::oob_keep
    ) +
    labs(title = mytitle) +
    scale_y_discrete(labels = myfit_df_all$variablefullname) +
    scale_color_manual(values = c("red", "blue", "black"),breaks = c("red", "blue", "black")) +
    xlab("Odds Ratio") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=6,colour = "black"),
          axis.text.x = element_text(size=6,colour = "black",angle=0,hjust = 0.5),
          axis.text.y = element_text(size=6,colour = "black"),
          strip.text.x = element_text(size=6,colour = "black"),
          axis.ticks = element_line(color = "black",size=0.1),
          axis.ticks.length=unit(0.02,"inch"),
          legend.position = "none",
          plot.title = element_text(size = 6, colour = "black"),
          panel.background =element_blank(),
          panel.border=element_blank(), 
          axis.line=element_line(color = "black",size=0.1))
  return(p)
}


logistic_gene_test <- function(gene, data) {
  # Make a copy of the dataset
  df <- data %>% filter(.data[[gene]] %in% c("NA")==FALSE)
  
  # Binary transform of gene column
  df[[gene]] <- ifelse(is.na(df[[gene]]), 0, 1)
  
  # Choose the formula
  if (gene == "EGFR") {
    formula <- as.formula(paste0("`", gene, "` ~ Smoking + log10(mut_load + 1) + log10(indel_load + 1)"))
  } else if (gene %in% c("fusion", cnv_values, loss_values, gain_values)) {
    formula <- as.formula(paste0("`", gene, "` ~ Smoking + svburden"))
  } else {
    formula <- as.formula(paste0("`", gene, "` ~ Smoking + log10(mut_load + 1)"))
  }
  
  
  # Fit the model
  model <- glm(formula, data = df, family = binomial)
  
  # Get tidy results with OR and CI
  results <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::select(term, estimate, conf.low, conf.high, p.value)
  colnames(results) <- c("Variable", "OR", "CI_2.5%", "CI_97.5%", "pvalue")
  
  # add significant star
  star_p <- c(0.05,0.01,0.001)
  results$star <- apply(results,1,function(x) paste0(rep("*",sum(as.numeric(x["pvalue"]) < star_p,na.rm = T)),collapse = ""))
  
  # add number
  results$wt.smoking = df %>% filter(Smoking=="Smoker") %>% filter(.data[[gene]] == 0) %>% nrow()
  results$mut.smoking = df %>% filter(Smoking=="Smoker") %>% filter(.data[[gene]] == 1) %>% nrow()
  results$mutpct.smoking = results$mut.smoking / (results$wt.smoking + results$mut.smoking)
  results$wt.nonsmoking = df %>% filter(Smoking=="Non-Smoker") %>% filter(.data[[gene]] == 0) %>% nrow()
  results$mut.nonsmoking = df %>% filter(Smoking=="Non-Smoker") %>% filter(.data[[gene]] == 1) %>% nrow()
  results$mutpct.nonsmoking = results$mut.nonsmoking / (results$wt.nonsmoking + results$mut.nonsmoking)
  
  # add gene
  results$gene = gene
  return(results %>% as.data.frame())
}

boxplot.function <- function(mydf=df.8grp.consistent,mycompasison=comparisons.8grp,myn="n_driver",
                             myy="# of driver genes per tumor altered by SNVs/indels",mytitle){
  mydf[,"n"] <- mydf[,myn]
  max_y <- max(mydf$n, na.rm = TRUE)
  test_results <- mydf %>%
    t_test(n ~ group, comparisons = mycompasison) 
  # smoking_test
  smoking_result <- mydf %>%
    t_test(n~Smoking, comparisons = c("Smoker","Non-Smoker"))
  
  # Step 1: Make sure both have a column named 'p'
  test_results$p_raw <- test_results$p
  smoking_result$p_raw <- smoking_result$p
  # Step 2: Combine the two result tables
  combined_results <- bind_rows(
    test_results
    # smoking_result
    )
  # Step 3: Apply FDR correction
  combined_results$p.adj <- p.adjust(combined_results$p_raw, method = "fdr")
  
  combined_results <- combined_results %>%
    mutate(p.adj = p.adjust(p, method = "fdr")) %>%
    mutate(p.adj.label = ifelse(p.adj<0.01,
                                formatC(p.adj, format = "e", digits = 1),
                                signif(p.adj, 2)
                                ) ) %>%
    # add_significance(p.col = "p.adj") %>%
    arrange(group1, group2) %>%
    mutate(y.position = max_y + 1 + row_number() * 0.5)
  
  print(combined_results)
  
  grp.num <- table(mydf$group) %>% as.data.frame() %>% setnames(.,old="Var1",new="group")
  p <- ggplot() +
    geom_violin(data=mydf, aes(x = group, y = n,fill = Smoking),trim = TRUE, color = "black") +
    geom_boxplot(data=mydf, aes(x = group, y = n),width = 0.1, outlier.shape = NA) +
    geom_text(data=grp.num,aes(x = group, y = 0,label = Freq),size = 2,vjust=2) +
    #geom_jitter(width = 0.1, size = 0.01, alpha = 0.2) +
    stat_pvalue_manual(
      combined_results,
      label = "p.adj.label",
      size = 2,
      tip.length = 0.01
    ) +
    scale_y_continuous(breaks = c(seq(0,100,1)))+
    scale_fill_manual(values = c("#006400", "#FF6EC7"),
                      labels=c("Never smoker","Smoker"),
                      breaks = c("Non-Smoker","Smoker")) +
    labs(x = "Group", y = myy,title = mytitle) +
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
  return(p)
}
# fisher_function <- function(df,genei){
#   df$status <- ifelse(is.na(df[,genei]),0,1)
#   df$smoking_status <- ifelse(df$Smoking=="Non-Smoker",0,1)
#   status0_smoking0 <- df %>% filter(status==0 & smoking_status==0) %>% nrow()
#   status0_smoking1 <- df %>% filter(status==0 & smoking_status==1) %>% nrow()
#   status1_smoking0 <- df %>% filter(status==1 & smoking_status==0) %>% nrow()
#   status1_smoking1 <- df %>% filter(status==1 & smoking_status==1) %>% nrow()
#   myfisher <- fisher.test(matrix(c(status0_smoking0,status0_smoking1,status1_smoking0,status1_smoking1),nrow = 2))
#   p <- myfisher$p.value
#   Estimate <- myfisher$estimate
#   n0=sum(df$status==0&df$smoking_status==0)
#   n1=sum(df$status==1&df$smoking_status==0)
#   s0=sum(df$status==0&df$smoking_status==1)
#   s1=sum(df$status==1&df$smoking_status==1)
#   return(c(p,Estimate,n0,n1,s0,s1))
# }
# 
# p.adjust_function <- function(df,col=3,p.bfrn,p.fdr){
#   df[,p.bfrn] <- p.adjust(df[,col],method = "bonferroni",n=nrow(df))
#   df[,p.fdr] <- p.adjust(df[,col],method = "BH")
#   return(df)
# }
# 
# point_plot_function <- function(result){
#   ggplot(result %>% 
#            mutate(sigfdr=ifelse(fdr<0.05,ifelse(est>1,"red","blue"),"grey"))%>%
#            mutate(label=ifelse(fdr<0.05,gene,"")))+
#     geom_point(aes(x=nonsmokerpct,y=smokerpct,size=-log10(fdr),color=sigfdr))+
#     geom_text_repel(aes(x=nonsmokerpct,y=smokerpct,label = label,color=sigfdr),
#                     size=1,          # Label font size
#                     segment.size = 0.2, # Line size (thickness of connecting lines)
#                     max.overlaps = 20  )+
#     scale_color_manual(values=c("grey","red","blue"),
#                        breaks=c("grey","red","blue"),
#                        labels=c("NA","smokers","never smokers"))+
#     scale_size_binned(
#       limits = c(0,60),
#       n.breaks = 4,  # Set the number of size bins to 3
#       range = c(0.1, 1),  # Define the range of point sizes
#       name = "-log10(FDR) from Fisher's test"  # Add a legend title for size
#     ) +
#     scale_x_continuous(limits = c(0,0.65))+
#     scale_y_continuous(limits = c(0,0.65))+
#     theme_minimal() +
#     labs(
#       x = "Percentage of mut-carrying nonsmokers",
#       y = "Percentage of mut-carrying smokers"
#     )+
#     theme(
#       text = element_text(size = 6, color = "black"),         # All text size 6 and black
#       axis.text.x = element_text(size = 6, color = "black"),    # Axis labels
#       axis.text.y = element_text(size = 6, color = "black"),    # Axis labels
#       axis.title = element_text(size = 6, color = "black"),   # Axis titles
#       plot.title = element_text(size = 6, color = "black", face = "bold"),  # Title
#       legend.text = element_text(size = 6, color = "black"),  # Legend text
#       legend.title = element_text(size = 6, color = "black"),  # Legend title
#       panel.background = element_blank(),
#       panel.spacing.x = unit(1, "mm"),
#       panel.grid = element_blank(),
#       panel.border = element_blank(),
#       legend.key.height = unit(3, "mm"),
#       legend.key.width  = unit(2, "mm"),
#       axis.ticks.length.x =unit(0.8, "mm"),
#       axis.ticks =  element_line(linewidth = 0.234,colour = "black"),
#       axis.line = element_line(linewidth = 0.234,colour = "black")
#     )
# }
# 
# barplot_plot_function <- function(result,mytitle,p,myx){
#   ggplot(data=result)+
#     geom_histogram(aes(x=total_mutated),binwidth =0.5,width = 0.5)+
#     facet_grid(Smoking~.,  # Facet by smoking across columns
#                labeller = as_labeller(c(`Non-Smoker` = "Never smoker", Smoker = "Smoker")))+
#     theme_minimal() +
#     labs(
#       x = myx,
#       y = "Number of tumors",
#       title = mytitle,
#     )+
#     geom_text(aes(x=10,
#                   y=200,
#                   label=paste0("P = ",signif(p,2),"\n","Linear regression adjusted for mutation load")),
#               size=2)+
#     theme(
#       text = element_text(size = 6, color = "black"),         # All text size 6 and black
#       axis.text.x = element_text(size = 6, color = "black"),    # Axis labels
#       axis.text.y = element_text(size = 6, color = "black"),    # Axis labels
#       axis.title = element_text(size = 6, color = "black"),   # Axis titles
#       plot.title = element_text(size = 6, color = "black"),  # Title
#       legend.text = element_text(size = 6, color = "black"),  # Legend text
#       legend.title = element_text(size = 6, color = "black"),  # Legend title
#       panel.background = element_blank(),
#       panel.spacing.x = unit(1, "mm"),
#       panel.grid = element_blank(),
#       panel.border = element_blank(),
#       legend.key.height = unit(3, "mm"),
#       legend.key.width  = unit(2, "mm"),
#       axis.ticks.length.x =unit(0.8, "mm"),
#       axis.ticks =  element_line(linewidth = 0.234,colour = "black"),
#       axis.line = element_line(linewidth = 0.234,colour = "black")
#     )
# }
# 
# percentage_plot_function <- function(result, mytitle, p, myx) {
#   # Compute proportions
#   result <- result %>%
#     group_by(Smoking, total_mutated) %>%
#     summarise(count = n(), .groups = "drop") %>%
#     mutate(percentage = count / sum(count) * 100)  # Convert count to percentage
#   
#   # Adjust sign for Smokers (negative for downward bars)
#   result <- result %>%
#     mutate(percentage = ifelse(Smoking == "Non-Smoker", percentage, -percentage))
#   
#   ggplot(data = result) +
#     geom_bar(aes(x = total_mutated, y = percentage, fill = Smoking), 
#              stat = "identity", 
#              width = 0.5) +
#     scale_fill_manual(values = c("Non-Smoker" = "#006400", "Smoker" = "#FF6EC7")) +  # Custom colors
#     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Reference line at y=0
#     theme_minimal() +
#     labs(
#       x = myx,
#       y = "Percentage of tumors",
#       title = mytitle,
#       fill = "Smoking Status"
#     ) +
#     geom_text(aes(x = max(result$total_mutated, na.rm = TRUE) * 0.9, 
#                   y = max(result$percentage, na.rm = TRUE) * 0.8,  # Adjusted for visibility
#                   label = paste0("P = ", signif(p, 2), "\nLinear regression adjusted for mutation load")),
#               size = 2) +
#     theme(
#       text = element_text(size = 6, color = "black"),         
#       axis.text.x = element_text(size = 6, color = "black"),  
#       axis.text.y = element_text(size = 6, color = "black"),  
#       axis.title = element_text(size = 6, color = "black"),   
#       plot.title = element_text(size = 6, color = "black"),   
#       legend.text = element_text(size = 6, color = "black"),  
#       legend.title = element_text(size = 6, color = "black"), 
#       panel.background = element_blank(),
#       panel.grid = element_blank(),
#       axis.ticks = element_line(linewidth = 0.234, colour = "black"),
#       axis.line = element_line(linewidth = 0.234, colour = "black")
#     )
# }



### path
sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
mysys <- "/Users/yangyang/Library/CloudStorage/"
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.suppl.txt")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.main.txt")
outpath <- file.path(scratch.path,"22-comprehensive-oncoplot","plot-test")
sbs4.path <- file.path(sherlock.path,"Data/SBS4/SBS4_annotation.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")

if (!dir.exists(outpath)){
  dir.create(outpath,showWarnings = F,recursive = T)
}

### necessary variable
driver_set <- c("EGFR","TP53","RBM10","KRAS",
                "PIK3CA","ERBB2","SETD2","CTNNB1","ARID1A",
                "HUWE1","MET","SMAD4","GNAS","RB1",
                "PTEN","CDKN2A","ARID2","KEAP1","NKX2-1",
                "NF1","SMARCA4","ATM","BRAF","ATRX",
                "STK11","UBA1","ARID1B","PDGFRA","PIK3R1",
                "ALK","TERT","CHD2","NFE2L2","RET",
                "MEN1","AKT1","EIF1AX","MDM2","ROS1",
                "MAP2K1","CCNE1","CDK4","NRAS","RIT1",
                "CUL3","U2AF1","HRAS","CCND1")
driver_set <- c("EGFR","TP53","RBM10","KRAS",
                "PIK3CA","SETD2","ARID1A","HUWE1","GNAS",
                "RB1","PTEN","CDKN2A","ARID2","KEAP1",
                "NF1","SMARCA4","ATM","STK11","ARID1B",
                "PDGFRA","ALK")

cnv_set <- c("CDKN2A","TERT","MDM2","EGFR","PTPRD",
             "MYC","CTNND2","NKX2-1","PTK6","MCL1",
             "FGFR1",
             # "GRIN2A",
             "CCNE1","STK11","MET","BRAF")
loss_set <- c("CDKN2A","PTPRD","STK11")
gain_set <- cnv_set[cnv_set %in% loss_set==FALSE]

fusion_set <- c("ALK","RET","ROS1","MET","NRG1",
                "FGFR1","FGFR2","FGFR3","FGFR4",
                "NTRK1","NTRK2","NTRK3")

cnv_values <- paste0("cnv.",cnv_set)
gistic_values <- paste0("gistic.",cnv_set)
loss_values <- paste0(loss_set,".loss")
gain_values <- paste0(gain_set,".gain")
loss_gistic <- paste0("gistic.",loss_set)
gain_gistic <- paste0("gistic.",gain_set)
### read
sbs4 <- read.delim2(sbs4.path)
oncoprint_table <- read.delim(oncoprint_table_path,check.names = F) %>%
  mutate(across(all_of(cnv_values), ~ ifelse((.x >= 5 & baseline_cn <= 2.7)|(.x >= 9 & baseline_cn > 2.7), "gain",
                                             ifelse((.x == 0 & baseline_cn <= 2.7)|(.x < (baseline_cn-2.7) & baseline_cn > 2.7),
                                                    "loss",
                                                    NA)))) %>%
  mutate(across(all_of(gain_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==3,"gain",NA)))) %>%
  mutate(across(all_of(loss_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==-3,"loss",NA))))

oncoprint_table <- oncoprint_table %>% setnames(.,old=c(loss_gistic,gain_gistic),new=c(loss_values,gain_values))
oncoprint_table <- merge(oncoprint_table,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
oncoprint_table <- oncoprint_table %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
histology <- read.csv(histology.path)
oncoprint_table <- merge(oncoprint_table,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

# #################### logistic test for driver genes/fusion/cnv genes ########################
# oncoprint_regression <- oncoprint_table
# allresult <- data_frame()
# for (mygene in c(driver_set,"fusion",loss_values,gain_values)){
#   #result <- logistic_gene_test(mygene, oncoprint_regression %>% filter(Smoking_SBS4 %in% c("Non-Smoker-SBS4","Smoker-noSBS4")))
#   result <- logistic_gene_test(mygene, oncoprint_regression %>% filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")))
#   result <- result %>% filter(Variable=="SmokingSmoker")
#   allresult <- rbind(allresult,result)
# }
# 
# # Categorize genes
# allresult <- allresult %>%
#   mutate(category = case_when(
#     gene %in% driver_set ~ "driver",
#     gene == "fusion" ~ "fusion",
#     gene %in% c(loss_values,gain_values) ~ "gistic",
#     TRUE ~ "other"
#   )) %>%
#   mutate(gistic_type = case_when(
#     category == "gistic" & gene %in% gain_values ~ "gain",
#     category == "gistic" & gene %in% loss_values ~ "loss",
#     TRUE ~ NA_character_
#   ))
# 
# # Create gene order
# ordered_genes <- allresult %>%
#   mutate(cat_rank = case_when(
#     category == "driver" ~ 1,
#     category == "fusion" ~ 2,
#     category == "gistic" ~ 3,
#     TRUE ~ 4
#   )) %>%
#   arrange(cat_rank, desc(mutpct.nonsmoking)) %>%
#   pull(gene)
# 
# allresult$gene <- factor(allresult$gene, levels = ordered_genes)
# 
# # Add star label positions
# allresult <- allresult %>%
#   mutate(
#     star_y = ifelse(OR < 1, mutpct.nonsmoking + 0.03, -mutpct.smoking - 0.03),
#     star_hjust = ifelse(OR < 1, 0, 1)
#   )
# 
# # Reshape data for plotting
# plot_df <- allresult %>%
#   dplyr::select(gene, category, gistic_type, mutpct.smoking, mutpct.nonsmoking, OR, pvalue, star) %>%
#   pivot_longer(
#     cols = c(mutpct.nonsmoking, mutpct.smoking),
#     names_to = "group",
#     values_to = "percent"
#   ) %>%
#   mutate(
#     direction = ifelse(group == "mutpct.smoking", "down", "up"),
#     percent = ifelse(direction == "down", -percent, percent)
#   )
# 
# 
# # Plot
# ggplot(plot_df, aes(x = gene, y = percent, fill = direction)) +
#   geom_col(width = 0.6, color = "black",linewidth=0) +
#   geom_text(data = allresult,
#             aes(x = gene, y = star_y, label = star, hjust = star_hjust),
#             inherit.aes = FALSE, size = 2, angle = 90, vjust = 0.75) +
#   scale_fill_manual(values = c("#006400", "#FF6EC7"),
#                     labels=c("Never smoker","Smoker"),
#                     breaks = c("up","down")) +
#   labs(x = "Gene", y = "Mutation %", title = "Mutation % in Smoking vs Non-Smoking") +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(size=6, colour = "black", angle = 90),
#     axis.text.y = element_text(size=6, colour = "black"),
#     axis.title.y = element_text(size=6, colour = "black"),
#     axis.title.x = element_blank(),
#     axis.line.y = element_line(size=0.1, colour = "black"),
#     axis.line.x.top = element_line(size=0.1, colour = "black"),
#     axis.line.x.bottom = element_line(size=0.1, colour = "black"),
#     axis.ticks.x = element_line(linewidth = 0.234),
#     axis.ticks.length.y = unit(0.8, "mm"),
#     axis.ticks.y = element_line(linewidth = 0.234),
#     plot.title = element_text(size=6, colour = "black"),
#     panel.background = element_blank(),
#     panel.spacing.x = unit(1, "mm"),
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     strip.text.x = element_text(size=6, angle = 0, colour = "black", hjust = 1, vjust = 0.5),
#     strip.text.y.left = element_text(size=6, angle = 0, colour = "black", hjust = 1, vjust = 0),
#     strip.background = element_blank(),
#     strip.placement = "top",
#     legend.key.height = unit(2, "mm"),
#     legend.key.width = unit(2, "mm"),
#     legend.title = element_blank(),
#     legend.text = element_text(size=6),
#     legend.position = "bottom",
#     plot.margin = grid::unit(c(1,2,0,2), "mm")
#   )
# ggsave(filename = file.path(outpath,"barplot.logistic.regression.correct.mutload.main.pdf"),width = 4,height = 2)
# ggsave(filename = file.path(outpath,"barplot.logistic.regression.correct.mutload.suppl.pdf"),width = 6,height = 2)


#################### Barplot Fisher's test for driver genes/fusion/cnv genes ########################
fisher_results <- map_dfr(c(driver_set,gain_values,loss_values), function(gene) {
  df <- oncoprint_table %>%
    filter(Smoking_SBS4 %in% c("Non-Smoker-noSBS4","Smoker-SBS4")) %>%
    dplyr::select(Smoking, gene = all_of(gene))

  if (gene %in% driver_set || gene %in% "fusion") {
    df <- df %>%
      mutate(mut = ifelse(is.na(gene), 0, 1))
  } else if (gene %in% gain_values) {
    df <- df %>%
      filter(gene %in% "NA"==FALSE) %>%
      mutate(mut = ifelse(gene %in% "gain", 1, 0))
  } else if (gene %in% loss_values) {
    df <- df %>%
      filter(gene %in% "NA"==FALSE) %>%
      mutate(mut = ifelse(gene %in% "loss", 1, 0))
  }

  # Create contingency table
  tab <- table(df$Smoking, df$mut)
  if (all(dim(tab) == c(2, 2))) {
    fisher_test <- fisher.test(tab)
    tibble(
      gene = gene,
      smoker_mut = tab["Smoker", "1"],
      smoker_wt = tab["Smoker", "0"],
      nonsmoker_mut = tab["Non-Smoker", "1"],
      nonsmoker_wt = tab["Non-Smoker", "0"],
      smoker_pct = tab["Smoker", "1"] / sum(tab["Smoker", ]),
      nonsmoker_pct = tab["Non-Smoker", "1"] / sum(tab["Non-Smoker", ]),
      OR = fisher_test$estimate,
      p_value = fisher_test$p.value
    )
  } else {
    tibble(
      gene = gene,
      smoker_mut = NA, smoker_wt = NA,
      nonsmoker_mut = NA, nonsmoker_wt = NA,
      smoker_pct = NA, nonsmoker_pct = NA,
      OR = NA, p_value = NA
    )
  }
})

# Adjust p-values
fisher_results <- fisher_results %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  mutate(star = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01 ~ "**",
    FDR < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Prepare data for plot
plot_df <- fisher_results %>%
  pivot_longer(cols = c(nonsmoker_pct, smoker_pct), names_to = "group", values_to = "percent") %>%
  mutate(direction = ifelse(group == "smoker_pct", "down", "up"),
         percent = ifelse(direction == "down", -percent, percent))

# Add star positions
fisher_results <- fisher_results %>%
  mutate(star_y = ifelse(OR < 1, nonsmoker_pct + 0.03, -smoker_pct - 0.03),
         star_hjust = ifelse(OR < 1, 0, 1))

# Add category info
fisher_results <- fisher_results %>%
  mutate(
    category = case_when(
      gene %in% driver_set ~ "driver",
      gene == "fusion" ~ "fusion",
      gene %in% gain_values ~ "gain",
      gene %in% loss_values ~ "loss",
      TRUE ~ "other"
    ),
    category_rank = case_when(
      category == "driver" ~ 1,
      category == "fusion" ~ 2,
      category == "loss" ~ 3,
      category == "gain" ~ 4,
      TRUE ~ 4
    )
  )

# Sort and apply as factor levels
ordered_genes <- fisher_results %>%
  arrange(category_rank, desc(nonsmoker_pct)) %>%
  pull(gene)

fisher_results$gene <- factor(fisher_results$gene, levels = ordered_genes)


ggplot(plot_df, aes(x = gene, y = percent, fill = direction)) +
  geom_col(width = 0.6, color = "black", linewidth = 0) +
  geom_text(data = fisher_results,
            aes(x = gene, y = star_y, label = star, hjust = star_hjust),
            inherit.aes = FALSE, size = 2, angle = 90, vjust = 0.75) +
  scale_fill_manual(values = c("#006400", "#FF6EC7"),
                    labels=c("Never smoker","Smoker"),
                    breaks = c("up","down")) +
  labs(x = "Gene", y = "Mutation %", title = "Mutation % in Smoking vs Non-Smoking") +
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
ggsave(filename = file.path(outpath,"barplot.fishertest.main.pdf"),width = 3.5,height = 2)
ggsave(filename = file.path(outpath,"barplot.fishertest.suppl.pdf"),width = 4,height = 2)

#################### Boxplot, T-test # of drivers per tumor ##################
driver_set <- c("EGFR","TP53","RBM10","KRAS",
                "PIK3CA","ERBB2","SETD2","CTNNB1","ARID1A",
                "HUWE1","MET","SMAD4","GNAS","RB1",
                "PTEN","CDKN2A","ARID2","KEAP1","NKX2-1",
                "NF1","SMARCA4","ATM","BRAF","ATRX",
                "STK11","UBA1","ARID1B","PDGFRA","PIK3R1",
                "ALK","TERT","CHD2","NFE2L2","RET",
                "MEN1","AKT1","EIF1AX","MDM2","ROS1",
                "MAP2K1","CCNE1","CDK4","NRAS","RIT1",
                "CUL3","U2AF1","HRAS","CCND1")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.suppl.txt")
oncoprint_table <- read.delim(oncoprint_table_path,check.names = F) %>%
  mutate(across(all_of(cnv_values), ~ ifelse((.x >= 5 & baseline_cn <= 2.7)|(.x >= 9 & baseline_cn > 2.7), "gain",
                                             ifelse((.x == 0 & baseline_cn <= 2.7)|(.x < (baseline_cn-2.7) & baseline_cn > 2.7),
                                                    "loss",
                                                    NA)))) %>%
  mutate(across(all_of(gain_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==3,"gain",NA)))) %>%
  mutate(across(all_of(loss_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==-3,"loss",NA))))
oncoprint_table <- merge(oncoprint_table,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

oncoprint_table <- oncoprint_table %>% setnames(.,old=c(loss_gistic,gain_gistic),new=c(loss_values,gain_values))
oncoprint_table <- merge(oncoprint_table,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
oncoprint_table <- oncoprint_table %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))
oncoprint_table_pub <- oncoprint_table %>%
  select(Tumor_Barcode,driver_set,loss_values,gain_values) 
oncoprint_table_pub[oncoprint_table_pub == "NA"] <- "Unavailable"
oncoprint_table_pub <- oncoprint_table_pub %>% setnames(.,old="Tumor_Barcode",new="Sample_ID")


write.csv(oncoprint_table_pub,
          file.path(scratch.path,"22-comprehensive-oncoplot","plot","snv_indel_focal_change_table.csv"),
          row.names = F,,quote = F)

# Step 1: Filter consistent smokers and non-smokers
df <- oncoprint_table %>%
  mutate(
    EGFR_status = ifelse(!is.na(EGFR), "EGFRmut", NA),
    KRAS_status = ifelse(!is.na(KRAS), "KRASmut", NA),
    TP53_status = ifelse(!is.na(TP53), "TP53mut", NA),
    Fusion_status = ifelse(!is.na(fusion), "Fusionmut", "NA")
  ) %>%
  rowwise() %>%
  mutate(
    n_driver = sum(c_across(all_of(driver_set)) %>% mutation_score()),
    n_cnv = sum(c_across(c("fusion", gain_values, loss_values)) %>% mutation_score()),
    n_driver_clean = sum(c_across(all_of(setdiff(driver_set, c("EGFR", "KRAS", "TP53")))) %>% mutation_score()),
    n_cnv_clean = sum(c_across(c(gain_values, loss_values)) %>% mutation_score())
  ) %>%
  ungroup() 


# df.2grp.consistent <- df %>%
#   filter(Smoking_SBS4 %in% c("Smoker-SBS4", "Non-Smoker-noSBS4")) %>%
#   mutate(
#     group = case_when(
#       Smoking_SBS4 == "Smoker-SBS4" ~ "Smoker",
#       Smoking_SBS4 == "Non-Smoker-noSBS4" ~ "Non-Smoker"
#     )
#   )
# 
# df.4grp <- df %>%
#   filter(Smoking_SBS4 %in% c("Smoker-SBS4","Non-Smoker-noSBS4","Smoker-noSBS4","Non-Smoker-SBS4")) %>%
#   mutate(
#     group = case_when(
#       Smoking_SBS4 == "Smoker-SBS4" ~ "Smoker-SBS4",
#       Smoking_SBS4 == "Non-Smoker-noSBS4" ~ "Non-Smoker-noSBS4",
#       Smoking_SBS4 == "Smoker-noSBS4" ~ "Smoker-noSBS4",
#       Smoking_SBS4 == "Non-Smoker-SBS4" ~ "Non-Smoker-SBS4"
#     )
#   )
# 
# df.6grp.consistent <-  df %>%
#   filter(Smoking_SBS4 %in% c("Smoker-SBS4", "Non-Smoker-noSBS4")) %>%
#   filter(!(EGFR_status %in% "EGFRmut" & KRAS_status %in% "KRASmut")) %>%
#   filter(Histology=="Adenocarcinoma") %>%
#   mutate(
#     mut_group = case_when(
#       EGFR_status == "EGFRmut" ~ "EGFRmut",
#       KRAS_status == "KRASmut" ~ "KRASmut",
#       TRUE ~ "Other"
#     ),
#     smoking_group = case_when(
#       Smoking_SBS4 == "Smoker-SBS4" ~ "Smoker",
#       Smoking_SBS4 == "Non-Smoker-noSBS4" ~ "Non-Smoker"
#     ),
#     TP53_group = case_when(
#       TP53_status == "TP53mut" ~ "TP53mut",
#       TRUE ~ "TP53wt"
#     ),
#     Fusion_group = case_when(
#       Fusion_status == "Fusionmut" ~ "Fusionmut",
#       TRUE ~ "Fusionwt"
#     ),
#     group = paste(smoking_group, mut_group, sep = "_")
#   )

df.8grp.consistent <-  df %>%
  filter(Smoking_SBS4 %in% c("Smoker-SBS4", "Non-Smoker-noSBS4")) %>%
  filter(Histology=="Adenocarcinoma") %>%
  #filter(!(EGFR_status %in% "EGFRmut" & KRAS_status %in% "KRASmut")) %>%
  mutate(
    mut_group = case_when(
      EGFR_status == "EGFRmut" ~ "EGFRmut",
      KRAS_status == "KRASmut" ~ "KRASmut",
      Fusion_status == "Fusionmut" ~ "Fusionmut",
      TRUE ~ "Other"
    ),
    smoking_group = case_when(
      Smoking_SBS4 == "Smoker-SBS4" ~ "Smoker",
      Smoking_SBS4 == "Non-Smoker-noSBS4" ~ "Non-Smoker"
    ),
    group = paste(smoking_group, mut_group, sep = "_")
  ) %>%
  mutate(group = factor(group,levels = c("Non-Smoker_EGFRmut","Non-Smoker_KRASmut","Non-Smoker_Fusionmut","Non-Smoker_Other",
                                            "Smoker_EGFRmut","Smoker_KRASmut","Smoker_Fusionmut","Smoker_Other")))

df.8grp.consistent %>% count(group)
df.8grp.consistent %>% nrow()
# Only valid combinations that share one factor
comparisons.8grp <- list(
  c("Smoker_KRASmut","Smoker_Other"),      #kras
  c("Smoker_KRASmut","Smoker_EGFRmut"),   #kras-egfr
  c("Smoker_KRASmut","Non-Smoker_KRASmut"), #smoker
  # c("Smoker_KRASmut","Smoker_Fusionmut"),   #kras-fusion
  c("Smoker_Other","Non-Smoker_Other"),  #smoker
  c("Smoker_Other","Smoker_EGFRmut"),    #egfr
  # c("Smoker_Other","Smoker_Fusionmut"),  #fusion
  c("Non-Smoker_Fusionmut","Non-Smoker_Other"), #fusion 
  c("Non-Smoker_Fusionmut","Non-Smoker_EGFRmut"), #egfr-fusion
  c("Non-Smoker_Fusionmut","Non-Smoker_KRASmut"), #kras-fuison
  # c("Non-Smoker_Fusionmut","Smoker_Fusionmut"),  #smoker
  c("Non-Smoker_Other","Non-Smoker_EGFRmut"), #egfr
  c("Non-Smoker_Other","Non-Smoker_KRASmut"), #kras
  c("Smoker_EGFRmut","Non-Smoker_EGFRmut"), #smoker
  # c("Smoker_EGFRmut","Smoker_Fusionmut"),   #egfr-fusion
  c("Non-Smoker_EGFRmut","Non-Smoker_KRASmut") #kras-egfr
)
# comparisons.6grp <- list(
#   c("Smoker_EGFRmut", "Non-Smoker_EGFRmut"),
#   c("Smoker_KRASmut", "Non-Smoker_KRASmut"),
#   c("Smoker_Other", "Non-Smoker_Other"),
#   c("Smoker_EGFRmut", "Smoker_KRASmut"),
#   c("Non-Smoker_EGFRmut", "Non-Smoker_KRASmut"),
#   c("Smoker_EGFRmut", "Smoker_Other"),
#   c("Non-Smoker_EGFRmut", "Non-Smoker_Other"),
#   c("Smoker_KRASmut", "Smoker_Other"),
#   c("Non-Smoker_KRASmut", "Non-Smoker_Other")
# )
# comparisons.2grp <- list(
#   c("Smoker", "Non-Smoker")
# )
# comparisons.4grp <- list(
#   c("Non-Smoker-noSBS4","Non-Smoker-SBS4"),
#   c("Smoker-noSBS4","Smoker-SBS4"),
#   c("Non-Smoker-noSBS4","Smoker-noSBS4"),
#   c("Non-Smoker-SBS4","Smoker-SBS4")
# )
# all
# boxplot.function(mydf=df.4grp,mycompasison=comparisons.4grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "number.of.snvdriver.4grp")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.4grp.pdf"),width = 2,height = 4) #full genelist
# 
# boxplot.function(mydf=df.4grp,mycompasison=comparisons.4grp,myn="n_cnv",myy="# of driver genes per tumor altered by SVs",mytitle = "number.of.svdriver.4grp")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.svdriver.4grp.pdf"),width = 2,height = 4) #full genelist
# 
# boxplot.function(mydf=df.2grp.consistent,mycompasison=comparisons.2grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "snvdriver.2grp.consistentsmk")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.2grp.consistentsmk.pdf"),width = 2,height = 4) #full genelist
# 
# boxplot.function(mydf=df.2grp.consistent,mycompasison=comparisons.2grp,myn="n_cnv",myy="# of driver genes per tumor altered by SVs",mytitle = "svdriver.2grp.consistentsmk")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.svdriver.2grp.consistentsmk.pdf"),width = 2,height = 4)
# 
# boxplot.function(mydf=df.6grp.consistent,mycompasison=comparisons.6grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "snvdriver.6grp.consistentsmk.luad")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.6grp.consistentsmk.luad.pdf"),width = 3,height = 5) #full genelist
# 
# boxplot.function(mydf=df.6grp.consistent,mycompasison=comparisons.6grp,myn="n_cnv",myy="# of driver genes per tumor altered by SVs",mytitle = "svdriver.6grp.consistentsmk.luad")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.svdriver.6grp.consistentsmk.luad.pdf"),width = 3,height = 5)

boxplot.function(mydf=df.8grp.consistent,mycompasison=comparisons.8grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "snvdriver.8grp.consistentsmk.luad")
ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.8grp.consistentsmk.luad.pdf"),width = 3.5,height = 7) #full genelist

boxplot.function(mydf=df.8grp.consistent,mycompasison=comparisons.8grp,myn="n_cnv_clean",myy="# of driver genes per tumor altered by SVs",mytitle = "svdriver.8grp.consistentsmk.luad")
ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.cleansvdriver.8grp.consistentsmk.luad.pdf"),width = 3.5,height = 7)

# boxplot.function(mydf=df.8grp.consistent,mycompasison=comparisons.8grp,myn="Age",myy="Age",mytitle = "Age.8grp.consistentsmk.luad")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.age.8grp.consistentsmk.luad.pdf"),width = 4,height = 7)
# high quality
# boxplot.function(mydf=df.4grp %>% filter(high_quality==1),mycompasison=comparisons.4grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "number.of.snvdriver.4grp.hq")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.4grp.hq.pdf"),width = 2,height = 4) #full genelist
# 
# boxplot.function(mydf=df.4grp %>% filter(high_quality==1),mycompasison=comparisons.4grp,myn="n_cnv",myy="# of driver genes per tumor altered by SVs",mytitle = "number.of.svdriver.4grp.hq")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.svdriver.4grp.hq.pdf"),width = 2,height = 4) #full genelist
# 
# boxplot.function(mydf=df.2grp.consistent %>% filter(high_quality==1),mycompasison=comparisons.2grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "snvdriver.2grp.consistentsmk.hq")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.2grp.consistentsmk.hq.pdf"),width = 2,height = 4) #full genelist
# 
# boxplot.function(mydf=df.2grp.consistent %>% filter(high_quality==1),mycompasison=comparisons.2grp,myn="n_cnv",myy="# of driver genes per tumor altered by SVs",mytitle = "svdriver.2grp.consistentsmk.hq")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.svdriver.2grp.consistentsmk.hq.pdf"),width = 2,height = 4)
# 
# boxplot.function(mydf=df.6grp.consistent %>% filter(high_quality==1),mycompasison=comparisons.6grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "snvdriver.6grp.consistentsmk.luad.hq")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.6grp.consistentsmk.luad.hq.pdf"),width = 3,height = 5) #full genelist
# 
# boxplot.function(mydf=df.6grp.consistent %>% filter(high_quality==1),mycompasison=comparisons.6grp,myn="n_cnv",myy="# of driver genes per tumor altered by SVs",mytitle = "svdriver.6grp.consistentsmk.luad.hq")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.svdriver.6grp.consistentsmk.luad.hq.pdf"),width = 3,height = 5)

boxplot.function(mydf=df.8grp.consistent %>% filter(high_quality==1),mycompasison=comparisons.8grp,myn="n_driver",myy="# of driver genes per tumor altered by SNVs/indels",mytitle = "snvdriver.8grp.consistentsmk.luad.hq")
ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.snvdriver.8grp.consistentsmk.luad.hq.pdf"),width = 3.5,height = 7) #full genelist

boxplot.function(mydf=df.8grp.consistent %>% filter(high_quality==1),mycompasison=comparisons.8grp,myn="n_cnv_clean",myy="# of driver genes per tumor altered by SVs",mytitle = "svdriver.8grp.consistentsmk.luad.hq")
ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.cleansvdriver.8grp.consistentsmk.luad.hq.pdf"),width = 3.5,height = 7)

# boxplot.function(mydf=df.8grp.consistent %>% filter(high_quality==1),mycompasison=comparisons.8grp,myn="Age",myy="Age",mytitle = "Age.8grp.consistentsmk.luad.hq")
# ggsave(filename = file.path(outpath,"boxplot.ttest.number.of.age.8grp.consistentsmk.luad.hq.pdf"),width = 4,height = 7)
######################## linear regression for # of drivers per sample #################
driver_set <- c("EGFR","TP53","RBM10","KRAS",
                "PIK3CA","ERBB2","SETD2","CTNNB1","ARID1A",
                "HUWE1","MET","SMAD4","GNAS","RB1",
                "PTEN","CDKN2A","ARID2","KEAP1","NKX2-1",
                "NF1","SMARCA4","ATM","BRAF","ATRX",
                "STK11","UBA1","ARID1B","PDGFRA","PIK3R1",
                "ALK","TERT","CHD2","NFE2L2","RET",
                "MEN1","AKT1","EIF1AX","MDM2","ROS1",
                "MAP2K1","CCNE1","CDK4","NRAS","RIT1",
                "CUL3","U2AF1","HRAS","CCND1")
oncoprint_table_path <- file.path(scratch.path,"22-comprehensive-oncoplot","plot","plot.table.suppl.txt")
oncoprint_table <- read.delim(oncoprint_table_path,check.names = F) %>%
  mutate(across(all_of(cnv_values), ~ ifelse((.x >= 5 & baseline_cn <= 2.7)|(.x >= 9 & baseline_cn > 2.7), "gain",
                                             ifelse((.x == 0 & baseline_cn <= 2.7)|(.x < (baseline_cn-2.7) & baseline_cn > 2.7),
                                                    "loss",
                                                    NA)))) %>%
  mutate(across(all_of(gain_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==3,"gain",NA)))) %>%
  mutate(across(all_of(loss_gistic), ~ ifelse(is.na(.x),"NA",ifelse(.x ==-3,"loss",NA))))
oncoprint_table <- merge(oncoprint_table,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

oncoprint_table <- oncoprint_table %>% setnames(.,old=c(loss_gistic,gain_gistic),new=c(loss_values,gain_values))
oncoprint_table <- merge(oncoprint_table,sbs4,by.x = "Tumor_Barcode",by.y = "Tumor_Barcode",all.x = T)
oncoprint_table <- oncoprint_table %>% mutate(Smoking_SBS4=paste0(Smoking,"-",ifelse(SBS4>0,"SBS4","noSBS4")))

df <- oncoprint_table %>%
  filter(Smoking_SBS4 %in% c("Smoker-SBS4", "Non-Smoker-noSBS4")) %>%
  filter(Assigned_Population %in% c("EAS","EUR")) %>%
  mutate(
    EGFR_status = ifelse(!is.na(EGFR), "Mut", "Wild"),
    KRAS_status = ifelse(!is.na(KRAS), "Mut", "Wild"),
    TP53_status = ifelse(!is.na(TP53), "Mut", "Wild"),
    Fusion_status = ifelse(!is.na(fusion), "Mut", "Wild"),
    Stage=ifelse(Stage %in% c("I","IA","IA1","IA2","IA3","IB"),"I",
                 ifelse(Stage %in% c("II","IIA","IIB"),"II",
                        ifelse(Stage %in% c("III","IIIA","IIIB"),"III",
                               ifelse(Stage %in% c("IV","IVA"),"IV",NA))))
  ) %>%
  mutate(Assigned_Population=factor(Assigned_Population,levels = c("EAS","EUR")),
         Gender=factor(Gender,levels = c("Female","Male")),
         Smoking=factor(Smoking,levels = c("Non-Smoker","Smoker","Unknown")),
         wgd=factor(wgd,levels = c("nWGD","WGD")),
         scna_group=factor(scna_group,levels = c("Piano","Mezzo-forte","Forte")),
         across(all_of(c("EGFR_status","KRAS_status","TP53_status")), ~ factor(.,levels = c("Wild","Mut"))),
         Fusion_status=factor(Fusion_status,levels = c("Wild","Mut")),
         Histology=factor(Histology,levels = c("Adenocarcinoma",
                                               "Squamous cell carcinoma",
                                               "Carcinoid tumor",
                                               "Adenosquamous carcinoma",
                                               "Others")),
         Stage=factor(Stage,levels=c("I","II","III","IV"))
  ) %>%
  rowwise() %>%
  mutate(
    n_driver = sum(c_across(all_of(driver_set)) %>% mutation_score()),
    n_cnv = sum(c_across(c("fusion", gain_values, loss_values)) %>% mutation_score()),
    n_driver_clean = sum(c_across(all_of(setdiff(driver_set, c("EGFR", "KRAS", "TP53")))) %>% mutation_score()),
    n_cnv_clean = sum(c_across(c(gain_values, loss_values)) %>% mutation_score())
  ) %>%
  ungroup()

# SNV linear regression in all
# myformula.snv.gene.load <- as.formula(paste0("n_driver","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status+log10(mut_load+1)+log10(indel_load+1)"))
# model.snv.gene.load <- glm(myformula.snv.gene.load, data = df, family = gaussian())
# forest.plot(model=model.snv.gene.load,mytitle = "snvdriver.gene.load")
# ggsave(filename = file.path(outpath,"linear.regression.number.of.snvdriver.gene.load.pdf"),width = 2.5,height = 3)

myformula.snv.gene <- as.formula(paste0("n_driver","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status"))
model.snv.gene <- glm(myformula.snv.gene, data = df, family = gaussian())
forest.plot(model=model.snv.gene,mytitle = "snvdriver.gene")
ggsave(filename = file.path(outpath,"linear.regression.number.of.snvdriver.gene.pdf"),width = 2,height = 2)

# SV linear regression in all
# myformula.sv.gene.load <- as.formula(paste0("n_cnv_clean","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status+svburden"))
# model.sv.gene.load <- glm(myformula.sv.gene.load, data = df, family = gaussian())
# forest.plot(model=model.sv.gene.load,mytitle = "clean.svdriver.gene.load")
# ggsave(filename = file.path(outpath,"linear.regression.number.of.cleansvdriver.gene.load.pdf"),width = 2.5,height = 3)

myformula.sv.gene <- as.formula(paste0("n_cnv_clean","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status"))
model.sv.gene <- glm(myformula.sv.gene, data = df, family = gaussian())
forest.plot(model=model.sv.gene,mytitle = "clean.svdriver.gene")
ggsave(filename = file.path(outpath,"linear.regression.number.of.cleansvdriver.gene.pdf"),width = 2,height = 2)

# SNV linear regression high quality
# myformula.snv.gene.load <- as.formula(paste0("n_driver","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status+log10(mut_load+1)+log10(indel_load+1)"))
# model.snv.gene.load <- glm(myformula.snv.gene.load, data = df %>% filter(high_quality==1), family = gaussian())
# forest.plot(model=model.snv.gene.load,mytitle = "snvdriver.gene.load.hq")
# ggsave(filename = file.path(outpath,"linear.regression.number.of.snvdriver.gene.load.hq.pdf"),width = 2.5,height = 3)

myformula.snv.gene <- as.formula(paste0("n_driver","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status"))
model.snv.gene <- glm(myformula.snv.gene, data = df %>% filter(high_quality==1), family = gaussian())
forest.plot(model=model.snv.gene,mytitle = "snvdriver.gene.hq")
ggsave(filename = file.path(outpath,"linear.regression.number.of.snvdriver.gene.hq.pdf"),width = 2,height = 2)

# SV linear regressionhigh quality
# myformula.sv.gene.load <- as.formula(paste0("n_cnv_clean","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status+svburden"))
# model.sv.gene.load <- glm(myformula.sv.gene.load, data = df %>% filter(high_quality==1), family = gaussian())
# forest.plot(model=model.sv.gene.load,mytitle = "cleansvdriver.gene.load.hq")
# ggsave(filename = file.path(outpath,"linear.regression.number.of.cleansvdriver.gene.load.hq.pdf"),width = 2.5,height = 3)

myformula.sv.gene <- as.formula(paste0("n_cnv_clean","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+TP53_status+EGFR_status+KRAS_status+Fusion_status"))
model.sv.gene <- glm(myformula.sv.gene, data = df %>% filter(high_quality==1), family = gaussian())
forest.plot(model=model.sv.gene,mytitle = "cleansvdriver.gene.hq")
ggsave(filename = file.path(outpath,"linear.regression.number.of.cleansvdriver.gene.hq.pdf"),width = 2,height = 2)

myformula.snv.gene <- as.formula(paste0("n_driver","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+TP53_status+EGFR_status+KRAS_status+Fusion_status"))
model.snv.gene <- glm(myformula.snv.gene, data = df %>% filter(high_quality==1), family = gaussian())
forest.plot(model=model.snv.gene,mytitle = "snvdriver.gene.hq")
ggsave(filename = file.path(outpath,"linear.regression.number.of.snvdriver.gene.addstage.hq.pdf"),width = 2,height = 2)

myformula.sv.gene <- as.formula(paste0("n_cnv_clean","~Assigned_Population+Smoking+Gender+Histology+purity+wgd+scna_group+Stage+TP53_status+EGFR_status+KRAS_status+Fusion_status"))
model.sv.gene <- glm(myformula.sv.gene, data = df %>% filter(high_quality==1), family = gaussian())
forest.plot(model=model.sv.gene,mytitle = "cleansvdriver.gene.hq")
ggsave(filename = file.path(outpath,"linear.regression.number.of.cleansvdriver.gene.addstage.hq.pdf"),width = 2,height = 2)
# #####################################################
# ### test1-1: count % of tumors carrying driver SNVs/indels, Fishers' test
# df <- oncoprint_table %>% filter(Smoking %in% c("Non-Smoker","Smoker"))
# result <- data.frame()
# for (genei in driver_set){
#   mylist <- fisher_function(df,genei)
#   p <- mylist[1]
#   est <- mylist[2]
#   n0 <- mylist[3]
#   n1 <- mylist[4]
#   s0 <- mylist[5]
#   s1 <- mylist[6]
#   result <- rbind(result,
#                             data.frame("gene"=genei,
#                                        "p"=p,
#                                        "est"=est,
#                                        "nonsmokerwild"=n0,
#                                        "nonsmokermut"=n1,
#                                        "smokerwild"=s0,
#                                        "smokermut"=s1,
#                                        "nonsmokerpct"=n1/(n0+n1),
#                                        "smokerpct"=s1/(s0+s1)))
# }
# result <- p.adjust_function(result,2,"bfrn","fdr")
# p.snv <- point_plot_function(result)
# 
# 
# ### test1-2: count % of tumors carrying CNVs, Fishers' test
# cnv_set_loss <- c("CDKN2A","STK11","PTPRD")
# cnv_set_gain <- cnv_set[cnv_set %in% cnv_set_loss==FALSE]
# cnv_set_loss <- paste0("cnv.",cnv_set_loss)
# cnv_set_gain <- paste0("cnv.",cnv_set_gain)
# df <- oncoprint_table %>% 
#   filter(Smoking %in% c("Non-Smoker","Smoker")) %>%
#   mutate(across(all_of(cnv_set_gain), ~ ifelse(. %in% c("gain"), "gain", NA))) %>% 
#   mutate(across(all_of(cnv_set_loss), ~ ifelse(. %in% c("loss"), "loss", NA)))
# result <- data.frame()
# for (genei in c("fusion",cnv_values)){
#   mylist <- fisher_function(df,genei)
#   p <- mylist[1]
#   est <- mylist[2]
#   n0 <- mylist[3]
#   n1 <- mylist[4]
#   s0 <- mylist[5]
#   s1 <- mylist[6]
#   result <- rbind(result,
#                   data.frame("gene"=genei,
#                              "p"=p,
#                              "est"=est,
#                              "nonsmokerwild"=n0,
#                              "nonsmokermut"=n1,
#                              "smokerwild"=s0,
#                              "smokermut"=s1,
#                              "nonsmokerpct"=n1/(n0+n1),
#                              "smokerpct"=s1/(s0+s1)))
# }
# result <- p.adjust_function(result,2,"bfrn","fdr")
# result$gene <- gsub("cnv.", "", result$gene)
# p.cnv <- point_plot_function(result)+
#   scale_x_continuous(limits = c(0,0.3))+
#   scale_y_continuous(limits = c(0,0.3))
# 
# ggarrange(p.snv,p.cnv,nrow=1)
# ggsave(filename = file.path(outpath,"fishers.test.pdf"),width = 8,height = 3)
# 
# ### test2-1: count number of driver SNVs/indels per tumor
# df <- oncoprint_table %>% filter(Smoking %in% c("Non-Smoker","Smoker"))
# # Calculate the total number of mutated genes for each sample
# result <- df %>%
#   rowwise() %>%  # Perform row-wise operations
#   mutate(total_mutated = sum(!is.na(c_across(all_of(driver_set))))) %>%  # Count non-NA entries for each sample
#   ungroup() %>%  # Remove row-wise grouping
#   select(Tumor_Barcode,Smoking, total_mutated,mut_load)  # Keep only relevant columns
# # Perform linear regression
# model <- lm(total_mutated ~ Smoking + mut_load, data = result)
# # Summarize the model to see the results
# p <- summary(model)$coefficients["SmokingSmoker","Pr(>|t|)"]
# # plot
# p.snv <- barplot_plot_function(result,mytitle="Driver mutations",p,myx="Number of mutated driver genes")+
#   scale_x_continuous(breaks = seq(0,14,1))
# pct.snv <- percentage_plot_function(result,mytitle="Driver mutations",p,myx="Number of mutated driver genes")+
#   scale_x_continuous(breaks = seq(0,14,1))+
#   scale_y_continuous(breaks = c(30,15,0,-15,-30),limits = c(-15,30))
# 
# ### test2-2: count number of CNV per tumor
# cnv_set_loss <- c("CDKN2A","STK11","PTPRD")
# cnv_set_gain <- cnv_set[cnv_set %in% cnv_set_loss==FALSE]
# cnv_set_loss <- paste0("cnv-",cnv_set_loss)
# cnv_set_gain <- paste0("cnv-",cnv_set_gain)
# df <- oncoprint_table %>% 
#   filter(Smoking %in% c("Non-Smoker","Smoker")) %>%
#   mutate(across(all_of(cnv_set_gain), ~ ifelse(. %in% c("gain"), "gain", NA))) %>% 
#   mutate(across(all_of(cnv_set_loss), ~ ifelse(. %in% c("loss"), "loss", NA)))
# # Calculate the total number of mutated genes for each sample
# result <- df %>%
#   rowwise() %>%  # Perform row-wise operations
#   mutate(total_mutated = sum(!is.na(c_across(all_of(c("fusion",cnv_values)))))) %>%  # Count non-NA entries for each sample
#   ungroup() %>%  # Remove row-wise grouping
#   select(Tumor_Barcode,Smoking, total_mutated,svburden)  # Keep only relevant columns
# # Perform linear regression
# model <- lm(total_mutated ~ Smoking + svburden, data = result)
# # Summarize the model to see the results
# p <- summary(model)$coefficients["SmokingSmoker","Pr(>|t|)"]
# # plot
# p.cnv <- barplot_plot_function(result,mytitle="Fusions/CNVs",p,myx="Number of CN change genes/fusions+")+
#   scale_x_continuous(breaks = seq(0,13,1))
# pct.cnv <- percentage_plot_function(result,mytitle="Fusions/CNVs",p,myx="Number of CN change genes/fusions+")+
#   scale_x_continuous(breaks = seq(0,13,1))+
#   scale_y_continuous(breaks = c(30,15,0,-15,-30),limits = c(-15,30))
# 
# p.out <- ggarrange(p.snv,p.cnv,nrow=1)
# ggsave(plot=p.out,filename = file.path(outpath,"linear.regression.pdf"),width = 6,height = 3)
# 
# pct.out <- ggarrange(pct.snv,pct.cnv,nrow=1)
# ggsave(plot=pct.out,filename = file.path(outpath,"linear.regression.pct.pdf"),width = 6,height = 2)
# 
# ############################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# test.df <- df %>% mutate(across(all_of(driver_set), ~ ifelse(is.na(.), 0, 1))) %>% 
#   gather(.,gene,status,driver_set) %>%
#   group_by(Smoking,gene,status) %>%
#   count() %>%
#   group_by(Smoking,gene) %>%
#   mutate(pct=n/sum(n)) %>%
#   filter(status==0) %>% 
#   mutate(pct.mut=1-pct) %>%
#   as.data.frame()
# mywilcox <- wilcox.test(test.df %>% filter(Smoking=="Smoker") %>% pull(pct.mut), 
#             test.df %>% filter(Smoking=="Non-Smoker") %>% pull(pct.mut), 
#             alternative = "greater", paired = TRUE)
# mywilcox
# 
# ### test1-2: count % of tumors carrying CNV/fusions and test the same way as test genes
# df <- oncoprint_table %>% filter(Smoking %in% c("Non-Smoker","Smoker"))
# cnv_set_loss <- c("CDKN2A","STK11","PTPRD")
# cnv_set_gain <- cnv_set[cnv_set %in% cnv_set_loss==FALSE]
# cnv_set_loss <- paste0("cnv-",cnv_set_loss)
# cnv_set_gain <- paste0("cnv-",cnv_set_gain)
# 
# test.df <- df %>% 
#   mutate(across(all_of(c("fusion")), ~ ifelse(is.na(.), 0, 1))) %>% 
#   mutate(across(all_of(cnv_set_gain), ~ ifelse(. %in% c("gain"), 1, 0))) %>% 
#   mutate(across(all_of(cnv_set_loss), ~ ifelse(. %in% c("loss"), 1, 0))) %>% 
#   gather(.,gene,status,c("fusion",cnv_values)) %>%
#   group_by(Smoking,gene,status) %>%
#   count() %>%
#   group_by(Smoking,gene) %>%
#   mutate(pct=n/sum(n)) %>%
#   filter(status==0) %>% 
#   mutate(pct.mut=1-pct) %>%
#   as.data.frame()
# mywilcox <- wilcox.test(test.df %>% filter(Smoking=="Smoker") %>% pull(pct), 
#                         test.df %>% filter(Smoking=="Non-Smoker") %>% pull(pct), 
#                         alternative = "greater", paired = TRUE)
# mywilcox
# t.test(test.df %>% filter(Smoking=="Smoker") %>% pull(pct), 
#             test.df %>% filter(Smoking=="Non-Smoker") %>% pull(pct), 
#             alternative = "greater", paired = TRUE)
# 
# ggplot(test.df %>% mutate(gene=factor(gene,levels=c("fusion",cnv_values))))+
#   geom_bar(aes(x=gene,y=pct.mut),stat = "identity")+facet_grid(.~Smoking)+
#   theme(axis.text.x = element_text(angle = 90))
# 
