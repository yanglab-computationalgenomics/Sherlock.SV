library(ConsensusClusterPlus)
library(clusterCrit)
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)

# args <- c("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/05-sigprofiler/1217_samples_hg38/Sherlock_SV49_simple/CH49/All_Solutions/SBS49_9_Signatures/Activities/SBS49_S9_NMF_Activities.txt",
#           "15",
#           "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/scratch/07-cluster/1217_samples_hg38/Sherlock_SV49_simple/9signature",
#           "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/SV_Signatures_Clusters/Sherlock_SV49_simple_9signature/evaluation",
#          "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/Sherlock/Data/SAMPLE_INFO/sherlock_sampleinfo_merged.txt")
# sample_sig_path <- args[1]
# cluster.num  <- as.numeric(args[2])
# output_path <-  args[3]
# setwd_path  <- args[4]
# info_path <- args[5]

sherlock.path <- "/Users/yangyang/Library/CloudStorage/OneDrive-TheUniversityofChicago/office/GitHub/Sherlock"
scratch.path <- file.path(sherlock.path,"scratch")
sample_sig_path <- file.path(scratch.path,"04-shatterseek","1217_samples_hg38","sherlock_simple_sv_with_signature.tsv")
cluster.num  <- 15
output_path <-  file.path(scratch.path,"07-cluster/1217_samples_hg38/Sherlock_SV49_simple/8signature_myannotation")
setwd_path  <- file.path(sherlock.path,"SV_Signatures_Clusters/Sherlock_SV49_simple_8signature_15cluste_myannotation")
info_path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_info.20250128.txt")
histology.path <- file.path(sherlock.path,"Data/SAMPLE_INFO/wgs_1217_to_1209_histology.csv")


### set wd
if (!dir.exists(setwd_path)){
  dir.create(setwd_path,showWarnings = F,recursive = T)
}
setwd(setwd_path)


### build dir
if (!dir.exists(output_path)){
  dir.create(output_path,recursive = T,showWarnings = F)
}

### read files
# sample_sig  <- read.delim(sample_sig_path,row.names = 1)
sample_info <- read.delim(info_path)
histology <- read.csv(histology.path)
sample_info <- merge(sample_info,histology,all.x = T) %>%
  # Remove 'Histology' and 'Histology_Detail' columns
  select(-Histology, -Histology_Detail) %>%
  # Rename 'New_tumor_type' to 'Histology'
  dplyr::rename(Histology = New_tumor_type) %>%
  # Filter rows where Inclusion == "In"
  filter(Inclusion == "In")

sample_sig  <- read.delim(sample_sig_path) %>% count(SAMPLE,Signature) %>% spread(.,Signature,n) %>% replace(is.na(.), 0)
row.names(sample_sig) <- sample_sig$SAMPLE
sample_sig <- sample_sig %>% select("Del1","Del2","Del3","TD1","TD2","Fb inv","Large intra","Tra")

sample_sig <- sample_sig[rownames(sample_sig) %in% sample_info$Tumor_Barcode,]
sample_sig_matrix <- as.matrix(t(sample_sig))
sample_sig_matrix_numeric <- matrix(as.numeric( as.matrix(sample_sig)),
                                    ncol = ncol( as.matrix(sample_sig)))
### run clusters
run_cluster_function <- function(input,mytitle,myclusterAlg,mydistance,mycluster.num){
  results_sig_cluster <-  ConsensusClusterPlus(input,maxK = mycluster.num,reps = 50,pItem=0.8,pFeature = 1,
                                               title = mytitle, clusterAlg = myclusterAlg,distance = mydistance, seed = 123456,plot = "png")
  
  for (i in  2:mycluster.num){
    # add  consensusClass
    sv_signature <- rbind(sample_sig_matrix,results_sig_cluster[[i]][["consensusClass"]])
    row.names(sv_signature)[nrow(sv_signature)] <- "consensusClass"
    # add order
    sv_signature <- rbind(sv_signature,results_sig_cluster[[i]]$consensusTree[["order"]])
    row.names(sv_signature)[nrow(sv_signature)] <- "order"
    # add color
    coli <- results_sig_cluster[[i]]$clrs[[1]]
    sv_signature <- rbind(sv_signature,coli)
    row.names(sv_signature)[nrow(sv_signature)] <- "col"
    
    sv_signature_t <- t(sv_signature)
    
    if  (!dir.exists(file.path(output_path,mytitle))){
      dir.create(file.path(output_path,mytitle),recursive = T,showWarnings = F)
    }
    write.table(sv_signature_t,file.path(output_path,mytitle,paste0("sig_cluster",i,".tsv")),quote = F,sep = "\t")
  }
}

run_cluster_function(input=sample_sig_matrix,mytitle = "hc_pearson",myclusterAlg = "hc",mydistance = "pearson",mycluster.num=cluster.num)
run_cluster_function(input=sample_sig_matrix,mytitle = "hc_euclidean",myclusterAlg = "hc",mydistance = "euclidean",mycluster.num=cluster.num)
run_cluster_function(input=sample_sig_matrix,mytitle = "km",myclusterAlg = "km",mydistance = "euclidean",mycluster.num=cluster.num)
run_cluster_function(input=sample_sig_matrix,mytitle = "pam_euclidean",myclusterAlg = "pam",mydistance = "euclidean",mycluster.num=cluster.num)
run_cluster_function(input=sample_sig_matrix,mytitle = "pam_pearson",myclusterAlg = "pam",mydistance = "pearson",mycluster.num=cluster.num)




### evaluation
cluter_evaluation <- data.frame()
cluster_list <- c("Hierarchical_pearson"="hc_pearson",
                     "Hierarchical_euclidean"="hc_euclidean",
                     "Kmeans"="km",
                     "PAM_pearson"="pam_pearson",
                     "PAM_euclidean"="pam_euclidean")
for (j in 1:length(cluster_list)){
  clusteri_name <- names(cluster_list[j])
  for (i in 2:cluster.num){
    filepath <- file.path(output_path,cluster_list[j],paste0("sig_cluster",i,".tsv"))
    cl <- read.delim(filepath)
    cl <- cl$consensusClass
    index <- intCriteria(sample_sig_matrix_numeric,cl,"all")
    dunn <- index$dunn
    silhouette <- index$silhouette
    cindex <- index$c_index
    calinskiharabasz <- index$calinski_harabasz
    daviesbouldin <- index$davies_bouldin
    cluter_evaluation <- rbind(cluter_evaluation,
                               data.frame(
                                 "method"=clusteri_name,
                                 "num.cluster"=i,
                                 "silhouette"=silhouette,
                                 "dunn"=dunn,
                                 "calinskiharabasz"=calinskiharabasz,
                                 "cindex"=cindex,
                                 "daviesbouldin"=daviesbouldin
                               ))
  }
}


### plot
cluster_evaluation_long <- melt(cluter_evaluation,id.vars = c("method","num.cluster"),
                                measure.vars = c("silhouette","dunn","calinskiharabasz",
                                                 "cindex","daviesbouldin"))

measures <- c("Silhouette score",
              "Dunn-index",
              "Calinski-Harabasz score",
              "C-indx",
              "Davies-Bouldin index")
names(measures) <- c("silhouette","dunn","calinskiharabasz",
                       "cindex","daviesbouldin")
ggplot(data=cluster_evaluation_long,
       aes(x=num.cluster,
           y=value,
           color=method))+
  geom_point() + 
  geom_line() + 
  facet_grid(variable ~.,
             scales = "free",
             labeller = labeller(variable = measures))+
  scale_x_continuous(breaks = seq(2,cluster.num,1),
                     labels = seq(2,cluster.num,1))+
  xlab("Num. of cluster")+
  theme(axis.title.y = element_blank())
ggsave(file.path(output_path,"evaluation_of_arg_and_clusternum.png"))
ggsave(file.path(output_path,"evaluation_of_arg_and_clusternum.pdf"),width = 6,height = 7)

# The value of the silhouette coefﬁcient is between [-1, 1]. A score of 1 denotes the best meaning that the data point i is very compact within the cluster to which it belongs and far away from the other clusters.
#Higher value of calinskiharabasz index means the clusters are dense and well separated
#The Dunn Index is the ratio of the smallest distance between observations not in the same cluster to the largest intra-cluster distance. The Dunn Index has a value between zero and infinity, and should be maximized.
#The C-index is limited to the interval [0, 1] and should be minimized.
#Consequently, the number of clusters that minimizes DB is taken as the optimal number of clusters.







