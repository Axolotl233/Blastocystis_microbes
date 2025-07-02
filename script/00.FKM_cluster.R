c_value_cluster <- function(x,k,prefix){
  
  dir.create(prefix)
  out_dir = paste(prefix,k,sep = "/")
  dir.create(out_dir)
  distanceMat<-as.dist((1-cor(t(x),method="pearson"))/2)
  initClust <- fanny(distanceMat, k, diss = TRUE, memb.exp= 1.02,maxit=1000)
  initClustMem <- as.data.frame(initClust$membership)
  initClustMem$sample = row.names(initClustMem)
  
  #clusterCutoff <- initClust$membership[sort(initClust$membership,index.return=TRUE,decreasing=TRUE)$ix[nrow(initClust$membership)]]
  clusterCutoff <- 0.9
  cantAssignSample<-(apply(initClust$membership,1,max)<clusterCutoff)
  goodSamples<-(apply(initClust $membership,1,max)>=clusterCutoff)
  initClust$cluster[cantAssignSample]<-0
  initClust$membership[cantAssignSample,]<-0
  
  clustProfiles <- getClustProfiles(x,initClust,clusterCutoff)
  whichClusters <- clustProfiles$whichClust
  clustProfiles <- clustProfiles$profiles
  profileCor <- cor(t(clustProfiles))
  hClust <- hclust(as.dist(1-profileCor),method="single")
  #plot(hClust, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "", cex = 0.8)
  collapseGrp <- cutree(hClust,h=0.2)
  finalProfiles<-makeCollapsedProfiles(x,
                                       clustProfiles, 
                                       clusterCutoff,
                                       collapseGrp,
                                       whichClusters,
                                       initClust)
  rownames(finalProfiles)<- paste("Pattern_",1:nrow(finalProfiles),sep="")
  colnames(finalProfiles)<-colnames(x)
  
  SampleToPatternCor <- correlateSampleToProfiles(x,finalProfiles)
  finalRes <- as.data.frame(GetMemberList(SampleToPatternCor,0.85))
  finalProfiles <- as.data.frame(finalProfiles)
  
  membershipData <- initClust$membership
  PC <- mean(rowSums(membershipData^2))
  PE <- -mean(rowSums(membershipData * log(membershipData)))
  
  clusterAssignments <- apply(initClust$membership, 1, which.max)
  sil <- silhouette(clusterAssignments, distanceMat)
  #plot(sil, main = "Silhouette plot for FKM clustering")
  mean_sil_width <- mean(sil[, 3])
  
  n <- nrow(x)
  n_iter <- 100
  prop_sample <- 0.8
  consensus_matrix <- matrix(0, n, n)
  consensus_matrix <- as.data.frame(consensus_matrix)
  rownames(consensus_matrix) <- initClustMem$sample
  colnames(consensus_matrix) <- initClustMem$sample
  for (i in 1:n_iter) {
    #idx <- sort(sample(1:n, size = floor(n * prop_sample), replace = FALSE))
    idx <- sample(1:n, n, replace = FALSE)
    data_shuffled <- x[idx, ]
    distanceMatShuffled<-as.dist((1-cor(t(data_shuffled),method="pearson"))/2)
    ShuffledClust <- fanny(distanceMatShuffled, k, diss = TRUE, memb.exp= 1.02,maxit=1000)
    #ShuffledClustMem <- as.data.frame(ShuffledClust$membership)
    #ShuffledClustMem$sample <- row.names(ShuffledClustMem)
    hard_labels <- apply(ShuffledClust$membership, 1, which.max)
    for (a in 1:length(idx)) {
      for (b in 1:length(idx)) {
        if (hard_labels[a] == hard_labels[b]) {
          consensus_matrix[names(hard_labels[a]), names(hard_labels[b])] <- 
            consensus_matrix[names(hard_labels[a]), names(hard_labels[b])] + 1
        }
      }
    }
  }
  consensus_matrix <- consensus_matrix / n_iter
  #pheatmap(consensus_matrix, main = "Consensus Matrix")
  mean_consensus <- mean(consensus_matrix[lower.tri(consensus_matrix)])
  #cat("Mean Consensus Index:", mean_consensus, "\n")
  
  out_file1 <- paste(out_dir,"patternOutputFile.txt",sep = "/")
  write_tsv(file=out_file1,x=finalProfiles)
  out_file2 <-  paste(out_dir,"patternsample.txt",sep = "/")
  write_tsv(file=out_file2,x=finalRes)
  
  list_return <- list(finalRes,PC,PE,mean_sil_width,consensus_matrix)
  return(list_return)
}
batch_cvalue <- function(d,k,value_col,pre){
  data_return = data.frame(
    k = "",
    clustered_sample = "",
    PC = "",
    PE = "",
    sil_width = ""
  )
  n = 1
  data_cdf = data.frame()
  for(kcluster in k){
    outdir_o <- pre
    d_tmp <- d[,value_col]
    tmp_return <-  c_value_cluster(d_tmp,kcluster,outdir_o)
    tmp_pattern <- tmp_return[[1]]
    colnames(tmp_pattern) <- c("ID","Pattern")
    d_tmp_tern <- left_join(d,tmp_pattern,by = "ID")
    d_tmp_tern[is.na(d_tmp_tern$Pattern),"Pattern"] <- "Unclass"
    d_tmp_tern$parasite_load <- log2(d_tmp_tern$parasite_load)
    p <- ggtern(d_tmp_tern,aes(x = ST1, y = ST2, z = ST3))+
      geom_mask()+
      geom_point(aes(shape = Pattern,color = Pattern),size = 5,alpha=0.8)+
      theme_rgbw()+
      scale_color_manual(
        values = color0_1
      )
    tern_plot_n <- paste(outdir_o,kcluster,"tern.pdf",sep = "/")
    ggsave(tern_plot_n,plot = p,width = 7,height = 5)
    tern_table_n <- paste(outdir_o,kcluster,"tern.csv",sep = "/")
    write_csv(x = d_tmp_tern,file = tern_table_n)
    heatmap_plot <- paste(outdir_o,kcluster,"consensus_heatmap.pdf",sep = "/")
    pdf(heatmap_plot,width = 15,height = 15)
    pheatmap(tmp_return[[5]], main = paste("Consensus Matrix (k =", kcluster, ")"),
             show_rownames = T, 
             show_colnames = T)
    dev.off()
    
    consensus_vals <- tmp_return[[5]][lower.tri(tmp_return[[5]])]
    x_vals <- seq(0, 1, 0.01)
    y_vals <- ecdf(consensus_vals)(x_vals)
    tmp_data_cdf = data.frame(
      x = x_vals,
      y = y_vals,
      k = kcluster
    )
    data_cdf = rbind(data_cdf,tmp_data_cdf)
    cdf_plot <- paste(outdir_o,kcluster,"CDF.pdf",sep = "/")
    pdf(cdf_plot,width = 8,height = 8)
    plot(x_vals, y_vals, type = "l", lwd = 2,
         xlab = "Consensus Index", ylab = "Cumulative Proportion",
         main = "Consensus CDF Curve")
    dev.off()
    
    data_return[n,] = c(kcluster,nrow(tmp_pattern),tmp_return[[2]],tmp_return[[3]],tmp_return[[4]])
    n = n + 1
  }
  data_return <- as.data.frame(apply(data_return, c(1,2), as.numeric))
  return(list(data_return,data_cdf))
}

library(cluster)
library(NbClust)
library(ggtern)
library(clue)
library(factoextra)
library(vegan)
library(pheatmap)
source("script/00.FKM_cluster.function.R")
out_dir = "output/subtype_class"
dir.create(out_dir,showWarnings = F)

d1 <- read_csv(in_meta)
d2 <- as.data.frame(t(read.csv(in_amp,head = F,row.names = 1)))

d_1 <- left_join(d1,d2,by = "ID") %>% 
  select(ID,species,parasite_load,ST1,ST2,ST3,ST5)
d_1[,4:7] <- apply(d_1[,4:7],2,function(x){as.numeric(x)})

## filter_ST5 sample
d_1 %<>% drop_na(parasite_load,ST1) %>% 
   select(ID,parasite_load,species,ST1,ST2,ST3,ST5)
d_1[,4:7] <- sweep(d_1[,4:7],1,apply(d_1[,4:7],1,sum),"/")
d_1[,4:7] = 100 * d_1[,4:7]

#d_1[,4:6] <- sweep(d_1[,4:6],1,apply(d_1[,4:6],1,sum),"/")
#d_1[,4:6] = 100 * d_1[,4:6]
#data_maf <- d_12 %>% filter(species == "mfas")
#data_mmu <- d_12 %>% filter(species == "mmul")
#data_mmu %<>% select(!ST5)
d_12 <- as.data.frame(d_1)
row.names(d_12) <- d_12$ID
d_cluster_res <- batch_cvalue(d_12,c(2:7),c(4:7),paste(out_dir,"cvalue",sep = "/"))
save(d_cluster_res,file = "workfile/FKM_cluster.Rdata")
d_cdf <- d_cluster_res[[2]]
d_cdf$k <- as.character(d_cdf$k)
ggplot(d_cdf)+
  geom_line(aes(x = x,y = y,color = k),size = 1)+
  theme_cowplot()+
  labs(x = "Consensus Index", y = "Cumulative Proportion",title = "Consensus CDF Curve" )
ggsave(paste(out_dir,"Consensus_CDF.pdf",sep = "/"),height = 5,width = 7)

k_used <- unique(d_cdf$k)
d_area = matrix(NA,nrow = length(k_used),ncol = 2)
for(i in 1:length(k_used)){
  d_cdf_tmp <- d_cdf[d_cdf$k == k_used[i],]
  d_area_tmp <- sum(d_cdf_tmp$y) / nrow(d_cdf_tmp)
  d_area[i,] = c(k_used[i],d_area_tmp)
}
d_area <- as.data.frame(d_area)
d_area[,2] <- as.numeric(as.vector(d_area[,2]))
d_area_d <- diff(d_area$V2)
pdf(paste(out_dir,"Consensus_DAP.pdf",sep = "/"),height = 5,width = 7)
plot(c(3,4,5,6,7), d_area_d, type = "b", pch = 16, col = "blue",
     xlab = "Number of Clusters (k)", ylab = "Delta Area",
     main = "Delta Area Plot (Consensus CDF)")
dev.off()
write_csv(d_cluster_res[[1]],file = paste(out_dir,"cluster.csv",sep = "/"))
