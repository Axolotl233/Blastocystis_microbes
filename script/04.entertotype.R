dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x<=0.0000001,pseudocount,x))
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
pam.medoids=function(x,k) {
  require(cluster)
  medoids = as.vector(pam(as.dist(x), k, diss=TRUE)$medoids)
  return(medoids)
}

library(ade4)
library(cluster)
library(clusterSim)
library(reshape2)
library(vegan)
library(ggstatsplot)
library(class)
library(caret)
library(ape)
library(ggsankey)

set.seed(12345)
out_dir <- paste("output/enterotype/",sp[3],sep = "")
dir.create(out_dir,showWarnings = F,recursive = T)
load("workfile/z.phyloseq.Rdata")

tmp_meta <- metadatadf %>% filter(species == sp[3]) %>% filter(Pattern_class %in% pattern_type[c(1:3)])
tmp_otu <- read.table(in_otu_ab,header = T,row.names = 1)
tmp_otu <- tmp_otu[,tmp_meta$ID]

identical(tmp_meta$ID,colnames(tmp_otu))
tmp_otu_genus <- tmp_otu[grep("g__",row.names(tmp_otu),perl = TRUE),]
tmp_otu_genus <- tmp_otu_genus[grep("s__",row.names(tmp_otu_genus),perl = TRUE,invert = T),]
row.names(tmp_otu_genus) <- gsub(".*\\|","",row.names(tmp_otu_genus))
tmp_otu_list <- read_csv("input/tmp.csv")
tmp_otu_genus <- tmp_otu_genus[row.names(tmp_otu_genus) %in% tmp_otu_list$genus,]
#tmp_otu_genus <- noise_removal(tmp_otu_genus,percent = 0.2,low = 0)
#tmp_otu_genus <- noise_removal(tmp_otu_genus,low = 0.05,method = "mean_cut")

tmp_jsd = dist.JSD(tmp_otu_genus)
tmp_jsd_data <- as.data.frame(as.matrix(tmp_jsd))

#tmp_train_sample <- tmp_meta[tmp_meta$groups %in% group[c(1,2)],2]
#tmp_ex_sample <- tmp_meta[tmp_meta$group1 %in% groupl[c(3)],2]
#tmp_meta_train <- tmp_meta[tmp_meta$sample %in% tmp_train_sample,]
#tmp_jsd_train <- as.dist(as.matrix(tmp_jsd_data[tmp_train_sample,tmp_train_sample]))
#tmp_otu_genus_train <- tmp_otu_genus[,tmp_train_sample]
tmp_nclusters=NULL
#k = 2
for (k in 1:20) { 
  if (k==1) {
    tmp_nclusters[k]=NA 
  } else {
    tmp_data_clusters=pam.clustering(tmp_jsd, k)
    tmp_nclusters[k]=index.G1(t(tmp_otu_genus),tmp_data_clusters,  d = tmp_jsd,
                              centrotypes = "medoids")
  }
}
tmp_nclusters[1] = 0
tmp_ncluster_d <- data.frame(n = 1:20,ch = tmp_nclusters)
ggplot(tmp_ncluster_d,aes(x = n,y = ch))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(1,20,1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x = "k clusters", y = "Calinski-Harabasz index", title = "Optimal number of clusters")
ggsave(paste(out_dir,"opti_clusters_pam.pdf",sep = '/'),width = 6,height = 4)

k_best = which(tmp_nclusters == max(tmp_nclusters), arr.ind = TRUE)
tmp_cluster=pam.clustering(tmp_jsd, k = k_best)
tmp_medoids=pam.medoids(tmp_jsd, k = k_best)
tmp_silhouette=mean(silhouette(tmp_cluster, tmp_jsd)[,3])
cat(tmp_silhouette)

tmp_meta$enterotype <- tmp_cluster
tmp_meta$enterotype <- ifelse(tmp_meta$enterotype == 1,"ET1","ET2")
#tmp_jsd_data <- as.data.frame(as.matrix(tmp_jsd))
identical(row.names(tmp_jsd_data),tmp_meta$ID)
tmp_jsd_data$class <- tmp_meta$enterotype

tmp_pca=dudi.pca(as.data.frame(t(tmp_otu_genus)), scannf=F, nf=10)
tmp_bet=bca(tmp_pca, fac=as.factor(tmp_cluster), scannf=F, nf=2) 
tmp_bet_res <- as.data.frame(t(tmp_bet$tab))
tmp_bet_res$taxon = row.names(tmp_bet_res)

tmp_et1 <- tmp_bet_res[which(tmp_bet_res[,1] == max(tmp_bet_res[,1])),3]
tmp_et2 <- tmp_bet_res[which(tmp_bet_res[,2] == max(tmp_bet_res[,2])),3]
res_et <- c(tmp_et1,tmp_et2) 

tmp_pcoa_res <- pcoa(tmp_jsd)
tmp_pcoa_res_matrix <- tmp_pcoa_res$vectors[,c(1:3)]
colnames(tmp_pcoa_res_matrix) <- c("dim1","dim2","dim3")
identical(tmp_meta$ID,row.names(tmp_pcoa_res_matrix))

tmp_res <- tmp_meta %>% mutate(PC1 = tmp_pcoa_res_matrix[,1],PC2 = tmp_pcoa_res_matrix[,2]) %>% 
  dplyr::select(PC1,PC2,enterotype,Pattern_class,ID)
tmp_res$enterotype <- as.character(tmp_res$enterotype)
write_csv(tmp_meta,paste(out_dir,"entertype.csv",sep = '/'))
tmp_res_f <- tmp_res[tmp_res$Pattern_class %in% pattern_type[c(1:3)],]
tmp_meta_f <- tmp_meta[tmp_meta$Pattern_class %in% pattern_type[c(1:3)],]
tmp_jsd_f <- tmp_jsd_data[tmp_meta_f$ID,tmp_meta_f$ID]
tmp_adonis_res <- adonis2(tmp_jsd_f~enterotype+Pattern_class+age+sex+bmi+parasite_load_log, tmp_meta_f,by = "term")
tmp_res_f$Pattern_class <- factor(tmp_res_f$Pattern_class,levels = pattern_type[c(1:3)])
ggplot(data=tmp_res_f)+
  stat_ellipse(data=tmp_res,geom = "polygon",aes(x=PC1,y=PC2,fill=enterotype),level = 0.9,alpha=0.05)+
  geom_point(aes(x=PC1,y=PC2,color=enterotype,shape = Pattern_class),size=2.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1),
  )+
  labs(x=c("PCoA1"),
       y=c("PCoA2")
  )+
  scale_fill_manual(values = color1)+
  scale_color_manual(values = color1)+
  scale_x_continuous(breaks = seq(-4,4,2))+
  scale_y_continuous(breaks = seq(-4,4,2))
  
ggsave(paste(out_dir,"pcoa.enterotype.pdf",sep = '/'),height = 5,width = 6)


plot_l <- list()
for(i in 1:length(res_et)){
  tmp_d <- as.data.frame(t(tmp_otu_genus))
  tmp_d$sample = row.names(tmp_d)
  tmp_d <- tmp_d[,c("sample",res_et[i])]
  colnames(tmp_d) <- c("ID","value")
  tmp_plot <- inner_join(tmp_res,tmp_d,by = "ID")
  p1 <- ggplot(tmp_plot, aes(x=enterotype, y =value))+
    geom_boxplot(aes(color=enterotype),outlier.shape = NA)+
    geom_jitter(aes(color=enterotype),width = 0.15,height = 0,size = 3,alpha = 0.7)+
    stat_compare_means(method="wilcox.test",
                       comparisons = list(c("ET1","ET2")),
                       label = "p.format"
    )+
    theme_light()+
    theme(legend.position="none", 
          axis.title.x = element_blank(),
          panel.grid = element_blank()
    )+
    scale_color_manual(
      values = color1[c(1,2)]
    )+
    labs(y = "Relative Abdunance",title = res_et[i])
  plot_l[[i]] <- p1
}
plot_grid(plot_l[[1]],plot_l[[2]],nrow = 1)
ggsave(paste(out_dir,"entertotype.sp.pdf",sep = '/'),width = 6,height = 6)

ggstatsplot::ggbarstats(
  data = tmp_meta_f,
  x = Pattern_class,
  y = enterotype
)
ggsave(paste(out_dir,"entertotype.ratio2.pdf",sep = '/'),width =6,height = 8)
tmp_sankey_p <- tmp_meta_f %>% dplyr::select(Pattern_class,enterotype) %>% make_long(enterotype,Pattern_class)
ggplot(tmp_sankey_p, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = 0.3,width = 0.2)+
  geom_alluvial_text(size = 3, color = "white")+
  theme_sankey()+
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank())+
  scale_fill_manual(
    values = c(color1[c(1,2)],color0[c(1,2,3)])
  )+
  labs(x = "")
ggsave(paste(out_dir,"entertotype.ratio.pdf",sep = '/'),width = 7,height = 5)

ggplot(tmp_meta_f, aes(x=enterotype, y =parasite_load_log)) +
  geom_boxplot(aes(color=enterotype,),outlier.shape = NA)+
  geom_jitter(aes(color=enterotype),width = 0.15,height = 0,size = 3,alpha = 0.7)+
  stat_compare_means(method="wilcox.test",
                     comparisons = list(c("ET1","ET2")),
                     label = "p.format"
  )+
  facet_wrap(.~Pattern_class)+
  theme_light()+
  theme(legend.position="none", 
        axis.title.x = element_blank(),
        panel.grid = element_blank()
  )+
  scale_color_manual(
    values = color1[c(1,2)]
  )
ggsave(paste(out_dir,"entertotype.parasite_load_pattern.pdf",sep = '/'),width = 6,height = 4)

ggplot(tmp_meta_f, aes(x=enterotype, y =parasite_load_log)) +
  geom_boxplot(aes(color=enterotype,),outlier.shape = NA,width = 0.7)+
  geom_jitter(aes(color=enterotype),width = 0.15,height = 0,size = 3,alpha = 0.7)+
  stat_compare_means(method="wilcox.test",
                     comparisons = list(c("ET1","ET2")),
                     label = "p.format"
  )+
  theme_light()+
  theme(legend.position="none", 
        axis.title.x = element_blank(),
        panel.grid = element_blank()
  )+
  scale_color_manual(
    values = color1[c(1,2)]
  )
ggsave(paste(out_dir,"entertotype.parasite_load.pdf",sep = '/'),width = 3,height = 4)
tmp_res_f %>% group_by(enterotype,Pattern_class) %>% summarise(n())
