library(vegan)
library(aPCoA)

set.seed(10001)
out_dir <- c("output/sample_stratum_gr1")
dir.create(out_dir,showWarnings = F,recursive = T)
load("workfile/z.phyloseq.Rdata")

phylo_used <- subset_samples(phyloseqin, species == sp[3])
phylo_used <- subset_samples(phylo_used,Pattern_class %in% pattern_type[c(1:3)])
meta_used <- metadatadf %>% filter(species == sp[3],Pattern_class %in% pattern_type[c(1:3)])
meta_used$Pattern_class <- factor(meta_used$Pattern_class,levels = pattern_type[c(1:3)])
meta_used %>% group_by(parasite_load_gr3) %>% summarise(n())

phylo_used_f <- filter_taxa(phylo_used, function(x) sum(x > 0.005) > (0.2*length(x)), TRUE)
data_otu <- as.data.frame(t(phylo_used_f@otu_table))
data_dist = as.matrix(vegdist(data_otu, method = 'bray'))

tmp_apcoa_res <- aPCoA(data_dist ~ age + sex + bmi,
                       meta_used,maincov = Pattern_class)
tmp_apcoa_matrix <- data.frame(tmp_apcoa_res$plotMatrix) %>% select(X1,X2)
tmp_apcoa_matrix$ID <- row.names(tmp_apcoa_matrix)
tmp_apcoa_matrix_p <- inner_join(tmp_apcoa_matrix,meta_used,by = "ID")
tmp_apcoa_matrix_p$Pattern_class <- factor(tmp_apcoa_matrix_p$Pattern_class,levels = pattern_type)
tmp_list <- combine_list(tmp_apcoa_matrix_p$Pattern_class)
a <- ggplot(tmp_apcoa_matrix_p,aes(x=Pattern_class,y=X1)) + 
  geom_boxplot(aes(color=Pattern_class),outlier.shape = NA)+
  geom_jitter(aes(color=Pattern_class),size = 2.5,height=0,width = 0.2,alpha = 0.7)+
  theme_light()+
  theme(axis.title.x=element_blank(), 
        legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_color_manual(values = color0)+
  labs(y="PCoA1")+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")+
  facet_wrap(.~parasite_load_gr3)

tmp_formula = y ~ x
tmp_apcoa_matrix_p$parasite_load_gr3 <- gsub("Q","",tmp_apcoa_matrix_p$parasite_load_gr3)
tmp_apcoa_matrix_p$parasite_load_gr3 <- as.numeric(tmp_apcoa_matrix_p$parasite_load_gr3)
b <- ggplot(tmp_apcoa_matrix_p) +
  stat_smooth(aes(x=parasite_load_gr3,y=X1),method='lm',formula = tmp_formula,color = "navy",alpha = 0.35)+
  geom_point(aes(x=parasite_load_gr3,y=X1,color = Pattern_class),size=3,alpha =0.7)+
  stat_cor(aes(x=parasite_load_gr3,y=X1), 
           cor.coef.name = "rho",method = "spearman",label.y.npc = 0.1,label.x.npc = 0.05)+
  theme_light()+
  theme(panel.grid = element_blank(),
        legend.position = "none"
  )+
  scale_x_continuous(breaks = c(1,2,3),labels = c("Q1","Q2","Q3"))+
  scale_color_manual(values = color0)+
  labs(x = "Parasite_load Group", y = "PCoA1")
plot_grid(a,b,nrow = 1,rel_widths = c(1.3,1))
ggsave(paste(out_dir,"stratum_PCoA1.pdf",sep = '/'),width = 11,height = 5)

model <- lm(X1~parasite_load_gr3 + sex + age + bmi, data = tmp_apcoa_matrix_p)
summary(model)

tmp_apcoa_matrix_p$parasite_load_gr3 <- paste("Q",tmp_apcoa_matrix_p$parasite_load_gr3,sep = "")
tmp_gr <- sort(unique(tmp_apcoa_matrix_p$parasite_load_gr3))
tmp_res_adonis_all <- data.frame()
for(i in 1:length(tmp_gr)){
  tmp_d <- tmp_apcoa_matrix_p %>% filter(parasite_load_gr3 == tmp_gr[i])
  tmp_dist <- data_dist[tmp_d$ID,tmp_d$ID]
  tmp_res_adonis = adonis2(tmp_dist~Pattern_class + age + sex + bmi + 
                             parasite_load_log ,tmp_d,by = "term",permutations = 999)
  tmp_res_adonis$class <- row.names(tmp_res_adonis)
  tmp_res_adonis$gr <- tmp_gr[i]
  tmp_res_adonis_all <- bind_rows(tmp_res_adonis,tmp_res_adonis_all)
}
tmp_res_adonis_all <- tmp_res_adonis_all %>% filter(class == "Pattern_class")
ggplot(tmp_res_adonis_all)+
  geom_bar(aes(x = gr, y = R2,fill = gr),stat = "identity",width = 0.6)+
  geom_text(aes(x = gr, y = R2+0.005,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = color2[c(3:5)])+
  theme(legend.position = "none")+
  labs(x = "")
ggsave(paste(out_dir,"stratum_adonis2.pdf",sep = '/'),width = 5,height = 5)

#data_dist_pair <- data_dist
#data_dist_pair[upper.tri(data_dist_pair)] <- NA
sample_list <- combine_list(meta_used$ID)
data_dist_pair <- data.frame(matrix(NA,nrow = length(sample_list),ncol = 4))
colnames(data_dist_pair) <- c("sample1","sample2","dist","delta")
for(i in 1:length(sample_list)){
  tmp_para1 <- meta_used[meta_used$ID == sample_list[[i]][1],13]
  tmp_para2 <- meta_used[meta_used$ID == sample_list[[i]][2],13]
  tmp_delta <- abs(tmp_para1-tmp_para2)
  data_dist_pair[i,] = c(sample_list[[i]][1],sample_list[[i]][2],
                         data_dist[sample_list[[i]][1],sample_list[[i]][2]],
                         tmp_delta
  )
}
data_dist_pair$dist <- as.numeric(data_dist_pair$dist)
data_dist_pair$delta <- as.numeric(data_dist_pair$delta)
ggplot(data_dist_pair) + 
  stat_smooth(aes(x=dist,y=delta),method='lm',formula = tmp_formula,color = "navy",alpha = 0.35)+
  ggrastr::geom_point_rast(aes(x=dist,y=delta),color="grey80",size=2,alpha =0.3)+
  stat_cor(aes(x=dist,y=delta),cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.35)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid = element_blank(),
        legend.position = "none"
  )+
  scale_color_manual(values = color0)
ggsave(paste(out_dir,"cor.bc_to_paraload_diff.pdf",sep = '/'),width = 4,height = 4)

plot_l <- list()
data_dist_pair_gr <- data_frame()
for(i in 1:length(tmp_gr)){
  tmp_meta <- meta_used %>% filter(parasite_load_gr3 == tmp_gr[i])
  tmp_dist <- data_dist[tmp_d$ID,tmp_d$ID]
  tmp_sample_list <- combine_list(tmp_meta$ID)
  tmp_data_dist_pair <- data.frame(matrix(NA,nrow = length(tmp_sample_list),ncol = 4))
  colnames(tmp_data_dist_pair) <- c("sample1","sample2","dist","delta")
  for(j in 1:length(tmp_sample_list)){
    tmp_para1 <- tmp_meta[tmp_meta$ID == tmp_sample_list[[j]][1],13]
    tmp_para2 <- tmp_meta[tmp_meta$ID == tmp_sample_list[[j]][2],13]
    tmp_delta <- abs(tmp_para1-tmp_para2)
    tmp_data_dist_pair[j,] = c(tmp_sample_list[[j]][1],tmp_sample_list[[j]][2],
                           data_dist[tmp_sample_list[[j]][1],tmp_sample_list[[j]][2]],
                           tmp_delta
    )
  }
  tmp_data_dist_pair$dist <- as.numeric(tmp_data_dist_pair$dist)
  tmp_data_dist_pair$delta <- as.numeric(tmp_data_dist_pair$delta)
  p <- ggplot(tmp_data_dist_pair) + 
    stat_smooth(aes(x=dist,y=delta),method='lm',formula = tmp_formula,color = "navy",alpha = 0.35)+
    ggrastr::geom_point_rast(aes(x=dist,y=delta),color="grey80",size=2,alpha =0.3)+
    stat_cor(aes(x=dist,y=delta),cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0.35)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
          legend.position = "none"
    )+
    labs(title = tmp_gr[i])
  tmp_data_dist_pair$gr <- tmp_gr[i]
  data_dist_pair_gr <- bind_rows(tmp_data_dist_pair,data_dist_pair_gr)
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],plot_l[[3]],nrow = 1)
ggsave(paste(out_dir,"cor.bc_to_paraload_diff.stratum.pdf",sep = '/'),width = 11,height = 4)

data_dist_pair$class = NA
for(i in 1:nrow(data_dist_pair)){
  tmp_s1 <- data_dist_pair[i,1]
  tmp_s2 <- data_dist_pair[i,2]
  tmp_gr1 <- meta_used[meta_used$ID == tmp_s1,"parasite_load_gr3"]
  tmp_gr2 <- meta_used[meta_used$ID == tmp_s2,"parasite_load_gr3"]
  data_dist_pair[i,5] <- paste(sort(c(tmp_gr1,tmp_gr2)),collapse = "-")
}

group_list <- list(
  c("Q1-Q1","Q1-Q2","Q1-Q3"),
  c("Q2-Q2","Q1-Q2","Q2-Q3"),
  c("Q3-Q3","Q1-Q3","Q2-Q3")
)
tmp_sample <- meta_used[meta_used$Pattern_class %in% pattern_type[c(1,3)],1]
data_dist_pair <- data_dist_pair %>% 
  filter(sample1 %in% tmp_sample) %>%
  filter(sample2 %in% tmp_sample)
tmp_gr <- sort(unique(meta_used$parasite_load_gr3))
plot_l <- list()
for(i in 1:length(group_list)){
  tmp_d <- data_dist_pair %>% filter(class %in% group_list[[i]])
  tmp_d$class <- factor(tmp_d$class,levels = group_list[[i]])
  tmp_lll <- list(c(group_list[[i]][1],group_list[[i]][2]),
                  c(group_list[[i]][1],group_list[[i]][3]))
  p <- ggplot(data = tmp_d,aes(x=class,y=dist))+
    geom_violin(aes(fill = class),color = "white")+
    geom_boxplot(width = 0.3)+
    stat_compare_means(method="wilcox.test",
                       comparisons = tmp_lll,
                       label = "p.format"
    )+
    labs(y = "Bray-Curtis distance")+
    theme_cowplot()+
    theme(legend.position="none", 
          strip.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    scale_fill_manual(
      values = c("#9abf82","grey85","grey85")
    )+
    labs(x = "",title = tmp_gr[i])+
    scale_y_continuous(limits = c(0,1.2))
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],plot_l[[3]],align = "hv",nrow = 1)
data_dist_pair_stat <- summary_stat(data_dist_pair,3,5)
ggsave(paste(out_dir,"distance.stratum.pdf",sep = '/'),width = 12.5,height = 4)
write_csv(data_dist_pair_stat,paste(out_dir,"distance.stratum_stat.csv",sep = '/'))
