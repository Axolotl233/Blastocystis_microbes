merge_samples_mean <- function(physeq, group){
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  merged <- merge_samples(physeq, group)
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}
plot_correalation_alpha <- function(dd,class_custom){
  tmp_formula <- y ~ x
  ddd <- dd %>% dplyr::select(type,value,key,class_custom)
  colnames(ddd) <- c("type","value","key","class")
  p <- ggplot(ddd) + 
    stat_smooth(aes(x=class,y=value),
                method='lm',formula = tmp_formula,color = "grey70")+
    geom_point(aes(x=class,y=value,color=type),size=3.5,alpha = 0.8)+
    scale_color_manual(
      values = color0
    )+
    stat_cor(aes(x=class,y=value), 
             cor.coef.name = "rho",method = "spearman",label.y.npc = 0.1,label.x.npc = 0.05)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()
    )+
    facet_wrap(.~key,scales = "free")+
    labs(title = class_custom)
  return(p)
}

library(vegan)
library(aPCoA)
library(factoextra)
library(microbiome)
library(microbiomeSeq)
library(pheatmap)
library(Maaslin2)
library(compositions)
library(mediation)
set.seed(10001)

out_dir <- paste("output/base_analysis/",sp[3],sep = "")
dir.create(out_dir,showWarnings = F,recursive = T)
load("workfile/z.phyloseq.Rdata")

phylo_used <- subset_samples(phyloseqin, species == sp[3])
phylo_used <- subset_samples(phylo_used,Pattern_class %in% pattern_type[c(1:3)])
meta_used <- metadatadf %>% filter(species == sp[3],Pattern_class %in% pattern_type[c(1:3)])
meta_used$Pattern_class <- factor(meta_used$Pattern_class,levels = pattern_type[c(1:3)])
phylo_used <- filter_taxa(phylo_used, function(x) sum(x > 0) > 0, TRUE)

data_bar <- merge_samples_mean(phylo_used, "Pattern_class")
data_bar_count <- 
  transform_sample_counts(data_bar, function(x) x/sum(x) * 100)
data_bar_count <-tax_glom(data_bar_count, 'Phylum', NArm = T)
data_bar_count@sam_data$Pattern_class <- factor(data_bar_count@sam_data$Pattern_class,
                                                levels = pattern_type[c(1:3)])

plot_bar(data_bar_count, fill='Phylum')+
  geom_bar(stat='identity',position='stack')+
  #scale_fill_manual(values = color19)+
  scale_fill_manual(values = color_bar)+
  ylim(c(0,101)) + theme_cowplot()+ theme(legend.key=element_blank())
ggsave(paste(out_dir,"phylum_bar.pdf",sep = '/'),width = 6,height = 6)

phylo_phylum <- tax_glom(phylo_used,taxrank = "Phylum")
data_otu_p <- as.data.frame(otu_table(phylo_phylum))
data_otu_p <- sweep(data_otu_p,2,colSums(data_otu_p),'/')*100
data_tax <- as.data.frame(tax_table(phylo_phylum))
identical(row.names(data_otu_p),row.names(data_otu_p))
row.names(data_otu_p) <- data_tax$Phylum

tmp_tax_p <- c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Spirochaetes")
data_otu_p <- data_otu_p %>% mutate(taxon = row.names(data_otu_p)) %>%
  filter(taxon %in% tmp_tax_p) %>%
  dplyr::select(-taxon)
data_otu_p <- as.data.frame(t(data_otu_p))
data_otu_p$ID <- row.names(data_otu_p)
meta_tmp <- meta_used %>% dplyr::select(ID,Pattern_class)
data_otu_pl <- inner_join(data_otu_p,meta_tmp,by = "ID") %>% 
  gather(key = "taxon",value = "value", tmp_tax_p)
data_otu_pl$Pattern_class <- factor(data_otu_pl$Pattern_class ,levels = pattern_type[c(1:3)])

plot_l = list() 
for(i in 1:length(tmp_tax_p)){
  tmp_data <- data_otu_pl %>% filter(taxon == tmp_tax_p[i])
  tmp_list <- combine_list(tmp_data$Pattern_class)
  p <- ggplot(tmp_data,aes(x=Pattern_class,y=value,color = Pattern_class))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(size = 2,width = 0.2,height = 0,alpha = 0.7)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = tmp_tax_p[i])+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color0)
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],plot_l[[5]],
          plot_l[[3]],plot_l[[4]],nrow = 2)
ggsave(paste(out_dir,"phylum_top5.pdf",sep = '/'),width = 11,height = 8)
data_otu_p$BF_ratio <- data_otu_p$Bacteroidetes/data_otu_p$Firmicutes
data_otu_p2 <- data_otu_p %>% dplyr::select(ID,BF_ratio)
data_otu_p2 <- inner_join(data_otu_p2,meta_tmp,by = "ID")
ggplot(data_otu_p2,aes(x=Pattern_class,y=BF_ratio,color = Pattern_class))+
  geom_boxplot(outlier.shape = NA,)+
  geom_jitter(size = 4,width = 0.15,height = 0,alpha = 0.7)+
  theme_cowplot()+
  labs(x="Group",y = "Values",title="B/F_ratio")+
  stat_compare_means(comparisons = tmp_list,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(values = color0)
ggsave(paste(out_dir,"phylum_bf_ratio.pdf",sep = '/'),width = 5,height = 5)

#=====> diversity analysis
data_otu <- as.data.frame(phylo_used@otu_table)
tmp_sample <- phylo_used@sam_data$ID
tmp_group <- factor(phylo_used@sam_data$Pattern_class,levels = pattern_type[c(2,1,3)])
tmp_sex <- phylo_used@sam_data$sex
tmp_age <- phylo_used@sam_data$age
tmp_bmi <- phylo_used@sam_data$bmi
tmp_paraload <- phylo_used@sam_data$parasite_load_log
tmp_paraload_ST1 <- phylo_used@sam_data$ST1_abs_log
tmp_paraload_ST2 <- phylo_used@sam_data$ST2_abs_log
tmp_paraload_ST3 <- phylo_used@sam_data$ST3_abs_log
tmp_ST1 <- phylo_used@sam_data$ST1
tmp_ST2 <- phylo_used@sam_data$ST2
tmp_ST3 <- phylo_used@sam_data$ST3

tmp_shannon <- vegan::diversity(t(data_otu),index = "shannon")
tmp_simpson <-  vegan::diversity(t(data_otu),index = "simpson")
tmp_invsimpson <- 1/tmp_simpson

data_alpha_raw <- data.frame(sample=tmp_sample
                            ,type=tmp_group
                            ,sex=tmp_sex
                            ,age=tmp_age
                            ,bmi=tmp_bmi
                            ,paraload_total = tmp_paraload
                            ,paraload_ST1 = tmp_paraload_ST1
                            ,paraload_ST2 = tmp_paraload_ST2
                            ,paraload_ST3 = tmp_paraload_ST3
                            ,ST1 = tmp_ST1
                            ,ST2 = tmp_ST2
                            ,ST3 = tmp_ST3
                            ,shannon=tmp_shannon
                            ,simpson=tmp_simpson
                            ,invsimpson=tmp_invsimpson
)

tmp_observe <- data.frame(sample = "",observed ="")
for(i in 1:ncol(data_otu)){
  tmp_sam = colnames(data_otu)[i]
  tmp_ob = length(which(data_otu[,i] != 0))
  tmp_observe[i,] = c(tmp_sam,tmp_ob)
}
data_alpha_raw <- inner_join(data_alpha_raw,tmp_observe,by = "sample")
data_alpha_raw$observed <- as.numeric(data_alpha_raw$observed)
write_csv(x = data_alpha_raw,paste(out_dir,"alpha_diversity.csv",sep = '/'))

phylo_used2 <- subset_samples(phyloseqin_count, Pattern_class %in% pattern_type[c(2,1,3)])
richness_index <- c("Observed","Shannon","Simpson","InvSimpson")
tmp_other <- estimate_richness(phylo_used2,measures = richness_index)
tmp_other$sample <- row.names(tmp_other)
data_alpha_raw2 <- data.frame(sample=tmp_sample
                             ,type=tmp_group
                             ,sex=tmp_sex
                             ,age=tmp_age
                             ,bmi=tmp_bmi
                             ,paraload_total = tmp_paraload
                             ,paraload_ST1 = tmp_paraload_ST1
                             ,paraload_ST2 = tmp_paraload_ST2
                             ,paraload_ST3 = tmp_paraload_ST3
                             ,ST1 = tmp_ST1
                             ,ST2 = tmp_ST2
                             ,ST3 = tmp_ST3
)
data_alpha_raw2 <- left_join(data_alpha_raw2,tmp_other,by = "sample")
write_csv(x = data_alpha_raw2,paste(out_dir,"alpha_diversity_phyloseq.csv",sep = '/'))

tmp_alpha_raw_plot <- data_alpha_raw %>% dplyr::select(sample,type,
                                                       paraload_total,paraload_ST1,paraload_ST2,paraload_ST3,
                                                       ST1,ST2,ST3,
                                                       observed,shannon,simpson) %>%
  gather(key = "key",value = "value",observed,shannon,simpson)
tmp_alpha_raw_plot$type <- factor(tmp_alpha_raw_plot$type,levels = pattern_type)
tmp_list <- combine_list(tmp_alpha_raw_plot$type)
ggplot(tmp_alpha_raw_plot,aes(x=type,y=value)) + 
  geom_boxplot(aes(color=type),outlier.shape = NA)+
  geom_jitter(aes(color=type),size = 2.5,height=0,width = 0.2,alpha = 0.7)+
  theme_light()+
  theme(axis.title.x=element_blank(), 
        legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_color_manual(values = color0)+
  labs(y="Observed value")+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")+
  facet_wrap(.~key,scales = "free")
ggsave(paste(out_dir,"alpha_raw.pdf",sep = '/'),width = 8,height = 5)

row.names(data_alpha_raw) <- data_alpha_raw$sample
masslin_data <- as.data.frame(t(data_alpha_raw[,c(13:16)]))
Maaslin2(
  input_data = masslin_data, 
  input_metadata = data_alpha_raw, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM" ,
  output = paste(out_dir,"alpha_cov",sep = '/'), 
  fixed_effects = c("type","sex","age","bmi","paraload_total"),
  reference = c("type,ST3_dom")
)

plot_correalation_alpha(tmp_alpha_raw_plot,"paraload_total")
ggsave(paste(out_dir,"cor.paraload_total_alpha.pdf",sep = '/'),width = 12.5,height = 4)

plot_correalation_alpha(tmp_alpha_raw_plot,"paraload_ST3")
ggsave(paste(out_dir,"cor.paraload_ST3_alpha.pdf",sep = '/'),width = 12.5,height = 4)

plot_correalation_alpha(tmp_alpha_raw_plot,"paraload_ST2")
ggsave(paste(out_dir,"cor.paraload_ST3_alpha.pdf",sep = '/'),width = 12.5,height = 4)

plot_correalation_alpha(tmp_alpha_raw_plot,"paraload_ST1")
ggsave(paste(out_dir,"cor.paraload_ST3_alpha.pdf",sep = '/'),width = 12.5,height = 4)

plot_correalation_alpha(tmp_alpha_raw_plot,"ST3")
ggsave(paste(out_dir,"cor.ST3_alpha.pdf",sep = '/'),width = 12.5,height = 4)

plot_correalation_alpha(tmp_alpha_raw_plot,"ST2")
ggsave(paste(out_dir,"cor.ST3_alpha.pdf",sep = '/'),width = 12.5,height = 4)

plot_correalation_alpha(tmp_alpha_raw_plot,"ST1")
ggsave(paste(out_dir,"cor.ST3_alpha.pdf",sep = '/'),width = 12.5,height = 4)

## beta
#tmp_otu2 <- noise_removal(tmp_otu,percent = 0.1,low = 0.001)
phylo_pca_f <- filter_taxa(phylo_used, function(x) sum(x > 0.005) > (0.2*length(x)), TRUE)
data_pca <- as.data.frame(t(phylo_pca_f@otu_table))
tmp_list <- combine_list(meta_used$Pattern_class)
res_adonis_list <- list()
for(i in 1:length(tmp_list)){
  tmp_meta2 = meta_used %>% filter(Pattern_class %in% tmp_list[[i]])
  tmp_data2 <- data_pca[row.names(data_pca) %in% tmp_meta2$ID,]
  tmp_dist2 <- as.matrix(vegdist(tmp_data2, method = 'bray'))
  tmp_adonis_result_dis2 = adonis2(tmp_dist2 ~ Pattern_class, tmp_meta2,by = "term")
  res_adonis_list[[i]] <- tmp_adonis_result_dis2
}
data_list_adonis <- data.frame()
for(i in 1:length(tmp_list)){
  tmp_d <- as.data.frame(res_adonis_list[[i]])
  tmp_d$item <- row.names(tmp_d)
  tmp_d$class <- paste(tmp_list[[i]],collapse = "-")
  data_list_adonis <- bind_rows(data_list_adonis,tmp_d[1,])
}
ggplot(data_list_adonis)+
  geom_bar(aes(x = class, y = `F`,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = class, y = `F`+0.2,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = color2)+
  theme(legend.position = "none")
ggsave(paste(out_dir,"pairpattern_adonis.pdf",sep = '/'),width = 6,height = 4)

all(row.names(data_pca) == phylo_used@sam_data$ID)
res_pca <- prcomp(data_pca, scale = T)
data_dist = as.matrix(vegdist(data_pca, method = 'bray'))
res_adonis_all = adonis2(data_dist~Pattern_class + age + sex + bmi + 
                           parasite_load_log + ST1_abs_log + ST2_abs_log + ST3_abs_log +
                           ST1 + ST2 + ST3 ,meta_used,by = "term",permutations = 999)
fviz_pca_ind(res_pca, label="none",axes = c(1, 2), alpha.ind =1,
             habillage=tmp_group,invisible="quali",pointsize = 2.5, addEllipses = TRUE,
             ellipse.level=0.95,palette = color0)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.1)
  )

ggsave(paste(out_dir,"beta_pca.pdf",sep = '/'),width = 6,height = 5)

tmp_apcoa_res <- aPCoA(data_dist ~ age + sex + bmi,
                       meta_used,maincov = Pattern_class)
tmp_apcoa_matrix <- data.frame(tmp_apcoa_res$plotMatrix)
identical(row.names(tmp_apcoa_matrix),meta_used$ID)
tmp_apcoa_matrix$group <- meta_used$Pattern_class
tmp_apcoa_matrix$sample = row.names(tmp_apcoa_matrix)
tmp_apcoa_matrix$parasite_log <- meta_used$parasite_load_log
tmp_list <- combine_list(tmp_apcoa_matrix$group)
a <- ggplot()+
  geom_point(data = tmp_apcoa_matrix,
             mapping = aes(x=X1,y=X2,color=group,shape=group),
             size=2.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.15,0.15),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PCoA1"),y=c("PCoA2"))+
  stat_ellipse(data=tmp_apcoa_matrix,geom = "polygon",
               aes(x=X1,y=X2,fill=group),
               color =NA,alpha=0.1,level = 0.9)+
  scale_fill_manual(values = color0)+
  scale_color_manual(values = color0)+
  scale_x_continuous(limits = c(-0.55,0.55))+
  scale_y_continuous(limits = c(-0.55,0.55))

b <- ggplot(data=tmp_apcoa_matrix,aes(x=group,y=X1,fill=group))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_boxplot(width = 0.8,alpha = 0.7)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  coord_flip()+
  scale_fill_manual(values = color0)+
  scale_y_continuous(limits = c(-0.55,0.55))
c <- ggplot(data=tmp_apcoa_matrix,aes(x=group,y=X2,fill=group))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_boxplot(width = 0.8,alpha = 0.7)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  scale_fill_manual(values = color0)+
  scale_y_continuous(limits = c(-0.55,0.55))

res_adonis_all$item <- row.names(res_adonis_all)
tmp_adonis_result_dis_p <- res_adonis_all[!is.na(res_adonis_all$F),]
tmp_adonis_result_dis_p$item <- factor(tmp_adonis_result_dis_p$item,levels = tmp_adonis_result_dis_p$item)
tmp_adonis_result_dis_p$class <- ifelse(tmp_adonis_result_dis_p$`Pr(>F)`>0.05,"no_sig","sig")
p1 <- ggplot(tmp_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = `F`,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = `F`+ 0.2,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p2 <- ggplot(tmp_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = R2,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2+0.002,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
ggsave(plot = p2,filename = paste(out_dir,"beta_adonis_R2.pdf",sep = '/'),width = 10,height = 5)
tmp_adonis_result_dis_p <- as.data.frame(t(tmp_adonis_result_dis_p))
tmp_adonis_result_dis_p <- tmp_adonis_result_dis_p[3:5,]
tmp_adonis_result_dis_p$posy <- c(0.25,0.5,0.75)
tmp_adonis_result_dis_p$posx <- 1
tmp_adonis_result_dis_p$labels <- paste(row.names(tmp_adonis_result_dis_p),round(as.numeric(tmp_adonis_result_dis_p[,1]),3),sep = ':')

d <- ggplot()+
  geom_bar(data = as.data.frame(matrix(c(1),nrow = 1)),aes(x = V1,y = V1),stat ="identity",fill="white")+
  geom_text(data = tmp_adonis_result_dis_p,mapping = aes(x = posx,y = posy,label = labels))+
  theme_bw()+theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
#plot_grid(b,d,a,c,align = "hv",rel_widths = c(1,0.225),rel_heights = c(0.225,1))
aplot::plot_list(b,d,a,c,widths = c(1,0.225),heights = c(0.225,1))
ggsave(paste(out_dir,"beta_pcoa_cov.pdf",sep = '/'),width = 7,height = 7)

tmp_adonis_result_dis_p <- res_adonis_all[!is.na(res_adonis_all$F),]
tmp_adonis_result_dis_p$item <- factor(tmp_adonis_result_dis_p$item,levels = tmp_adonis_result_dis_p$item)
tmp_adonis_result_dis_p$class <- ifelse(tmp_adonis_result_dis_p$`Pr(>F)`>0.05,"no_sig","sig")

tmp_adonis_result_dis_p_2 <- tmp_adonis_result_dis_p %>% filter(item %in% c("parasite_load_log","ST1","ST2","ST3"))
ggplot(tmp_adonis_result_dis_p_2)+
  geom_bar(aes(x = item, y = R2,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2+0.002,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
ggsave(filename = paste(out_dir,"beta_adonis_R2_2.pdf",sep = '/'),width = 5,height = 4)

tmp_formula <- y ~ x
ggplot(tmp_apcoa_matrix) + 
  stat_smooth(aes(x=X1,y=parasite_log),method='lm',formula = tmp_formula,color = "navy",alpha = 0.35)+
  geom_point(aes(x=X1,y=parasite_log,color = group),size=3,alpha =0.7)+
  stat_cor(aes(x=X1,y=parasite_log),cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.35)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid = element_blank(),
        legend.position = "none"
  )+
  scale_color_manual(values = color0)+
  labs(y = "Log2(Parasite load + 1)", x = "PCoA1")
ggsave(paste(out_dir,"cor.Pco1_to_paraload.pdf",sep = '/'),width = 4,height = 4)

masslin_data <- as.data.frame(t(tmp_apcoa_matrix %>% dplyr::select(X1)))
Maaslin2(
  input_data = masslin_data, 
  input_metadata = meta_used, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM" ,
  output = paste(out_dir,"PCoA_cov",sep = '/'), 
  fixed_effects = c("sex","age","bmi","parasite_load_log"),
  reference = c("sex,male")
)

## variance partitioning
meta_used2 <- meta_used %>% filter(Pattern_class %in% pattern_type[c(1,3)])
meta_used2$Pattern_class <- factor(meta_used2$Pattern_class,levels = pattern_type[c(1,3)])
tmp_dist <- as.data.frame(data_dist)
tmp_dist <- as.data.frame(data_dist)[meta_used2$ID,meta_used2$ID]
tmp_load_matrix <- meta_used2["parasite_load_log"]
tmp_ST_matrix <- meta_used2[c("ST1","ST3")]
tmp_ST_pattern_matrix <- meta_used2["Pattern_class"]
identical(row.names(tmp_dist),row.names(tmp_load_matrix))
identical(row.names(tmp_ST_matrix),row.names(tmp_load_matrix))

res_vp <- varpart(tmp_dist, tmp_ST_pattern_matrix , tmp_load_matrix)
print(res_vp)
plot(res_vp)

mod_ST_1 <- dbrda(tmp_dist ~ Pattern_class, data = meta_used2)
print(anova(mod_ST_1))
mod_load_1 <- dbrda(tmp_dist ~ parasite_load_log, data = meta_used2)
anova(mod_load_1)

mod_ST <- dbrda(tmp_dist ~ Pattern_class + Condition(parasite_load_log), data = meta_used2)
anova(mod_ST)
mod_load <- dbrda(tmp_dist ~ parasite_load_log + Condition(Pattern_class), data = meta_used2)
anova(mod_load)

r2_ST_total <- RsquareAdj(mod_ST_1)$adj.r.squared       # [a + c]
r2_load_total <- RsquareAdj(mod_load_1)$adj.r.squared   # [b + c]
r2_combined <- RsquareAdj(dbrda(tmp_dist ~ Pattern_class + parasite_load_log, data = meta_used2))$adj.r.squared  # [a + b + c]

## mediation
tmp_apcoa_matrix2 <- tmp_apcoa_matrix %>% filter(group %in% pattern_type[c(1,3)])
med_model <- lm(parasite_log ~ group, data = tmp_apcoa_matrix2)
out_model <- lm(X1 ~ group + parasite_log, data = tmp_apcoa_matrix2)

med_result <- mediate(med_model, out_model, treat = "group", mediator = "parasite_log",
                      boot = TRUE, sims = 1000)
summary(med_result)

## distance
data_list_pattern <- data.frame()
for(i in 1:length(tmp_list)){
  tmp_meta <- meta_used %>% filter(Pattern_class %in% tmp_list[[i]]) %>% dplyr::select(ID,Pattern_class)
  tmp_dist <- data_dist[row.names(data_dist) %in% tmp_meta$ID,]
  tmp_dist <- as.data.frame(tmp_dist)
  tmp_sam1 <- tmp_meta[tmp_meta$Pattern_class == tmp_list[[i]][1],1]
  tmp_sam2 <- tmp_meta[tmp_meta$Pattern_class == tmp_list[[i]][2],1]
  tmp_dist <- tmp_dist[tmp_sam1,tmp_sam2]
  tmp_dist$sample1 <- row.names(tmp_dist)
  tmp_dist_l <- gather(tmp_dist,key = "sample2",value = "value",-sample1)
  tmp_dist_l$class <- paste(tmp_list[[i]],collapse = "-")
  data_list_pattern <- bind_rows(data_list_pattern,tmp_dist_l)
}
tmp_list2 <- combine_list(data_list_pattern$class)
ggplot(data_list_pattern,aes(x=class,y=value))+
  geom_violin(aes(fill = class),color = "white")+
  geom_boxplot(width = 0.3)+
  stat_compare_means(method="wilcox.test",
                     comparisons = tmp_list2,
                     label = "p.format"
  )+
  labs(y = "Bray-Curtis distance")+
  theme_cowplot()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color2[c(2:4)])
ggsave(paste(out_dir,"pairpattern_distance_bray.pdf",sep = '/'),width = 6,height = 4)
i = 1
data_list_pattern$para_diff <- NA
for(i in 1:nrow(data_list_pattern)){
  sam1 <- data_list_pattern[i,1]
  sam2 <- data_list_pattern[i,2]
  para1 <- meta_used[meta_used$ID == sam1,"parasite_load_log"]
  para2 <- meta_used[meta_used$ID == sam2,"parasite_load_log"]
  tmp_diff <- abs(para1-para2)
  data_list_pattern[i,5] <- tmp_diff
}

formula_tmp = y~x 
ggplot(data_list_pattern) + 
  stat_smooth(aes(x=value,y=para_diff),method='lm',formula = tmp_formula,color = "navy",alpha = 0.35)+
  ggrastr::geom_point_rast(aes(x=value,y=para_diff),color="grey80",size=2,alpha =0.3)+
  stat_cor(aes(x=value,y=para_diff),cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.35)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid = element_blank(),
        legend.position = "none"
  )+
  scale_color_manual(values = color0)
ggsave(paste(out_dir,"cor.bc_to_paraload_diff.pdf",sep = '/'),width = 4,height = 4)

out_dir <- c("output/pattern_sampling")
phylo_used3 <- subset_samples(phyloseqin, species == sp[3])
phylo_used3 <- subset_samples(phylo_used,Pattern_class %in% pattern_type[c(1,3)])
meta_used3 <- metadatadf %>% filter(species == sp[3],Pattern_class %in% pattern_type[c(1,3)])
meta_used3$Pattern_class <- factor(meta_used3$Pattern_class,levels = pattern_type[c(1,3)])
phylo_used3 <- filter_taxa(phylo_used3, function(x) sum(x > 0) > 0, TRUE)

ggplot(meta_used3)+
  geom_histogram(aes(x = parasite_load_log,fill = Pattern_class
                     ,color = Pattern_class),alpha = 0.3,boundary = 0,binwidth = 0.5)+
  scale_color_manual(values = color0[c(1,3)])+
  scale_fill_manual(values = color0[c(1,3)])+
  scale_x_continuous(breaks = seq(0,15,1),limits = c(0,15))+
  theme_cowplot()+
  facet_wrap(Pattern_class~.,nrow = 2)

ggsave(paste(out_dir,"ST1_ST3.paraload.log_distribution.pdf",sep = '/'),width = 6,height = 6)

n_successful = 0
n_time = 1000
n_size = 15
d_res_sampling <- matrix(nrow = n_time,ncol = 3)
while (n_successful < n_time) {
  d_tmp_1 <- meta_used %>% filter(Pattern_class == pattern_type[c(1)]) %>% 
    filter(parasite_load_log > 2 & parasite_load_log < 9)
  d_tmp_2 <- meta_used %>% filter(Pattern_class == pattern_type[c(3)]) %>%
    filter(parasite_load_log > 2 & parasite_load_log < 9)
  d_tmp_1 <- as.data.frame(d_tmp_1)
  d_tmp_2 <- as.data.frame(d_tmp_2)
  
  s_1 = sample(d_tmp_1$ID,size = n_size)
  s_2 = sample(d_tmp_2$ID,size = n_size)
  
  d_tmp_1 <- d_tmp_1[d_tmp_1$ID %in% s_1,]
  d_tmp_2 <- d_tmp_2[d_tmp_2$ID %in% s_2,]
  
  t_tmp <- wilcox.test(d_tmp_1$parasite_load_log,d_tmp_2$parasite_load_log)
  if(t_tmp$p.value > 0.05){
    phylo_tmp <- subset_samples(phylo_used, ID %in% c(s_1,s_2))
    phylo_tmp <- filter_taxa(phylo_tmp, function(x) sum(x > 0) > 0, TRUE)
    d_otu <- as.data.frame(t(phylo_tmp@otu_table))
    d_meta_tmp <- rbind(d_tmp_1,d_tmp_2)
    dist_tmp <- as.matrix(vegdist(d_otu, method = 'bray'))
    r_adonis = adonis2(dist_tmp ~ Pattern_class + sex + age + bmi, d_meta_tmp,by = "term")
    n_successful = n_successful + 1
    print(paste(n_successful,r_adonis$`Pr(>F)`[1],sep = ":"))
    d_res_sampling[n_successful,] <- c(t_tmp$p.value,r_adonis$`Pr(>F)`[1],mean(d_meta_tmp$parasite_load_log))
  }
  #wilcox.test(d_tmp_1$age,d_tmp_2$age)
  #wilcox.test(d_tmp_1$bmi,d_tmp_2$bmi)
  #chisq.test(matrix(c(5,7,10,8),nrow = 2))
}
d_res_sampling <- as.data.frame(d_res_sampling)
colnames(d_res_sampling) <- c("wilcox","adonis","parasite_load_log")
d_res_sampling$tmp <- "n"
length(which(d_res_sampling[,2] < 0.05))

ggplot(d_res_sampling,aes(x = tmp,y = adonis))+
  geom_hline(yintercept = 0.05,linetype = "dashed")+
  geom_violin(fill = "grey80",color = "white")+
  geom_boxplot(width=0.2,fill = "white")+
  theme_cowplot()+
  scale_y_continuous(
    breaks = c(0,0.05,0.25,0.5,0.75,1),limits = c(0,1)
  )+
  labs(x = "", y = "P-value of adonis")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
ggsave(paste(out_dir,"ST1_ST3.paraload.adonis.pdf",sep = '/'),width = 4,height = 6)

meta_used4 <- meta_used3 %>% filter(parasite_load_log < 8 & parasite_load_log > 2)
phylo_used4 <- subset_samples(phylo_used3,ID %in% meta_used4$ID)
phylo_used4 <- filter_taxa(phylo_used4, function(x) sum(x > 0) > 0, TRUE)
d_otu <- as.data.frame(t(phylo_used4@otu_table))
d_dist = as.matrix(vegdist(d_otu, method = 'bray'))
r_adonis_all = adonis2(d_dist~Pattern_class + age + sex + bmi,
                       meta_used4,by = "term")

r_apcoa <- aPCoA(d_dist ~ age + sex + bmi + parasite_load_log,
                 meta_used4,maincov = Pattern_class)
d_apcoa_matrix <- data.frame(r_apcoa$plotMatrix)
identical(row.names(d_apcoa_matrix),meta_used4$ID)
d_apcoa_matrix$group <- meta_used4$Pattern_class
d_apcoa_matrix$sample = row.names(d_apcoa_matrix)
l_tmp <- combine_list(d_apcoa_matrix$group)

a <- ggplot()+
  geom_point(data = d_apcoa_matrix,
             mapping = aes(x=X1,y=X2,color=group,shape=group),
             size=2.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.15,0.15),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PCoA1"),y=c("PCoA2"))+
  stat_ellipse(data=d_apcoa_matrix,geom = "polygon",
               aes(x=X1,y=X2,fill=group),
               color =NA,alpha=0.1,level = 0.9)+
  scale_fill_manual(values = color0[c(1,3)])+
  scale_color_manual(values = color0[c(1,3)])+
  scale_x_continuous(limits = c(-0.45,0.6))+
  scale_y_continuous(limits = c(-0.4,0.35))
b <- ggplot(data=d_apcoa_matrix,aes(x=group,y=X1,fill=group))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_boxplot(width = 0.8,alpha = 0.7)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  coord_flip()+
  scale_fill_manual(values = color0[c(1,3)])+
  scale_y_continuous(limits = c(-0.45,0.6))
c <- ggplot(data=d_apcoa_matrix,aes(x=group,y=X2,fill=group))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_boxplot(width = 0.8,alpha = 0.7)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  scale_fill_manual(values = color0[c(1,3)])+
  scale_y_continuous(limits = c(-0.4,0.35))

r_adonis_all$item <- row.names(r_adonis_all)
d_adonis_result_dis_p <- r_adonis_all[!is.na(r_adonis_all$F),]
d_adonis_result_dis_p$item <- factor(d_adonis_result_dis_p$item,levels = d_adonis_result_dis_p$item)
d_adonis_result_dis_p$class <- ifelse(d_adonis_result_dis_p$`Pr(>F)`>0.05,"no_sig","sig")
p1 <- ggplot(d_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = `F`,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = `F`+ 0.2,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
p2 <- ggplot(d_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = R2,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2+0.002,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

ggsave(plot = p2,filename = paste(out_dir,"beta_adonis_R2.pdf",sep = '/'),width = 10,height = 5)
d_adonis_result_dis_p <- as.data.frame(t(d_adonis_result_dis_p))
d_adonis_result_dis_p <- d_adonis_result_dis_p[3:5,]
d_adonis_result_dis_p$posy <- c(0.25,0.5,0.75)
d_adonis_result_dis_p$posx <- 1
d_adonis_result_dis_p$labels <- paste(row.names(d_adonis_result_dis_p),round(as.numeric(d_adonis_result_dis_p[,1]),3),sep = ':')

d <- ggplot()+
  geom_bar(data = as.data.frame(matrix(c(1),nrow = 1)),aes(x = V1,y = V1),stat ="identity",fill="white")+
  geom_text(data = d_adonis_result_dis_p,mapping = aes(x = posx,y = posy,label = labels))+
  theme_bw()+theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
aplot::plot_list(b,d,a,c,widths = c(1,0.225),heights = c(0.225,1))
ggsave(paste(out_dir,"ST1_ST3.beta_pcoa_cov.pdf",sep = '/'),width = 7,height = 7)

meta_used5 <- meta_used4 %>% select(ID,Pattern_class)
c_sam3 <- meta_used5[meta_used5$Pattern_class == pattern_type[1],1]
c_sam4 <- meta_used5[meta_used5$Pattern_class == pattern_type[3],1]
d_distn2 <- as.data.frame(as.matrix(data_dist[c_sam3,c_sam4]))
d_distn2$ID <- row.names(d_distn2)
d_dist_l2 <- d_distn2 %>% gather(key = "ID2",value = "value",-ID)
d_dist_l2$class <- "same"
colnames(d_dist_l2) <- c("sample1","sample2","value","class")
d_dist_l2_c <- rbind(d_dist_l2,as.data.frame(data_list_pattern%>% filter(class == "ST1_dom-ST3_dom")))

l_tmp3 <- combine_list(d_dist_l2_c$class)
ggplot(d_dist_l2_c,aes(x=class,y=value))+
  geom_violin(aes(fill = class),color = "white")+
  geom_boxplot(width = 0.3)+
  stat_compare_means(method="wilcox.test",
                     comparisons = l_tmp3,
                     label = "p.format"
  )+
  labs(y = "Bray-Curtis distance")+
  theme_cowplot()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = color2[c(4,5)])
ggsave(paste(out_dir,"ST1_ST3.all_vs_no_paraload.pdf",sep = '/'),width = 3.5,height = 5)
ggstatsplot::ggbetweenstats(
  data = d_dist_l2_c,
  x = class,
  y = value,
  type = "parametric",
  )

rm(list=ls(pattern = "tmp"))
rm(list=ls(pattern = "dist"))
rm(list=ls(pattern = "ann"))
rm(list=ls(pattern = "phyloseqin"))


##useless
phylo_used2 <- subset_samples(phyloseqin,class %in% "control")
phylo_used2 <- subset_samples(phylo_used2,species == sp[3])
meta_used2 <- as.data.frame(as.matrix(sample_data(phylo_used2)))
meta_used2$Pattern_class <- "Control"
meta_used2 <- rbind(meta_used,meta_used2)
phylo_used2 <- subset_samples(phyloseqin,ID %in% meta_used2$ID)
phylo_used2 <- filter_taxa(phylo_used2, function(x) sum(x > 0) > 0, TRUE)
d_otu <- as.data.frame(t(phylo_used2@otu_table))
d_dist = as.matrix(vegdist(d_otu, method = 'bray'))
c_sam1 <- meta_used2[meta_used2$Pattern_class == "Control",1]
c_sam2 <- meta_used2[meta_used2$Pattern_class != "Control",1]
d_dist <- as.data.frame(d_dist[c_sam2,c_sam1])
d_dist_m<- as.data.frame(as.matrix(apply(d_dist,1,mean)))
d_dist_m$ID = row.names(d_dist_m)
colnames(d_dist_m) <- c("mean","ID")
#d_dist_m <- d_dist %>% select(ID,mean)
d_dist_m <- left_join(d_dist_m,meta_used2 %>% dplyr::select(ID,parasite_load_log,Pattern_class),by = "ID")
d_dist_m$mean <- as.numeric(d_dist_m$mean)
d_dist_m$parasite_load_log <- as.numeric(d_dist_m$parasite_load_log)
formula_tmp = y~x 
ggplot(d_dist_m,aes(x=mean,y=parasite_load_log)) + 
  stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
  geom_point(color = "gray40",size=3)+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid = element_blank(),
  )+
  facet_wrap(~Pattern_class,scales = "free")

d_distn = as.data.frame(as.matrix(vegdist(d_otu, method = 'bray')))
d_distn$ID <- rownames(d_distn)
d_distnl <- left_join(d_distn,meta_used2 %>% dplyr::select(ID),by = "ID") %>%
  gather(key = "ID2",value = "value",rownames(d_distn))
d_distnl <- d_distnl %>% filter(ID != ID2)
d_distnl$class <- NA
for(i in 1:nrow(d_distnl)){
  t_sam <- c(d_distnl[i,1],d_distnl[i,2])
  t_meta <- meta_used2 %>% filter(ID %in% t_sam) %>% dplyr::select (Pattern_class)
  t_t <- paste(t_meta$Pattern_class,collapse = "-")
  d_distnl[i,4] <- t_t
}
d_distnl2 <- d_distnl[grep("Control$",d_distnl$class),]
l_tmp <- combine_list(d_distnl2$class)

ggplot(d_distnl2,aes(x=class,y=value))+
  geom_violin(aes(fill = class),color = "white")+
  geom_boxplot(width = 0.3)+
  labs(y = "Bray-Curtis distance")+
  theme_cowplot()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())