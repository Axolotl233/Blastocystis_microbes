plot_sp_box <- function(otu,meta,sp,pre_cut = 0.2){
  plot_r <- list()
  data_r <- list() 
  j = 1
  for(i in 1:length(sp)){
    if(! sp[i] %in% row.names(otu)){
      print(paste(c("sp is not exist:",sp[i]),collapse = " "))
      next
    }
    tmp_data <- as.data.frame(t(otu[which(row.names(otu)==sp[i]),meta$sample]))
    #print(sp[i])
    tmp_data$sample <- row.names(tmp_data)
    colnames(tmp_data) <- c("value","sample")
    tmp_data_p <- inner_join(meta,tmp_data,by = "sample")
    if(length(which(tmp_data_p$value != 0)) < nrow(tmp_data_p) * pre_cut ){
      print(paste(c("sp have lower prevalence than threshold:",sp[i]),collapse = " "))
      next
    }
    data_r[[j]] <- tmp_data_p
    tmp_c=unique(tmp_data_p$groups)
    tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
      mutate(tmp = str_c(V1,",",V2)) %>% 
      select(tmp)
    tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
    plot_r[[j]] <- ggplot(tmp_data_p,aes(x=groups,y=value,color = groups))+
      geom_boxplot(outlier.shape = NA,)+
      geom_jitter(width = 0.2,height = 0)+
      theme_cowplot()+
      theme(
        legend.position = "none"
      )+
      labs(x="Group",y="",title = sp[i])+
      stat_compare_means(comparisons = tmp_list,
                         method = "wilcox.test", label = "p.format")+
      scale_color_manual(values = color0)
    j = j + 1
  }
  res_r <-list(plot_r,data_r) 
  return(res_r)
}
pathway_beta <- function(x,y,z){
  tmp_meta_pca <- y %>% dplyr::select(sample,groups,age,sex)
  out1 <- paste("output/cross",z,"pca.pdf",sep = ".")
  out2 <- paste("output/cross",z,"pcoa.pdf",sep = ".")
  data_path_pca=x
  data_path_pca <- noise_removal(data_path_pca,percent = 0.1,low = 0.001)
  data_path_pca <- as.data.frame(t(data_path_pca))
  data_path_pca$sample <- row.names(data_path_pca)
  data_path_pca <- data_path_pca %>% filter(sample %in% tmp_sam)
  data_path_pca <- inner_join(data_path_pca,tmp_meta_pca,by = "sample")
  row.names(data_path_pca) <- data_path_pca$sample
  
  tmp_meta_pca <- data_path_pca %>% dplyr::select(sample,groups,age,sex)
  data_path_pca <- data_path_pca %>% dplyr::select(!all_of(c("sample","groups","age","sex")))
  
  tmp_dist = vegdist(as.matrix(data_path_pca), method = 'bray')
  tmp_site = as.data.frame(data.frame(sample = tmp_meta_pca$sample,
                                      group = tmp_meta_pca$groups,
                                      age = tmp_meta_pca$age,
                                      sex = tmp_meta_pca$sex))
  tmp_adonis_result_dis = adonis2(data_path_pca~group, tmp_site)
  res_pca <- prcomp(data_path_pca, scale = T)
  
  p1 <- fviz_pca_ind(res_pca, label="none",axes = c(1, 2), alpha.ind =1,
                     habillage=tmp_meta_pca$groups,invisible="quali",pointsize = 2.5, addEllipses = TRUE,
                     ellipse.level=0.66,palette = color0)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(vjust = 0.5, hjust = 0.1)
    )+
    annotate("text", x=max(res_pca$x[,1])-25, y=max(res_pca$x[,2])-5,
             label=paste("p-value(adonis)=",tmp_adonis_result_dis$`Pr(>F)`[1]))+
    annotate("text", x=max(res_pca$x[,1])-25, y=max(res_pca$x[,2])-12, 
             label=paste("R2(adonis)=",round(tmp_adonis_result_dis$R2[1],3)))+
    labs(title = paste(z,"Indiviuals pathway PCA",sep = " "))
  ggsave(out1,p1,width = 6,height = 4)
  
  row.names(tmp_site) <- tmp_site$sample
  tmp_apcoa_res <- aPCoA(tmp_dist ~ age,tmp_site,maincov = group)
  tmp_apcoa_matrix <- data.frame(tmp_apcoa_res$plotMatrix)
  identical(row.names(tmp_apcoa_matrix),tmp_site$sample)
  tmp_apcoa_matrix$group <- tmp_site$group
  p2 <- ggplot(data=tmp_apcoa_matrix,aes(x=X1,y=X2,
                                         color=group,shape=group))+
    geom_point(size=2.5)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.title = element_blank())+
    geom_vline(xintercept = 0,lty="dashed")+
    geom_hline(yintercept = 0,lty="dashed")+
    labs(x=c("PCoA1"),
         y=c("PCoA2"),
         title = paste(z,"Indiviuals pathway PCoA",sep = " ")
    )+
    stat_ellipse(data=tmp_apcoa_matrix,
                 geom = "polygon",
                 aes(fill=group),
                 alpha=0.1)+
    scale_fill_manual(values = color0)+
    scale_color_manual(values = color0)
  ggsave(out2,p2,width = 6,height = 4)
}

library(factoextra)
library(vegan)
library(microeco)
library(file2meco)
library(magrittr)
library(aplot)
library(compositions)
library(Maaslin2)
library(maaslin3)

out_dir <- paste("output/pathway/",sp[3],sep = "")
dir.create(out_dir,showWarnings = F,recursive = T)
load("workfile/z.phyloseq.Rdata")
data_pathway <- read_csv("input/humann.pathabundance.ab.csv")

meta_used <- metadatadf %>% filter(species == sp[3],Pattern_class %in% pattern_type[c(1:3)])

tmp_data_pathway <- data_pathway[grep("\\|",data_pathway$Pathway,invert = T),]
tmp_data_pathway <- as.data.frame(tmp_data_pathway)
row.names(tmp_data_pathway) <- tmp_data_pathway$Pathway
tmp_data_pathway <- tmp_data_pathway %>% select (! Pathway) %>% select(meta_used$ID)
tmp_data_pathway <- noise_removal(tmp_data_pathway,method = "max_cut",low = 0)

tmp_data_pathway_s <- data_pathway[grep("s__",data_pathway$Pathway),]
tmp_data_pathway_s <- as.data.frame(tmp_data_pathway_s)
row.names(tmp_data_pathway_s) <- tmp_data_pathway_s$Pathway
tmp_data_pathway_s <- tmp_data_pathway_s %>% select (! Pathway)
tmp_data_pathway_s <- tmp_data_pathway_s %>% select(meta_used$ID)
tmp_data_pathway_s <- noise_removal(tmp_data_pathway_s,method = "max_cut",low = 0)

tmp_data_run1 <- as.data.frame(t(tmp_data_pathway))
tmp_data_run2 <- as.data.frame(t(tmp_data_pathway_s))

tmp_list <- combine_list(meta_used$Pattern_class)
for(i in 1:length(tmp_list)){
  tmp_list[[i]] = sort(tmp_list[[i]])
  tmp_meta <- meta_used %>% filter(Pattern_class %in% tmp_list[[i]])
  tmp_d1 <- tmp_data_run1[row.names(tmp_data_run1) %in% tmp_meta$ID,]
  tmp_d2 <- tmp_data_run2[row.names(tmp_data_run2) %in% tmp_meta$ID,]
  tmp_fn <- paste(tmp_list[[i]][1],tmp_list[[i]][2],sep = "_vs_")
  tmp_ref <- paste("Pattern_class",tmp_list[[i]][1],sep = ",")
  maaslin3(
    input_data = tmp_d1,
    input_metadata = tmp_meta,
    min_prevalence = 0,
    min_abundance = 0,
    
    formula = '~ Pattern_class + age + sex + bmi',
    normalization = "TSS",
    transform = 'LOG',
    
    output = paste(out_dir,paste(tmp_fn,"maaslin3_all",sep = "_"),sep = '/'),
    augment = TRUE,
    standardize = F,
    max_significance = 0.1,
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,
    max_pngs = 250,
    cores = 1
  )
  Maaslin2(
    input_data = tmp_d1, 
    input_metadata = tmp_meta, 
    min_prevalence = 0,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    output = paste(out_dir,paste(tmp_fn,"maaslin2_all",sep = "_"),sep = '/'),
    fixed_effects = c("Pattern_class","age","sex","bmi"),
    reference = tmp_ref,
  )
  maaslin3(
    input_data = tmp_d2,
    input_metadata = tmp_meta,
    min_prevalence = 0,
    min_abundance = 0,
    
    formula = '~ Pattern_class + age + sex + bmi',
    normalization = "TSS",
    transform = 'LOG',
    
    output = paste(out_dir,paste(tmp_fn,"maaslin3_sp",sep = "_"),sep = '/'),
    augment = TRUE,
    standardize = F,
    max_significance = 0.1,
    median_comparison_abundance = TRUE,
    median_comparison_prevalence = FALSE,
    max_pngs = 250,
    cores = 1
  )
  Maaslin2(
    input_data = tmp_d2, 
    input_metadata = tmp_meta, 
    min_prevalence = 0,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    output = paste(out_dir,paste(tmp_fn,"maaslin2_sp",sep = "_"),sep = '/'),
    fixed_effects = c("Pattern_class","age","sex","bmi"),
    reference = tmp_ref,
  )
}

Maaslin2(
  input_data = tmp_data_run1, 
  input_metadata = meta_used, 
  min_prevalence = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = paste(out_dir,"all",sep = '/'), 
  fixed_effects = c("ST1_abs_log","ST2_abs_log","ST3_abs_log")
)

maaslin3(
  input_data = tmp_data_run1,
  input_metadata = meta_used,
  min_prevalence = 0,
  min_abundance = 0,
  
  formula = '~ ST1_abs_log + ST2_abs_log + ST3_abs_log +
  ST1_abs_log:ST2_abs_log +  ST1_abs_log:ST3_abs_log +  ST3_abs_log:ST2_abs_log +
  ST1_abs_log:ST2_abs_log:ST3_abs_log',
  normalization = "TSS",
  transform = 'LOG',
  
  output = paste(out_dir,"all_3",sep = '/'),
  augment = TRUE,
  standardize = F,
  median_comparison_abundance = TRUE,
  median_comparison_prevalence = FALSE,
  max_pngs = 250,
  cores = 1
)

for(i in 1:length(pattern_type[c(1:3)])){
  tmp_meta <- meta_used %>% filter(Pattern_class %in% pattern_type[i])
  tmp_d3 <- tmp_data_run1[row.names(tmp_data_run1) %in% tmp_meta$ID,]
  #tmp_fn <- paste(tmp_list[[i]][1],tmp_list[[i]][2],sep = "_vs_")
  #tmp_ref <- paste("Pattern_class",tmp_list[[i]][1],sep = ",")
  Maaslin2(
    input_data = tmp_d3, 
    input_metadata = tmp_meta, 
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    output = paste(out_dir,pattern_type[i],sep = '/'), 
    fixed_effects = c("ST1_abs_log","ST2_abs_log","ST3_abs_log")
  )
}

tmp_maaslin <- read_table(paste(out_dir,"ST1_dom_vs_ST3_dom_maaslin2_all","significant_results.tsv",sep = '/'))
tmp_maaslin <- tmp_maaslin %>% filter(metadata == 'Pattern_class') %>% arrange( desc(coef))
tmp_maaslin$feature <- gsub("\\.\\."," ",tmp_maaslin$feature)
tmp_maaslin$feature <- gsub("\\."," ",tmp_maaslin$feature)
tmp_maaslin$feature <- factor(tmp_maaslin$feature,levels = tmp_maaslin$feature)
tmp_maaslin$log_pval <- -log10(tmp_maaslin$pval)

ggplot(tmp_maaslin)+
  geom_bar(aes(x = feature, y = coef, fill = qval),stat = "identity")+
  geom_errorbar(aes(x = feature, ymin = coef - stderr,ymax = coef + stderr),width = 0.2,stat = "identity")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  scale_fill_gradient2(high = "white",low = "#8856a7",mid = "#9ebcda",midpoint = 0.16)+
  coord_flip()+
  labs(title = "ST1_dom_vs_ST3_dom")
ggsave(paste(out_dir,"ST1_dom_vs_ST3_dom_maaslin2_all.pathway_diff.pdf",sep = '/'),width = 8,height = 5)

meco_pathway_metacyc <- humann2meco("input/humann.pathabundance.ab.tsv", 
                                    db = "MetaCyc", 
                                    sample_table = meta_used 
)
meco_pathway_metacyc$tidy_dataset()
meco_pathway_metacyc$cal_abund(rel = T)
meco_pathway_metacyc$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Superclass2 %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$pathway %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Species %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Genus %<>% .[!grepl("unclass", rownames(.)), ]

meco_pathway_metacyc_used <- trans_abund$new(meco_pathway_metacyc, taxrank = "Superclass2", 
                                             ntaxa = 10, use_percentage = FALSE)
meco_pathway_metacyc_used$ylabname <- "Relative Abundace"
meco_pathway_metacyc_used$plot_bar(facet = "Pattern_class", bar_type = "notfull")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

ggsave(paste(out_dir,"pathway_metacyc_bar.all.pdf",sep = '/'),width = 20,height = 6)
tmp_plot_data <- meco_pathway_metacyc_used$data_abund %>% arrange(desc(all_mean_abund))
tmp_plot_data$Taxonomy <- factor(tmp_plot_data$Taxonomy,levels = rev(unique(tmp_plot_data$Taxonomy)))
tmp_plot_pathway <- unique(tmp_plot_data$Taxonomy)[1:15]

ggplot(tmp_plot_data %>% filter(Taxonomy %in% tmp_plot_pathway),
       aes(x=Pattern_class,y=Abundance,fill = Taxonomy)) +
  stat_summary(fun=median, geom="bar" ,width = 0.6,position = "stack")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = rev(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                   "#A6761D","#666666","#8DCEBB","#ECAF80","#BAB7D9","#F394C4",
                   "#B2D28E","#F2D480","#D2BA8E"))
  )+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
ggsave((paste(out_dir,"pathway_metacyc_bar.mean.pdf",sep = '/')),width = 8,height = 6)

tmp_data_stat <- data.frame()
for(i in 1:length(tmp_plot_pathway)){
  tmp_data <- tmp_plot_data %>% filter(Taxonomy == tmp_plot_pathway[i])
  tmp <- summary_stat(tmp_data,3,13) %>% mutate(Taxaonomy = tmp_plot_pathway[i]) 
  tmp_data_stat <- rbind(tmp_data_stat,tmp)
}


meco_pathway_metacyc$cal_betadiv(method = "bray")
meco_pathway_metacyc_distance <- as.dist(meco_pathway_metacyc$beta_diversity$bray)
meco_pathway_metacyc_beta <- trans_beta$new(dataset = meco_pathway_metacyc, 
                                            group = "Pattern_class", measure = "bray")

meco_pathway_metacyc_beta$cal_ordination(method = "PCoA")

meco_pathway_metacyc_beta_tmp <- meco_pathway_metacyc_beta$res_ordination$scores
meco_pathway_metacyc_beta_tmp2 <- trans_env$new(dataset = meco_pathway_metacyc, 
                                                add_data = meco_pathway_metacyc_beta_tmp[, 1:2])
meco_pathway_metacyc_beta$plot_ordination(plot_color = "Pattern_class", 
                                                plot_type = c("point", "ellipse"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(
    values = color0[c(2,1,3)]
  )+
  scale_fill_manual(
    values = color0[c(2,1,3)]
  )

ggplot()+
  geom_point(data = meco_pathway_metacyc_beta_tmp,
             mapping = aes(x=PCo1,y=PCo2,color=Pattern_class,shape=Pattern_class),
             size=2.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        #legend.position = c(0.88,0.88),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  labs(x=c("PCoA1"),y=c("PCoA2"))+
  stat_ellipse(data=meco_pathway_metacyc_beta_tmp,geom = "polygon",
               aes(x=PCo1,y=PCo2,fill=Pattern_class),
               color =NA,alpha=0.1)+
  scale_fill_manual(values = color0[c(2,1,3)])+
  scale_color_manual(values = color0[c(2,1,3)])

ggsave(paste(out_dir,"pathway_metacyc.pcoa.pdf",sep = '/'),width = 7,height = 5)

res_adonis_list <- list()
for(i in 1:length(tmp_list)){
  tmp_meta2 = meta_used %>% filter(Pattern_class %in% tmp_list[[i]])
  tmp_dist2 <- as.matrix(meco_pathway_metacyc_distance)
  tmp_dist2 <- as.dist(tmp_dist2[tmp_meta2$ID,tmp_meta2$ID])
  tmp_adonis_result_dis2 = adonis2(tmp_dist2 ~ Pattern_class + age + sex + bmi + parasite_load_log, tmp_meta2)
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

