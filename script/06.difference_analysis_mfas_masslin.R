phase_ancom_add_tax <- function(tab,tax,target_class = c("Phylum","Order"),ref_class = "Species"){
  tab$taxon_id <- gsub("s__","",tab$taxon_id)
  tax_use = tax[tax$Species %in% tab$taxon_id,]
  tax_use = tax_use[,c(target_class,ref_class)]
  colnames(tax_use) = c(target_class,"taxon_id")
  res_tab <- merge(tab,tax_use,by.x = "taxon_id")
  return(res_tab)
}
phase_tax_file <- function(tax){
  tax[,c(1:6)] <- apply(tax[,c(1:6)], 2, function(x){gsub("Candidatus_","",x)})
  tax[,c(1:6)] <- apply(tax[,c(1:6)], 2, function(x){ifelse(grepl("unclassified",x),"Unclassified",x)})
  tax <- as.data.frame(tax)
  tax$Class <- ifelse(grepl("FGB",tax$Class),"Unclassified",tax$Class)
  tax$Order <- ifelse(grepl("FGB",tax$Order),"Unclassified",tax$Order)
  tax$Family <- ifelse(grepl("FGB",tax$Family),"Unclassified",tax$Family)
  tax$Genus <- ifelse(grepl("GGB",tax$Genus),"Unclassified",tax$Genus)
  return(tax)
}

library(Maaslin2)
library(ggsankey)

load("workfile/z.phyloseq.Rdata")

out_dir <- paste("output/different_analysis/",sp[3],sep = "")
dir.create(out_dir,showWarnings = F,recursive = T)

phylo_used <- subset_samples(phyloseqin, species == sp[3])
phylo_used <- subset_samples(phylo_used,Pattern_class %in% pattern_type[c(1:3)])
meta_used <- metadatadf %>% filter(species == sp[3],Pattern_class %in% pattern_type[c(1:3)])
meta_used$Pattern_class <- factor(meta_used$Pattern_class,levels = pattern_type[c(1:3)])
phylo_used <- filter_taxa(phylo_used, function(x) sum(x > 0) > 0, TRUE)
tmp_tax <- as.data.frame(phylo_used@tax_table)
tmp_tax <- phase_tax_file(tmp_tax)
d_otu <- as.data.frame(t(otu_table(phylo_used)))

Maaslin2(
  input_data = d_otu, 
  input_metadata = meta_used, 
  min_prevalence = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = paste(out_dir,"all",sep = '/'), 
  fixed_effects = c("parasite_load_log","sex","age","bmi")
)

Maaslin2(
  input_data = d_otu, 
  input_metadata = meta_used, 
  min_prevalence = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = paste(out_dir,"subtype",sep = '/'), 
  fixed_effects = c("ST1_abs_log","ST2_abs_log","ST3_abs_log","sex","age","bmi")
)

diff_res <- read_table(paste(out_dir,"all/significant_results.tsv",sep = '/'))
diff_res <- diff_res %>% filter(metadata == "parasite_load_log")
diff_res <- diff_res %>% filter(qval < 0.1)
tmp_tax$feature <- row.names(tmp_tax)
diff_res_tax <- left_join(diff_res,tmp_tax,by = "feature") %>% arrange(coef)
diff_res_tax$Species2 <- 100:(nrow(diff_res_tax)+99)

diff_res_tax_p <- diff_res_tax %>%
  make_long(Phylum,Class,Order,Family,Genus,Species2)

ggplot(diff_res_tax_p, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
   geom_alluvial(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "black", fill = "white") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 18) +
  labs(x = NULL)+
  theme(legend.position = "none")

ggplot(diff_res_tax_p, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6,
                node.color = "gray30") +
  geom_alluvial_text(size = 3, color = "black", fill = "white") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 18) +
  labs(x = NULL)+
  theme(legend.position = "none")
ggsave(paste(out_dir,"sankey_0.1.pdf",sep = '/'),width = 10,height = 10)
diff_res_tax$log_qval <- -log10(diff_res_tax$qval)
diff_res_tax$feature <- factor(diff_res_tax$feature,levels = diff_res_tax$feature)
diff_res_tax$Species <- factor(diff_res_tax$Species,levels = diff_res_tax$Species)
ggplot(diff_res_tax)+
  geom_bar(aes(x = Species, y = coef, fill = log_qval),stat = "identity")+
  geom_errorbar(aes(x = Species, ymin = coef - stderr,ymax = coef + stderr),width = 0.2,stat = "identity")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
  )+
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(5,"BuPu"))(50))+
  coord_flip()
ggsave(paste(out_dir,"coef_0.1.pdf",sep = '/'),width = 6,height = 10)
d_otu_tmp <- data.frame(ID=row.names(d_otu),sp=d_otu[,"s__Limosilactobacillus_reuteri"])
d_otu_tmp_p <- inner_join(d_otu_tmp,meta_used %>% select(ID,parasite_load_log),by = "ID")
tmp_formula <- y ~ x

ggplot(d_otu_tmp_p,aes(x=sp,y=parasite_load_log)) + 
  stat_smooth(method='lm',formula = tmp_formula,color = "navy",level = 0.95)+
  geom_point(color = "gray40",size=3)+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.35)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid = element_blank(),
  )+labs(title = "Limosilactobacillus reuteri",x = "Abundance (%)",y = "Parasite")
ggsave(paste(out_dir,"sp_Limosilactobacillus_reuteri.pdf",sep = '/'),width = 4,height = 4)

Maaslin2(
  input_data = d_otu, 
  input_metadata = meta_used, 
  min_prevalence = 0,
  normalization = "TSS",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = paste(out_dir,"Pattern_class",sep = '/'), 
  fixed_effects = c("Pattern_class","sex","age","bmi")
)

tmp_list <- combine_list(meta_used$Pattern_class)
pattern_res <- data.frame()
pattern_res_sig <- c()
for(i in 1:length(tmp_list)){
  tmp_list[[i]] = sort(tmp_list[[i]])
  tmp_meta <- meta_used %>% filter(Pattern_class %in% tmp_list[[i]])
  tmp_d <- d_otu[rownames(d_otu) %in% tmp_meta$ID,]
  tmp_fn <- paste(tmp_list[[i]][1],tmp_list[[i]][2],sep = "-")
  tmp_ref <- paste("Pattern_class",tmp_list[[i]][1],sep = ",")
  Maaslin2(
    input_data = tmp_d, 
    input_metadata = tmp_meta, 
    min_prevalence = 0,
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    output = paste(out_dir,tmp_fn,sep = '/'), 
    fixed_effects = c("Pattern_class","sex","age","bmi"),
    reference = tmp_ref,
  )
  tmp_res <-  read_tsv(paste(c(out_dir,tmp_fn,"all_results.tsv"),collapse = '/'))
  tmp_res <- tmp_res %>% 
    filter(value == tmp_list[[i]][2]) %>% 
    dplyr::select(feature,coef,qval) %>% mutate(group = tmp_fn)
  pattern_res <- bind_rows(pattern_res,tmp_res)
  tmp_sig <- tmp_res %>% filter(qval < 0.25)
  tmp_sig <- tmp_sig$feature
  pattern_res_sig <- c(pattern_res_sig,tmp_sig)
}
pattern_res_sig <- unique(pattern_res_sig)
pattern_res <- pattern_res %>% filter(feature %in% pattern_res_sig) 
pattern_res$sig <- ifelse(pattern_res$qval<0.25,"aquamarine3", "black")

tmp_lo = floor(min(pattern_res$coef))
tmp_up = ceiling(max(pattern_res$coef))
tmp_mid = (tmp_lo + tmp_up)/2

ggplot(pattern_res, aes(x = group, y = feature, fill = coef)) + 
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = tmp_mid, limit = c(tmp_lo,tmp_up),
                       name = NULL) +
  geom_text(aes(group, feature, label = round(coef,2), color = sig), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "coefficient as compared to co-current patterns") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste(out_dir,"massling.pdf",sep = '/'),width = 6,height = 6)

phylo_used <- subset_samples(phyloseqin_all, species == sp[3])
phylo_used <- subset_samples(phylo_used,Pattern_class %in% pattern_type[c(1:3)])
phylo_used <- filter_taxa(phylo_used, function(x) sum(x > 0) > 0, TRUE)
d_otu <- as.data.frame(t(otu_table(phylo_used)))
d_otu <- d_otu[,grep("g__",colnames(d_otu))]

tmp_list <- combine_list(meta_used$Pattern_class)
pattern_res <- data.frame()
pattern_res_sig <- c()
for(i in 1:length(tmp_list)){
  tmp_list[[i]] = sort(tmp_list[[i]])
  tmp_meta <- meta_used %>% filter(Pattern_class %in% tmp_list[[i]])
  tmp_d <- d_otu[rownames(d_otu) %in% tmp_meta$ID,]
  tmp_fn <- paste(tmp_list[[i]][1],tmp_list[[i]][2],sep = "-")
  tmp_fn <- paste(tmp_fn,"genus",sep = "_")
  tmp_ref <- paste("Pattern_class",tmp_list[[i]][1],sep = ",")
  Maaslin2(
    input_data = tmp_d, 
    input_metadata = tmp_meta, 
    min_prevalence = 0,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    output = paste(out_dir,tmp_fn,sep = '/'), 
    fixed_effects = c("Pattern_class","sex","age","bmi"),
    reference = tmp_ref,
  )
  tmp_res <-  read_tsv(paste(c(out_dir,tmp_fn,"all_results.tsv"),collapse = '/'))
  tmp_res <- tmp_res %>% 
    filter(value == tmp_list[[i]][2]) %>% 
    dplyr::select(feature,coef,qval) %>% mutate(group = tmp_fn)
  pattern_res <- bind_rows(pattern_res,tmp_res)
  tmp_sig <- tmp_res %>% filter(qval < 0.25)
  tmp_sig <- tmp_sig$feature
  pattern_res_sig <- c(pattern_res_sig,tmp_sig)
}
pattern_res_sig <- unique(pattern_res_sig)
pattern_res <- pattern_res %>% filter(feature %in% pattern_res_sig) 
pattern_res$sig <- ifelse(pattern_res$qval<0.25,"aquamarine3", "black")

tmp_lo = floor(min(pattern_res$coef))
tmp_up = ceiling(max(pattern_res$coef))
tmp_mid = (tmp_lo + tmp_up)/2

ggplot(pattern_res, aes(x = group, y = feature, fill = coef)) + 
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = tmp_mid, limit = c(tmp_lo,tmp_up),
                       name = NULL) +
  geom_text(aes(group, feature, label = round(coef,2), color = sig), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "coefficient as compared to co-current patterns") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste(out_dir,"massling_genus.pdf",sep = '/'),width = 6,height = 6)

