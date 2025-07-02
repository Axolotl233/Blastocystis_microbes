out_dir = "output/subtype_class"
dir.create(out_dir,showWarnings = F)
d_pattern <- read.csv(paste(out_dir,"cvalue/5/tern.csv",sep = "/")) %>% select(ID,Pattern)
#data_base_maf <- d_1 %>% filter(species == "mfas")
d_pattern <- left_join(d_1,d_pattern,by = "ID")
d_pattern$Pattern <- ifelse(is.na(d_pattern$Pattern),"Unclass",d_pattern$Pattern)
d_pattern <- d_pattern %>%
  mutate(
    Pattern_class = case_when(
      Pattern == "Pattern_1" ~ "ST1_dom",
      Pattern == "Pattern_2" ~ "ST3_dom",
      Pattern == "Pattern_3" ~ "Balance",
      Pattern == "Pattern_4" ~ "ST2_dom",
      TRUE ~ "Unclass"
    )
  )

d_pattern <- as.data.frame(d_pattern)
row.names(d_pattern) <- d_pattern$ID
d_pattern$Pattern_class <- factor(d_pattern$Pattern_class, levels = pattern_type)

d_tmp <- d_pattern[,c(4:7)]
tmp_pca_res <- prcomp(d_tmp, scale = F)
tmp_pca_plot <- as.data.frame(tmp_pca_res$x)
tmp_pca_plot$ID <- row.names(tmp_pca_plot)
tmp_pca_plot <- left_join(tmp_pca_plot,d_pattern,by = "ID")

ggplot()+
  geom_point(data = tmp_pca_plot,
             mapping = aes(x=PC1,y=PC2,color=Pattern_class,shape=species),
             size=3,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PC1"),y=c("PC2"),title = "PCA")+
  stat_ellipse(data=tmp_pca_plot,geom = "polygon",
               aes(x=PC1,y=PC2,fill=Pattern_class),
               color =NA,alpha=0.1,level = 0.9)+
  scale_fill_manual(values = color0)+
  scale_color_manual(values = color0)+
  scale_shape_manual(values = shape0)

#pca_plot_n <- paste(outdir_o,kcluster,"pca.pdf",sep = "/")
ggsave(paste(out_dir,"cvalue.pca.pdf",sep = "/"),width = 6,height = 5)

d_pattern_stat <- d_pattern %>% group_by(Pattern_class,species) %>% summarise(Count = n())
ggplot(d_pattern_stat) +
  geom_bar(aes(x = Pattern_class,y=Count,fill = species),stat = "identity",position = "dodge",
           width = 0.6)+
  geom_text(aes(x = Pattern_class,y=Count + 1 ,label  = Count))+
  scale_fill_manual(
    values = color0
  )+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 30,vjust = 1,hjust = 1),
    legend.position = "none",
    panel.grid = element_blank()
  )
ggsave(paste(out_dir,"cvalue.count.pdf",sep = "/"),width = 5,height = 3)

d_pattern_l <- d_pattern %>% select(!parasite_load)%>%select(!Pattern) %>% select(!ST5) %>%
  gather(key = "Subtype",value = "RA",ST1,ST2,ST3)
tmp_list1 <- combine_list(d_pattern_l$Subtype)
d_pattern_l$Pattern_class <- factor(d_pattern_l$Pattern_class, levels = pattern_type)
d_pattern_l$Pattern_merge <- paste(d_pattern_l$Pattern_class,d_pattern_l$Subtype,sep = '-')

d_pattern_l_p <- d_pattern_l  %>% filter (Pattern_class != "Unclass")  %>% filter (Pattern_class != "ST2_dom")
d_pattern_l_p_stat <- summary_stat(d_pattern_l_p,5,6)
d_pattern_l_p_stat$Pattern_class <- gsub("-.*","",d_pattern_l_p_stat$Group)
d_pattern_l_p_stat$Subtype <- gsub(".*-","",d_pattern_l_p_stat$Group)

ggplot(d_pattern_l_p,
       aes(x = Subtype, y = RA,color = Pattern_class)) +
  geom_line(aes(group = ID),linetype = "dashed",alpha = 0.3)+
  geom_point(size = 3,alpha = 0.5)+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid = element_blank())+
  facet_grid(species~Pattern_class)+
  scale_color_manual(
    values = color0
  )+
  labs(x="Subtype",y="Abundance")
ggsave(paste(out_dir,"cvalue.line.pdf",sep = "/"),width = 10,height = 6)

ggplot(d_pattern_l_p_stat,aes(x = Subtype, y = Mean))+
  geom_ribbon(aes(ymin = Mean-Sd,ymax = Mean+Sd,group = Pattern_class),fill = "grey90"
  )+
  geom_point(aes(color = Pattern_class),size = 3)+
  geom_line(aes(group = Pattern_class,color = Pattern_class),linewidth = 0.5)+
  theme_cowplot()+
  theme(legend.position = "none")+
  facet_wrap(Pattern_class~.,nrow = 1)+
  scale_color_manual(
    values = color0
  )
ggsave(paste(out_dir,"cvalue.line_ribbon.pdf",sep = "/"),width = 8,height = 6)

d_pattern_l2 <- d_pattern_l %>% filter (Pattern_class != "Unclass") %>% 
  filter (Pattern_class != "ST2_dom")
tmp_list2 <- combine_list(d_pattern_l2$Pattern_class)
ggplot(d_pattern_l2,
       aes(x = Pattern_class, y = RA,color = Pattern_class)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size = 2.5,height = 0,width = 0.15,alpha = 0.7) +
  theme_bw()+
  theme( panel.grid = element_blank(),legend.position = "none")+
  facet_grid(species~Subtype)+
  #stat_compare_means(comparisons = tmp_list2,
  #                  method = "wilcox.test", label = "p.format")+
  scale_color_manual(
    values = color0
  )
ggsave(paste(out_dir,"cvalue.class2type.pdf",sep = "/"),height = 6,width = 8)

ggstatsplot::ggbarstats(
  data = d_pattern  %>% filter (Pattern_class != "Unclass")  %>% filter (Pattern_class != "ST2_dom"),
  y = species,
  x = Pattern_class
)+ coord_flip()+
  scale_fill_manual(
    values = rev(color0[c(1:3)])
  )
ggsave(paste(out_dir,"cvalue.class2sp.pdf",sep = "/"),height = 3,width = 10)

d_3 <- d1 %>% select(ID,sex,age,bmi)
d_pattern_all <-  inner_join(d_pattern,d_3,by = "ID")
d_pattern_mfas <- d_pattern_all %>% filter(species == "mfas")
d_pattern_mmul <- d_pattern_all %>% filter(species == "mmul")

#dist_subtype <- vegdist(d_pattern_mfas[,4:7],method = "euclidean")
#dist_subtype_mfas <- as.dist((1-cor(t(d_pattern_mfas[,4:7]),method="pearson"))/2)
dist_subtype_mfas <- vegdist(d_pattern_mfas[,4:7],method = "euclidean")
r_subtype_permanova_mfas <-  adonis2(dist_subtype_mfas ~ sex + bmi ++age , 
                                     d_pattern_mfas,by = "terms",permutations = 999)
r_subtype_permanova_mfas$item <- row.names(r_subtype_permanova_mfas)
r_subtype_permanova_mfas_p <- r_subtype_permanova_mfas[!is.na(r_subtype_permanova_mfas$F),]
r_subtype_permanova_mfas_p$item <- factor(r_subtype_permanova_mfas_p$item,levels = r_subtype_permanova_mfas_p$item)
r_subtype_permanova_mfas_p$class <- ifelse(r_subtype_permanova_mfas_p$`Pr(>F)`>0.05,"no_sig","sig")
r_subtype_permanova_mfas_p$species <- sp[3]

p_a <- ggplot(r_subtype_permanova_mfas_p)+
  geom_bar(aes(x = item, y = R2,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2+ 0.02,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none")+
  labs(x = "", y = "R2",title = sp[1])

dist_subtype_mmul <- vegdist(d_pattern_mmul[,4:7],method = "euclidean")
r_subtype_permanova_mmul <-  adonis2(dist_subtype_mmul ~ sex  + bmi  +age, 
                                     d_pattern_mmul,by = "terms",permutations = 999)
r_subtype_permanova_mmul$item <- row.names(r_subtype_permanova_mmul)
r_subtype_permanova_mmul_p <- r_subtype_permanova_mmul[!is.na(r_subtype_permanova_mmul$F),]
r_subtype_permanova_mmul_p$item <- factor(r_subtype_permanova_mmul_p$item,levels = r_subtype_permanova_mmul_p$item)
r_subtype_permanova_mmul_p$class <- ifelse(r_subtype_permanova_mmul_p$`Pr(>F)`>0.05,"no_sig","sig")
r_subtype_permanova_mmul_p$species <- sp[4]

p_b <- ggplot(r_subtype_permanova_mmul_p)+
  geom_bar(aes(x = item, y = R2,fill = class),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2+ 0.005,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none")+
  labs(x = "", y = "R2",title = sp[2])
plot_grid(p_a,p_b,nrow = 2)
ggsave(paste(out_dir,"cvalue_cluster.adonis2.pdf",sep = "/"),height = 6,width = 5)

r_subtype_permanova_p <- rbind(r_subtype_permanova_mfas_p,r_subtype_permanova_mmul_p)
ggplot(r_subtype_permanova_p)+
  geom_bar(aes(x = item, y = R2,fill = species),stat = "identity",position = "dodge",width = 0.6)+
  geom_text(aes(x = item, y = R2+ 0.002,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_light()+
  scale_fill_manual(values = c("grey80","#5698c4"))+
  theme(legend.position = "none",
        panel.grid = element_blank())+
  labs(x = "", y = "R2")
ggsave(paste(out_dir,"cvalue_cluster.adonis2.pdf",sep = "/"),height = 3,width = 5)

d_pattern_all$parasite_load_log <- log2(d_pattern_all$parasite_load + 1)
d_pattern_all_p <- d_pattern_all %>% select(ID,Pattern_class,species,parasite_load_log,age,bmi) %>% 
  gather(key = "class",value = "value",age,bmi,parasite_load_log) %>%
  filter (Pattern_class %in% pattern_type[c(1:3)])
tmp_list2 <- combine_list(d_pattern_all_p$Pattern_class)

ggplot(d_pattern_all_p %>% filter(species == "mfas"),
       aes(x = Pattern_class, y = value,color = Pattern_class)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size = 2.5,height = 0,width = 0.15,alpha = 0.7) +
  theme_cowplot()+
  theme(panel.grid = element_blank(),legend.position = "none")+
  facet_wrap(.~class,scales = "free",nrow = 1)+
  stat_compare_means(comparisons = tmp_list2,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(
    values = color0
  )
ggsave(paste(out_dir,"cvalue.var2pattern.mfas.pdf",sep = "/"),height = 5,width = 8)

ggplot(d_pattern_all_p %>% filter(species == "mfas",class == "parasite_load_log"),
       aes(x = Pattern_class, y = value,color = Pattern_class)) +
  geom_boxplot(outlier.shape = NA,width = 0.7)+
  geom_jitter(size = 3.5,height = 0,width = 0.15,alpha = 0.7) +
  theme_cowplot()+
  theme(panel.grid = element_blank(),legend.position = "none")+
  stat_compare_means(comparisons = tmp_list2,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(
    values = color0
  )+ canvas(5, 4, units = "in")
ggsave(paste(out_dir,"cvalue.para2pattern.mfas.pdf",sep = "/"),height = 4,width = 5)

ggplot(d_pattern_all_p %>% filter(species == "mmul"),
       aes(x = Pattern_class, y = value,color = Pattern_class)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size = 2.5,height = 0,width = 0.15,alpha = 0.7) +
  theme_cowplot()+
  theme(panel.grid = element_blank(),legend.position = "none")+
  facet_wrap(.~class,scales = "free",nrow = 1)+
  stat_compare_means(comparisons = tmp_list2,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(
    values = color0
  )
ggsave(paste(out_dir,"cvalue.var2pattern.mmul.pdf",sep = "/"),height = 5,width = 8)

ggstatsplot::ggbarstats(
  data = d_pattern_all  %>% 
    filter (Pattern_class != "Unclass")%>% 
    filter (Pattern_class != "ST2_dom") %>%
    filter (species == "mfas"),
  y = sex,
  x = Pattern_class
)+ coord_flip()+
  scale_fill_manual(
    values = rev(color0[c(1:3)])
  )
ggsave(paste(out_dir,"cvalue.sex2pattern.mfas.pdf",sep = "/"),height = 4,width = 10)

ggstatsplot::ggbarstats(
  data = d_pattern_all  %>% 
    filter (Pattern_class != "Unclass")%>% 
    filter (Pattern_class != "ST2_dom") %>%
    filter (species == "mmul"),
  y = sex,
  x = Pattern_class
)+ coord_flip()+
  scale_fill_manual(
    values = rev(color0[c(1:3)])
  )
ggsave(paste(out_dir,"cvalue.sex2pattern.mmul.pdf",sep = "/"),height = 4,width = 10)
write_csv(x = d_pattern_all,paste(out_dir,"subtype_pattern.csv",sep = "/"))
