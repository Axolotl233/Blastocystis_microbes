table_pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    p <- kruskal.test(y ~ g)$p.value
  } else {
    p <- fisher.test(table(y, g))$p.value
  }
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

library(table1)
library(UpSetR)
library(ape)
library(treeio)
library(tidytree)
library(ggtree)
library(ggtreeExtra)
library(ggtext)
library(ggnewscale)
library(ggsci)

out_dir = "output/metadata"
dir.create(out_dir,showWarnings = F)
d <- read_csv("input/Blastocystis_meta_final.fix.csv")

tree <- read.nexus(in_tree)
tree <- tidytree::as.treedata(tree)
tree_c <- data.frame(OTU = tree@phylo$tip.label)
tree_anno <- read.csv(in_amp_taxon,row.names = 1,header = T)
tree_anno <- tree_anno[row.names(tree_anno) %in% tree_c$OTU,colnames(tree_anno) %in% d$ID]
tree_otu <- row.names(tree_anno)
tree_anno <- as.data.frame(t(tree_anno))
tree_anno$ID <- row.names(tree_anno)
tree_anno <- left_join(tree_anno,d %>% select(ID,species),by = "ID")
tree_anno_l <- gather(tree_anno,key = "OTU",value = "num",tree_otu)
tree_anno_l_stat_1 <- tree_anno_l %>% group_by(OTU,species) %>% summarise(n_sample = length(which(num != 0)))
#tree_c <- left_join(tree_c,tree_anno_l_stat_1,by = "ID")
#tree_c[is.na(tree_c$n_sample),2] <-  120
ggtree(tree,layout = "circular",size=0.1,branch.length = "none")+
  geom_tiplab(align = T,size=2)+
  geom_fruit(
    data=tree_anno_l_stat_1,
    geom=geom_bar,
    mapping=aes(y=OTU, x=n_sample,fill = species),
#    position = "stack",
    pwidth=0.6,
    stat="identity",
    offset = 0.3
  )
ggsave(paste(out_dir,"tree.pdf",sep = "/"),width = 7,height = 5)
tree_anno_l_stat_2 <- tree_anno_l %>% group_by(ID,species) %>% summarise(n_sample = length(which(num != 0)))

ggplot(tree_anno_l_stat_2)+
  geom_histogram(aes(x = n_sample),binwidth = 1,boundary = 0)+
  scale_x_continuous(breaks = seq(0,40,10),limits = c(0,40))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  facet_wrap(.~species)+
  labs(x="Number of OTU in single sample",y="Number of sample")
ggsave(paste(out_dir,"otu_sample_stat.pdf",sep = "/"),width = 10,height = 3)

d <- inner_join(d,tree_anno_l_stat_2 %>% select(!species),by = "ID")
dd <- d %>% filter(Pattern_class %in% pattern_type[c(1:3)])
dd$Pattern_class <- factor(dd$Pattern_class,levels = pattern_type[c(1:3)])
formula_tmp <- y ~ x 
ggplot(dd,aes(x=parasite_load_log,y=n_sample)) + 
  stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
  geom_point(aes(color = Pattern_class),size=3,alpha = 0.7)+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid = element_blank(),
  )+
  facet_wrap(.~species)+
  scale_color_manual(
    values = color0
  )+
  labs(x = "Log2(Parasite load + 1)",y = "OTU number")
ggsave(paste(out_dir,"cor.otu_sample_parasite_load.pdf",sep = "/"),width = 10,height = 3)
d_1 <- d %>% select(ID,species,parasite_load,ST1,ST2,ST3,ST5)
d_2 <- apply(d_1[,4:7],c(1,2),function(x){if(x != 0){x = 1}else{x = 0}})

pdf(paste(out_dir,"subtype_upsets.pdf",sep = "/"),width = 6,height = 5)
upset(as.data.frame(d_2), order.by = c('freq'),
      ,sets.bar.color = "#56B4E9")
dev.off()

d_2_l <- gather(d_1[,c(2,4:7)],key = "Subtype",value = "value",ST1,ST2,ST3,ST5)
#p_a <- ggplot(d_2_l,aes(x=Subtype,y=value,color=species,fill = species)) +
#  stat_summary(fun=mean, geom="bar",position = position_dodge(0.7),width=0.65)+
#  stat_summary(fun.data = mean_se, geom = "errorbar",color = "black", 
#               position = position_dodge(0.7),width = 0.2)+
#  scale_y_continuous(breaks = seq(0,100,10),)+
#  theme_bw()+
#  theme(panel.grid = element_blank(),legend.position = "none")+
#  labs(x="Subtype",y="Abundance")

p_a <- ggplot(d_2_l,aes(x=Subtype,y=value,color=species)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75,
                                              jitter.width = 0.15,
                                              jitter.height = 0),
              size = 2.5,alpha = 0.7
              )+
  scale_y_continuous(breaks = seq(0,100,10),)+
  theme_bw()+
  theme(panel.grid = element_blank()
        ,legend.position = "none"
        )+
  labs(x="Subtype",y="Abundance")

p_b <- ggplot(d_2_l %>% filter(value > 0) %>% group_by(species,Subtype) %>% summarise(n = n()))+
  geom_bar(aes(x=Subtype,y = n,color=species,fill = species),position = position_dodge(0.7),
           stat = "identity",width=0.65)+
  geom_text(aes(x=Subtype,y = n+3,label = n),position = position_dodge(0.7))+
  theme_bw()+
  theme(panel.grid = element_blank(),legend.position = "none")+
  labs(x="Subtype",y="Sample number")

plot_grid(p_a,p_b)
ggsave(paste(out_dir,"subtype_stat.pdf",sep = "/"),width = 9,height = 3)

d_2 <- as.data.frame(d_2)
c_colname_d_2 <- colnames(d_2)
c_tmp <- c()
for(i in 1:nrow(d_2)){
  c_tmp[i] <- paste(c_colname_d_2[d_2[i,] != 0],collapse = "+")
}
d_2$type <- c_tmp
d_2$species <- d_1$species
d_2_stat <- d_2 %>% group_by(type) %>% summarise(n = n())
d_2_stat$fraction <- d_2_stat$n/126
d_2_stat = d_2_stat[order(d_2_stat$fraction), ]
d_2_stat$ymax = cumsum(d_2_stat$fraction)
d_2_stat$ymin = c(0, head(d_2_stat$ymax, n=-1))
c_fill = pal_npg(palette = c("nrc"), alpha = 1)(nrow(d_2_stat))
d_2_stat$type <- factor(d_2_stat$type,levels = rev(unique(d_2_stat$type)))
ggplot(d_2_stat, aes(fill=type, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +  
  geom_rect(colour="White") +
  coord_polar(theta="y")+
  scale_fill_manual(labels = paste("<span style='color:",
                                   c_fill,
                                   "'>",
                                   rev(unique(d_2_stat$type)),
                                   "</span>"),
                    values = c_fill)+
  theme_void()+
  geom_text_repel(aes(label=paste(round(d_2_stat$fraction*100,2),"%"),x=4,y=
                        (ymin+ymax)/2),inherit.aes = F)+
  xlim(c(0, 4)) +
  theme(legend.title = element_blank())+
  theme(legend.text=element_markdown(size=12))
ggsave(paste(out_dir,"subtype_class_ratio.pdf",sep = "/"),width = 7,height = 5)
d_2_stat4 <- d_2 %>% group_by(type,species) %>% summarise(n = n())

d_2_stat3 <-  d_2_stat4 %>% filter(species == "mfas")
d_2_stat3$fraction <- d_2_stat3$n/100
d_2_stat3 = d_2_stat3[order(d_2_stat3$fraction), ]
d_2_stat3$ymax = cumsum(d_2_stat3$fraction)
d_2_stat3$ymin = c(0, head(d_2_stat3$ymax, n=-1))
c_fill = pal_npg(palette = c("nrc"), alpha = 1)(nrow(d_2_stat3))
d_2_stat3$type <- factor(d_2_stat$type,levels = rev(unique(d_2_stat3$type)))
ggplot(d_2_stat3, aes(fill=type, ymax=ymax, ymin=ymin, xmax=1, xmin=0)) +  
  geom_rect(colour="White") +
  coord_polar(theta="y")+
  scale_fill_manual(values = c_fill)+
  theme_void()+
  geom_text_repel(aes(label=paste(round(d_2_stat3$fraction*100,2),"%"),x=1,y=
                   (ymin+ymax)/2),inherit.aes = F)+
  xlim(c(0, 4)) +
  theme(legend.title = element_blank())+
  theme(legend.text=element_markdown(size=12))

ggsave(paste(out_dir,"subtype_class_ratio_mfas.pdf",sep = "/"),width = 7,height = 5)

d_2_stat2 <- d_2_stat4 %>% filter(species == "mmul")
d_2_stat2$fraction <- d_2_stat2$n/26
d_2_stat2 = d_2_stat2[order(d_2_stat2$fraction), ]
d_2_stat2$ymax = cumsum(d_2_stat2$fraction)
d_2_stat2$ymin = c(0, head(d_2_stat2$ymax, n=-1))
c_fill = pal_npg(palette = c("nrc"), alpha = 1)(nrow(d_2_stat2))
d_2_stat2$type <- factor(d_2_stat$type,levels = rev(unique(d_2_stat2$type)))
ggplot(d_2_stat2, aes(fill=type, ymax=ymax, ymin=ymin, xmax=1, xmin=0)) +  
  geom_rect(colour="White") +
  coord_polar(theta="y")+
  scale_fill_manual(values = c_fill)+
  theme_void()+
  geom_text_repel(aes(label=paste(round(d_2_stat2$fraction*100,2),"%"),x=1,y=
                        (ymin+ymax)/2),inherit.aes = F)+
  xlim(c(0, 4)) +
  theme(legend.title = element_blank())+
  theme(legend.text=element_markdown(size=12))
ggsave(paste(out_dir,"subtype_class_ratio_mmul.pdf",sep = "/"),width = 7,height = 5)
d_2_stat$species <- "all"

d_2_stat_all <- rbind(d_2_stat,d_2_stat2,d_2_stat3) %>% select(type,fraction,species)
d_2_stat_all[17,1] <- c("tmp")
d_2_stat_all[17,2] <- 1
d_2_stat_all[17,3] <- "un"
d_2_stat_all$species <- factor(d_2_stat_all$species,levels = c("un","mmul","mfas","all"))

c_fill = pal_npg(palette = c("nrc"), alpha = 1)(length(unique(d_2_stat_all$type)))
ggplot(d_2_stat_all, aes(fill=type, y = fraction, x = species))+
  geom_bar(stat = "identity",width = 0.65)+
  coord_polar(theta = "y")+
  theme_void()+
  scale_fill_manual(values = c_fill)

ggsave(paste(out_dir,"subtype_class_ratio_all.pdf",sep = "/"),width = 7,height = 5)


c_tmp1 <- c("ID","species","age","sex","bmi","ST1","ST2","ST3","ST1_abs_log","ST2_abs_log","ST3_abs_log","parasite_load_log")
d_3 <- d %>% select(all_of(c_tmp1))
d_3_mfas <- d_3 %>% filter(species == "mfas")
d_3_mmul <- d_3 %>% filter(species == "mmul")

wilcox.test(ST1_abs_log~sex,d_3_mfas)
wilcox.test(ST2_abs_log~sex,d_3_mfas)
wilcox.test(ST3_abs_log~sex,d_3_mfas)
wilcox.test(age~sex,d_3_mfas)

wilcox.test(ST1_abs_log~sex,d_3_mmul)
wilcox.test(ST2_abs_log~sex,d_3_mmul)
wilcox.test(ST3_abs_log~sex,d_3_mmul)
wilcox.test(age~sex,d_3_mmul)

c_tmp2 <- c("age","ST1_abs_log","ST2_abs_log","ST3_abs_log","parasite_load_log")
p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d_3_mfas[,c_tmp2[c(1,i)]]
  colnames(d_tmp1) <- c("v1","v2")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )+
    labs(x = "Age",y =c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[4]],nrow = 2)
ggsave(paste(out_dir,"cor.abs_parasite2age.mfas.pdf",sep = "/"),width = 7,height = 7)

p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d_3_mmul[,c_tmp2[c(1,i)]]
  colnames(d_tmp1) <- c("v1","v2")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )+
    labs(x = "Age",y =c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[4]],nrow = 2)
ggsave(paste(out_dir,"cor.abs_parasite2age.mmul.pdf",sep = "/"),width = 7,height = 7)

p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d_3[,c_tmp2[c(1,i)]]
  colnames(d_tmp1) <- c("v1","v2")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )+
    labs(x = "Age",y =c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[4]],nrow = 2)
ggsave(paste(out_dir,"cor.abs_parasite2age.all.pdf",sep = "/"),width = 7,height = 7)

c_tmp2 <- c("age","ST1","ST2","ST3")
p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d_3_mfas[,c_tmp2[c(1,i)]]
  colnames(d_tmp1) <- c("v1","v2")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )+
    labs(x = "Age",y =c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],nrow = 2)
ggsave(paste(out_dir,"cor.ra_parasite2age.mfas.pdf",sep = "/"),width = 7,height = 7)

p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d_3_mmul[,c_tmp2[c(1,i)]]
  colnames(d_tmp1) <- c("v1","v2")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )+
    labs(x = "Age",y =c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],nrow = 2)
ggsave(paste(out_dir,"cor.ra_parasite2age.mmul.pdf",sep = "/"),width = 7,height = 7)

p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d_3[,c_tmp2[c(1,i)]]
  colnames(d_tmp1) <- c("v1","v2")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )+
    labs(x = "Age",y =c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],nrow = 2)
ggsave(paste(out_dir,"cor.ra_parasite2age.all.pdf",sep = "/"),width = 7,height = 7)

p_list <- list()
for(i in 2:length(c_tmp2)){
  #d_tmp1 <- d %>% filter(species == sp[3])
  d_tmp1 <- d[,c("parasite_load_log",c_tmp2[i],"Pattern_class")]
  d_tmp1 <- d_tmp1 %>% filter(Pattern_class %in% pattern_type[c(1:3)])
  d_tmp1$Pattern_class <- factor(d_tmp1$Pattern_class, levels = pattern_type[c(1:3)])
  colnames(d_tmp1) <- c("v1","v2","V3")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy")+
    geom_point(aes(color = V3),size=3,alpha = 0.7)+
    scale_y_continuous(breaks = seq(0,100,20),limits = c(-20,120))+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
          legend.title = element_blank())+
    scale_color_manual(
      values = color0
    )+
    labs(x = "log2(parasite_load+1)",y = c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[3]],p_list[[2]],nrow = 2)
ggsave(paste(out_dir,"cor.ra_parasite2load.all.pdf",sep = "/"),width = 10,height = 7)

p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d[d$species == sp[3],c("parasite_load_log",c_tmp2[i],"Pattern_class")]
  d_tmp1 <- d_tmp1 %>% filter(Pattern_class %in% pattern_type[c(1:3)])
  d_tmp1$Pattern_class <- factor(d_tmp1$Pattern_class, levels = pattern_type[c(1:3)])
  colnames(d_tmp1) <- c("v1","v2","V3")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy",alpha = 0.1)+
    geom_point(aes(color = V3),size=3,alpha = 0.7)+
    scale_y_continuous(breaks = seq(0,100,20),limits = c(-20,120))+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
          legend.title = element_blank())+
    scale_color_manual(
      values = color0
    )+
    labs(x = "log2(parasite_load+1)",y = c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[3]],p_list[[2]],nrow = 2)
ggsave(paste(out_dir,"cor.ra_parasite2load.mfas.pdf",sep = "/"),width = 9.5,height = 7)

p_list <- list()
for(i in 2:length(c_tmp2)){
  d_tmp1 <- d[d$species == sp[4],c("parasite_load_log",c_tmp2[i],"Pattern_class")]
  d_tmp1 <- d_tmp1 %>% filter(Pattern_class %in% pattern_type[c(1:3)])
  d_tmp1$Pattern_class <- factor(d_tmp1$Pattern_class, levels = pattern_type[c(1:3)])
  colnames(d_tmp1) <- c("v1","v2","V3")
  formula_tmp <- y ~ x 
  p <- ggplot(d_tmp1,aes(x=v1,y=v2)) + 
    stat_smooth(method='lm',formula = formula_tmp,color = "navy",alpha = 0.1)+
    geom_point(aes(color = V3),size=3,alpha = 0.7)+
    scale_y_continuous(breaks = seq(0,100,20),limits = c(-20,120))+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
          legend.title = element_blank())+
    scale_color_manual(
      values = color0
    )+
    labs(x = "log2(parasite_load+1)",y = c_tmp2[i])
  p_list[[i-1]] <- p
}
plot_grid(p_list[[1]],p_list[[2]],p_list[[3]],nrow = 2)
ggsave(paste(out_dir,"cor.ra_parasite2load.mmul.pdf",sep = "/"),width = 9.5,height = 7)


