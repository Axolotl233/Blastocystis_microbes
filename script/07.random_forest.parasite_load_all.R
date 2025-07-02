library(randomForest)
library(rfUtilities)
library(rfPermute)
library(VSURF)
library(Boruta)
library(caret)
library(pROC)

set.seed(19970519)
out_dir <- paste("output/randomforest","all",sep = "/")
dir.create(out_dir,showWarnings = F,recursive = T)

load("workfile/z.phyloseq.Rdata")
phyloseqin_used <-  subset_samples(phyloseqin, species %in% c(sp[3],sp[4]))
phyloseqin_used <- subset_samples(phyloseqin_used,Pattern_class %in% pattern_type[c(1:3)])
tmp_otu <- as.data.frame(otu_table(phyloseqin_used))
tmp_otu <- noise_removal(tmp_otu,percent=0.2,low = 0.005, method = "pre_cut")
#tmp_otu <- sweep(tmp_otu,2,colSums(tmp_otu),"/") * 100
tmp_meta <- metadatadf %>% dplyr::select(ID,parasite_load_log)
tmp_meta <- tmp_meta[tmp_meta$ID %in% colnames(tmp_otu), ]
ggplot(tmp_meta) +
  geom_histogram(aes(x = parasite_load_log),boundary = 0,binwidth = 0.5)+
  scale_x_continuous(breaks = seq(0,15,1),limits = c(0,15))+
  scale_y_continuous(breaks = seq(0,12,1))+
  theme_classic()
#tmp_meta <- tmp_meta %>% filter(parasite_load_log < 9) 
tmp_otu <- tmp_otu[,colnames(tmp_otu)%in%tmp_meta$ID]

identical(tmp_meta$ID,colnames(tmp_otu))
tmp_data <- as.data.frame(t(tmp_otu))
tmp_data$ID <- row.names(tmp_data)
tmp_data <- inner_join(tmp_data,tmp_meta,by="ID") 
row.names(tmp_data) <- tmp_data$ID
tmp_data <- tmp_data %>% dplyr::select(!ID) 

model_classify2 <- randomForest(x=tmp_data[,1:(ncol(tmp_data)-1)] , 
                               y=tmp_data[,"parasite_load_log"],
                               ntree=1000, 
                               importance=TRUE, 
                               proximity =TRUE, 
                               replace = TRUE,
#                               mtry=which.min(tmp_errrate),
                               na.action=na.omit,
                               oob_score=TRUE
)
plot(model_classify2)
model_imp <- importance(model_classify2)
model_imp <- data.frame(predictors = rownames(model_imp), model_imp)
model_imp_sort <- arrange(model_imp, desc(X.IncMSE))
model_classify_pm <- rfPermute(x=tmp_data[,1:(ncol(tmp_data)-1)] , 
                               y=tmp_data[,"parasite_load_log"],
                               ntree=1000, 
                               importance=TRUE, 
                               proximity =TRUE, 
                               replace = TRUE,
                               #mtry=which.min(tmp_errrate),
                               na.action=na.omit,
                               oob_score=TRUE
)
model_imp_pvalue <- as.data.frame(model_classify_pm$pval[,,2]) %>% dplyr::select(`%IncMSE`)
model_imp_pvalue$predictors <- row.names(model_imp_pvalue)
colnames(model_imp_pvalue) <- c("Pvalue","predictors")
model_imp_sort <- inner_join(model_imp_sort,model_imp_pvalue,by = "predictors")
model_imp_sort <- model_imp_sort %>% mutate(
  sig = case_when(
    Pvalue < 0.01 ~ "**",
    Pvalue < 0.05 & Pvalue >= 0.01 ~ "*",
    TRUE ~ "ns"
  )
)
colnames(model_imp_sort) <- c("predictors","IncMSE","IncNodePurity","Pvalue","sig")
boruta_result <- Boruta(x=tmp_data[,1:(ncol(tmp_data)-1)] , 
                        y=tmp_data[,"parasite_load_log"],
                        doTrace = 2)

boruta_result_table <- as.data.frame(boruta_result$finalDecision)
colnames(boruta_result_table) <- "Class"
boruta_result_table$predictors <- row.names(boruta_result_table)
model_imp_sort <- inner_join(model_imp_sort,boruta_result_table,by = "predictors")
write_csv(model_imp_sort,file = paste(out_dir,"model_imp.csv",sep = '/'))

model_imp_sort_top <- model_imp_sort %>% filter(Pvalue < 0.05) %>% 
  filter(Class != "Rejected") 
#model_imp_sort_top$predictors <- factor(model_imp_sort_top$predictors,levels = rev(model_imp_sort_top$predictors))
model_imp_sort_top$Species <- gsub("s__","",model_imp_sort_top$predictors)
model_imp_sort_top$Species <- factor(model_imp_sort_top$Species,levels = rev(model_imp_sort_top$Species))

ggplot(model_imp_sort_top) +
  geom_bar(aes(x = Species, y = IncMSE,fill = Class),stat = "identity") +
  geom_text(aes(x = Species, y = IncMSE+0.25,label = sig),angle = 270)+
  coord_flip()+
  labs(title= "The important species",x = "Species")+
  theme_bw()+
  theme(strip.background = element_blank(),panel.grid = element_blank())+
  scale_fill_manual(
    values = c("grey80","grey50")
  )
ggsave(paste(out_dir,"rf_imp.pdf",sep = '/'),width = 8,height = 8)

plot_l1 <- list()
tmp_meta2 <- metadatadf %>% filter(ID %in% tmp_meta$ID)
for(i in 1:nrow(model_imp_sort_top)){
  tmp_t <- as.data.frame(t(tmp_otu))
  tmp_t$ID <- row.names(tmp_t)
  tmp_d <- inner_join(tmp_meta2,tmp_t,by ="ID") %>% dplyr::select(model_imp_sort_top[i,1],ID,Pattern_class)
  
  tmp_c=unique(tmp_d$Pattern_class)
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    mutate(tmp = str_c(V1,",",V2)) %>% 
    dplyr::select(tmp)
  tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
  
  colnames(tmp_d) <- c("value","sample","groups")
  tmp_d$value <- log10(tmp_d$value + 1e-6)
  tmp_d$groups <- factor(tmp_d$groups,levels = pattern_type[c(1:3)])
  p <- ggplot(tmp_d,aes(x=groups,y=value,color = groups))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(width = 0.2,height = 0)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="log10(Ab + 1e-6)",title = model_imp_sort_top[i,1])+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color0)
  plot_l1[[i]] <- p
}
plot_grid(plot_l1[[1]],plot_l1[[2]],plot_l1[[3]],plot_l1[[4]],
          plot_l1[[5]],plot_l1[[6]],plot_l1[[7]],plot_l1[[8]],
          plot_l1[[9]],plot_l1[[10]],plot_l1[[11]],plot_l1[[12]],
          plot_l1[[13]],plot_l1[[14]],plot_l1[[15]],plot_l1[[16]],
          plot_l1[[16]],plot_l1[[18]],plot_l1[[19]],plot_l1[[20]],
          plot_l1[[21]],plot_l1[[22]],plot_l1[[23]],plot_l1[[24]],
          plot_l1[[25]],plot_l1[[26]],plot_l1[[27]],plot_l1[[28]],
          plot_l1[[29]],plot_l1[[30]],plot_l1[[31]],plot_l1[[32]],
          nrow = )
ggsave(paste(out_dir,"rf_imp_sp.pdf",sep = '/'),width = 25,height = 25)

tmp_train_use2 <- sample(nrow(tmp_data), nrow(tmp_data)*0.7)
tmp_data_train_f <- tmp_data[tmp_train_use2,]
tmp_data_train_f <- tmp_data_train_f[,c(model_imp_sort_top$predictors,"parasite_load_log")]
tmp_data_pre <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train_f),]
tmp_data_pre_f <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train_f),c(model_imp_sort_top$predictors,"parasite_load_log")]

model_train_2 <- randomForest(x = tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],
                              y = tmp_data_train_f[,"parasite_load_log"],
                              ntree=1000, 
                              importance=TRUE, 
                              proximity=TRUE
                              )

train_predict <- predict(model_train_2,tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)])
test_predict <- predict(model_train_2,tmp_data_pre_f[,1:(ncol(tmp_data_pre_f)-1)])

defaultSummary(data.frame(obs = tmp_data_train_f$parasite_load_log,
                          pred = train_predict))

defaultSummary(data.frame(obs = tmp_data_pre_f$parasite_load_log,
                          pred = test_predict))

rf_estimation <- function(obs,pred){
  tmp_d <- data.frame(obs = obs,pred = pred)
  tmp_d$ID <- row.names(tmp_d)
  tmp_d[,c(1:2)] <- as.data.frame(apply(tmp_d[,c(1:2)], 2,as.numeric))
  tmp_d$residuals <- tmp_d$obs - tmp_d$pred
  tmp_list <- list()
  p1 <- ggplot(data = tmp_d, aes(x = pred, y = residuals)) +
    geom_point(alpha = 0.8,size = 3) +
    geom_hline(yintercept = 0, color = "red") +
    labs(x = "Predicted", y = "Residuals",
         title = "Residuals vs. Predicted") +
    theme_cowplot()
  
  p2 <- ggplot(data = tmp_d, aes(x = residuals)) +
    geom_histogram(bins = 20, fill = "skyblue", color = "black") +
    labs(title = "Histogram of Residuals", x = "Residuals") +
    theme_cowplot()
  
  tmp_d$rss <- sum((tmp_d$obs - tmp_d$pred)^2)
  tmp_d$tss <- sum((tmp_d$obs - mean(tmp_d$obs))^2)
  r_squared <- 1 - tmp_d$rss/tmp_d$tss
  
  p3 <- ggplot(data = NULL, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.8,size = 3) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    #annotate("text", x = min(tmp_d$obs), y = max(tmp_d$pred), 
     #        label = paste0("R² = ", round(r_squared, 3)),
     #        hjust = 0, vjust = 0, size = 5,) +
    labs(x = "Actual", y = "Predicted", 
         title = "Predicted vs. Actual (with R²)",
         subtitle = paste0("R² = ", round(r_squared, 3))) +
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_x_continuous(breaks = seq(0,14,1),limits = c(0,14))+
    scale_y_continuous(breaks = seq(0,14,1),limits = c(0,14))
  tmp_list[[1]] = p1
  tmp_list[[2]] = p2
  tmp_list[[3]] = p3
  return(tmp_list)
}
r_list_train <- rf_estimation(tmp_data_train_f$parasite_load_log,train_predict)
r_list_pre <- rf_estimation(tmp_data_pre_f$parasite_load_log,test_predict)

plot_grid(r_list_train[[1]],r_list_train[[2]],r_list_train[[3]],nrow = 2,align = "hv")
ggsave(paste(out_dir,"rf_train.pdf",sep = '/'),width = 8,height = 8.5)

plot_grid(r_list_pre[[1]],r_list_pre[[2]],r_list_pre[[3]],nrow = 2,align = "hv")
ggsave(paste(out_dir,"rf_pre.pdf",sep = '/'),width = 8,height = 8.5)

save(model_train_2,model_classify2,tmp_data_train_f,tmp_data_pre_f,
     file = "workfile/z.rf_load.Rdata")

