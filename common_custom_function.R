noise_removal <- function(data, percent=0.1,low = 0, method = "pre_cut",dim = 1){
  matrix <- data
  tmp_col_len <- ncol(matrix)
  if(method == "max_cut"){
    bigones <- apply(matrix,dim,function(x){max(x) > low})
    print(c(method,low))
  }
  if(method == "pre_cut"){
    bigones <- apply(matrix,dim,function(x){length(x[x > low]) >= percent*tmp_col_len})
    print(c(method,percent,low))
  }
  if(method == "mean_cut"){
    bigones <- apply(matrix,dim,function(x){sum(x,na.rm = T)/sum(!is.na(x)) > low})
    print(c(method,low))
  }
  if(method == "na_cut"){
    bigones <- apply(matrix,dim,function(x){sum(!is.na(x)) >= percent*tmp_col_len})
    print(c(method,percent))
  }
  matrix_return <- matrix[bigones,]
  return(matrix_return)
}
combine_list <- function(x,y = "none"){
  x = as.character(x)
  tmp_c=unique(x)
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    dplyr::mutate(tmp = str_c(V1,",",V2)) %>% 
    dplyr::select(tmp)
  if(y != "none"){
    tmp_a <- as.data.frame(tmp_a[grepl(y,tmp_a[,1]),])
    colnames(tmp_a) <- "tmp"
  }
  t_list <- (str_split(tmp_a$tmp,pattern = ","))
  return(t_list)
}
summary_stat <- function(data,value_col = 1,class_col = 2){
  
  zero_count <- function(x){
    c = length(which(x == 0) )
    return(c)
  }
  
  data <- as.data.frame(data)
  
  median <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=median)
  mean <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=mean)
  sd <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=sd)
  len <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=length)
  zero <- aggregate(data[,value_col], by=list(data[,class_col]), FUN=zero_count)
  
  data_stat <- data.frame(mean,median=median$x, sd=sd$x, len=len$x,zero=zero$x)
  colnames(data_stat) = c("Group","Mean","Median", "Sd", "Count","ZeroCount")
  data_stat$Se <- data_stat$Sd/sqrt(data_stat$Count)
  
  data_stat$q25 = 0
  data_stat$q75 = 0
  
  for(i_t in 1:nrow(data_stat)){
    tmp_class = data_stat[i_t,1]
    tmp_data <- data[data[,class_col] == tmp_class,value_col]
    tmp_q25 <- quantile(tmp_data,probs = 0.25) 
    tmp_q75 <- quantile(tmp_data,probs = 0.75)
    data_stat[i_t,8] = tmp_q25
    data_stat[i_t,9] = tmp_q75
  }
  return(data_stat)
}
zero_one_norm <- function(x){
  max = max(x)
  min = min(x)
  t = (x-min)/(max-min)
  return(t)
}
plot_specify_sp <- function(otu,meta,sp,color_m,sample_col = 1,group_col = 8){
  
  if(! sp %in% row.names(otu)){
    c <- print(paste(c("sp is not exist:",sp),collapse = " "))
    return(c)
  }
  
  meta <- meta[,c(sample_col,group_col)]
  colnames(meta) <- c("sample","group")
  tmp_data <- as.data.frame(t(otu[which(row.names(otu)==sp),meta$sample]))
  tmp_data$sample <- row.names(tmp_data)
  colnames(tmp_data) <- c("value","sample")
  tmp_data_p <- inner_join(meta,tmp_data,by = "sample")
  tmp_list <- combine_list(tmp_data_p$group)
  #tmp_data_p$value2 <- log10(tmp_data_p$value + 1e-6)
  res_l <- list()
  p2 <- ggplot(tmp_data)+
    geom_bar(aes(x=sample,y=value),
             stat="identity", position="stack")+
    theme_cowplot()+
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
    )+
    labs(title = sp)
  p <- ggplot(tmp_data_p,aes(x=group,y=value2,color = group))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(width = 0.2,height = 0)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = sp)+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color_m)
  res_l[[1]] <- tmp_data_p
  res_l[[2]] <- p
  res_l[[3]] <- p2
  return(res_l)
}
plot_correalation <- function(d,col_1 = 1,col_2 = 2, meth = "pearson"){
  dd <- d[,c(col_1,col_2)]
  colnames(dd) <- c("g1","g2")
  res_l <- list()
  tmp_formula <- y ~ x
  p <- ggplot(dd,aes(x=g1,y=g2)) + 
    stat_smooth(method='lm',formula = tmp_formula,color = "navy")+
    geom_point(color = "gray40",size=3)+
    stat_cor(cor.coef.name = "rho",method = meth,
             label.y.npc = 1,label.x.npc = 0.35)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid = element_blank(),
    )
  t <- cor.test(dd$g1,dd$g2)
  m <- lm(g2~g1,dd)
  res_l[[1]] <- p
  res_l[[2]] <- t
  res_l[[3]] <- m
  return(res_l)
}
list_pair_test <- function(t_data,t_list,value_col,group_col,method = "np" ){
  t_res <- list()
  t_data <- t_data[,c(group_col,value_col)]
  colnames(t_data) <- c("group","value")
  for(i in 1:length(tmp_list)){
    i = 1
    t_g <- tmp_list[[i]]
    t_data1 <- t_data[t_data$group == t_g[1],]
    t_data2 <- t_data[t_data$group == t_g[2],]
    if(method == "np"){
      t_res[[i]] <- wilcox.test(t_data1$value,t_data2$value)
    }
    if(method == "p"){
      t_res[[i]] <- t.test(t_data1$value,t_data2$value)
    }
  }
  return(t)
}


