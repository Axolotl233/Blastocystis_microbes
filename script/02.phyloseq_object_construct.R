metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ## if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

library(microbiomeSeq)
library(microbiome)
out_dir = "output/metadata"
dir.create(out_dir,showWarnings = F)

data_m <- read_csv("output/subtype_class/subtype_pattern.csv")
data_m <- as.data.frame(data_m)
data_m$ST1_abs <- (data_m$ST1 * data_m$parasite_load)/100
data_m$ST2_abs <- (data_m$ST2 * data_m$parasite_load)/100
data_m$ST3_abs <- (data_m$ST3 * data_m$parasite_load)/100
data_m$parasite_load_log <- log2(data_m$parasite_load + 1)
data_m$ST1_abs_log <- log2(data_m$ST1_abs + 1)
data_m$ST2_abs_log <- log2(data_m$ST2_abs + 1)
data_m$ST3_abs_log <- log2(data_m$ST3_abs + 1)
data_m <- data_m %>%
  mutate(parasite_load_gr3 = ntile(parasite_load_log, 4))
data_m$parasite_load_gr3 <- ifelse(data_m$parasite_load_gr3 == 3,2,data_m$parasite_load_gr3)
data_m$parasite_load_gr3 <- ifelse(data_m$parasite_load_gr3 == 4,3,data_m$parasite_load_gr3)
data_m$parasite_load_gr3 <- paste("Q",data_m$parasite_load_gr3,sep = "")

#data_m %>% group_by(parasite_load_gr1,Pattern_class) %>%
#  summarise(count = n())

data_gr <- read_csv("input/Abs_abundance_6.24.csv")
data_m <- inner_join(data_m,data_gr,by = "ID")
write_csv(data_m,"input/Blastocystis_meta_final.fix.csv")

data_m2 <- read_csv(in_meta) %>% filter(class == "control")
data_m$class = "positive"
tmp_a <- (nrow(data_m)+1)
tmp_b <- (nrow(data_m)+nrow(data_m2))
data_m[tmp_a:tmp_b,"ID"] <- data_m2$ID
data_m[tmp_a:tmp_b,"sex"] <- data_m2$sex
data_m[tmp_a:tmp_b,"bmi"] <- data_m2$bmi
data_m[tmp_a:tmp_b,"species"] <- data_m2$species
data_m[tmp_a:tmp_b,"class"] <- data_m2$class

row.names(data_m) <- data_m$ID
metadatadf <- data_m
sampledata <- sample_data(metadatadf)

mphlanin <- read.csv(in_otu_ab, sep = "\t", strip.white = T, stringsAsFactors = F,row.names=1)
mphlanin <- mphlanin[,colnames(mphlanin) %in% metadatadf$ID]
mphlanin <- noise_removal(mphlanin,method = "max_cut",percent = 0,low = 0)
which(rowSums(mphlanin) == 0)
mphlanin_used_tmp <- mphlanin
#mphlanin_used_tmp<-data.frame(apply(mphlanin, 2, function(x) as.numeric(as.character(x))))
#rownames(mphlanin_used_tmp) <- row.names(mphlanin)
# only Bacteria included
mphlanin_used_tmp2 <- mphlanin_used_tmp[(grepl("k__Bacteria", rownames(mphlanin_used_tmp),perl = TRUE)),]
#filtering # keep only the species level
mphlanin_used <- mphlanin_used_tmp2[grep("s__",rownames(mphlanin_used_tmp2),perl =TRUE,invert=FALSE),]
mphlanin_used <- mphlanin_used[grep("t__",rownames(mphlanin_used),perl =TRUE,invert=TRUE),]
#mphlanin_used <- mphlanin_used[-(unlist(lapply(as.list(tmp_sp), function(x) {grep(x,row.names(mphlanin_used))}))),]
mphlanin_used[is.na(mphlanin_used)] <- 0

all(sampledata$ID %in% colnames(mphlanin_used))
mphlanin_used <- mphlanin_used[,sampledata$ID]
identical(sort(rownames(sampledata)), sort(colnames(mphlanin_used)))
phyloseqin <- metaphlanToPhyloseq(mphlanin_used, metadat = sampledata)
phyloseqin_all <- metaphlanToPhyloseq(mphlanin_used_tmp2, metadat = sampledata)

### count
mphlanin_count <- read.csv(in_otu_count, sep = "\t", strip.white = T, stringsAsFactors = F,row.names=1)
mphlanin_count <- mphlanin_count[,-1]
mphlanin_count <- mphlanin_count[,colnames(mphlanin_count) %in% metadatadf$ID]
#mphlanin_count <- noise_removal(mphlanin_count,method = "max_cut",percent = 0,low = 0)
#which(rowSums(mphlanin_count) == 0)
mphlanin_count_used_tmp<-data.frame(apply(mphlanin_count, 2, function(x) as.numeric(as.character(x))))
rownames(mphlanin_count_used_tmp) <- row.names(mphlanin_count)
mphlanin_count_used_tmp <- noise_removal(mphlanin_count_used_tmp,method = "max_cut",percent = 0,low = 0)
which(rowSums(mphlanin_count_used_tmp) == 0)
# only Bacteria included
mphlanin_count_used_tmp2 <- mphlanin_count_used_tmp[(grepl("k__Bacteria", rownames(mphlanin_count_used_tmp),perl = TRUE)),]
#filtering # keep only the species level
mphlanin_count_used <- mphlanin_count_used_tmp2[grep("s__",rownames(mphlanin_count_used_tmp2),perl =TRUE,invert=FALSE),]
mphlanin_count_used <- mphlanin_count_used[grep("t__",rownames(mphlanin_count_used),perl =TRUE,invert=TRUE),]
#mphlanin_count_used <- mphlanin_count_used[-(unlist(lapply(as.list(tmp_sp), function(x) {grep(x,row.names(mphlanin_count_used))}))),]
mphlanin_count_used[is.na(mphlanin_count_used)] <- 0

all(sampledata$ID %in% colnames(mphlanin_count_used))
mphlanin_count_used <- mphlanin_count_used[,sampledata$ID]
identical(sort(rownames(sampledata)), sort(colnames(mphlanin_count_used)))
phyloseqin_count <- metaphlanToPhyloseq(mphlanin_count_used, metadat = sampledata)
phyloseqin_count_all <- metaphlanToPhyloseq(mphlanin_count_used_tmp2, metadat = sampledata)
mphlanin_count_used <- noise_removal(mphlanin_count_used,percent = 0.1,low = 0,method = "pre_cut")
mphlanin_count_used <- mphlanin_count_used + 1
phyloseqin_count_ancombc <-  metaphlanToPhyloseq(mphlanin_count_used, metadat = sampledata)

save(phyloseqin,phyloseqin_all,phyloseqin_count,phyloseqin_count_all,sampledata,metadatadf,file="workfile/z.phyloseq.Rdata")
save(phyloseqin_count_ancombc,file="workfile/z.phyloseq.ancombc.Rdata")
rm(list = ls(pattern = "tmp"))
rm(list = ls(pattern = "mphlanin"))
rm(data_m)
