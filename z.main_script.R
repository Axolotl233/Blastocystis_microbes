rm(list=ls())

set.seed(12345)
dir.create("output")
dir.create("workfile")

source("common_custom_function.R")
library(tidyverse)
library(magrittr)
library(phyloseq)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(ggpmisc)
library(ggview)
library(RColorBrewer)
library(Maaslin2)
#library(ggtext)
#library(ggnewscale)
#library(MoMAColors)
#library(ggsci)

#=====> global variant

in_meta <- "input/Blastocystis_meta_2310.csv"
in_meta_fix <- "input/Blastocystis_meta_2310.fix.csv"
in_amp <- "input/Blastocystis_depth.csv"
in_amp_taxon <- "input/subtype_taxon_depth.csv"
in_tree <- "input/filter.treefile2"
in_otu_ab <- "input/metaphlan4.RA.txt"
in_otu_count <- "input/metaphlan4.readscount.txt"
subtype <- c("ST1","ST2","ST3","ST5")
pattern_type <- c("ST1_dom","Balance","ST3_dom","ST2_dom","Unclass")
sp <- c("Macaca fascicularis","Macaca mulatta","mfas","mmul")

color0 <- c("#94ad90","#447593","#be6f7c","#66446f","grey70")
color0_1 <- c("#94ad90","#be6f7c","#447593","grey70","#d4d7c0","#66446f")
color1 <- c("#2A7185","#725CA5")
color2 <- c("#dbc0af","#eac47c","#a5c2cd","#9fba95","#9b8e95")
color_bar <- c("#65387d","#7c9fb0","#9163b6","#be5168","#447c69","#e9d78e","#e0598b","#5698c4",
               "#e2975d","#c94a53","#51574a","#f19670","#993767","#9abf88","#a34974","#4e2472",
               "#8e8c6d","#e4bf80","#e279a3")
shape0 <- c(16,15)
