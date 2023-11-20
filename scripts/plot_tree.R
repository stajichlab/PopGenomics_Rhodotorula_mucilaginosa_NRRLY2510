#!/usr/bin/env R
library(ggtree)
library(dplyr)
library(tidyverse)
library(treedataverse)
library(ape)
library(treeio)
library(tidytree)
library(ggtreeExtra)
library(ggstar)
library(RColorBrewer)
tree<-read.tree("strain_tree/Rmuc_v7.All.SNP.fasttree.tre")
#tree<-read.tree("strain_tree/Rmuc_v7.All.SNP.poppr.nj.tre")
meta <- read_csv("strain_metdata.csv",col_names=TRUE) %>% 
  mutate(label=Strain) %>% select(label,SimpleEnv,culture_collection) %>% 
  filter(!is.na(SimpleEnv))

unique(meta$SimpleEnv)
unique(meta$culture_collection)
p <- ggtree(tree, layout="circular", size=0.3, branch.length = "none")
p

colourCount = length(unique(meta$SimpleEnv))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p2<-p + geom_fruit(
  data=meta,
  geom=geom_star,
  mapping=aes(y=label, fill=SimpleEnv, starshape=culture_collection),
  starstroke=0.2
) +   scale_fill_manual(values = getPalette(colourCount))

ggsave("phylogram.pdf",p2)

p <- ggtree(tree, size=0.5,layout="circular") +
    geom_tippoint(aes(color=SimpleEnv), 
                size=1.5,
                show.legend=FALSE) +
  scale_fill_manual(values = getPalette(colourCount))
p

p2<-p + geom_fruit(
  data=meta,
  geom=geom_star,
  mapping=aes(y=label, fill=SimpleEnv, group=culture_collection),
  outlier.size=0.5,
  outlier.stroke=0.08,
  outlier.shape=21,
  starstroke=0.2
) +   scale_fill_manual(values = getPalette(colourCount))
  
