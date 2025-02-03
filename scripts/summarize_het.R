#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)

hetdata <- read_tsv("heterozygocity.tsv") %>% filter(MIN_HET > -100)
hist(hetdata$MIN_HET,100)
