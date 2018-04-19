library(tidyverse)

platinum_snps <- readr::read_tsv('./data/platinum_list.tsv',col_types=cols(.default = 'c'))
save(platinum_snps,file='./data/platinum_snps.RDA')
gold_snps <- readr::read_tsv('./data/gold_snps.tsv',col_types=cols(.default = 'c'))
save(gold_snps,file='./data/gold_snps.RDA')
# Load and clean complete marker set and readcounts
markers <- readr::read_tsv('./data/Supplemental_Dataset_2-SomaticSnvsTargetedForValidation-ReadCounts.tsv',col_types=cols(.default= 'c')) %>%
  mutate(tumor_vaf=as.numeric(`tumor vaf`),
         relapse_vaf=as.numeric(`relapse vaf`),
         normal_vaf=as.numeric(`normal vaf`),
         normal_depth=as.numeric(`normal depth`),
         tumor_depth=as.numeric(`tumor depth`),
         relapse_depth=as.numeric(`relapse depth`)) %>%
  dplyr::select(tumor_vaf,relapse_vaf,normal_vaf,normal_depth,tumor_depth,relapse_depth,everything())
save(markers,file='markers.RDA')
