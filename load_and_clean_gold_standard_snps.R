library(tidyverse)

platinum_snps <- readr::read_tsv('./data/platinum_list.tsv',col_types=cols(.default = 'c'))
save(platinum_snps,file='./data/platinum_snps.RDA')
gold_snps <- readr::read_tsv('./data/gold_snps.tsv',col_types=cols(.default = 'c'))
save(gold_snps,file='./data/gold_snps.RDA')
