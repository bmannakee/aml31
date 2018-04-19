library(tidyverse)
library(smartcallr)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(precrec)
load('./data/gold_snps.RDA')
load('./markers.RDA')
markers$start <- markers$end # they are off by one, for reasons I can not explain!
fr <- smartcallr::compute_odds('./data/aml31P_mutect2.filtered.vcf',BSgenome.Hsapiens.1000genomes.hs37d5,sample_name="aml31P",prior_lod=50)
fr$chrom <- as.character(fr$seqnames)
fr$start <- as.character(fr$start)
fr$end <- as.character(fr$end)
 
joined <- fr %>% dplyr::full_join(gold_snps,by=c("chrom"="chromosome_name","start"="start","end"="stop","ref"="reference","alt"="variant"))
joined <- joined %>% dplyr::full_join(markers,by=c("chrom"="chr","start"="start","ref"="ref","alt"="var"))
joined <- joined %>% dplyr::select("tumor_capture_vaf"=TUMOR_CAPTURE_VAF,'relapse_capture_vaf'='RELAPSE_CAPTURE_VAF',tumor_vaf,relapse_vaf,normal_vaf,default_gene_name,TLOD,AF,mutect_odds,prior_odds,softFilterMatrix.clustered_events,pass_all,tlod_only,chrom,start,everything())

# both_pass_and_found <- joined %>% dplyr::filter(mutect_odds > 2 & prior_odds > 2 & !is.na(RELAPSE_WGS_VAF))
# both_fail_and_found <- joined %>% dplyr::filter(mutect_odds < 2 & prior_odds < 2 & !is.na(RELAPSE_WGS_VAF))
# mutect_pass_prior_fail_and_found <- joined %>% dplyr::filter(mutect_odds >= 2 & prior_odds < 2 & !is.na(RELAPSE_WGS_VAF))
# mutect_fail_prior_pass_and_found <- joined %>% dplyr::filter(mutect_odds < 2 & prior_odds >= 2 & !is.na(RELAPSE_WGS_VAF))
# both_pass_and_not_found <- joined %>% dplyr::filter(mutect_odds > 2 & prior_odds > 2 & is.na(RELAPSE_WGS_VAF))

# called_by_them_not_seen_mutect <- joined %>% dplyr::filter(is.na(sampleNames) & as.numeric(RELAPSE_WGS_VAF) > 0.0) # n=33
# called_by_mutect_not_seen_them <- joined %>% dplyr::filter(pass_all & is.na(RELAPSE_WGS_VAF)) %>% dplyr::select(start) # n=1108. No idea why. for the first one, 1:5100397 I don't see any reads, but I might have the wrong file!
# mutect_tlod_only_called_them <- joined %>% dplyr::filter(tlod_only & as.numeric(RELAPSE_ALLDNA_VAF) > 0.0) # n=36 none elevated
# mutect_pass_all_called_them <- joined %>% dplyr::filter(pass_all & as.numeric(RELAPSE_ALLDNA_VAF) > 0.0) # n=36 none elevated
# 
# prior_elevated <- joined %>% dplyr::filter(tlod_only & prior_odds > 2) # n=0
# prior_reduced <- joined %>% dplyr::filter(pass_all & prior_odds < 2)
# prior_wtf <- joined %>% dplyr::filter(pass_all & prior_odds > 2 & is.na(RELAPSE_WGS_VAF))
truth_table <- joined %>% dplyr::filter(!is.na(tumor_vaf) & (pass_all | tlod_only)) %>%
  dplyr::mutate(present=(!is.na(tumor_capture_vaf) & tumor_capture_vaf > 0.0),
                                        not_present=(is.na(tumor_capture_vaf) | tumor_capture_vaf==0.0)) %>%
  dplyr::select(present,not_present,everything())

sscurves_mutect <- precrec::evalmod(scores=truth_table$mutect_odds,labels=truth_table$present)
sscurves_prior <-  precrec::evalmod(scores=truth_table$prior_odds,labels=truth_table$present)

p1 <- ggplot(truth_table,aes(x=log10(mutect_odds),y=log10(prior_odds))) + geom_point() + xlim(c(-3,3)) + ylim(c(-3,3)) + theme_bw() + ggtitle("Comparison of odds for called variants only")
p2 <- ggplot(joined,aes(x=log10(mutect_odds),y=log10(prior_odds))) + geom_point() + xlim(c(-3,3)) + ylim(c(-3,3)) + theme_bw() + ggtitle("Comparison of odds for all variants")

prior_density_plot <- ggplot(joined,aes(logit_prior)) + geom_density() +  geom_vline(xintercept=-6) + theme_classic() + xlab('Log prior - line at -6 is MuTect') + ggtitle('Observed distribution of\nprior values in primary aml')
prior_table <- unique(joined$joint_prior)
prior_table <- prior_table[!is.na(prior_table)]
prior_fr <- data_frame(prior=prior_table)
prior_plot <- ggplot(prior_fr,aes(prior)) + geom_density() + geom_vline(xintercept=1e-6) + theme_classic() + xlab('Log prior - line at -6 is MuTect') + ggtitle("Density of the prior over all contexts. (n=96)")
prior_fr <- prior_fr %>% left_join(joined,by=c('prior'='joint_prior')) %>% dplyr::select(prior,context,cref,calt) %>% unique()
