library(tidyverse)
library(smartcallr)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(precrec)
load('./data/gold_snps.RDA')
load('./markers.RDA')
markers$start <- markers$end # they are off by one, for reasons I can not explain!
fr <- smartcallr::compute_odds('./data/aml31R_mutect2.no.downsample.filtered.vcf',BSgenome.Hsapiens.1000genomes.hs37d5,sample_name="aml31R",prior_lod=10)
fr$chrom <- as.character(fr$seqnames)
fr$start <- as.character(fr$start)
fr$end <- as.character(fr$end)
 
joined <- fr %>% dplyr::full_join(gold_snps,by=c("chrom"="chromosome_name","start"="start","end"="stop","ref"="reference","alt"="variant"))
joined <- joined %>% dplyr::full_join(markers,by=c("chrom"="chr","start"="start","ref"="ref","alt"="var"))
joined <- joined %>% dplyr::select("tumor_capture_vaf"=TUMOR_CAPTURE_VAF,'tumor_deep_var'=`tumor var_count`,'relapse_deep_var'=`relapse var_count`,'relapse_capture_vaf'=RELAPSE_CAPTURE_VAF,tumor_vaf,relapse_vaf,normal_vaf,default_gene_name,TLOD,AF,mutect_odds,prior_odds,softFilterMatrix.clustered_events,pass_all,tlod_only,chrom,start,everything())

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
truth_table <- joined %>% dplyr::filter((pass_all | tlod_only)) %>%
  dplyr::mutate(present=!is.na(tumor_vaf) & tumor_vaf > 0 & normal_vaf ==0) %>%
  dplyr::select(present,everything())

sscurves_mutect <- precrec::evalmod(scores=truth_table$mutect_odds,labels=truth_table$present)
sscurves_prior <-  precrec::evalmod(scores=truth_table$prior_odds,labels=truth_table$present)
# 
# p1 <- ggplot(truth_table,aes(x=log10(mutect_odds),y=log10(prior_odds))) + geom_point() + xlim(c(-3,3)) + ylim(c(-3,3)) + theme_bw() + ggtitle("Comparison of odds for called variants only")
# p2 <- ggplot(joined,aes(x=log10(mutect_odds),y=log10(prior_odds))) + geom_point() + xlim(c(-3,3)) + ylim(c(-3,3)) + theme_bw() + ggtitle("Comparison of odds for all variants")
# 
# prior_density_plot <- ggplot(joined,aes(logit_prior)) + geom_density() +  geom_vline(xintercept=-6) + theme_classic() + xlab('Log prior - line at -6 is MuTect') + ggtitle('Observed distribution of\nprior values in primary aml')
# prior_table <- unique(joined$joint_prior)
# prior_table <- prior_table[!is.na(prior_table)]
# prior_fr <- data_frame(prior=prior_table)
# prior_plot <- ggplot(prior_fr,aes(prior)) + geom_density() + geom_vline(xintercept=1e-6) + theme_classic() + xlab('Log prior - line at -6 is MuTect') + ggtitle("Density of the prior over all contexts. (n=96)")
# prior_fr <- prior_fr %>% left_join(joined,by=c('prior'='joint_prior')) %>% dplyr::select(prior,context,cref,calt) %>% unique()

relapse_fr <- joined %>% dplyr::mutate(posterior_log_odds_diff=(TLOD+logit_prior)-(TLOD-6)) %>% dplyr::select(posterior_log_odds_diff,everything())
save(relapse_fr,file='relapse_frame.RDA')
#primary_fr <- joined %>% dplyr::mutate(posterior_log_odds_diff=(TLOD+logit_prior)-(TLOD-6)) %>% dplyr::select(posterior_log_odds_diff,everything())
#save(primary_fr,file='primary_frame.RDA')

ss_mutect <- as.data.frame(sscurves_mutect) %>% dplyr::filter(type=='PRC')
ss_prior <- as.data.frame(sscurves_prior) %>% dplyr::filter(type=='PRC')

prc_mutect <- ggplot(ss_mutect,aes(x=x,y=y)) + geom_line() + theme_bw() + xlab("recall") + ylab("precision") + ggtitle("PRC for MuTect2")
prc_prior <- ggplot(ss_prior,aes(x=x,y=y)) + geom_line() + theme_bw() + xlab("recall") + ylab("precision") + ggtitle("PRC for Prior")

# compute the odds at major recall inflections (RELAPSE)
ordered_mutect <- truth_table %>% dplyr::select(mutect_odds) %>% dplyr::arrange(desc(mutect_odds))
ordered_prior <- truth_table %>% dplyr::select(prior_odds) %>% dplyr::arrange(desc(prior_odds))
# Mutect2 at .75 recall
slicemutectat75 <- nrow(ss_mutect) - nrow(ss_mutect %>% dplyr::filter(x > 0.75))
mutectoddsat75 <- ordered_mutect %>% dplyr::slice(slicemutectat75)
# Prior at .75 recall
slicepriorat75 <- nrow(ss_prior) - nrow(ss_prior %>% dplyr::filter(x > 0.75))
prioroddsat75 <- ordered_prior %>% dplyr::slice(slicepriorat75)
# Mutect2 at .8 recall
slicemutectat8 <- nrow(ss_mutect) - nrow(ss_mutect %>% dplyr::filter(x > 0.80))
mutectoddsat8 <- ordered_mutect %>% dplyr::slice(slicemutectat8)
# Prior at .8 recall prior odds of 2 is between .79 and .8
slicepriorat8 <- nrow(ss_prior) - nrow(ss_prior %>% dplyr::filter(x > 0.80))
prioroddsat8 <- ordered_prior %>% dplyr::slice(slicepriorat8)

# compute the odds at major recall inflections (PRIMARY)
ordered_mutect <- truth_table %>% dplyr::select(mutect_odds) %>% dplyr::arrange(desc(mutect_odds))
ordered_prior <- truth_table %>% dplyr::select(prior_odds) %>% dplyr::arrange(desc(prior_odds))
# Mutect2 at .78 recall
slicemutectat78 <- nrow(ss_mutect) - nrow(ss_mutect %>% dplyr::filter(x > 0.78))
mutectoddsat78 <- ordered_mutect %>% dplyr::slice(slicemutectat78)
# Prior at .78 recall
slicepriorat78 <- nrow(ss_prior) - nrow(ss_prior %>% dplyr::filter(x > 0.78))
prioroddsat78 <- ordered_prior %>% dplyr::slice(slicepriorat78)
# Mutect2 at .98 recall
slicemutectat98 <- nrow(ss_mutect) - nrow(ss_mutect %>% dplyr::filter(x > 0.98))
mutectoddsat98 <- ordered_mutect %>% dplyr::slice(slicemutectat98)
# Prior at .98 recall prior odds of 2 is between .79 and .8
slicepriorat98 <- nrow(ss_prior) - nrow(ss_prior %>% dplyr::filter(x > 0.98))
prioroddsat98 <- ordered_prior %>% dplyr::slice(slicepriorat98)