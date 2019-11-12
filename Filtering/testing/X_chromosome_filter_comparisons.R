library(tidyverse)
library(VariantAnnotation)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/vcfFilterTree_functions.R")
#This script was used to help determine how the X chromosome filters should be set

vcf = readVcf("~/hpc/pmc_vanboxtel/projects/Freek_trees/test/ACHSC19_ACMSCBULK_Q50_CGQ10_SGQ10_PASS_10X_X_nonBlacklist_noEvidenceCon.vcf")

sample = "ACHSC19"
samples_bulk = samples(header(vcf))
sample_col = grep(sample, samples_bulk)

gt = geno(vcf)$GT[,sample_col]
homo_gt = ifelse(gt == "1/1", "Homozygous", "Heterozygous")
vaf = get_vaf(vcf, sample)
gq = geno(vcf)$GQ[,sample_col]
mq = info(vcf)$MQ


tb = tibble(gt, vaf, gq, mq, homo_gt)
gq_gt_fig = ggplot(tb, aes(x = homo_gt, y = gq)) +
    geom_boxplot(aes(color = homo_gt), outlier.shape = NA, fill = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 1, dotsize = 1) +
    scale_color_discrete(guide=FALSE) +
    labs(title = "GQ between heterozygous and homozygous calls on the X chromosome (pre-final SNVFI)", x = "", y = "GQ") +
    theme_bw() +
    theme(text = element_text(size = 20))

mq_gt_fig = ggplot(tb, aes(x = homo_gt, y = mq)) +
    geom_boxplot(aes(color = homo_gt), outlier.shape = NA, fill = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 1, dotsize = 0.2) +
    scale_color_discrete(guide=FALSE) +
    labs(title = "MQ between heterozygous and homozygous calls on the X chromosome (pre-final SNVFI)", x = "", y = "MQ") +
    theme_bw() +
    theme(text = element_text(size = 20))

gq_vaf_fig = ggplot(tb, aes(x = vaf, y = gq)) +
    geom_jitter(width = 0.015, height = 1) +
    labs(title = "GQ comparted to vaf on the X chromosome (pre-final SNVFI) (jitter)", x = "VAF", y = "GQ") +
    theme_bw() +
    theme(text = element_text(size = 20))

mq_vaf_fig = ggplot(tb, aes(x = vaf, y = mq)) +
    geom_jitter(width = 0.015, height = 0.5) +
    labs(title = "MQ comparted to vaf on the X chromosome (pre-final SNVFI) (jitter)", x = "VAF", y = "MQ") +
    theme_bw() +
    theme(text = element_text(size = 20))

pdf(paste0("~/surfdrive/Shared/Projects/Freek_Smallprojects/Xfiltering.pdf"))
gq_gt_fig
mq_gt_fig
gq_vaf_fig
mq_vaf_fig
dev.off()
