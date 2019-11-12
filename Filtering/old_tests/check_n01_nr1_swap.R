library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(VariantAnnotation)
#Script to check for sample swaps

#Function from vcfFilterTree.R
vcf_gt_inf = function(vcf){
    gt = geno(vcf)$GT %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_GT"))
    gq = geno(vcf)$GQ %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_GQ")) %>% replace(., is.na(.), 0)
    dp = geno(vcf)$DP %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_DP")) %>% replace(., is.na(.), 0)
    ad = geno(vcf)$AD %>% as_tibble()
    ad_alt = apply(ad, 1:2, function(x) x[[1]][2]) %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_ADALT"))
    vaf = apply(ad, 1:2, function(x) x[[1]][2] / sum(x[[1]])) %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_vaf")) %>% replace(., is.na(.), 0)
    geno_tb = bind_cols(gt, gq, dp, ad_alt, vaf)
    return(geno_tb)
}

filter_row3 = function(row, cols, nsamples){
    gt_f = row[cols] == "0/1"
    gq_f = row[cols + nsamples] >= 99
    dp_f = row[cols + 2*nsamples] >= 20
    vaf_f = row[cols + 4*nsamples] >= 0.3 & row[cols + 4*nsamples] <= 0.7
    
    col_f = gt_f & gq_f & dp_f & vaf_f
    keep_row = sum(col_f) == length(cols)
    return(keep_row)
}

#N01
vcf_n01 = readVcf("~/hpc/pmc_vanboxtel/projects/Freek_trees/N01/shared/N01_somatic_filtered_noblacklist_MQ.vcf")
geno_tb_n01 = vcf_gt_inf(vcf_n01)
nsamples = 6

keep_rows_n01_variants = sapply(seq(nrow(geno_tb_n01)), function(i) filter_row3(geno_tb_n01[i,], cols = c(2,3,4,5,6), nsamples))
n01_variants = granges(vcf_n01[keep_rows_n01_variants])
keep_rows_n01_bm = sapply(seq(nrow(geno_tb_n01)), function(i) filter_row3(geno_tb_n01[i,], cols = c(1), nsamples))
n01_bm = granges(vcf_n01[keep_rows_n01_bm])

#NR1
vcf_nr1 = readVcf("~/hpc/pmc_vanboxtel/projects/Freek_trees/NR1/shared/NR1_somatic_filtered_noblacklist_MQ.vcf")
geno_tb_nr1 = vcf_gt_inf(vcf_nr1)
nsamples = 4

keep_rows_nr1_variants = sapply(seq(nrow(geno_tb_nr1)), function(i) filter_row3(geno_tb_nr1[i,], cols = c(2,3,4), nsamples))
nr1_variants = granges(vcf_nr1[keep_rows_nr1_variants])
keep_rows_nr1_bm = sapply(seq(nrow(geno_tb_nr1)), function(i) filter_row3(geno_tb_nr1[i,], cols = c(1), nsamples))
nr1_bm = granges(vcf_nr1[keep_rows_nr1_bm])


n01_var_bm = findOverlaps(n01_variants, n01_bm) %>% length()
nr1_var_bm = findOverlaps(nr1_variants, nr1_bm) %>% length()

n01var_nr1bm = findOverlaps(n01_variants, nr1_bm) %>% length()
nr1var_n01bm = findOverlaps(nr1_variants, n01_bm) %>% length()

