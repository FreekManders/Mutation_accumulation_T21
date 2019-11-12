library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(gplots)
library(MutationalPatterns)
library(ComplexHeatmap)
library(ggeffects)

set.seed(42)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/vcfFilterTree_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/Do_mutationalpatterns.R")
genome = "BSgenome.Hsapiens.UCSC.hg19"

out_dir_base = "~/hpc/pmc_vanboxtel/projects/Freek_trees/"
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

cosmic_sig_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_working_signatures_incl_hspc.txt"
signatures = read.table(cosmic_sig_fname, sep = "\t", header = T)
signatures = as.matrix(signatures[,-c(1,2)])

fetuses_tb = read_tsv(paste0(out_dir_base, "fetuses_inclewart.txt"), col_types = cols(.default = "c"))
fetuses = split(fetuses_tb, seq(1,nrow(fetuses_tb)))
fetuses = lapply(fetuses, unlist)
overview_samples = create_overview_samples(fetuses, out_dir_base)

chromosomes = paste0("chr", c(1:22, "X"))
overview_high = overview_samples[overview_samples$sample %in% c("OS1LIHSPCC1D23", "OS1LIHSPCC2B21", "OS1LIHSPCC2D6"),]
gr_list = apply(overview_high, 1, function(x) get_all_snvs(x, chromosomes))
grl = GRangesList(gr_list)
names(grl) = overview_high$sample
grl = keepSeqlevels(grl, chromosomes, pruning.mode = "fine")
mut_mat = mut_matrix(grl, ref_genome)
mut_mat_pool = pool_mut_mat(mut_mat, rep("High mutational load", 3))
mut_mat_vaf = get_vafsplit(overview_high$unique_vcf, samples = overview_high$sample, chromosomes)
mut_mat_vaf_pool = pool_mut_mat(mut_mat_vaf, rep(c("low vaf", "high vaf"), 3))

refit_per_sample = function(mut_mat, max_delta = 0.05){
    figs_l = vector("list", ncol(mut_mat))
    for (i in 1:ncol(mut_mat)){
        mut_mat_sample = mut_mat[,i, drop = F]
        refit_out = fit_to_signatures_selection(mut_mat_sample, signatures)
        sim_decay_fig = refit_out$sim_decay_fig #Show how the cosine similarity is decreased as signatures are removed
        refit_heatmap_fig = plot_contribution_heatmap(refit_out$fit_res$contribution, cluster_samples = F) #Show a heatmap of the refitting.
        
        contri_boots = fit_to_signatures_bootstrapped(mut_mat_sample, signatures, method = "selection", max_delta = max_delta)
        nr_selected_fig = plot_nr_signatures_selected(contri_boots) #Plot how many signatures were selected during the different bootstrapping iterations.
        fraction_fig = plot_fraction_contri(contri_boots) #Plot how often each signature was selected.
        contri_boots_fig = plot_boots_contri(contri_boots) #Plot the contribution of bootstrapped selected signatures
        ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat_sample, contri_boots, signatures) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
        sig_cor_figs = plot_sig_contri_cor(contri_boots) #Plot the signature contribution correlations between the signatures.
        figs = list(sim_decay_fig, refit_heatmap_fig, nr_selected_fig, fraction_fig, contri_boots_fig, ori_vs_rec_fig, sig_cor_figs)
        figs_l[[i]] = figs
    }
    return(figs_l)
}

figs_l = refit_per_sample(mut_mat)
figs_l_pool = refit_per_sample(mut_mat_pool)
figs_l_vaf = refit_per_sample(mut_mat_vaf)
figs_l_vaf_pool = refit_per_sample(mut_mat_vaf_pool)



pdf("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/OS1_refitting_individually.pdf", width = 10)
figs_l
figs_l_pool
figs_l_vaf
figs_l_vaf_pool
dev.off()

chromosomes = paste0("chr", c(1:22, "X"))
nr_muts_tb_raw = read_tsv(paste0(out_dir_base, "count_muts.txt"))
nr_muts_tb_raw = left_join(nr_muts_tb_raw, overview_samples[,c("fetus", "sample", "celltype", "surveyed_region_pct", "origin", "gender", "trisomy", "age_year", "age_weeks")], by = c("fetus_name" = "fetus", "sample")) %>% 
    filter(celltype != "Liver")
samples_high = nr_muts_tb_raw %>% dplyr::filter(nr_total_snv > 50) %>% pull(sample)
overview_high = overview_samples %>% dplyr::filter(sample %in% samples_high)
gr_list = apply(overview_high, 1, function(x) get_all_snvs(x, chromosomes))
grl = GRangesList(gr_list)
names(grl) = overview_high$sample
grl = keepSeqlevels(grl, chromosomes, pruning.mode = "fine")
mut_mat = mut_matrix(grl, ref_genome)
figs_l = refit_per_sample(mut_mat)
pdf("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/high_refitting_individually.pdf", width = 10)
figs_l
dev.off()

