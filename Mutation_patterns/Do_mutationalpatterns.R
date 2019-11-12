library(plyr)
library(tidyverse)
library(MutationalPatterns)
library(Homo.sapiens)
library(gridExtra)
library(biomaRt)
library(ggpubr)
library(PCDimension)
library(ggbiplot)
library(ComplexHeatmap)
library(ggdendro)
library(cowplot)
library(Rtsne)
library(broom)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/Run_mut_signature_analysis/Better_refitting.R")

ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


library(reshape2)
# plot_inner_cosheat = function(signatures, col_order){
#     cos_sim_matrix = cos_sim_matrix(signatures, signatures)
#     sample_order = col_order
#     Cosine.sim = NULL
#     Signature = NULL
#     Sample = NULL
#     x = NULL
#     y = NULL
#     xend = NULL
#     yend = NULL
#     cos_sim_matrix.m = melt(cos_sim_matrix)
#     colnames(cos_sim_matrix.m) = c("Sample", "Signature", "Cosine.sim")
#     cos_sim_matrix.m$Signature = factor(cos_sim_matrix.m$Signature, 
#                                         levels = col_order)
#     cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample, 
#                                      levels = sample_order)
#     heatmap = ggplot(cos_sim_matrix.m, aes(x = Signature, y = Sample, 
#                                            fill = Cosine.sim, order = Sample)) + geom_tile(color = "white") + 
#         scale_fill_distiller(palette = "YlGnBu", direction = 1, 
#                              name = "Cosine \nsimilarity", limits = c(0, 1.000001)) + 
#         theme_bw() + theme(axis.text.x = element_text(angle = 90, 
#                                                       hjust = 1)) + labs(x = NULL, y = NULL)
#     return(heatmap)
# }
plot_inner_cosheat = function(mat){
    cos_sim_matrix = cos_sim_matrix(mat, mat)
    
    clust = as.dist(1-cos_sim_matrix) %>% hclust()
    #clust = dist(t(mat)) %>% hclust()
    order = colnames(mat)[clust$order]
    
    
    dhc = as.dendrogram(clust)
    ddata = dendro_data(dhc, type = "rectangle")
    sig_dendrogram = ggplot(segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
        theme_dendro()
    
    sample_order = order
    Cosine.sim = NULL
    Signature = NULL
    Sample = NULL
    x = NULL
    y = NULL
    xend = NULL
    yend = NULL
    cos_sim_matrix.m = cos_sim_matrix %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather(-Sample, key = "Sample2", value = "Cosine.sim")
    cos_sim_matrix.m$Sample2 = factor(cos_sim_matrix.m$Sample2, levels = order)
    cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample, levels = sample_order)
    heatmap = ggplot(cos_sim_matrix.m, aes(x = Sample2, y = Sample, 
                                           fill = Cosine.sim, order = Sample)) + geom_tile(color = "white") + 
        scale_fill_distiller(palette = "YlGnBu", direction = 1, 
                             name = "Cosine \nsimilarity", limits = c(0, 1.000001)) + 
        theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                      hjust = 1)) + labs(x = NULL, y = NULL)
    
    heatmap = plot_grid(sig_dendrogram, heatmap, align = "v", rel_heights = c(0.2, 1), axis = "lr", nrow = 2, scale = c(1,1))
    return(heatmap)
}

binomial_test = function(p, n, x){
    expected = p * n
    if (x < expected) {
        pval = pbinom(x, n, p, lower.tail = TRUE)
        effect = "depletion"
    }
    else {
        pval = pbinom(x - 1, n, p, lower.tail = FALSE)
        effect = "enrichment"
    }
    pval = pval * 2 #make test two sided.
    #if (pval < 0.05) 
    #    significant = "*"
    #else significant = ""
    res = data.frame(effect, pval)
    return(res)
}

enrichment_depletion_test = function (x, by = c()){
    if (length(by) > 0) {
        x$by = by
        res2 = stats::aggregate(cbind(n_muts, surveyed_length, 
                                      surveyed_region_length, observed) ~ by + region, 
                                data = x, sum)
    }
    else {
        res2 = x
        res2$by = res2$sample
        res2 = res2[, c(9, 1, 3, 4, 6, 8)]
    }
    res2$prob = res2$n_muts/res2$surveyed_length
    res2$expected = res2$prob * res2$surveyed_region_length
    res3 = data.frame()
    for (i in 1:nrow(res2)) {
        x = res2[i, ]
        res3 = rbind(res3, binomial_test(x$prob, x$surveyed_region_length, 
                                         x$observed))
    }
    df = cbind(res2, res3)
    df$fdr = df$pval %>% p.adjust(method = "fdr")
    df$significant = ifelse(df$fdr < 0.1, "*", "")
    return(df)
}

plot_enrichment_depletion = function (df) 
{
    df2 = melt(df[, c(1, 2, 6, 8)], id = c("by", "region"))
    value = NULL
    variable = NULL
    observed = NULL
    expected = NULL
    significant = NULL
    plot1 = ggplot(df2, aes(x = by, y = value, fill = by, group = variable, 
                            alpha = variable)) + geom_bar(colour = "black", stat = "identity", 
                                                          position = position_dodge()) + facet_grid(~region) + 
        theme_classic() + theme(axis.ticks = element_blank(), axis.text.x = element_blank(), 
                           legend.title = element_blank()) + xlab("") + ylab("No. mutations") + 
        scale_x_discrete(breaks = NULL)
    max = round(max(abs(log2((df$observed + 0.1)/(df$expected + 
                                                      0.1)))), digits = 1) + 0.1
    plot2 = ggplot(data = df, aes(x = by, y = log2((observed + 
                                                        0.1)/(expected + 0.1)), fill = by)) + geom_bar(colour = "black", 
                                                                                                       stat = "identity", position = position_dodge()) + scale_y_continuous(limits = c(-max, 
                                                                                                                                                                                       max)) + geom_text(aes(x = by, y = log2((observed + 0.1)/(expected + 
                                                                                                                                                                                                                                                    0.1)), label = significant, vjust = ifelse(sign(log2((observed + 
                                                                                                                                                                                                                                                                                                              0.1)/(expected + 0.1))) > 0, 0.5, 1)), size = 8, position = position_dodge(width = 1)) + 
        facet_grid(~region) + theme_classic() + theme(axis.ticks = element_blank(), 
                                                 axis.text.x = element_blank(), legend.title = element_blank()) + 
        xlab("") + ylab("log2(observed/expected)") + scale_x_discrete(breaks = NULL)
    output <- cowplot::plot_grid(plot1, plot2, ncol = 1, nrow = 2, 
                                 rel_heights = c(2, 1.2))
    return(output)
}


strand_bias_test = function (strand_occurrences){
    group = NULL
    type = NULL
    strand = NULL
    variable = NULL
    df_strand = reshape2::dcast(melt(strand_occurrences), group + 
                                    type ~ strand, sum, subset = plyr::.(variable == "no_mutations"))
    df_strand$total = df_strand[, 3] + df_strand[, 4]
    df_strand$ratio = df_strand[, 3]/df_strand[, 4]
    df_strand$p_poisson = apply(df_strand, 1, function(x) poisson.test(c(as.numeric(x[3]), 
                                                                         as.numeric(x[4])), r = 1)$p.value)
    df_strand$fdr = df_strand$p_poisson %>% p.adjust(method = "fdr")
    df_strand$significant = ifelse(df_strand$fdr < 0.1, "*", " ")
    return(df_strand)
}

enrichment_depletion_between_test = function(distr, by = experiment_order){
    distr2 = cbind(distr, experiment_order)
    distr_group_betweenexps = distr2 %>% 
        dplyr::select(region, experiment_order, expected, observed) %>% 
        dplyr::group_by(region, experiment_order) %>%
        dplyr::summarise(expected = sum(expected), observed = sum(observed)) %>% 
        dplyr::group_map(~ tidy(chisq.test(cbind(.$observed, .$expected), simulate.p.value = T))) %>%
        dplyr::ungroup() %>% 
        dplyr::mutate(fdr = p.adjust(p.value, method = "fdr"))
    return(distr_group_betweenexps)
}

COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")

plot_strand_bias = function(strand_bias, colors){
    var_names = colnames(strand_bias)[3:4]
    if (missing(colors)) 
        colors = COLORS6
    log2_ratio = log2((strand_bias[, 3] + 0.1)/(strand_bias[, 
                                                            4] + 0.1))
    max = round(max(abs(log2_ratio)), digits = 1) + 0.1
    
    pos_stars = abs(log2((strand_bias[, 3])/(strand_bias[, 4] + 0.1)))
    max_pos_star = round(max(pos_stars[is.finite(pos_stars)]), digits = 1) + 0.1
    if(max < max_pos_star){
        max = max_pos_star
    }
    
    type = NULL
    significant = NULL
    label2 = log2(strand_bias$ratio)
    select = which(is.finite(label2))
    label2[select] = " "
    plot = ggplot(strand_bias, aes(x = type, y = log2((strand_bias[, 3] + 0.1)/(strand_bias[, 4] + 0.1)), fill = type)) + 
        scale_fill_manual(values = COLORS6) + 
        geom_bar(colour = "black", stat = "identity", position = "identity") + 
        scale_y_continuous(limits = c(-max,max)) + 
        geom_text(aes(x = type, y = log2((strand_bias[,3])/(strand_bias[, 4] + 0.1)), ymax = log2((strand_bias[,3])/(strand_bias[, 4] + 0.1)), label = significant, vjust = ifelse(sign(log2((strand_bias[, 3])/(strand_bias[, 4] + 0.1))) > 0, 0.5, 1)), size = 8, position = ggplot2::position_dodge(width = 1)) + 
        facet_grid(. ~ group) + theme_bw() + 
        theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.title = element_blank()) + 
        xlab("") + 
        ylab(paste("log2(", var_names[1], "/", var_names[2],")", sep = "")) + 
        scale_x_discrete(breaks = NULL)
    return(plot)
}

#Faster mut_matrix function
C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
TRIPLETS_96 = paste(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3), sep = "")
mut_96_occurrences2 = function(type_context, gr_sizes){
    cats = tibble("categories" = factor(TRIPLETS_96, levels = TRIPLETS_96))
    full_context = paste0(substr(type_context$context, 1, 1), "[", type_context$types, "]", substr(type_context$context, 3, 3)) %>% factor(levels = TRIPLETS_96)
    
    if (is.null(names(gr_sizes))){
        n = length(gr_sizes)
        names(gr_sizes) = seq(n)
    }
    
    sample_vector = rep(names(gr_sizes), gr_sizes) %>% factor(levels = names(gr_sizes))
    
    #Count the mutations per type and per sample
    counts = tibble("categories" = full_context, "sample" = sample_vector) %>% dplyr::filter(!is.na(categories)) %>% dplyr::group_by(categories, sample) %>% dplyr::summarise(count = n())
    counts = left_join(cats, counts, by = "categories")
    
    #Transform the data into a mutation matrix
    counts = spread(counts, key = sample, value = count, fill = 0)
    unnecesary_cols = which(colnames(counts) == "<NA>")
    mut_mat = as.matrix(counts[,-c(1, unnecesary_cols)])
    rownames(mut_mat) = counts$categories
    return(mut_mat)
    
    # cats = tibble("categories" = TRIPLETS_96)
    # full_context = paste0(substr(type_context$context, 1, 1), "[", type_context$types, "]", substr(type_context$context, 3, 3))
    #counts = table(full_context) %>% enframe(name = "categories", value = "count") %>% as.data.frame()
    #counts = left_join(cats, counts, by = "categories")
    #counts[is.na(counts$count), "count"] = 0
    #vector = counts$count
    #names(vector) = counts$categories
    #return(vector)
}

mut_matrix2 = function (grl, ref_genome) {
    if (class(grl)[[1]] == "CompressedGRangesList"){
        gr_sizes = elementNROWS(grl)
        gr = unlist(grl)
    } else if (class(grl)[[1]] == "GRanges"){
        gr = grl
        gr_sizes = length(gr)
        names(gr_sizes) = "My_sample"
    } else{
        stop("This function requires either a GRanges or a CompressedGRangesList object as its first argument.")
    }
    type_context = type_context(gr, ref_genome)
    mut_mat = mut_96_occurrences2(type_context, gr_sizes)
    return(mut_mat)
}

refit_per_sample = function(mut_mat, signatures, max_delta = 0.05){
    figs_l = vector("list", ncol(mut_mat))
    for (i in 1:ncol(mut_mat)){
        mut_mat_sample = mut_mat[,i, drop = F]
        refit_out = fit_to_signatures_selection(mut_mat_sample, signatures)
        sim_decay_fig = refit_out$sim_decay_fig #Show how the cosine similarity is decreased as signatures are removed
        refit_heatmap_fig = plot_contribution_heatmap(refit_out$fit_res$contribution, cluster_samples = F) #Show a heatmap of the refitting.
        
        #plot similarity of signatures with mut_mat
        clust = dist(t(signatures)) %>% hclust()
        sig_order = colnames(signatures)[clust$order]
        dhc = as.dendrogram(clust)
        ddata = dendro_data(dhc, type = "rectangle")
        sig_dendrogram = ggplot(segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
            theme_dendro()
        cos_sim_samples_signatures = cos_sim_matrix(mut_mat_sample, signatures)
        sig_cosim_heat_fig = plot_cosine_heatmap(cos_sim_samples_signatures, cluster_rows = F, col_order = sig_order)
        sig_cosim_heat_fig = plot_grid(sig_dendrogram, sig_cosim_heat_fig, align = "v", rel_heights = c(0.2, 1), axis = "lr", nrow = 2, scale = c(1.062,1))
        
        
        contri_boots = fit_to_signatures_bootstrapped(mut_mat_sample, signatures, method = "selection", max_delta = max_delta)
        nr_selected_fig = plot_nr_signatures_selected(contri_boots) #Plot how many signatures were selected during the different bootstrapping iterations.
        fraction_fig = plot_fraction_contri(contri_boots) #Plot how often each signature was selected.
        contri_boots_fig = plot_boots_contri(contri_boots) #Plot the contribution of bootstrapped selected signatures
        ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat_sample, contri_boots, signatures) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
        sig_cor_figs = plot_sig_contri_cor(contri_boots) #Plot the signature contribution correlations between the signatures.
        figs = list(sig_cosim_heat_fig, sim_decay_fig, refit_heatmap_fig, nr_selected_fig, fraction_fig, contri_boots_fig, ori_vs_rec_fig, sig_cor_figs)
        figs_l[[i]] = figs
    }
    return(figs_l)
}
#Script to perform mutational patterns analyses

####__________Helper function to write a config file and run rurika's mutsigpipe analysis_______####
do_mutsigpipe = function(my_dir, exp_name, type_input = "mut_mat", mutsig_ref = "hg19", mutsig_tresh = 0.9, mutsig = "~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/MutSigPipe/MutSigPipe_v8.R"){
    
    out_dir_mutsigpipe = paste0(my_dir, exp_name, "_mutsig8")
    if (dir.exists(out_dir_mutsigpipe)){
        print(paste0("Mut sig pip was already run for experiment: ", exp_name))
        return(0)
    }
    
    #Create config file
    config_fname = paste0(my_dir, exp_name, "_config_v8.txt")
    
    if (type_input == "mut_mat"){
        invar = "inMutMat"
        input_path = paste0(my_dir, exp_name,  "_mut_mat.txt")
        
    } else if (type_input == "vcf"){
        invar = "inmy_dir"
        input_path = my_dir
        
    }
    line_input = paste0("input_fmt <- '", type_input, "'\n")
    line_input_path = paste0(invar, " <- '", input_path, "'\n")
    line_outmy_dir = paste0("wkDIR <- '", my_dir, "'\n")
    line_refgenome = paste0("refGenome <- '", mutsig_ref, "'\n")
    line_treshold = paste0("threshold <- '", mutsig_tresh, "'\n")
    
    lines_all = c(line_input, line_input_path, line_outmy_dir, line_refgenome, line_treshold)
    filecon = file(config_fname)
    writeLines(lines_all, filecon)
    close(filecon)
    
    #Run mutsigpipe
    system(paste("Rscript", mutsig, config_fname, sep = " "))
    file.rename(paste0(my_dir, "MutSigPipe_v8"), out_dir_mutsigpipe)
    print("Ran mutsigpipe")
}

#Helper function to read and reduce bed files.
get_bed_as_gr = function(bed_fname, chromosomes){
    bed = read_tsv(bed_fname, col_names = F, col_types = cols(X1 = "c"))
    colnames(bed)[1:3] = c("CHROM", "START", "END")
    bed_gr = makeGRangesFromDataFrame(bed, starts.in.df.are.0based = T) %>% sort() %>% GenomicRanges::reduce()
    seqlevelsStyle(bed_gr) = "UCSC"
    seqlevels(bed_gr, pruning.mode = "coarse") = chromosomes
    return(bed_gr)
}


####__________________Main functions__________________####
do_mutationalpatterns = function(grl, exp_name, factors, out_dir_mut, chromosomes){
    cosmic_sig_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_working_signatures.txt"
    slopes_bed_fname = "~/surfdrive/Shared/Boxtel_General/Data/Replication_timing/ENCODE/blood/B_lymphocytes_direction.bed"
    
    ####____________________add exp info to the grl. Also create a grl for the groups and a gr with everything.___________####
    
    
    grl = keepSeqlevels(grl, chromosomes, pruning.mode = "fine")
    for (i in seq(1, nrow(factors))){
        exp = factors$experiment[i]
        grl[[i]]$experiment = exp
    }
    gr_all = unlist(grl) %>% sort()
    grl_group = split(gr_all, gr_all$experiment)
    
    factors$experiment = as.factor(factors$experiment)
    
    
    
    ####________________________________spectra______________________________####
    #Per sample
    type_occur = mut_type_occurrences(grl, ref_genome)
    type_occur_ungrouped = type_occur %>% dplyr::group_by("experiment" = factors$experiment)
    write_tsv(type_occur_ungrouped, paste0(out_dir_mut, exp_name,  "_type_occur_ungrouped.txt"))
    type_occur_all = type_occur %>% dplyr::group_by("experiment" = factors$experiment) %>% dplyr::summarise_all(sum)
    write_tsv(type_occur_all, paste0(out_dir_mut, exp_name,  "_type_occur.txt"))
    
    spectrum_fig = plot_spectrum(type_occur, by = names(grl), CT = T) + theme_classic()
    
    mut_mat = mut_matrix(grl, ref_genome)
    write.table(mut_mat, paste0(out_dir_mut, exp_name,  "_mut_mat.txt"), quote = F, sep = "\t")
    
    index_list = seq(1, ncol(mut_mat)) %>% split(., ceiling(seq_along(.)/8))
    spectrum96_figs = lapply(index_list, function(x){
        mut_mat_sub = mut_mat[, x, drop = F]
        spectrum96_fig = plot_96_profile(mut_mat_sub) + theme_classic()
        return(spectrum96_fig)
    })
    
    
    
    #Per group
    spectrum_group_fig = plot_spectrum(type_occur, by = factors$experiment, CT = T) + theme_classic()
    
    #Calculate if 6 profile is different between the two groups
    if (length(unique(factors$experiment)) > 1){
        type_occur_noct = type_occur %>% dplyr::select(-`C>T`)#Use the C>T cols that are split over cpg instead.
        type_occur_norm = apply(type_occur_noct, 1, function(i) i/sum(i)) %>% t() %>% as_tibble()
        pvals = apply(type_occur_norm, 2, function(i){#Compare each feature between groups
            wil_res = wilcox.test(i ~ factors$experiment, exact = F)
            pval = wil_res$p.value
            return(pval)
            })
        fdrs = pvals %>% p.adjust(method = "fdr")

        type_occur_group = type_occur %>% group_by(factors$experiment) %>% summarise_all(list(~sum))
        type_occur_group_m = type_occur_group[,-c(1,4)] %>% as.matrix()
        rownames(type_occur_group_m) = type_occur_group$`factors$experiment`
        chisq_res = chisq.test(type_occur_group_m, simulate.p.value = T)#perform chi-square to see if there is any difference between the groups
        
        vals = c(fdrs, "chisq_pval" = chisq_res$p.value) %>% enframe(name = "Mutation_class", value = "fdr/pval")
        write_tsv(vals, paste0(out_dir_mut, exp_name, "_spectra6_pvals.txt"))
    }
    
    #96 profile per group
    mut_mat_group = mut_mat %>% t(.) %>% as_tibble() %>% dplyr::mutate(factor = factors$experiment) %>% group_by(factor) %>% summarise_all(funs(sum)) %>% dplyr::select(-factor) %>% t(.)
    colnames(mut_mat_group) = levels(factors$experiment)
    write.table(mut_mat_group, paste0(out_dir_mut, exp_name,  "_combined_mut_mat.txt"), quote = F, sep = "\t")
    spectrum96_group_fig = plot_96_profile(mut_mat_group, condensed = T) + theme_classic()
    
    if (length(unique(factors$experiment)) > 1){
        mut_mat_norm = mut_mat %>% apply(2, function(i) i/sum(i))
        pvals = apply(mut_mat, 1, function(i){
            wil_res = wilcox.test(i ~ factors$experiment, exact = F)
            pval = wil_res$p.value
            return(pval)
        })
        fdrs = pvals %>% p.adjust(method = "fdr")
        
        mut_mat_group_exist_levels = mut_mat_group[rowSums(mut_mat_group) != 0,] #Remove groups that have no mutations, because the expected number of mutations here is 0, which breaks chisq.
        chisq_res = chisq.test(mut_mat_group_exist_levels, simulate.p.value = T)
        
        vals = c(fdrs, "chisq_pval" = chisq_res$p.value) %>% enframe(name = "Mutation_class", value = "fdr/pval")
        write_tsv(vals, paste0(out_dir_mut, exp_name, "_spectra96_pvals.txt"))
    }
    
    spectrum96_group_fig2 = plot_compare_profiles(mut_mat_group[,1], mut_mat_group[,2], profile_names = colnames(mut_mat_group)[1:2], condensed = T) +
        facet_grid(variable ~ substitution) + theme_classic()
    
    print("Plotted spectra")
    
    
    ###PCA of spectra
    signatures = read.table(cosmic_sig_fname, sep = "\t", header = T)
    signatures = as.matrix(signatures[,-c(1,2)])
    rownames(signatures) = rownames(mut_mat)
    
    #Get features with enough mutations
    good_features = rowSums(mut_mat) > 5
    mut_mat_goodfeatures = mut_mat[good_features,]
    signatures_goodfeatures = signatures[good_features,]
    
    #Perform pca and determine number of significant pcs.
    prcom_96 = mut_mat_goodfeatures %>% scale() %>% t() %>% prcomp(scale = T, center = T)
    prcom_var = prcom_96$sdev * prcom_96$sdev
    total_variance = sum(prcom_var)
    nr_sig_pcs = bsDimension(prcom_var)
    nr_pcs = prcom_96$sdev %>% length()
    expected_var_pc1 = brokenStick(1, nr_pcs)
    expected_var_pc2 = brokenStick(2, nr_pcs)
    
    #cosine similarities between pcs and signatures
    pc_sig_cosim = cos_sim_matrix(prcom_96$rotation, signatures_goodfeatures)
    pc_sig_cosim_long = as.data.frame(pc_sig_cosim) %>% rownames_to_column("PC") %>% gather(-PC, key = "Signature", value = "cossim")
    pc_sig_cosim_long$PC = factor(pc_sig_cosim_long$PC, levels = rev(rownames(pc_sig_cosim)))
    pc_sig_cosim_long$Signature = factor(pc_sig_cosim_long$Signature, levels = colnames(pc_sig_cosim))
    
    pc_sig_fig = ggplot(pc_sig_cosim_long, aes(x = Signature, fill = cossim, y = PC), order = NULL) +
        geom_tile(color = "white") +
        scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Difference") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = NULL, y = NULL)
    
    
    pc1_positive = prcom_96$rotation[,1] > 0
    pc_sig_cosim = cos_sim_matrix(prcom_96$rotation[pc1_positive, 1, drop = F], signatures_goodfeatures[pc1_positive,])
    pc_sig_fig_pos = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
        labs(title = "Positive PC1")
    pc_sig_cosim = cos_sim_matrix(abs(prcom_96$rotation[!pc1_positive, 1, drop = F]), signatures_goodfeatures[!pc1_positive,])
    pc_sig_fig_neg = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
        labs(title = "Negative PC1")

    pc2_positive = prcom_96$rotation[,2] > 0
    pc_sig_cosim = cos_sim_matrix(prcom_96$rotation[pc2_positive, 2, drop = F], signatures_goodfeatures[pc2_positive,])
    pc_sig2_fig_pos = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
        labs(title = "Positive PC2")
    pc_sig_cosim = cos_sim_matrix(abs(prcom_96$rotation[!pc2_positive, 2, drop = F]), signatures_goodfeatures[!pc2_positive,])
    pc_sig2_fig_neg = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
        labs(title = "Negative PC2")
    

    #Plot pca
    pca_axis_fig = ggbiplot(prcom_96, obs.scale = 1, ellipse = T, var.axes = T,
             groups = as.vector(factors$experiment), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic() +
        annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
    
    #Plot pca no axis
    pca_fig = ggbiplot(prcom_96, obs.scale = 1, ellipse = T, var.axes = F,
             groups = as.vector(factors$experiment), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic() +
        annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
    
    #Plot pca with fetus as group
    if(as.vector(factors$fetus) %>% table() %>% max() >= 3){
        pca_fetus_fig = ggbiplot(prcom_96, obs.scale = 1, ellipse = T, var.axes = F,
                           groups = as.vector(factors$fetus), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
            theme_classic() +
            annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
    } else{
        pca_fetus_fig = ggbiplot(prcom_96, obs.scale = 1, ellipse = T, var.axes = F, choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
            theme_classic() +
            annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
    }
    
    #Plot pca 2 and 3
    pca_23_fig = ggbiplot(prcom_96, obs.scale = 1, ellipse = T, var.axes = F,
                       groups = as.vector(factors$experiment), choices = c(2,3) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic()

    #Plot pca 3 and 4
    pca_34_fig = ggbiplot(prcom_96, obs.scale = 1, ellipse = T, var.axes = F,
                          groups = as.vector(factors$experiment), choices = c(3,4) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic()
    
    #overlay signature 1 to existing pca model
    sigs_pca_input = signatures_goodfeatures[, "SBS1", drop = F] %>% scale() %>% t()
    predicted = predict(prcom_96, newdata = sigs_pca_input)
    prcom_withsig1 = prcom_96
    prcom_withsig1$x = rbind(prcom_withsig1$x, predicted)
    pca_sig1_fig = ggbiplot(prcom_withsig1, obs.scale = 1, ellipse = T, var.axes = F,
             groups = c(as.vector(factors$experiment), rownames(sigs_pca_input)), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic() +
        annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
    
    
    #Overlay all signatures to existing pca model
    nrsigs = ncol(signatures_goodfeatures)
    sig_l = seq(1, nrsigs) %>% split(., cut(seq_along(.), 6, labels = FALSE))
    pca_sigs_fig_l = lapply(sig_l, function(sigs){
        sigs_pca_input = signatures_goodfeatures[, sigs, drop = F] %>% scale() %>% t()
        predicted = predict(prcom_96, newdata = sigs_pca_input)
        prcom_withsigs = prcom_96
        prcom_withsigs$x = rbind(prcom_withsigs$x, predicted)
        pca_part_sigs_fig = ggbiplot(prcom_withsigs, obs.scale = 1, ellipse = T, var.axes = F,
                                groups = c(as.vector(factors$experiment), rownames(sigs_pca_input)), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
            theme_classic() +
            annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
        return(pca_part_sigs_fig)
    })
    
    
    #For 6 way spectrum
    prcom_6 = type_occur %>% t() %>% scale() %>% t() %>% prcomp(scale = T, center = T)
    prcom_var = prcom_6$sdev * prcom_6$sdev
    nr_sig_pcs = bsDimension(prcom_var)
    nr_pcs = prcom_6$sdev %>% length()
    expected_var_pc1 = brokenStick(1, nr_pcs)
    expected_var_pc2 = brokenStick(2, nr_pcs)
    pca_sixway_fig = ggbiplot(prcom_6, obs.scale = 1, var.axes = F, ellipse = T, 
                              groups = factors$experiment, choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic() +
        annotate("text", y = 3, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant. (Based on 7 base substitutions.)"))
    
    pca_sixway_axis_fig = ggbiplot(prcom_6, obs.scale = 1, var.axes = T, ellipse = T, 
                                   groups = factors$experiment, choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic() +
        annotate("text", y = 3, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant. (Based on 7 base substitutions.)"))
    
    prcom_6 = type_occur %>% prcomp(scale = T, center = T)
    prcom_var = prcom_6$sdev * prcom_6$sdev
    nr_sig_pcs = bsDimension(prcom_var)
    nr_pcs = prcom_6$sdev %>% length()
    expected_var_pc1 = brokenStick(1, nr_pcs)
    expected_var_pc2 = brokenStick(2, nr_pcs)
    pca_sixway_lessscale_fig = ggbiplot(prcom_6, obs.scale = 1, var.axes = F, ellipse = T, 
                              groups = factors$experiment, choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic() +
        annotate("text", y = 3, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant. (Based on 7 base substitutions. Not scaled per sample)"))
    
    pca_sixway_axis_lessscale_fig = ggbiplot(prcom_6, obs.scale = 1, var.axes = T, ellipse = T, 
                                   groups = factors$experiment, choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
        theme_classic() +
        annotate("text", y = 3, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant. (Based on 7 base substitutions. Not scaled per sample)"))
    
    #Also make a heatmap of it
    ha = HeatmapAnnotation(as.data.frame(factors$experiment))
    heat_fig = mut_mat %>% scale() %>% Heatmap(top_annotation = ha)
    
    #t-SNE
    mut_mat_scaled = mut_mat %>% scale() %>% t()
    
    if (nrow(mut_mat_scaled) > 10){
        KL_div = Inf
        for (i in 1:10){ #Repeat t-SNE ten times. Pick the one, that has the lowest KL divergence.
            tsne_loop = Rtsne(mut_mat_scaled, perplexity = 5, max_iter = 5000, pca_center = T, pca_scale = F)
            KL_div_loop = tsne_loop$itercosts[length(tsne_loop$itercosts)]
            if (KL_div_loop < KL_div){
                tsne = tsne_loop
                KL_div = KL_div_loop
            }
        }
        
        d_tsne = as.data.frame(tsne$Y) %>% dplyr::mutate("Group" = factors$experiment)
        tsne_fig = ggplot(d_tsne, aes(x = V1, y = V2, colour = Group)) +
            geom_point(size = 1) +
            labs(x = "", y = "", title = "t-SNE") +
            theme_bw() +
            theme(text = element_text(size = 24),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()
            )
    } else{
        tsne_fig = NULL
    }
    ####_______________________Run mutsigpipe_____________________________####
    do_mutsigpipe(out_dir_mut, paste0(exp_name, "_combined"))
    
    
    ###_____________________Signatures manually. (Without rurikas extensive filtering)___________####
    signatures = read.table(cosmic_sig_fname, sep = "\t", header = T)
    signatures = as.matrix(signatures[,-c(1,2)])
    
    clust = dist(t(signatures)) %>% hclust()
    sig_order = colnames(signatures)[clust$order]
    dhc = as.dendrogram(clust)
    ddata = dendro_data(dhc, type = "rectangle")
    sig_dendrogram = ggplot(segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
        theme_dendro()
    
    #sig contribution
    fit_res = fit_to_signatures(mut_mat_group, signatures)
    select = which(rowSums(fit_res$contribution) > 10)
    sigcontri_fig = plot_contribution(fit_res$contribution[select,,drop = F], signatures[,select], coord_flip = F) +
        theme(text = element_text(size = 20))
    sigcontri_abs_fig = plot_contribution(fit_res$contribution[select,,drop = F], signatures[,select], coord_flip = F, mode = "absolute") +
        theme(text = element_text(size = 20))
    sigcontri_heat_fig = plot_contribution_heatmap(fit_res$contribution, cluster_samples = F, method = "complete")
    
    rel_contr = prop.table(fit_res$contribution, 2) #Divides by colsums
    diff_contr = rel_contr[,1, drop = F] - rel_contr[,2, drop = F]
    colnames(diff_contr) = "Difference"
    diff_contr %<>% as.data.frame() %>% rownames_to_column("Signature") %>% mutate(Signature = factor(Signature, levels = Signature))
    name_diff = paste0("Difference (", colnames(rel_contr)[1], " - ", colnames(rel_contr)[2], ")")
    diff_contr$sample = name_diff
    sigcontri_diff_heat_fig = ggplot(diff_contr, aes(x = Signature, fill = Difference, y = sample), order = NULL) +
        geom_tile(color = "white") +
        scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Difference") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = NULL, y = NULL)

    
    #Cos simularity
    cos_sim_samples_signatures = cos_sim_matrix(mut_mat_group, signatures)
    sig_cosim_heat_fig = plot_cosine_heatmap(cos_sim_samples_signatures, cluster_rows = F, col_order = sig_order)
    sig_cosim_heat_fig = plot_grid(sig_dendrogram, sig_cosim_heat_fig, align = "v", rel_heights = c(0.2, 1), axis = "lr", nrow = 2, scale = c(1.062,1))
    
    cos_sim_diff = cos_sim_samples_signatures[1,, drop = F] - cos_sim_samples_signatures[2,, drop = F]
    rownames(cos_sim_diff) = "Difference"
    cos_sim_diff %<>% t() %>% as.data.frame() %>% rownames_to_column("Signature") %>% mutate(Signature = factor(Signature, levels = sig_order))
    name_diff = paste0("Difference (", rownames(cos_sim_samples_signatures)[1], " - ", rownames(cos_sim_samples_signatures)[2], ")")
    cos_sim_diff$sample = name_diff
    sig_cossim_diff_heat_fig = ggplot(cos_sim_diff, aes(x = Signature, fill = Difference, y = sample), order = NULL) +
        geom_tile(color = "white") +
        scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Difference") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = NULL, y = NULL)
    sig_cossim_diff_heat_fig = plot_grid(sig_dendrogram, sig_cossim_diff_heat_fig, align = "v", rel_heights = c(0.2, 1), axis = "lr", nrow = 2, scale = c(1.062,1))
    
    
    
    #Look at similarity between the reconstructed spectra (based on the cosmic signatures) and the original spectra
    cos_sim_ori_rec = cos_sim_matrix(mut_mat_group, fit_res$reconstructed)
    cos_sim_ori_rec = as.data.frame(diag(cos_sim_ori_rec))
    colnames(cos_sim_ori_rec) = "cos_sim"
    cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
    orivsrec_fig = ggplot(cos_sim_ori_rec, aes(y = cos_sim, x = sample)) +
        geom_bar(stat = "identity", fill = "skyblue3") +
        coord_cartesian(ylim = c(0.6, 1), expand = F) +
        labs(y = "Cosine similarity\n original VS reconstructed", x = "") +
        theme_bw() +
        theme(text = element_text(size = 20)) +
        geom_hline(aes(yintercept = 0.95))
    
    
    #Perform cosine similarity for diff off the the mut mats.
    mut_mat_rel = prop.table(mut_mat_group, 2)
    mut_mat_diff = mut_mat_rel[,1, drop = F] - mut_mat_rel[,2, drop = F]
    cossime_sig_diff = cos_sim_matrix(mut_mat_diff, signatures)
    name_exp = paste0(colnames(mut_mat_rel)[1], "-", colnames(mut_mat_rel)[2], "\nall features")
    rownames(cossime_sig_diff) = name_exp
    cossime_sig_diff_long = cossime_sig_diff %>% as.data.frame() %>% rownames_to_column("sample") %>% gather(-sample, key = "Signature", value = "cossim")
    cossime_sig_diff_long$Signature = factor(cossime_sig_diff_long$Signature, levels = colnames(cossime_sig_diff))
    cossim_sig_diff_fig = ggplot(cossime_sig_diff_long, aes(x = Signature, fill = cossim, y = sample), order = NULL) +
        geom_tile(color = "white") +
        scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Difference") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(x = NULL, y = NULL)
    
    
    #Take either the negative or the positive features from the diff.
    pos_features = mut_mat_diff > 0
    cossime_sig_diff_pos = cos_sim_matrix(mut_mat_diff[pos_features,, drop = F], signatures[pos_features,, drop = F])
    name_exp_pos = paste0(colnames(mut_mat_rel)[1], "-", colnames(mut_mat_rel)[2])
    rownames(cossime_sig_diff_pos) = name_exp_pos
    cossim_sig_diff_pos_fig = plot_cosine_heatmap(cossime_sig_diff_pos, cluster_rows = F)
    cossime_sig_diff_neg = cos_sim_matrix(abs(mut_mat_diff[!pos_features,, drop = F]), signatures[!pos_features,, drop = F])
    name_exp_neg = paste0(colnames(mut_mat_rel)[2], "-", colnames(mut_mat_rel)[1])
    rownames(cossime_sig_diff_neg) = name_exp_neg
    cossim_sig_diff_neg_fig = plot_cosine_heatmap(cossime_sig_diff_neg, cluster_rows = F)
    
    
    
    print("Performed manual signature refitting")
    
   
    ####_________________________Rainfall plot_____________________________####
    #All samples combined
    rain_fig = plot_rainfall(gr_all, chromosomes = chromosomes, title = "All samples", ylim = 1e09) + coord_cartesian(ylim = c(1,1e09))
    
    #Per sample
    n = length(grl)
    rain_fig_l = vector("list", n)
    for (i in seq(1, n)){
        name = names(grl)[i]
        if (seqnames(grl[[name]]) %>% as.vector() %>% duplicated() %>% sum() == 0){
            print(paste0("No mutations in the same chromosome for sample: ", name))
            next
        }
        rain_sample_fig = plot_rainfall(grl[[name]], chromosomes = chromosomes, title = name, ylim = 1e09) + coord_cartesian(ylim = c(1,1e09))
        rain_fig_l[[i]] = rain_sample_fig
    }
    
    #Per group
    n = nlevels(factors$experiment)
    rain_fig_group_l = vector("list", n)
    for (i in seq(1, n)){
        exp = levels(factors$experiment)[i]
        gr_exp = gr_all[gr_all$experiment == exp]
        rain_group_fig = plot_rainfall(gr_exp, chromosomes = chromosomes, title = exp, ylim = 1e09) + coord_cartesian(ylim = c(1,1e09))
        rain_fig_group_l[[i]] = rain_group_fig
    }
    
    print("Created rainfall plots")
    
    ####_________________________________Create transcription strand bias figures_________####
    genes = genes(Homo.sapiens)
    
    #Per sample
    mut_mat_s = mut_matrix_stranded(grl, ref_genome, genes)
    strand_counts = strand_occurrences(mut_mat_s, by = names(grl))
    strand_bias = strand_bias_test(strand_counts)
    p6_trans_fig = plot_strand(strand_counts, mode = "relative")
    p6_trans2_fig = plot_strand_bias(strand_bias)
    p6_trans_abs_fig = plot_strand(strand_counts, mode = "absolute")
    
    p192_trans_fig = plot_192_profile(mut_mat_s, ymax = 0.25)
    
    #Per group
    strand_counts_group = strand_occurrences(mut_mat_s, by = factors$experiment)
    strand_bias_group = strand_bias_test(strand_counts_group)
    p6_trans_group_fig = plot_strand(strand_counts_group, mode = "relative")
    p6_trans2_group_fig = plot_strand_bias(strand_bias_group)
    p6_trans_abs_group_fig = plot_strand(strand_counts_group, mode = "absolute")
    
    mut_mat_s_group = mut_matrix_stranded(grl_group, ref_genome, genes)
    p192_trans_group_fig = plot_192_profile(mut_mat_s_group, ymax = 0.25)
    
    #Check if there is a bias in all genes
    seqlevels(genes, pruning.mode = "coarse") = chromosomes
    genes_count = GenomicRanges::reduce(genes)
    freq = getSeq(get(ref_genome), genes_count) %>% alphabetFrequency(., baseOnly = T)
    nuc_total = freq %>% colSums() %>% enframe() %>% mutate(name = factor(name, levels = c("A", "T", "C", "G", "other")))
    nuc_total_fig = ggplot(data = nuc_total, aes(x = name, y = value, fill = name)) +
        geom_bar(stat = "identity") +
        labs(title = "Nucleotide counts in genes", x = "Count", y = "Nucleotide") +
        theme_bw() +
        theme(text = element_text(size = 20), legend.position = "none")
    
    
    ####_____________Create replication bias figures____________####
    #Take only regions with replication direction
    slopes_bed = read_tsv(slopes_bed_fname)
    direction = slopes_bed %>% filter(strand_info != "No direction") %>% dplyr::select(-slope) %>% makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = T)
    direction$strand_info = as.factor(direction$strand_info)
    
    #Create figures per sample
    mut_mat_rep = mut_matrix_stranded(grl, ref_genome, direction, mode = "replication")
    strand_counts_rep = strand_occurrences(mut_mat_rep, by = names(grl))
    strand_bias_rep = strand_bias_test(strand_counts_rep)
    p6_rep_fig = plot_strand(strand_counts_rep, mode = "relative")
    p6_rep_abs_fig = plot_strand(strand_counts_rep, mode = "absolute")
    p6_rep2_fig = plot_strand_bias(strand_bias_rep)
    
    mut_mat_rep = mut_matrix_stranded(grl, ref_genome, direction, mode = "replication")
    p192_rep_fig = plot_192_profile(mut_mat_rep, ymax = 0.25)
    
    #Create figures per group
    strand_counts_rep_group = strand_occurrences(mut_mat_rep, by = factors$experiment)
    strand_bias_rep_group = strand_bias_test(strand_counts_rep_group)
    p6_rep_group_fig = plot_strand(strand_counts_rep_group, mode = "relative")
    p6_rep_abs_group_fig = plot_strand(strand_counts_rep_group, mode = "absolute")
    p6_rep2_group_fig = plot_strand_bias(strand_bias_rep_group)
    
    mut_mat_rep_group = mut_matrix_stranded(grl_group, ref_genome, direction, mode = "replication")
    p192_rep_group_fig = plot_192_profile(mut_mat_rep_group, ymax = 0.25)
    
    print("Created transcription and replication strand bias figures")
    
    ####__________________Write out all plots everything_________________####
    pdf(paste0(out_dir_mut, exp_name, "mutpatterns.pdf"), width = 10)
    #Group figs
    print(spectrum_group_fig)
    print(spectrum96_group_fig)
    print(spectrum96_group_fig2)
    print(rain_fig_group_l)
    grid.arrange(p6_trans_group_fig, p6_trans2_group_fig)
    grid.arrange(p6_trans_abs_group_fig, p6_trans2_group_fig)
    print(p192_trans_group_fig)
    grid.arrange(p6_rep_group_fig, p6_rep2_group_fig)
    grid.arrange(p6_rep_abs_group_fig, p6_rep2_group_fig)
    print(p192_rep_group_fig)
    print(pca_axis_fig)
    print(pca_fig)
    print(pca_fetus_fig)
    print(pca_23_fig)
    print(pca_34_fig)
    print(pca_sig1_fig)
    print(pca_sigs_fig_l)
    print(pca_sixway_axis_fig)
    print(pca_sixway_fig)
    print(pca_sixway_axis_lessscale_fig)
    print(pca_sixway_lessscale_fig)
    print(tsne_fig)
    print(heat_fig)
    print(pc_sig_fig)
    print(pc_sig_fig_pos)
    print(pc_sig_fig_neg)
    print(pc_sig2_fig_pos)
    print(pc_sig2_fig_neg)
    print(sigcontri_fig)
    print(sigcontri_abs_fig)
    print(sigcontri_heat_fig)
    print(sigcontri_diff_heat_fig)
    print(sig_cosim_heat_fig)
    print(sig_cossim_diff_heat_fig)
    print(cossim_sig_diff_fig)
    print(cossim_sig_diff_pos_fig)
    print(cossim_sig_diff_neg_fig)
    print(orivsrec_fig)
    
    #everything figs
    print(rain_fig)
    print(nuc_total_fig)
    
    #per sample figs
    print(spectrum_fig)
    print(spectrum96_figs)
    print(rain_fig_l)
    grid.arrange(p6_trans_fig, p6_trans2_fig)
    grid.arrange(p6_trans_abs_fig, p6_trans2_fig)
    print(p192_trans_fig)
    grid.arrange(p6_rep_fig, p6_rep2_fig)
    grid.arrange(p6_rep_abs_fig, p6_rep2_fig)
    print(p192_rep_fig)

    dev.off()
}

#This function is very slow, because it calculates the distribution for everything.
do_distribution = function(grl, exp_name, factors, out_dir_mut, chromosomes){
    repli_fname = "~/surfdrive/Shared/Boxtel_General/Data/Replication_timing/ENCODE/blood/B_lymphocytes.bed"
    gata1_chip_fname = "~/hpc/pmc_vanboxtel/data/Chip_seq/GATA1/ENCFF957CWW.bed"
    hypermeth_fname = "~/hpc/pmc_vanboxtel/data/Bisulphite/G199.hyper_meth.bs_call.GRCh37.20150707_ranges.bed"
    hypometh_fname = "~/hpc/pmc_vanboxtel/data/Bisulphite/G199.hypo_meth.bs_call.GRCh37.20150707_ranges.bed"
    ####____________________add exp info to the grl. Also create a grl for the groups and a gr with everything.___________####
    
    
    grl = keepSeqlevels(grl, chromosomes)
    for (i in seq(1, nrow(factors))){
        exp = factors$experiment[i]
        grl[[i]]$experiment = exp
    }
    gr_all = unlist(grl) %>% sort()
    grl_group = split(gr_all, gr_all$experiment)
    
    factors$experiment = as.factor(factors$experiment)
    
    ####_____________________Read in regulatory regions_____________________####
    regulatory = useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature", GRCh = 37)
    regsites = getBM(attributes = c("chromosome_name", "chromosome_start",  "chromosome_end", "feature_type_name"), mart = regulatory)
    features = unique(regsites$feature_type_name)
    regregions_list = lapply(features, function(feature){
        regsitesfeature = dplyr::filter(regsites, feature_type_name == feature)
        Grangeobj = GenomicRanges::reduce(GRanges(regsitesfeature$chromosome_name, IRanges(regsitesfeature$chromosome_start, regsitesfeature$chromosome_end)))
        return(Grangeobj)
    })
    regregions_list = GRangesList(regregions_list)
    names(regregions_list) = features
    seqlevelsStyle(regregions_list) = "UCSC"
    seqlevels(regregions_list, pruning.mode = "fine") = chromosomes
    
    #exons
    exons_gr = exons(Homo.sapiens)
    strand(exons_gr) = "*"
    seqlevels(exons_gr, pruning.mode = "coarse") = chromosomes
    exons_gr %<>% sort() %>% GenomicRanges::reduce()
    regregions_list$exons = exons_gr
    
    #gene bodies
    genes_gr = genes(Homo.sapiens)
    strand(genes_gr) = "*"
    seqlevels(genes_gr, pruning.mode = "coarse") = chromosomes
    genes_gr %<>% sort() %>% GenomicRanges::reduce()
    regregions_list$genes = genes_gr
    
    #Gata chip-seq on erythroblast male (encode)
    gata_gr = get_bed_as_gr(gata1_chip_fname, chromosomes)
    regregions_list$gata1_chip = gata_gr
    
    #hyper and hypo methylation for a hematopoietic multipotent progenitor cell (blueprint)
    hypermeth_gr = get_bed_as_gr(hypermeth_fname, chromosomes)
    regregions_list$hypermethylation = hypermeth_gr
    hypometh_gr = get_bed_as_gr(hypometh_fname, chromosomes)
    regregions_list$hypomethylation = hypometh_gr
    
    #replication timing data
    repli_bed = read_tsv(repli_fname, col_types = c("ciid"))
    colnames(repli_bed) = c("CHROM", "START", "END", "B_Lymphocytes")
    repli_bed %<>% dplyr::filter(B_Lymphocytes > 0.5) #Remove sites with very low values, because these are generally missing values.
    repli_gr = makeGRangesFromDataFrame(repli_bed, keep.extra.columns = T, starts.in.df.are.0based = T) %>% sort()
    seqlevels(repli_gr, pruning.mode = "coarse") = chromosomes
    
    early_gr = repli_gr[repli_gr$B_Lymphocytes >= 60] %>% sort() %>% GenomicRanges::reduce()
    regregions_list$early = early_gr
    intermediate_gr = repli_gr[repli_gr$B_Lymphocytes < 60 & repli_gr$B_Lymphocytes > 33] %>% sort() %>% GenomicRanges::reduce()
    regregions_list$intermediate = intermediate_gr
    late_gr = repli_gr[repli_gr$B_Lymphocytes <= 33] %>% sort() %>% GenomicRanges::reduce()
    regregions_list$late = late_gr
    
    print("Read in regulatory regions")
    
    ####_________________________________Look at repli-seq data_______________________####
    median_timing = median(repli_gr$B_Lymphocytes)
    
    
    mtch = findOverlaps(gr_all, repli_gr)
    muts_repli = gr_all[queryHits(mtch)] %>% mcols()
    repli_time = repli_gr[subjectHits(mtch)] %>% mcols()
    repli_df = cbind(muts_repli, repli_time) %>% as.data.frame() %>% mutate(experiment = as.vector(experiment))
    
    my_comparisons = combn(unique(repli_df$experiment), 2, simplify = F)
    replitime_fig = ggplot(data = repli_df, aes(x = experiment, y = B_Lymphocytes)) +
        geom_boxplot(outlier.shape = NA) +
        geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 1) +
        geom_hline(aes(yintercept = median_timing)) +
        labs(title = "Difference in replication time", x = "", y = "Replication time") +
        theme_bw() +
        theme(text = element_text(size = 20)) +
        stat_compare_means(comparisons = my_comparisons, label = "p.adj") +
        stat_compare_means(size = 6, label.y = 85)
    
    
    
    ####________________________Look at genomic distribution________________####
    #factors$callableloci_merged_surf = gsub("~/hpc/pmc_vanboxtel/projects/Freek_trees/.*/callableloci", "~/surfdrive/Shared/Projects/Freek/Freek_trees/CallableLoci", factors$callableloci_merged)
    #factors %<>% dplyr::mutate(callableloci_merged = ifelse(file.exists(callableloci_merged_surf), callableloci_merged_surf, callableloci_merged))
    factors_distr = factors %>% dplyr::filter(callableloci_merged_true)
    n = nrow(factors_distr)
    distr_l = vector("list", n)
    for (i in seq(1, n)){
        bed = get_bed_as_gr(factors_distr$callableloci_merged[i], chromosomes)
        distr = genomic_distribution(grl[i], list(bed), regregions_list)
        distr_l[[i]] = distr
        print(paste0("Calculated the distribution for sample: ", names(grl[i])))
    }
    distr = bind_rows(distr_l)
    
    #per group
    factors$experiment = factor(factors$experiment, levels = sort(levels(factors$experiment))) #Reorder factor levels, so they are alphabetically displayed.
    sample_2_exp = tibble("sample" = names(grl), "experiment" = factors$experiment)
    experiment_order = left_join(distr[,"sample", drop = F], sample_2_exp, by = "sample") %>% pull(experiment)
    distr_group_test = enrichment_depletion_test(distr, by = experiment_order)
    write_tsv(distr_group_test, paste0(out_dir_mut, exp_name, "_distr_vals.txt"))

    distr_group_betweenexps = enrichment_depletion_between_test(distr, by = experiment_order)
    write_tsv(distr_group_betweenexps, paste0(out_dir_mut, exp_name, "_distr_betweenexps_vals.txt"))

    distr_group_fig1 = plot_enrichment_depletion(distr_group_test[1:12,])
    distr_group_fig2 = plot_enrichment_depletion(distr_group_test[13:16,])
    distr_group_fig3 = plot_enrichment_depletion(distr_group_test[17:18,])
    distr_group_fig4 = plot_enrichment_depletion(distr_group_test[19:22,])
    distr_group_fig5 = plot_enrichment_depletion(distr_group_test[23:28,])
    
    #per sample
    distr_figs = lapply(distr_l, function(x){
        distr_test = enrichment_depletion_test(x)
        #distr_group_test$fdr = distr_group_test$pval %>% p.adjust(method = "fdr")
        #distr_group_test$significant = ifelse(distr_group_test$fdr < 0.1, "*", "")
        fig = plot_enrichment_depletion(distr_test)
        return(fig)
    })
    
    print("Looked at the genomic distribution")
    
    
    ####__________________Write out all plots everything_________________####
    
    pdf(paste0(out_dir_mut, exp_name, "mutpatterns_distribution.pdf"), width = 10)
    #Group figs
    print(distr_group_fig1)
    print(distr_group_fig2)
    print(distr_group_fig3)
    print(distr_group_fig4)
    print(distr_group_fig5)
    print(replitime_fig)

    #per sample figs
    print(distr_figs)
    
    dev.off()
    
}

#This function runs the better signature refitting.
do_better_refitting = function(grl, exp_name, factors, out_dir_mut, chromosomes){
    
    cosmic_sig_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_working_signatures_incl_hspc.txt"
    signatures = read.table(cosmic_sig_fname, sep = "\t", header = T)
    signatures = as.matrix(signatures[,-c(1,2)])
    
    signatures_incl_hspc_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_working_signatures_incl_hspc.txt"
    signatures_incl_hspc = read.table(signatures_incl_hspc_fname, sep = "\t", header = T)
    signatures_incl_hspc = as.matrix(signatures_incl_hspc[,-c(1,2)])
    interest_sigs_f = colnames(signatures_incl_hspc) %in% c("SBS1", "SBS5", "SBS18", "HSPC")
    signatures_incl_hspc = signatures_incl_hspc[,interest_sigs_f, drop = F]
    signatures_1_5 = signatures[,c("SBS1", "SBS5")]
    
    
    grl = keepSeqlevels(grl, chromosomes, pruning.mode = "fine")
    mut_mat = mut_matrix(grl, ref_genome)
    mut_mat_group = pool_mut_mat(mut_mat, factors$experiment)
    
    
    ####________________________Better_sig_refitting_____________________####
    fit_to_signatures(mut_mat_group, signatures)
    
    refit_out = fit_to_signatures_selection(mut_mat_group, signatures)
    write_rds(refit_out, paste0(out_dir_mut, exp_name, "refit_out.rds"))
    
    pc1 = plot_contribution_nonmf(refit_out$fit_res$contribution, mode = "relative") +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    pc2 = plot_contribution_nonmf(refit_out$fit_res$contribution, mode = "absolute")
    refit_contri_fig = plot_grid(pc1, pc2, align = "v", nrow = 2)
    
    sigs_fit = rownames(refit_out$fit_res$contribution) %>% unique()
    signatures_fit = signatures[,sigs_fit, drop = F]
    
    contri_tb = refit_out$fit_res$contribution %>% 
        as.data.frame() %>% 
        rownames_to_column("Signature") %>% 
        tidyr::gather(key = "Sample", value = "Contribution", -Signature)
    non_stacked_bar_fig = ggplot(contri_tb, aes(x = fct_rev(Signature), y = Contribution, fill = Signature)) +
        geom_bar(stat = "identity") +
        facet_grid(Sample ~ .) +
        theme_classic() +
        labs(x = "Signature", y = "Contribution") +
        theme(text = element_text(size = 20)) +
        coord_flip()
    
    contri = refit_out$fit_res$contribution
    contri[is.na(contri)] = 0
    contri = prop.table(contri, 2)
    contri_tb = contri %>% 
        as.data.frame() %>% 
        rownames_to_column("Signature") %>% 
        tidyr::gather(key = "Sample", value = "Rel_Contribution", -Signature)
    
    rel_non_stacked_bar_fig = ggplot(contri_tb, aes(x = fct_rev(Signature), y = Rel_Contribution, fill = Signature)) +
        geom_bar(stat = "identity") +
        facet_grid(Sample ~ .) +
        theme_classic() +
        labs(x = "Signature", y = "Contribution") +
        theme(text = element_text(size = 20)) +
        coord_flip(ylim = c(0,1)) +
        scale_y_continuous(labels = scales::percent)
    
    sim_decay_fig = refit_out$sim_decay_fig #Show how the cosine similarity is decreased as signatures are removed
    refit_heatmap_fig = plot_contribution_heatmap(refit_out$fit_res$contribution, cluster_samples = F) #Show a heatmap of the refitting.
    
    #Perform the signature refitting with both selection of signatures and bootstrapping.
    contri_boots = fit_to_signatures_bootstrapped(mut_mat_group, signatures, method = "selection", max_delta = 0.05)
    write_rds(contri_boots, paste0(out_dir_mut, exp_name, "contri_boots.rds"))
    
    nr_selected_fig = plot_nr_signatures_selected(contri_boots) #Plot how many signatures were selected during the different bootstrapping iterations.
    fraction_fig = plot_fraction_contri(contri_boots) #Plot how often each signature was selected.
    contri_boots_fig = plot_boots_contri(contri_boots) #Plot the contribution of bootstrapped selected signatures
    rel_contri_boots_fig = plot_boots_contri(contri_boots, mode = "relative") #Plot the contribution of bootstrapped selected signatures
    ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat_group, contri_boots, signatures) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
    sig_cor_figs = plot_sig_contri_cor(contri_boots) #Plot the signature contribution correlations between the signatures.
    
    #Perform the signature refitting with bootstrapping, on the signatures, that were most often selected during previous bootstrapping.
    common_sigs = get_common_signatures(contri_boots, signatures, frac_sel_lim = 0.5)
    contri_boots_common = fit_to_signatures_bootstrapped(mut_mat_group, common_sigs, method = "mutpatterns_default")
    common_contri_boots_fig = plot_boots_contri(contri_boots_common, mode = "relative") #Plot the contribution of bootstrapped common signatures
    common_ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat_group, contri_boots_common, common_sigs) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
    
    
    #If you have an issue with your signatures being too similar, you can merge them, before the refitting.
    merged_sigs = merge_signatures(signatures)
    contri_boots_merge = fit_to_signatures_bootstrapped(mut_mat_group, merged_sigs, method = "selection")
    merge_nr_selected_fig = plot_nr_signatures_selected(contri_boots_merge) #Plot how many signatures were selected during the different bootstrapping iterations.
    merge_fraction_fig = plot_fraction_contri(contri_boots_merge) #Plot how often each signature was selected.
    merge_contri_boots_fig = plot_boots_contri(contri_boots_merge) #Plot the contribution of bootstrapped selected signatures
    merge_ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat_group, contri_boots_merge, merged_sigs) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
    merge_sig_cor_figs = plot_sig_contri_cor(contri_boots_merge) #Plot the signature contribution correlations between the signatures.
    
    
    
    #Difference between signatures. This can again be done on common or merged signatures if desired.
    contri_boots_perm = fit_to_signatures_permuted(mut_mat_group, merged_sigs, method = "selection", n_boots = 2000, max_delta = 0.05) #For pvalue calculation a higher number of boots is recommended
    diff_sigs = calc_diff_sigs_perm(contri_boots_perm, mut_mat_group, merged_sigs, max_delta = 0.05)
    perm_contri_boots_figs = plot_perm_contri(contri_boots_perm, diff_sigs, mode = "absolute")
    perm_contri_boots_rel_figs = plot_perm_contri(contri_boots_perm, diff_sigs, mode = "relative")
    
    contri_boots_perm2 = fit_to_signatures_permuted(mut_mat_group, common_sigs, method = "mutpatterns_default", n_boots = 2000) #For pvalue calculation a higher number of boots is recommended
    diff_sigs2 = calc_diff_sigs_perm(contri_boots_perm2, mut_mat_group, common_sigs, method = "mutpatterns_default")
    perm_contri_boots_figs2 = plot_perm_contri(contri_boots_perm2, diff_sigs2, mode = "absolute")
    perm_contri_boots_rel_figs2 = plot_perm_contri(contri_boots_perm2, diff_sigs2, mode = "relative")
    
    #Signatures of special interest
    contri_boots_interest = fit_to_signatures_bootstrapped(mut_mat_group, signatures_incl_hspc, method = "mutpatterns_default")
    interest_contri_boots_fig = plot_boots_contri(contri_boots_interest, mode = "relative") #Plot the contribution of bootstrapped common signatures
    interest_ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat_group, contri_boots_interest, signatures_incl_hspc) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions

    contri_boots_perm3 = fit_to_signatures_permuted(mut_mat_group, signatures_incl_hspc, method = "mutpatterns_default", n_boots = 2000) #For pvalue calculation a higher number of boots is recommended
    diff_sigs3 = calc_diff_sigs_perm(contri_boots_perm3, mut_mat_group, signatures_incl_hspc, method = "mutpatterns_default")
    perm_contri_boots_figs3 = plot_perm_contri(contri_boots_perm3, diff_sigs3, mode = "absolute")
    perm_contri_boots_rel_figs3 = plot_perm_contri(contri_boots_perm3, diff_sigs3, mode = "relative")
 
    contri_boots_perm4 = fit_to_signatures_permuted(mut_mat_group, signatures_1_5, method = "mutpatterns_default", n_boots = 2000) #For pvalue calculation a higher number of boots is recommended
    diff_sigs4 = calc_diff_sigs_perm(contri_boots_perm4, mut_mat_group, signatures_1_5, method = "mutpatterns_default")
    perm_contri_boots_figs4 = plot_perm_contri(contri_boots_perm4, diff_sigs4, mode = "absolute")
    perm_contri_boots_rel_figs4 = plot_perm_contri(contri_boots_perm4, diff_sigs4, mode = "relative")
    
    contri_boots_perm5 = fit_to_signatures_permuted(mut_mat_group, signatures_fit, method = "mutpatterns_default", n_boots = 2000) #For pvalue calculation a higher number of boots is recommended
    diff_sigs5 = calc_diff_sigs_perm(contri_boots_perm5, mut_mat_group, signatures_fit, method = "mutpatterns_default")
    perm_contri_boots_figs5 = plot_perm_contri(contri_boots_perm5, diff_sigs5, mode = "absolute")
    perm_contri_boots_rel_figs5 = plot_perm_contri(contri_boots_perm5, diff_sigs5, mode = "relative")
    
    
    
    pdf(paste0(out_dir_mut, exp_name, "better_refitting.pdf"))
    print(sim_decay_fig)
    print(refit_contri_fig)
    print(non_stacked_bar_fig)
    print(rel_non_stacked_bar_fig)
    print(refit_heatmap_fig)
    print(nr_selected_fig)
    print(fraction_fig)
    print(contri_boots_fig)
    print(rel_contri_boots_fig)
    print(ori_vs_rec_fig)
    print(sig_cor_figs)
    print(common_contri_boots_fig)
    print(common_ori_vs_rec_fig)
    print(merge_nr_selected_fig)
    print(merge_fraction_fig)
    print(merge_contri_boots_fig)
    print(merge_ori_vs_rec_fig)
    print(merge_sig_cor_figs)
    print(perm_contri_boots_figs)
    print(perm_contri_boots_rel_figs)
    print(perm_contri_boots_figs2)
    print(perm_contri_boots_rel_figs2)
    print(interest_contri_boots_fig)
    print(interest_ori_vs_rec_fig)
    print(perm_contri_boots_figs3)
    print(perm_contri_boots_rel_figs3)
    print(perm_contri_boots_figs4)
    print(perm_contri_boots_rel_figs4)
    print(perm_contri_boots_figs5)
    print(perm_contri_boots_rel_figs5)
    
    dev.off()
    
    
    ####SPARSE
    sparse_sig_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sparse_signatures.txt"
    sparse_signatures = read.table(sparse_sig_fname, sep = "\t", header = T)
    sparse_signatures = as.matrix(sparse_signatures[,-c(1,2)])
    
    #Perform signature refitting with selection of signatures. Here this is done without refitting, to show how the signature selection process is functioning.
    refit_out = fit_to_signatures_selection(mut_mat_group, sparse_signatures)
    sim_decay_fig = refit_out$sim_decay_fig #Show how the cosine similarity is decreased as signatures are removed
    refit_heatmap_fig = plot_contribution_heatmap(refit_out$fit_res$contribution, cluster_samples = F) #Show a heatmap of the refitting.
    
    #Perform the signature refitting with both selection of signatures and bootstrapping.
    contri_boots = fit_to_signatures_bootstrapped(mut_mat_group, sparse_signatures, method = "selection")
    nr_selected_fig = plot_nr_signatures_selected(contri_boots) #Plot how many signatures were selected during the different bootstrapping iterations.
    fraction_fig = plot_fraction_contri(contri_boots) #Plot how often each signature was selected.
    contri_boots_fig = plot_boots_contri(contri_boots, mode = "relative") #Plot the contribution of bootstrapped selected signatures
    ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat_group, contri_boots, sparse_signatures) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
    sig_cor_figs = plot_sig_contri_cor(contri_boots) #Plot the signature contribution correlations between the signatures.
    
    pdf(paste0(out_dir_mut, exp_name, "sparse_sigs_better_refitting.pdf"))
    print(sim_decay_fig)
    print(refit_heatmap_fig)
    print(nr_selected_fig)
    print(fraction_fig)
    print(contri_boots_fig)
    print(ori_vs_rec_fig)
    print(sig_cor_figs)
    dev.off()
    
    print("Performed better_refitting for the signatures")
}



do_indel_refitting = function(indel_counts, out_dir_mut, exp_name){
    id_sigs = read_csv("~/surfdrive/Shared/Boxtel_General/Data/Sigs/ID/sigProfiler_ID_signatures.csv")
    id_sigs = as.matrix(id_sigs[,-1])
    indel_counts_m = indel_counts %>% dplyr::select(-muttype, -muttype_sub) %>% as.matrix()
    
    ####________________________Better_sig_refitting_____________________####
    refit_out = fit_to_signatures_selection(indel_counts_m, id_sigs)
    sim_decay_fig = refit_out$sim_decay_fig #Show how the cosine similarity is decreased as signatures are removed
    refit_heatmap_fig = plot_contribution_heatmap(refit_out$fit_res$contribution, cluster_samples = F) #Show a heatmap of the refitting.
    
    
    #Perform the signature refitting with both selection of signatures and bootstrapping.
    contri_boots = fit_to_signatures_bootstrapped(indel_counts_m, id_sigs, method = "selection")
    nr_selected_fig = plot_nr_signatures_selected(contri_boots) #Plot how many signatures were selected during the different bootstrapping iterations.
    fraction_fig = plot_fraction_contri(contri_boots) #Plot how often each signature was selected.
    contri_boots_fig = plot_boots_contri(contri_boots) #Plot the contribution of bootstrapped selected signatures
    ori_vs_rec_fig = plot_cosine_bootstrapped(indel_counts_m, contri_boots, id_sigs) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
    sig_cor_figs = plot_sig_contri_cor(contri_boots) #Plot the signature contribution correlations between the signatures.
    
    #Perform the signature refitting with bootstrapping, on the signatures, that were most often selected during previous bootstrapping.
    common_sigs = get_common_signatures(contri_boots, id_sigs, frac_sel_lim = 0.5)
    contri_boots_common = fit_to_signatures_bootstrapped(indel_counts_m, common_sigs, method = "mutpatterns_default")
    common_contri_boots_fig = plot_boots_contri(contri_boots_common, mode = "relative") #Plot the contribution of bootstrapped common signatures
    common_ori_vs_rec_fig = plot_cosine_bootstrapped(indel_counts_m, contri_boots_common, common_sigs) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
    
    
    #Difference between signatures. This can again be done on common or merged signatures if desired.
    contri_boots_perm = fit_to_signatures_permuted(indel_counts_m, id_sigs, method = "selection", n_boots = 2000, max_delta = 0.05) #For pvalue calculation a higher number of boots is recommended
    diff_sigs = calc_diff_sigs_perm(contri_boots_perm, indel_counts_m, id_sigs, max_delta = 0.05)
    perm_contri_boots_figs = plot_perm_contri(contri_boots_perm, diff_sigs, mode = "absolute")
    
    contri_boots_perm2 = fit_to_signatures_permuted(indel_counts_m, common_sigs, method = "mutpatterns_default", n_boots = 2000) #For pvalue calculation a higher number of boots is recommended
    diff_sigs2 = calc_diff_sigs_perm(contri_boots_perm2, indel_counts_m, common_sigs, method = "mutpatterns_default")
    perm_contri_boots_figs2 = plot_perm_contri(contri_boots_perm2, diff_sigs2, mode = "absolute")
    
    pdf(paste0(out_dir_mut, exp_name, "_indel_refitting.pdf"))
    print(sim_decay_fig)
    print(refit_heatmap_fig)
    print(nr_selected_fig)
    print(fraction_fig)
    print(contri_boots_fig)
    print(ori_vs_rec_fig)
    print(sig_cor_figs)
    print(common_contri_boots_fig)
    print(common_ori_vs_rec_fig)
    print(perm_contri_boots_figs)
    print(perm_contri_boots_figs2)
    dev.off()
}
    
