#Get the fixed effects from a lmm fitted with lme.
get_fixed = function(mixed_model){
    lme_sum = summary(mixed_model)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    ci = intervals(mixed_model, which = "fixed")$fixed %>% as_tibble()
    vals = cbind(vals, ci[,c("lower", "upper")])
    call = lme_sum$call %>% as.character() %>% paste0(collapse = " ")
    vals %<>% dplyr::mutate(logLik = mixed_model$logLik, BIC = lme_sum$BIC, AIC = lme_sum$AIC, call = call)
    return(vals)
}



get_vaf = function(vcf, sample){
    ad = geno(vcf)$AD[,sample]
    vaf = sapply(ad, function(x) x[[2]] / sum(x))
    vaf[is.na(vaf)] = 0 #For sites with 0 reads
    return(vaf)
}


sort_vcf = function(vcf){

    if (class(vcf)[1] != "CollapsedVCF"){
            warning(paste0("This function was tested on the Collapsed VCF class. It may not work on the supplied input vcf, which is a ", class(vcf[1]), " object."))
    }
    
    if (length(vcf) == 0){
        warning("vcf is empty. Returning vcf as is.")
        return(vcf)
    }
    gr = granges(vcf)
    gr$id = seq(1, length(gr))
    gr_sorted = sort(gr)
    vcf = vcf[gr_sorted$id,]
    return(vcf)
}

contri_pie = function(vaf_tb, name){
    blank_theme <- theme_minimal()+
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold")
        )
    
    contri_fig = ggplot(vaf_tb, aes(y = contribution, x = "", fill = group)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = c("grey", "red"), guide = F) +
        #scale_fill_grey() +  
        blank_theme +
        labs(title = name) +
        geom_text(aes(y = contribution[1]/3, x = 1.1, 
                      label = percent(contribution[1])), size=18)
    return(contri_fig)
}

get_allbulks = function(other_bulks, bulk){
    if (!is.na(other_bulks)){
        other_bulks_v = strsplit(other_bulks, ";")[[1]]
        all_bulks = c(other_bulks_v, bulk) %>% paste0(collapse = "|")
    } else{
        all_bulks = bulk
    }
    return(all_bulks)
}

get_surveyed_region = function(surveyed_fname){
    if (file.exists(surveyed_fname)){
        surveyed_region = scan(surveyed_fname, quiet = T, n = 1)
        return(surveyed_region)
    } else{
        return(NA)
    }
}

#Get all snvs for a sample
get_all_snvs = function(overview_sample, chroms){
    sample = overview_sample[["sample"]]
    vcf_fnames = c(overview_sample[["shared_vcf"]], overview_sample[["unique_vcf"]])
    vcf_fnames = vcf_fnames[file.exists(vcf_fnames)]
    gr_list = lapply(vcf_fnames, function(x) {
        vcf = readVcf(x, genome = "hg19")
        gt = geno(vcf)$GT
        sample_col = grep(sample, colnames(gt))
        presence_f = gt[,sample_col] == "0/1" | gt[,sample_col] == "1/1"
        vcf = vcf[presence_f,]
        gr = granges(vcf)
        print(paste0("Read in vcf: ", x))
        return(gr)
    })
    gr = do.call(c, gr_list)
    gr %<>% sort()
    seqlevelsStyle(gr) = "UCSC"
    seqlevels(gr, pruning.mode = "coarse") = chroms
    
    return(gr)
}

#Create a mutation matrix for the high vaf and the low vaf muts from samples seperately. (Run on the all unique snv vcf files.)
get_vafsplit = function(vcfs, samples, chromosomes){
    mut_mat_l = vector("list", length(vcfs))
    for (i in 1:length(vcfs)){
        vcf = readVcf(vcfs[i], genome = "hg19")
        sample = samples[i]
        vaf = get_vaf(vcf, sample)
        low_vaf = vaf <= median(vaf, na.rm = T)
        gr = granges(vcf)
        seqlevelsStyle(gr) = "UCSC"
        seqlevels(gr, pruning.mode = "coarse") = chromosomes
        gr_lowvaf = gr[low_vaf]
        gr_highvaf = gr[!low_vaf]
        grl = GRangesList(gr_lowvaf, gr_highvaf)
        names(grl) = paste0(sample, c("_lowvaf", "_highvaf"))
        mut_mat_sample = mut_matrix(grl, ref_genome)
        mut_mat_l[[i]] = mut_mat_sample
    }
    mut_mat = do.call(cbind, mut_mat_l)
    return(mut_mat)
}

#Get the subclonal mutations.
write_subclonal_vcfs = function(overview_samples, out_dir_base){
    apply(overview_samples, 1, function(row){
        fetus = row[["fetus"]]
        sample = row[["sample"]]
        bulk = row[["bulk"]]
        
        snvfi_out = paste0(out_dir_base, fetus, "/", sample, "/", sample, "_", bulk, "_Q50_CGQ10_SGQ99_PASS_20X_VAF0.1_NoY_nonBlacklist_final.vcf")
        vcf = readVcf(snvfi_out, genome = "hg19")
        vaf = get_vaf(vcf, sample)
        vaf_f = vaf < 0.3
        
        sample_col = grep(sample, samples(header(vcf)))
        bulk_col = grep(bulk, samples(header(vcf)))
        dp_f = geno(vcf)$DP[,sample_col] >= 30 & geno(vcf)$DP[,bulk_col] >= 30
        
        ad = geno(vcf)$AD
        alt = apply(ad, 1:2, function(x) x[[1]][2])
        ad_f = alt[,sample_col] >= 5 & rowSums(alt[,-sample_col, drop = F]) == 0 #At least 5 reads in the clone. 0 reads in the rest.
        
        gq_f = geno(vcf)$GQ[,bulk_col] >= 99
        
        subclonal_vcf = vcf[vaf_f & dp_f & ad_f & gq_f]
        if (nrow(subclonal_vcf) >= 1){
            vcf_out_fname = row[["subclonal_vcf"]]
            bed_validated_fname = row[["subclonal_validated"]]
            bed_tovalidate_fname = row[["subclonal_tovalidate"]]
            vcf_out(subclonal_vcf, vcf_out_fname, bed_validated_fname, bed_tovalidate_fname)
        }
    })
}


#Get all indels for a sample
get_all_indels = function(overview_sample, chroms){
    sample = overview_sample[["sample"]]
    vcf_fnames = c(overview_sample[["indel_shared"]], overview_sample[["indel_uniq"]])
    vcf_fnames = vcf_fnames[file.exists(vcf_fnames)]
    
    if (isEmpty(vcf_fnames)){
        gr = GRanges()
        genome(gr) = "hg19"
        return(gr)
    }
    gr_list = lapply(vcf_fnames, function(x) {
        vcf = readVcf(x, genome = "hg19")
        gt = geno(vcf)$GT
        sample_col = grep(sample, colnames(gt))
        presence_f = gt[,sample_col] == "0/1" | gt[,sample_col] == "1/1"
        vcf = vcf[presence_f,]
        gr = granges(vcf)
        print(paste0("Read in vcf: ", x))
        return(gr)
    })
    gr = do.call(c, gr_list)
    gr %<>% sort()
    seqlevelsStyle(gr) = "UCSC"
    seqlevels(gr, pruning.mode = "coarse") = chroms
    
    return(gr)
}

#Get all unique indels for a sample
get_uniq_indels = function(overview_sample, chroms){
    sample = overview_sample[["sample"]]
    vcf_fname = overview_sample[["indel_uniq"]]
    if (!file.exists(vcf_fname)){
        gr = GRanges()
        genome(gr) = "hg19"
        return(gr)
    }
    
    if (isEmpty(vcf_fname)){
        gr = GRanges()
        genome(gr) = "hg19"
        return(gr)
    }
    vcf = readVcf(vcf_fname, genome = "hg19")
    gr = granges(vcf)
    print(paste0("Read in vcf: ", vcf_fname))
    seqlevelsStyle(gr) = "UCSC"
    seqlevels(gr, pruning.mode = "coarse") = chroms
    return(gr)
}

#Get all indels per fetus. Each indel is only counted once.
get_all_indels_fetus = function(fetus, overview_samples, chromosomes){
    print(fetus)
    overview_fetus = overview_samples %>% dplyr::filter(fetus == !!fetus)
    shared_vcf = overview_fetus$indel_shared[overview_fetus$indel_shared_true] %>% unique()
    unique_vcfs = overview_fetus$indel_uniq[overview_fetus$indel_uniq_true]
    vcfs = c(shared_vcf, unique_vcfs)
    if (length(vcfs) == 0){
        gr = GRanges()
        seqlevels(gr, pruning.mode = "coarse") = chromosomes
        return(gr)
    }
    grl = lapply(vcfs, function(x){
        vcf = readVcf(x, genome = "hg19")
        gr = granges(vcf)
        return(gr)
    })
    gr = do.call(c, grl)
    gr %<>% sort()
    seqlevelsStyle(gr) = "UCSC"
    seqlevels(gr, pruning.mode = "coarse") = chromosomes
    return(gr)
}

get_snvs_indels_vcf = function(overview_sample, chroms){
    sample = overview_sample[["sample"]]
    vcf_fnames = c(overview_sample[["indel_shared"]], overview_sample[["indel_uniq"]], overview_sample[["shared_vcf"]], overview_sample[["unique_vcf"]])
    vcf_fnames = vcf_fnames[file.exists(vcf_fnames)]
    
    if (isEmpty(vcf_fnames)){
        gr = GRanges()
        genome(gr) = "hg19"
        return(gr)
    }
    vcf_list = lapply(vcf_fnames, function(x) {
        vcf = readVcf(x, genome = "hg19")
        gt = geno(vcf)$GT
        sample_col = grep(sample, colnames(gt))
        presence_f = gt[,sample_col] == "0/1" | gt[,sample_col] == "1/1"
        vcf = vcf[presence_f,]
        return(vcf)
    })
    vcf = do.call(rbind, vcf_list) %>% sort_vcf()
    return(vcf)
}

#Function to identify the vcf annotated with a gene in a specified cancer gene symbol vector and predicted to possibly have an effect.
get_possible_drivers = function(vcf, cancer_gene_symbols){
    cosm_f = grepl("COSM", names(vcf))
    
    ann = info(vcf)$ANN
    nmuts = nrow(vcf)
    possible_driver = vector("list", nmuts)
    for (i in 1:nmuts){
        x = ann[i]
        x = unlist(x)
        list_of_lists = strsplit(x, "\\|")
        gene_list = lapply(list_of_lists, function(x) x[[4]]) %>% unique()
        genes = unlist(gene_list)
        nr_cancer_genes =  genes %in% cancer_gene_symbols %>% sum()
        
        eff = lapply(list_of_lists, function(x) x[[3]]) %>% unique() %>% unlist()
        effect_f = eff %in% c("MODERATE", "HIGH") %>% sum()
        
        possible_driver_mut = (nr_cancer_genes > 0 | cosm_f[i]) & effect_f > 0
        possible_driver[[i]] = possible_driver_mut
    }
    possible_driver = unlist(possible_driver)
    vcf_posdriver = vcf[possible_driver]
    return(vcf_posdriver)
}

#Identify and write out potential driver mutations.
create_possible_driver_vcfs = function(overview_sample, chromosomes, out_dir_base, cancer_gene_symbols){
    fetus = overview_sample[["fetus"]]
    sample = overview_sample[["sample"]]
    vcf = get_snvs_indels_vcf(overview_sample, chromosomes)
    
    vcf_pos_driver = get_possible_drivers(vcf, cancer_gene_symbols)
    
    #Check if any pos driver muts were found
    if (length(vcf_pos_driver) == 0){
        print(paste0("No drivers in sample: ", sample, " in fetus: ", fetus))
        return(0)
    }
    nr_muts = length(vcf_pos_driver)
    print(paste0(nr_muts, " Possible drivers were found in sample: ", sample, " in fetus: ", fetus))
    
    #Write a vcf with the identified possible driver mutations.
    somatic_dir = paste0(out_dir_base, "/", fetus, "/somatic_drivers/")
    if (!dir.exists(somatic_dir)){
        dir.create(somatic_dir)
    }
    writeVcf(vcf_pos_driver, paste0(somatic_dir, sample, "_possible_drivers.vcf"))
    return(0)
}


#Function to pool the granges in a grl that belong to the same exp
pool_grl_exp = function(grl, factor){
    
    #Split muts based on the experiment
    factor = as.factor(factor)
    n_exp = nlevels(factor)
    gr_list_exps = vector("list", n_exp)
    for (i in 1:n_exp){
        exp = levels(factor)[i]
        index = grep(exp, factor)
        grl_exp = grl[index]
        gr_exp = unlist(grl_exp)
        gr_list_exps[[i]] = gr_exp
    }
    grl_exp = GRangesList(gr_list_exps)
    names(grl_exp) = levels(factor)
    return(grl_exp)
}

#Create overview of all samples, with the location of their vcfs. Samples are mined from the _diffclones.vcf, so this needs to exist.
get_overview_fetus = function(fetus, out_dir_base){
    fetus_name = fetus[[1]]
    bulk = fetus[[2]]
    gender = fetus[[3]]
    trisomy = fetus[[4]]
    samples = strsplit(fetus[[5]], ";")[[1]]
    other_bulks = fetus[[6]]
    age_weeks = fetus[[7]]
    #vcf_shared = readVcf(paste0(out_dir_base, fetus_name, "/shared/", fetus_name, "_somatic_filtered_noblacklist_MQ_diffclones.vcf"))
    #samples_bulk = colnames(geno(vcf_shared)$GT)
    #samples = samples_bulk[-grep(bulk, samples_bulk)]
    overview_samples = tibble("fetus" = fetus_name, "sample" = samples, "bulk" = bulk, "gender" = gender, "trisomy" = trisomy, "other_bulks" = other_bulks, "age_weeks" = age_weeks, "uniq_vcf_inbulk_unfiltered" = paste0(out_dir_base, fetus_name, "/dirichlet/", samples, "_inbulk.vcf"), "uniq_vcf_inbulk" = paste0(out_dir_base, fetus_name, "/dirichlet/", samples, "_inbulk_true.vcf"), "uniq_vcf_notbulk_unfiltered" = paste0(out_dir_base, fetus_name, "/dirichlet/", samples, "_notbulk.vcf"), "uniq_vcf_notbulk" = paste0(out_dir_base, fetus_name, "/dirichlet/", samples, "_notbulk_true.vcf"), "shared_vcf" = paste0(out_dir_base, fetus_name, "/shared/", fetus_name, "_shared_true.vcf"), "unique_vcf" = paste0(out_dir_base, fetus_name, "/dirichlet/", samples, "_true.vcf"), "shared_inbulk_vcf" = paste0(out_dir_base, fetus_name, "/shared/", fetus_name, "_shared_inbulk.vcf"), "shared_notbulk_vcf" = paste0(out_dir_base, fetus_name, "/shared/", fetus_name, "_shared_notbulk.vcf"))
    return(overview_samples)
}

create_overview_samples = function(fetuses, out_dir_base){
    overview_samples = lapply(fetuses, function(x) get_overview_fetus(x, out_dir_base)) %>% bind_rows() %>% mutate(origin = ifelse(grepl("BM", sample), "BM", ifelse(grepl("SI", sample, ignore.case = T), "SI", ifelse(grepl("-I-", sample), "SI", ifelse(grepl("INTESTINE", sample), "SI", "Liver")))))
    overview_samples %<>% mutate(callableloci_merged = paste0(out_dir_base, fetus, "/callableloci/", sample, "_", bulk, "_CallableLoci_merged_autosomal_X.bed"), callableloci_merged_auto = paste0(out_dir_base, fetus, "/callableloci/", sample, "_", bulk, "_CallableLoci_merged_autosomal.bed"))
    overview_samples$celltype = ifelse(grepl("SI|Si|-I-", overview_samples$sample), "SI", ifelse(grepl("-L-", overview_samples$sample), "Liver", "HSC"))
    
    overview_samples$shared_vcf_true = file.exists(overview_samples$shared_vcf)
    overview_samples$shared_vcf_inbulk_true = file.exists(overview_samples$shared_inbulk_vcf)
    overview_samples$shared_vcf_notbulk_true = file.exists(overview_samples$shared_notbulk_vcf)
    overview_samples$callableloci_merged_true = file.exists(overview_samples$callableloci_merged)
    overview_samples$callableloci_merged_auto_true = file.exists(overview_samples$callableloci_merged_auto)
    
    overview_samples$uniq_vcf_inbulk_unfiltered_true = file.exists(overview_samples$uniq_vcf_inbulk_unfiltered)
    overview_samples$uniq_vcf_inbulk_true = file.exists(overview_samples$uniq_vcf_inbulk)
    overview_samples$uniq_vcf_notbulk_unfiltered_true = file.exists(overview_samples$uniq_vcf_notbulk_unfiltered)
    overview_samples$uniq_vcf_notbulk_true = file.exists(overview_samples$uniq_vcf_notbulk)
    overview_samples$unique_vcf_true = file.exists(overview_samples$unique_vcf)
    
    #Present in bulk, not in clones
    overview_samples %<>% mutate(vcf_inbulk_not_clones_unfiltered = paste0(out_dir_base, fetus, "/shared/", fetus, "_subclonalbulks_notclones.vcf"))
    overview_samples %<>% mutate(vcf_inbulk_not_clones_unfiltered_true = file.exists(overview_samples$vcf_inbulk_not_clones_unfiltered))
    overview_samples %<>% mutate(vcf_inbulk_not_clones = paste0(out_dir_base, fetus, "/shared/", fetus, "_subclonalbulks_notclones_true.vcf"))
    
    #Indels
    overview_samples %<>% mutate(indel_uniq_notbulk_indelfi = paste0(out_dir_base, fetus, "/", sample, "_indel/", sample, "_", bulk, "_QUAL250_RD10_VAF0.1_GQ10_flank100_MQ60_INDELs_autosomal_X_noEvidenceControl_noBlacklist.vcf"), indel_uniq_notbulk_indelfi_true = file.exists(indel_uniq_notbulk_indelfi))
    overview_samples %<>% mutate(indel_uniq_dir = paste0(out_dir_base, fetus, "/indels_unique/"))
    overview_samples %<>% mutate(indel_uniq_notbulk_filtered = paste0(indel_uniq_dir, sample, "_notbulk.vcf"), indel_uniq_notbulk_filtered_true = file.exists(indel_uniq_notbulk_filtered))
    overview_samples %<>% mutate(indel_uniq_notbulk = paste0(indel_uniq_dir, sample, "_notbulk_true.vcf"), indel_uniq_notbulk_true = file.exists(indel_uniq_notbulk))
    overview_samples %<>% mutate(indel_uniq_notbulk_bed_tovalidate = paste0(indel_uniq_dir, sample, "_notbulk_tovalidate.bed"), indel_uniq_notbulk_bed_validated = paste0(indel_uniq_dir, sample, "_notbulk_validated.bed"))
    overview_samples %<>% mutate(indel_uniq_inbulk_filtered = paste0(indel_uniq_dir, sample, "_inbulk.vcf"), indel_uniq_inbulk_filtered_true = file.exists(indel_uniq_inbulk_filtered))
    overview_samples %<>% mutate(indel_uniq_inbulk = paste0(indel_uniq_dir, sample, "_inbulk_true.vcf"), indel_uniq_inbulk_true = file.exists(indel_uniq_inbulk))
    overview_samples %<>% mutate(indel_uniq_inbulk_manualadded = paste0(indel_uniq_dir, sample, "_inbulk_manualadded.vcf"), indel_uniq_inbulk_manualadded_true = file.exists(indel_uniq_inbulk_manualadded))
    overview_samples %<>% mutate(indel_uniq_inbulk_bed_tovalidate = paste0(indel_uniq_dir, sample, "_inbulk_tovalidate.bed"), indel_uniq_inbulk_bed_validated = paste0(indel_uniq_dir, sample, "_inbulk_validated.bed"))
    overview_samples %<>% mutate(indel_uniq = paste0(indel_uniq_dir, sample, "_true.vcf"), indel_uniq_true = file.exists(indel_uniq))
    
    overview_samples %<>% mutate(indel_somatic = paste0(out_dir_base, fetus, "/indels/", fetus, "_somatic_filtered_noblacklist_MQ.vcf"), indel_shared_dir = paste0(out_dir_base, fetus, "/indels_shared/"), indel_shared_filtered = paste0(indel_shared_dir, fetus, "_shared.vcf"), indel_shared_filtered_true = file.exists(indel_shared_filtered), indel_shared = paste0(indel_shared_dir, fetus, "_shared_true.vcf"), indel_shared_true = file.exists(indel_shared), indel_shared_tovalidate = paste0(indel_shared_dir, fetus, "_tovalidate.bed"), indel_shared_validated = paste0(indel_shared_dir, fetus, "_validated.bed"))
    overview_samples %<>% mutate(indel_shared_inbulk = paste0(indel_shared_dir, fetus, "_shared_inbulk.vcf"), indel_shared_inbulk_true = file.exists(indel_shared_inbulk), indel_shared_notbulk = paste0(indel_shared_dir, fetus, "_shared_notbulk.vcf"), indel_shared_notbulk_true = file.exists(indel_shared_notbulk))
    overview_samples %<>% mutate(indel_shared_manualadded = paste0(indel_shared_dir, fetus, "_shared_true_manualadded.vcf"), indel_shared_manualadded_true = file.exists(indel_shared_manualadded))
    
    #Subclonal snvs
    overview_samples %<>% mutate(subclonal_vcf = paste0(out_dir_base, fetus, "/subclonal/", fetus, "_", sample, "_subclonal.vcf"), subclonal_tovalidate = paste0(out_dir_base, fetus, "/subclonal/", fetus, "_", sample, "_subclonal_tovalidate.bed"), subclonal_validated = paste0(out_dir_base, fetus, "/subclonal/", fetus, "_", sample, "_subclonal_validated.bed"))
    
    #other
    overview_samples$surveyed_fname = paste0(out_dir_base, overview_samples$fetus, "/callableloci/", overview_samples$sample, "_", overview_samples$bulk, "_surveyed_autosomal_X.txt")
    overview_samples$surveyed_region = sapply(overview_samples$surveyed_fname, get_surveyed_region)
    chr = paste0("chr", c(1:22,"X"))
    length_chr = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chr] %>% sum()
    overview_samples %<>% mutate(surveyed_region_pct = surveyed_region / length_chr)
    overview_samples$surveyed_fname_auto = paste0(out_dir_base, overview_samples$fetus, "/callableloci/", overview_samples$sample, "_", overview_samples$bulk, "_surveyed.txt")
    overview_samples$surveyed_region_auto = sapply(overview_samples$surveyed_fname_auto, get_surveyed_region)
    chr = paste0("chr", c(1:22))
    length_chr_auto = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chr] %>% sum()
    overview_samples %<>% mutate(surveyed_region_auto_pct = surveyed_region_auto / length_chr_auto)
    
    weeksbyyear = 52.177457
    overview_samples %<>% mutate(age_weeks = as.numeric(age_weeks), age_year = age_weeks / weeksbyyear)
    
    
    write_tsv(overview_samples, paste0(out_dir_base, "overview_samples.txt"))
    return(overview_samples)
}

#function to make a co-occurence heatmap, to see how often mutations co-occur between samples
co_occurence_map = function(gt_clones){
    
    if (nrow(gt_clones) < 1){
        return("There were no mutations in this category")
    }
    
    gt_clones[gt_clones == "0/0"] = 0
    gt_clones[gt_clones == "0/1" | gt_clones == "1/1"] = 1
    gt_clones %<>% mutate_all(list(~as.integer)) %>% as.matrix()
    gt_clones_t = t(gt_clones)
    co_occurence = gt_clones_t %*% gt_clones
    
    heat_fig = Heatmap(co_occurence, col = colorRamp2(seq(0,70,7.777778), brewer.pal(9, "OrRd")), column_dend_reorder = F, row_dend_reorder = F, cluster_rows = T, cluster_columns = T,
                       row_names_side = "left", column_names_side = "top", show_row_names = T, 
                       column_title = "Shared mutations", row_title = "",
                       heatmap_legend_param = list(title = "Number of shared mutations"),
                       row_names_gp = gpar(fontsize = 18), column_names_gp = gpar(fontsize = 18), column_title_gp = gpar(fontsize = 20), row_title_gp = gpar(fontsize = 20))
    return(heat_fig)
}

#Function to count how often each combination of mutations occurs, both for all informative mutations and for mutations present in bulk
count_combis = function(vcf, out_name, bulk, all_bulks){
    
    gt = geno(vcf)$GT
    bulk_columns = grep(all_bulks, colnames(gt))
    
    gt_combis_all = gt[, -bulk_columns] %>% as_tibble()
    colnames(gt_combis_all) = gsub("-", "_", colnames(gt_combis_all))
    combi_count = plyr::count(gt_combis_all, vars = colnames(gt_combis_all))
    write_tsv(combi_count, paste0(out_name, "_freq.txt"))
    
    vaf = get_vaf(vcf, bulk)
    combi_vaf = as_tibble(gt)
    combi_vaf$vaf = vaf
    
    #gt_combis_inbulk = gt %>% as.tibble() %>% dplyr::filter(UQ(as.name(bulk)) == "0/1") %>% dplyr::select(-bulk)
    if (nrow(combi_vaf) >= 1){
        write_tsv(combi_vaf, paste0(out_name, "_inbulk_vaf.txt"))
    }
}

#Function to get only sites from a vcf that are different between the clones (based on gt calls)
vcf_diffclones = function(sample_name, out_dir, vcf_diffclones_fname, bulk, all_bulks){
    vcf_fname = paste0(out_dir, sample_name, "_somatic_filtered_noblacklist_MQ.vcf")
    vcf = readVcf(vcf_fname)
    
    gt = geno(vcf)$GT
    bulk_columns = grep(all_bulks, colnames(gt))
    gt_clones = gt[, -bulk_columns]
    
    diff_clones = apply(gt_clones, 1, function(x) "0/0" %in% x & ("0/1" %in% x | "1/1" %in% x))
    
    vcf = vcf[diff_clones]
    writeVcf(vcf, vcf_diffclones_fname)
    return(vcf)
}

filter_row = function(row, nsamples, min_vaf, min_absent, min_present, min_dp = 20, min_pres_gq, mut_type){
    
    if (mut_type == "snv"){
        min_abs_gq = 10
    } else if (mut_type == "indel"){
        min_abs_gq = 99
    }
    
    absent = grep("0/0", row)
    present = grep("0/1|1/1", row)
    absent_goodgq = row[absent + nsamples] >= min_abs_gq
    present_goodgq = row[present + nsamples] >= min_pres_gq
    absent_gooddp = row[absent + 2*nsamples] >= min_dp
    absent_goodad = row[absent + 3*nsamples] <= 0
    present_gooddp = row[present + 2*nsamples] >= min_dp
    present_goodvaf = row[present + 4*nsamples] >= min_vaf
    
    good_absents = absent_goodgq & absent_gooddp & absent_goodad
    good_presents = present_goodgq & present_gooddp & present_goodvaf
    
    keep_row = sum(good_absents) >= min_absent & sum(good_presents) >= min_present
    return(keep_row)
}

vcf_gt_inf = function(vcf, gt_clones_cols){
    gt = geno(vcf)$GT %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_GT")) %>% dplyr::select(gt_clones_cols)
    gq = geno(vcf)$GQ %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_GQ")) %>% replace(., is.na(.), 0) %>% dplyr::select(gt_clones_cols)
    dp = geno(vcf)$DP %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_DP")) %>% replace(., is.na(.), 0) %>% dplyr::select(gt_clones_cols)
    ad = geno(vcf)$AD %>% as_tibble() %>% dplyr::select(gt_clones_cols)
    
    ad_alt = apply(ad, 1:2, function(x) x[[1]][2]) %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_ADALT"))
    vaf = apply(ad, 1:2, function(x) x[[1]][2] / sum(x[[1]])) %>% as_tibble() %>% rename_all(., function(x) paste0(names(.), "_vaf")) %>% replace(., is.na(.), 0)
    geno_tb = bind_cols(gt, gq, dp, ad_alt, vaf)
    return(geno_tb)
}

find_tree_sites = function(vcf, fetus, bulk, min_absent, min_present, trisomy, gender, mut_type, all_bulks){
    
    #Filter on gt
    gt = geno(vcf)$GT
    bulk_columns = grep(all_bulks, colnames(gt))
    gt_clones = gt[, -bulk_columns]
    
    shared_sites_f = apply(gt_clones, 1, function(x) table(x)["0/0"] >= min_absent & (sum("0/1" == x) + sum("1/1" == x)) >= min_present) %>% replace(., is.na(.), FALSE)
    vcf_shared = vcf[shared_sites_f]
    
    if (nrow(vcf_shared) < 1){
        print(paste0("No shared mutations for sample: ", fetus, ". The function is returning empty"))
        return(0)
    }
    
    samples = samples(header(vcf_shared))
    nsamples_bulk = length(samples)
    gt_cols = seq(1, nsamples_bulk)
    bulk_gt_cols = grep(all_bulks, samples)
    gt_clones_cols = gt_cols[!gt_cols %in% bulk_gt_cols]
    nsamples = nsamples_bulk - length(bulk_gt_cols)
    
    chroms_vcf = rowRanges(vcf_shared) %>% seqnames() %>% as.vector()
    other_chroms = unique(chroms_vcf)
    
    #Trisomy filter
    if (trisomy != "FALSE"){
        chr_tri = chroms_vcf == trisomy
        vcf_tri = vcf_shared[chr_tri,]
        geno_tb_tri = vcf_gt_inf(vcf_tri, gt_clones_cols)
        nsites_tri = nrow(geno_tb_tri)
        keep_rows_tri = sapply(seq(1, nsites_tri), function(i) filter_row(geno_tb_tri[i,], nsamples, min_vaf = 0.2, min_absent, min_present, min_dp = 20, min_pres_gq = 99, mut_type))
        vcf_tri = vcf_tri[keep_rows_tri]
        other_chroms = other_chroms[!other_chroms == trisomy]
    } else{
        vcf_tri = vcf_shared[0]
    }
    
    #X chromosome filter for men
    if (gender == "M"){
        chr_x = chroms_vcf == "X"
        vcf_x = vcf_shared[chr_x,]
        geno_tb_x = vcf_gt_inf(vcf_x, gt_clones_cols)
        nsites_x = nrow(geno_tb_x)
        keep_rows_x = sapply(seq(nsites_x), function(i) filter_row(geno_tb_x[i,], nsamples, min_vaf = 0.99, min_absent, min_present, min_dp = 10, min_pres_gq = 10, mut_type))
        vcf_x = vcf_x[keep_rows_x]
        other_chroms = other_chroms[!other_chroms == "X"]
    } else{
        vcf_x = vcf_shared[0]
    }
    
    #Filter for the other chromosomes
    chr_rest = chroms_vcf %in% other_chroms
    vcf_rest = vcf_shared[chr_rest,]
    geno_tb_rest = vcf_gt_inf(vcf_rest, gt_clones_cols)
    nsites_rest = nrow(geno_tb_rest)
    keep_rows_and_good_samples_rest = lapply(seq(nsites_rest), function(i) filter_row(geno_tb_rest[i,], nsamples, min_vaf = 0.3, min_absent, min_present, min_dp = 20, min_pres_gq = 99, mut_type))
    keep_rows_rest = lapply(keep_rows_and_good_samples_rest, function(x) x[[1]]) %>% unlist()
    vcf_rest = vcf_rest[keep_rows_rest]
    
    #Combine the output from the different filters
    vcf_shared_filtered = rbind(vcf_rest, vcf_tri, vcf_x)
    vcf_shared_filtered = sort_vcf(vcf_shared_filtered)
    
    return(vcf_shared_filtered)
}

vcf_out = function(vcf, vcf_out_fname, bed_validated_fname, bed_tovalidate_fname){
    
    writeVcf(vcf, vcf_out_fname)
    print(paste0("Wrote out the vcf: ", vcf_out_fname))
    
    #Create a bed to use in validating muts in igv
    ranges = rowRanges(vcf)
    bed = tibble("CHROM" = as.vector(seqnames(ranges)), "START" = start(ranges)-1, "END" = end(ranges), "Real_mut" = TRUE)
    
    #Filter tovalidate bed file based on sites that were already validated or falsified.
    if (file.exists(bed_validated_fname)){
        bed_validated = read_tsv(bed_validated_fname, col_types = cols(CHROM = "c"))
        bed = anti_join(bed, bed_validated, by = c("CHROM", "START", "END"))
    }
    if (nrow(bed) >= 1){
        write_tsv(bed, bed_tovalidate_fname)
    }
    return(0)
}

#Filter on igv validation
igv_filter = function(vcf_fname, bed_fname, vcf_out_fname){
    if (file.exists(vcf_fname) & file.exists(bed_fname)){
        vcf = readVcf(vcf_fname, genome = "hg19")
        
        #Check if any muts exist
        bed = read_tsv(bed_fname, col_types = cols(CHROM = "c"))
        if (nrow(bed) == 0){
            print(paste0("bed file  ", bed_fname, " was empty"))
            return(0)
        }
        
        #Check if any TRUE muts exist
        bed %<>% dplyr::filter(Real_mut)
        if (nrow(bed) == 0){
            print(paste0("bed file ", bed_fname,  " did not contain any true mutations"))
            return(0)
        }
        bed %<>% makeGRangesFromDataFrame(starts.in.df.are.0based = T) %>% sort()
        vcf_ranges = granges(vcf)
        true_index = findOverlaps(vcf_ranges, bed, type = "equal") %>% queryHits()
        vcf = vcf[true_index]
        if (nrow(vcf) > 0){
            writeVcf(vcf, vcf_out_fname)
        }
    }
    return(0)
}

