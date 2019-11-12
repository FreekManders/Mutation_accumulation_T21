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
library(scales)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/vcfFilterTree_functions.R")

out_dir_base = "~/hpc/pmc_vanboxtel/projects/Freek_trees/"
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

#Get info on samples
fetuses_tb = read_tsv(paste0(out_dir_base, "fetuses_inclewart.txt"), col_types = cols(.default = "c"))
fetuses = split(fetuses_tb, seq(1,nrow(fetuses_tb)))
fetuses = lapply(fetuses, unlist)
overview_samples = create_overview_samples(fetuses, out_dir_base)

####____________Identify shared mutations. these need to be validated in igv.______________####
for (fetus in fetuses){
    sample_name = fetus[[1]]
    bulk = fetus[[2]]
    gender = fetus[[3]]
    trisomy = fetus[[4]]
    other_bulks = fetus[[6]]
    out_dir = paste0(out_dir_base, sample_name, "/shared/")

    all_bulks = get_allbulks(other_bulks, bulk)
    
    #Filter for muts different between the clones
    vcf_diffclones_fname = paste0(out_dir, sample_name, "_somatic_filtered_noblacklist_MQ_diffclones.vcf")
    if (!file.exists(vcf_diffclones_fname)){
        vcf = vcf_diffclones(sample_name, out_dir, vcf_diffclones_fname, bulk, all_bulks)
        print(paste0("Made diffclones vcf for sample: ", sample_name))
    } else{
        vcf = readVcf(vcf_diffclones_fname)
        print(paste0("diffclones vcf for sample: ", sample_name, " already exists."))
    }

    #Perform filtering for mutations present in at least two, but absent in at least one clones.
    vcf_shared_filtered = find_tree_sites(vcf, sample_name, bulk, min_absent = 1, min_present = 2, trisomy, gender, mut_type = "snv", all_bulks)
    print(paste0("Performed filtering for sample: ", sample_name))
    
    #Write the output
    if (class(vcf_shared_filtered) != "numeric"){
        nr_shared = nrow(vcf_shared_filtered)
        print(paste0(nr_shared, " mutations were left after filtering"))
        vcf_shared_filtered_fname = paste0(out_dir, sample_name, "_shared.vcf")
        bed_validated_fname = paste0(out_dir, sample_name, "_shared_validated.bed")
        bed_tovalidate_fname = paste0(out_dir, sample_name, "_shared_tovalidate.bed")
        trash = vcf_out(vcf_shared_filtered, vcf_shared_filtered_fname, bed_validated_fname, bed_tovalidate_fname)
    }
    
    
    #Filter for unique mutations present in the bulk
    vcf_one_filtered = find_tree_sites(vcf, sample_name, bulk, min_absent = 1, min_present = 1, trisomy, gender, "snv", all_bulks)
    samples = samples(header(vcf_one_filtered))
    nsamples_bulk = length(samples)
    gt_cols = seq(1, nsamples_bulk)
    bulk_cols = grep(all_bulks, samples)
    ref_bulk_col = grep(bulk, samples)
    gt_clones_cols = gt_cols[!gt_cols %in% bulk_cols]
    nsamples = nsamples_bulk - length(bulk_cols)

    #Filter on whether mutations are unique and present in the bulk
    gt = geno(vcf_one_filtered)$GT
    #inbulk_f = gt[,ref_bulk_col] %in% c("0/1", "1/1")
    clones_called = gt[ ,gt_clones_cols] == "0/1" | gt[ ,gt_clones_cols] == "1/1"
    uniq_clone_f = rowSums(clones_called) == 1
    vaf_bulk = get_vaf(vcf_one_filtered, bulk)
    #ad_bulk = geno(vcf_one_filtered)$AD[,ref_bulk_col]
    #vaf_bulk = sapply(ad_bulk, function(x) x[[2]]/sum(x))
    vaf_f = vaf_bulk > 0
    vcf_uniq_inbulk = vcf_one_filtered[vaf_f & uniq_clone_f, ]
    gt = geno(vcf_uniq_inbulk)$GT
    
    
    #Write out unique mutations per clone
    samples_notbulk = colnames(gt)[gt_clones_cols]
    
    dirichlet_dir = paste0(out_dir_base, sample_name, "/dirichlet/")
    if (!dir.exists(dirichlet_dir)){
        dir.create(dirichlet_dir)
    }
    
    for (clone in samples_notbulk){
        called_clone_f = gt[,clone] %in% c("0/1", "1/1")
        vcf_uniq_clone = vcf_uniq_inbulk[called_clone_f,]
        if (nrow(vcf_uniq_clone) == 0){
            next
        }
        vcf_uniq_clone_fname = paste0(out_dir_base, sample_name, "/dirichlet/", clone, "_inbulk.vcf")
        bed_validated_fname = paste0(out_dir_base, sample_name, "/dirichlet/", clone, "_inbulk_validated.bed")
        bed_tovalidate_fname =  paste0(out_dir_base, sample_name, "/dirichlet/", clone, "_inbulk_tovalidate.bed")
        vcf_out(vcf_uniq_clone, vcf_uniq_clone_fname, bed_validated_fname, bed_tovalidate_fname)
    }
    
    
    
}
#update overview samples
overview_samples = create_overview_samples(fetuses, out_dir_base)

####______Check if there are any mutations (sub)clonally present in (both) bulks, but not in the clones. Also check for muts that are present in all clones, but subclonally present in the bulk_____________####
for (fetus in fetuses){
    sample_name = fetus[[1]]
    bulk = fetus[[2]]
    gender = fetus[[3]]
    trisomy = fetus[[4]]
    other_bulks = fetus[[6]]
    out_dir = paste0(out_dir_base, sample_name, "/shared/")
    
    all_bulks = get_allbulks(other_bulks, bulk)
    
    nr_bulks = length(strsplit(all_bulks, "\\|")[[1]])
    vcf = readVcf(paste0(out_dir, sample_name, "_somatic_filtered_noblacklist_MQ.vcf"))
    gt = geno(vcf)$GT
    bulk_cols = grep(all_bulks, colnames(gt))
    
    gt_bulks =  gt[,bulk_cols, drop = F]
    presence_f_tb = apply(gt_bulks, 2, function(x) x == "0/1" | x == "1/1")
    presence_f = rowSums(presence_f_tb) == ncol(presence_f_tb)
    gt_clones = gt[,-bulk_cols]
    absence_f_tb = apply(gt_clones, 2, function(x) x == "0/0")
    absence_f = rowSums(absence_f_tb) == ncol(absence_f_tb)
    combined_f = presence_f & absence_f
    
    vcf_filtered_gt = vcf[combined_f]
    nsites = nrow(vcf_filtered_gt)
    print(paste0(nsites, " SNVs were present in the bulks, but not in the clones, based on the GTs."))
    
    nsamples = ncol(gt)
    #min_absent = round((nsamples - nr_bulks + 0.01) / 2) #Ensures .5 is consistently rounded up
    min_absent = nsamples - nr_bulks
    geno_tb = vcf_gt_inf(vcf_filtered_gt, seq(1, nsamples))
    keep_rows = sapply(seq(1, nsites), function(i) filter_row(geno_tb[i,], nsamples, min_vaf = 0.01, min_absent = min_absent, min_present = nr_bulks, min_dp = 20, min_pres_gq = 99, mut_type = "snv"))
    vcf_filtered = vcf_filtered_gt[keep_rows,]
    nmuts = nrow(vcf_filtered)
    print(paste0(nmuts, " SNVs were present in the bulks, but not in the clones after filtering."))
    if (nmuts > 0){
        vcf_filtered_fname = paste0(out_dir, sample_name, "_subclonalbulks_notclones.vcf")
        bed_validated_fname = paste0(out_dir, sample_name, "_subclonalbulks_notclones_validated.bed")
        bed_tovalidate_fname = paste0(out_dir, sample_name, "_subclonalbulks_notclones_tovalidate.bed")
        trash = vcf_out(vcf_filtered, vcf_filtered_fname, bed_validated_fname, bed_tovalidate_fname)
    }
    
    
    
    #check for muts that are present in all clones, but subclonally present in the bulk
    presence_f = rowSums(absence_f_tb) == 0
    bulk_col = grep(bulk, colnames(gt))
    gt_bulk = gt[,bulk_col]
    presence_bulk_f = gt_bulk == "0/1" | gt_bulk == "1/1"
    vaf = get_vaf(vcf, bulk)
    vaf_f = vaf < 0.25
    total_f = presence_f & presence_bulk_f & vaf_f
    vcf_pres = vcf[total_f,]
    nsites = nrow(vcf_pres)
    print(paste0(nsites, " SNVs were subclonally present in the bulks, and in all clones, based on the GTs."))
    
    
    #filter_bulk
    geno_tb = vcf_gt_inf(vcf_pres, bulk_col)
    keep_rows_bulk = sapply(seq(1, nsites), function(i) filter_row(geno_tb[i,], nsamples = 1, min_vaf = 0.01, min_absent = 0, min_present = 1, min_dp = 20, min_pres_gq = 99, mut_type = "snv"))
    
    #filter_clones
    cols = seq(1, nsamples)
    sample_cols = cols[!cols %in% bulk_cols]
    geno_tb = vcf_gt_inf(vcf_pres, sample_cols)
    nclones = length(sample_cols)
    keep_rows_clones = sapply(seq(1, nsites), function(i) filter_row(geno_tb[i,], nsamples = nclones, min_vaf = 0.3, min_absent = 0, min_present = nclones, min_dp = 20, min_pres_gq = 99, mut_type = "snv"))
    
    keep_rows = keep_rows_bulk & keep_rows_clones
    
    vcf_filt = vcf_pres[keep_rows]
    nmuts = nrow(vcf_filt)
    print(paste0(nmuts, " SNVs were clonally present in all clones and subclonally in the bulk."))
    if (nmuts > 0){
        vcf_filtered_fname = paste0(out_dir, sample_name, "_subclonalbulks_allclones.vcf")
        bed_validated_fname = paste0(out_dir, sample_name, "_subclonalbulks_allclones_validated.bed")
        bed_tovalidate_fname = paste0(out_dir, sample_name, "_subclonalbulks_allclones_tovalidate.bed")
        trash = vcf_out(vcf_filt, vcf_filtered_fname, bed_validated_fname, bed_tovalidate_fname)
        
    }
}

####________Perform dirichlet modeling on all the mutations found by snvfi_________####
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/Perform_dirichlet.R")
overview_samples = read_tsv(paste0(out_dir_base, "overview_samples.txt"))
for (fetus in fetuses){
    fetus_name = fetus[[1]]
    bulk = fetus[[2]]
    gender = fetus[[3]]
    trisomy = fetus[[4]]
    fetus_dir = paste0(out_dir_base, fetus_name, "/")

    dirichlet_dir = paste0(out_dir_base, fetus_name, "/dirichlet/")
    if (!dir.exists(dirichlet_dir)){
        dir.create(dirichlet_dir)
    }
    setwd(dirichlet_dir)
    
    #Determine names of vcf files
    samples = overview_samples %>% dplyr::filter(fetus == fetus_name) %>% pull(sample)
    vcf_files = paste0(out_dir_base, fetus_name, "/", samples, "/", samples, "_", bulk, "_Q50_CGQ10_SGQ99_PASS_20X_VAF0.1_NoY_nonBlacklist_final.vcf")
    normal_size_f = file.info(vcf_files)$size <= 10000000
    too_large = vcf_files[!normal_size_f]
    if (length(too_large) > 0){
        print(paste0("File ", too_large, " was not used for dirichlet modeling because of its large size"))
    }
    vcf_files = vcf_files[normal_size_f]
    #samples = samples[normal_size_f]
    print(paste0("Determined the location of the vcf files for fetus: ", fetus_name))
    
    perform_dirichlet(vcf_files, fetus_dir = fetus_dir, fetus_name, trisomy, gender, bulk)
}
#update overview samples
overview_samples = create_overview_samples(fetuses, out_dir_base)


for (indel_dir in unique(c(overview_samples$indel_uniq_dir, overview_samples$indel_shared_dir))){
    if (!dir.exists(indel_dir)){
        dir.create(indel_dir)
    }
}
###______________look for shared indels and unique indels present in the bulk_________####
overview_shared = overview_samples %>% dplyr::filter(!duplicated(indel_somatic))
trash = apply(overview_shared, 1, function(row){
    
    vcf_fname = row[["indel_somatic"]]
    fetus = row[["fetus"]]
    bulk = row[["bulk"]]
    gender = row[["gender"]]
    trisomy = row[["trisomy"]]
    indel_shared_dir = row[["indel_shared_dir"]]
    vcf_filtered_fname = row[["indel_shared_filtered"]]
    bed_tovalidate_fname = row[["indel_shared_tovalidate"]]
    bed_validated_fname = row[["indel_shared_validated"]]
    other_bulks = row[["other_bulks"]]
    
    all_bulks = get_allbulks(other_bulks, bulk)
    
    vcf = readVcf(vcf_fname, genome = "hg19")
    
    #Identify shared mutations
    vcf_shared_filtered = find_tree_sites(vcf, fetus, bulk, min_absent = 1, min_present = 2, trisomy, gender, "indel", all_bulks)
    
    #Write the output
    if (class(vcf_shared_filtered) != "numeric"){
        nr_shared = nrow(vcf_shared_filtered)
        print(paste0(nr_shared, " mutations were left after filtering"))
        if (nr_shared > 0){
            trash = vcf_out(vcf_shared_filtered, vcf_filtered_fname, bed_validated_fname, bed_tovalidate_fname)
        }
    }
    
    #Identify uniq mutations present in the bulk
    vcf_one_filtered = find_tree_sites(vcf, fetus, bulk, min_absent = 1, min_present = 1, trisomy, gender, "indel", all_bulks)
    samples = samples(header(vcf_one_filtered))
    nsamples_bulk = length(samples)
    gt_cols = seq(1, nsamples_bulk)
    bulk_cols = grep(all_bulks, samples)
    gt_clones_cols = gt_cols[!gt_cols %in% bulk_cols]
    nsamples = nsamples_bulk - length(bulk_cols)
    
    #Filter on whether mutations are unique and present in the bulk
    gt = geno(vcf_one_filtered)$GT
    #inbulk_f = gt[,ref_bulk_col] %in% c("0/1", "1/1")
    clones_called = gt[ ,gt_clones_cols, drop = F] == "0/1" | gt[ ,gt_clones_cols, drop = F] == "1/1"
    uniq_clone_f = rowSums(clones_called) == 1
    vaf_bulk = get_vaf(vcf_one_filtered, bulk)
    #ad_bulk = geno(vcf_one_filtered)$AD[,ref_bulk_col]
    #vaf_bulk = sapply(ad_bulk, function(x) x[[2]]/sum(x))
    vaf_f = vaf_bulk > 0
    vcf_uniq_inbulk = vcf_one_filtered[vaf_f & uniq_clone_f, ]
    gt = geno(vcf_uniq_inbulk)$GT
    
    
    #write the output per clone
    samples_notbulk = colnames(gt)[gt_clones_cols]
    for (clone in samples_notbulk){
        overview_sample = overview_samples %>% dplyr::filter(sample == UQ(clone))
        if (nrow(overview_sample) == 0){#Check if clone should be filtered, or whether it was removed from analysis.
            next
        }
        called_clone_f = gt[,clone] %in% c("0/1", "1/1")
        vcf_uniq_clone = vcf_uniq_inbulk[called_clone_f,]
        if (nrow(vcf_uniq_clone) == 0){
            next
        }
        vcf_filtered_fname = overview_sample$indel_uniq_inbulk_filtered
        bed_tovalidate_fname = overview_sample$indel_uniq_inbulk_bed_tovalidate
        bed_validated_fname = overview_sample$indel_uniq_inbulk_bed_validated
        
        trash = vcf_out(vcf_uniq_clone, vcf_filtered_fname, bed_validated_fname, bed_tovalidate_fname)
    }
    return(0)
})



####_____________Look for subclonal mutations______________________________####

#Create subclonal dirs.
for (fetus in unique(overview_samples$fetus)){
    dir = paste0(out_dir_base, fetus, "/subclonal/")
    if (!dir.exists(dir)){
        dir.create(dir)
    }
}
trash = write_subclonal_vcfs(overview_samples, out_dir_base)

####_____________Filter indelfi output based on vaf and depth and gq________####

overview_samples_indelfi = overview_samples %>% dplyr::filter(indel_uniq_notbulk_indelfi_true)
trash = apply(overview_samples_indelfi, 1, function(row){
    vcf_fname = row[["indel_uniq_notbulk_indelfi"]]
    vcf_out_fname = row[["indel_uniq_notbulk_filtered"]]
    bed_tovalidate_fname = row[["indel_uniq_notbulk_bed_tovalidate"]]
    bed_validated_fname = row[["indel_uniq_notbulk_bed_validated"]]
    vcf_shared_filtered_fname = row[["indel_shared_filtered"]]
    gender = row[["gender"]]
    trisomy = row[["trisomy"]]
    bulk = row[["bulk"]]
    sample = row[["sample"]]
    
    vcf = readVcf(vcf_fname, genome = "hg19")
    vcf = vcf[!isSNV(vcf)]
    
    chroms_vcf = rowRanges(vcf) %>% seqnames() %>% as.vector()
    other_chroms = unique(chroms_vcf)
    #Filter on vaf, dp and gq
    #Trisomy filter
    if (trisomy != "FALSE"){
        chr_tri = chroms_vcf == trisomy
        vcf_tri = vcf[chr_tri,]
        ad = geno(vcf_tri)$AD[,sample]
        vaf = sapply(ad, function(x) x[[2]] / sum(x))
        dp_sample = geno(vcf_tri)$DP[,sample]
        dp_bulk = geno(vcf_tri)$DP[,bulk]
        gq_sample = geno(vcf_tri)$GQ[,sample]
        gq_bulk = geno(vcf_tri)$GQ[,bulk]
        keep_f = vaf >= 0.2 & dp_sample >= 20 & dp_bulk >= 20 & gq_sample >= 99 & gq_bulk >= 99
        vcf_tri = vcf_tri[keep_f]
        other_chroms = other_chroms[!other_chroms == trisomy]
    } else{
        vcf_tri = vcf[0]
    }
    #seperate filter for x
    if (gender == "M"){
        chr_x = chroms_vcf == "X"
        vcf_x = vcf[chr_x,]
        ad = geno(vcf_x)$AD[,sample]
        vaf = sapply(ad, function(x) x[[2]] / sum(x))
        dp_sample = geno(vcf_x)$DP[,sample]
        dp_bulk = geno(vcf_x)$DP[,bulk]
        gq_sample = geno(vcf_x)$GQ[,sample]
        gq_bulk = geno(vcf_x)$GQ[,bulk]
        keep_f = vaf >= 0.99 & dp_sample >= 10 & dp_bulk >= 10 & gq_sample >= 10 & gq_bulk >= 99
        vcf_x = vcf_x[keep_f]
        other_chroms = other_chroms[!other_chroms == "X"]
    } else{
        vcf_x = vcf[0]
    }
    #filter the other chroms
    chr_rest = chroms_vcf %in% other_chroms
    vcf_rest = vcf[chr_rest,]
    ad = geno(vcf_rest)$AD[,sample]
    vaf = sapply(ad, function(x) x[[2]] / sum(x))
    dp_sample = geno(vcf_rest)$DP[,sample]
    dp_bulk = geno(vcf_rest)$DP[,bulk]
    gq_sample = geno(vcf_rest)$GQ[,sample]
    gq_bulk = geno(vcf_rest)$GQ[,bulk]
    keep_f = vaf >= 0.3 & dp_sample >= 20 & dp_bulk >= 20 & gq_sample >= 99 & gq_bulk >= 99
    vcf_rest = vcf_rest[keep_f]
    
    vcf = rbind(vcf_tri, vcf_x, vcf_rest)
    
    #Remove muts that are in the shared part of the tree
    if (file.exists(vcf_shared_filtered_fname)){
        vcf_shared = readVcf(vcf_shared_filtered_fname, genome = "hg19")
        not_unique_clone = findOverlaps(granges(vcf), granges(vcf_shared)) %>% queryHits() %>% unique()
        if (length(not_unique_clone) > 0){
            vcf = vcf[-not_unique_clone]
        }
    }
    
    if(nrow(vcf) > 0){
        vcf = sort_vcf(vcf)
        trash = vcf_out(vcf, vcf_out_fname, bed_validated_fname, bed_tovalidate_fname)
    }
    return(0)
})







# ####__________________Filter muts based on igv validation____________####
# #Shared (Can't be simplified, because some muts need to be manually added)
# for (fetus in fetuses){
#     sample_name = fetus[[1]]
#     bulk = fetus[[2]]
#     gender = fetus[[3]]
#     out_dir = paste0(out_dir_base, sample_name, "/shared/")
#     
#     vcf_shared_filtered_fname = paste0(out_dir, sample_name, "_shared.vcf")
#     bed_fname = paste0(out_dir, sample_name, "_shared_validated.bed")
#     
#     if (!file.exists(vcf_shared_filtered_fname) | !file.exists(bed_fname)){
#         print(paste0("A shared.vcf or shared_validated.bed file does not exist for sample: ", sample_name))
#         next
#     }
#     vcf_shared_filtered = readVcf(vcf_shared_filtered_fname)
#     bed = read_tsv(bed_fname, col_types = cols(CHROM = "c"))
#     bed %<>% dplyr::filter(Real_mut) %>% makeGRangesFromDataFrame()
#     vcf_ranges = granges(vcf_shared_filtered)
#     true_index = findOverlaps(vcf_ranges, bed) %>% queryHits()
#     vcf_shared_filtered_true = vcf_shared_filtered[true_index]
#     
#     #Some mutations were put in a clone specific vcf, but are actually shared. These are added here.
#     vcf_manual_addition_fname = paste0(out_dir, sample_name, "_shared_manualadded.vcf")
#     if (file.exists(vcf_manual_addition_fname)){
#         vcf_manual_add = readVcf(vcf_manual_addition_fname)
#         
#         already_in_vcf = findOverlaps(granges(vcf_manual_add), granges(vcf_shared_filtered_true)) %>% queryHits() %>% unique() #Check if the mutations that are being added havent been added by a change in the filtering.
#         if (length(already_in_vcf) > 0){
#             vcf_manual_add = vcf_manual_add[-already_in_vcf]
#         }
#         vcf_shared_filtered_true = rbind(vcf_shared_filtered_true, vcf_manual_add)
#         vcf_shared_filtered_true = sort_vcf(vcf_shared_filtered_true)
#         
#     }
#     writeVcf(vcf_shared_filtered_true, paste0(out_dir, sample_name, "_shared_true.vcf"))
# }
# 
# #Unique in bulk
# overview_samples_tofilter = overview_samples[overview_samples$uniq_vcf_inbulk_unfiltered_true, ]
# vcf_fnames = overview_samples_tofilter$uniq_vcf_inbulk_unfiltered
# bed_fnames = paste0(out_dir_base, overview_samples_tofilter$fetus, "/dirichlet/", overview_samples_tofilter$sample, "_inbulk_validated.bed")
# vcf_filtered_fnames = overview_samples_tofilter$uniq_vcf_inbulk
# 
# for(i in seq(length(vcf_fnames))){
#     igv_filter(vcf_fnames[i], bed_fnames[i], vcf_filtered_fnames[i])
# }
# 
# #Unique notinbulk
# overview_samples_tofilter = overview_samples[overview_samples$uniq_vcf_notbulk_unfiltered_true, ]
# vcf_fnames = overview_samples_tofilter$uniq_vcf_notbulk_unfiltered
# bed_fnames = paste0(out_dir_base, overview_samples_tofilter$fetus, "/dirichlet/", overview_samples_tofilter$sample, "_notbulk_validated.bed")
# vcf_filtered_fnames = overview_samples_tofilter$uniq_vcf_notbulk
# for(i in seq(length(vcf_fnames))){
#     igv_filter(vcf_fnames[i], bed_fnames[i], vcf_filtered_fnames[i])
# }
# 
# #update overview samples
# overview_samples = create_overview_samples(fetuses, out_dir_base)
# 
# 
# ####_____________________Combine the bulk and not bulk mutations, for the unique mutations_____####
# for (i in seq(1, nrow(overview_samples))){
#     overview_sample = overview_samples[i,]
#     if(overview_sample$uniq_vcf_inbulk_true & overview_sample$uniq_vcf_notbulk_true){
#         inbulk = readVcf(overview_sample$uniq_vcf_inbulk, genome = "hg19")
#         notbulk = readVcf(overview_sample$uniq_vcf_notbulk, genome = "hg19")
#         vcf = rbind(inbulk, notbulk)
#         vcf = sort_vcf(vcf)
#     } else if (overview_sample$uniq_vcf_inbulk_true){
#         vcf = readVcf(overview_sample$uniq_vcf_inbulk, genome = "hg19")
#     } else if (overview_sample$uniq_vcf_notbulk_true){
#         vcf = readVcf(overview_sample$uniq_vcf_notbulk, genome = "hg19")
#     } else{
#         print(paste0("No unique mutations exist for sample: ", overview_sample$sample, " in fetus ", overview_sample$fetus))
#         next
#     }
#     writeVcf(vcf, overview_sample$unique_vcf)
# }
# #update overview samples
# overview_samples = create_overview_samples(fetuses, out_dir_base)
# 
# ####___________Split true shared mutation between bulk and not bulk.___________####
# for (fetus in fetuses){
#     sample_name = fetus[[1]]
#     bulk = fetus[[2]]
#     gender = fetus[[3]]
#     out_dir = paste0(out_dir_base, sample_name, "/shared/")
#     
#     vcf_fname = paste0(out_dir, sample_name, "_shared_true.vcf")
#     if (!file.exists(vcf_fname)){
#         next
#     }
#     vcf = readVcf(vcf_fname)
#     
#     #Write a vcf with all informative sites, that are present in the bulk
#     vaf = get_vaf(vcf, bulk)
#     vaf_f = vaf > 0
#     vcf_inbulk = vcf[vaf_f, ]
#     writeVcf(vcf_inbulk, paste0(out_dir, sample_name, "_shared_inbulk.vcf"))
#     
#     #Write a vcf with shared sites that are not present in the bulk
#     vcf_notbulk = vcf[!vaf_f,]
#     writeVcf(vcf_notbulk, paste0(out_dir, sample_name, "_shared_notbulk.vcf"))
# }
# #update overview samples
# overview_samples = create_overview_samples(fetuses, out_dir_base)


