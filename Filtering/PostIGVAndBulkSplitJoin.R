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
fetuses_tb = read_tsv(paste0(out_dir_base, "fetuses_ewart.txt"), col_types = cols(.default = "c"))
fetuses = split(fetuses_tb, seq(1,nrow(fetuses_tb)))
fetuses = lapply(fetuses, unlist)
overview_samples = create_overview_samples(fetuses, out_dir_base)



#################################SNVS##################################
####__________________Filter snvs based on igv validation____________####
#Shared (Can't be simplified, because some muts need to be manually added)
for (fetus in fetuses){
    sample_name = fetus[[1]]
    bulk = fetus[[2]]
    gender = fetus[[3]]
    out_dir = paste0(out_dir_base, sample_name, "/shared/")
    
    vcf_shared_filtered_fname = paste0(out_dir, sample_name, "_shared.vcf")
    bed_fname = paste0(out_dir, sample_name, "_shared_validated.bed")
    
    if (!file.exists(vcf_shared_filtered_fname) | !file.exists(bed_fname)){
        print(paste0("A shared.vcf or shared_validated.bed file does not exist for sample: ", sample_name))
        next
    }
    vcf_shared_filtered = readVcf(vcf_shared_filtered_fname)
    bed = read_tsv(bed_fname, col_types = cols(CHROM = "c"))
    bed %<>% dplyr::filter(Real_mut) %>% makeGRangesFromDataFrame(starts.in.df.are.0based = T)
    vcf_ranges = granges(vcf_shared_filtered)
    true_index = findOverlaps(vcf_ranges, bed) %>% queryHits()
    vcf_shared_filtered_true = vcf_shared_filtered[true_index]
    
    #Some mutations were put in a clone specific vcf, but are actually shared. These are added here.
    vcf_manual_addition_fname = paste0(out_dir, sample_name, "_shared_manualadded.vcf")
    if (file.exists(vcf_manual_addition_fname)){
        vcf_manual_add = readVcf(vcf_manual_addition_fname)
        
        already_in_vcf = findOverlaps(granges(vcf_manual_add), granges(vcf_shared_filtered_true)) %>% queryHits() %>% unique() #Check if the mutations that are being added havent been added by a change in the filtering.
        if (length(already_in_vcf) > 0){
            vcf_manual_add = vcf_manual_add[-already_in_vcf]
        }
        vcf_shared_filtered_true = rbind(vcf_shared_filtered_true, vcf_manual_add)
        vcf_shared_filtered_true = sort_vcf(vcf_shared_filtered_true)
        
    }
    writeVcf(vcf_shared_filtered_true, paste0(out_dir, sample_name, "_shared_true.vcf"))
}

#Unique in bulk
overview_samples_tofilter = overview_samples[overview_samples$uniq_vcf_inbulk_unfiltered_true, ]
vcf_fnames = overview_samples_tofilter$uniq_vcf_inbulk_unfiltered
bed_fnames = paste0(out_dir_base, overview_samples_tofilter$fetus, "/dirichlet/", overview_samples_tofilter$sample, "_inbulk_validated.bed")
vcf_filtered_fnames = overview_samples_tofilter$uniq_vcf_inbulk

for(i in seq(length(vcf_fnames))){
    igv_filter(vcf_fnames[i], bed_fnames[i], vcf_filtered_fnames[i])
}

#Unique notinbulk
overview_samples_tofilter = overview_samples[overview_samples$uniq_vcf_notbulk_unfiltered_true, ]
vcf_fnames = overview_samples_tofilter$uniq_vcf_notbulk_unfiltered
bed_fnames = paste0(out_dir_base, overview_samples_tofilter$fetus, "/dirichlet/", overview_samples_tofilter$sample, "_notbulk_validated.bed")
vcf_filtered_fnames = overview_samples_tofilter$uniq_vcf_notbulk
for(i in seq(length(vcf_fnames))){
    igv_filter(vcf_fnames[i], bed_fnames[i], vcf_filtered_fnames[i])
}

#In bulk not in clones
overview_samples_tofilter = overview_samples[overview_samples$vcf_inbulk_not_clones_unfiltered_true, ] %>% dplyr::filter(!duplicated(fetus))
vcf_fnames = overview_samples_tofilter$vcf_inbulk_not_clones_unfiltered
bed_fnames = paste0(out_dir_base, overview_samples_tofilter$fetus, "/shared/", overview_samples_tofilter$fetus, "_subclonalbulks_notclones_validated.bed")
vcf_filtered_fnames = overview_samples_tofilter$vcf_inbulk_not_clones
for(i in seq(length(vcf_fnames))){
    igv_filter(vcf_fnames[i], bed_fnames[i], vcf_filtered_fnames[i])
    
    #Get the location and vaf of these muts, so they can be used to make primers.
    if (file.exists(vcf_filtered_fnames[i])){
        vcf = readVcf(vcf_filtered_fnames[i], genome = "hg19")
        bulk = overview_samples_tofilter$bulk[i]
        vaf = get_vaf(vcf, bulk)
        out = tibble("Chrom" = as.vector(seqnames(vcf)), "Coord" = start(vcf), "VAF" = vaf)
        write_tsv(out, paste0(out_dir_base, overview_samples_tofilter$fetus[i], "/shared/", overview_samples_tofilter$fetus[i], "_subclonalbulks_notclones_loc.txt"))
    }
}


#update overview samples
overview_samples = create_overview_samples(fetuses, out_dir_base)


####_____________________Combine the bulk and not bulk mutations, for the unique snvs_____####
for (i in seq(1, nrow(overview_samples))){
    overview_sample = overview_samples[i,]
    print(overview_sample$sample)
    if(overview_sample$uniq_vcf_inbulk_true & overview_sample$uniq_vcf_notbulk_true){
        inbulk = readVcf(overview_sample$uniq_vcf_inbulk, genome = "hg19")
        notbulk = readVcf(overview_sample$uniq_vcf_notbulk, genome = "hg19")
        vcf = rbind(inbulk, notbulk)
        vcf = sort_vcf(vcf)
    } else if (overview_sample$uniq_vcf_inbulk_true){
        vcf = readVcf(overview_sample$uniq_vcf_inbulk, genome = "hg19")
    } else if (overview_sample$uniq_vcf_notbulk_true){
        vcf = readVcf(overview_sample$uniq_vcf_notbulk, genome = "hg19")
    } else{
        print(paste0("No unique mutations exist for sample: ", overview_sample$sample, " in fetus ", overview_sample$fetus))
        next
    }
    writeVcf(vcf, overview_sample$unique_vcf)
}
#update overview samples
overview_samples = create_overview_samples(fetuses, out_dir_base)

####___________Split true shared snvs between bulk and not bulk.___________####
for (fetus in fetuses){
    sample_name = fetus[[1]]
    bulk = fetus[[2]]
    gender = fetus[[3]]
    out_dir = paste0(out_dir_base, sample_name, "/shared/")
    
    vcf_fname = paste0(out_dir, sample_name, "_shared_true.vcf")
    if (!file.exists(vcf_fname)){
        next
    }
    vcf = readVcf(vcf_fname)
    
    #Write a vcf with all informative sites, that are present in the bulk
    vaf = get_vaf(vcf, bulk)
    vaf_f = vaf > 0
    vcf_inbulk = vcf[vaf_f, ]
    writeVcf(vcf_inbulk, paste0(out_dir, sample_name, "_shared_inbulk.vcf"))
    
    #Write a vcf with shared sites that are not present in the bulk
    vcf_notbulk = vcf[!vaf_f,]
    writeVcf(vcf_notbulk, paste0(out_dir, sample_name, "_shared_notbulk.vcf"))
}
#update overview samples
overview_samples = create_overview_samples(fetuses, out_dir_base)

################################INDELS######################
####____________________Filter indels based on IGV______________####
#shared muts
overview_fetuses = overview_samples %>% dplyr::filter(!duplicated(fetus))
trash = apply(overview_fetuses, 1, function(row){
    bed_validated = row[["indel_shared_validated"]]
    vcf_fname = row[["indel_shared_filtered"]]
    vcf_out_fname = row[["indel_shared"]]
    igv_filter(vcf_fname, bed_validated, vcf_out_fname)
    
})

#Add shared indels that were manually identified
overview_fetuses = overview_samples %>% dplyr::filter(!duplicated(fetus) & indel_shared_manualadded_true)
trash = apply(overview_fetuses, 1, function(row){
    print(row[["indel_shared_manualadded"]])
    add_manual = readVcf(row[["indel_shared_manualadded"]], genome = "hg19")
    vcf_shared = readVcf(row[["indel_shared"]], genome = "hg19")
    vcf = rbind(add_manual, vcf_shared) %>% sort_vcf()
    writeVcf(vcf, row[["indel_shared"]])
})

#unique muts
trash = apply(overview_samples, 1, function(row){
    bed_notbulk_validated = row[["indel_uniq_notbulk_bed_validated"]]
    vcf_notbulk_fname = row[["indel_uniq_notbulk_filtered"]]
    vcf_notbulk_out = row[["indel_uniq_notbulk"]]
    print(vcf_notbulk_fname)
    igv_filter(vcf_notbulk_fname, bed_notbulk_validated, vcf_notbulk_out)
    
    bed_inbulk_validated = row[["indel_uniq_inbulk_bed_validated"]]
    vcf_inbulk_fname = row[["indel_uniq_inbulk_filtered"]]
    vcf_inbulk_out = row[["indel_uniq_inbulk"]]
    print(vcf_inbulk_fname)
    igv_filter(vcf_inbulk_fname, bed_inbulk_validated, vcf_inbulk_out)
    
    #Add manual added for inbulk indels
    if (row[["indel_uniq_inbulk_manualadded_true"]]){
        
        if (file.exists(row[["indel_uniq_inbulk"]])){
            vcf_inbulk = readVcf(row[["indel_uniq_inbulk"]], genome = "hg19")
            add_manual = readVcf(row[["indel_uniq_inbulk_manualadded"]], genome = "hg19")
            vcf = rbind(vcf_inbulk, add_manual) %>% sort_vcf()
        } else{
            vcf = readVcf(row[["indel_uniq_inbulk_manualadded"]], genome = "hg19")
        }
        writeVcf(vcf, row[["indel_uniq_inbulk"]])
    }
    
})
#update overview samples
overview_samples = create_overview_samples(fetuses, out_dir_base)

####_____________________Combine the bulk and not bulk mutations, for the unique mutations_____####
for (i in seq(1, nrow(overview_samples))){
    overview_sample = overview_samples[i,]
    if(overview_sample$indel_uniq_inbulk_true & overview_sample$indel_uniq_notbulk_true){
        inbulk = readVcf(overview_sample$indel_uniq_inbulk, genome = "hg19")
        notbulk = readVcf(overview_sample$indel_uniq_notbulk, genome = "hg19")
        vcf = rbind(inbulk, notbulk)
        vcf = sort_vcf(vcf)
    } else if (overview_sample$indel_uniq_inbulk_true){
        vcf = readVcf(overview_sample$indel_uniq_inbulk, genome = "hg19")
    } else if (overview_sample$indel_uniq_notbulk_true){
        vcf = readVcf(overview_sample$indel_uniq_notbulk, genome = "hg19")
    } else{
        print(paste0("No unique indels exist for sample: ", overview_sample$sample, " in fetus ", overview_sample$fetus))
        next
    }
    writeVcf(vcf, overview_sample$indel_uniq)
}
#update overview samples
overview_samples = create_overview_samples(fetuses, out_dir_base)

####___________Split the shared mutations between inbulk and notbulk_____####
overview_fetuses = overview_samples %>% dplyr::filter(!duplicated(fetus) & indel_shared_true)
trash = apply(overview_fetuses, 1, function(row){
    vcf_fname = row[["indel_shared"]]
    bulk = row[["bulk"]]
    vcf_out_inbulk = row[["indel_shared_inbulk"]]
    vcf_out_notbulk = row[["indel_shared_notbulk"]]
    vcf = readVcf(vcf_fname, genome = "hg19")
    vaf = get_vaf(vcf, bulk)
    vaf_f = vaf > 0
    vcf_inbulk = vcf[vaf_f,]
    if (nrow(vcf_inbulk) > 0)
        writeVcf(vcf_inbulk, vcf_out_inbulk)
    vcf_notbulk = vcf[!vaf_f,]
    if (nrow(vcf_notbulk) > 0)
        writeVcf(vcf_notbulk, vcf_out_notbulk)
    
})



