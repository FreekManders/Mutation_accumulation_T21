library(gdata)
library(tidyverse)
library(VariantAnnotation)

wdir = "~/hpc/pmc_vanboxtel/projects/Freek_germline_AMLALL/"
pre_dir = "~/surfdrive/Shared/Boxtel_General/Data/Predisposition/"
samples = c("HMFreg0312_PMC22813AML", "ALL")

#Create different gene lists. The Zhang gene list contains general cancer genes. The ad and ar lists are subsets of this list, containing genes that are predisposition associated.
Zhang_genes_file = read.xls(paste0(pre_dir, "germline_mutations_predisposition.xlsx"), sheet = "TableS2", stringsAsFactors = F, skip = 1)
Zhang_genes = Zhang_genes_file %>% dplyr::select(Gene) %>% unique()

genes_ad = read_tsv(paste0(pre_dir, "Predisposition_genes_AD.txt"))
genes_ad_query = paste("\\|", genes_ad$pre_gene_ad, "\\|", collapse = "|", sep = "")

genes_ar = read_tsv(paste0(pre_dir, "Predisposition_genes_AR.txt"))
genes_ar_query = paste("\\|", genes_ar$pre_gene_ar, "\\|", collapse = "|", sep = "")

other_genes = Zhang_genes %>% dplyr::filter(!Gene %in% genes_ad$pre_gene_ad & !Gene %in% genes_ar$pre_gene_ar)
other_genes_query = paste("\\|", other_genes$Gene, "\\|", collapse = "|", sep = "")

DGR_genes_file = read.xls(paste0(pre_dir, "2018_DGR_paper_predisposition_genes.xlsx"), stringsAsFactors = F)
DGR_genes_query = paste("\\|", DGR_genes_file$Gene, "\\|", collapse = "|", sep = "")

#Combine Zhang and DGR lists
all_genes = c(DGR_genes_file$Gene, Zhang_genes$Gene) %>% unique()
all_genes_query = paste("\\|", all_genes, "\\|", collapse = "|", sep = "")

for (sample in samples){
    #Read in vcf and collapse the annotation field, so it can be grepped.
    vcf_fname = paste0(wdir, sample, "_highimpact.vcf")
    vcf = readVcf(vcf_fname)
    ANN = info(vcf)$ANN %>% paste(., collapse = "__")
    
    #Filter vcf on all cancer genes
    muts_index = grepl(all_genes_query, ANN)
    muts = vcf[muts_index]
    ANN_cancer = info(muts)$ANN %>% paste(., collapse = "__")
    
    #Add column to show which part of the Zhang list a mut was linked to
    header_add = DataFrame("Number" = 1, "Type" = "String", "Description" = "The group of predisposition genes this snv was linked to", row.names = "Matching_genelist_Zhang")
    info(header(muts)) = rbind(info(header(muts)), header_add)
    info(muts)$Matching_genelist_Zhang = "NA"
    
    other_genes_index = grepl(other_genes_query, ANN_cancer)
    info(muts)$Matching_genelist_Zhang[other_genes_index] = "Cancer_related"
    genes_ad_index = grepl(genes_ad_query, ANN_cancer)
    info(muts)$Matching_genelist_Zhang[genes_ad_index] = "Predisposition_AD"
    genes_ar_index = grepl(genes_ar_query, ANN_cancer)
    info(muts)$Matching_genelist_Zhang[genes_ar_index] = "Predisposition_AR"
    
    ##Add column to show whether a muts was linked to the DGR list
    header_add = DataFrame("Number" = 1, "Type" = "String", "Description" = "The group of predisposition genes this snv was linked to", row.names = "Matching_genelist_DGR")
    info(header(muts)) = rbind(info(header(muts)), header_add)
    info(muts)$Matching_genelist_DGR = "NA"
    
    DGR_genes_index = grepl(DGR_genes_query, ANN_cancer)
    info(muts)$Matching_genelist_DGR[DGR_genes_index] = "Predisposition"
    
    meta_head = meta(header(muts))
    first_line = grep("fileformat", names(meta_head))
    meta(header(muts)) = c(meta_head[first_line], meta_head[-first_line])
    
    #write out vcf
    writeVcf(muts, paste0(wdir, sample, "_predisposition.vcf"))
    print(paste0("Wrote out the vcf for sample: ", sample))
}