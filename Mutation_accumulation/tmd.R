library(tidyverse)
library(VariantAnnotation)
library(GenomicRanges)
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg19)
library(magrittr)
library(nlme)
library(ggdendro)
library(lmvar)
library(cowplot)
library(ggbeeswarm)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/vcfFilterTree_functions.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


owndata_dir = "~/hpc/pmc_vanboxtel/projects/Freek_TMD/"

#Read in published data and pre-processing
tmd_df = read_delim("~/surfdrive/Shared/Projects/Ruben/Trisomie21/muts2.txt", delim = " ", col_names = F)
colnames(tmd_df) <- c("patient","disease","gene","location","chr","pos","REF","ALT")

tmd_df %<>% dplyr::filter(REF != "-" & ALT != "-" & nchar(REF) == 1 & nchar(ALT) == 1) #Remove indels
tmd_df %<>% dplyr::filter(disease != "TAM" | patient == "003") #get clonal in tam, by taking muts occuring in both tam and amkl
tmd_df %<>% mutate(disease = gsub("TAM/AMKL", "TAM", disease), chr = gsub("Chr", "chr", chr)) #Set muts occuring in both to the tam, as they originated there.

gr_published = makeGRangesFromDataFrame(tmd_df, keep.extra.columns = T, start.field = "pos", end.field = "pos")
gr_published$REF = DNAStringSet(gr_published$REF)
gr_published$ALT = lapply(gr_published$ALT, DNAStringSet) %>% DNAStringSetList()
chroms = paste0("chr", c(1:22, "X"))
seqlevels(gr_published, pruning.mode = "coarse") = chroms

gr_tam = gr_published[gr_published$disease == "TAM"]
grl_tam = split(gr_tam, gr_tam$patient)


#Read in own data
fetuses_tb = read_tsv(paste0(owndata_dir, "fetuses.txt"), col_types = cols(.default = "c"))
fetuses = split(fetuses_tb, seq(1,nrow(fetuses_tb)))
fetuses = lapply(fetuses, unlist)
overview_samples = create_overview_samples(fetuses, owndata_dir)
overview_samples %<>% dplyr::filter(uniq_vcf_notbulk_true)
grl_owndata = read_vcfs_as_granges(overview_samples$uniq_vcf_notbulk, overview_samples$sample, ref_genome)

for (i in 1:length(grl_owndata)){
    grl_owndata[[i]]$patient = names(grl_owndata)[i]
}

gr_owndata = unlist(grl_owndata)
gr_owndata$disease = "TAM"


#combine the data
gr = c(gr_owndata[,c("patient", "disease", "REF", "ALT")], gr_published[,c("patient", "disease", "REF", "ALT")])
grl = split(gr, gr$disease)

gr_tmd = grl$TAM
saveRDS(gr_tmd, "~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/tmd_gr.rds")

#Count muts
counts_tmd = grl$TAM$patient %>% table() %>% enframe(name = "patient", value = "total_muts")

ages_patients = tibble("patient" = paste0("00", 1:4), "age_diag_tam" = rep(0.728, 4), "age_diag_amkl" = c(2.89, 3.23, 3.31, 1.23)) %>% mutate(diag_diff = age_diag_amkl - age_diag_tam)

chroms_auto = paste0("chr", 1:22)
grl_auto = grl
seqlevels(grl_auto, pruning.mode = "fine") = chroms_auto
counts_tmd_auto = grl_auto$TAM$patient %>% table() %>% enframe(name = "patient", value = "total_muts")
counts_tmd_auto %<>% dplyr::mutate(total_muts = ifelse(patient == "003", total_muts - 7, total_muts)) #remove subclones in patient 003
#counts = tmd_df %>% dplyr::group_by(patient, disease) %>% dplyr::summarise(count = n())



nr_muts_tb = read_tsv("~/hpc/pmc_vanboxtel/projects/Freek_trees/count_muts_extended.txt")
counts = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si", sample)) %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2")) %>% dplyr::select(Donor = fetus_name, total_muts = nr_total_snv_notbulk_auto, Age = age_year)
weeksbyyear = 52.177457
age_at_birth = 38 / weeksbyyear
#muts_counts_axel = read_tsv(mut_counts_axel_fname) %>% unique() %>% dplyr::filter(Tissue == "Blood") %>% dplyr::select(person = Donor, Age, norm_snvs = norm_muts) %>% mutate(Age = Age + age_at_birth)
combined_counts_all = counts
combined_counts_all %<>% mutate("age_cat" = ifelse(Age < 0.6, "fetal", "Post-infant"))
combined_counts_all$age_cat = factor(combined_counts_all$age_cat, levels = c("Post-infant", "fetal"))


combined_counts = combined_counts_all


counts_tmd = grl$TAM$patient %>% table() %>% enframe(name = "Donor", value = "total_muts") %>% dplyr::mutate(Type = "TMD")
counts_fetuses = nr_muts_tb %>% dplyr::filter(trisomy == "21" & celltype == "HSC") %>% dplyr::select(Donor = fetus_name, total_muts = nr_total_snv_notbulk) %>% dplyr::mutate(Type = "T21 Fetal clones")
all_counts = rbind(counts_tmd, counts_fetuses) %>% 
    dplyr::mutate(total_muts = as.integer(total_muts))

m1 = lme(total_muts ~ Type, random = ~ 1 | Donor, data = all_counts)
pval = m1 %>% 
    get_fixed() %>% 
    dplyr::filter(variable == "TypeTMD") %>% 
    pull(`p-value`)

tmd_vs_t21_fig = ggplot(all_counts, aes(y = total_muts, x = Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_quasirandom(aes(colour = Donor), size = 3) +
    labs(title = "TMD vs T21 fetal clones", x = "", y = "Nr. of corrected mutations", colour = "Donor") +
    theme_bw() +
    theme(text = element_text(size = 24)) +
    annotate("text", x = 1, y = 120, label = paste0("P: ", round(pval, 3)), size = 6)

pdf(paste0(owndata_dir, "tmd_vs_t21.pdf"), width = 10)
tmd_vs_t21_fig
dev.off()


#Create complete vcf and bed
for (fetus in unique(overview_samples$fetus)){
    overview_fetus = overview_samples %>% dplyr::filter(fetus == !!fetus)
    vcf_fnames1 = overview_fetus$indel_uniq
    vcf_fnames2 = overview_fetus$indel_shared
    vcf_fnames3 = overview_fetus$unique_vcf
    vcf_fnames4 = overview_fetus$shared_vcf
    vcf_fnames = c(vcf_fnames1, vcf_fnames2, vcf_fnames3, vcf_fnames4) %>% unique
    vcf_fnames = vcf_fnames[file.exists(vcf_fnames)]
    
    vcf_l = lapply(vcf_fnames, function(x) readVcf(x, genome = "hg19"))
    
    vcf = do.call(rbind, vcf_l) %>%
        sort() %>%
        unique()
    writeVcf(vcf, paste0(owndata_dir, fetus, "/", fetus, "_complete.vcf"))
    gr = granges(vcf)
    bed = tibble("CHROM" = as.vector(seqnames(gr)), "START" = start(gr)-1, "END" = end(gr))
    bed_fname = paste0(owndata_dir, fetus, "/", fetus, "_complete.bed")
    writeLines(paste0("#", paste0(colnames(bed), collapse = "\t")), bed_fname)
    write.table(bed, bed_fname, append = T, quote = F, sep = "\t", col.names = F, row.names = F)
}

#Create an exon table
get_exon_table = function(vcf){
    ann = info(vcf)$ANN %>%
        as.list()
    ann[elementNROWS(ann) == 0] = ""
    ann = purrr::map_chr(ann, str_c, collapse = ";")
    coding_f = str_detect(ann, "synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant")
    ann_exon = ann[coding_f]
    ann_l = str_split(ann_exon, pattern = "\\|")
    effect_type = purrr::map_chr(ann_l, 2) #Missense, nonnsense ect.
    effect = purrr::map_chr(ann_l, 3)
    gene = purrr::map_chr(ann_l, 4)
    gr = granges(vcf)
    gr$ALT = gr$ALT %>% unlist() %>% as.vector()
    tb = gr %>% 
        as.data.frame() %>% 
        as_tibble() %>% 
        dplyr::filter(coding_f) %>% 
        dplyr::select(seqnames, start, end, REF, ALT) %>% 
        dplyr::mutate(gene = gene, effect = effect, effect_type = effect_type)
    return(tb)
}

fetuses = unique(overview_samples$fetus)
vcf_fnames = str_c("~/hpc/pmc_vanboxtel/projects/Freek_TMD/", fetuses, "/", fetuses, "_complete.vcf")
vcf_l = purrr::map(vcf_fnames, readVcf)
exon_tb = purrr::map(vcf_l, get_exon_table) %>% 
    set_names(fetuses) %>% 
    bind_rows(.id = "fetus")
write_tsv(exon_tb, str_c(owndata_dir, "exon_muts.txt"))
