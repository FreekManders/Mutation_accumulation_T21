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
library(ggpubr)
library(nlme)
library(lsmeans)
library(ape)
library(PCDimension)
library(ComplexHeatmap)
library(ggsignif)
library(lme4)
library(ggeffects)
library(lmvar)
library(extrafont)
loadfonts()

set.seed(42)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/vcfFilterTree_functions.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/Do_mutationalpatterns.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/id_context.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/dbs_context.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/lmm_leave_one_out_leverage.R")
genome = "BSgenome.Hsapiens.UCSC.hg19"

#mut_counts_axel_fname = "~/surfdrive/Shared/Projects/Zzz_sleeping_projects/Healthy_BMseq/Frequencies/Table_mutation_counts_total.txt"
muts_counts_axel_fname = "~/surfdrive/Shared/Projects/AgeLine/data/Healthy_BM_Freek.txt"
muts_counts_arianne_fname = "~/surfdrive/Shared/Projects/AgeLine/data/Arianne_AML_ALL_HMFreg0372.txt"
muts_counts_mirjam_fname = "~/surfdrive/Shared/Projects/Mirjam_HSCT/Data/Clone_data.txt"
muts_counts_si_fname = "~/surfdrive/Shared/Projects/AgeLine/data/SI_data_blokzijl.txt"
vcfs_axel_dir = "~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/SNVs/Healthy_bone_marrow/"
vcfs_blokzijl_dir = "~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/SNVs/Published_data/Co_SI_Li_ASCs_Blokzijl_et_al_2016_Nature/intestine/"
callable_axel_dir = "~/surfdrive/Shared/Boxtel_General/Data/callableLoci/Healthy_bone_marrow/"
cosmic_sig_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_working_signatures.txt"
telomere_dir_fname = "~/hpc/pmc_vanboxtel/projects/Freek_trees/telomeres/telomerecat/"

out_dir_base = "~/hpc/pmc_vanboxtel/projects/Freek_trees/"
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

fetuses_tb = read_tsv(paste0(out_dir_base, "fetuses_inclewart.txt"), col_types = cols(.default = "c"))
fetuses = split(fetuses_tb, seq(1,nrow(fetuses_tb)))
fetuses = lapply(fetuses, unlist)
overview_samples = create_overview_samples(fetuses, out_dir_base)



#overview_samples = read_tsv(paste0(out_dir_base, "overview_samples.txt"))


####___________Create figures and usefull files for validated shared mutations.___________####

for (fetus in fetuses){
    sample_name = fetus[[1]]
    bulk = fetus[[2]]
    gender = fetus[[3]]
    out_dir = paste0(out_dir_base, sample_name, "/shared/")
    other_bulks = fetus[[6]]
    
    vcf_fname = paste0(out_dir, sample_name, "_shared_true.vcf")
    if (!file.exists(vcf_fname)){
        next
    }
    vcf = readVcf(vcf_fname)
    
    all_bulks = get_allbulks(other_bulks, bulk)
    
    #Create heatmaps of how often mutations cooccur between samples
    gt = geno(vcf)$GT
    bulkcols = grep(all_bulks, colnames(gt))
    gt_clones = gt %>% as_tibble() %>% dplyr::select(-bulkcols)
    co_occur_clones_fig = co_occurence_map(gt_clones)
    
    vaf = get_vaf(vcf, bulk)
    vaf_f = vaf > 0
    gt_notbulk = gt %>% as_tibble() %>% dplyr::filter(!vaf_f) %>% dplyr::select(-bulkcols)
    co_occur_notbulk_fig = co_occurence_map(gt_notbulk)
    gt_inbulk = gt %>% as_tibble() %>% dplyr::filter(vaf_f) %>% dplyr::select(-bulkcols)
    co_occur_inbulk_fig = co_occurence_map(gt_inbulk)
    
    pdf(paste0(out_dir, sample_name, "_co_occurence.pdf"), width = 11.5, height = 8)
    print(co_occur_clones_fig)
    print(co_occur_notbulk_fig)
    print(co_occur_inbulk_fig)
    dev.off()
    
    
    #Count how often combinations of mutations occur
    count_combis(vcf, paste0(out_dir, sample_name, "_shared"), bulk, all_bulks)
    
    #Create a preliminary tree based on the shared mutations.
    co_occur_m = gt_clones
    co_occur_m[co_occur_m == "0/0"] = 0
    co_occur_m[co_occur_m == "0/1" | co_occur_m == "1/1"] = 1
    co_occur_m %<>% mutate_all(list(~as.integer)) %>% as.matrix()
    tree = co_occur_m %>% t() %>% dist(method = "binary") %>% hclust(method = "ward.D2")
    
    pdf(paste0(out_dir, sample_name, "_tree.pdf"))
    plot(tree)
    dev.off()
}


####_______Look at the surveyed region of the different samples (callable loci)_______#####

surveyed_tb = overview_samples %>% dplyr::filter(!is.na(surveyed_region))

surveyed_fig = ggplot(surveyed_tb, aes(x = fetus, y = surveyed_region_pct, fill = origin)) +
    geom_boxplot(aes(color = fetus), outlier.shape = NA, fill = NA) +
    geom_quasirandom(shape = 21, size = 3) +
    #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.001, dotsize = 2) +
    scale_color_discrete(guide=FALSE) +
    labs(title = "Surveyed region of the genome", x = "Fetus", y = "Surveyed region of the genome") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"), axis.text.x = element_text(angle = 90))

surveyed_fullaxis_fig = ggplot(surveyed_tb, aes(x = fetus, y = surveyed_region_pct, fill = origin)) +
    geom_boxplot(aes(color = fetus), outlier.shape = NA, fill = NA) +
    #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.01, dotsize = 2) +
    geom_quasirandom(shape = 21, size = 3) +
    scale_color_discrete(guide=FALSE) +
    labs(title = "Surveyed region of the genome", x = "Fetus", y = "Surveyed region of the genome") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"), axis.text.x = element_text(angle = 90)) +
    coord_cartesian(ylim = c(0,1), expand = F)

pdf(paste0(out_dir_base, "surveyed_region.pdf"))
surveyed_fig
surveyed_fullaxis_fig
dev.off()

####______________Count mutations hscs fetuses_______________####
count_muts = function(row){
    sample = row[["sample"]]
    fetus_name = row[["fetus"]]
    shared_vcf_fname = row[["shared_vcf"]]
    bulk = row[["bulk"]]
    if (file.exists(shared_vcf_fname)){
        shared_vcf = readVcf(shared_vcf_fname)
        gt = geno(shared_vcf)$GT
        vaf_bulk = get_vaf(shared_vcf, bulk)
        nr_shared = gt[,sample] %in% c("0/1", "1/1") %>% sum()
        #nr_shared_inbulk_f = gt[,sample] %in% c("0/1", "1/1") & gt[,bulk] %in% c("0/1", "1/1")
        nr_shared_inbulk_f = gt[,sample] %in% c("0/1", "1/1") & vaf_bulk > 0
        nr_shared_inbulk = sum(nr_shared_inbulk_f)
        
        seqs = granges(shared_vcf) %>% seqnames() %>% as.vector()
        shared_auto_vcf = shared_vcf[seqs %in% 1:22]
        gt_auto = geno(shared_auto_vcf)$GT
        vaf_auto_bulk = get_vaf(shared_auto_vcf, bulk)
        nr_shared_auto = gt_auto[,sample] %in% c("0/1", "1/1") %>% sum()
        nr_shared_inbulk_f_auto = gt_auto[,sample] %in% c("0/1", "1/1") & vaf_auto_bulk > 0
        nr_shared_inbulk_auto = sum(nr_shared_inbulk_f_auto)
        
    } else{
        nr_shared = 0
        nr_shared_inbulk = 0
        nr_shared_auto = 0
        nr_shared_inbulk_auto = 0
        
    }
    unique_vcf_fname = row[["unique_vcf"]]
    if (file.exists(unique_vcf_fname)){
        unique_vcf = readVcf(unique_vcf_fname)
        nr_unique = nrow(unique_vcf)
        vaf_bulk = get_vaf(unique_vcf, bulk)
        nr_unique_inbulk = sum(vaf_bulk > 0)
        
        seqs = granges(unique_vcf) %>% seqnames() %>% as.vector()
        unique_auto_vcf = unique_vcf[seqs %in% 1:22]
        nr_unique_auto = nrow(unique_auto_vcf)
        vaf_auto_bulk = get_vaf(unique_auto_vcf, bulk)
        nr_unique_inbulk_auto = sum(vaf_auto_bulk > 0)
        
    } else{
        nr_unique = 0
        nr_unique_inbulk = 0
        nr_unique_auto = 0
        nr_unique_inbulk_auto = 0
        
    }
    
    indel_shared_fname = row[["indel_shared"]]
    if (file.exists(indel_shared_fname)){
        shared_vcf = readVcf(indel_shared_fname)
        gt = geno(shared_vcf)$GT
        vaf_bulk = get_vaf(shared_vcf, bulk)
        nr_shared_indel = gt[,sample] %in% c("0/1", "1/1") %>% sum()
        nr_shared_inbulk_indel_f = gt[,sample] %in% c("0/1", "1/1") & vaf_bulk > 0
        nr_shared_inbulk_indel = sum(nr_shared_inbulk_indel_f)
        
        seqs = granges(shared_vcf) %>% seqnames() %>% as.vector()
        shared_auto_vcf = shared_vcf[seqs %in% 1:22]
        gt_auto = geno(shared_auto_vcf)$GT
        vaf_auto_bulk = get_vaf(shared_auto_vcf, bulk)
        nr_shared_indel_auto = gt_auto[,sample] %in% c("0/1", "1/1") %>% sum()
        nr_shared_inbulk_indel_f_auto = gt_auto[,sample] %in% c("0/1", "1/1") & vaf_auto_bulk > 0
        nr_shared_inbulk_indel_auto = sum(nr_shared_inbulk_indel_f_auto)
        
    } else{
        nr_shared_indel = 0
        nr_shared_inbulk_indel = 0
        nr_shared_indel_auto = 0
        nr_shared_inbulk_indel_auto = 0
        
    }
    indel_uniq_fname = row[["indel_uniq"]]
    if (file.exists(indel_uniq_fname)){
        unique_vcf = readVcf(indel_uniq_fname)
        nr_unique_indel = nrow(unique_vcf)
        vaf_bulk = get_vaf(unique_vcf, bulk)
        nr_unique_inbulk_indel = sum(vaf_bulk > 0)
        #nr_unique_inbulk = geno(unique_vcf)$GT[,bulk] %in% c("0/1", "1/1") %>% sum()
        
        seqs = granges(unique_vcf) %>% seqnames() %>% as.vector()
        unique_auto_vcf = unique_vcf[seqs %in% 1:22]
        nr_unique_indel_auto = nrow(unique_auto_vcf)
        vaf_auto_bulk = get_vaf(unique_auto_vcf, bulk)
        nr_unique_inbulk_indel_auto = sum(vaf_auto_bulk > 0)
        
    } else{
        nr_unique_indel = 0
        nr_unique_inbulk_indel = 0
        nr_unique_indel_auto = 0
        nr_unique_inbulk_indel_auto = 0
        
    }
    
    nr_total_snv = nr_shared + nr_unique
    nr_total_indel = nr_shared_indel + nr_unique_indel
    nr_total = nr_total_snv + nr_total_indel
    nr_total_snv_inbulk = nr_shared_inbulk + nr_unique_inbulk
    nr_total_snv_notbulk = nr_total_snv - nr_total_snv_inbulk
    nr_total_indel_inbulk = nr_shared_inbulk_indel + nr_unique_inbulk_indel
    nr_total_indel_notbulk = nr_total_indel - nr_total_indel_inbulk
    nr_total_inbulk = nr_total_snv_inbulk + nr_total_indel_inbulk
    nr_total_notbulk = nr_total_snv_notbulk + nr_total_indel_notbulk
    
    nr_total_snv_auto = nr_shared_auto + nr_unique_auto
    nr_total_indel_auto = nr_shared_indel_auto + nr_unique_indel_auto
    nr_total_auto = nr_total_snv_auto + nr_total_indel_auto
    nr_muts_tb = tibble(fetus_name, sample, nr_shared, nr_shared_inbulk, nr_unique, nr_unique_inbulk, nr_total_snv, nr_shared_indel, nr_shared_inbulk_indel, nr_unique_indel, nr_unique_inbulk_indel, nr_total_indel, nr_total, nr_shared_auto, nr_shared_inbulk_auto, nr_unique_auto, nr_unique_inbulk_auto, nr_total_snv_auto, nr_shared_indel_auto, nr_shared_inbulk_indel_auto, nr_unique_indel_auto, nr_unique_inbulk_indel_auto, nr_total_indel_auto, nr_total_auto, nr_total_snv_inbulk, nr_total_snv_notbulk, nr_total_indel_inbulk, nr_total_indel_notbulk, nr_total_inbulk, nr_total_notbulk)
    print(paste0("Counted mutations for sample: ", sample, " from fetus: ", fetus_name))
    return(nr_muts_tb)
}
nr_muts_tb = apply(overview_samples, 1, function(x) count_muts(x)) %>% bind_rows()
write_tsv(nr_muts_tb, paste0(out_dir_base, "count_muts.txt"))

####_______Make plots based on the mutation counts. age lines___________________####
nr_muts_tb = read_tsv(paste0(out_dir_base, "count_muts.txt"))
nr_muts_tb = left_join(nr_muts_tb, overview_samples[,c("fetus", "sample", "surveyed_region_pct", "origin", "celltype", "gender", "trisomy", "age_year", "age_weeks", "surveyed_region_auto_pct")], by = c("fetus_name" = "fetus", "sample")) %>% mutate(nr_total_corrected_snv = nr_total_snv / surveyed_region_pct)
#ages = tibble("fetus_name" = c("N01", "NR1", "NR2", "MH2"), "age_weeks" = c(10, 14, 12, 13)) #I substracted two weeks from the ages, because they were likely counted from start of pregnancy and not from conception
#weeksbyyear = 52.177457
nr_muts_tb %<>% mutate(muts_per_year = nr_total_snv / age_year,  muts_per_year_corrected = muts_per_year / surveyed_region_pct)
nr_muts_tb %<>% mutate(nr_total_snv_notbulk_auto = nr_total_snv_auto - (nr_shared_inbulk_auto + nr_unique_inbulk_auto), nr_total_snv_notbulk_auto_corrected = nr_total_snv_notbulk_auto / surveyed_region_auto_pct)
#nr_muts_tb$celltype = ifelse(grepl("SI|Si|-I-", nr_muts_tb$origin), "SI", "HSC")
write_tsv(nr_muts_tb, paste0(out_dir_base, "count_muts_extended.txt"))

#Create combined counts all, that is used later.
counts = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si|-I-", sample)) %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2")) %>% dplyr::select(Donor = fetus_name, norm_muts = nr_total_corrected_snv, Age = age_year)
counts_notbulk_auto = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si|-I-", sample)) %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2")) %>% dplyr::select(Donor = fetus_name, norm_muts = nr_total_snv_notbulk_auto_corrected, Age = age_year)
weeksbyyear = 52.177457
age_at_birth = 38 / weeksbyyear
#muts_counts_axel = read_tsv(mut_counts_axel_fname) %>% unique() %>% dplyr::filter(Tissue == "Blood") %>% dplyr::select(person = Donor, Age, norm_snvs = norm_muts) %>% mutate(Age = Age + age_at_birth)
muts_counts_axel = read_tsv(muts_counts_axel_fname) %>% dplyr::mutate(norm_muts = total_muts / percent_covered) %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth)
muts_counts_arianne = read_tsv(muts_counts_arianne_fname) %>% dplyr::filter(celltype == "HSC" | celltype == "MPP") %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth)
muts_counts_mirjam = read.table(muts_counts_mirjam_fname, sep = "\t", header = T) %>% as_tibble() %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth)
combined_counts_all = rbind(counts_notbulk_auto, muts_counts_axel, muts_counts_arianne, muts_counts_mirjam)
combined_counts_all %<>% mutate("age_cat" = ifelse(Age < 0.6, "fetal", "Post-infant"))
combined_counts_all$age_cat = factor(combined_counts_all$age_cat, levels = c("Post-infant", "fetal"))


###Linear mixed model age line for the fetuses
counts = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si|-I-", sample)) %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2"))
lme_res_fetuses = lme(nr_total_corrected_snv ~ age_year, random = ~ -1 + age_year | fetus_name, data = counts) #This models a random slope. As the mutation rate might be different between samples. The number of muts at conception is 0, so the intercept should always be 0 and is thus not a random variable.
lme_sum = summary(lme_res_fetuses)
vals_fetuses = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% mutate(test = "fetuses_snv_normalized")
pval_slope = lme_sum$tTable[[10]]

fetuses_fig = ggplot(counts, aes(x = age_year, y = nr_total_corrected_snv)) +
    geom_point(size = 4, shape = 21, fill = "red", color = "black") +
    geom_line(aes(y = predict(lme_res_fetuses, level = 0)), color = "red", size = 1.5) +
    labs(title = "Age line SNVs (fetuses)", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 0.235, y = 35, label = paste0("P = ", round(pval_slope, 3)), size = 6)


###Age line (snvs) for the fetuses, combined with adult data
#All blood data age lines
combined_counts = combined_counts_all
lme_all = lme(norm_muts ~ Age*age_cat, random = ~ -1 + Age |Donor, data = combined_counts) #Check if the slope is different between young and fetal. Not sure if this is correct formula
lme_sum = summary(lme_all)
vals_all = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% mutate(test = "fetuses_all_snv_normalized")

#Other calculations, result in the same pvalue as already present in the lme.
m.lst = lstrends(lme_all, "age_cat", var="Age")
comparisons = pairs(m.lst)
write_tsv(as_tibble(m.lst), paste0(out_dir_base, "values_youngvsold_agelines.txt"))
write_tsv(as_tibble(comparisons), paste0(out_dir_base, "values_youngvsold_contrast_agelines.txt"))

combined_counts$predict = predict(lme_all, level = 0)

pval_diff = lme_sum$tTable["Age:age_catfetal", "p-value"]
all_blood_log_fig = ggplot(combined_counts, aes(x = Age, y = norm_muts)) +
    geom_point(aes(fill = age_cat), size = 4, shape = 21, color = "black") +
    geom_line(aes(y = predict, color = age_cat), size = 1.5) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = "Age line SNVs (all blood data)", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1, y = 1100, label = paste0("P (diff between fetal and adult): ", round(pval_diff, 3)), size = 6)

all_blood_fig = ggplot(combined_counts, aes(x = Age, y = norm_muts)) +
    geom_point(aes(fill = age_cat), size = 4, shape = 21, color = "black") +
    geom_line(aes(y = predict, color = age_cat), size = 1.5) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    labs(title = "Age line SNVs (all blood data)", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 20, y = 1100, label = paste0("P (diff between fetal and adult): ", round(pval_diff, 3)), size = 6)

all_blood_zoom_fig = ggplot(combined_counts, aes(x = Age, y = norm_muts)) +
    geom_point(aes(fill = age_cat), size = 4, shape = 21, color = "black") +
    geom_line(aes(y = predict, color = age_cat), size = 1.5) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    labs(title = "Age line SNVs (all blood data)", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1, y = 1100, label = paste0("P (diff between fetal and adult): ", round(pval_diff, 3)), size = 6) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0,100))

#Still all data, but Age and age_cat are not crossed
combined_counts = combined_counts_all
lme_all2 = lme(norm_muts ~ Age+age_cat, random = ~ -1 + Age |Donor, data = combined_counts) #Check if the slope is different between young and fetal. Not sure if this is correct formula
lme_sum = summary(lme_all2)
vals_all2 = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% mutate(test = "fetuses_all_snv_normalized")
combined_counts$predict = predict(lme_all2, level = 0)

pval_diff = lme_sum$tTable["age_catfetal", "p-value"]

all_blood_zoom_fig2 = ggplot(combined_counts, aes(x = Age, y = norm_muts)) +
    geom_point(aes(fill = age_cat), size = 4, shape = 21, color = "black") +
    geom_line(aes(y = predict, color = age_cat), size = 1.5) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    labs(title = "Age line SNVs (all blood data)", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1, y = 100, label = paste0("P (diff between fetal and adult): ", round(pval_diff, 3)), size = 6) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0,100))




###Use only adult data for model
combined_counts_adult = combined_counts_all %>% dplyr::filter(age_cat != "fetal")
lme_adult = lme(norm_muts ~ Age, random = ~ -1 + Age |Donor, data = combined_counts_adult)
lme_sum_adult = summary(lme_adult)
vals_all_adult = lme_sum_adult$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% mutate(test = "Adult_all_snv_normalized")
pval = vals_all_adult[2, "p-value"]
combined_counts$predict_adult = predict(lme_adult, newdata = combined_counts, level = 0)
#combined_counts %<>% dplyr::mutate(resid = norm_muts - predict_adult)

all_blood_zoom_adult_fig = ggplot(combined_counts, aes(x = Age, y = norm_muts)) +
    geom_point(size = 4, shape = 21, color = "black", aes(fill = age_cat)) +
    geom_line(aes(y = predict_adult), size = 1.5) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    labs(title = "Model based on post-natal data", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1, y = 100, label = paste0("P = ", round(pval, 3)), size = 6) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0,100))


#Create a prediction interval. The interval is slightly extended. (Such an extension is not super clean)
lme_adult2 = lmer(norm_muts ~ Age + (-1 + Age |Donor), data = combined_counts_adult)
pi = ggpredict(lme_adult2, "Age", type = "re")
x_vals = tibble("x" = c(0, pi$x))
conf.low = lm(conf.low ~ x, data = pi) %>% predict(., newdata = x_vals)
predicted = lm(predicted ~ x, data = pi) %>% predict(., newdata = x_vals)
conf.high = lm(conf.high ~ x, data = pi) %>% predict(., newdata = x_vals)
pi = cbind("x" = x_vals$x, conf.low, predicted, conf.high) %>% as.data.frame()

all_blood_adult_pi_fig = ggplot(data = pi, aes(x = x, y = predicted)) +
    geom_line(color = "red") +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
    geom_point(data = combined_counts, size = 4, shape = 21, color = "black", aes(fill = age_cat, y = norm_muts, x = Age)) +
    labs(title = "Model based on post-natal data", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"))

#lmvar prediction interval
X = model.matrix(~ Age - 1, combined_counts_adult)
fit = lmvar(combined_counts_adult$norm_muts, X_mu = X, X_sigma = X)
intervals = fitted(fit, interval = "prediction", level = 0.95)
combined_counts_adult_lmvar = cbind(combined_counts_adult, intervals[,c("mu", "lwr", "upr")])

x_vals = tibble("Age" = c(0, combined_counts$Age))
lwr = lm(lwr ~ poly(Age, 4), data = combined_counts_adult_lmvar) %>% predict(., newdata = x_vals)
mu = lm(mu ~ poly(Age, 4), data = combined_counts_adult_lmvar) %>% predict(., newdata = x_vals)
upr = lm(upr ~ poly(Age, 4), data = combined_counts_adult_lmvar) %>% predict(., newdata = x_vals)
extended_pi = cbind(lwr, mu, upr, x_vals) %>% as.data.frame()

all_blood_adult_pi_fig2 = ggplot(data = combined_counts_adult_lmvar, aes(x = Age, y = mu)) +
    geom_line(color = "black") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), data = extended_pi, alpha = 0.05, fill = "red") +
    geom_line(aes(x = Age, y = mu), data = extended_pi, color  = "red", linetype = "dashed") +
    #geom_point(aes(y = norm_muts)) +
    geom_point(data = combined_counts, size = 4, shape = 21, color = "black", aes(fill = age_cat, y = norm_muts, x = Age)) +
    labs(title = "Model based on post-natal data", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"))

all_blood_adult_zoom_pi_fig2 = ggplot(data = combined_counts_adult_lmvar, aes(x = Age, y = mu)) +
    geom_line(color = "black") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = .1) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), data = extended_pi, alpha = 0.05, fill = "red") +
    geom_line(aes(x = Age, y = mu), data = extended_pi, color  = "red", linetype = "dashed") +
    #geom_point(aes(y = norm_muts)) +
    geom_point(data = combined_counts, size = 4, shape = 21, color = "black", aes(fill = age_cat, y = norm_muts, x = Age)) +
    labs(title = "Model based on post-natal data", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0,100))


###use all data with no seperate categories.
combined_counts = combined_counts_all
lme_all_nocat = lme(norm_muts ~ Age, random = ~ -1 + Age |Donor, data = combined_counts)
lme_sum_nocat = summary(lme_all_nocat)
vals_all_nocat = lme_sum_nocat$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% mutate(test = "fetuses_all_snv_normalized")
pval = vals_all_nocat[2, "p-value"]
combined_counts$predict_nocat = predict(lme_all_nocat, level = 0)


all_blood_zoom_nocat_fig = ggplot(combined_counts, aes(x = Age, y = norm_muts)) +
    geom_point(size = 4, shape = 21, color = "black", fill = "red") +
    geom_line(aes(y = predict_nocat), size = 1.5) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    labs(title = "Model based on all data", x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1, y = 100, label = paste0("P = ", round(pval, 3)), size = 6) +
    coord_cartesian(xlim = c(0, 2), ylim = c(0,100))

#Compare models
combined_counts = combined_counts_all
lme_all = lme(norm_muts ~ Age*age_cat, random = ~ -1 + Age |Donor, data = combined_counts, method = "ML")
lme_all2 = lme(norm_muts ~ Age+age_cat, random = ~ -1 + Age |Donor, data = combined_counts, method = "ML")
lme_all3 = lme(norm_muts ~ Age + Age:age_cat, random = ~ -1 + Age |Donor, data = combined_counts, method = "ML")
lme_all_nocat = lme(norm_muts ~ Age, random = ~ -1 + Age |Donor, data = combined_counts, method = "ML")


an1 = anova(lme_all, lme_all_nocat) %>% as.data.frame() %>% rownames_to_column("model_name")
an2 = anova(lme_all2, lme_all_nocat) %>% as.data.frame() %>% rownames_to_column("model_name")
an3 = anova(lme_all3, lme_all_nocat) %>% as.data.frame() %>% rownames_to_column("model_name")
ans = rbind(an1, an2, an3)
write_tsv(ans, paste0(out_dir_base, "age_lines_anovas.txt"))


outlier_fig = plot(lme_all_nocat, id = 0.05, idLabels = combined_counts$Donor)

lm_all = lm(norm_muts ~ Age, data = combined_counts) #Check if model performs much worse without the lme random effects
anova(lme_all_nocat, lm_all)

#explain age line issue
pdf(paste0(out_dir_base, "age_line_issue.pdf"))
print(all_blood_fig)
print(all_blood_zoom_adult_fig)
print(all_blood_zoom_fig)
print(all_blood_zoom_fig2)
print(all_blood_zoom_nocat_fig)
print(outlier_fig)
print(all_blood_adult_pi_fig)
print(all_blood_adult_pi_fig2)
print(all_blood_adult_zoom_pi_fig2)
dev.off()

###Look how much points deviate from the age line
combined_counts_adult = combined_counts_all %>% dplyr::filter(age_cat != "fetal")
lme_adult = lme(norm_muts ~ Age, random = ~ -1 + Age |Donor, data = combined_counts_adult)
lme_sum_adult = summary(lme_adult)
vals_all_adult = lme_sum_adult$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% mutate(test = "adult_snv_normalized")
pval = vals_all_adult[2, "p-value"]
combined_counts$predict_adult = predict(lme_adult, newdata = combined_counts, level = 0)

var_dev_tb = combined_counts %>% dplyr::group_by(Donor) %>% dplyr::summarise(mean_muts = mean(norm_muts), Age = Age[1], n = n(), predict_adult = predict_adult[1], sd = sd(norm_muts)) %>% mutate(dev_from_pred = mean_muts - predict_adult, dev_by_sd = abs(dev_from_pred / sd), norm_dev = abs(dev_from_pred / predict_adult))
var_dev_tb %<>% mutate(age_cat = factor(ifelse(Age < 0.6, "Fetal", ifelse(Age < 1, "Cord blood", "Post-infant")), levels = c("Fetal", "Cord blood", "Post-infant")))

my_comparisons = combn(levels(var_dev_tb$age_cat), 2, simplify = F)
dev_fig = ggplot(data = var_dev_tb, aes(x = age_cat, y = abs(dev_from_pred), color = age_cat)) +
    geom_point(size = 3) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    labs(title = "Deviation from expected based on age line (Post-infant)", x = "", y = "Abs. Mean deviation from expected") +
    stat_compare_means(label.y = 150, size = 6) + 
    stat_compare_means(comparisons = my_comparisons, size = 6) +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"))


rel_dev_fig = ggplot(data = var_dev_tb, aes(x = age_cat, y = norm_dev, color = age_cat)) +
    geom_point(size = 3) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    labs(title = "Deviation from expected based on age line (Post-infant)", x = "", y = "Abs. Mean percent deviation from expected") +
    stat_compare_means(label.y = 0.5, size = 6) + 
    stat_compare_means(comparisons = my_comparisons, size = 6) +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"))


var_dev_tb2 = var_dev_tb %>% dplyr::filter(n >= 3)
z_dev_fig = ggplot(data = var_dev_tb2, aes(x = age_cat, y = dev_by_sd, color = age_cat)) +
    geom_point(size = 3) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    labs(title = "Deviation from expected based on age line (Post-infant)", x = "", y = "Abs. Mean deviation / sd") +
    stat_compare_means(label.y = 3, size = 6) + 
    stat_compare_means(comparisons = my_comparisons, size = 6) +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"))


pdf(paste0(out_dir_base, "age_line_deviations.pdf"), width = 8)
dev_fig
rel_dev_fig
z_dev_fig
dev.off()


###Look at within-group standard error vs random effect size for adult data and both fetal and adult data. NOTE!! The random effect of the fetuses is here modeled as a random intercept
rndm_effect_sizes = lme(norm_muts ~ Age, random = ~ 1 |Donor, data = combined_counts_adult) %>% intervals()
rndm_effect_sizes_df = rbind("Random effect size" = rndm_effect_sizes$reStruct$Donor, "Within-group standard error" = rndm_effect_sizes$sigma) %>% rownames_to_column("type")
rndm_effect_sizes_fig = ggplot(data = rndm_effect_sizes_df, aes(x = type)) +
    geom_point(aes(y = est.), size = 4) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, size = 1) +
    labs(title = "non-fetal", x = "", y = "95% CI") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"))

rndm_effect_sizes = lme(norm_muts ~ Age, random = ~ 1 |Donor, data = combined_counts) %>% intervals()
rndm_effect_sizes_df = rbind("Random effect size" = rndm_effect_sizes$reStruct$Donor, "Within-group standard error" = rndm_effect_sizes$sigma) %>% rownames_to_column("type")
rndm_effect_sizes_fig2 = ggplot(data = rndm_effect_sizes_df, aes(x = type)) +
    geom_point(aes(y = est.), size = 4) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, size = 1) +
    labs(title = "all healthy data", x = "", y = "95% CI") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"))


pdf(paste0(out_dir_base, "within_error_vs_random_effect_size.pdf"))
rndm_effect_sizes_fig
rndm_effect_sizes_fig2
dev.off()


###Look at leukemia clones. The adult lme model is used to draw the age line
counts_fetal = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si|-I-", sample)) %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2")) %>% dplyr::select(Donor = fetus_name, norm_muts = nr_total_corrected_snv, Age = age_year) %>% dplyr::mutate(type = "Healthy")
muts_counts_axel = read_tsv(muts_counts_axel_fname) %>% dplyr::mutate(norm_muts = total_muts / percent_covered) %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth, type = "Healthy")
muts_counts_arianne = read_tsv(muts_counts_arianne_fname) %>% mutate(type = ifelse(celltype == "HSC" | celltype == "MPP", "Healthy", "Leukemia")) %>% dplyr::select(Donor, Age, norm_muts, type) %>% dplyr::mutate(Age = Age + age_at_birth)
muts_counts_mirjam = read.table(muts_counts_mirjam_fname, sep = "\t", header = T) %>% as_tibble() %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth, type = "Healthy")
combined_counts_sick = rbind(counts_fetal, muts_counts_axel, muts_counts_arianne, muts_counts_mirjam)
combined_counts_sick$predict_adult = predict(lme_adult, newdata = combined_counts_sick, level = 0)

zoom_abline_fig = ggplot(combined_counts_sick, aes(x = Age, y = norm_muts)) +
    geom_point(aes(fill = type), size = 4, shape = 21, color = "black") +
    geom_line(aes(y = predict_adult), size = 1.5, color = "red") +
    geom_abline(slope = 115) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    labs(x = "Age (from conception)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1, y = 1100, label = paste0("P (diff between fetal and adult): ", round(pval_diff, 3)), size = 6) +
    coord_cartesian(xlim = c(0, 5), ylim = c(0,300), expand = F)


###Only cord blood and fetuses
counts_forcord = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si|-I-", sample)) %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2")) %>% dplyr::select(Donor = fetus_name, norm_muts = nr_total_corrected_snv, Age = age_year)
muts_counts_axel = read_tsv(muts_counts_axel_fname) %>% dplyr::mutate(norm_muts = total_muts / percent_covered) %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth)
muts_counts_axel_cord = muts_counts_axel %>% dplyr::filter(Age == age_at_birth)
combined_counts_cord = rbind(counts_forcord, muts_counts_axel_cord)
lme_res_cord = lme(norm_muts ~ Age, random = ~ -1 + Age | Donor, data = combined_counts_cord) #This models a random slope. As the mutation rate might be different between samples. The number of muts at conception is 0, so the intercept should always be 0 and is thus not a random variable.
lme_sum = summary(lme_res_cord)
pval_slope = lme_sum$tTable["Age","p-value"]
vals_cord = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% mutate(test = "fetuses_cord_snv_normalized")

fetuses_cord_fig = ggplot(combined_counts_cord, aes(x = Age, y = norm_muts)) +
    geom_point(size = 5, shape = 21, fill = "red", color = "black") +
    geom_line(aes(y = predict(lme_res_cord, level = 0)), color = "red", size = 1.5) +
    labs(title = "Age line SNVs (fetuses + cord)", x = "Age (years)", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 0.3, y = 50, label = paste0("P = ", round(pval_slope, 3)), size = 6)
    

#Write out the age lines
pdf(paste0(out_dir_base, "age_lines.pdf"))
fetuses_fig
all_blood_fig
all_blood_log_fig
all_blood_zoom_fig
all_blood_zoom_fig2
fetuses_cord_fig
zoom_abline_fig
all_blood_zoom_adult_fig
all_blood_zoom_nocat_fig
all_blood_adult_pi_fig
all_blood_adult_pi_fig2
all_blood_adult_zoom_pi_fig2
dev.off()

#Write out the pvals, std errors and values of all the slopes and intercepts
vals_all = rbind(vals_fetuses, vals_all, vals_all2, vals_cord, vals_all_adult)
write_tsv(vals_all, paste0(out_dir_base, "values_agelines.txt"))


####___________________Mutation count fetus vs adult comparisons_____________________####
combined_mut_counts = combined_counts_all %>% mutate(norm_muts_byage = norm_muts / Age)
combined_mut_counts$age_cat = factor(combined_mut_counts$age_cat, levels = c("fetal", "Post-infant"))
combined_mut_counts$age_cat = revalue(combined_mut_counts$age_cat, c("fetal" = "Fetal"))
lme_divstri = lme(norm_muts_byage ~ age_cat, random = ~ 1 | Donor, data = combined_mut_counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
lme_sum = summary(lme_divstri)
vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
pval = vals[2, "p-value"]
mutcount_fetalvsadult_lmm_fig = ggplot(combined_mut_counts, aes(x = age_cat, y = norm_muts_byage)) +
    geom_boxplot(aes(color = age_cat), outlier.shape = NA, fill = NA) +
    geom_quasirandom(aes(fill = Donor), shape = 21, size = 3) +
    #geom_dotplot(aes(fill = Donor), binaxis = "y", stackdir = "center", binwidth = 2, stackgroups = T, binpositions = "all", dotsize = 2) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(title = "Nr of substitutions in HSPC clones per year corrected for surveyed region", x = "", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1.5, y = 220, label = paste0("P = ", round(pval, 3), " (lmm)"), size = 6) +
    coord_cartesian(ylim = c(0,240))


#Same as before. However now only autosomal non bulk mutations are used. As other mutations would not be identified in Axels data. The comparison is also only made to axels data, as the rest is not yet published and might not be healthy.
counts = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si", sample)) %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2")) %>% dplyr::select(Donor = fetus_name, norm_muts = nr_total_snv_notbulk_auto_corrected, Age = age_year)
age_at_birth = 38 / weeksbyyear
#muts_counts_axel = read_tsv(mut_counts_axel_fname) %>% unique() %>% dplyr::filter(Tissue == "Blood") %>% dplyr::select(person = Donor, Age, norm_snvs = norm_muts) %>% mutate(Age = Age + age_at_birth)
muts_counts_axel = read_tsv(muts_counts_axel_fname) %>% dplyr::mutate(norm_muts = total_muts / percent_covered) %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth)
#muts_counts_mirjam = read.table(muts_counts_mirjam_fname, sep = "\t", header = T) %>% as_tibble() %>% dplyr::select(Donor, Age, norm_muts) %>% dplyr::mutate(Age = Age + age_at_birth)
combined_counts = rbind(counts, muts_counts_axel)
combined_counts %<>% mutate("age_cat" = ifelse(Age < 0.6, "Fetal", "Post-infant"))
combined_counts$age_cat = factor(combined_counts$age_cat, levels = c("Fetal", "Post-infant"))
combined_mut_counts = combined_counts %>% mutate(norm_muts_byage = norm_muts / Age)
lme_divstri = lme(norm_muts_byage ~ age_cat, random = ~ 1 | Donor, data = combined_mut_counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
lme_sum = summary(lme_divstri)
vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
pval = vals[2, "p-value"]
mutcount_fetalvsadult_lmm_2_fig = ggplot(combined_mut_counts, aes(x = age_cat, y = norm_muts_byage)) +
    geom_boxplot(aes(color = age_cat), outlier.shape = NA, fill = NA) +
    geom_quasirandom(aes(fill = Donor), shape = 21, size = 3) +
    #geom_dotplot(aes(fill = Donor), binaxis = "y", stackdir = "center", binwidth = 2, stackgroups = T, binpositions = "all", dotsize = 2) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(title = "Nr of autosomal non-bulk substitutions in HSPC clones per year corrected for surveyed region", x = "", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1.5, y = 220, label = paste0("P = ", round(pval, 3), " (lmm)"), size = 6) +
    coord_cartesian(ylim = c(0,240))

write_tsv(vals, paste0(out_dir_base, "nrmuts_fetalvsnonfetal_vals.txt"))


#Use cord blood as a seperate age category.
levels(combined_mut_counts$age_cat) = c(levels(combined_mut_counts$age_cat), "Cord blood")
combined_mut_counts$age_cat = factor(combined_mut_counts$age_cat, levels = c("Fetal", "Cord blood", "Post-infant"))

combined_mut_counts$age_cat[combined_mut_counts$Age > 0.6 & combined_mut_counts$Age < 1] = "Cord blood"

comp1 = combined_mut_counts %>% dplyr::filter(age_cat != "Post-infant")
lme_comp1 = lme(norm_muts_byage ~ age_cat, random = ~ 1 | Donor, data = comp1) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
#lme_sum1 = summary(lme_comp1)
#vals1 = lme_sum1$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
vals1 = get_fixed(lme_comp1)
comp2 = combined_mut_counts %>% dplyr::filter(age_cat != "Cord blood")
lme_comp2 = lme(norm_muts_byage ~ age_cat, random = ~ 1 | Donor, data = comp2) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
#lme_sum2 = summary(lme_comp2)
#vals2 = lme_sum2$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
vals2 = get_fixed(lme_comp2)
comp3 = combined_mut_counts %>% dplyr::filter(age_cat != "Fetal")
lme_comp3 = lme(norm_muts_byage ~ age_cat, random = ~ 1 | Donor, data = comp3) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
#lme_sum3 = summary(lme_comp3)
#vals3 = lme_sum3$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
vals3 = get_fixed(lme_comp3)
vals = rbind(vals1, vals2, vals3)
vals$category = c(rep("CordbloodvsFetal", 2), rep("AdultvsFetal", 2), rep("AdultvsCordblood", 2))
write_tsv(vals, paste0(out_dir_base, "nrmuts_fetalvsnonfetalvscordblood_vals.txt"))

CN = combn(levels(combined_mut_counts$age_cat), 2, simplify = FALSE)
pvals = vals[c(2,4,6), "p-value"] %>% round(3) #Order in the same way as the CN
hspc_mean_rates = combined_mut_counts %>% dplyr::group_by(age_cat) %>% dplyr::summarise("mean_hspc" = mean(norm_muts_byage)) %>% dplyr::mutate(age_cat = as.character(age_cat))
mutcount_fetalvsadult_lmm_3_fig = ggplot(combined_mut_counts, aes(x = age_cat, y = norm_muts_byage)) +
    geom_boxplot(aes(color = age_cat), outlier.shape = NA, fill = NA) +
    geom_quasirandom(aes(fill = Donor), shape = 21, size = 3) +
    #geom_dotplot(aes(fill = Donor), binaxis = "y", stackdir = "center", binwidth = 2, stackgroups = T, binpositions = "all", dotsize = 2) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(title = "Nr of autosomal non-bulk substitutions in HSPC clones per year corrected for surveyed region", x = "", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    coord_cartesian(ylim = c(0,240)) +
    geom_signif(comparisons = CN, annotations = pvals, y_position = c(180, 200, 220), size = 1, textsize = 6)


#Look at SI
counts = nr_muts_tb %>% dplyr::filter(celltype == "SI") %>% dplyr::filter(fetus_name %in% c("NR1", "NR2", "MH2", "F100916W15", "F100916W17")) %>% dplyr::select(Donor = fetus_name, norm_muts = nr_total_snv_notbulk_auto_corrected, Age = age_year)
age_at_birth = 38 / weeksbyyear
chr = paste0("chr", c(1:22))
length_chr_auto = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chr] %>% sum()
muts_counts_si = read_tsv(muts_counts_si_fname) %>% dplyr::mutate(percent_covered = Percent_covered_old * 2682655440 / length_chr_auto / 100, norm_muts = Total_muts / percent_covered) %>% dplyr::select(Donor, norm_muts, Age)
combined_counts = rbind(counts, muts_counts_si)
combined_counts %<>% mutate("age_cat" = ifelse(Age < 0.6, "Fetal", "Post-infant"))
combined_counts$age_cat = factor(combined_counts$age_cat, levels = c("Fetal", "Post-infant"))
combined_mut_counts = combined_counts %>% mutate(norm_muts_byage = norm_muts / Age)
lme_divstri = lme(norm_muts_byage ~ age_cat, random = ~ 1 | Donor, data = combined_mut_counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
lme_sum = summary(lme_divstri)
vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
pval = vals[2, "p-value"]
si_mean_rates = combined_mut_counts %>% dplyr::group_by(age_cat) %>% dplyr::summarise("mean_si" = mean(norm_muts_byage)) %>% dplyr::mutate(age_cat = as.character(age_cat))
mutcount_fetalvsadult_lmm_4_fig = ggplot(combined_mut_counts, aes(x = age_cat, y = norm_muts_byage)) +
    geom_boxplot(aes(color = age_cat), outlier.shape = NA, fill = NA) +
    geom_quasirandom(aes(fill = Donor), shape = 21, size = 3) +
    #geom_dotplot(aes(fill = Donor), binaxis = "y", stackdir = "center", binwidth = 2, stackgroups = T, binpositions = "all", dotsize = 2) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(title = "Nr of autosomal non-bulk substitutions in SI clones per year corrected for surveyed region", x = "", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 1.5, y = 220, label = paste0("P = ", round(pval, 3), " (lmm)"), size = 6) +
    coord_cartesian(ylim = c(0,240))

write_tsv(vals, paste0(out_dir_base, "nrmuts_fetalvsnonfetal_si_vals.txt"))

mean_rates = full_join(hspc_mean_rates, si_mean_rates, by = "age_cat")
write_tsv(mean_rates, paste0(out_dir_base, "mean_mut_rate_agecats.txt"))

pdf(paste0(out_dir_base, "muts_per_fetus_fetalvsadult.pdf"), width = 10)
mutcount_fetalvsadult_lmm_fig
mutcount_fetalvsadult_lmm_2_fig
mutcount_fetalvsadult_lmm_3_fig
mutcount_fetalvsadult_lmm_4_fig
dev.off()

#ggplot(combined_mut_counts, aes(x = age_cat, y = norm_muts_byage)) +
#     geom_boxplot(aes(color = age_cat), outlier.shape = NA, fill = NA) +
#     geom_jitter(width = 0.1, height = 0, aes(color = Donor), size = 3) +
#     labs(title = "Nr of mutations in HSC clones per year corrected for surveyed region", x = "", y = "Nr. of mutations per year") +
#     theme_classic() +
#     theme(text = element_text(size = 20)) +
#     annotate("text", x = 1.5, y = 160, label = paste0("P = ", round(pval, 3), " (lmm)"), size = 6)




####_______________________Mutation count comparisons di vs tri________________####
make_mut_count_figs = function(counts, exp_name, lims, muttype, out_dir_base){
    muts_fig = ggplot(counts, aes(x = fetus_name, y = nr_muts, fill = factor(origin))) +
        geom_boxplot(aes(color = fetus_name), outlier.shape = NA, fill = NA) +
        geom_jitter(width = 0.1, height = 1, size = 2, shape = 21) +
        #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2, stackgroups = T, binpositions = "all", dotsize = 1) +
        scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = paste0("Nr of substitutions in ", exp_name, " clones per fetus"), x = "Fetus", y = "Nr. of substitutions") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        coord_cartesian(ylim = c(0, lims[1]))
    
    muts_corrected_fig = ggplot(counts, aes(x = fetus_name, y = nr_total_corrected, fill = origin)) +
        geom_boxplot(aes(color = fetus_name), outlier.shape = NA, fill = NA) +
        geom_jitter(width = 0.1, height = 1, size = 2, shape = 21) +
        #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 2, stackgroups = T, binpositions = "all", dotsize = 1) +
        scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = paste0("Nr of substitutions in ", exp_name, " clones per fetus corrected for surveyed region"), x = "Fetus", y = "Nr. of substitutions") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        coord_cartesian(ylim = c(0, lims[1]))
    
    
    muts_year_fig = ggplot(counts, aes(x = fetus_name, y = muts_per_year, fill = origin)) +
        geom_boxplot(aes(color = fetus_name), outlier.shape = NA, fill = NA) +
        geom_jitter(width = 0.1, height = 1, size = 2, shape = 21) +
        #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 3, stackgroups = T, binpositions = "all", dotsize = 3) +
        scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = paste0("Nr of substitutions in ", exp_name, " clones per fetus per year"), x = "Fetus", y = "Number of somatic base substitutions\n per genome per year") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        coord_cartesian(ylim = c(0, lims[2]))
    
    
    muts_year_corrected_fig = ggplot(counts, aes(x = fetus_name, y = muts_per_year_corrected, fill = origin)) +
        geom_boxplot(aes(color = fetus_name), outlier.shape = NA, fill = NA) +
        geom_jitter(width = 0.1, height = 1, size = 2, shape = 21) +
        #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 3, stackgroups = T, binpositions = "all", dotsize = 3) +
        scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = paste0("Nr of substitutions in ", exp_name, " clones per fetus per year corrected for surveyed region"), x = "Fetus", y = "Number of somatic base substitutions\n per genome per year") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        coord_cartesian(ylim = c(0, lims[2]))
    
    
    #Compare disomy to trisomy (wilcoxon)
    disomy_muts_year_corrected_fig = ggplot(counts, aes(x = trisomy, y = muts_per_year_corrected)) +
        geom_boxplot(aes(color = trisomy), outlier.shape = NA, fill = NA) +
        geom_jitter(aes(fill = fetus_name), width = 0.1, height = 1, size = 2, shape = 21) +
        #geom_dotplot(aes(fill = fetus_name), binaxis = "y", stackdir = "center", binwidth = 3, stackgroups = T, binpositions = "all", dotsize = 3) +
        scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = paste0("Nr of substitutions in ", exp_name, " clones per year corrected for surveyed region"), x = "Trisomy", y = "Number of somatic base substitutions\n per genome per year") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        stat_compare_means(cex = 8, label.x.npc = 0.4) +
        coord_cartesian(ylim = c(0, lims[2]))
    
        
    #Compare disomy to trisomy (lmm)
    lme_divstri = lme(muts_per_year_corrected ~ trisomy, random = ~ 1 | fetus_name, data = counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
    lme_sum = summary(lme_divstri)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    pval = vals[2, "p-value"]
    disomy_muts_year_corrected_lmm_fig = ggplot(counts, aes(x = trisomy, y = muts_per_year_corrected)) +
        geom_boxplot(aes(color = trisomy), outlier.shape = NA, fill = NA) +
        geom_jitter(aes(fill = fetus_name), width = 0.1, height = 1, size = 2, shape = 21) +
        #geom_dotplot(aes(fill = fetus_name), binaxis = "y", stackdir = "center", binwidth = 3, stackgroups = T, binpositions = "all", dotsize = 3) +
        scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = paste0("Nr of substitutions in ", exp_name, " clones per year corrected for surveyed region"), x = "Trisomy", y = "Number of somatic base substitutions\n per genome per year") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        annotate("text", x = 1.5, y = (lims[2] - 0.2*lims[2]), label = paste0("P = ", round(pval, 3), " (lmm)"), size = 6) +
        coord_cartesian(ylim = c(0, lims[2]))
    
    
    write_tsv(vals, paste0(out_dir_base, exp_name, "_", muttype, "_nrmuts_divstri_vals.txt"))
    figs = list(muts_fig, muts_corrected_fig, muts_year_fig, muts_year_corrected_fig, disomy_muts_year_corrected_fig, disomy_muts_year_corrected_lmm_fig)
    return(figs)
}


plot_divstri_plots = function(nr_muts_tb_raw, mut, muttype, lims, out_dir_base){
    
    #Calculate corrected mutations counts and mutation counts per year
    nr_muts_tb = nr_muts_tb_raw
    nr_muts_tb$nr_muts = nr_muts_tb %>% pull(!!mut)
    nr_muts_tb %<>% mutate(nr_total_corrected = nr_muts / surveyed_region_pct)
    nr_muts_tb %<>% mutate(muts_per_year = nr_muts / age_year,  muts_per_year_corrected = muts_per_year / surveyed_region_pct)
    
    #Create figures for different celltypes
    counts = nr_muts_tb %>% dplyr::filter(!grepl("SI|Si|-I-", sample)) %>% mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), to = c("T21", "D21")))
    hsc_mut_count_figs = make_mut_count_figs(counts, "HSPC", lims, muttype, out_dir_base)
    counts = nr_muts_tb %>% mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), to = c("T21", "D21")))
    all_mut_count_figs = make_mut_count_figs(counts, "ALL", lims, muttype, out_dir_base)
    counts = nr_muts_tb %>% dplyr::filter(grepl("SI|Si|-I-", sample)) %>% mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), to = c("T21", "D21")))
    si_mut_count_figs = make_mut_count_figs(counts, "SI", lims, muttype, out_dir_base)
    
    
    
    #Compare SI vs HSC
    counts = nr_muts_tb %>% dplyr::filter(trisomy == "21")
    lme_hscvssi_tri = lme(muts_per_year_corrected ~ celltype, random = ~ 1 | fetus_name, data = counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
    lme_sum = summary(lme_hscvssi_tri)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    pval_tri = vals[2, "p-value"]
    write_tsv(vals, paste0(out_dir_base, muttype, "_nrmuts_hsvvssi_tri_vals.txt"))
    
    counts = nr_muts_tb %>% dplyr::filter(trisomy == "FALSE")
    lme_hscvssi_di = lme(muts_per_year_corrected ~ celltype, random = ~ 1 | fetus_name, data = counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
    lme_sum = summary(lme_hscvssi_di)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    pval_di = vals[2, "p-value"]
    write_tsv(vals, paste0(out_dir_base, muttype, "_nrmuts_hsvvssi_di_vals.txt"))
    
    counts = nr_muts_tb
    counts$celltype = revalue(counts$celltype, c("HSC" = "HSPC"))
    pvals = c(pval_tri, pval_di)
    dat_text = tibble(trisomy = c("21", "FALSE"), label = paste0("P = ", round(pvals, 3), " (lmm)"))
    hscvssi_fig = ggplot(counts, aes(x = celltype, y = muts_per_year_corrected)) +
        geom_boxplot(outlier.shape = NA, fill = NA) + 
        geom_jitter(aes(color = fetus_name), width = 0.1, height = 1, size = 2) +
        #geom_dotplot(aes(fill = fetus_name), binaxis = "y", stackdir = "center", binwidth = 3, dotsize = 3) +
        facet_grid(. ~ trisomy) +
        #scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = "Nr of substitutions per year between HSPC and SI", x = "Cell type", y = "Nr. of substitutions") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        geom_text(data = dat_text, aes(label = label, x = 1.5, y = (lims[2] - 0.2*lims[2])), size = 6) +
        coord_cartesian(ylim = c(0,lims[2]))
    
    
    #Compare D21 with T21 for both SI and HSC in a single figure
    counts = nr_muts_tb %>% dplyr::filter(celltype == "HSC")
    lme_res = lme(muts_per_year_corrected ~ trisomy, random = ~ 1 | fetus_name, data = counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
    lme_sum = summary(lme_res)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    pval_hsc = vals[2, "p-value"]
    
    counts = nr_muts_tb %>% dplyr::filter(celltype == "SI")
    lme_res = lme(muts_per_year_corrected ~ trisomy, random = ~ 1 | fetus_name, data = counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.
    lme_sum = summary(lme_res)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    pval_si = vals[2, "p-value"]
    
    counts = nr_muts_tb
    counts$celltype = revalue(counts$celltype, c("HSC" = "HSPC"))
    counts$trisomy = revalue(counts$trisomy, c("21" = "T21", "FALSE" = "D21"))
    pvals = c(pval_hsc, pval_si)
    dat_text = tibble(celltype = c("HSPC", "SI"), label = paste0("P = ", round(pvals, 3), " (lmm)"))
    trivsdi_fig = ggplot(counts, aes(x = trisomy, y = muts_per_year_corrected)) +
        geom_boxplot(outlier.shape = NA, fill = NA) + 
        geom_jitter(aes(color = fetus_name), width = 0.1, height = 1, size = 2) +
        #geom_dotplot(aes(fill = fetus_name), binaxis = "y", stackdir = "center", binwidth = 40, dotsize = 0.5) +
        facet_grid(. ~ celltype) +
        #scale_color_brewer(palette = "Dark2", guide = F) +
        labs(title = "Nr of substitutions per year between HSPC and SI", x = "Type", y = "Number of somatic base substitutions\n per genome per year") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        geom_text(data = dat_text, aes(label = label, x = 1.5, y = (lims[2] - 0.2*lims[2])), size = 6) +
        coord_cartesian(ylim = c(0,lims[2]))
    
    ###Look at combination of trisomy state with celltype.
    counts = nr_muts_tb
    lme_divstri_combi = lme(muts_per_year_corrected ~ trisomy*celltype, random = ~ 1 | fetus_name, data = counts) #Random intercept per fetus. Some fetuses might have a higher base mutation rate by chance.

    #Plot outliers 
    outlier_fig = plot(lme_divstri_combi, id = 0.05, idLabels = counts$sample)
    
    #Allow the variance to be different between disomy and trisomy.
    lme_divstri_combi_diffvar = lme(muts_per_year_corrected ~ trisomy*celltype, random = ~ 1 | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy))
    anova_res = anova(lme_divstri_combi, lme_divstri_combi_diffvar)
    write_tsv(anova_res, paste0(out_dir_base, muttype, "_diffvar_anova.txt"))
    outlier_fig2 = plot(lme_divstri_combi_diffvar, id = 0.05, idLabels = counts$sample)
    
    m1_fixed = get_fixed(lme_divstri_combi)
    m2_fixed = get_fixed(lme_divstri_combi_diffvar)
    m_fixed = rbind(m1_fixed, m2_fixed)
    write_tsv(m_fixed, paste0(out_dir_base, muttype, "_divstri_combi_models_.txt"))
    
    
    #Use age as an explanatory variable. Trisomy and celltype are crossed explanatory variables. Also allow the variance to be different between trisomy and disomy.
    counts$trisomy = factor(mapvalues(counts$trisomy, from = c("21", "FALSE"), to = c("T21", "D21")), levels = c("D21", "T21"))
    counts$celltype = factor(counts$celltype, levels = c("HSC", "SI"))
    #counts$celltype = factor(counts$celltype, levels = c("SI", "HSC"))
    lme_divstri_age_vardiff = lme(nr_total_corrected ~ age_year + trisomy * celltype, random = ~ -1 + age_year | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy))
    resid(lme_divstri_age_vardiff, type = "pearson") %>% abs() %>% multiply_by(-1) %>% pnorm() %>% multiply_by(2) %>% p.adjust(method = "fdr")
    lme_divstri_age = lme(nr_total_corrected ~ age_year + trisomy * celltype, random = ~ -1 + age_year | fetus_name, data = counts)
    resid(lme_divstri_age, type = "pearson") %>% abs() %>% multiply_by(-1) %>% pnorm() %>% multiply_by(2) %>% p.adjust(method = "fdr")
    anova_res = anova(lme_divstri_age_vardiff, lme_divstri_age)
    write_tsv(anova_res, paste0(out_dir_base, muttype, "_diffvar_anova_withage.txt"))
    
    #lme_divstri_age_vardiff = lme(nr_total_corrected ~ age_year + trisomy * celltype, random = ~ -1 + age_year | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy), method = "ML")
    #lme_divstri_age_vardiff_notrisomy = lme(nr_total_corrected ~ age_year + celltype, random = ~ -1 + age_year | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy), method = "ML")
    #lme_divstri_age_vardiff_nocross = lme(nr_total_corrected ~ age_year + trisomy + celltype, random = ~ -1 + age_year | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy), method = "ML")
    #anova(lme_divstri_age_vardiff_nocross, lme_divstri_age_vardiff_notrisomy)
    
    #Sensitivity to removing n points test
    my_vars = c("p-value", "Value")
    sample_col = "sample"
    vals_1 = leave_n_out(lme_divstri_age_vardiff, my_vars, nr_removed = 1, sample = sample_col)
    leave_out_1_p_fig = plot_n_out_pval(vals_1, lme_divstri_age_vardiff)
    leave_out_1_fig = plot_n_out(vals_1, lme_divstri_age_vardiff)
    max_pval_1 = vals_1 %>% dplyr::group_by(variable) %>% dplyr::summarise(max_1 = max(`p-value`))
    
    vals_2 = leave_n_out(lme_divstri_age_vardiff, my_vars, nr_removed = 2, sample = sample_col)
    leave_out_2_p_fig = plot_n_out_pval(vals_2, lme_divstri_age_vardiff)
    leave_out_2_fig = plot_n_out(vals_2, lme_divstri_age_vardiff)
    max_pval_2 = vals_2 %>% dplyr::group_by(variable) %>% dplyr::summarise(max_2 = max(`p-value`))
    
    #vals_3 = leave_n_out(lme_divstri_age_vardiff, my_vars, nr_removed = 3, sample = sample_col)
    #leave_out_3_p_fig = plot_n_out_pval(vals_3, lme_divstri_age_vardiff)
    #leave_out_3_fig = plot_n_out(vals_3, lme_divstri_age_vardiff)
    #max_pval_3 = vals_3 %>% dplyr::group_by(variable) %>% dplyr::summarise(max_3 = max(`p-value`))
    
    max_pvals = cbind(max_pval_1, max_pval_2[,"max_2"])
    write_tsv(max_pvals, paste0(out_dir_base, muttype, "_leave_n_out_maxp.txt"))
    
    pdf(paste0(out_dir_base, muttype, "_leave_n_out.pdf"))
    print(leave_out_1_p_fig)
    print(leave_out_1_fig)
    print(leave_out_2_p_fig)
    print(leave_out_2_fig)
    #print(leave_out_3_p_fig)
    #print(leave_out_3_fig)
    dev.off()
    
    #Compare multiple models
    m1_fixed = get_fixed(lme_divstri_age_vardiff)
    m2_fixed = get_fixed(lme_divstri_age)
    m3_fixed = lme(nr_total_corrected ~ age_year + trisomy, random = ~ -1 + age_year | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy)) %>% get_fixed()
    m4_fixed = lme(nr_total_corrected ~ age_year + trisomy, random = ~ -1 + age_year | fetus_name, data = counts) %>% get_fixed()
    m5_fixed = lme(nr_total_corrected ~ age_year:trisomy, random = ~ -1 + age_year | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy)) %>% get_fixed()
    m6_fixed = lme(nr_total_corrected ~ age_year:trisomy, random = ~ -1 + age_year | fetus_name, data = counts) %>% get_fixed()
    m7_fixed = lme(nr_total_corrected ~ age_year + age_year:trisomy, random = ~ -1 + age_year | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy)) %>% get_fixed()
    m8_fixed = lme(nr_total_corrected ~ age_year + age_year:trisomy, random = ~ -1 + age_year | fetus_name, data = counts) %>% get_fixed()
    
    m_fixed = rbind(m1_fixed, m2_fixed, m3_fixed, m4_fixed, m5_fixed, m6_fixed, m7_fixed, m8_fixed)
    write_tsv(m_fixed, paste0(out_dir_base, muttype, "_divstri_age_models.txt"))
    
    outlier_fig3 = plot(lme_divstri_age, id = 0.05, idLabels = counts$sample)
    outlier_fig4 = plot(lme_divstri_age_vardiff, id = 0.05, idLabels = counts$sample)
    
    #plot the model
    
    if (length(grep("indel", muttype)) == 1){
        ymax = 20
        y_annotation = 12
    } else{
        ymax = 200
        y_annotation = 120
    }
    
    
    counts_plot = counts
    counts$predict = predict(lme_divstri_age_vardiff, level = 0)
    pval = m1_fixed %>% dplyr::filter(variable == "trisomyT21") %>% pull("p-value")
    age_line_fig = ggplot(counts, aes(x = age_year, y = nr_total_corrected)) +
        geom_point(aes(color = trisomy, shape = celltype), size = 4) +
        geom_line(aes(y = predict, color = trisomy, linetype = celltype), size = 1.5) +
        #facet_wrap(. ~ age_cat, scales = "free") +
        labs(title = "LMM containing all data", x = "Age (from conception)", y = "Nr. of corrected mutations") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        annotate("text", x = 0.2, y = y_annotation, label = paste0("P (D21 vs T21): ", round(pval, 3)), size = 6) +
        coord_cartesian(xlim = c(0, 0.3), ylim = c(0, ymax))
    
    #Plot confidence intervals
    t21_ci = intervals(lme_divstri_age_vardiff, which = "fixed")$fixed[3,]
    t21_ci %<>% as.list() %>% as_tibble()
    extra_muts_fig = ggplot(t21_ci, aes(x = "", y = est.)) +
        geom_col(fill = "darkred") +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        labs(title = "Extra mutations from T21", subtitle = "95% CI", x = "", y = "Extra SNVs") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial"))
    
    m_ci = intervals(lme_divstri_age_vardiff, which = "fixed")$fixed
    m_ci %<>% as.data.frame() %>% rownames_to_column("Effect")
    extra_muts_fig2 = ggplot(m_ci, aes(x = Effect, y = est.)) +
        geom_col(fill = "darkred") +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
        labs(title = "Fixed effects", subtitle = "95% CI", x = "Fixed effects", y = "Extra SNVs") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        coord_flip()
    
    #Dont use age at all. This is necessary for the inbulk muts
    lme_divstri_noage_vardiff = lme(nr_total_corrected ~ trisomy*celltype, random = ~ 1 | fetus_name, data = counts, weights = varIdent(form = ~1|trisomy))
    lme_divstri_noage = lme(nr_total_corrected ~ trisomy*celltype, random = ~ 1 | fetus_name, data = counts)
    anova_res = anova(lme_divstri_noage_vardiff, lme_divstri_noage)
    write_tsv(anova_res, paste0(out_dir_base, muttype, "_diffvar_anova_noage.txt"))
    
    lme_sum = summary(lme_divstri_noage_vardiff)
    vals = lme_sum$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable")
    ci = intervals(lme_divstri_noage_vardiff, which = "fixed")$fixed[,c("lower", "upper")]
    vals = cbind(vals, ci) %>% as_tibble()
    
    write_tsv(vals, paste0(out_dir_base, muttype, "_divstri_combi_diffvar_nrmuts_noage.txt"))
    outlier_fig5 = plot(lme_divstri_noage_vardiff, id = 0.05, idLabels = counts$sample)
    
    #Plot no age model
    pval = vals %>% dplyr::filter(variable == "trisomyT21") %>% pull("p-value")
    no_age_fig = ggplot(counts, aes(y = nr_total_corrected, x = trisomy, colour = fetus_name, shape = celltype)) +
        geom_boxplot(outlier.shape = NA, colour = "black", shape = "") +
        geom_quasirandom(size = 3) +
        labs(title = "DivsTri no age", x = "", y = "Nr. of corrected mutations", colour = "Fetus", shape = "Celltype") +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial")) +
        annotate("text", x = 1, y = 20, label = paste0("P (D21 vs T21): ", round(pval, 3)), size = 6)
        
    
    pdf(paste0(out_dir_base, muttype, "_muts_per_fetus.pdf"), width = 10)
    print(hsc_mut_count_figs)
    print(all_mut_count_figs)
    print(si_mut_count_figs)
    print(hscvssi_fig)
    print(trivsdi_fig)
    print(age_line_fig)
    print(extra_muts_fig)
    print(extra_muts_fig2)
    print(no_age_fig)
    print(outlier_fig)
    print(outlier_fig2)
    print(outlier_fig3)
    print(outlier_fig4)
    print(outlier_fig5)
    dev.off()
}

#Read in mutation counts
nr_muts_tb_raw = read_tsv(paste0(out_dir_base, "count_muts.txt"))
nr_muts_tb_raw = left_join(nr_muts_tb_raw, overview_samples[,c("fetus", "sample", "celltype", "surveyed_region_pct", "origin", "gender", "trisomy", "age_year", "age_weeks")], by = c("fetus_name" = "fetus", "sample")) %>% 
    filter(celltype != "Liver")

#Remove liver clones #For ewarts data
#liver_clones_f = grepl("-L-", nr_muts_tb_raw$sample)
#nr_muts_tb_raw = nr_muts_tb_raw[!liver_clones_f,]
#nr_muts_tb_raw$celltype = ifelse(grepl("SI|Si|-I-", nr_muts_tb_raw$sample), "SI", "HSC")
#nr_muts_tb_raw %<>% dplyr::mutate(surveyed_region_pct = ifelse(is.na(surveyed_region_pct), 0.8637, surveyed_region_pct)) #Guess the surveyed region, since callable loci was not yet run on them.

#Calculate correlations between indel and snv mutation counts
nr_muts_tb_raw %<>% mutate(group = paste0(trisomy, "_", celltype))
cor_group = nr_muts_tb_raw %>% dplyr::group_by(group) %>% dplyr::summarise(cor = cor(nr_total_snv, nr_total_indel))
cor_fetus = nr_muts_tb_raw %>% dplyr::group_by(fetus_name) %>% dplyr::summarise(cor = cor(nr_total_snv, nr_total_indel)) %>% dplyr::select(group = fetus_name, cor)
cor_all = nr_muts_tb_raw %>% dplyr::summarise(cor = cor(nr_total_snv, nr_total_indel)) %>% mutate(group = "all")
cors = rbind(cor_all, cor_group, cor_fetus)
write_tsv(cors, paste0(out_dir_base, "snv_indel_cors.txt"))

#Create mutation count plots for different mutation types.
out_dir_divstri = paste0(out_dir_base, "divstri_counts/")
if (!dir.exists(out_dir_divstri)){
    dir.create(out_dir_divstri)
}
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_snv", muttype = "snv", lims = c(150, 700), out_dir_divstri)
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_indel", muttype = "indel", lims = c(15, 50), out_dir_divstri)
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total", muttype = "all", lims = c(150, 700), out_dir_divstri)
# 
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_snv_inbulk", muttype = "snv_inbulk", lims = c(50, 200), out_dir_divstri)
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_indel_inbulk", muttype = "indel_inbulk", lims = c(15, 50), out_dir_divstri)
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_inbulk", muttype = "all_inbulk", lims = c(50, 200), out_dir_divstri)
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_snv_notbulk", muttype = "snv_notbulk", lims = c(150, 700), out_dir_divstri)
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_indel_notbulk", muttype = "indel_notbulk", lims = c(15, 50), out_dir_divstri)
# plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_notbulk", muttype = "all_notbulk", lims = c(150, 700), out_dir_divstri)

#For ewarts data
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_snv", muttype = "snv_inclewart", lims = c(150, 700), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_snv_notbulk", muttype = "snv_notbulk_inclewart", lims = c(150, 700), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_snv_inbulk", muttype = "snv_inbulk_inclewart", lims = c(50, 200), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_indel", muttype = "indel_inclewart", lims = c(15, 50), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_indel_notbulk", muttype = "indel_notbulk_inclewart", lims = c(15, 50), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_indel_inbulk", muttype = "indel_inbulk_inclewart", lims = c(15, 50), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total", muttype = "all_inclewart", lims = c(150, 700), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_notbulk", muttype = "all_notbulk_inclewart", lims = c(150, 700), out_dir_divstri)
plot_divstri_plots(nr_muts_tb_raw, mut = "nr_total_inbulk", muttype = "all_inbulk_inclewart", lims = c(50, 200), out_dir_divstri)


#Compare different clones of n01
nr_muts_tb = nr_muts_tb_raw
nr_muts_tb$nr_muts = nr_muts_tb$nr_total_snv
nr_muts_tb %<>% mutate(nr_total_corrected = nr_muts / surveyed_region_pct)
nr_muts_tb %<>% mutate(muts_per_year = nr_muts / age_year,  muts_per_year_corrected = muts_per_year / surveyed_region_pct)
counts = nr_muts_tb %>% dplyr::filter(fetus_name == "N01")
muts_n01_fig = ggplot(counts, aes(x = origin, y = nr_muts, fill = origin)) +
    geom_boxplot(aes(color = fetus_name), outlier.shape = NA, fill = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 1, dotsize = 1) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(title = "Nr of substitutions in the N01 clones", x = "Clone type", y = "Nr. of substitutions") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    coord_cartesian(ylim = c(0,50))

muts_n01_corrected_fig = ggplot(counts, aes(x = origin, y = nr_total_corrected, fill = origin)) +
    geom_boxplot(aes(color = fetus_name), outlier.shape = NA, fill = NA) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 1, dotsize = 1) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(title = "Nr of substitutions in the N01 clones", x = "Clone type", y = "Nr. of substitutions") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    coord_cartesian(ylim = c(0,50))


#Compare di with tri for both HSC and SI seperately, divided further by mut type
chromosomes = paste0("chr", c(1:22, "X"))
gr_list = apply(overview_samples, 1, function(x) get_all_snvs(x, chromosomes))
grl = GRangesList(gr_list)
factors = overview_samples %>% dplyr::select(fetus, sample, trisomy, celltype, age_year, surveyed_region_pct)
names(grl) = factors$sample
type_occur = mut_type_occurrences(grl, ref_genome)
mut_counts_type = type_occur %>% as.data.frame() %>% rownames_to_column("sample") %>% inner_join(factors, by = "sample") %>% gather(-fetus, -trisomy, -celltype, - sample, -age_year, -surveyed_region_pct, key = "muttype", value = "nr_muts")
mut_counts_type %<>% mutate(muts_per_year = nr_muts / age_year, muts_per_year_corrected = muts_per_year / surveyed_region_pct)
mut_counts_type$celltype = revalue(mut_counts_type$celltype, c("HSC" = "HSPC"))
mut_counts_type$trisomy = revalue(mut_counts_type$trisomy, c("21" = "T21", "FALSE" = "D21"))

maxy = max(mut_counts_type$muts_per_year_corrected)
trivsdi_bytype_fig = ggplot(mut_counts_type, aes(x = trisomy, y = muts_per_year_corrected)) +
    geom_boxplot(outlier.shape = NA, fill = NA) + 
    geom_dotplot(aes(fill = fetus), binaxis = "y", stackdir = "center", binwidth = 4, dotsize = 2) +
    facet_grid(muttype ~ celltype) +
    scale_color_brewer(palette = "Dark2", guide = F) +
    labs(title = "Nr of substitutions per year between HSPC and SI", x = "Cell type", y = "Number of somatic base substitutions\n per genome per year") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    coord_cartesian(ylim = c(0,maxy))

pdf(paste0(out_dir_base, "muts_per_fetus_other.pdf"), width = 10)
trivsdi_bytype_fig
muts_n01_fig
muts_n01_corrected_fig
dev.off()


####_____________Quality of data_____________####
#flagstats_fnames1 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0335/", c("N01", "NC1", "NR1", "NR2"), "/QCStats/flagstat_summary.txt")
#wgsmetrics_fnames1 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0335/", c("N01", "NC1", "NR1", "NR2"), "/QCStats/WGSMetrics_summary.txt")
#flagstats_fnames2 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0349/", c("NR1BMMPPCM16", "N01BMHSPCCB8"), "/QCStats/flagstat_summary.txt")
#wgsmetrics_fnames2 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0349/", c("NR1BMMPPCM16", "N01BMHSPCCB8"), "/QCStats/WGSMetrics_summary.txt")
#flagstats_fnames3 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0367/", c("MH2", "MH3", "N01", "NR1", "NR2"), "/QCStats/flagstat_summary.txt")
#wgsmetrics_fnames3 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0367/", c("MH2", "MH3", "N01", "NR1", "NR2"), "/QCStats/WGSMetrics_summary.txt")
#flagstats_fnames = c(flagstats_fnames1, flagstats_fnames3)
#wgsmetrics_fnames = c(wgsmetrics_fnames1, wgsmetrics_fnames3)

flagstats_fnames1 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0367/", c("MH2"), "/QCStats/flagstat_summary.txt")
wgsmetrics_fnames1 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0367/", c("MH2"), "/QCStats/WGSMetrics_summary.txt")
flagstats_fnames2 = paste0("~/hpc/pmc_vanboxtel/processed/fetal_clones/", c("N01", "NR1", "NR2", "MH3", "OS1"), "/output/QCStats/flagstat_summary.txt")
wgsmetrics_fnames2 = paste0("~/hpc/pmc_vanboxtel/processed/fetal_clones/", c("N01", "NR1", "NR2", "MH3", "OS1"), "/output/QCStats/WGSMetrics_summary.txt")
flagstats_fnames3 = paste0("~/hpc/pmc_vanboxtel/processed/Fetal_Body_Map/", c("E080416", "F100916W17", "F100916W15"), "/QCStats/flagstat_summary.txt")
wgsmetrics_fnames3 = paste0("~/hpc/pmc_vanboxtel/processed/Fetal_Body_Map/", c("E080416", "F100916W17", "F100916W15"), "/QCStats/WGSMetrics_summary.txt")

#flagstats_fnames3 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0409/", c("OS1"), "/QCStats/flagstat_summary.txt")
#wgsmetrics_fnames3 = paste0("~/hpc/pmc_vanboxtel/processed/HMFreg0409/", c("OS1"), "/QCStats/WGSMetrics_summary.txt")
flagstats_fnames = c(flagstats_fnames1, flagstats_fnames2, flagstats_fnames3)
wgsmetrics_fnames = c(wgsmetrics_fnames1, wgsmetrics_fnames2, wgsmetrics_fnames3)


wgsmetrics_fnames = wgsmetrics_fnames[file.exists(wgsmetrics_fnames)]
flagstats_fnames = flagstats_fnames[file.exists(flagstats_fnames)]
wgsmetrics = lapply(wgsmetrics_fnames, function(x)read_tsv(x, col_types = cols(sample = "c"))) %>% bind_rows()
flagstats = lapply(flagstats_fnames, function(x){
    flagstat = read_tsv(x, col_types = cols(.default = "c")) %>% t() %>% as.data.frame(stringsAsFactors = F) %>% rownames_to_column()
    colnames(flagstat) = gsub(" ", ".", flagstat[1,])
    flagstat = flagstat[-1,]
    return(flagstat)
})
flagstats %<>% bind_rows()
qctable = inner_join(wgsmetrics, flagstats, by = "sample") %>% mutate(sample = gsub(".realigned", "", gsub("_dedup", "", sample)), Percentage.reads.mapped = as.numeric(gsub("%", "", Percentage.reads.mapped))/100, Total.number.of.reads = as.numeric(Total.number.of.reads))

qctable$fetus = sapply(qctable$sample, function(i){
    index = grep(i, paste(fetuses_tb$bulk, fetuses_tb$samples, fetuses_tb$other_bulks))
    fetus_name = fetuses_tb$sample_name[index]
    return(fetus_name)
    })

sample_not_used = isEmpty(qctable$fetus)
qctable %<>% dplyr::filter(!sample_not_used) %>% mutate(fetus = unlist(fetus))
#qctable = inner_join(overview_samples, intermed_table, by = "sample")

qctable$mapped.reads = qctable$Total.number.of.reads * qctable$Percentage.reads.mapped
qctable$unmapped.reads = qctable$Total.number.of.reads - qctable$mapped.reads
nrcolors = length(unique(qctable$fetus))
set.seed(001)
colors = sample(colorRampPalette(brewer.pal(8, "Accent"))(nrcolors))
maxyfig1 = 1.05 * max(qctable$MEAN_COVERAGE + qctable$SD_COVERAGE)
MeanCovFig = ggplot(qctable, aes(x = sample, y = MEAN_COVERAGE)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(x = sample, ymin = MEAN_COVERAGE - SD_COVERAGE, ymax = MEAN_COVERAGE + SD_COVERAGE)) + 
    geom_hline(aes(yintercept = 30)) + 
    facet_grid(. ~ fetus, scales = "free", space = "free") + 
    theme_classic() + 
    labs(x = "Sample", y = "Mean coverage +/- SD") + 
    coord_cartesian(ylim = c(0, maxyfig1), expand = F) + 
    theme(axis.text.x = element_text(angle = 90, size = 9, margin = margin(t = 5)),  text = element_text(size=20)) + 
    guides(fill=guide_legend("Family")) +
    scale_fill_manual(values = colors)

qctable2 = qctable %>% dplyr::select(sample, fetus, unmapped.reads, mapped.reads) %>% gather(key = "Legend", value = "nrreads", unmapped.reads, mapped.reads)
qctable2$Legend = gsub("\\.", " ", qctable2$Legend)
maxyfig2 = 1.1 * max(qctable2$nrreads)
NrReadsFig = ggplot(qctable2, aes(x = sample, y = nrreads, fill = Legend)) + 
    facet_grid(. ~ fetus, scales = "free", space = "free") + 
    geom_bar(stat = "identity") + 
    theme_classic() + 
    labs(y = "Number of reads", x = "Sample") + 
    coord_cartesian(ylim = c(0, maxyfig2), expand = F) + 
    theme(axis.text.x = element_text(angle = 90, size = 9, margin = margin(t = 5)),  text = element_text(size=20))

cols_qctable3 = grep("PCT_.*X$", colnames(qctable))
qctable3 = qctable %>% dplyr::select(sample, cols_qctable3, fetus) %>% gather(key = "Xaxis", value = "Coverage", -sample, -fetus)
qctable3$Xaxis = gsub("PCT_", "", qctable3$Xaxis) %>% gsub("X", "", .) %>% as.numeric()
qctable3$Coverage = qctable3$Coverage * 100
n_fetuses = qctable3 %>% dplyr::select(fetus, sample) %>% unique() %>% group_by(fetus) %>% dplyr::summarise(n = n())
nrcolors = length(unique(qctable3$sample))

nrcolors = max(n_fetuses$n)
nrcolors = sum(n_fetuses$n)

colors = sample(colorRampPalette(brewer.pal(8, "Accent"))(nrcolors))

colors_ordered = lapply(n_fetuses$n, function(n){
    cols = colors[1:n]
    return(cols)
}) %>% do.call(c, .)

MinReadCovFig = ggplot(qctable3, aes(x = Xaxis, y=Coverage, color = sample)) + 
    geom_line(aes(linetype = fetus)) + 
    labs(x = "Minimal read coverage", y = "% of genome") + 
    theme_classic() + 
    coord_cartesian(ylim = c(0,100), xlim = c(0,100), expand = F) + 
    theme(legend.text=element_text(size=6), text = element_text(size=20)) + 
    guides(color=guide_legend(title="sample")) +
    scale_color_manual(values = colors_ordered)

qctable4 = qctable
cols_qctable4 = grep("^PCT_EXC_(?!.*TOTAL)", colnames(qctable), perl = T)
#qctable4[,cols_qctable4] = qctable4[,cols_qctable4] * qctable4$Total.number.of.reads
qctable4 = qctable4 %>% dplyr::select(sample, fetus, cols_qctable4) %>% gather(key = "Legend", value = "nrexcluded", -sample, -fetus)
qctable4$Legend = gsub("PCT_EXC_", "", qctable4$Legend)
qctable4$nrexcluded = qctable4$nrexcluded * 100
maxyfig4 = 1.1 * qctable4 %>% dplyr::group_by(sample) %>% dplyr::summarise(sum_sample = sum(nrexcluded)) %>% pull(sum_sample) %>% max()
ReadExclusionFig = ggplot(qctable4, aes(x = sample, y = nrexcluded, fill = Legend)) + 
    geom_bar(stat = "identity") + 
    facet_grid(. ~ fetus, scales = "free", space = "free") + 
    theme_classic() + 
    labs(x = "Sample", y = "% reads excluded") + 
    coord_cartesian(ylim = c(0, maxyfig4), expand = F) + 
    theme(axis.text.x = element_text(angle = 90, size = 9, margin = margin(t = 5)),  text = element_text(size=20)) + 
    guides(fill=guide_legend("Exclusion filter"))

pdf(paste0(out_dir_base, "sequencing_qc_inclewart.pdf"), width = 10)
MeanCovFig
NrReadsFig
MinReadCovFig
ReadExclusionFig
dev.off()


####___________Tree all snvs and indels_______________________####
plot_tree = function(co_occur_m, bulkcols, overview_fetus, out_dir_base, muttype, fetus_name){
    bulkcols = grep(all_bulks, colnames(co_occur_m))
    co_occur_m %<>% as_tibble() %>% dplyr::select(-bulkcols) %>% mutate_all(list(~as.integer)) %>% as.matrix()
    clone_names = colnames(co_occur_m)
    co_occur_m = cbind(co_occur_m, "root" = c(rep(0, nrow(co_occur_m))))
    max_muts = colSums(co_occur_m) %>% max()
    
    #Create tree
    tree = co_occur_m %>% t() %>% dist.gene() %>% nj() #neighbour joining tree construction
    rooted_tree = root(tree, outgroup = "root", resolve.root = T) %>% drop.tip("root", trim.internal = T)
    
    #Color edges based on tissue
    color_map = tibble("origin" = c("Liver", "SI", "BM", "Intestine"), "color" = c("Red", "Darkblue", "Darkred", "Blue"))
    colors_tb = left_join(overview_fetus[,c("origin", "sample")], color_map, by = "origin")
    colors_tb_ordered = enframe(clone_names, name = "nr", value = "sample") %>% inner_join(colors_tb, by = "sample") #Switch color ordering from the order in overview_samples to the order of cols in the vcf.
    color_index = match(colors_tb_ordered$nr, rooted_tree$edge[,2])
    tree_edge_cols = rep("black", nrow(rooted_tree$edge))
    tree_edge_cols[color_index] = colors_tb_ordered$color
    
    #Plot tree
    pdf(paste0(out_dir_base, fetus_name, "/", fetus_name, "_", muttype, "_fulltree.pdf"))
    my_tree_plot = plot(rooted_tree, use.edge.length = T, label.offset = 1, edge.color = tree_edge_cols, edge.width = 2, tip.color = colors_tb_ordered$color)
    axis(1,at =seq(from =0, to = max_muts, by = 1),cex.axis = 1)
    #nodelabels() # we id node label
    edgelabels(rooted_tree$edge.length, bg="black", col="white", font=2)
    dev.off()
}

for (fetus_name in unique(overview_samples$fetus)){
    overview_fetus = overview_samples %>% dplyr::filter(fetus == fetus_name & (shared_vcf_true | indel_shared_true))
    if (nrow(overview_fetus) == 0){
        next
    }
    bulk = overview_fetus$bulk[1]
    other_bulks = overview_fetus$other_bulks[1]
    all_bulks = get_allbulks(other_bulks, bulk)
    
    
    #Get unique snvs The genotype is set to 1 in the called sample and 0 in the not called sample.
    gts_unique = apply(overview_fetus, 1, function(row){
        vcf_fname = row[["unique_vcf"]]
        sample = row[["sample"]]
        if (!file.exists(vcf_fname)){
            return()
        }
        vcf = readVcf(vcf_fname, genome = "hg19")
        gt = geno(vcf)$GT
        sample_col = grep(sample, colnames(gt))
        gt[,sample_col] = 1
        gt[,-sample_col] = 0
        return(gt)
    })
    gt_unique = do.call(rbind, gts_unique)
    
    #Get shared snvs
    if (file.exists(overview_fetus$shared_vcf[1])){
    shared_vcf = readVcf(overview_fetus$shared_vcf[1], genome = "hg19")
    gt_shared = geno(shared_vcf)$GT
    gt_shared[gt_shared == "0/0"] = 0
    gt_shared[gt_shared == "0/1" | gt_shared == "1/1"] = 1
    } else{
        gt_shared = NULL
    }
    #Get unique indels
    if (sum(overview_fetus$indel_uniq_true) != 0){
        gts_unique_indel = apply(overview_fetus, 1, function(row){
            vcf_fname = row[["indel_uniq"]]
            sample = row[["sample"]]
            if (!file.exists(vcf_fname)){
                return()
            }
            vcf = readVcf(vcf_fname, genome = "hg19")
            gt = geno(vcf)$GT
            sample_col = grep(sample, colnames(gt))
            gt[,sample_col] = 1
            gt[,-sample_col] = 0
            return(gt)
        })
        gt_unique_indel = do.call(rbind, gts_unique_indel)
    } else{
        gt_unique_indel = NULL
    }
    
    #Get shared indels
    if (file.exists(overview_fetus$indel_shared[1])){
        shared_vcf_indel = readVcf(overview_fetus$indel_shared[1], genome = "hg19")
        gt_shared_indel = geno(shared_vcf_indel)$GT
        gt_shared_indel[gt_shared_indel == "0/0"] = 0
        gt_shared_indel[gt_shared_indel == "0/1" | gt_shared_indel == "1/1"] = 1
    } else{
        gt_shared_indel = NULL
        }
    
    #Combine all muts in the fetus
    co_occur_m_snv = rbind(gt_unique, gt_shared)
    co_occur_m_indel = rbind(gt_unique_indel, gt_shared_indel)
    co_occur_m_all = rbind(gt_unique, gt_shared, gt_unique_indel, gt_shared_indel)
    plot_tree(co_occur_m_all, bulkcols, overview_fetus, out_dir_base, "all", fetus_name)
    plot_tree(co_occur_m_snv, bulkcols, overview_fetus, out_dir_base, "snv", fetus_name)
    plot_tree(co_occur_m_indel, bulkcols, overview_fetus, out_dir_base, "indel", fetus_name)
    
    
}


# Create vcf and bed with all muts ----------------------------------------
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
    writeVcf(vcf, paste0(out_dir_base, fetus, "/", fetus, "_complete.vcf"))
    gr = granges(vcf)
    bed = tibble("CHROM" = as.vector(seqnames(gr)), "START" = start(gr)-1, "END" = end(gr))
    bed_fname = paste0(out_dir_base, fetus, "/", fetus, "_complete.bed")
    writeLines(paste0("#", paste0(colnames(bed), collapse = "\t")), bed_fname)
    write.table(bed, bed_fname, append = T, quote = F, sep = "\t", col.names = F, row.names = F)
}

####___________Make tree like in axels paper, but with all snvs that are in the bulk_________####
for (fetus_name in unique(overview_samples$fetus)){
    overview_fetus = overview_samples %>% dplyr::filter(fetus == fetus_name & unique_vcf_true & shared_vcf_true)
    if (nrow(overview_fetus) == 0){
        next
    }
    bulk = overview_fetus$bulk[1]
    other_bulks = overview_fetus$other_bulks[1]
    all_bulks = get_allbulks(other_bulks, bulk)
    
    
    
    #Get mutations that are unique to a clone, but occur in the bulk
    vcfs_uniq = lapply(overview_fetus$unique_vcf, function(x) readVcf(x, genome = "hg19"))
    vcf_uniq = do.call(rbind, vcfs_uniq)
   
    vaf_bulk = get_vaf(vcf_uniq, bulk)
    inbulk = vaf_bulk > 0
    #inbulk = gt_uniq[,bulk] == "0/1" | gt_uniq[,bulk] == "1/1"
    vcf_uniq_inbulk = vcf_uniq[inbulk,]
    
    #Combine with the shared mutations
    shared_vcf = readVcf(overview_fetus$shared_vcf[1], genome = "hg19")
    vcf = rbind(shared_vcf, vcf_uniq_inbulk)
    
    gt = geno(vcf)$GT
    bulkcols = grep(all_bulks, colnames(gt))
    co_occur_m = gt %>% as_tibble() %>% dplyr::select(-bulkcols)
    co_occur_m[co_occur_m == "0/0"] = 0
    co_occur_m[co_occur_m == "0/1" | co_occur_m == "1/1"] = 1
    co_occur_m %<>% mutate_all(funs(as.integer)) %>% as.matrix()
    
    #Get vaf
    vaf = get_vaf(vcf, bulk)
    vaf_order = order(vaf, decreasing = T)
    if (!is.na(other_bulks)){
        other_bulks_v = strsplit(other_bulks, ";")[[1]]
        vafs = vector("list", length(other_bulks_v))
        vafs = lapply(other_bulks_v, function(other_bulk){
            vaf_other_bulk = get_vaf(vcf, other_bulk)
        })
        vafs = do.call(cbind, vafs)
        vafs_all = cbind(vaf, vafs)
        colnames(vafs_all) = c(bulk, other_bulks_v)
    } else{
        vafs_all = data.frame(vaf, vaf)
        colnames(vafs_all) = c(bulk, bulk)
    }
    
    
    
    #Order muts based on vaf
    co_occur_ordered = co_occur_m[vaf_order,]
    vaf_ordered = vafs_all[vaf_order,, drop = F]
    
    #Name muts based on letters or numbers if there are too many muts
    if (nrow(co_occur_ordered) > 26){
        rownames(co_occur_ordered) = seq(1, nrow(co_occur_ordered))
        rownames(vaf_ordered) = rownames(co_occur_ordered)
    } else{
        rownames(co_occur_ordered) = letters[seq(1, nrow(co_occur_ordered))]
        rownames(vaf_ordered) = rownames(co_occur_ordered)
    }
    
    #Count unique mutations that are not in the bulk
    vcf_uniq_notbulk = vcf_uniq[!inbulk,]
    gt_uniq_notbulk = geno(vcf_uniq_notbulk)$GT[,-bulkcols]
    uniq_notbulk = gt_uniq_notbulk == "0/1" | gt_uniq_notbulk == "1/1"
    nr_uniq_notbulk = colSums(uniq_notbulk)
    clust = t(co_occur_ordered) %>% dist(method = "binary") %>% hclust(method = "ward.D2")
    nr_uniq_notbulk_sort = nr_uniq_notbulk[clust$order] %>% enframe()
    nr_uniq_notbulk_sort$name = factor(nr_uniq_notbulk_sort$name, levels = nr_uniq_notbulk_sort$name)
    maxy = max(nr_uniq_notbulk_sort$value) * 1.2
    
    
    pdf(paste0(out_dir_base, fetus_name, "/axeltree.pdf"))
    my_palette <- colorRampPalette(c("#eeeeee", "#333333"))(n = 3)
    col_breaks <- c(-0.5, 0.5, 1.5, 2.5)
    heatmap.2(co_occur_ordered, distfun=function(x) dist(x,method = 'binary'),
              hclustfun=function(x) hclust(x,method = 'ward.D2'),
              Rowv = F, labCol = FALSE,
              key=FALSE, # remove color key
              col=my_palette, breaks=col_breaks,
              colsep=c(21, 44, 57, 81, 101), sepcolor="black",
              trace="none", density.info="none",
              lwid=c(0.1,2), lhei=c(0.5,5), dendrogram = "column")
    
    
    my_palette <- colorRampPalette(c("white", "#0068b7"))(n = 299)
    heatmap.2(as.matrix(vaf_ordered),
              Rowv = F,  dendrogram = "column",
              col=my_palette,
              trace="none", density.info="none")
    
    nr_uniq_notbulk_fig = ggplot(nr_uniq_notbulk_sort, aes(x = name, ymin = 0, ymax = value)) +
        geom_linerange() +
        scale_y_reverse() +
        theme_classic() + 
        labs(x = "Sample", y = "Unique muts not in vaf") + 
        theme(text = element_text(size=20))
    print(nr_uniq_notbulk_fig)
    
    dev.off()
}

####________Look at early mutations (shared + unique, but in bulk)_________####



#Heplper function that takes a shared vcf and returns the mutations that are unique to the specified samples. It returns the mutations as a granges object.
shared_muts_only_samples = function(shared_vcf, samples){
    samples_grep = paste0(samples, collapse = "$|^") %>% paste0("^", ., "$")
    gt = geno(shared_vcf)$GT
    col_samples = grep(samples_grep, colnames(gt))
    muts_othersamples = gt[, -col_samples, drop = F] == "0/1" | gt[, -col_samples, drop = F] == "1/1"
    muts_only_selected_samples_f = rowSums(muts_othersamples) == 0
    shared_vcf = shared_vcf[muts_only_selected_samples_f]
    shared_gr = granges(shared_vcf)
    return(shared_gr)
}


#Function to get muts that are not in the bulk. It gets all unique muts. It also gets shared muts if they are only shared by the specified samples.
get_notbulk_gr = function(overview_fetus){
    uniq_vcf_fnames = overview_fetus$uniq_vcf_inbulk[overview_fetus$uniq_vcf_inbulk_true]
    if (length(uniq_vcf_fnames > 0)){
        uniq_vcfs = lapply(uniq_vcf_fnames, function(x) readVcf(x, genome = "hg19"))
        uniq_grl = lapply(uniq_vcfs, granges)
        uniq_gr = do.call(c, uniq_grl) %>% sort()
    } else{
        uniq_gr = GRanges()
    }
    samples = overview_fetus$sample
    shared_vcf_fnames = overview_fetus$shared_notbulk_vcf[overview_fetus$shared_vcf_notbulk_true] %>% unique()
    if (length(shared_vcf_fnames > 0)){
        shared_vcfs = lapply(shared_vcf_fnames, function(x) readVcf(x, genome = "hg19"))
        shared_grl = lapply(shared_vcfs, function(x) shared_muts_only_samples(x, samples))
        shared_gr = do.call(c, shared_grl) %>% sort()
    } else{
        shared_gr = GRanges()
    }
    gr = c(uniq_gr, shared_gr) %>% sort()
    return(gr)
}

####_______________Mut patterns____________________####

#Create directory for mutational patterns analyses
out_dir_mut = paste0(out_dir_base, "mutpatterns/")
if (!dir.exists(out_dir_mut)){
    dir.create(out_dir_mut)
}
overview_samples = overview_samples %>% dplyr::filter(celltype != "Liver")

chromosomes = paste0("chr", c(1:22, "X"))
overview_samples_bm = overview_samples %>% dplyr::filter(origin == "BM")
grl = apply(overview_samples_bm, 1, function(x) get_all_snvs(x, chromosomes)) %>% 
    GRangesList()

mut_mat = mut_matrix2(grl, ref_genome) %>% pool_mut_mat(rep("BM", 2))
cosmic_sig_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_working_signatures_incl_hspc.txt"
signatures = read.table(cosmic_sig_fname, sep = "\t", header = T)
signatures = as.matrix(signatures[,-c(1,2)])

figs_l = refit_per_sample(mut_mat)
pdf(paste0(out_dir_mut, "BM_betterrefit.pdf"))
print(figs_l)
dev.off()

#Function to get the unique muts for the specified samples, that are not present in the bulk. The function returns a grl object.
vcf_to_gr = function(vcf_fname){
    vcf = readVcf(vcf_fname, genome = "hg19")
    gr = granges(vcf)
    seqlevelsStyle(gr) = "UCSC"
    return(gr)
}

exp_name = "not_bulk_uniq_blood_di_vs_tri"
chromosomes = paste0("chr", c(1:22, "X"))
overview_blood = overview_samples %>% dplyr::filter(celltype == "HSC" & uniq_vcf_notbulk_true)
gr_list = lapply(overview_blood$uniq_vcf_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
names(grl) = overview_blood$sample
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = overview_blood %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

#Get unique not bulk indels
overview_blood = overview_samples %>% dplyr::filter(celltype == "HSC" & indel_uniq_notbulk_true)
gr_list = lapply(overview_blood$indel_uniq_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
names(grl) = overview_blood$sample
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = overview_blood %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
grl_exp = pool_grl_exp(grl, factors$experiment)

grl_exp = get_indel_context(grl_exp, ref_genome)
indel_counts = count_indel_contexts(grl_exp)
main_indel_context_fig = plot_main_indel_contexts(indel_counts)
indel_context_fig = plot_indel_contexts(indel_counts)

pdf(paste0(out_dir_mut, exp_name, "_indel_contexts.pdf"), width = 12)
main_indel_context_fig
indel_context_fig
dev.off()

exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "_indels_uniq_chisq.txt"))

do_indel_refitting(indel_counts, out_dir_mut, exp_name)

exp_name = "not_bulk_uniq_blood_di_vs_tri_noOS1B21"
chromosomes = paste0("chr", c(1:22, "X"))
overview_blood = overview_samples %>% dplyr::filter(celltype == "HSC" & uniq_vcf_notbulk_true & sample != "OS1LIHSPCC2B21")
gr_list = lapply(overview_blood$uniq_vcf_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
names(grl) = overview_blood$sample
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = overview_blood %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

exp_name = "not_bulk_uniq_blood_tmd_vs_tri"
chromosomes = paste0("chr", c(1:22, "X"))
overview_blood = overview_samples %>% dplyr::filter(celltype == "HSC" & uniq_vcf_notbulk_true & trisomy == "21")
gr_list = lapply(overview_blood$uniq_vcf_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
gr = unlist(grl)
mcols(gr) = mcols(gr)[c("REF", "ALT")]
gr_tmd = readRDS("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/tmd_gr.rds")
mcols(gr_tmd) = mcols(gr_tmd)[c("REF", "ALT")]
grl = GRangesList(gr, gr_tmd)
names(grl) = c("T21", "TMD")
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = tibble("fetus" = c("T21", "TMD"), "callableloci_merged" = c(NA, NA), "callableloci_merged_true" = c(F, F), experiment = c("T21", "TMD"))

do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

exp_name = "not_bulk_uniq_blood_tmd_vs_tri_noOS1B21"
chromosomes = paste0("chr", c(1:22, "X"))
overview_blood = overview_samples %>% dplyr::filter(celltype == "HSC" & uniq_vcf_notbulk_true & trisomy == "21" & sample != "OS1LIHSPCC2B21")
gr_list = lapply(overview_blood$uniq_vcf_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
gr = unlist(grl)
mcols(gr) = mcols(gr)[c("REF", "ALT")]
gr_tmd = readRDS("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/Before_final_version/tmd_gr.rds")
mcols(gr_tmd) = mcols(gr_tmd)[c("REF", "ALT")]
grl = GRangesList(gr, gr_tmd)
names(grl) = c("T21", "TMD")
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = tibble("fetus" = c("T21", "TMD"), "callableloci_merged" = c(NA, NA), "callableloci_merged_true" = c(F, F), experiment = c("T21", "TMD"))

do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)


tmd_type = read_tsv("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/not_bulk_uniq_blood_tmd_vs_tri_noOS1B21_type_occur.txt") %>% dplyr::filter(experiment == "TMD")
divstri_type = read_tsv("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/di_vs_tri_inclewart_noOS1B21_type_occur.txt")
all_type = rbind(divstri_type, tmd_type)
spectra_fig = all_type %>% column_to_rownames("experiment") %>% plot_spectrum(., by = rownames(.)) + theme_classic()
pdf("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/d21vst21vstmd_spectrum.pdf", width = 10)
spectra_fig
dev.off()

exp_name = "not_bulk_uniq_SI_di_vs_tri"
chromosomes = paste0("chr", c(1:22, "X"))
overview_si = overview_samples %>% dplyr::filter(celltype == "SI" & uniq_vcf_notbulk_true)
gr_list = lapply(overview_si$uniq_vcf_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
names(grl) = overview_si$sample
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = overview_si %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

#Get unique not bulk indels
overview_si = overview_samples %>% dplyr::filter(celltype == "SI" & indel_uniq_notbulk_true)
gr_list = lapply(overview_si$indel_uniq_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
names(grl) = overview_si$sample
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = overview_si %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
grl_exp = pool_grl_exp(grl, factors$experiment)

grl_exp = get_indel_context(grl_exp, ref_genome)
indel_counts = count_indel_contexts(grl_exp)
main_indel_context_fig = plot_main_indel_contexts(indel_counts)
indel_context_fig = plot_indel_contexts(indel_counts)

pdf(paste0(out_dir_mut, exp_name, "_indel_contexts.pdf"), width = 12)
main_indel_context_fig
indel_context_fig
dev.off()

exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "_indels_uniq_chisq.txt"))

do_indel_refitting(indel_counts, out_dir_mut, exp_name)

exp_name = "not_bulk_uniq_SI_di_vs_tri_noOS1B21"
chromosomes = paste0("chr", c(1:22, "X"))
overview_si = overview_samples %>% dplyr::filter(celltype == "SI" & uniq_vcf_notbulk_true & sample != "OS1LIHSPCC2B21")
gr_list = lapply(overview_si$uniq_vcf_notbulk, vcf_to_gr)
grl = GRangesList(gr_list)
names(grl) = overview_si$sample
seqlevels(grl, pruning.mode = "fine") = chromosomes
factors = overview_si %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

# exp_name = "not_bulk_uniq_di_vs_tri_inclewart"
# chromosomes = paste0("chr", c(1:22, "X"))
# overview_fetus1 = overview_samples %>% dplyr::filter(celltype %in% c("HSC", "SI") & uniq_vcf_notbulk_true)
# grl_blood = get_unique_muts(overview_fetus1)
# overview_fetus2 = overview_samples %>% dplyr::filter(celltype %in% c("HSC", "SI") & uniq_vcf_notbulk_true)
# grl_blood_tri = get_unique_muts(overview_fetus2)
# overview_fetus_all = bind_rows(overview_fetus1, overview_fetus2)
# factors = tibble("fetus" = overview_fetus_all$fetus, "callableloci_merged" = overview_fetus_all$callableloci_merged, "callableloci_merged_true" = overview_fetus_all$callableloci_merged_true, "experiment" = as.factor(c(rep("Disomy_blood", nrow(overview_fetus1)), rep("Trisomy_blood", nrow(overview_fetus2)))))
# grl = c(grl_blood, grl_blood_tri)
# do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
# do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
# do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)


#Function to get all the muts for the specified fetuses. This includes shared muts ect. This way each mutation is only counted ONCE
get_all_muts = function(overview_fetus){
    shared_vcfs = overview_fetus$shared_vcf[overview_fetus$shared_vcf_true] %>% unique()
    unique_vcfs = overview_fetus$unique_vcf[overview_fetus$unique_vcf_true]
    all_vcfs = c(shared_vcfs, unique_vcfs)
    grl = lapply(all_vcfs, function(x){
        vcf = readVcf(x, genome = "hg19")
        gr = granges(vcf)
        return(gr)
    })
    gr = do.call(c, grl)
    gr %<>% sort()
    seqlevelsStyle(gr) = "UCSC"
    return(gr)
}

exp_name = "di_vs_tri_inclewart"
chromosomes = paste0("chr", c(1:22, "X"))
gr_list = apply(overview_samples, 1, function(x) get_all_snvs(x, chromosomes))
fetuses_v = unique(overview_samples$fetus)
grl_fetal = vector("list", length(fetuses_v))
for (i in 1:length(fetuses_v)){#Pool per fetus. Only count shared muts once.
    fetus = fetuses_v[i]
    j = grep(fetus, overview_samples$fetus)
    gr_fetus = do.call("c", gr_list[j]) %>% unique()
    grl_fetal[[i]] = gr_fetus
}
grl = GRangesList(grl_fetal)
names(grl) = fetuses_v

# gr_list = lapply(unique(overview_samples$fetus), function(x){
#     overview_fetus = overview_samples %>% dplyr::filter(fetus == !!x)
#     gr = get_all_muts(overview_fetus)
#     return(gr)
# })
# grl = GRangesList(gr_list)
# names(grl) = unique(overview_samples$fetus)

overview_fetuses = overview_samples %>% dplyr::filter(!duplicated(fetus)) %>% mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), c("T21", "D21")))
factors = tibble("fetus" = overview_fetuses$fetus, "callableloci_merged" = paste0(out_dir_base, overview_fetuses$fetus, "/callableloci/", overview_fetuses$bulk, "CallableLoci_autosomal_X.bed"), "callableloci_merged_true" = TRUE, "experiment" = as.factor(overview_fetuses$trisomy))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

#Get all indels per fetus. Each indel is only counted once.
gr_list = apply(overview_samples, 1, function(x) get_all_indels(x, chromosomes))
fetuses_v = unique(overview_samples$fetus)
grl_fetal = vector("list", length(fetuses_v))
for (i in 1:length(fetuses_v)){#Pool per fetus. Only count shared muts once.
    fetus = fetuses_v[i]
    j = grep(fetus, overview_samples$fetus)
    gr_fetus = do.call("c", gr_list[j]) %>% unique()
    grl_fetal[[i]] = gr_fetus
}

grl = GRangesList(grl_fetal)
names(grl) = fetuses_v

# gr_list = lapply(factors$fetus, function(x) get_all_indels_fetus(x, overview_samples, chromosomes))
# grl = GRangesList(gr_list)
# names(grl) = factors$fetus
grl_exp = pool_grl_exp(grl, factors$experiment)

grl_exp = get_indel_context(grl_exp, ref_genome)
indel_counts = count_indel_contexts(grl_exp)
main_indel_context_fig = plot_main_indel_contexts(indel_counts)
indel_context_fig = plot_indel_contexts(indel_counts)

pdf(paste0(out_dir_mut, exp_name, "_indel_contexts.pdf"), width = 12)
main_indel_context_fig
indel_context_fig
dev.off()

exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "_indels_uniq_chisq.txt"))

do_indel_refitting(indel_counts, out_dir_mut, exp_name)


exp_name = "di_vs_tri_inclewart_noOS1B21"
chromosomes = paste0("chr", c(1:22, "X"))
overview_samples_noOSB21 = overview_samples %>% dplyr::filter(sample != "OS1LIHSPCC2B21")
gr_list = apply(overview_samples_noOSB21, 1, function(x) get_all_snvs(x, chromosomes))
fetuses_v = unique(overview_samples_noOSB21$fetus)
grl_fetal = vector("list", length(fetuses_v))
for (i in 1:length(fetuses_v)){#Pool per fetus. Only count shared muts once.
    fetus = fetuses_v[i]
    j = grep(fetus, overview_samples_noOSB21$fetus)
    gr_fetus = do.call("c", gr_list[j]) %>% unique()
    grl_fetal[[i]] = gr_fetus
}
grl = GRangesList(grl_fetal)
names(grl) = fetuses_v


overview_fetuses = overview_samples_noOSB21 %>% dplyr::filter(!duplicated(fetus)) %>% mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), c("T21", "D21")))
factors = tibble("fetus" = overview_fetuses$fetus, "callableloci_merged" = paste0(out_dir_base, overview_fetuses$fetus, "/callableloci/", overview_fetuses$bulk, "CallableLoci_autosomal_X.bed"), "callableloci_merged_true" = TRUE, "experiment" = as.factor(overview_fetuses$trisomy))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)





# exp_name = "blood_di_vs_tri_inclewart"
# overview_samples_blood = overview_samples %>% dplyr::filter(celltype == "HSC")
# chromosomes = paste0("chr", c(1:22, "X"))
# gr_list = lapply(unique(overview_samples_blood$fetus), function(x){
#     overview_fetus = overview_samples_blood %>% dplyr::filter(fetus == !!x)
#     gr = get_all_muts(overview_fetus)
#     return(gr)
# })
# grl = GRangesList(gr_list)
# names(grl) = unique(overview_samples_blood$fetus)
# 
# overview_fetuses = overview_samples_blood %>% dplyr::filter(!duplicated(fetus)) %>% mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), c("T21", "D21")))
# factors = tibble("fetus" = overview_fetuses$fetus, "callableloci_merged" = paste0(out_dir_base, overview_fetuses$fetus, "/callableloci/", overview_fetuses$bulk, "CallableLoci_CALLABLE.bed"), "callableloci_merged_true" = TRUE, "experiment" = as.factor(overview_fetuses$trisomy))
# do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
# do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
# do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)
# 
# #Get all indels per fetus. Each indel is only counted once.
# gr_list = lapply(factors$fetus, function(x) get_all_indels_fetus(x, overview_samples_blood, chromosomes))
# grl = GRangesList(gr_list)
# names(grl) = factors$fetus
# grl_exp = pool_grl_exp(grl, factors$experiment)
# 
# grl_exp = get_indel_context(grl_exp, ref_genome)
# indel_counts = count_indel_contexts(grl_exp)
# main_indel_context_fig = plot_main_indel_contexts(indel_counts)
# indel_context_fig = plot_indel_contexts(indel_counts)
# 
# pdf(paste0(out_dir_mut, exp_name, "_indel_contexts.pdf"), width = 12)
# main_indel_context_fig
# indel_context_fig
# dev.off()
# 
# exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
# indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
# chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
# chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
# write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "_indels_uniq_chisq.txt"))
# 
# do_indel_refitting(indel_counts, out_dir_mut, exp_name)
# 
# exp_name = "si_di_vs_tri_inclewart"
# overview_samples_si = overview_samples %>% dplyr::filter(celltype == "SI")
# chromosomes = paste0("chr", c(1:22, "X"))
# gr_list = lapply(unique(overview_samples_si$fetus), function(x){
#     overview_fetus = overview_samples_si %>% dplyr::filter(fetus == !!x)
#     gr = get_all_muts(overview_fetus)
#     return(gr)
# })
# grl = GRangesList(gr_list)
# names(grl) = unique(overview_samples_si$fetus)
# 
# overview_fetuses = overview_samples_si %>% dplyr::filter(!duplicated(fetus)) %>% mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), c("T21", "D21")))
# factors = tibble("fetus" = overview_fetuses$fetus, "callableloci_merged" = paste0(out_dir_base, overview_fetuses$fetus, "/callableloci/", overview_fetuses$bulk, "CallableLoci_CALLABLE.bed"), "callableloci_merged_true" = TRUE, "experiment" = as.factor(overview_fetuses$trisomy))
# do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
# do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
# do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)
# 
# #Get all indels per fetus. Each indel is only counted once.
# gr_list = lapply(factors$fetus, function(x) get_all_indels_fetus(x, overview_samples_si, chromosomes))
# grl = GRangesList(gr_list)
# names(grl) = factors$fetus
# grl_exp = pool_grl_exp(grl, factors$experiment)
# 
# grl_exp = get_indel_context(grl_exp, ref_genome)
# indel_counts = count_indel_contexts(grl_exp)
# main_indel_context_fig = plot_main_indel_contexts(indel_counts)
# indel_context_fig = plot_indel_contexts(indel_counts)
# 
# pdf(paste0(out_dir_mut, exp_name, "_indel_contexts.pdf"), width = 12)
# main_indel_context_fig
# indel_context_fig
# dev.off()
# 
# exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
# indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
# chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
# chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
# write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "_indels_uniq_chisq.txt"))
# 
# do_indel_refitting(indel_counts, out_dir_mut, exp_name)


exp_name = "hsc_all_muts_per_cell_divstri"
chromosomes = paste0("chr", c(1:22, "X"))
overview_hsc = overview_samples %>% dplyr::filter(origin %in% c("Liver", "BM"))
gr_list = apply(overview_hsc, 1, function(x) get_all_snvs(x, chromosomes))
grl = GRangesList(gr_list)
names(grl) = overview_hsc$sample
factors = overview_hsc %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

#Get_dbs
dbs_grl = lapply(grl, get_dbs) %>% GRangesList()
dbs_grl_exp = pool_grl_exp(dbs_grl, factors$experiment)
dbs_grl_exp = lapply(dbs_grl_exp, set_context_dbs) %>% GRangesList()
dbs_counts = count_dbs_contexts(dbs_grl_exp)
main_dbs_context_fig = plot_main_dbs_contexts(dbs_counts)
dbs_context_fig = plot_dbs_contexts(dbs_counts)

pdf(paste0(out_dir_mut, exp_name, "_dbs_contexts.pdf"), width = 12)
main_dbs_context_fig
dbs_context_fig
dev.off()

#dbs only count each dbs once
dbs_grl_fetus = pool_grl_exp(dbs_grl, factors$fetus)
dbs_grl_uniq = lapply(dbs_grl_fetus, unique) %>% GRangesList()
names_grl = names(dbs_grl_uniq)
factors_dbs = factors %>% dplyr::filter(!duplicated(fetus)) %>% arrange(fetus, names_grl)
dbs_grl_exp = pool_grl_exp(dbs_grl_uniq, factors_dbs$experiment)

dbs_grl_exp = lapply(dbs_grl_exp, set_context_dbs) %>% GRangesList()
dbs_counts = count_dbs_contexts(dbs_grl_exp)
main_dbs_context_fig = plot_main_dbs_contexts(dbs_counts)
dbs_context_fig = plot_dbs_contexts(dbs_counts)

pdf(paste0(out_dir_mut, exp_name, "_dbs_contexts_nondupli.pdf"), width = 12)
main_dbs_context_fig
dbs_context_fig
dev.off()


#Get indels
gr_list = apply(overview_hsc, 1, function(x) get_all_indels(x, chromosomes))
grl = GRangesList(gr_list)
names(grl) = overview_hsc$sample
grl_exp = pool_grl_exp(grl, factors$experiment)

grl_exp = get_indel_context(grl_exp, ref_genome)
indel_counts = count_indel_contexts(grl_exp)
main_indel_context_fig = plot_main_indel_contexts(indel_counts)
indel_context_fig = plot_indel_contexts(indel_counts)

pdf(paste0(out_dir_mut, exp_name, "_indel_contexts.pdf"), width = 12)
main_indel_context_fig
indel_context_fig
dev.off()

exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "indels_chisq.txt"))

do_indel_refitting(indel_counts, out_dir_mut, exp_name)




exp_name = "all_muts_per_cell_divstri"
chromosomes = paste0("chr", c(1:22, "X"))
gr_list = apply(overview_samples, 1, function(x) get_all_snvs(x, chromosomes))
grl = GRangesList(gr_list)
names(grl) = overview_samples$sample
factors = overview_samples %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy) %>% mutate(experiment = mapvalues(as.factor(experiment), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)


#Get_dbs
dbs_grl = lapply(grl, get_dbs) %>% GRangesList()
dbs_grl_exp = pool_grl_exp(dbs_grl, factors$experiment)
dbs_grl_exp = lapply(dbs_grl_exp, set_context_dbs) %>% GRangesList()
dbs_counts = count_dbs_contexts(dbs_grl_exp)
main_dbs_context_fig = plot_main_dbs_contexts(dbs_counts)
dbs_context_fig = plot_dbs_contexts(dbs_counts)

pdf(paste0(out_dir_mut, exp_name, "_dbs_contexts.pdf"), width = 12)
main_dbs_context_fig
dbs_context_fig
dev.off()

counts = as.list(dbs_grl) %>% sapply(length)
wilcox.test(counts ~ overview_samples$trisomy, exact = F)
glm(counts ~ overview_samples$trisomy, family="poisson") %>% summary()

#dbs only count each dbs once
dbs_grl_fetus = pool_grl_exp(dbs_grl, factors$fetus)
dbs_grl_uniq = lapply(dbs_grl_fetus, unique) %>% GRangesList()
names_grl = names(dbs_grl_uniq)
factors_dbs = factors %>% dplyr::filter(!duplicated(fetus)) %>% arrange(fetus, names_grl)
dbs_grl_exp = pool_grl_exp(dbs_grl_uniq, factors_dbs$experiment)

dbs_grl_exp = lapply(dbs_grl_exp, set_context_dbs) %>% GRangesList()
dbs_counts = count_dbs_contexts(dbs_grl_exp)
main_dbs_context_fig = plot_main_dbs_contexts(dbs_counts)
dbs_context_fig = plot_dbs_contexts(dbs_counts)

pdf(paste0(out_dir_mut, exp_name, "_dbs_contexts_nondupli.pdf"), width = 12)
main_dbs_context_fig
dbs_context_fig
dev.off()

dbs_gr_uniq = unlist(dbs_grl_uniq)
mut_names = names(dbs_gr_uniq)
sample_names = str_replace(mut_names, "^.*\\.(.*)\\..*$", "\\1")
counts2 = table("sample" = sample_names) %>% as.data.frame()
count_table = overview_samples[,c("sample", "trisomy")]
counts2 = left_join(count_table, counts2, by = "sample", fill = 0)
counts2$Freq = replace_na(counts2$Freq, 0)
wilcox.test(counts2$Freq ~ counts2$trisomy, exact = F)


#Get indels
gr_list = apply(overview_samples, 1, function(x) get_all_indels(x, chromosomes))
grl = GRangesList(gr_list)
names(grl) = overview_samples$sample
grl_exp = pool_grl_exp(grl, factors$experiment)

grl_exp = get_indel_context(grl_exp, ref_genome)
indel_counts = count_indel_contexts(grl_exp)
main_indel_context_fig = plot_main_indel_contexts(indel_counts)
indel_context_fig = plot_indel_contexts(indel_counts)

pdf(paste0(out_dir_mut, exp_name, "_indel_contexts.pdf"), width = 12)
main_indel_context_fig
indel_context_fig
dev.off()

exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "indels_chisq.txt"))

do_indel_refitting(indel_counts, out_dir_mut, exp_name)

# #Get indels uniq
# gr_list = apply(overview_samples, 1, function(x) get_uniq_indels(x, chromosomes))
# grl = GRangesList(gr_list)
# names(grl) = overview_samples$sample
# grl_exp = pool_grl_exp(grl, factors$experiment)
# 
# grl_exp = get_indel_context(grl_exp, ref_genome)
# indel_counts = count_indel_contexts(grl_exp)
# main_indel_context_fig = plot_main_indel_contexts(indel_counts)
# indel_context_fig = plot_indel_contexts(indel_counts)
# 
# pdf(paste0(out_dir_mut, exp_name, "_uniq_indel_contexts.pdf"), width = 12)
# main_indel_context_fig
# indel_context_fig
# dev.off()
# 
# exist_f = indel_counts[,-c(1:2)] %>% rowSums() != 0
# indel_counts_exist = indel_counts %>% dplyr::filter(exist_f)
# chisq_res = indel_counts_exist %>% dplyr::select(-muttype, -muttype_sub) %>% chisq.test(simulate.p.value = T)
# chisq_tb = tibble("pvalue" = chisq_res$p.value, "x_squared" = chisq_res$statistic)
# write_tsv(chisq_tb, paste0(out_dir_mut, exp_name, "_indels_uniq_chisq.txt"))
# 
# do_indel_refitting(indel_counts, out_dir_mut, paste0(exp_name, "_uniq"))



exp_name = "all_muts_per_cell_hscvssi_healthy"
chromosomes = paste0("chr", c(1:22, "X"))
overview_healthy = overview_samples %>% dplyr::filter(trisomy == "FALSE")
gr_list = apply(overview_healthy, 1, function(x) get_all_snvs(x, chromosomes))
grl = GRangesList(gr_list)
names(grl) = overview_healthy$sample
factors = overview_healthy %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = celltype)
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)


#Function to get the unique muts for the specified samples, that are not present in the bulk. The function returns a grl object.
get_unique_muts = function(overview_fetus){
    gr_list = lapply(overview_fetus$uniq_vcf_notbulk, function(x) {
        vcf = readVcf(x, genome = "hg19")
        gr = granges(vcf)
        return(gr)
    })
    grl = GRangesList(gr_list)
    names(grl) = overview_fetus$sample
    seqlevelsStyle(grl) = "UCSC"
    return(grl)
}

exp_name = "hscvssi_healthy"
overview_samples_healthy = overview_samples %>% dplyr::filter(trisomy == "FALSE")
chromosomes = paste0("chr", c(1:22, "X"))
grl = get_unique_muts(overview_samples_healthy)

factors = overview_samples_healthy %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = celltype)
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)

exp_name = "fetalvsadult"
chromosomes = paste0("chr", 1:22)
overview_samples_blood = overview_samples %>% dplyr::filter(celltype == "HSC" & trisomy == "FALSE")
gr_list = apply(overview_samples_blood, 1, function(x) get_all_snvs(x, chromosomes))
fetuses_v = unique(overview_samples_blood$fetus)
grl_fetal = vector("list", length(fetuses_v))

for (i in 1:length(fetuses_v)){#Pool per fetus. Only count shared muts once.
    fetus = fetuses_v[i]
    j = grep(fetus, overview_samples_blood$fetus)
    gr_fetus = do.call("c", gr_list[j]) %>% unique()
    grl_fetal[[i]] = gr_fetus
}

grl_fetal = GRangesList(grl_fetal)
names(grl_fetal) = fetuses_v

overview_fetuses = overview_samples_blood %>% dplyr::filter(!duplicated(fetus))
factors_fetal = tibble("fetus" = overview_fetuses$fetus, "callableloci_merged" = paste0(out_dir_base, overview_fetuses$fetus, "/callableloci/", overview_fetuses$bulk, "CallableLoci_autosomal.bed"), "callableloci_merged_true" = TRUE, "experiment" = "Fetal")

#Get vcf files from axels data
files_axel =list.files(vcfs_axel_dir, full.names = T)
vcfs_axel = files_axel[grepl("_MQ60.vcf", files_axel)]
samples_axel = gsub(".*/(.*)_.*_Q100.*", "\\1", vcfs_axel)

#Get callable loci files from axels data and ensure its sorted in the same way as the vcfs.
callable_files_axel_all = list.files(callable_axel_dir, full.names = T)
callable_files_axel = callable_files_axel_all[grepl("_merged_autosomal.bed", callable_files_axel_all)]
samples_callable_axel = gsub(".*/(.*)_.*_CallableLoci.*", "\\1", callable_files_axel)
donor_callable_axel = gsub(".*/.*_(.*)MSC.*_CallableLoci.*", "\\1", callable_files_axel)
callable_tb = tibble(fetus = donor_callable_axel, callableloci_merged = callable_files_axel, callableloci_merged_true = TRUE, experiment = "Adult", samples_callable_axel) #Create tibble, that can be used as the factors for Do_mutational patterns.
factors_adult = enframe(samples_axel) %>% left_join(callable_tb, by = c("value" = "samples_callable_axel")) %>% dplyr::select(-name, -value) #Sort the factors tibble based on the sample order from the vcfs


cb_samples = c("CB112MPPC17", "CB112MPPC9", "CB112MPPO6", "10817MPPA4")
cb_index = samples_axel %in% cb_samples
grl_axel = read_vcfs_as_granges(vcfs_axel[!cb_index], samples_axel[!cb_index], genome)
factors_adult = factors_adult[!cb_index,]

grl = c(grl_fetal, grl_axel)
factors = rbind(factors_fetal, factors_adult)
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)


exp_name = "all_muts_per_cell_fetalvsadult"
chromosomes = paste0("chr", 1:22)
overview_hsc_di = overview_samples %>% dplyr::filter(trisomy == "FALSE" & origin %in% c("Liver", "BM"))
gr_list = apply(overview_hsc_di, 1, function(x) get_all_snvs(x, chromosomes))
grl_fetal = GRangesList(gr_list)
names(grl_fetal) = overview_hsc_di$sample
factors_fetal = overview_hsc_di %>% dplyr::select(fetus, callableloci_merged = callableloci_merged_auto, callableloci_merged_true = callableloci_merged_auto_true) %>% mutate(experiment = "Fetal")

#Get vcf files from axels data
files_axel =list.files(vcfs_axel_dir, full.names = T)
vcfs_axel = files_axel[grepl("_MQ60.vcf", files_axel)]
samples_axel = gsub(".*/(.*)_.*_Q100.*", "\\1", vcfs_axel)

#Get callable loci files from axels data and ensure its sorted in the same way as the vcfs.
callable_files_axel_all = list.files(callable_axel_dir, full.names = T)
callable_files_axel = callable_files_axel_all[grepl("_merged_autosomal.bed", callable_files_axel_all)]
samples_callable_axel = gsub(".*/(.*)_.*_CallableLoci.*", "\\1", callable_files_axel)
donor_callable_axel = gsub(".*/.*_(.*)MSC.*_CallableLoci.*", "\\1", callable_files_axel)
callable_tb = tibble(fetus = donor_callable_axel, callableloci_merged = callable_files_axel, callableloci_merged_true = TRUE, experiment = "Adult", samples_callable_axel) #Create tibble, that can be used as the factors for Do_mutational patterns.
factors_adult = enframe(samples_axel) %>% left_join(callable_tb, by = c("value" = "samples_callable_axel")) %>% dplyr::select(-name, -value) #Sort the factors tibble based on the sample order from the vcfs


cb_samples = c("CB112MPPC17", "CB112MPPC9", "CB112MPPO6", "10817MPPA4")
cb_index = samples_axel %in% cb_samples
grl_axel = read_vcfs_as_granges(vcfs_axel[!cb_index], samples_axel[!cb_index], genome)
factors_adult = factors_adult[!cb_index,]

grl = c(grl_fetal, grl_axel)
factors = rbind(factors_fetal, factors_adult)
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)



exp_name = "fetalvsadult_si"
chromosomes = paste0("chr", 1:22)
overview_samples_si = overview_samples %>% dplyr::filter(celltype == "SI" & trisomy == "FALSE")
gr_list = apply(overview_samples_si, 1, function(x) get_all_snvs(x, chromosomes))
fetuses_v = unique(overview_samples_si$fetus)
grl_fetal = vector("list", length(fetuses_v))

for (i in 1:length(fetuses_v)){#Pool per fetus. Only count shared muts once.
    fetus = fetuses_v[i]
    j = grep(fetus, overview_samples_si$fetus)
    gr_fetus = do.call("c", gr_list[j]) %>% unique()
    grl_fetal[[i]] = gr_fetus
}

grl_fetal = GRangesList(grl_fetal)
names(grl_fetal) = fetuses_v

overview_fetuses = overview_samples_si %>% dplyr::filter(!duplicated(fetus))
factors_fetal = tibble("fetus" = overview_fetuses$fetus, "callableloci_merged" = paste0(out_dir_base, overview_fetuses$fetus, "/callableloci/", overview_fetuses$bulk, "CallableLoci_autosomal.bed"), "callableloci_merged_true" = TRUE, "experiment" = "Fetal")

#Get vcf files from blokzijl data
vcfs_blokzijl =list.files(vcfs_blokzijl_dir, full.names = T)
samples_blokzijl = gsub(".*/(.*)_.*_Q100.*", "\\1", vcfs_blokzijl)
grl_blokzijl = read_vcfs_as_granges(vcfs_blokzijl, samples_blokzijl, genome)

donor_blokzijl = gsub(".*/.*_(.*)_Q100.*", "\\1", vcfs_blokzijl)
factors_adult = tibble("fetus" = donor_blokzijl, "callableloci_merged" = NA, "callableloci_merged_true" = F, "experiment" = "Post-infant")

grl = c(grl_fetal, grl_blokzijl)
factors = rbind(factors_fetal, factors_adult)
do_mutationalpatterns(grl, exp_name, factors, out_dir_mut, chromosomes)
do_better_refitting(grl, exp_name, factors, out_dir_mut, chromosomes)
#do_distribution(grl, exp_name, factors, out_dir_mut, chromosomes)





###Subclonal muts
chromosomes = paste0("chr", c(1:22, "X"))
exist_f = file.exists(overview_samples$subclonal_vcf)
overview_subclonal = overview_samples[exist_f,]
grl = read_vcfs_as_granges(overview_subclonal$subclonal_vcf, overview_subclonal$sample, ref_genome)
factors = overview_subclonal %>% dplyr::select(fetus, callableloci_merged, callableloci_merged_true, experiment = trisomy)

#subclone_mut_mat = mut_matrix(grl, ref_genome)
#subclone_mut_mat = pool_mut_mat(subclone_mut_mat, overview_subclonal$trisomy)
do_better_refitting(grl, "subclonal_inclewart", factors, out_dir_mut, chromosomes)


#low vs high vaf
chromosomes = paste0("chr", c(1:22, "X"))
overview_samples %<>% dplyr::mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), to = c("T21", "D21")))

mut_mat_vaf = get_vafsplit(overview_samples$unique_vcf, samples = overview_samples$sample, chromosomes)
col_names = colnames(mut_mat_vaf)
highmut_mut_mat_vaf_i = grep("OS1LIHSPCC1D23|OS1LIHSPCC2B21|OS1LIHSPCC2D6", col_names)
highmut_mut_mat_vaf = mut_mat_vaf[,highmut_mut_mat_vaf_i]

group = paste0(rep(overview_samples$trisomy, each = 2), "_", rep(overview_samples$celltype, each = 2), c("_lowvaf", "_highvaf"))
mut_mat_vaf_pooled = pool_mut_mat(mut_mat_vaf[,-highmut_mut_mat_vaf_i], group[-highmut_mut_mat_vaf_i])
combi_mut_mat_vaf = cbind(highmut_mut_mat_vaf, mut_mat_vaf_pooled)

cos_heat_fig = plot_inner_cosheat(combi_mut_mat_vaf)

#add subclonal
exist_f = file.exists(overview_samples$subclonal_vcf)
overview_subclonal = overview_samples[exist_f,]
grl = read_vcfs_as_granges(overview_subclonal$subclonal_vcf, overview_subclonal$sample, ref_genome)
mut_mat_sub = mut_matrix(grl, ref_genome)

#split in highly mutated and normal.
high_mut_i = colnames(mut_mat_sub) %in% c("OS1LIHSPCC1D23", "OS1LIHSPCC2B21", "OS1LIHSPCC2D6")
highmut_mut_mat_sub = mut_mat_sub[,high_mut_i, drop = F]
colnames(highmut_mut_mat_sub) = paste0(colnames(highmut_mut_mat_sub), "_subclonal")
group = paste0(overview_subclonal$trisomy, "_", overview_subclonal$celltype, "_subclonal")
mut_mat_sub_pooled = pool_mut_mat(mut_mat_sub[, -high_mut_i, drop = F], group[-high_mut_i])

combi_mut_mat_sub = cbind(highmut_mut_mat_sub, mut_mat_sub_pooled)
combi_mut_mat = cbind(combi_mut_mat_vaf, combi_mut_mat_sub)

cos_heat_fig2 = plot_inner_cosheat(combi_mut_mat)

pdf(paste0(out_dir_mut, "lowvshighvaf_heatfig.pdf"))
print(cos_heat_fig)
print(cos_heat_fig2)
dev.off()

# #Make figure 3a of manuscript
# divstri_types = read_tsv("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/di_vs_tri_inclewart_noOS1B21_type_occur_ungrouped.txt")
# fetalvsadult_types = read_tsv("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/fetalvsadult_type_occur_ungrouped.txt")
# types = rbind(divstri_types, fetalvsadult_types) %>% dplyr::select(-experiment)
# group = factor(c("T21", rep("D21", 3), "T21", "T21", "D21", "D21", "T21", rep("Fetal HSPC", 3), rep("Post-infant HSPC", 18)), levels = c("Post-infant HSPC", "Fetal HSPC", "D21", "T21"))
# fig = plot_spectrum(types, CT = T, by = group) + theme_classic()
# pdf("~/hpc/pmc_vanboxtel/projects/Freek_trees/mutpatterns/fig3a.pdf")
# fig
# dev.off()

#Make figure 5b of revision
fetal_adult_mat = read.table("~/surfdrive/Shared/Projects/Freek/Freek_trees/mutpatterns/fetalvsadult_combined_mut_mat.txt", 
                             sep = "\t",
                             stringsAsFactors = F)
di_tri_mat = read.table("~/surfdrive/Shared/Projects/Freek/Freek_trees/mutpatterns/not_bulk_uniq_blood_di_vs_tri_noOS1B21_combined_mut_mat.txt",
                        sep = "\t",
                        stringsAsFactors = F)
mut_mat = cbind(fetal_adult_mat$Adult, di_tri_mat)
colnames(mut_mat) = c("post_infant", "d21_hspc", "t21_hspc")



#
####__________________Make pca of all the data_____________####
#Create directory for mutational patterns analyses
out_dir_mut = paste0(out_dir_base, "mutpatterns/")
if (!dir.exists(out_dir_mut)){
    dir.create(out_dir_mut)
}


chromosomes = paste0("chr", 1:22)
gr_list = apply(overview_samples, 1, function(x) get_all_snvs(x, chromosomes))
grl_fetal = GRangesList(gr_list)
names(grl_fetal) = overview_samples$sample
factors_fetal = overview_samples %>% dplyr::select(donor = fetus, trisomy, celltype) %>% mutate(age = "Fetal")

#Get vcf files from axels data
files_axel =list.files(vcfs_axel_dir, full.names = T)
vcfs_axel = files_axel[grepl("_MQ60.vcf", files_axel)]
samples_axel = gsub(".*/(.*)_.*_Q100.*", "\\1", vcfs_axel)
cb_samples = c("CB112MPPC17", "CB112MPPC9", "CB112MPPO6", "10817MPPA4")
cb_index = samples_axel %in% cb_samples

grl_axel = read_vcfs_as_granges(vcfs_axel, samples_axel, genome)
grl = c(grl_fetal, grl_axel)
factors_adult = tibble(donor = samples_axel, trisomy = "FALSE", celltype = "HSC") %>% mutate(age = ifelse(cb_index, "Cord_blood",  "Adult"))
factors = rbind(factors_fetal, factors_adult) %>% mutate(trisomy = mapvalues(as.factor(trisomy), from = c("21", "FALSE"), to = c("Trisomy 21", "Disomy")))
factors %<>% mutate(group = paste(age, celltype, trisomy))

mut_mat = mut_matrix(grl, ref_genome)
pca_res = mut_mat %>% scale() %>% t() %>% prcomp(scale = T, center = T)

prcom_var = pca_res$sdev * pca_res$sdev
total_variance = sum(prcom_var)
nr_sig_pcs = bsDimension(prcom_var)
nr_pcs = pca_res$sdev %>% length()
expected_var_pc1 = brokenStick(1, nr_pcs)
expected_var_pc2 = brokenStick(2, nr_pcs)

pca_axis_fig = ggbiplot(pca_res, obs.scale = 1, ellipse = T, var.axes = T,
                        groups = as.vector(factors$group), choices = c(1,2) ,varname.size = 1, point_size = 2.5) +
    theme_classic() +
    annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
pca_fig = ggbiplot(pca_res, obs.scale = 1, ellipse = T, var.axes = F,
                        groups = as.vector(factors$group), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
    theme_classic() +
    annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))
pca_donor_fig = ggbiplot(pca_res, obs.scale = 1, ellipse = T, var.axes = F,
                   groups = as.vector(factors$donor), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
    theme_classic() +
    annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))

#Add signatures
signatures = read.table(cosmic_sig_fname, sep = "\t", header = T)
signatures = as.matrix(signatures[,-c(1,2)])
rownames(signatures) = rownames(mut_mat)

#Sig 1
sigs_pca_input = signatures[, "SBS1", drop = F] %>% scale() %>% t()
predicted = predict(pca_res, newdata = sigs_pca_input)
prcom_withsig1 = pca_res
prcom_withsig1$x = rbind(prcom_withsig1$x, predicted)
pca_sig1_fig = ggbiplot(prcom_withsig1, obs.scale = 1, ellipse = T, var.axes = F,
                        groups = c(as.vector(factors$group), rownames(sigs_pca_input)), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
    theme_classic() +
    annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))


#Overlay all signatures to existing pca model
sigs_pca_input = signatures %>% scale() %>% t()
predicted = predict(pca_res, newdata = sigs_pca_input)
prcom_withsigs = pca_res
prcom_withsigs$x = rbind(prcom_withsigs$x, predicted)
pca_sigs_fig = ggbiplot(prcom_withsigs, obs.scale = 1, ellipse = T, var.axes = F,
                        groups = c(as.vector(factors$group), rownames(sigs_pca_input)), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
    theme_classic() +
    annotate("text", y = 4, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant"))


#cosine similarities between pcs and signatures
pc_sig_cosim = cos_sim_matrix(pca_res$rotation, signatures) %>% abs()
pc_sig_fig = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F)

pc1_positive = pca_res$rotation[,1] > 0
pc_sig_cosim = cos_sim_matrix(pca_res$rotation[pc1_positive,], signatures[pc1_positive,]) %>% abs()
pc_sig_fig_pos = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
    labs(title = "Positive PC1")
pc_sig_cosim = cos_sim_matrix(pca_res$rotation[!pc1_positive,], signatures[!pc1_positive,]) %>% abs()
pc_sig_fig_neg = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
    labs(title = "Negative PC1")

pc2_positive = pca_res$rotation[,2] > 0
pc_sig_cosim = cos_sim_matrix(pca_res$rotation[pc2_positive,], signatures[pc2_positive,]) %>% abs()
pc_sig2_fig_pos = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
    labs(title = "Positive PC2")
pc_sig_cosim = cos_sim_matrix(pca_res$rotation[!pc2_positive,], signatures[!pc2_positive,]) %>% abs()
pc_sig2_fig_neg = plot_cosine_heatmap(pc_sig_cosim, cluster_rows = F) +
    labs(title = "Negative PC2")




#Sum up all fetal stuf, and see if the increased number of mutations (and thus reduced noise) causes it to be placed somewhere else
mut_mat2 = mut_mat
fetal_summed = mut_mat[,factors$age == "Fetal", drop = F] %>% rowSums()
mut_mat2 = cbind(mut_mat2, fetal_summed)
pca_res2 = mut_mat2 %>% scale() %>% t() %>% prcomp(scale = T, center = T)
pca_fetal_summed_fig = ggbiplot(pca_res2, obs.scale = 1, ellipse = T, var.axes = F,
                   groups = c(as.vector(factors$group), "fetal_summed"), choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
    theme_classic()

#6 way
type_occur = mut_type_occurrences(grl, ref_genome)
prcom_6 = type_occur %>% t() %>% scale() %>% t() %>% prcomp(scale = T, center = T)
prcom_var = prcom_6$sdev * prcom_6$sdev
nr_sig_pcs = bsDimension(prcom_var)
nr_pcs = prcom_6$sdev %>% length()
expected_var_pc1 = brokenStick(1, nr_pcs)
expected_var_pc2 = brokenStick(2, nr_pcs)
pca_sixway_fig = ggbiplot(prcom_6, obs.scale = 1, var.axes = F, ellipse = T, 
                          groups = factors$group, choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
    theme_classic() +
    annotate("text", y = 3, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant. (Based on 7 base substitutions.)"))

pca_sixway_axis_fig = ggbiplot(prcom_6, obs.scale = 1, var.axes = T, ellipse = T, 
                               groups = factors$group, choices = c(1,2) ,varname.size = 2.5, point_size = 2.5) +
    theme_classic() +
    annotate("text", y = 3, x = 0, label = paste0("First ", nr_sig_pcs, " PCs are significant. (Based on 7 base substitutions.)"))

ha = HeatmapAnnotation(as.data.frame(factors))
mut_mat_fract = t(t(mut_mat) / colSums(mut_mat))
heat_fig = Heatmap(top_annotation = ha, mut_mat_fract, clustering_method_columns = "ward.D2", clustering_distance_columns = "euclidean",
                   row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 7),
                   heatmap_legend_param = list(title = "Fraction"),
                   col = colorRamp2(seq(0, 0.20, 0.025), brewer.pal(9, "Blues")))
#heat_fig = mut_mat %>% scale() %>% Heatmap(top_annotation = ha)


#t-SNE
mut_mat_scaled = mut_mat %>% scale() %>% t()

KL_div = 1e6
for (i in 1:10){ #Repeat t-SNE ten times. Pick the one, that has the lowest KL divergence.
    tsne_loop = Rtsne(mut_mat_scaled, perplexity = 5, max_iter = 5000, pca_center = T, pca_scale = F, theta = 0.5)
    KL_div_loop = tsne$itercosts[length(tsne$itercosts)]
    if (KL_div_loop < KL_div){
        tsne = tsne_loop
        KL_div = KL_div_loop
    }
}

d_tsne = as.data.frame(tsne$Y) %>% dplyr::mutate("Group" = factors$group)
tsne_fig = ggplot(d_tsne, aes(x = V1, y = V2, colour = Group)) +
    geom_point(size = 1) +
    labs(x = "", y = "", title = "t-SNE") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    )


pdf(paste0(out_dir_mut, "pca_everything.pdf"), height = 8.5)
pca_axis_fig
pca_fig
pca_donor_fig
pca_sig1_fig
pca_sigs_fig
pc_sig_fig
pc_sig_fig_pos
pc_sig_fig_neg
pc_sig2_fig_pos
pc_sig2_fig_neg
pca_fetal_summed_fig
pca_sixway_fig
pca_sixway_axis_fig
heat_fig
tsne_fig
dev.off()


#Function that returns a granges object for all the mutations that are in bulk for the specified fetuses.
get_inbulk_gr = function(fetuses_interest){
    gr_list = vector("list", length(fetuses_interest))
    for (i in seq(1, length(fetuses_interest))){
        fetus_name = fetuses_interest[i]
        overview_fetus = overview_samples %>% dplyr::filter(fetus == fetus_name)
        bulk = overview_fetus$bulk[1]
        unique_vcf_fnames = overview_fetus$uniq_vcf_inbulk[overview_fetus$uniq_vcf_inbulk_true]
        
        #Check if any muts in the bulk exist for this fetus
        if (length(unique_vcf_fnames) == 0 & !overview_fetus$shared_vcf_inbulk_true[1]){
            print(paste0("No mutations were found for fetus ", fetus_name))
            gr_list[[i]] = GRanges()
            next
        }
        
        #Get mutations that are unique to a clone
        if (length(unique_vcf_fnames > 0)){
            vcfs_uniq = lapply(unique_vcf_fnames, function(x) readVcf(x, genome = "hg19"))
            vcf_uniq = do.call(rbind, vcfs_uniq)
            vcf_uniq_exists = T
        } else{
            vcf_uniq_exists = F
        }
        
        #Combine with the shared mutations
        if (overview_fetus$shared_vcf_inbulk_true[1]){
            shared_vcf = readVcf(overview_fetus$shared_inbulk_vcf[1], genome = "hg19")
            if (vcf_uniq_exists == T){
                vcf = rbind(shared_vcf, vcf_uniq)
            } else{
                vcf = shared_vcf
            }
        } else{
            vcf = vcf_uniq
        }
        
        #Transform vcf into granges object and store in list.
        gr = granges(vcf)
        gr_list[[i]] = gr
    }
    gr = do.call(c, gr_list) %>% sort()
    return(gr)
}

gr = get_inbulk_gr(c("N01", "OS1", "MH3", "E080416"))
grl = GRangesList(gr)
seqlevelsStyle(grl) = "UCSC"
mut_type = mut_type_occurrences(grl, ref_genome = genome)
spectrum_early_fig = plot_spectrum(mut_type, CT = F)
pdf(paste0(out_dir_mut, "inbulk_trisomy_mutprofile.pdf"))
spectrum_early_fig
dev.off()

shared_vcf_fnames = overview_samples %>% dplyr::filter((fetus %in% c("N01", "OS1", "MH3", "E080416") & shared_vcf_true == "TRUE")) %>% pull(shared_vcf) %>% unique()
gr_list = lapply(shared_vcf_fnames, function(x) {
    vcf = readVcf(x, genome = "hg19")
    gr = granges(vcf)
    print(paste0("Read in vcf: ", x))
    return(gr)
})
grl = GRangesList(gr_list)
seqlevelsStyle(grl) = "UCSC"
mut_type = mut_type_occurrences(grl, ref_genome = genome)
spectrum_early_fig = plot_spectrum(mut_type, CT = F)
pdf(paste0(out_dir_mut, "shared_trisomy_mutprofile.pdf"))
spectrum_early_fig
dev.off()


ff
####_______________Calculate whether early mutations occur in the same cell division. Also calculate median vaf and create contribution pie charts._____####

#Colour function to use for plotting vafs
contri_cols = RColorBrewer::brewer.pal(9, "RdYlBu") %>% colorRampPalette()
contri_col_tb = tibble("value" = seq(0,1, by = 0.01), "color" = rev(contri_cols(101)))

overview_samples = create_overview_samples(fetuses, out_dir_base)
for (fetus in unique(overview_samples$fetus)){
    overview_fetus = overview_samples %>% dplyr::filter(fetus == UQ(fetus))
    bulk = overview_fetus$bulk[1]
    other_bulks = overview_fetus$other_bulks[1]
    uniq = overview_fetus$uniq_vcf_inbulk[overview_fetus$uniq_vcf_inbulk_true]
    shared = overview_fetus$shared_inbulk_vcf[overview_fetus$shared_vcf_inbulk_true]
    indel_uniq = overview_fetus$indel_uniq_inbulk[overview_fetus$indel_uniq_inbulk_true]
    indel_shared = overview_fetus$indel_shared_inbulk[overview_fetus$indel_shared_inbulk_true]
    inbulk_vcfs_fnames = c(uniq, shared, indel_uniq, indel_shared) %>% unique()
    if (length(inbulk_vcfs_fnames) == 0){
        next
    }
    vcfs = lapply(inbulk_vcfs_fnames, function(x) readVcf(x, genome = "hg19"))
    vcf = do.call(rbind, vcfs)
    if (nrow(vcf) == 0){
        next
    }
    
    all_bulks = get_allbulks(other_bulks, bulk)
    bulkcols = grep(all_bulks, samples(header(vcf)))
    gt = geno(vcf)$GT[,-bulkcols, drop = F]
    gt[gt == "1/1"] = "0/1" #Makes sure homozygous mutations on the x chromosome are placed in the same combi as muts on the autosomes.
    samples = colnames(gt)
    
    
    
    #Identify combinations
    combis = apply(gt, 1, paste, collapse = "\t")
    combi_order = order(combis)
    combi_ids = factor(combis)
    levels(combi_ids) = seq(1, nlevels(combi_ids))
    unique_combis = unique(combis)
    
    
    for (bulk_selected in strsplit(all_bulks, "\\|")[[1]]){
        bulkcol = grep(bulk_selected, samples(header(vcf)))
        
        #Create table with the vaf of all inbulk mutations
        vaf = get_vaf(vcf, bulk_selected) %>% unname()
        gt_tb = gt %>% as_tibble() %>% mutate(vaf = vaf, id = combi_ids, rounded_vaf = round(vaf, digits = 2))
        gt_tb = left_join(gt_tb, contri_col_tb, by = c("rounded_vaf" = "value")) %>% dplyr::select(-rounded_vaf, contribution_col = color)
        gt_tb = gt_tb[combi_order,]
        write_tsv(gt_tb, paste0(out_dir_base, fetus, "/", bulk_selected,"_bulk_vafs.txt"))
        
        #Perform calculations for muts that occured in the same combi of clones
        rows_contrifigs = lapply(unique_combis, function(combi){
            index = grep(combi, combis)
            vcf_combi = vcf[index,]
            nr_muts = length(index)
            id = combi_ids[index][[1]]
            
            #Calculate whether muts occured in the same cell division
            if(nr_muts > 1){
                ad = geno(vcf_combi)$AD[,bulkcol]
                ad = do.call(rbind, ad)
                chi_res = chisq.test(ad, simulate.p.value = T)
                pval = chi_res$p.value
            } else{
                pval = NA
            }
            vaf = get_vaf(vcf_combi, bulk_selected)
            vaf_med = median(vaf)
            combi_split = strsplit(combi, "\t")
            row = c(unlist(combi_split), pval, vaf_med, nr_muts, id)
            
            #Create contribution figure
            contri = vaf_med * 2
            vaf_tb = tibble("contribution" = c(contri, (1-contri)), "group" = c("vaf", "rest"))
            contri_fig = contri_pie(vaf_tb, paste0(combi, "\tid: ", id))
            
            contri_round = round(contri, digits = 2)
            my_contri_col = contri_col_tb %>% dplyr::filter(dplyr::near(value, contri_round)) %>% pull(color)
            row = c(row, my_contri_col)
            
            output = list(row, contri_fig)
            return(output)
        })
        
        rows = lapply(rows_contrifigs, function(x) x[[1]])
        contrifigs = lapply(rows_contrifigs, function(x) x[[2]])
        gt_pvals = do.call(rbind, rows)
        colnames(gt_pvals) = c(samples, "Pval", "Median_vaf", "nr_muts", "id", "vaf_col")
        write_tsv(as_tibble(gt_pvals), paste0(out_dir_base, fetus, "/", bulk_selected,"_bulk_samedivision.txt"))
        
        pdf(paste0(out_dir_base, fetus, "/", bulk_selected,"_contri_figs.pdf"))
        print(contrifigs)
        dev.off()
        
        #Calculate whether different pairs of combinations. (So different parts of the tree) have a significantly different bulk)
        if (length(unique_combis) < 2){
            next
        }
        combi_pairs = combn(unique_combis, 2, simplify = F)
        rows_combi_pairs = lapply(combi_pairs, function(combi_pair){
            index1 = grep(combi_pair[[1]], combis)
            vcf_combi1 = vcf[index1,]
            id1 = combi_ids[index1][[1]]
            ad1 = geno(vcf_combi1)$AD[,bulkcol]
            ad1_sum = do.call(rbind, ad1) %>% colSums()
            
            index2 = grep(combi_pair[[2]], combis)
            vcf_combi2 = vcf[index2,]
            id2 = combi_ids[index2][[1]]
            ad2 = geno(vcf_combi2)$AD[,bulkcol]
            ad2_sum = do.call(rbind, ad2) %>% colSums()
            
            ad_both = rbind(ad1_sum, ad2_sum)
            fisher_res = fisher.test(ad_both)
            pval = fisher_res$p.value
            row = c(id1, id2, pval)
        })
        combi_pairs_tb = do.call(rbind, rows_combi_pairs) %>% as_tibble()
        colnames(combi_pairs_tb) = c("id1", "id2", "pval")
        write_tsv(combi_pairs_tb, paste0(out_dir_base, fetus, "/", bulk_selected,"_combi_pairs.txt"))
    }
    
    #Compare vafs of combinations between the different bulks
    all_bulks_v = strsplit(all_bulks, "\\|")[[1]]
    if (length(all_bulks_v) < 2){
        next
    }
    bulk_pairs = combn(all_bulks_v, 2, simplify = F)
    for (bulk_pair in bulk_pairs){
        rows_combi_pairs = lapply(unique_combis, function(combi){
            index = grep(combi, combis)
            vcf_combi = vcf[index,]
            id = combi_ids[index][[1]]
            ad1 = geno(vcf_combi)$AD[,bulk_pair[1]]
            ad1_sum = do.call(rbind, ad1) %>% colSums()
            
            ad2 = geno(vcf_combi)$AD[,bulk_pair[2]]
            ad2_sum = do.call(rbind, ad2) %>% colSums()
            
            ad_both = rbind(ad1_sum, ad2_sum)
            fisher_res = fisher.test(ad_both)
            pval = fisher_res$p.value
            row = c(id, pval)
        })
        combi_pairs_tb = do.call(rbind, rows_combi_pairs) %>% as_tibble()
        colnames(combi_pairs_tb) = c("id", "pval")
        write_tsv(combi_pairs_tb, paste0(out_dir_base, fetus, "/", bulk_pair[1], "_vs_", bulk_pair[2],"_combi_pairs.txt"))
    }
}

col_legend_fig = ggplot(contri_col_tb, aes(x = value, y = 1, color = value)) +
    geom_line() +
    scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(9, "RdYlBu")), limits = c(0,1))

pdf(paste0(out_dir_base, "tree_vafcolor_legend.pdf"))
print(col_legend_fig)
dev.off()

empty_tb = tibble("contribution" = c(0,1), "group" = c("vaf", "rest"))
empty_pie_fig = contri_pie(empty_tb, "Empty")
full_tb = tibble("contribution" = c(1,0), "group" = c("vaf", "rest"))
full_pie_fig = contri_pie(full_tb, "Full")

pdf(paste0(out_dir_base, "pie_full_empty.pdf"))
print(empty_pie_fig)
print(full_pie_fig)
dev.off()



####__________________Create vaf vs depth plot___________####
for (fetus in unique(overview_samples$fetus)){
    overview_fetus = overview_samples %>% dplyr::filter(fetus == UQ(fetus))
    bulk = overview_fetus$bulk[1]
    other_bulks = overview_fetus$other_bulks[1]
    shared = overview_fetus$shared_vcf[overview_fetus$shared_vcf_true] %>% unique()
    shared_indel = overview_fetus$indel_shared[overview_fetus$indel_shared_true] %>% unique()
    vcf_fnames = c(shared, shared_indel)

    if (length(vcf_fnames) == 0){
        next
    }
    vcfs = lapply(vcf_fnames, function(x) readVcf(x, genome = "hg19"))
    vcf = do.call(rbind, vcfs)
    if (nrow(vcf) == 0){
        next
    }
    all_bulks = get_allbulks(other_bulks, bulk)
    bulkcols = grep(all_bulks, samples(header(vcf)))
    gt = geno(vcf)$GT[,-bulkcols, drop = F]
    gt[gt == "1/1"] = "0/1" #Makes sure homozygous mutations on the x chromosome are placed in the same combi as muts on the autosomes.
    samples = colnames(gt)
    
    #Identify combinations
    combis = apply(gt, 1, paste, collapse = "\t")
    combi_order = order(combis)
    combi_ids = factor(combis)
    levels(combi_ids) = seq(1, nlevels(combi_ids))
    unique_combis = unique(combis)
    bulk_v = strsplit(all_bulks, "\\|")[[1]]
    nrbulks = length(bulk_v)
    vaf_dp_tb = vector("list", nrbulks)
    for (i in 1:nrbulks){
        bulk_selected = bulk_v[i]
        bulkcol = grep(bulk_selected, samples(header(vcf)))
        
        #Create table with the vaf of all inbulk mutations
        vaf = get_vaf(vcf, bulk_selected)
        dp = geno(vcf)$DP[,bulkcol]
        mut_nr = seq(1, nrow(vcf))
        vaf_dp_tb_bulk = tibble(vaf, dp, bulk_selected, combis)
        vaf_dp_tb_bulk %<>% group_by(combis) %>%dplyr::mutate(nr_mut = dplyr::row_number())
        vaf_dp_tb[[i]] = vaf_dp_tb_bulk
    }
    vaf_dp_tb = do.call(rbind, vaf_dp_tb)
    vaf_dp_tb %<>% mutate(nr_mut = as.factor(nr_mut))
    
    vaf_tb_fig = ggplot(vaf_dp_tb, aes(x = vaf, y = dp, color = bulk_selected, shape = nr_mut)) +
        geom_point() +
        facet_grid(. ~ combis) +
        scale_shape_manual(values=1:nlevels(vaf_dp_tb$nr_mut)) +
        theme_classic() +
        theme(text = element_text(size = 24, family = "Arial"), strip.text.x = element_text(size = 12))

    pdf(paste0(out_dir_base, fetus, "/shared_dpvsvaf.pdf"), width = 14)
    print(vaf_tb_fig)
    dev.off()
    
}




####___________________Look for somatic drivers________________####
chromosomes = paste0("chr", c(1:22, "X"))
cancer_genes = read_tsv("~/hpc/pmc_vanboxtel/data/Cancer_genes/Cosmic_cancer_gene_census_09052019.txt")
cancer_genes_somatic = cancer_genes %>% dplyr::filter(Somatic == "yes")

trash_can = apply(overview_samples, 1, function(x) create_possible_driver_vcfs(x, chromosomes, out_dir_base, cancer_genes_somatic$`Gene Symbol`))






# Telomeres ---------------------------------------------------------------

tel_fnames = list.files(telomere_dir_fname, full.names = T) %>% grep(".csv$", ., value = T)
tel_tb = purrr::map(tel_fnames, read_csv, col_types = c("cddddddddddd")) %>% 
    do.call("rbind", .) %>% 
    dplyr::mutate(Sample = gsub("_dedup.realigned.bam|_combined.bam", "", Sample))
write_tsv(tel_tb, paste0(out_dir_base, "telomeres/telomere_lengths.txt"))

ggplot(tel_tb, aes(x = "", y = Length)) +
    geom_quasirandom() +
    theme_classic()

overview_samples = overview_samples %>% dplyr::filter(celltype != "Liver")
tel_tb_clean = tel_tb %>% dplyr::filter(!grepl("skin|bulk", Sample, ignore.case = T) & Length != 0)
tel_tb_clean = inner_join(tel_tb_clean, overview_samples, by = c("Sample" = "sample")) %>% dplyr::mutate(trisomy = mapvalues(trisomy, from = c("21", "FALSE"), to = c("T21", "D21")))

tel_dot_fig = ggplot(tel_tb_clean, aes(x = "", y = Length, colour = trisomy, shape = fetus)) +
    geom_quasirandom(size = 2) +
    theme_classic()

lme_tel = lme(Length ~ age_year + trisomy * celltype, random = ~ -1 + age_year | fetus, data = tel_tb_clean)
lme_tel_tb = get_fixed(lme_tel)
write_tsv(lme_tel_tb, paste0(out_dir_base, "telomeres/lme_vars.txt"))

tel_tb_clean$predict = predict(lme_tel, level = 0)
pval = lme_tel_tb %>% dplyr::filter(variable == "trisomyT21") %>% pull("p-value")
tel_line_fig = ggplot(tel_tb_clean, aes(x = age_year, y = Length)) +
    geom_point(aes(color = trisomy, shape = celltype), size = 4) +
    geom_line(aes(y = predict, color = trisomy, linetype = celltype), size = 1.5) +
    #facet_wrap(. ~ age_cat, scales = "free") +
    labs(title = "LMM containing all data", x = "Age (years)", y = "Telomere length") +
    theme_classic() +
    theme(text = element_text(size = 24, family = "Arial")) +
    annotate("text", x = 0.2, y = 15000, label = paste0("P (D21 vs T21): ", round(pval, 3)), size = 6) +
    coord_cartesian(xlim = c(0.15, 0.3))

pdf(paste0(out_dir_base, "telomeres/lme_divstri.pdf"))
tel_dot_fig
tel_line_fig
dev.off()


# Make table of exonic genes ----------------------------------------------

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
vcf_fnames = str_c("~/hpc/pmc_vanboxtel/projects/Freek_trees/", fetuses, "/", fetuses, "_complete.vcf")
vcf_l = purrr::map(vcf_fnames, readVcf)
exon_tb = purrr::map(vcf_l, get_exon_table) %>% 
    set_names(fetuses) %>% 
    bind_rows(.id = fetus)
write_tsv(exon_tb, str_c(out_dir_base, "exon_muts.txt"))
