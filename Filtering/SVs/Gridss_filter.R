library(tidyverse)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(optparse)
library(ggbio)
library(Homo.sapiens)



sort_vcf = function(vcf){
    
    if (class(vcf)[1] != "CollapsedVCF"){
        warning(paste0("This function was hg19_gred on the Collapsed VCF class. It may not work on the supplied input vcf, which is a ", class(vcf[1]), " object."))
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

#Add genes that are near the breakpoints.
add_nearby_genes = function(vcf){
    gr = granges(vcf)
    seqlevelsStyle(gr) = "UCSC"
    
    #Get the genes
    genes_gr = suppressMessages(genes(Homo.sapiens, columns = c("SYMBOL")))
    genes_gr = keepSeqlevels(genes_gr, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")
    
    #Overlap the SVs and the genes.
    overlap_index = findOverlaps(gr, genes_gr) %>% as.data.frame()
    gene_symbols = genes_gr[overlap_index$subjectHits]$SYMBOL %>% unlist()
    
    #Add the gene name to the vcf
    header_add = DataFrame("Number" = 1, "Type" = "Character", "Description" = "Gene located at breakpoint", row.names = "Gene")
    info(header(vcf)) = rbind(info(header(vcf)), header_add)
    info(vcf)$Gene = "NA"
    info(vcf)$Gene[overlap_index$queryHits] = gene_symbols
    
    #gr$Gene = NA
    #gr[overlap_index$queryHits]$Gene = gene_symbols
    
    #Find genes near the SVs.
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        gr_ext = gr + 5000
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    gr_ext = gr_ext %>% trim()
    overlap_index = findOverlaps(gr_ext, genes_gr) %>% as.data.frame()
    overlap_index$gene_symbols = genes_gr[overlap_index$subjectHits]$SYMBOL %>% unlist()
    near_genes = overlap_index %>% dplyr::group_by(queryHits) %>% dplyr::summarise(near_genes = paste(gene_symbols, collapse = ";"), n = n())
    
    header_add = DataFrame("Number" = 1, "Type" = "Character", "Description" = "Genes located within 5kb", row.names = "Near_gene_5kb")
    info(header(vcf)) = rbind(info(header(vcf)), header_add)
    header_add = DataFrame("Number" = 1, "Type" = "Character", "Description" = "Nr of genes located within 5kb", row.names = "Nr_near_gene_5kb")
    info(header(vcf)) = rbind(info(header(vcf)), header_add)
    info(vcf)$Near_gene_5kb = "NA"
    info(vcf)$Nr_near_gene_5kb = 0
    info(vcf)$Near_gene_5kb[near_genes$queryHits] = near_genes$near_genes
    info(vcf)$Nr_near_gene_5kb[near_genes$queryHits] = near_genes$n
    # gr$near_genes_5kb = NA
    # gr$n_near_genes_5kb = 0
    # gr[near_genes$queryHits]$near_genes_5kb = near_genes$near_genes
    # gr[near_genes$queryHits]$n_near_genes_5kb = near_genes$n
    
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        gr_ext = gr + 100000
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    gr_ext = gr_ext %>% trim()
    overlap_index = findOverlaps(gr_ext, genes_gr) %>% as.data.frame()
    overlap_index$gene_symbols = genes_gr[overlap_index$subjectHits]$SYMBOL %>% unlist()
    near_genes = overlap_index %>% dplyr::group_by(queryHits) %>% dplyr::summarise(near_genes = paste(gene_symbols, collapse = ";"), n = n())
    
    header_add = DataFrame("Number" = 1, "Type" = "Character", "Description" = "Genes located within 100kb", row.names = "Near_gene_100kb")
    info(header(vcf)) = rbind(info(header(vcf)), header_add)
    header_add = DataFrame("Number" = 1, "Type" = "Character", "Description" = "Nr of genes located within 100kb", row.names = "Nr_near_gene_100kb")
    info(header(vcf)) = rbind(info(header(vcf)), header_add)
    info(vcf)$Near_gene_100kb = "NA"
    info(vcf)$Nr_near_gene_100kb = 0
    info(vcf)$Near_gene_100kb[near_genes$queryHits] = near_genes$near_genes
    info(vcf)$Nr_near_gene_100kb[near_genes$queryHits] = near_genes$n
    
    # gr$near_genes_100kb = NA
    # gr$n_near_genes_100kb = 0
    # gr[near_genes$queryHits]$near_genes_100kb = near_genes$near_genes
    # gr[near_genes$queryHits]$n_near_genes_100kb = near_genes$n
    return(vcf)
}

#Check whether the identified nearby genes are associated with cancer.
check_nearby_genes_cancer = function(vcf, cancer_genes_symbols){
    near_genes = info(vcf)$Near_gene_100kb
    near_genes_l = strsplit(near_genes, ";")
    nearby_cancer_genes = lapply(near_genes_l, function(x){
        index = x %in% cancer_genes_symbols
        nearby_cancer_genes = x[index]
        if (length(nearby_cancer_genes) == 0){
            return(NA)
        }
        nearby_cancer_genes = paste(nearby_cancer_genes, collapse = ";")
        return(nearby_cancer_genes)
    })
    
    header_add = DataFrame("Number" = 1, "Type" = "Character", "Description" = "Genes located within 100kb that are in the cosmic cancer census list and which are in the category somatic", row.names = "Cancer_gene_within_100kb")
    info(header(vcf)) = rbind(info(header(vcf)), header_add)
    info(vcf)$Cancer_gene_within_100kb = unlist(nearby_cancer_genes)

    return(vcf)
}

#Function copied from GRIDSS example.
simpleEventType <- function(gr) {
    return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
                  ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                         ifelse(strand(gr) == strand(partner(gr)), "INV",
                                ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
                                       "DUP")))))
}

#Function to annotate a vcf with the SV type, using the results from the simpleEvenType function from GRIDSS. This may not work well for complex rearrangements.
get_sv_type = function(vcf){
    bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)
    sv_type = simpleEventType(bpgr)
    bp_f = names(vcf) %in% names(bpgr)
    
    header_add = DataFrame("Number" = 1, "Type" = "Character", "Description" = "The SV sub-category", row.names = "SV_type")
    info(header(vcf)) = rbind(info(header(vcf)), header_add)
    
    info(vcf)$SV_type = "NA"
    info(vcf)$SV_type[bp_f] = sv_type
    return(vcf)
}

####_________________________________Parse arguments___________________________________####
option_list = list(
    make_option(c("-s", "--SAMPLE_NAME"), type="character", help="The name of the sample", default = NA),
    make_option(c("-d", "--DIR"), type="character", help="The working directory", default = NA),
    make_option(c("--pon_dir"), type = "character", help = "The directory containing the HMF PON files", default = NA),
    make_option(c("--vaf"), type = "double", help = "The vaf, that will be used for filtering.", default = NA),
    make_option(c("--ref_cell"), type = "character", help = "The name of the sample that will be used as the control (Name of bam before first underscore). If ALL is chosen, then the script will hg19_gr for all cells wether they lack the SV."),
    make_option(c("--cancer_genes"), type = "character", help = "A cosmic cancer gene census list", default = NA)

); 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


sample = opt$SAMPLE_NAME
wdir = opt$DIR
pon_dir = opt$pon_dir
vaf_lim = opt$vaf
ref_cell = opt$ref_cell
cancer_genes_fname = opt$cancer_genes

nr_missing_args = c(sample, wdir, pon_dir, vaf_lim, ref_cell, cancer_genes_fname) %>% is.na() %>% sum()
if(nr_missing_args > 0){
    stop("All arguments must be provided. See script usage (--help)")
}


####____________________________________________Filter vcf on PON file__________________________________________####
if (substring(wdir, nchar(wdir)) != "/"){
    wdir = paste0(wdir, "/")
}

vcf = readVcf(paste0(wdir, sample, ".sv.filtered.vcf"))



#Functions copied and slightly modified from hmf pipeline
gridss_overlaps_breakpoint_pon = function(gr, pon_dir=NULL, pongr=read_gridss_breakpoint_pon(paste(pon_dir, "gridss_pon_breakpoint.bedpe", sep="/")), ...) {
    hasHit = rep(FALSE, length(gr))
    if (!is.null(pongr)) {
        hasHit[findBreakpointOverlaps(gr, pongr[pongr$score >= 1], sizemargin=NULL, restrictMarginToSizeMultiple=NULL, ...)$queryHits] = TRUE
    }
    return(hasHit)
}
read_gridss_breakpoint_pon = function(file) {
    df = read_tsv(file,
                  col_names=c("chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2"),
                  col_types="ciiciicccc")
    gro = GRanges(
        seqnames=df$chr1,
        ranges=IRanges(
            start=df$start1 + 1,
            end=df$end1),
        strand=df$strand1,
        partner=paste0(seq_len(nrow(df)), "h"),
        score=df$score)
    names(gro) = paste0(seq_len(nrow(df)), "o")
    grh = GRanges(
        seqnames=df$chr2,
        ranges=IRanges(
            start=df$start2 + 1,
            end=df$end2),
        strand=df$strand2,
        partner=paste0(seq_len(nrow(df)), "o"),
        score=df$score)
    names(grh) = paste0(seq_len(nrow(df)), "h")
    return(c(gro, grh))
}

#Get SVs that have a breakpoint and determine whether they are in the hmf pon. This does not work for breakend sites.
full_bpgr = breakpointRanges(vcf, unpartneredBreakends=FALSE)
overlap_bpgr = gridss_overlaps_breakpoint_pon(full_bpgr, pon_dir)

bp_f = names(vcf) %in% names(full_bpgr)
vcf_pon = vcf[bp_f,][!overlap_bpgr,]

#For breakpoints check if the partner was kept. If not then remove the breakpoint.
bpgr = breakpointRanges(vcf_pon)
bp_f = names(vcf_pon) %in% names(bpgr)
vcf_pon = vcf_pon[bp_f,]

####____________________________________________Annotate breakpoints with nearby genes and check if they are in the cosmic gene census. Also add the SV type if possible_________####
cancer_genes = read_tsv(cancer_genes_fname)
cancer_genes_symbols = cancer_genes %>% dplyr::filter(Somatic == "yes") %>% pull(`Gene Symbol`)

vcf_pon = add_nearby_genes(vcf_pon)
vcf_pon = check_nearby_genes_cancer(vcf_pon, cancer_genes_symbols)
vcf_pon = get_sv_type(vcf_pon)


writeVcf(vcf_pon, paste0(wdir, sample, ".sv.filtered_PON.vcf"))

####________________________________________________________Determine whether SVs are somatic._________________________####
quals = geno(vcf_pon)$QUAL
if (ref_cell == "ALL"){
    absence = rowSums(quals == 0) > 0
} else{
    cells = gsub("_.*", "", colnames(quals))
    absence = quals[,cells %in% ref_cell] == 0
}

absence[is.na(absence)] = F
vcf_somatic = vcf_pon[absence,]

#For breakpoints check if the partner was kept. If not then remove the breakpoint.
bpgr = breakpointRanges(vcf_somatic)
bp_f = names(vcf_somatic) %in% names(bpgr)
vcf_somatic = vcf_somatic[bp_f,]


writeVcf(vcf_somatic, paste0(wdir, sample, ".sv.somatic.vcf"))

####____________________________________________________Determine vaf of SVs_____________________________________####
#Functions copied from hmf pipeline
.genosum <- function(genotypeField, columns) {
    rowSums(genotypeField[,columns, drop=FALSE])
}

is_short_deldup = function(gr) {
    is_deldup = rep(FALSE, length(gr))
    if (!is.null(gr$partner)) {
        isbp <- gr$partner %in% names(gr)
        bpgr <- gr[isbp]
        bp_short_deldup = strand(bpgr) != strand(partner(bpgr)) &
            seqnames(bpgr) == seqnames(partner(bpgr)) &
            abs(start(bpgr)-start(partner(bpgr))) < 1000
        is_deldup[isbp] = bp_short_deldup
    }
    return(is_deldup)
}

#For short deletions, you can get refpairs even on a mut allele. Therefore the refpairs should not be counted.
gridss_bp_af = function(gr, vcf, ordinal) {
    return(.gridss_af(gr, vcf, ordinal, !is_short_deldup(gr), includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE))
}
gridss_be_af = function(gr, vcf, ordinal) {
    return(.gridss_af(gr, vcf, ordinal, includeRefPair=rep(TRUE, length(gr)), includeBreakpointSupport=FALSE, includeBreakendSupport=TRUE))
}

.gridss_af = function(gr, vcf, ordinal, includeRefPair, no_coverage_af=0, includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE) {
    assertthat::are_equal(length(gr), length(includeRefPair))
    genotype = geno(vcf[names(gr)])
    g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], ordinal) } else { genotype[[field]] } })
    names(g) <- names(genotype)
    support = rep(0, length(gr))
    if (includeBreakpointSupport) {
        support = support + g$VF
    }
    if (includeBreakendSupport) {
        support = support + g$BVF
    }
    vf_af = support / (support + g$REF + ifelse(!includeRefPair, 0, g$REFPAIR))
    vf_af[is.nan(vf_af)] = no_coverage_af
    return(vf_af)
}

#Add vaf name to header
header_add = DataFrame("Number" = 1, "Type" = "Float", "Description" = "Variant allele frequency", row.names = "VAF")
geno(header(vcf_somatic)) = rbind(geno(header(vcf_somatic)), header_add)



#Add vafs to the vcfs.
bpgr = breakpointRanges(vcf_somatic, unpartneredBreakends=FALSE)
n_samples = ncol(vcf_somatic)
vafs_bp = lapply(1:n_samples, function(i) gridss_bp_af(bpgr, vcf_somatic, i)) #i (ordinal) gives you the vaf of that column
vafs_bp = do.call(cbind, vafs_bp)
vafs_bp = round(vafs_bp, 3)
colnames(vafs_bp) = samples(header(vcf_somatic))
geno(vcf_somatic)$VAF = vafs_bp
#info(vcf_bp)$VAF =  apply(vafs_bp, 1, paste, collapse = ";")

vcf_vaf = vcf_somatic

#Filter on vaf
vafs = geno(vcf_vaf)$VAF
max_vaf = apply(vafs, 1, max)
vaf_f = max_vaf >= vaf_lim

vcf_clonal = vcf_vaf[vaf_f]

#For breakpoints check if the partner was kept. If not then remove the breakpoint.
bpgr = breakpointRanges(vcf_clonal, unpartneredBreakends=FALSE)
bp_f = names(vcf_clonal) %in% names(bpgr)
vcf_clonal = vcf_clonal[bp_f]


writeVcf(vcf_clonal, paste0(wdir, sample, ".sv.vaf.vcf"))



####_________________________________________Create circos plots________________________________________####
#Set sequence lengths
hg19_gr = GRanges(seqnames = c(1:22, "X"), ranges = IRanges(start = rep(1, 23), end = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560)))
seqlengths(hg19_gr) = end(hg19_gr)

bpgr = breakpointRanges(vcf_clonal, unpartneredBreakends=FALSE)
present = geno(vcf_clonal)$QUAL != 0

for (i in 1:ncol(present)){
    pres_f = present[,i]
    cell = gsub("_.*", "", colnames(present)[i])
    
    if (sum(pres_f) == 0){
        print(paste0("No somatic SVs for sample: ", cell))
        next
    }
    bpgr_cell = bpgr[pres_f]
    bpgr_cell <- bpgr_cell[seqnames(bpgr_cell) %in% seqlevels(hg19_gr)] # Remove breakends, outside the specified genome region.
    bpgr_cell = bpgr_cell[bpgr_cell$partner %in% names(bpgr_cell)] #Remove breakpoints for which the partner was considered not present.

    if (length(bpgr_cell) == 0){
        print(paste0("No somatic SVs for sample: ", cell))
        next
    }
    gr.circos = bpgr_cell
    seqlevels(gr.circos) <- seqlevels(hg19_gr)
    mcols(gr.circos)$to.gr <- granges(partner(gr.circos))
    
    circos_fig = ggbio() +
        circle(gr.circos, geom="link", linked.to="to.gr") +
        circle(hg19_gr, geom='ideo', fill='gray70') +
        circle(hg19_gr, geom='scale', size=2) +
        circle(hg19_gr, geom='text', aes(label=seqnames), vjust=0, size=3)
    pdf(paste0(wdir, sample, "_", cell, "_circos.pdf"))
    print(circos_fig)
    dev.off()
    
}
