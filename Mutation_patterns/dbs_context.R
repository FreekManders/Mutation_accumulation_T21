

#Identify double base substitutions. Takes a granges object as its input. The input should contain only snvs and have been filtered on quality.
get_dbs = function(gr){
    
    #Identify location of each mut and its subsequent mut.
    chroms = seqnames(gr) %>% as.vector()
    first_chrom = chroms[-length(chroms)]
    second_chrom = chroms[-1]
    coords = start(gr)
    first_coords = coords[-length(coords)]
    second_coords = coords[-1]
    
    #Compare location of mut and subsequent mut. If they are only one base apart then they are a dbs.
    dbs_locs = which(first_chrom == second_chrom & first_coords + 1 == second_coords)
    
    if (length(dbs_locs) == 0){
        gr = GRanges()
        genome(gr) = "hg19"
        return(gr)
    }
    
    #Create seperate granges for the first and second base of each dbs. Then incorporate the info from the second base into the first grange.
    gr_first = gr[dbs_locs]
    gr_second = gr[dbs_locs+1]
    end(gr_first) = end(gr_second)
    
    first_ref = gr_first$REF %>% as.vector()
    second_ref = gr_second$REF %>% as.vector()
    gr_first$REF = paste0(first_ref, second_ref) %>% DNAStringSet()
    
    first_alt = gr_first$ALT %>% unlist() %>% as.vector()
    second_alt = gr_second$ALT %>% unlist() %>% as.vector()
    gr_first$ALT = paste0(first_alt, second_alt) %>% as.list() %>% DNAStringSetList()
    
    gr_first$QUAL = cbind(gr_first$QUAL, gr_second$QUAL) %>% rowMeans()
    return(gr_first)
}

#Set the correct features for the dbs muts. AA>NN is for examle changed to TT>NN
set_context_dbs = function(gr){
    ref = gr$REF
    rev_main = ref %in% c("AA", "GG", "AG", "CA", "GA", "GT")
    ref[rev_main] = reverseComplement(ref[rev_main])
    
    alt = gr$ALT %>% unlist()
    alt[rev_main] = reverseComplement(alt[rev_main])
    
    #If the reference is identifcal to its reverse complement then the number of possible alts is reduced to 6. Since A CG>TT is the same as a CG>AA.
    #In these cases the rev complement of the alt has to be taken in some cases.
    rev_alt1 = which(ref == "TA" & alt %in% c("CC", "AC", "AG"))
    rev_alt2 = which(ref == "AT" & alt %in% c("GG", "TG", "TC"))
    rev_alt3 = which(ref == "GC" & alt %in% c("CT", "TT", "TG"))
    rev_alt4 = which(ref == "CG" & alt %in% c("AA", "AC", "GA"))
    rev_alt = c(rev_alt1, rev_alt2, rev_alt3, rev_alt4)
    alt[rev_alt] = reverseComplement(alt[rev_alt])
    alt %<>% as.vector() %>% as.list() %>% DNAStringSetList()
    gr$REF = ref
    gr$ALT = alt
    return(gr)
}

#Helper function to count how often different dbs features occur for a single gr.
count_dbs_contexts_gr = function(gr){
    categories = tibble("REF" = c(rep("AC", 9), rep("AT", 6), rep("CC", 9), rep("CG", 6), rep("CT", 9), rep("GC", 6), rep("TA", 6), rep("TC", 9), rep("TG", 9), rep("TT", 9)), "ALT" = c("CA", "CG", "CT", "GA", "GG", "GT", "TA", "TG", "TT", "CA", "CC", "CG", "GA", "GC", "TA", "AA", "AG", "AT", "GA", "GG", "GT", "TA", "TG", "TT", "AT", "GC", "GT", "TA", "TC", "TT", "AA", "AC", "AG", "GA", "GC", "GG", "TA", "TC", "TG", "AA", "AG", "AT", "CA", "CG", "TA", "AT", "CG", "CT", "GC", "GG", "GT", "AA", "AG", "AT", "CA", "CG", "CT", "GA", "GG", "GT", "AA", "AC", "AT", "CA", "CC", "CT", "GA", "GC", "GT", "AA", "AC", "AG", "CA", "CC", "CG", "GA", "GC", "GG"))
    
    context = cbind("REF" = as.vector(gr$REF), "ALT" = as.vector(unlist(gr$ALT)))
    counts = context %>% as_tibble() %>% group_by(REF, ALT) %>% dplyr::summarise(count = n())
    counts_full = left_join(categories, counts, by = c("REF", "ALT")) %>% dplyr::select(-REF, -ALT)
    return(counts_full)
}

#Count how often different dbs features occur.
count_dbs_contexts = function(grl){
    if (class(grl)[[1]] == "CompressedGRangesList"){
        
        counts_l = lapply(grl, count_dbs_contexts_gr)
        counts = do.call(cbind, counts_l)
        colnames(counts) = names(grl)
        
    } else if (class(grl)[[1]] == "GRanges"){
        counts = count_indel_contexts_gr(gr)
        colnames(counts) = "My_sample"
    }
    categories = tibble("REF" = c(rep("AC", 9), rep("AT", 6), rep("CC", 9), rep("CG", 6), rep("CT", 9), rep("GC", 6), rep("TA", 6), rep("TC", 9), rep("TG", 9), rep("TT", 9)), "ALT" = c("CA", "CG", "CT", "GA", "GG", "GT", "TA", "TG", "TT", "CA", "CC", "CG", "GA", "GC", "TA", "AA", "AG", "AT", "GA", "GG", "GT", "TA", "TG", "TT", "AT", "GC", "GT", "TA", "TC", "TT", "AA", "AC", "AG", "GA", "GC", "GG", "TA", "TC", "TG", "AA", "AG", "AT", "CA", "CG", "TA", "AT", "CG", "CT", "GC", "GG", "GT", "AA", "AG", "AT", "CA", "CG", "CT", "GA", "GG", "GT", "AA", "AC", "AT", "CA", "CC", "CT", "GA", "GC", "GT", "AA", "AC", "AG", "CA", "CC", "CG", "GA", "GC", "GG"))
    counts = cbind(categories, counts)
    counts[is.na(counts)] = 0
    counts = as_tibble(counts)
    counts$REF = factor(counts$REF, levels = unique(counts$REF))
    
    bases = c("A", "C", "G", "T")
    bases_combi = crossing(bases, bases)
    counts$ALT = factor(counts$ALT, levels = paste0(bases_combi$bases, bases_combi$bases1))
    return(counts)
}

#Plot the main dbs contexts. Comparable to a 6 feature snv profile
plot_main_dbs_contexts = function(counts, same_y = F){
    counts_main = counts %>% dplyr::select(-ALT) %>% group_by(REF) %>% summarise_all(list( ~sum(.)))
    counts_main %<>% gather(key = "sample", value = "count", -REF)
    nr_muts = counts_main %>% group_by(sample) %>% dplyr::summarise(nr_muts = sum(count))
    facet_labs_y = paste0(nr_muts$sample, " (n = ", nr_muts$nr_muts, ")")
    names(facet_labs_y) = nr_muts$sample
    colors = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
    if (same_y){
        facet_scale = "free_x"
    } else{
        facet_scale = "free"
    }
    fig = ggplot(counts_main, aes(x = REF, y = count, fill = REF)) +
        geom_bar(stat = "identity") +
        facet_grid(sample ~ ., labeller = labeller(sample = facet_labs_y), scale = facet_scale) +
        labs(x = "", y = "Nr of DBSs") +
        scale_fill_manual(guide=FALSE, values = colors) +
        theme(axis.text.x = element_text(angle = 90))
    return(fig)
}


#Plot the dbs contexts. Comparable to a 96 feature snv profile
plot_dbs_contexts = function(counts, same_y = F){
    counts %<>% gather(key = "sample", value = "count", -REF, -ALT)
    nr_muts = counts %>% group_by(sample) %>% dplyr::summarise(nr_muts = round(sum(count)))
    
    if (same_y){
        facet_scale = "free_x"
    } else{
        facet_scale = "free"
    }
    facet_labs_y = paste0(nr_muts$sample, " (n = ", nr_muts$nr_muts, ")")
    names(facet_labs_y) = nr_muts$sample
    #counts %<>% mutate(Category = ifelse(grepl("deletion_with", muttype), "Deletion With Microhomology", ifelse(grepl("deletion", muttype), "Deletion", "Insertion")))
    facet_labs_x = paste0(levels(counts$REF), ">NN")
    names(facet_labs_x) = levels(counts$REF)
    colors = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
    fig = ggplot(counts, aes(x = ALT, y = count, fill = REF)) +
        geom_bar(stat = "identity") +
        facet_grid(sample ~ REF, scale = facet_scale, space = "free_x", labeller = labeller(REF = facet_labs_x, sample = facet_labs_y)) +
        scale_fill_manual(guide = F, values = colors) +
        theme_minimal() +
        theme(panel.grid.major.x = element_blank(), strip.background =element_rect(fill="cadetblue"), axis.text.x = element_text(angle = 90)) +
        labs(fill = "Mutation type", title = "", y = "Nr of DBSs", x = "")
    
    return(fig)
}

plot_dbs_sigs = function(dbs_sigs, dbs_features_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/DBS/sig_features.txt"){
    dbs_features = read_tsv(dbs_features_fname)
    dbs_sigs = cbind(dbs_features, dbs_sigs)
    dbs_sigs$REF = factor(dbs_sigs$REF, levels = unique(dbs_sigs$REF))
    dbs_sigs$ALT = factor(dbs_sigs$ALT, levels = unique(dbs_sigs$ALT))
    sig_cols = 3:ncol(dbs_sigs)
    sig_figs = lapply(sig_cols, function(x){
        fig = plot_dbs_contexts(dbs_sigs[,c(1,2,x)]) +
            coord_cartesian(ylim = c(0,1), expand = F)
        return(fig)
    })
    return(sig_figs)
}
