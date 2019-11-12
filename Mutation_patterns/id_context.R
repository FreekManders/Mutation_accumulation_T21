# Written by Freek Manders, 2019
remove_snvs_gr = function(gr, sample_name = "sample"){
    snv_f = width(gr$REF) == 1 & width(unlist(gr$ALT)) == 1
    nr_snv = sum(snv_f)
    gr = gr[!snv_f]
    if (nr_snv >= 1){
        print(paste0("Removed ", nr_snv, " SNVs from sample: ", sample_name))
    }
    return(gr)
}


remove_snvs = function(grl){
    if (class(grl)[[1]] == "CompressedGRangesList"){
        gr_list = mapply(remove_snvs_gr, grl, names(grl))
        grl = GRangesList(gr_list)
        return(grl)
    } else if (class(grl)[[1]] == "GRanges"){
        gr = remove_snvs_gr(grl)
        return(gr)
    }
}


indel_check = function(grl){
    if (class(grl)[[1]] == "CompressedGRangesList"){
        gr = unlist(grl)
    } else if (class(grl)[[1]] == "GRanges"){
        gr = grl
    }
    trash = indel_check_gr(gr)
    return(0)
}

indel_check_gr = function(gr){
    snv_f = width(gr$REF) == 1 & width(unlist(gr$ALT)) == 1
    nr_snv = sum(snv_f)
    snv_present = nr_snv >= 1
    if (snv_present){
        stop(paste0("There seem to be ", nr_snv, " SNVs present in your data. Make sure to remove all SNVs with the remove_snvs function."))
    }
}


#1bp deletions
get_1bp_dels = function(gr, mut_size, ref_genome){
    gr = gr[mut_size == -1]
    if (length(gr) == 0){
        return(gr)
    }
    
    del_bases = gr$REF %>% as.character() %>% substring(2)
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        gr_extended = flank(gr, 19, start = F)
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    gr_extended %<>% trim() #Trim the ranges that are extended beyond the actual length of the chromosome.
    seq = getSeq(eval(as.symbol(ref_genome)), gr_extended)
    
    seq_z = str_replace_all(as.character(seq), del_bases, rep("Z", length(seq))) #For each mut replace the deleted basetype in the flanking sequence with Zs.
    homopolymer_length = gsub("[^Z].*", "", as.character(seq_z)) %>% nchar() + 1 #Remove all bases after the Zs and count how many bases are left. Add one for the deleted base itself.
    del_bases[del_bases == "A"] = "T"
    del_bases[del_bases == "G"] = "C"
    
    gr$muttype = paste0(del_bases, "_deletion")
    gr$muttype_sub = homopolymer_length
    return(gr)
}
#1bp insertions
get_1bp_ins = function(gr, mut_size, ref_genome){
    gr = gr[mut_size == 1]
    if (length(gr) == 0){
        return(gr)
    }
    
    
    ins_bases = gr$ALT %>% unlist() %>% as.character() %>% substring(2)
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        gr_extended = flank(gr, 20, start = F)
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    gr_extended %<>% trim() #Trim the ranges that are extended beyond the actual length of the chromosome.
    seq = getSeq(eval(as.symbol(ref_genome)), gr_extended)
    seq_z = str_replace_all(as.character(seq), ins_bases, rep("Z", length(seq))) #For each mut replace the inserted basetype in the flanking sequence with Zs.
    homopolymer_length = gsub("[^Z].*", "", as.character(seq_z)) %>% nchar() #Remove all bases after the Zs and count how many bases are left.
    ins_bases[ins_bases == "A"] = "T"
    ins_bases[ins_bases == "G"] = "C"
    
    gr$muttype = paste0(ins_bases, "_insertion")
    gr$muttype_sub = homopolymer_length
    return(gr)
}

#Bigger insertions
get_big_ins = function(gr, mut_size, ref_genome){
    gr = gr[mut_size > 1]
    if (length(gr) == 0){
        return(gr)
    }
    
    
    mut_size = mut_size[mut_size > 1]
    
    ins_bases = gr$ALT %>% unlist() %>% as.character() %>% substring(2)
    biggest_ins = ins_bases %>% nchar() %>% max()
    flank_dist = biggest_ins * 20
    
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        gr_extended = flank(gr, flank_dist, start = F)
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    gr_extended %<>% trim() #Trim the ranges that are extended beyond the actual length of the chromosome.
    
    seq = getSeq(eval(as.symbol(ref_genome)), gr_extended)
    seq_z = str_replace_all(as.character(seq), ins_bases, rep("Z", length(seq))) #For each mut replace the deleted basetype in the flanking sequence with Zs.
    n_repeats = gsub("[^Z].*", "", as.character(seq_z)) %>% nchar() #Remove all bases after the Zs and count how many bases are left.
    
    gr$muttype = paste0(mut_size, "bp_insertion")
    gr$muttype_sub = n_repeats
    return(gr)
}

#Bigger deletions
get_big_dels = function(gr, mut_size, ref_genome){
    gr = gr[mut_size < -1]
    if (length(gr) == 0){
        return(gr)
    }
    
    mut_size = mut_size[mut_size < -1]
    
    del_bases = gr$REF %>% as.character() %>% substring(2)
    biggest_dels = del_bases %>% nchar() %>% max()
    flank_dist = biggest_dels * 20
    
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        gr_extended = flank(gr, flank_dist, start = F)
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    gr_extended %<>% trim() #Trim the ranges that are extended beyond the actual length of the chromosome.
    
    seq = getSeq(eval(as.symbol(ref_genome)), gr_extended)
    seq_z = str_replace_all(as.character(seq), del_bases, rep("Z", length(seq))) #For each mut replace the deleted basetype in the flanking sequence with Zs.
    n_repeats = gsub("[^Z].*", "", as.character(seq_z)) %>% nchar() + 1 #Remove all bases after the Zs and count how many bases are left. Add +1 for the deleted bases itself
    
    gr$muttype = paste0(abs(mut_size), "bp_deletion")
    gr$muttype_sub = n_repeats
    
    
    #Determine if there is microhomology for deletions that are not in repeat regions.
    pos_mh = gr$muttype_sub == 1 #There is always at least 1 'repeat', because of the deleted bases themselves.
    
    gr_repeat = gr[!pos_mh]
    gr_mh = gr[pos_mh]
    
    if (length(gr_mh) == 0){
        return(gr_repeat)
    }
    
    mut_size_mh = mut_size[pos_mh]
    del_bases_mh = del_bases[pos_mh]
    del_bases_s = strsplit(del_bases_mh, "")
    seq_s = strsplit(as.character(seq[pos_mh]), "")
    
    
    #Also check for microhomology to the left of the deletion. For this take the reverse sequence to the left of the deletion and the reverse deleted bases.
    rev_del_bases = reverse(del_bases_mh)
    rev_l_del_bases_s = strsplit(rev_del_bases, "")
    
    withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
        l_flank = flank(gr_mh, biggest_dels)
    }, warning = function(w) {
        if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
            invokeRestart("muffleWarning")
    })
    l_flank %<>% trim() %>% shift(1) #Trim the ranges that are extended beyond the actual length of the chromosome. #Add 1 base, because the first base in the granges obj is not deleted and should be used in the flank.
    rev_l_seq = getSeq(eval(as.symbol(ref_genome)), l_flank) %>% reverse()
    rev_l_seq_s = strsplit(as.character(rev_l_seq), "")
    
    #For each mutation determine how many bases show hm
    nr_pos_mh = length(del_bases_s)
    nr_mh = vector("list", nr_pos_mh)
    for (i in 1:nr_pos_mh){
        del_bases_sample = del_bases_s[[i]]
        seq_s_sample = seq_s[[i]][1:length(del_bases_sample)]
        same = del_bases_sample == seq_s_sample
        r_nr_mh_sample = cumprod(same) %>% sum(na.rm = T) #Determine how many bases are the same before the first difference. na.rm is for when a sequence has been trimmed.
        
        l_del_bases_sample = rev_l_del_bases_s[[i]]
        l_seq_s_sample = rev_l_seq_s[[i]][1:length(l_del_bases_sample)]
        l_same = l_del_bases_sample == l_seq_s_sample
        l_nr_mh_sample = cumprod(l_same) %>% sum(na.rm = T)
        
        nr_mh_sample = max(r_nr_mh_sample, l_nr_mh_sample)
        
        nr_mh[[i]] = nr_mh_sample
    }
    nr_mh = unlist(nr_mh)
    
    #Update gr when mh is indeed present
    mh_f = nr_mh > 0
    gr_mh$muttype[mh_f] = paste0(abs(mut_size_mh[mh_f]), "bp_deletion_with_microhomology")
    gr_mh$muttype_sub[mh_f] = nr_mh[mh_f]
    
    #Combine muts with and without mh
    gr = c(gr_mh, gr_repeat) %>% sort()
    return(gr)
}


#Function to get all the indel contexts of a gr object. This gr object should contain only indels
get_indel_context_gr = function(gr, ref_genome){
    
    mult_alts = elementNROWS(gr$ALT) > 1
    nr_mult_alts = sum(mult_alts)
    if (nr_mult_alts > 0){
        gr = gr[!mult_alts]
        warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
    }
    
    #Calculate indel size to determine main categorie
    ref_sizes = gr$REF %>% width()
    alt_sizes = gr$ALT %>% unlist() %>% width()
    mut_size = alt_sizes - ref_sizes
    
    #For the main indel categories, determine their sub categories. (Also split the big deletion categorie into repeat and micro homology.)
    gr_1b_dels = get_1bp_dels(gr, mut_size, ref_genome)
    gr_1b_ins = get_1bp_ins(gr, mut_size, ref_genome)
    gr_big_dels = get_big_dels(gr, mut_size, ref_genome)
    gr_big_ins = get_big_ins(gr, mut_size, ref_genome)
    
    gr = c(gr_1b_dels, gr_1b_ins, gr_big_dels, gr_big_ins) %>% sort()
    return(gr)
}

#Convenience function to retreive the indel contexts of a GRangeslist. It can also be used on a GRanges object itself.
get_indel_context = function(grl, ref_genome){
    if (class(grl)[[1]] == "CompressedGRangesList"){
        gr_list = lapply(grl, function(x) get_indel_context_gr(x, ref_genome))
        grl = GRangesList(gr_list)
        return(grl)
    } else if (class(grl)[[1]] == "GRanges"){
        gr = get_indel_context_gr(grl, ref_genome)
        return(gr)
    }
}

count_indel_contexts_gr = function(gr){
    categories = tibble("muttype" = c(rep("C_deletion", 6), rep("T_deletion", 6), rep("C_insertion", 6), rep("T_insertion", 6), rep("2bp_deletion", 6), rep("3bp_deletion", 6), rep("4bp_deletion", 6), rep("5+bp_deletion", 6), rep("2bp_insertion", 6),rep("3bp_insertion", 6),rep("4bp_insertion", 6),rep("5+bp_insertion", 6), rep("2bp_deletion_with_microhomology", 1), rep("3bp_deletion_with_microhomology", 2), rep("4bp_deletion_with_microhomology", 3), rep("5+bp_deletion_with_microhomology", 5)), "muttype_sub" =  c(rep(c(1:5, "6+"), 2), rep(c(0:4, "5+"), 2), rep(c(1:5, "6+"), 4), rep(c(0:4, "5+"), 4), 1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
    if (length(gr) == 0){
        categories %<>% mutate(count = 0) %>% dplyr::select(-muttype, -muttype_sub)
        return(categories)
    }
    
    id_context = tibble("muttype" = gr$muttype, "muttype_sub" = gr$muttype_sub)
    
    #Classify the number of repeat units/ homopolymer length / microhomology length to either 5+ or 6+ depending on whether the indel is a ins or del.
    id_context[id_context$muttype_sub >= 6, "muttype_sub"] = "6+"
    id_context[grepl("insertion|microhomology", id_context$muttype) & id_context$muttype_sub >= 5, "muttype_sub"] = "5+"
    
    #Classify large indels as size 5+
    ref_sizes = gr$REF %>% width()
    alt_sizes = gr$ALT %>% unlist() %>% width()
    mut_size = abs(alt_sizes - ref_sizes)
    mut_size_f = mut_size >= 5
    id_context$muttype = ifelse(mut_size_f, gsub("[0-9]+bp", "5+bp", id_context$muttype, perl = T), id_context$muttype)
    
    id_context_count = id_context %>% group_by(muttype, muttype_sub) %>% dplyr::summarise(count = n())
    id_context_count_full = left_join(categories, id_context_count, by = c("muttype", "muttype_sub")) %>% dplyr::select(-muttype, -muttype_sub)
    #colnames(id_context_count_full)[3] = name
    return(id_context_count_full)
}

count_indel_contexts = function(grl){
    if (class(grl)[[1]] == "CompressedGRangesList"){
        
        counts_l = lapply(grl, count_indel_contexts_gr)
        counts = do.call(cbind, counts_l)
        colnames(counts) = names(grl)
        
    } else if (class(grl)[[1]] == "GRanges"){
        counts = count_indel_contexts_gr(grl)
        colnames(counts) = "My_sample"
    }
    categories = tibble("muttype" = c(rep("C_deletion", 6), rep("T_deletion", 6), rep("C_insertion", 6), rep("T_insertion", 6), rep("2bp_deletion", 6), rep("3bp_deletion", 6), rep("4bp_deletion", 6), rep("5+bp_deletion", 6), rep("2bp_insertion", 6),rep("3bp_insertion", 6),rep("4bp_insertion", 6),rep("5+bp_insertion", 6), rep("2bp_deletion_with_microhomology", 1), rep("3bp_deletion_with_microhomology", 2), rep("4bp_deletion_with_microhomology", 3), rep("5+bp_deletion_with_microhomology", 5)), "muttype_sub" =  c(rep(c(1:5, "6+"), 2), rep(c(0:4, "5+"), 2), rep(c(1:5, "6+"), 4), rep(c(0:4, "5+"), 4), 1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
    counts = cbind(categories, counts)
    counts[is.na(counts)] = 0
    counts = as_tibble(counts)
    counts$muttype = factor(counts$muttype, levels = unique(counts$muttype))
    return(counts)
}

plot_main_indel_contexts = function(counts, same_y = F){
    counts_main = counts %>% dplyr::select(-muttype_sub) %>% group_by(muttype) %>% summarise_all(list( ~sum(.)))
    counts_main %<>% gather(key = "sample", value = "count", -muttype)
    nr_muts = counts_main %>% group_by(sample) %>% dplyr::summarise(nr_muts = sum(count))
    facet_labs_y = paste0(nr_muts$sample, " (n = ", nr_muts$nr_muts, ")")
    names(facet_labs_y) = nr_muts$sample
    colors = c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B")
    if (same_y){
        facet_scale = "free_x"
    } else{
        facet_scale = "free"
    }
    fig = ggplot(counts_main, aes(x = muttype, y = count, fill = muttype)) +
        geom_bar(stat = "identity") +
        facet_grid(sample ~ ., labeller = labeller(sample = facet_labs_y), scale = facet_scale) +
        labs(x = "", y = "Nr of indels") +
        scale_fill_manual(guide=FALSE, values = colors) +
        theme(axis.text.x = element_text(angle = 90))
    return(fig)
}


plot_indel_contexts = function(counts, same_y = F){
    counts %<>% gather(key = "sample", value = "count", -muttype, -muttype_sub)
    nr_muts = counts %>% group_by(sample) %>% dplyr::summarise(nr_muts = round(sum(count)))
    
    if (same_y){
        facet_scale = "free_x"
    } else{
        facet_scale = "free"
    }
    
    facet_labs_y = paste0(nr_muts$sample, " (n = ", nr_muts$nr_muts, ")")
    names(facet_labs_y) = nr_muts$sample
    #counts %<>% mutate(Category = ifelse(grepl("deletion_with", muttype), "Deletion With Microhomology", ifelse(grepl("deletion", muttype), "Deletion", "Insertion")))
    facet_labs_x = c("1: C", "1: T", "1: C", "1: T", 2 ,3, 4, "5+", 2, 3, 4, "5+", 2, 3, 4, "5+")
    names(facet_labs_x) = levels(counts$muttype)
    colors = c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B")
    fig = ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype)) +
        geom_bar(stat = "identity") +
        facet_grid(sample ~ muttype, scale = facet_scale, space = "free_x", labeller = labeller(muttype = facet_labs_x, sample = facet_labs_y)) +
        theme_minimal() +
        scale_fill_manual(values = colors) +
        theme(panel.grid.major.x = element_blank(), strip.background =element_rect(fill="cadetblue")) +
        labs(fill = "Mutation type", title = "Deletion         Insertion        Deletion                              Insertion                              Deletion (MH)", y = "Nr of indels",
             x = "Homopolymer length                         Number of repeat units                                                                Microhomology length")
    
    return(fig)
}

plot_indel_sigs = function(id_sigs, sig_features_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/ID/sig_features.txt"){
    id_sig_features = read_tsv(sig_features_fname)
    id_sigs = cbind(id_sig_features, id_sigs)
    id_sigs$muttype = factor(id_sigs$muttype, levels = unique(id_sigs$muttype))
    sig_cols = 3:ncol(id_sigs)
    sig_figs = lapply(sig_cols, function(x){
        fig = plot_indel_contexts(id_sigs[,c(1,2,x)]) + 
            coord_cartesian(ylim = c(0,1), expand = F)
        return(fig)
    })
    return(sig_figs)
}

