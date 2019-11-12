library(VariantAnnotation)
library(GenomicRanges)
library(ggplot2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)
library(reshape2)

source("~/hpc/pmc_vanboxtel/tools/Wedge_Dirichlet_Gibbs_sampler/ClusterAssignment.R")
source("~/hpc/pmc_vanboxtel/tools/Wedge_Dirichlet_Gibbs_sampler/interconvertMutationBurdens.R")
source("~/hpc/pmc_vanboxtel/tools/Wedge_Dirichlet_Gibbs_sampler/subclone_Dirichlet_Gibbs_sampler_binomial.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/vcfFilterTree_functions.R")

perform_dirichlet = function(vcf_files, fetus_dir, fetus_name, trisomy, gender, bulk){

    df.all <- data.frame()
    for (i in vcf_files){
        sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
        #sample <- gsub("10817", "X10817", sample)
        
        vcf <- readVcf(i, "hg19")
        seqlevels(vcf) <- paste("chr", seqlevels(vcf), sep="")
        muts <- paste(ref(vcf), unlist(alt(vcf)), sep=">")
        
        df.AD <- data.frame(geno(vcf)$AD, check.names = F)
        df.DP <- data.frame(round(apply(df.AD, 1:2, function(x) sum(unlist(x))), 2), check.names = F)
        df.AA <- data.frame(round(apply(df.AD, 1:2, function(x) sum(unlist(x)) - unlist(x)[1]), 2), check.names = F)
        
        df.mut <- data.frame(sample = sample,
                             chr = seqnames(vcf),
                             pos = start(vcf),
                             mut = muts,
                             DP = df.DP[,sample],
                             AA = df.AA[,sample])
        df.all <- rbind(df.all, df.mut)
    }
    hist(df.all$DP, breaks = 1000) #filter to lowly covered positions
    df.all <- df.all[which(df.all$DP >= 20),]
    
    
    
    
    for (i in unique(df.all$sample)){
        
        totalReads <- df.all[which(df.all$sample == i),]$DP
        mutReads <- df.all[which(df.all$sample == i),]$AA
        
        GS.data <- subclone.dirichlet.gibbs(y=mutReads,N=totalReads)
        
        density.file = paste(i, "density.txt", sep="_")
        Gibbs.subclone.density.est(GS.data, paste(i, "DirichletProcessplotBinomial.png", sep="_"),
                                   density.file = density.file,
                                   post.burn.in.start = 300,
                                   post.burn.in.stop = 1000,
                                   y.max=50)
        
        cluster.assignment <- getClusterAssignments(GS.data, density.file = density.file)
        
        write.table(cbind(1:length(cluster.assignment$localOptima),cluster.assignment$localOptima),
                    paste(i, "clusterPositions.txt", sep="_"),
                    quote=F,
                    sep="\t",
                    col.names=c("cluster_number","cluster_position"),
                    row.names=F)
        
        write.table(cluster.assignment$most.likely.cluster,
                    paste(i, "clusterAssignments.txt", sep="_"),
                    quote=F,
                    sep="\t",
                    col.names="cluster_number",
                    row.names=F)
    }
    
    #process DPC
    
    df.data <- data.frame()
    df.density <- data.frame()
    df.polygon <- data.frame()
    
    for (i in unique(df.all$sample)){
        df.mut <- df.all[which(df.all$sample ==i),]  
        clusters <- data.frame(read.table(paste(i, "clusterAssignments.txt", sep="_"), header=T))
        df.mut <- cbind(df.mut, clusters)
        df.data <- rbind(df.data, df.mut) 
        
        density <- data.frame(read.table(paste(i, "density.txt", sep="_"), header=T))
        polygon <- data.frame(read.table(paste(i, "DirichletProcessplotBinomialpolygonData.txt", sep="_"), header=T))
        polygon <- data.frame(x=c(density$mutation.burden, rev(density$mutation.burden)),
                              y=polygon$x)
        
        density$sample <- i
        df.density <- rbind(df.density, density)
        
        polygon$sample <- i
        df.polygon <- rbind(df.polygon, polygon)
    }
    
    for (i in unique(df.data$sample)){
        plot.dis <- ggplot(df.data[which(df.data$sample == i),], aes(x = AA/DP)) +
            geom_histogram(aes(y = ..count..),
                           binwidth = .025,
                           color = "grey",
                           fill = "light grey",
                           size = .2) +
            geom_line(data = df.density[which(df.density$sample == i),], aes(x = mutation.burden , y = (nrow(df.data[which(df.data$sample == i),]) * 0.025) * median.density), color = "black", size = .5) +
            geom_polygon(data = df.polygon[which(df.polygon$sample == i),], aes(x = x, y = (nrow(df.data[which(df.data$sample == i),]) * 0.025) * y, alpha = .1), fill = "#FFCCFF", color = "#FFCCFF", size = .2) +
            xlab("Variant allele frequency") +
            ylab("Number of base substitutions") +
            coord_cartesian(ylim=c(0, 1000)) +
            scale_y_continuous(breaks=seq(0, 1000, 200)) +
            geom_vline(aes(xintercept=.3), col = "red", linetype="dashed") +
            theme_bw() +
            theme(axis.title.x = element_text(size = 8),
                  axis.text.x = element_text(size = 6),
                  axis.title.y = element_text(size = 8,vjust = 1),
                  axis.text.y = element_text(size = 6),
                  axis.ticks = element_line(size = .3),
                  axis.ticks.length = unit(.6, "mm"),
                  strip.text.x = element_text(size = 8, face = "bold"),
                  strip.text.y = element_text(size=8, face = "bold"),
                  #panel.border = element_rect(colour = "black", fill=NA, size=.2),
                  #panel.border = element_blank(),
                  #panel.grid.major = element_blank(), 
                  #panel.grid.minor = element_blank(),
                  #panel.background = element_blank(),
                  axis.line = element_line(colour = "black"),
                  legend.position = "none")
        ggsave(plot.dis, file= paste(i, "DPCplot.pdf", sep = "_"), dpi = 300, width = 2, height = 2, useDingbats=FALSE)
    } 
    
    #Create vcf with clonal muts, that are not in the shared.vcf (This part is written by Freek)
    for (i in vcf_files){
        sample <- unlist(strsplit(tail(unlist(strsplit(i, '/', fixed = T)), n = 1), '_', fixed = T))[1]
        #sample <- gsub("10817", "X10817", sample)
        
        vcf <- readVcf(i, "hg19")
        # clust_pos_fname = paste0(sample, "_clusterPositions.txt")
        # clust_pos = read_tsv(clust_pos_fname)
        # clonal_clusts = clust_pos %>% dplyr::filter(cluster_position > 0.4 & cluster_position < 0.6) %>% pull(cluster_number)
        # clust_assign_fname = paste0(sample, "_clusterAssignments.txt")
        # clust_assign = read_tsv(clust_assign_fname)
        # muts_clonal_f = clust_assign$cluster_number %in% clonal_clusts
        
        chroms_vcf = rowRanges(vcf) %>% seqnames() %>% as.vector()
        other_chroms = unique(chroms_vcf)
        
        if (trisomy != "FALSE"){
            chr_tri = chroms_vcf == trisomy
            vcf_tri = vcf[chr_tri,]
            ads = geno(vcf_tri)$AD[,sample]
            vaf = sapply(ads, function(x) x[[2]]/ sum(x))
            vaf_f = vaf > 0.2
            vcf_tri = vcf_tri[vaf_f,]
            other_chroms = other_chroms[!other_chroms == trisomy]
        } else{
            vcf_tri = vcf[0]
        }
        
        chr_rest = chroms_vcf %in% other_chroms
        vcf_rest = vcf[chr_rest,]
        ads = geno(vcf_rest)$AD[,sample]
        vaf = sapply(ads, function(x) x[[2]]/ sum(x))
        vaf_f = vaf > 0.3
        vcf_rest = vcf_rest[vaf_f,]
        
        vcf = rbind(vcf_tri, vcf_rest)
        vcf = sort_vcf(vcf)
        
        # if(trisomy != F){
        # chr_tri = rowRanges(vcf) %>% seqnames() %>% as.vector() == trisomy #trisomy muts are not removed because of their dirichlet filtering.
        # vcf = vcf[muts_clonal_f | chr_tri]
        # } else{
        #     vcf = vcf[muts_clonal_f]
        # }
        
        #Add muts on the x chromosome for males
        if (gender == "M"){
            vcf_x_fname = paste0(fetus_dir, sample, "_X/", sample, "_", bulk, "_Q50_CGQ10_SGQ10_PASS_10X_VAF0.99_X_nonBlacklist_final.vcf")
            vcf_x = readVcf(vcf_x_fname, genome = "hg19")
            vcf = rbind(vcf, vcf_x)
        }
        
        #Remove muts that are in the shared part of the tree
        vcf_shared_fname = paste0(fetus_dir, "shared/", fetus_name, "_shared.vcf")
        if (file.exists(vcf_shared_fname)){
            vcf_shared = readVcf(vcf_shared_fname, genome = "hg19")
            not_unique_clone = findOverlaps(granges(vcf), granges(vcf_shared)) %>% queryHits() %>% unique()
            if (length(not_unique_clone) > 0){
                vcf = vcf[-not_unique_clone]
            }
        }
        vcf_out_fname = paste0(sample, "_notbulk.vcf")
        bed_validated_fname = paste0(sample, "_notbulk_validated.bed")
        bed_tovalidate_fname = paste0(sample, "_notbulk_tovalidate.bed")
        trash = vcf_out(vcf, vcf_out_fname, bed_validated_fname, bed_tovalidate_fname)
    }
}
#     
#     #In vitro spectrum
#     df.all$cat <- "In vivo"
#     df.all[which(df.all$AA/df.all$DP < 0.3),]$cat <- "In vitro"
#     
#     df.all$sample <- gsub("X10817", "10817", df.all$sample)
#     
#     write.table(df.all, "Invitro_muts.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#     
#     for (i in 1:nrow(df.all)){
#         if(grepl("G", substring(df.all[i,"mut"], 1, 1)) == TRUE)
#             df.all[i,"strand"] <- "-"
#         else if(grepl("A", substring(df.all[i,"mut"], 1, 1)) == TRUE)
#             df.all[i,"strand"] <- "-"
#         else
#             df.all[i, "strand"] <- "+"
#     }
#     
#     df.all$triplet <- as.character(getSeq(Hsapiens,
#                                           df.all$chr,
#                                           start = df.all$pos - 1,
#                                           end = df.all$pos + 1,
#                                           strand=df.all$strand))
#     df.all$context <- paste(substring(df.all$triplet, 1,1),
#                             substring(df.all$triplet,3),
#                             sep=".")
#     
#     #--add mut types
#     df.all$type <- df.all$mut
#     df.all$type <- gsub('G>T', 'C>A', df.all$type)
#     df.all$type <- gsub('G>C', 'C>G', df.all$type)
#     df.all$type <- gsub('G>A', 'C>T', df.all$type)
#     df.all$type <- gsub('A>T', 'T>A', df.all$type)
#     df.all$type <- gsub('A>G', 'T>C', df.all$type)
#     df.all$type <- gsub('A>C', 'T>G', df.all$type)
#     
#     df.cluster <- melt(table(df.all$type, df.all$context, df.all$cat))
#     colnames(df.cluster) <- c("type", "context", "cat", "count")
#     
#     #df.cluster$norm_count <- df.cluster$count/sum(df.cluster[which(df.cluster$cluster == 3),]$count)
#     #df.cluster[which(df.cluster$cluster == 1),]$norm_count <- df.cluster[which(df.cluster$cluster == 1),]$count/sum(df.cluster[which(df.cluster$cluster == 1),]$count)
#     
#     df.cluster <- df.cluster[order(df.cluster$cat, df.cluster$type),]
#     
#     df.sig <- data.frame(type = unique(paste(df.cluster$type, df.cluster$context, sep = "_")),
#                          invivo = df.cluster[which(df.cluster$cat == "In vivo"),4],
#                          invitro = df.cluster[which(df.cluster$cat == "In vitro"),4])
#     
#     plot_96_profile(df.sig[,2:3])
#     
#     #Refit
#     
#     cancer_signatures <- read.table("~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_signatures.csv", sep = ",", header = T)
#     cancer_signatures <- cancer_signatures[,3:67]
#     
#     mut_mat <- as.matrix(df.sig[,2:3])
#     
#     fit_res <- fit_to_signatures(mut_mat, as.matrix(cancer_signatures[,-45]))
#     plot_contribution(fit_res$contribution, as.matrix(cancer_signatures))
#     cos_sim_sigs <- cos_sim_matrix(mut_mat, fit_res$reconstructed)
#     sim_all <-diag(cos_sim_sigs)
#     
#     sort(fit_res$contribution[,2])
#     abline(h = 0.10)
#     abline(h = 0.05)
#     
#     fit_res3 <- fit_to_signatures(mut_mat, as.matrix(cancer_signatures[,c("SBS5","SBS32","SBS1")]))
#     plot_contribution(fit_res3$contribution,
#                       as.matrix(cancer_signatures[,c("SBS5","SBS32","SBS1")]),
#                       coord_flip = FALSE,
#                       mode = "relative")
#     
#     cos_sim_sigs3 <- cos_sim_matrix(mut_mat, fit_res3$reconstructed)
#     sim3 <- diag(cos_sim_sigs3)
# }