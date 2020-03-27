library(tidyverse)

#Function to create a pie chart
create_pie = function(v, title, name_legend){
    
    #Count effects
    tb = tibble("variable" = v) %>%
        group_by(variable) %>% 
        dplyr::count()
    
    fig = ggplot(tb, aes(y = n, x = "", fill = variable)) +
        geom_bar(width = 1, stat = "identity") +
        geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
        coord_polar("y", start = 0) +
        scale_fill_discrete(drop=TRUE, limits = levels(tb$variable)) +
        labs(title = title, fill = name_legend) +
        theme_minimal() +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #panel.border = element_blank(),
            panel.grid=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks = element_blank(),
            plot.title=element_text(size=14, face="bold")
        )
    return(fig)
}

#Set wd
setwd("~/surfdrive/Shared/Projects/Freek/Freek_trees/")

#Read exon table and combine with sample overview
exon_clones_tb = read_tsv("exon_muts.txt")
overview_samples = read_tsv("overview_samples.txt") %>% 
    dplyr::select(fetus, trisomy) %>% 
    unique()
exon_clones_tb = exon_clones_tb %>% 
    left_join(overview_samples, by = "fetus") %>% 
    dplyr::mutate(trisomy = ifelse(trisomy == 21, "T21", "D21"))

#Add tmd
exon_tmd_tb = read_tsv("~/surfdrive/Shared/Projects/Freek/Freek_tmd/exon_muts.txt") %>% 
    dplyr::mutate(trisomy = "TMD")
exon_tb = rbind(exon_clones_tb, exon_tmd_tb)

#Make pretty
exon_tb = exon_tb %>% 
    dplyr::mutate(effect_type = str_remove(effect_type, "&.*"),
                  effect_type = str_replace_all(effect_type, "_", " "),
                  effect_type = factor(effect_type, levels = unique(effect_type)),
                  effect = factor(effect, levels = unique(effect)))

#Split based on trisomy state
exon_tb_l = split(exon_tb, exon_tb$trisomy)

#Create effect type pie charts
fig_l = exon_tb_l %>% 
    purrr::map("effect_type") %>% 
    purrr::imap(create_pie, "Effect type")
pdf("effect_type_piecharts.pdf", useDingbats = F)
print(fig_l)
dev.off()

#Create effect size figures
fig_l = exon_tb_l %>% 
    purrr::map("effect") %>% 
    purrr::imap(create_pie, "Effect size")
pdf("effect_size_piecharts.pdf", useDingbats = F)
print(fig_l)
dev.off()



