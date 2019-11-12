#mixed model cross validation
library(tidyverse)
library(ggbeeswarm)
library(stringr)
library(nlme)

#Helper function, which removes a row from a dataset and runs a formula on it. Results of interest are then selected and returned.
lme_sub = function(n, index, data, new_call, sample, my_vars){
    if (sample != F){
        sample = data %>% dplyr::slice(n) %>% dplyr::pull(!!sample)
    } else{
        sample = n
    }
    sample %<>% paste(collapse = ";")
    sub_index = index[-n]
    sub_data = data[sub_index,]
    new_m = eval(new_call)
    sum_m = summary(new_m)
    vals = sum_m$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% dplyr::select(variable, my_vars) %>% mutate(sample = sample)
    return(vals)
}

#Function to test the effect of removing points on a model.
leave_n_out = function(m, my_vars = c("p-value", "Value"), nr_removed = 1, sample = F){
    
    #Creates a new formula, that refers to data missing one sample.
    call = deparse(m$call) %>% paste0(collapse = "")
    new_call = sub("data = .*?, ", "data = sub_data, ", call) %>% parse(text = .)
    
    #Run the lme_sub function, which runs the formula.
    data = m$data %>% as_tibble()
    n_cells = nrow(data)
    index = seq_len(n_cells)
    to_remove = combn(index, nr_removed, simplify = F)
    vals_l = lapply(to_remove, function(x) lme_sub(x, index, data, new_call, sample = sample, my_vars))
    
    vals = do.call(rbind, vals_l) %>% as_tibble()
    return(vals)
}

#Plot p values of the model.
plot_n_out_pval = function(vals, m, p_name = "p-value"){
    #my_vars = vals %>% dplyr::select(-variable, -sample) %>% colnames()
    sum_m_main = summary(m)
    vals_main = sum_m_main$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% dplyr::select(variable, !!p_name) %>% mutate(sample = "original")
    
    n_points = m$data %>% nrow()
    nr_left_out = str_count(vals$sample, ";") %>% min() %>% add(1)
    #n_samples = vals$sample %>% unique() %>% length()
    p_value_fig = ggplot(data = vals, aes(x = variable, y = get(p_name))) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(size = 0.8) +
        geom_point(data = vals_main, color = "red") +
        geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
        labs(x = "Term", y = p_name, title = paste0("Leave ", nr_left_out, " out validation (n = ", n_points, ")"))
    return(p_value_fig)
}

#Plot a specified column of the vals tibble. The original values, are plotted as a red point.
plot_n_out = function(vals, m, value = "Value"){
    sum_m_main = summary(m)
    vals_main = sum_m_main$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% dplyr::select(variable, !!value) %>% mutate(sample = "original")
    
    n_points = m$data %>% nrow()
    nr_left_out = str_count(vals$sample, ";") %>% min() %>% add(1)
    value_fig = ggplot(data = vals, aes(x = variable, y= get(value))) +
        geom_boxplot(outlier.shape = NA) +
        geom_quasirandom(size = 0.8, groupOnX = T) +
        geom_point(data = vals_main, color = "red") +
        #geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
        labs(x = "Term", y = value, title = paste0("Leave ", nr_left_out, " out validation (n = ", n_points, ")")) +
        coord_flip()
    return(value_fig)
}



#Function to calculate, how many combinations can be made when leaving out r observations, from a total of n. If this number becomes really high, then your compute will take too long.
nr_combis = function(n, r){
    y = factorial(n) / (factorial(n-r) * factorial(r))
    return(y)
}

# combn(index, n_cells-1, function(x) lme_sub(x, data), simplify = F)
# 
# lme_sub = function(sub_index, data){
#     sub_data = data[sub_index,]
#     new_m = parse(text = new_call) %>% eval()
#     sum_m = summary(new_m)
#     vals = sum_m$tTable %>% as.data.frame() %>% rownames_to_column(var = "variable") %>% dplyr::select(variable, Value, pval = `p-value`)
#     return(vals)
# }
