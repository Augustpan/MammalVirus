
library(mgcv)
library(oddsratio)
library(caret)

tar_load("tb_virus_list_clean")
tar_load("mat_virus_fam_pre_mat")
tar_load("mat_host_ord_pre_mat")

{
    virus_list = tb_virus_list_clean
    virus_list$Genome_Type_Bin = ifelse(virus_list$Genome_Type == "DNA", "DNA", "RNA")
    
    virusFamPrecisionMat = mat_virus_fam_pre_mat %>%
        `[`(unique(virus_list$Viral_Family), unique(virus_list$Viral_Family)) %>%
        solve()
    virusFamRank = virusFamPrecisionMat %>%
        rankMatrix() %>%
        as.vector()
    hostOrdPrecisionMat = mat_host_ord_pre_mat %>%
        `[`(c("Insectivora", "Chiroptera", "Rodentia"), c("Insectivora", "Chiroptera", "Rodentia")) %>%
        solve()
    hostOrdRank = hostOrdPrecisionMat %>%
        rankMatrix() %>%
        as.vector()
    

}

{
    # the biological characteristics of virus
    effect_viral_family = "s(Viral_Family, bs = 're', xt = list(S = list(virusFamPrecisionMat), rank = virusFamRank))"
    effect_tissue_tropism_axis1 = "s(TT_CA1, bs = 'tp', k = 10)"
    effect_tissue_tropism_axis2 = "s(TT_CA2, bs = 'tp', k = 10)"
    effect_genome_type = "Genome_Type"
    effect_zoonosis = "Zoonosis"
    effect_novelty = "Novelty"
    
    # the ecological characteristics of virus
    effect_com_viral_richness = "s(Region_Viral_Richness, bs = 'tp', k = 3)"
    effect_pop_viral_richness = "s(Population_Viral_Richness, bs = 'tp', k = 9)"
    effect_viral_com_genus_axis1 = "s(VG_CA1, bs = 'tp', k = 10)"
    effect_viral_com_genus_axis2 = "s(VG_CA2, bs = 'tp', k = 10)"
    
    # the characteristics of host
    effect_host_order = "s(Host_Order, bs = 're')"
    effect_host_order_specific_viral_family = "s(Host_Order, Viral_Family, bs = 're')"
    effect_host_order_phylogeny = "s(Host_Order, bs = 're', xt = list(S = list(hostOrdPrecisionMat), rank = hostOrdRank))"
    effect_com_host_richness = "s(Region_Host_Richness, bs = 'tp', k = 3)"
}


fit = gam(
    Multi_Host ~
        Genome_Type +
        #s(TT_CA1, bs = 'tp', k = 10) + 
        s(TT_CA2, bs = 'tp', k = 10) +
        s(Viral_Family, bs = 're', xt = list(S = list(virusFamPrecisionMat), rank = virusFamRank)) +
        s(VG_CA1, bs = 'tp', k = 10) + 
        #s(VG_CA2, bs = 'tp', k = 10) +
        s(Population_Viral_Richness, bs = 'tp', k = 9) + 
        s(Region_Host_Richness, bs = 'tp', k = 3) +
        s(Host_Order, bs = 're'),
    data = virus_list,
    family = "binomial",
    select = FALSE
)


de = c(fit %>% dev_expl() - update(fit, . ~ . -Genome_Type) %>% dev_expl(),
#fit %>% dev_expl() - update(fit, . ~ . -s(TT_CA1, bs = 'tp', k = 10)) %>% dev_expl(),
fit %>% dev_expl() - update(fit, . ~ . -s(TT_CA2, bs = 'tp', k = 10)) %>% dev_expl(),
fit %>% dev_expl() - update(fit, . ~ . -s(Viral_Family, bs = 're', xt = list(S = list(virusFamPrecisionMat), rank = virusFamRank))) %>% dev_expl(),
fit %>% dev_expl() - update(fit, . ~ . -s(VG_CA1, bs = 'tp', k = 10)) %>% dev_expl(),
#fit %>% dev_expl() - update(fit, . ~ . -s(VG_CA2, bs = 'tp', k = 10)) %>% dev_expl(),
fit %>% dev_expl() - update(fit, . ~ . -s(Population_Viral_Richness, bs = 'tp', k = 9)) %>% dev_expl(),
fit %>% dev_expl() - update(fit, . ~ . -s(Region_Host_Richness, bs = 'tp', k = 3)) %>% dev_expl(),
fit %>% dev_expl() - update(fit, . ~ . -s(Host_Order, bs = 're')) %>% dev_expl())

x = predict(fit, terms="s(Viral_Family)", se.fit = TRUE, type="iterms")
tibble(effect_size = x$fit, se = x$se.fit, term = fit$model$Viral_Family) %>% distinct()


