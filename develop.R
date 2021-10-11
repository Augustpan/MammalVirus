
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

df = read_csv("data/host_occurence.csv")
meta = read_csv("data/host_metadata.csv")
mf = read_csv("data/Mammal_family.csv")

ndf = df %>% 
    left_join(meta) %>%
    left_join(mf)
ca = cca(ndf[,2:5])
sp_score = ca$CA$u[,1]
st_score = ca$CA$v[,1]

sp_ord = order(sp_score)
st_ord = order(st_score)

x = cbind(
    select(ndf, Mammal, Mammal_family, Mammal_genus, Mammal_species)[sp_ord,],
    ndf[,2:5][sp_ord, st_ord]
)


vir = read_tsv("data/VIR.tsv")
vcf = read_csv("variant_calls.csv")
df = filter(vcf, grepl(vir$VIR[1], VIR))


x = tb_virus_list %>%
    select(Host_Name, Virus_Name, Total_Abundance) %>%
    pivot_wider(names_from=Virus_Name, values_from=Total_Abundance, values_fill=0) %>%
    mutate(across(2:ncol(.), ~ as.integer(.x > 0))) %>%
    as.data.frame()
rownames(x) = x$Host_Name
x$Host_Name = NULL

tmp = str_replace(colnames(x), "\\(.+\\)", "")

colnames(x) = str_split(tmp, "\\s", simplify = T) %>% 
    plyr::aaply(1, function(x){paste(x, collapse ="_")}) %>%
    str_replace("_{2,}", "") %>%
    str_to_lower() %>%
    str_replace("_&.+", "")

for (tree_file in list.files("data/virus_tree")[-4]) {
    #tree_file = "data/virus_tree/5_my_Rhabdo_L_protein_alin PhyML Tree.nwk"
    print(tree_file)
    tre = phyloseq::read_tree(paste0("data/virus_tree/", tree_file))
    tre$tip.label = tre$tip.label %>%
        str_replace("_partial", "") %>%
        str_replace("_complete_genome", "") %>%
        str_replace("Niviventer_coninga", "Leopoldamys_edwardsi") %>%
        str_to_lower() %>%
        str_replace("jingshan", "jingmen")
    unk = setdiff(tre$tip.label, colnames(x))
    xs = select(x, all_of(intersect(colnames(x), tre$tip.label))) %>%
        filter(rowSums(.)>0)
    
    ot = otu_table(xs, taxa_are_rows = F)
    ps = phyloseq(ot, tre)
    dst = UniFrac(ps, normalized = T) %>%
        as.matrix() %>%
        as.data.frame() %>%
        mutate(s1 = rownames(.)) %>%
        pivot_longer(cols=1:ncol(.)-1, names_to="s2", values_to="uw_unifrac") %>%
        mutate(uw_unifrac = uw_unifrac * length(tre$edge))
}


