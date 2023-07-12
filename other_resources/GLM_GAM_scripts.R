
## GLM (Fig.1)
my_GLM_model <- glm(N_all_recept ~ log10_Nt_seqs, family= poisson(link=log), data=my_data)


## GAM models
library(mgcv)

## Fig. 2
# Predicting number of total receptors (protein and moieties) per virus
my_GAM_model <- mgcv::gam(N_receptors_and_moieties_per_virus ~ s(log10p1_Nseqs) + Enveloped +
                             s(Virion_diameter) + Transm_vector + s(V_family, bs="re"), 
                         data=my_data, family="poisson", method ="REML")
# Predicting number of protein receptors per virus
my_GAM_model <- mgcv::gam(N_prot_receptors_per_virus ~ s(log10p1_Nseqs) + Enveloped +
                             s(Virion_diameter) + Transm_vector + s(V_family, bs="re"), 
                         data=my_data, family="poisson", method ="REML")
# Predicting number of moieties per virus
my_GAM_model <- mgcv::gam(N_moieties_per_virus ~ s(log10p1_Nseqs) + Enveloped +
                             s(Virion_diameter) + Transm_vector + s(V_family, bs="re"), 
                         data=my_data, family="poisson", method ="REML")
# Predicting number of sufficient receptors per virus
my_GAM_model <- mgcv::gam(N_sufficient_prot_receptors_per_virus ~ s(log10p1_Nseqs) + Enveloped +
                             s(Virion_diameter) + Transm_vector + s(V_family, bs="re"), 
                         data=my_data, family="poisson", method ="REML")
# Predicting number of accessory receptors per virus
my_GAM_model <- mgcv::gam(N_accessory_per_virus ~ s(log10p1_Nseqs) + Enveloped +
                             s(Virion_diameter) + Transm_vector + s(V_family, bs="re"), 
                         data=my_data, family="poisson", method ="REML")

## Fig. 6
# Predicting number of host species (excluding humans)
my_GAM_model <- mgcv::gam(Nspecies_excl_humans ~ s(log10p1_Nseqs) + ds1_ss0 + Enveloped + Segmented +
Nucleus0Cytplasm1 + DNA0RNA1 + s(log10_genome_size) + s(Genome_GC.) + s(Virion_diameter) +                             Transm_vector + s(V_family, bs="re") + s(N_sufficient_prot_receptors_per_virus) + s(N_accessory_per_virus) + s(N_moieties_per_virus, k=3), data=my_data, family="poisson", method ="REML")

# Predicting number of host families
my_GAM_model <- mgcv::gam(Nfams ~ s(log10p1_Nseqs) + ds1_ss0 + Enveloped + Segmented +
Nucleus0Cytplasm1 + DNA0RNA1 + s(log10_genome_size) + s(Genome_GC.) + s(Virion_diameter) +                             Transm_vector + s(V_family, bs="re") + s(N_sufficient_prot_receptors_per_virus) + s(N_accessory_per_virus) + s(N_moieties_per_virus, k=3), data=my_data, family="poisson", method ="REML")


