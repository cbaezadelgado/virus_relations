## 1. VISUALIZING GAM PREDICTIONS

# model= my_gam_model
# original dataset = my_data
# Predictions in column "fit"


## 1. OBTAIN PREDICTIONS FROM MODEL
my_pred <- as.data.frame(predict(my_gam_model, type="response", se.fit=TRUE))

my_pred2 <- cbind(my_data, my_pred)

my_pred3 = cbind(my_pred2, CI_lower = my_pred2$fit-(1.96*my_pred2$se.fit), CI_upper = my_pred2$fit+(1.96*my_pred2$se.fit)) # calculating CI 95%


## 2. SETTING CLASSES
my_pred3$ds1_ss0 <- as.factor(my_pred3$ds1_ss0)
my_pred3$Enveloped <- as.factor(my_pred3$Enveloped)
my_pred3$Segmented <- as.factor(my_pred3$Segmented)
my_pred3$Nucleus0Cytplasm1 <- as.factor(my_pred3$Nucleus0Cytplasm1)
my_pred3$DNA0RNA1 <- as.factor(my_pred3$DNA0RNA1)
my_pred3$V_family <- as.factor(my_pred3$V_family)
my_pred3$Transm_vector <- as.factor(my_pred3$Transm_vector)
my_pred3$Virus_uses_multiple_prot_receptors <- as.factor(my_pred3$Virus_uses_multiple_prot_receptors)
my_pred3$Virus_uses_alternative_receptors <- as.factor(my_pred3$Virus_uses_alternative_receptors)
my_pred3$Virus_uses_accessory_receptors <- as.factor(my_pred3$Virus_uses_accessory_receptors)
my_pred3$Virus_uses_moieties <- as.factor(my_pred3$Virus_uses_moieties)


## 3. PLOT FIGURES
library(ggplot2)

## Figure 2
p1 <- my_pred3 %>%
  ggplot(aes(x=Enveloped, y=fit, fill=Enveloped, color=Enveloped)) +
  scale_color_manual(values=c("#F55050", "#86A3B8")) +
  scale_fill_manual(values=c("#F55050", "#86A3B8")) +
  geom_jitter(position=position_jitter(0.2), alpha=0.6, size=3) +
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="black", size=0.4) +
  theme_classic() +
  xlab("Enveloped") +
  theme(axis.text.y =element_text(size=10), axis.text.x =element_text(size=10), legend.position = "none")

p2 <- my_pred3 %>% ggplot(aes(Virion_diameter, fit)) +
  geom_smooth(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, fill="darkgrey") +
  geom_jitter(alpha=0.2, size=3, color="darkgrey") +
  xlab("Virion diameter") +
    theme_classic()

p3 <- my_pred3 %>%
  ggplot( aes(x=Transm_vector, y=fit, fill=Transm_vector, color=Transm_vector)) +
  scale_color_manual(values=c("#F55050", "#86A3B8")) +
  scale_fill_manual(values=c("#F55050", "#86A3B8")) +
  geom_jitter(position=position_jitter(0.2), alpha=0.6, size=3) +
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="black", size=0.4) +
  theme_classic() +
  xlab("Vector") +
  theme(axis.text.y =element_text(size=10), axis.text.x =element_text(size=10), legend.position = "none")


## Figure 6
p1 <- my_pred3 %>% ggplot(aes(N_sufficient_prot_receptors_per_virus, fit)) +
  geom_smooth(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, fill="darkgrey") +
  geom_jitter(aes(alpha=0.3, size=2, col=Enveloped)) +
  scale_color_manual(values=c("#F55050", "#86A3B8")) +
  scale_fill_manual(values=c("#F55050", "#86A3B8")) +
  xlab("N sufficient receptors/virus") +
    theme_classic() +
  theme(legend.position="none")

p2 <- my_pred3 %>% ggplot(aes(N_accessory_per_virus, fit)) +
  geom_smooth(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, fill="darkgrey") +
  geom_jitter(aes(alpha=0.3, size=2, col=Enveloped)) +
  scale_color_manual(values=c("#F55050", "#86A3B8")) +
  scale_fill_manual(values=c("#F55050", "#86A3B8")) +
  xlab("N accessory receptors/virus") +
  theme_classic() +
  theme(legend.position="none")

p3 <- my_pred3 %>% ggplot(aes(N_moieties_per_virus, fit)) +
  geom_smooth(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, fill="darkgrey") +
  geom_jitter(aes(alpha=0.3, size=2, col=Enveloped)) +
  scale_color_manual(values=c("#F55050", "#86A3B8")) +
  scale_fill_manual(values=c("#F55050", "#86A3B8")) +
  xlab("N moiety receptors/virus") +
  theme_classic() +
  theme(legend.position="none")







