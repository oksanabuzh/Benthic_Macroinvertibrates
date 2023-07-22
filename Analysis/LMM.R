ggplot(Benthic, aes(#log(Fish_abundance+9),
                    (Fish_abundance), 
                    Abund)) + geom_point()


#Waste.ef
#Gravel_ef
#Agric_m_log
#Highway_m_log


m1 <- lmer(log(Abund) ~    
             Waste.ef +
             Gravel_ef +
             Agric_m_log +
             Highway_m_log +
             (1|System_id:Year:Country),
           #  (1|Random_effects), 
           data = Benthic)
plot(m1)
check_heteroscedasticity(m1)

car::vif(m1) 
performance::check_collinearity(m1)

car::Anova(m1)


plot_model(m1, type = "pred", terms = c( "Highway_m_log"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway_m_log", "Abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Highway_m_log, y = Abund),  
             fill= "#0072B2", size=2, shape=21)


m1 <- lmer(log(Abund) ~    
             # System_type +
             Waste.ef +Gravel_ef +Agric_m_log +Highway_m_log +
             #  HI_PC1 + HI_PC2 + HI_PC3 +
              #  Dim.1 +Dim.2 + Dim.3 +
              Water_PC1 + Water_PC2 + Water_PC3 +        
             #   sqrt(Macrophyte_biomass) + 
            # Fish_abundance +
                  (1|System_id:Year:Country),
                #  (1|Random_effects), 
                data = Benthic)

car::Anova(m1)

plot_model(m1, type = "pred", terms = c( "Water_PC3"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC3", "Abund")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC3, y = Abund),  
             fill= "#0072B2", size=2, shape=21)

plot_model(m1, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Abund")) +
  geom_point(data = Benth, 
             mapping = aes(x = Macrophyte_biomass, y = Abund),  
             fill= "#0072B2", size=2, shape=21)

library(emmeans)
library(multcomp)

emmeans(m1, list(pairwise ~ Waste.ef), adjust = "tukey", 
        Letters = letters)

cld(emmeans(m1, list(pairwise ~ Waste.ef)),  
    #  type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benthic, aes(x = Waste.ef, y = log(Abund))) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Waste.ef", y="Abundance")+
  theme_bw()


####

# SR ----
m_SR1 <- glmer(SR ~ #Abund +  
                 #  System_type +
                 Waste.ef +
                 Gravel_ef +
                 Agric_m_log +
                 Highway_m_log +
                 (1|System_id:Year:Country), #  (1|Random_effects), 
               data = Benthic,
               family="poisson")

car::Anova(m_SR1)



performance::check_convergence(m_SR1)
performance::check_singularity(m_SR1)

car::Anova(m_SR1)
car::vif(m_SR1)

m_SR1 <- glmer(SR ~ #Abund +  
                 #  System_type +
                 #  HI_PC1 + HI_PC2 + HI_PC3 +
               #     Dim.1 +Dim.2 + Dim.3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +  
                 # Natural_areas + 
                 # sqrt(Macrophyte_biomass) + 
                 # Fish_abundance +
                 (1|System_id:Year:Country), #  (1|Random_effects), 
               data = Benthic,
               family="poisson")


performance::check_convergence(m_SR1)
performance::check_singularity(m_SR1)

car::Anova(m_SR1)



plot_model(m_SR1, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Species richness")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Macrophyte_biomass, y = SR),  
             fill= "#0072B2", size=2, shape=21)


###
# EVe----
m_FEve <- lmer(log(FEve) ~  # System_type +
                 Waste.ef +
                 Gravel_ef +
                 Agric_m_log +
                 Highway_m_log +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic )

plot(m_FEve)
check_heteroscedasticity(m_FEve) # no heteroscedasticity
check_outliers(m_FEve)

car::vif(m_FEve) 
performance::check_collinearity(m_FEve)
car::Anova(m_FEve)







Benthic$Fish_log <-  log(Benthic$Fish_abundance+9)

m_FEve <- lmer((FEve) ~  # System_type +
                #  HI_PC1 + HI_PC2 + HI_PC3 +
                #    Dim.1 +Dim.2 + Dim.3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                 # Natural_areas + 
                 # sqrt(Macrophyte_biomass) + 
               Fish_log +
                 (1|System_id:Year:System_type), # (1|Random_effects), 
               data = Benthic)

car::Anova(m_FEve)
summary(m_FEve)

car::vif(m_FEve) 
performance::check_collinearity(m_FEve)



plot_model(m_FEve, type = "pred", terms = c( "Highway_m_log"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway_m_log", "FEve")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Highway_m_log, y = FEve),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m_FEve, type = "pred", terms = c( "Water_PC1 "),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC1", "FEve")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC1 , y = FEve),  
             fill= "#0072B2", size=2, shape=21)

plot_model(m_FEve, type = "pred", terms = c( "Fish_abundance "),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Fish abundance ", "FEve")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Fish_abundance , y = FEve),  
             fill= "#0072B2", size=2, shape=21)



plot_model(m_FEve, type = "pred", terms = c( "Fish_log "),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Fish abundance ", "FEve")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Fish_log , y = FEve),  
             fill= "#0072B2", size=2, shape=21)










#FDiv-----
m_FDiv <- lmer(FDiv ~  # System_type +
                 Waste.ef +
                 Gravel_ef +
                 Agric_m_log +
                 Highway_m_log +
                # HI_PC1 + HI_PC2 + HI_PC3 +
                # Water_PC1 +  Water_PC2 + Water_PC3 +        
                 # Natural_areas + 
                # sqrt(Macrophyte_biomass) + 
               #  Fish_abundance  +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic)


car::Anova(m_FDiv)


m_FDiv <- lmer(FDiv ~  FEve+
                 # System_type +
                #  Waste.ef +Gravel_ef + #Agric_m_log + Highway_m_log +
                 #   HI_PC1 + HI_PC2 + HI_PC3 +
               Water_PC1 + Water_PC2 +  Water_PC3 +    
                 # Natural_areas + 
                 sqrt(Macrophyte_biomass) + 
                 Fish_abundance +
                # poly( Fish_abundance,2)  +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic)

car::Anova(m_FDiv)


performance::check_convergence(m_FDiv)

### singularity fit
performance::check_singularity(m_FDiv)

### multicolinearity
car::vif(m_FDiv) 
performance::check_collinearity(m_FDiv)

summary(m_FDiv)
car::Anova(m_FDiv)

plot_model(m_FDiv, type = "pred", terms = c( "Fish_abundance"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Fish_abundance", "FDiv")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Fish_abundance, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)


# FRic ----
m_FRic <- lmer(log(FRic) ~  # System_type +
                 Waste.ef +
                 Gravel_ef +
                 Agric_m_log +
                 Highway_m_log +
              #   HI_PC1 + HI_PC2 + HI_PC3 +
               #  Water_PC1 + Water_PC2 + Water_PC3 +        
                 # Natural_areas + 
                # sqrt(Macrophyte_biomass) + 
               #  Fish_abundance  +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic %>% 
                filter(!Sampling_site=="III1_2016")#, REML = F
              )

car::Anova(m_FRic)
plot(m_FRic)
check_heteroscedasticity(m_FRic) # no heteroscedasticity
check_outliers(m_FRic)



m_FRic <- lmer(FRic_log ~  poly(SR,2) +
                
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic %>% 
               filter(!Sampling_site=="III1_2016"))

car::Anova(m_FRic)

car::vif(m_FRic) 
performance::check_collinearity(m_FRic)



Benthic$FRic_log <- log(Benthic$FRic)


plot_model(m_FRic, type = "pred", terms = c( "SR"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("SR", "FRic")) +
  geom_point(data = Benthic %>% 
               filter(!Sampling_site=="III1_2016"), 
             mapping = aes(x = SR, y = FRic_log),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m_FDis, type = "pred", terms = c( "Water_PC2"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC2", "FDis")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC2, y = FDis),  
             fill= "#0072B2", size=2, shape=21)





m_FRic <- lmer((FEve) ~  poly(SR,2) +
                 
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic %>% 
                 filter(!Sampling_site=="III1_2016"))

car::Anova(m_FRic)


Benthic$FEve_log <- log(Benthic$FEve)


plot_model(m_FRic, type = "pred", terms = c( "SR"),  show.data=T, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("SR", "FEve")) +
  geom_point(data = Benthic ,
             mapping = aes(x = SR, y = FEve),  
             fill= "#0072B2", size=2, shape=21)


##


m_FRic <- lmer(FDiv ~  poly(SR,2) +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic %>% 
                 filter(!Sampling_site=="III1_2016"))

car::Anova(m_FRic)



plot_model(m_FRic, type = "pred", terms = c( "SR"),  show.data=T, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("SR", "FDiv")) +
  geom_point(data = Benthic ,
             mapping = aes(x = SR, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)





m_FRic <- lmer(FDis ~  
                 # poly(SR,2) +
                 SR +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic %>% 
                 filter(!Sampling_site=="III1_2016"))

car::Anova(m_FRic)



plot_model(m_FRic, type = "pred", terms = c( "SR"),  show.data=T, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("SR", "FDis")) +
  geom_point(data = Benthic ,
             mapping = aes(x = SR, y = FDis),  
             fill= "#0072B2", size=2, shape=21)




#####




m_FRic <- lmer(FDiv ~ #  poly(FEve,2) +
                 FEve +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benthic %>% 
                 filter(!Sampling_site=="III1_2016"))

car::Anova(m_FRic)



plot_model(m_FRic, type = "pred", terms = c( "FEve"),  show.data=T, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("FEve", "FDiv")) +
  geom_point(data = Benthic ,
             mapping = aes(x = FEve, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)
