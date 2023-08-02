# Mixed effect models (G)LMM

# Response variables: abundance, species richness, functional diversity

# Predictors

# Model 1: system type (pond, lake, channel) and human impact intensity 
## Waste.ef
## Gravel_ef
## Agric_m_log
## Highway_m_log

# Model 2: water properties, macrophyte biomass, fish predation

# packages----

library(tidyverse)
library(ggplot2)
library(car)
library(lme4) 
library(lmerTest)
library(performance)
library(sjPlot)
library(emmeans)
library(multcomp)

# set theme for the plots

set_theme(base = theme_bw(),
          axis.textsize.x = 0.9, axis.textsize.y = 0.9, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1,
          legend.pos = "top")


#Data----

Benth_dat <- read_csv ("Data/Benth_PCA2.csv")
str(Benth_dat)


FuncDiv <- read.csv ( "Data/FuncDiv.csv", header=T)%>% 
  dplyr::select(-qual.FRic, -sing.sp, -nbsp) %>% 
  rename(Sampling_site=X)

names(FuncDiv)

Benthic <- Benth_dat %>%
  full_join(FuncDiv, by="Sampling_site") %>%
  mutate(System_id = as_factor(System_id), Year = as_factor(Year), 
       Country=as_factor(Country),
       Waste.ef=as_factor(Waste.effluent),
       Gravel_ef=as_factor(Gravel_eff),
       Agric_m_log=log(Agriculture_m),
       Highway_m_log=log(Highway_m),
       Macroph_sqrt=sqrt(Macrophyte_biomass))

str(Benthic)


  

#--------------------------#
  
# 1) Abundance----

#--------------------------#
# Explore the relationships
ggplot(Benthic, aes(x = System_type, y = log(Abund))) + geom_boxplot()
ggplot(Benthic, aes(sqrt(Macrophyte_biomass), (Abund))) + geom_point()
ggplot(Benthic, aes(sqrt(Fish_abundance), (Abund))) + geom_point()


Benthic$Abund

### Mod1: HI----

m1 <- lmer(Abund ~  System_type +
              Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
                  (1|System_id:Year:Country),
                data = Benthic)
## assumptions 
plot(m1) 
check_heteroscedasticity(m1)
qqnorm(resid(m1))
qqline(resid(m1))

### ->  transform the response
Benthic$Abund_log <- log(Benthic$Abund)

m1.1 <- lmer(Abund_log ~  System_type +
             Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
             (1|System_id:Year:Country),
           data = Benthic)
plot(m1.1) 
check_heteroscedasticity(m1.1)
qqnorm(resid(m1.1))
qqline(resid(m1.1))

### singularity fit
performance::check_singularity(m1.1)
performance::check_convergence(m1.1)

### multicolinearity
car::vif(m1.1) 
performance::check_collinearity(m1.1)

summary(m1.1)
car::Anova(m1.1)

##### -> remove System_type----
# model selection:
m1.1 <- lmer(Abund_log ~  System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m1.2 <- lmer(Abund_log ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F)


library(MuMIn)
AICc(m1.1, m1.2)

anova(m1.1, m1.2)
# We select model m1.2 (without "System type")

car::Anova(m1.2)
# drop1(m1.2, test="Chi")

write.csv(car::Anova(m1.2),  file = "Results/Abund_LMM.csv")
# write.csv(drop1(m1.2, test="Chi"),  file = "Results/Abund_LMM.csv")

# write.csv(coef(summary(m1.2)),  file = "Results/Abund_summary.csv")


#### -> plot----
plot_model(m1.2, type = "pred", terms = c( "Highway_m_log"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway distance", "Abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Highway_m_log, y = Abund_log),  
             fill= "#0072B2", size=2, shape=21)

######  plot on a raw scale----

m1.2 <- lmer(log(Abund) ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
               data = Benthic, REML=F)

plot_model(m1.2, type = "pred", terms = c( "Highway_m_log"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway distance", "Abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Highway_m_log, y = Abund),  
             fill= "#0072B2", size=2, shape=21)

#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m1.2)
# Partial R2 for fixed effects

m1.2 <- lmer(Abund_log ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
              # (1|System_id:Year:Country),
               (1|Random_effects), # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
             data = Benthic, REML=F)
#Anova(m1.2)
R2_abund1 <- r2glmm::r2beta(m1.2,  partial = T, data = Benthic)
R2_abund1
plot(R2_abund1)

# save results
# write.csv(Anova(m1.2),  file = "Results/Abund_LMM.csv")
# write.csv(coef(summary(m1.2)),  file = "Results/Abund_summary.csv")


### Mod2: Water properties----

m1.3 <- lmer(Abund_log ~ 
               System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               sqrt(Macrophyte_biomass) + 
               Fish_abundance +
              (1|System_id:Year:Country), data = Benthic)


plot(m1.3) 
check_heteroscedasticity(m1.3)
qqnorm(resid(m1.3))
qqline(resid(m1.3))

### multicolinearity
performance::check_collinearity(m1.3)

summary(m1.3)
# car::Anova(m1.3)

##### -> remove System_type----
# model selection:

m1.3 <- lmer(Abund_log ~ 
               System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m1.4 <- lmer(Abund_log ~ 
              # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F)

library(MuMIn)
AICc(m1.3, m1.4)

anova(m1.3, m1.4)

# We select model m1.4 (without "System type")

car::Anova(m1.4)

#### -> plot---- 

plot_model(m1.4, type = "pred", terms = c( "Water_PC3"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water PC3", "Abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC3, y = Abund_log),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m1.4, type = "pred", terms = c( "Macroph_sqrt"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Macroph_sqrt, y = Abund_log),  
             fill= "#0072B2", size=2, shape=21)


######  plot on a raw scale----

m1.4 <- lmer(log(Abund) ~ 
               # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               sqrt(Macrophyte_biomass) + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F)


plot_model(m1.4, type = "pred", terms = c( "Water_PC3"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water PC3", "Abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC3, y = Abund),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m1.4, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Macrophyte_biomass, y = Abund),  
             fill= "#0072B2", size=2, shape=21)


#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m1.4)
# Partial R2 for fixed effects
# # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
m1.4 <- lmer(log(Abund) ~ 
               # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               sqrt(Macrophyte_biomass) + Fish_abundance +
              # (1|System_id:Year:Country),
             (1|Random_effects), data = Benthic, REML=F)

car::Anova(m1.4)
R2_abund2 <- r2glmm::r2beta(m1.4,  partial = T, data = Benthic)
R2_abund2
plot(R2_abund2)


### save partial R2 ----

R2_Abund <- rbind(as.data.frame(R2_abund1), as.data.frame(R2_abund2)) %>% 
  rename(Variables=Effect) %>% 
  dplyr::select(Variables, Rsq) %>% 
  filter(!Variables=="Model") %>% 
  mutate(Variables=recode_factor(Variables, 
         Highway_m_log="Highway_m", Agric_m_log ="Agric_m", 
         Waste.ef1="Waste", Gravel_eff="Gravel",
         "sqrt(Macrophyte_biomass)"="Macrophytes", Fish_abundance ="Fish")) %>% 
  mutate(Variables=fct_relevel(Variables, c("Agric_m", "Highway_m", "Gravel", "Waste", 
                                            "Water_PC1", "Water_PC2", "Water_PC3",
                                            "Macrophytes", "Fish"))) %>% 
  arrange(Variables)


R2_Abund
write.csv(R2_Abund, "Results/Abund_R2.csv")

#--------------------------#

# 2) SR - species richness ----
  
#--------------------------#



### Mod1: HI----

m2.1 <- glmer(SR ~ System_type +
                Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
                (1|System_id:Year:Country), data = Benthic, family="poisson")

performance::check_convergence(m2.1)
# Model failed to converge !

### multicolinearity
performance::check_collinearity(m2.1)

# remove System_type due to collinearity
plot(m2.1)

##### -> remove System_type----

m2.2 <- glmer(SR ~ # System_type +
                Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
                (1|System_id:Year:Country), data = Benthic, family="poisson")

performance::check_convergence(m2.2)
# model converged

### multicolinearity
performance::check_collinearity(m2.2)


# check overdispersion 
sum(residuals(m2.2, type = "pearson")^2) / df.residual(m2.2)

## assumptions 
plot(m2.2) 
qqnorm(resid(m2.2))
qqline(resid(m2.2))
check_outliers(m2.2)


anova(m2.1, m2.2)

summary(m2.2)
car::Anova(m2.2)

#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m2.2)
# Partial R2 for fixed effects
# # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 

m2.2 <- glmer(SR ~ Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
                (1|Random_effects), data = Benthic, family="poisson")

car::Anova(m2.2)
R2_SR1 <- r2glmm::r2beta(m2.2,  partial = T, data = Benthic, method = 'sgv')
R2_SR1
plot(R2_SR1)


### Mod2: Water properties----

m2.3 <- glmer(SR ~ System_type +
                Water_PC1 + Water_PC2 + Water_PC3 +
                Macroph_sqrt + 
                Fish_abundance +
                (1|System_id:Year:Country), data = Benthic, family="poisson")

performance::check_convergence(m2.3)
# model converged

plot(m2.3) 
qqnorm(resid(m2.3))
qqline(resid(m2.3))

### multicolinearity
performance::check_collinearity(m2.3)

summary(m2.3)
car::Anova(m2.3)

##### -> remove System_type----
# model selection:

m2.4 <- glmer(SR ~ # System_type +
                Water_PC1 + Water_PC2 + Water_PC3 +
                Macroph_sqrt + 
                Fish_abundance +
                (1|System_id:Year:Country), data = Benthic, family="poisson")

performance::check_convergence(m2.4)

library(MuMIn)
AICc(m2.3, m2.4)
anova(m2.3, m2.4)

# Model m2.4 did not differ significantly from model 2.3, but "System type" is significant, and the likelihood test is marginally significant, thus we select for model 2.3 

car::Anova(m2.3)

#### -> plot---- 
plot_model(update(m2.3, . ~ .-System_type), # for the plot we need to remove the categorical variable
           type = "pred", terms = c( "Macroph_sqrt"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "SR")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Macroph_sqrt, y = SR),  
             fill= "#0072B2", size=2, shape=21)

######  plot on a raw scale----

m2.3 <- glmer(SR ~  System_type +
                Water_PC1 + Water_PC2 + Water_PC3 +
                sqrt(Macrophyte_biomass) + 
                Fish_abundance +
                (1|System_id:Year:Country), data = Benthic, family="poisson")

plot_model(update(m2.3, . ~ .-System_type), # for the plot we need to remove the categorical variable
           type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "SR")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Macrophyte_biomass, y = SR),  
             fill= "#0072B2", size=2, shape=21)

# plot system type:
### Pairwise comparisons for System_type----

marg_means <- emmeans(m2.3, list(pairwise ~ System_type), type="response")
marg_means
cld(marg_means,Letters = letters, adjust = "tukey")

ggplot(Benthic, aes(x = System_type, y = SR)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Species richness")+
  theme_bw()

#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m2.4)
# Partial R2 for fixed effects
# # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 

m2.3 <- glmer(SR ~  System_type +
                Water_PC1 + Water_PC2 + Water_PC3 +
                Macroph_sqrt + Fish_abundance +
                (1|Random_effects), data = Benthic, family="poisson")

car::Anova(m2.3)
R2_SR2 <- r2glmm::r2beta(m2.3,  partial = T, data = Benthic, method = 'sgv')
R2_SR2
plot(R2_SR2)


### save partial R2 ----

R2_SR <- rbind(as.data.frame(R2_SR1), as.data.frame(R2_SR2)) %>% 
  rename(Variables=Effect) %>% 
  dplyr::select(Variables, Rsq) %>% 
  filter(!Variables=="Model", 
         !Variables=="System_typePond", !Variables=="System_typeLake") %>% 
  mutate(Variables=recode_factor(Variables, 
                                 Highway_m_log="Highway_m", Agric_m_log ="Agric_m", 
                                 Waste.ef1="Waste", Gravel_eff="Gravel",
                                 Macroph_sqrt="Macrophytes", Fish_abundance ="Fish")) %>% 
  mutate(Variables=fct_relevel(Variables, c("Agric_m", "Highway_m", "Gravel", "Waste", 
                                            "Water_PC1", "Water_PC2", "Water_PC3",
                                            "Macrophytes", "Fish"))) %>% 
  arrange(Variables)


R2_SR
write.csv(R2_SR, "Results/SR_R2.csv")

# 3) Functional diversity measures----

            
# FEve              
# FDiv              
# FDis             

# FRic and RaoQ are not considered, 
# because FRic correlates with SR, while RaoQ correlates with FDis



#--------------------------#

## 3.1) FEve ----

#--------------------------#

m3.1 <- lmer(FEve ~  System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic)

plot(m3.1) 
check_outliers(m3.1)
qqnorm(resid(m3.1))
qqline(resid(m3.1))

### singularity fit
performance::check_singularity(m3.1)
performance::check_convergence(m3.1)

### multicolinearity
performance::check_collinearity(m3.1)

summary(m3.1)
car::Anova(m3.1)

##### -> remove System_type----
# model selection:
m3.1 <- lmer(FEve ~  System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m3.2 <- lmer(FEve ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F)


library(MuMIn)
AICc(m3.1, m3.2)

anova(m3.1, m3.2)

# We select model m3.2 (without "System type")

car::Anova(m3.2)

# 
plot_model(m3.2, type = "pred", terms = c( "Highway_m_log"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway distance", "FEve")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Highway_m_log, y = FEve),  
             fill= "#0072B2", size=2, shape=21)


#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m3.2)
# Partial R2 for fixed effects

m3.2 <- lmer(FEve ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               # (1|System_id:Year:Country),
               (1|Random_effects), # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
             data = Benthic)

car::Anova(m3.2)
R2_FEve1 <- r2glmm::r2beta(m3.2,  partial = T, data = Benthic)
R2_FEve1
plot(R2_FEve1)


### Mod2: Water properties----

m3.3 <- lmer(FEve ~ System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic)


plot(m3.3) 
qqnorm(resid(m3.3))
qqline(resid(m3.3))

### multicolinearity
performance::check_collinearity(m3.3)

summary(m3.3)
car::Anova(m3.3)

##### -> remove System_type----
# model selection:

m3.3 <- lmer(FEve ~ System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m3.4 <- lmer(FEve ~ # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F)

library(MuMIn)
AICc(m3.3, m3.4)
anova(m3.3, m3.4)
# We select model m3.4 (without "System type")

car::Anova(m3.4)

#### -> plot---- 

plot_model(m3.4, type = "pred", terms = c( "Water_PC1"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water PC1", "FEve")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC1, y = FEve),  
             fill= "#0072B2", size=2, shape=21)




#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m3.4)
# Partial R2 for fixed effects
# # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
m3.4 <- lmer(FEve ~ # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               sqrt(Macrophyte_biomass) + Fish_abundance +
               # (1|System_id:Year:Country),
               (1|Random_effects), data = Benthic, REML=F)


car::Anova(m3.4)
R2_FEve2 <- r2glmm::r2beta(m3.4,  partial = T, data = Benthic)
R2_FEve2
plot(R2_FEve2)


### save partial R2 ----

R2_FEve <- rbind(as.data.frame(R2_FEve1), as.data.frame(R2_FEve2)) %>% 
  rename(Variables=Effect) %>% 
  dplyr::select(Variables, Rsq) %>% 
  filter(!Variables=="Model") %>% 
  mutate(Variables=recode_factor(Variables, 
                                 Highway_m_log="Highway_m", Agric_m_log ="Agric_m", 
                                 Waste.ef1="Waste", Gravel_eff="Gravel",
                                 "sqrt(Macrophyte_biomass)"="Macrophytes", Fish_abundance ="Fish")) %>% 
  mutate(Variables=fct_relevel(Variables, c("Agric_m", "Highway_m", "Gravel", "Waste", 
                                            "Water_PC1", "Water_PC2", "Water_PC3",
                                            "Macrophytes", "Fish"))) %>% 
  arrange(Variables)


R2_FEve
write.csv(R2_FEve, "Results/FEve_R2.csv")

#--------------------------#

## 3.2) FDiv ----

#--------------------------#

m4.1 <- lmer(FDiv ~  System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic)

plot(m4.1) 
check_outliers(m4.1)
qqnorm(resid(m4.1))
qqline(resid(m4.1))

### singularity fit
performance::check_singularity(m4.1)
performance::check_convergence(m4.1)

### multicolinearity
performance::check_collinearity(m4.1)

summary(m4.1)
car::Anova(m4.1)

##### -> remove System_type----
# model selection:
m4.1 <- lmer(FDiv ~  System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m4.2 <- lmer(FDiv ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F)


library(MuMIn)
AICc(m4.1, m4.2)

anova(m4.1, m4.2)

# We select model m4.2 (without "System type")

car::Anova(m4.2)



#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m4.2)
# Partial R2 for fixed effects

m4.2 <- lmer(FDiv ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               # (1|System_id:Year:Country),
               (1|Random_effects), # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
             data = Benthic)

car::Anova(m4.2)
R2_FDiv1 <- r2glmm::r2beta(m4.2,  partial = T, data = Benthic)
R2_FDiv1
plot(R2_FDiv1)


### Mod2: Water properties----

m4.3 <- lmer(FDiv ~ System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic)


plot(m4.3) 
qqnorm(resid(m4.3))
qqline(resid(m4.3))

### multicolinearity
performance::check_collinearity(m4.3)

summary(m4.3)
car::Anova(m4.3)

##### -> remove System_type----
# model selection:

m4.3 <- lmer(FDiv ~ System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m4.4 <- lmer(FDiv ~ # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F)


library(MuMIn)
AICc(m4.3, m4.4)
anova(m4.3, m4.4)
# We select model m4.4 (without "System type")


car::Anova(m4.4)

#### -> plot---- 

plot_model(m4.4, type = "pred", terms = c( "Water_PC1"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water PC1", "FDiv")) +
   geom_point(data = Benthic, 
             mapping = aes(x = Water_PC1, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m4.4, type = "pred", terms = c( "Water_PC2"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water PC2", "FDiv")) +
  xlim(-3,2)+ylim(0.54,1)+
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC2, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m4.4, type = "pred", terms = c( "Fish_abundance"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Fish abundance", "FDiv")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Fish_abundance, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)


#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m4.4)
# Partial R2 for fixed effects
# # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
m4.4 <- lmer(FDiv ~ # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               sqrt(Macrophyte_biomass) + Fish_abundance +
               # (1|System_id:Year:Country),
               (1|Random_effects), data = Benthic)


car::Anova(m4.4)
R2_FDiv2 <- r2glmm::r2beta(m4.4,  partial = T, data = Benthic)
R2_FDiv2
plot(R2_FDiv2)


### save partial R2 ----

R2_FDiv <- rbind(as.data.frame(R2_FDiv1), as.data.frame(R2_FDiv2)) %>% 
  rename(Variables=Effect) %>% 
  dplyr::select(Variables, Rsq) %>% 
  filter(!Variables=="Model") %>% 
  mutate(Variables=recode_factor(Variables, 
                                 Highway_m_log="Highway_m", Agric_m_log ="Agric_m", 
                                 Waste.ef1="Waste", Gravel_eff="Gravel",
                                 "sqrt(Macrophyte_biomass)"="Macrophytes", Fish_abundance ="Fish")) %>% 
  mutate(Variables=fct_relevel(Variables, c("Agric_m", "Highway_m", "Gravel", "Waste", 
                                            "Water_PC1", "Water_PC2", "Water_PC3",
                                            "Macrophytes", "Fish"))) %>% 
  arrange(Variables)


R2_FDiv
write.csv(R2_FDiv, "Results/FDiv_R2.csv")




#--------------------------#

## 3.3) FDis ----

#--------------------------#

m5.1 <- lmer(FDis ~  System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic)

plot(m5.1) 
check_homogeneity(m5.1)
check_outliers(m5.1)
qqnorm(resid(m5.1))
qqline(resid(m5.1))

### singularity fit
performance::check_singularity(m5.1)
performance::check_convergence(m5.1)

### multicolinearity
performance::check_collinearity(m5.1)

summary(m5.1)
car::Anova(m5.1)

##### -> remove System_type----
# model selection:
m5.1 <- lmer(FDis ~  System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m5.2 <- lmer(FDis ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               (1|System_id:Year:Country),
             data = Benthic, REML=F)


library(MuMIn)
AICc(m5.1, m5.2)

anova(m5.1, m5.2)

# We select model m5.2 (without "System type")

car::Anova(m5.2)

#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m5.2)
# Partial R2 for fixed effects

m5.2 <- lmer(FDis ~  # System_type +
               Waste.ef + Gravel_eff + Agric_m_log + Highway_m_log + 
               # (1|System_id:Year:Country),
               (1|Random_effects), # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
             data = Benthic)

car::Anova(m5.2)
R2_FDis1 <- r2glmm::r2beta(m5.2,  partial = T, data = Benthic)
R2_FDis1
plot(R2_FDis1)


### Mod2: Water properties----

m5.3 <- lmer(FDis ~ System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic)


plot(m5.3) 
qqnorm(resid(m5.3))
qqline(resid(m5.3))

### multicolinearity
performance::check_collinearity(m5.3)

summary(m5.3)
car::Anova(m5.3)

##### -> remove System_type----
# model selection:

m5.3 <- lmer(FDis ~ System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F) # use REML=F for comparison of fixed effects

m5.4 <- lmer(FDis ~ # System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=F)


library(MuMIn)
AICc(m5.3, m5.4)
anova(m5.3, m5.4)
# We select model m5.3 (with "System type") as the likelihood test is close to be significant


m5.3 <- lmer(FDis ~  System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macroph_sqrt + Fish_abundance +
               (1|System_id:Year:Country), data = Benthic, REML=T)

car::Anova(m5.3)

#### -> plot---- 

plot_model(update(m5.3, .~.-System_type), 
           type = "pred", terms = c( "Water_PC3"),  show.data=F, jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water PC3", "FDis")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC3, y = FDis),  
             fill= "#0072B2", size=2, shape=21)


#### -> R2 ----
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m5.4)
# Partial R2 for fixed effects
# # for r2beta the nested random effects do not work. Thus, I included manually assigned random effects, which are equal to "System_id:Year:Country" 
m5.4 <- lmer(FDis ~  System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               sqrt(Macrophyte_biomass) + Fish_abundance +
               # (1|System_id:Year:Country),
               (1|Random_effects), data = Benthic)


car::Anova(m5.4)
R2_FDis2 <- r2glmm::r2beta(m5.4,  partial = T, data = Benthic)
R2_FDis2
plot(R2_FDis2)

### save partial R2 ----

R2_FDis <- rbind(as.data.frame(R2_FDis1), as.data.frame(R2_FDis2)) %>% 
  rename(Variables=Effect) %>% 
  dplyr::select(Variables, Rsq) %>% 
  filter(!Variables=="Model",
         !Variables=="System_typePond", !Variables=="System_typeLake") %>% 
  mutate(Variables=recode_factor(Variables, 
                                 Highway_m_log="Highway_m", Agric_m_log ="Agric_m", 
                                 Waste.ef1="Waste", Gravel_eff="Gravel",
                                 "sqrt(Macrophyte_biomass)"="Macrophytes", Fish_abundance ="Fish")) %>% 
  mutate(Variables=fct_relevel(Variables, c("Agric_m", "Highway_m", "Gravel", "Waste", 
                                            "Water_PC1", "Water_PC2", "Water_PC3",
                                            "Macrophytes", "Fish"))) %>% 
  arrange(Variables)


R2_FDis
write.csv(R2_FDis, "Results/FDis_R2.csv")
