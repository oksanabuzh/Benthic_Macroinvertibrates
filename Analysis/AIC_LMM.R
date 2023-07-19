# Mixed effect models (G)LMM

# Response variables: abundance, species richness, functional diversity
# Predictors: System type (pond, lake, channel), human impact intensity (HI_PC1, HI_PC2, HI_PC3)

library(tidyverse)
library(ggplot2)
library(car)
library(lme4) 
library(lmerTest)
library(performance)
library(sjPlot)
library(emmeans)
library(multcomp)

#Data----

Benth_dat <- read_csv ("Data/Benth_PCA2.csv")
str(Benth_dat)


FuncDiv <- read.csv ( "Data/FuncDiv.csv", header=T)%>% 
  dplyr::select(-qual.FRic, -sing.sp, -nbsp) %>% 
  rename(Sampling_site=X)

names(FuncDiv)

Benth <- Benth_dat %>%
  full_join(FuncDiv, by="Sampling_site") %>% 
  mutate(System_id = as_factor(System_id), Year = as_factor(Year), Country=as_factor(Country)) 

str(Benth)

#--------------------------#

# Abundance----

#--------------------------#

# Explore the relationships
ggplot(Benth, aes(x = System_type, y = log(Abund))) + geom_boxplot()
ggplot(Benth, aes(sqrt(Macrophyte_biomass), log(Abund))) + geom_point()
ggplot(Benth, aes(Water_PC1, log(Abund))) + geom_point()
ggplot(Benth, aes(Water_PC3, log(Abund))) + geom_point()
ggplot(Benth, aes(Fish_abundance, log(Abund))) + geom_point()



# LMM: Abundance

m1 <- lmer(Abund ~    
             System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +        
             sqrt(Macrophyte_biomass) +  
             Fish_abundance +
             (1|Country/Year/System_id),
           #(1|Random_effects), 
           data = Benth)
## assumptions 
plot(m1) # heteroscedasticity  ->  transform response
qqnorm(resid(m1))
qqline(resid(m1))

### ->  transform the response

m2 <- lmer(log(Abund) ~    
             System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +        
             sqrt(Macrophyte_biomass) +  
             Fish_abundance +
             (1|System_id:Year:Country),
           #  (1|Random_effects), 
           data = Benth)

car::Anova(m2)
summary(m2)

## assumptions
#check_model(m2)

plot(m2) # homoscedasticity met
check_heteroscedasticity(m2)
qqnorm(resid(m2))
qqline(resid(m2))
check_outliers(m2)

### singularity fit
performance::check_singularity(m2)
performance::check_convergence(m2)

### multicolinearity
car::vif(m2) 
performance::check_collinearity(m2)


## Test random effects
# Models differ in fixed effects:  lmer (…, REML = F)
# Models differ in random effects: lmer (…, REML = T)
ranova(m2)


## Model estimates----

m2 <- lmer(log(Abund) ~    
             System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +        
             sqrt(Macrophyte_biomass) +  
             Fish_abundance +
             (1|System_id:Year:Country),
           #  (1|Random_effects), 
           data = Benth)

car::Anova(m2)

### Pairwise comparisons for System_type----

library(emmeans)
library(multcomp)

emmeans(m2, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters)

cld(emmeans(m2, list(pairwise ~ System_type)),  
    #  type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = log(Abund))) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Abundance")+
  theme_bw()



### check multicolinearity
car::vif(m2) 
performance::check_collinearity(m2) # System_type VIF is large  

# remove System_type

### AIC selection? ----
m_abund1 <- lmer(log(Abund) ~    
                  # System_type +
                  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  sqrt(Macrophyte_biomass) +  
                  Fish_abundance +
                  (1|System_id:Year:Country),
                #  (1|Random_effects), 
                data = Benth, REML = F)


m_abund2 <- lmer(log(Abund) ~    
                  # System_type +
                 # HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  sqrt(Macrophyte_biomass) +  
                  Fish_abundance +
                  (1|System_id:Year:Country),
                #  (1|Random_effects), 
                data = Benth, REML = F)

### AICc comparison----
library(MuMIn)
arrange(AICc(m_abund1, m_abund2), 
        +AICc)

car::Anova(m_abund2)

m_abund <- lmer(log(Abund) ~    
                  # System_type +
                  #  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  sqrt(Macrophyte_biomass) +  
                  #  Fish_abundance +
                  (1|System_id:Year:Country),
                #  (1|Random_effects), 
                data = Benth)

car::Anova(m_abund)

### Plot predictor effects ----

set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          legend.pos = "top")

plot_model(m_abund, type = "pred", terms = c( "Water_PC3"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water properties PC3", "Abundance")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC3, y = Abund),  
             fill= "#0072B2", size=2, shape=21)

# marginally significant

plot_model(m_abund, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Abundance")) +
  geom_point(data = Benth, 
             mapping = aes(x = Macrophyte_biomass, y = Abund),  
             fill= "#0072B2", size=2, shape=21)




# R2 for the entire model
# R2m and R2c are marginal (for fixed predictors) and 
## conditional (for fixed and random predictors) coefficients of determination
MuMIn::r.squaredGLMM(m_abund)

# Partial R2 for fixed effects

m_abund <- lmer(log(Abund) ~    
                  # System_type +
                  #  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  sqrt(Macrophyte_biomass) +  
                  # Fish_abundance +
                  #  (1|System_id:Year:Country),
                  (1|Random_effects), 
                data = Benth)

car::Anova(m_abund)
R2 <- r2glmm::r2beta(m_abund,  partial = T, data = Benth)
R2
plot(R2)

# save results

# write.csv(MuMIn::r.squaredGLMM(FE_m4),  file = "Results/Mod_R2_Abund.csv")
# write.csv(R2,  file = "Results/R2_part_Abund.csv")

#write.csv(Anova(m_abund),  file = "Results/LMM_Abund.csv")
#write.csv(coef(summary(m_abund)),  file = "Results/summary_Abund.csv")


#--------------------------#

# SR - species richness ----

#--------------------------#
m3 <- glmer(SR ~ System_type +
              HI_PC1 + HI_PC2 + HI_PC3 +
              Water_PC1 + Water_PC2 + Water_PC3 +        
              # Natural_areas + 
              sqrt(Macrophyte_biomass) + 
              Fish_abundance  +
              (1|System_id:Year:Country), # (1|Random_effects), 
            data = Benth,
            family="poisson")

performance::check_convergence(m3)

## assumptions 
plot(m3) 
qqnorm(resid(m3))
qqline(resid(m3))
check_outliers(m3)

### singularity fit
performance::check_singularity(m3)

### multicolinearity
car::vif(m3) 
performance::check_collinearity(m3)
# System_type correlates with HI_PC1      

summary(m3)
# car::Anova(m3)

# check overdispersion 
sum(residuals(m3, type = "pearson")^2) / df.residual(m3)
# or
# StatisticalModels::GLMEROverdispersion(model = m3)

### Pairwise comparisons for System_type----

emmeans(m3, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters)

cld(emmeans(m3, list(pairwise ~ System_type)),  
    #  type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = SR)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Species richness")+
  theme_bw()

# remove System_type

m_SR <- glmer(SR ~   
                # System_type +
                HI_PC1 + HI_PC2 + HI_PC3 +
                Water_PC1 + Water_PC2 + Water_PC3 +        
                # Natural_areas + 
                sqrt(Macrophyte_biomass) + 
                Fish_abundance +
                (1|System_id:Year:Country), #  (1|Random_effects), 
              data = Benth,
              family="poisson")

performance::check_convergence(m_SR)

### singularity fit
performance::check_singularity(m_SR)

### multicolinearity
car::vif(m_SR) 
performance::check_collinearity(m_SR)

# check overdispersion 
sum(residuals(m_SR, type = "pearson")^2) / df.residual(m_SR)

summary(m_SR)
car::Anova(m_SR)

### AICc comparison----

m_SR2 <- glmer(SR ~   
                # System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
                Water_PC1 + Water_PC2 + Water_PC3 +        
                # Natural_areas + 
                sqrt(Macrophyte_biomass) + 
                Fish_abundance +
                (1|System_id:Year:Country), #  (1|Random_effects), 
              data = Benth,
              family="poisson")
library(MuMIn)
arrange(AICc(m_SR, m_SR2), 
        +AICc)

car::Anova(m_SR2)

### Plot predictor effects ----

set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          legend.pos = "top")


plot_model(m_SR, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Species richness")) +
  geom_point(data = Benth, 
             mapping = aes(x = Macrophyte_biomass, y = SR),  
             fill= "#0072B2", size=2, shape=21)




# Functional diversity measures----


# FEve              
# FDiv              
# FDis             

# FRic and RaoQ are not considered, 
# becouse FRic correlates with SR, while RaoQ correlates with FDis



#--------------------------#

## FEve ----

#--------------------------#

m6 <- lmer(FEve ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +        
             # Natural_areas + 
             sqrt(Macrophyte_biomass) + 
             Fish_abundance  +
             (1|System_id:Year:Country), # (1|Random_effects), 
           data = Benth)

performance::check_convergence(m6)

## assumptions 
plot(m6) 
qqnorm(resid(m6))
qqline(resid(m6))
check_outliers(m6)

### singularity fit
performance::check_singularity(m6)

### multicolinearity
car::vif(m6) 
performance::check_collinearity(m6)
# System_type correlates with HI_PC1      

summary(m6)
car::Anova(m6)

### Pairwise comparisons for System_type----

emmeans(m6, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters)

cld(emmeans(m6, list(pairwise ~ System_type)),  
    #  type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = FEve)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="FEve")+
  theme_bw()

# remove System_type

m_FEve <- lmer(FEve ~ # System_type +
                 HI_PC1 + HI_PC2 + HI_PC3 +
                 Water_PC1 + Water_PC2 + Water_PC3 +        
                 # Natural_areas + 
                 sqrt(Macrophyte_biomass) + 
                 # sqrt(Fish_abundance)  +
                 Fish_abundance +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benth)

performance::check_convergence(m_FEve)

### multicolinearity
car::vif(m_FEve) 
performance::check_collinearity(m_FEve)

summary(m_FEve)
car::Anova(m_FEve)


# AICc comparisons -----
# REML = F

m_FEve1 <- lmer(FEve ~ # System_type +
                  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  # Natural_areas + 
                  sqrt(Macrophyte_biomass) + 
                  # sqrt(Fish_abundance)  +
                  Fish_abundance +
                  (1|System_id:Year:Country), # (1|Random_effects), 
                data = Benth, REML = F)

m_FEve2 <- lmer(FEve ~ # System_type +
                #  HI_PC1 + HI_PC2 + HI_PC3 +
                Water_PC1 + Water_PC2 + Water_PC3 +        
                 # Natural_areas + 
                 sqrt(Macrophyte_biomass) + 
                  # sqrt(Fish_abundance)  +
                    Fish_abundance +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benth, REML = F)







library(MuMIn)
arrange(AICc(m_FEve1, m_FEve2), 
        +AICc)



car::Anova(m_FEve2)


plot_model(m_FEve2, type = "pred", terms = c( "Water_PC1"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC1", "FEve")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC1, y = FEve),  
             fill= "#0072B2", size=2, shape=21)

# marginally significant

plot_model(m_FEve2, type = "pred", terms = c( "Water_PC3"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC3", "FEve")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC3, y = FEve),  
             fill= "#0072B2", size=2, shape=21)

plot_model(m_FEve2, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "FEve")) +
  geom_point(data = Benth, 
             mapping = aes(x = Macrophyte_biomass, y = FEve),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m_FEve2, type = "pred", terms = c( "Fish_abundance"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Fish abundance", "FEve")) +
  geom_point(data = Benth, 
             mapping = aes(x = Fish_abundance, y = FEve),  
             fill= "#0072B2", size=2, shape=21)

 #--------------------------#

## FDiv ----

#--------------------------#

m7 <- lmer(FDiv ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +        
             # Natural_areas + 
             sqrt(Macrophyte_biomass) + 
             Fish_abundance  +
             (1|System_id:Year:Country), # (1|Random_effects), 
           data = Benth)

performance::check_convergence(m7)

## assumptions 
plot(m7) 
qqnorm(resid(m7))
qqline(resid(m7))
check_outliers(m7)

### singularity fit
performance::check_singularity(m7)

### multicolinearity
car::vif(m7) 
performance::check_collinearity(m7)
# System_type correlates with HI_PC1      

summary(m7)
# car::Anova(m7)

### Pairwise comparisons for System_type----

library(emmeans)
library(multcomp)

emmeans(m7, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters)

cld(emmeans(m7, list(pairwise ~ System_type)),  
    #  type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = FDiv)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="FDiv")+
  theme_bw()

# remove System_type

m_FDiv <- lmer(FDiv ~  # System_type +
                 HI_PC1 + 
                 HI_PC2 + 
                 HI_PC3 +
                 Water_PC1 +  Water_PC2 + 
                 Water_PC3 +        
                 # Natural_areas + 
                 sqrt(Macrophyte_biomass) + 
                 Fish_abundance  +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benth)

performance::check_convergence(m_FDiv)

### singularity fit
performance::check_singularity(m_FDiv)

### multicolinearity
car::vif(m_FDiv) 
performance::check_collinearity(m_FDiv)

summary(m_FDiv)
car::Anova(m_FDiv)


# AICc comparisons -----
# REML = F

m_FDiv1 <- lmer(FDiv ~ # System_type +
                  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  # Natural_areas + 
                  sqrt(Macrophyte_biomass) + 
                  # sqrt(Fish_abundance)  +
                  Fish_abundance +
                  (1|System_id:Year:Country), # (1|Random_effects), 
                data = Benth, REML = F)

m_FDiv2 <- lmer(FDiv ~ # System_type +
                  #  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  # Natural_areas + 
                  sqrt(Macrophyte_biomass) + 
                  # sqrt(Fish_abundance)  +
                  Fish_abundance +
                  (1|System_id:Year:Country), # (1|Random_effects), 
                data = Benth, REML = F)







library(MuMIn)
arrange(AICc(m_FDiv1, m_FDiv2), 
        +AICc)


car::Anova(m_FDiv2)


plot_model(m_FDiv2, type = "pred", terms = c( "Water_PC1"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC1", "FDiv")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC1, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)

plot_model(m_FDiv2, type = "pred", terms = c( "Water_PC2"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC2", "FDiv")) +
  xlim(-2,3)+
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC2, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m_FDiv2, type = "pred", terms = c( "Fish_abundance"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Fish abundance", "FDiv")) +
  geom_point(data = Benth, 
             mapping = aes(x = Fish_abundance, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)

# marginally significant 

plot_model(m_FDiv, type = "pred", terms = c( "Water_PC2"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC2", "FDiv")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC2, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)



plot_model(m_FDiv2, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "FDiv")) +
  geom_point(data = Benth, 
             mapping = aes(x = Macrophyte_biomass, y = FDiv),  
             fill= "#0072B2", size=2, shape=21)




#--------------------------#

## FDis ----

#--------------------------#

m8 <- lmer(FDis ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +        
             # Natural_areas + 
             sqrt(Macrophyte_biomass) + 
             Fish_abundance  +
             (1|System_id:Year:Country), # (1|Random_effects), 
           data = Benth)

performance::check_convergence(m8)

## assumptions 
plot(m8) 
qqnorm(resid(m8))
qqline(resid(m8))
check_outliers(m8)


### ->  transform the response
m8 <- lmer(sqrt(FDis) ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +        
             # Natural_areas + 
             sqrt(Macrophyte_biomass) + 
             Fish_abundance  +
             (1|System_id:Year:Country), # (1|Random_effects), 
           data = Benth)

## assumptions 
plot(m8) 
qqnorm(resid(m8))
qqline(resid(m8))
check_outliers(m8)

### singularity fit
performance::check_singularity(m8)

### multicolinearity
car::vif(m8) 
performance::check_collinearity(m8)
# System_type correlates with HI_PC1      

summary(m8)
# car::Anova(m8)

### Pairwise comparisons for System_type----

library(emmeans)
library(multcomp)

emmeans(m8, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters)

cld(emmeans(m8, list(pairwise ~ System_type)),  
    #  type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = FDis)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="FDis")+
  theme_bw()

# remove System_type

m_FDis <- lmer(sqrt(FDis) ~  # System_type +
                 HI_PC1 +  
                 HI_PC2 + HI_PC3 +
                 Water_PC1 + Water_PC2 + Water_PC3 +        
                 # Natural_areas + 
                 sqrt(Macrophyte_biomass) + 
                 Fish_abundance  +
                 (1|System_id:Year:Country), # (1|Random_effects), 
               data = Benth, REML = F)

performance::check_convergence(m_FDis)

### singularity fit
performance::check_singularity(m_FDis)

### multicolinearity
car::vif(m_FDis) 
performance::check_collinearity(m_FDis)

summary(m_FDis)
car::Anova(m_FDis)


# AICc comparisons -----
# REML = F

m_FDis1 <- lmer(sqrt(FDis) ~ # System_type +
                  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  # Natural_areas + 
                  sqrt(Macrophyte_biomass) + 
                  # sqrt(Fish_abundance)  +
                  Fish_abundance +
                  (1|System_id:Year:Country), # (1|Random_effects), 
                data = Benth, REML = F)

m_FDis2 <- lmer(sqrt(FDis) ~ # System_type +
                  #  HI_PC1 + HI_PC2 + HI_PC3 +
                  Water_PC1 + Water_PC2 + Water_PC3 +        
                  # Natural_areas + 
                  sqrt(Macrophyte_biomass) + 
                  # sqrt(Fish_abundance)  +
                  Fish_abundance +
                  (1|System_id:Year:Country), # (1|Random_effects), 
                data = Benth, REML = F)







library(MuMIn)
arrange(AICc(m_FDis1, m_FDis2), 
        +AICc)


car::Anova(m_FDis2)

plot_model(m_FDis2, type = "pred", terms = c( "Water_PC3"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC3", "FDis")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC3, y = FDis),  
             fill= "#0072B2", size=2, shape=21)



plot_model(m_FDis2, type = "pred", terms = c( "Water_PC2"),  show.data=F, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC2", "FDis")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC2, y = FDis),  
             fill= "#0072B2", size=2, shape=21)
