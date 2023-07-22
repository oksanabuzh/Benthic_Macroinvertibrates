# Mixed effect models (G)LMM
# testing the effcets of human impact and system properties on the predictor variables

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


Benth <- Benth_dat %>%
  full_join(FuncDiv, by="Sampling_site") %>% 
  mutate(System_id = as_factor(System_id), Year = as_factor(Year), 
         Country=as_factor(Country),
         Waste.ef=as_factor(Waste.effluent),
         Gravel_ef=as_factor(Gravel_eff),
         Agric_m_log=log(Agriculture_m),
         Highway_m_log=log(Highway_m)) 

str(Benth)

# Human impact relations:

#Depth----
m1 <- lmer(Depth_cm ~ System_type +
  (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m1) # heteroscedasticity  ->  transform response
qqnorm(resid(m1))
qqline(resid(m1))

m1.1 <- lmer(log(Depth_cm) ~ System_type +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m1.1) # heteroscedasticity  ->  transform response
qqnorm(resid(m1.1))
qqline(resid(m1.1))

summary(m1.1)
car::Anova(m1.1)

### Pairwise comparisons
emmeans(m1.1, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters)

cld(emmeans(m1.1, list(pairwise ~ System_type)),  
    #  type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = (Depth_cm))) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="System size (depth, cm)")+
  theme_bw()


#HI----

##HI_PC1----
m2 <- lmer(HI_PC1 ~ System_type +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m2) # heteroscedasticity  ->  transform response
qqnorm(resid(m2))
qqline(resid(m2))

m2.1 <- lmer(log(HI_PC1+2) ~ System_type +
               (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m2.1) # heteroscedasticity  ->  transform response
qqnorm(resid(m2.1))
qqline(resid(m2.1))

summary(m2.1)
car::Anova(m2.1)

### Pairwise comparisons
emmeans(m2.1, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m2.1, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = log(HI_PC1+2))) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="HI_PC1")+
  theme_bw()


##HI_PC2----
m3 <- lmer(HI_PC2 ~ System_type +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m3) # heteroscedasticity ?
check_heteroscedasticity(m3) # no heteroscedasticity
qqnorm(resid(m3))
qqline(resid(m3))

check_outliers(m3)


summary(m3)
car::Anova(m3)

### Pairwise comparisons
emmeans(m3, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m3.1, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = HI_PC2)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="HI_PC2")+
  theme_bw()

##HI_PC3----
m4 <- lmer(HI_PC3 ~ System_type +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m4) # heteroscedasticity ?
check_heteroscedasticity(m4)
qqnorm(resid(m4))
qqline(resid(m4))


summary(m4)
car::Anova(m4)

### Pairwise comparisons
emmeans(m4, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m4, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = HI_PC3)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="HI_PC3")+
  theme_bw()



# Water properties----

## Water_PC1----

m5 <- lmer(Water_PC1 ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m5) # heteroscedasticity ?
check_heteroscedasticity(m5)
qqnorm(resid(m5))
qqline(resid(m5))

### singularity fit
performance::check_singularity(m5)
performance::check_convergence(m5)

### multicolinearity
car::vif(m5) 
performance::check_collinearity(m5)


summary(m5)
car::Anova(m5)

### Pairwise comparisons
emmeans(m5, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m5, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = Water_PC1)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Water_PC1")+
  theme_bw()



####  remove System_type----

m5.1 <- lmer(Water_PC1 ~ #System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
             (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m5.1) 
performance::check_collinearity(m5.1)


summary(m5.1)
car::Anova(m5.1)

set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          legend.pos = "top")

plot_model(m5.1, type = "pred", terms = c( "HI_PC1"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("HI_PC1", "Water_PC1")) +
  geom_point(data = Benth, 
             mapping = aes(x = HI_PC1, y = Water_PC1),  
             fill= "#0072B2", size=2, shape=21)

#### use raw HI data----

ggplot(Benth, aes(log(Highway_m), Water_PC1)) + geom_point()
ggplot(Benth, aes(log(Agriculture_m), Water_PC1)) + geom_point()


m5.2 <- lmer(Water_PC1 ~ #System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
                Waste.ef + 
               Gravel_ef + 
               Agric_m_log + Highway_m_log +
               (1|System_id:Year:Country), data = Benth)




### multicolinearity
car::vif(m5.2) 
performance::check_collinearity(m5.2)


summary(m5.2)
car::Anova(m5.2)

set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          legend.pos = "top")

plot_model(m5.2, type = "pred", terms = c( "Highway_m_log"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway distance", "Water_PC1")) +
  geom_point(data = Benth, 
             mapping = aes(x = Highway_m_log, y = Water_PC1),  
             fill= "#0072B2", size=2, shape=21)


plot_model(m5.2, type = "pred", terms = c( "Agric_m_log"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Agrofields distance", "Water_PC1")) +
  geom_point(data = Benth, 
             mapping = aes(x = Agric_m_log, y = Water_PC1),  
             fill= "#0072B2", size=2, shape=21)


### Pairwise comparisons 
# for Gravel_ef

emmeans(m5.2, list(pairwise ~ Gravel_ef), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m5.2, list(pairwise ~ Gravel_ef)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = Gravel_ef, y = Water_PC1)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Gravel exploitation", y="Water_PC1")+
  theme_bw()

# Waste.ef
emmeans(m5.2, list(pairwise ~ Waste.ef), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m5.2, list(pairwise ~ Waste.ef)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = Waste.ef, y = Water_PC1)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Waste effluent presence", y="Water_PC1")+
  theme_bw()




## Water_PC2----

m6 <- lmer(Water_PC2 ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m6) # heteroscedasticity ?
check_heteroscedasticity(m6)
qqnorm(resid(m6))
qqline(resid(m6))


min(Benth$Water_PC2)
m6.1 <- lmer(sqrt(Water_PC2+2.6) ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m6.1) # heteroscedasticity ?
check_heteroscedasticity(m6.1)
qqnorm(resid(m6.1))
qqline(resid(m6.1))





### singularity fit
performance::check_singularity(m6.1)
performance::check_convergence(m6.1)

### multicolinearity
car::vif(m6.1) 
performance::check_collinearity(m6.1)


summary(m6.1)
car::Anova(m6.1)

### Pairwise comparisons
emmeans(m6.1, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m6.1, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = Water_PC2)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Water_PC2")+
  theme_bw()



####  remove System_type----

m6.2 <- lmer(sqrt(Water_PC2+2.6) ~ #System_type +
               HI_PC1 + HI_PC2 + HI_PC3 +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m6.2) 
performance::check_collinearity(m6.2)


summary(m6.2)
car::Anova(m6.2)

set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          legend.pos = "top")

plot_model(m6.2, type = "pred", terms = c( "HI_PC2"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("HI_PC2", "Water_PC2")) +
  geom_point(data = Benth, 
             mapping = aes(x = HI_PC2, y = Water_PC2),  
             fill= "#0072B2", size=2, shape=21)

#### use raw HI data----

ggplot(Benth, aes(log(Highway_m), Water_PC2)) + geom_point()


m6.3 <- lmer(Water_PC2 ~ #System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
               Waste.ef + 
               Gravel_ef + 
               Agric_m_log + Highway_m_log +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m6.3) 
performance::check_collinearity(m6.3)


summary(m6.3)
car::Anova(m6.3)


### Pairwise comparisons 

# Waste.ef
emmeans(m6.3, list(pairwise ~ Waste.ef), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m6.3, list(pairwise ~ Waste.ef)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = Waste.ef, y = Water_PC2)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Waste effluent presence", y="Water_PC2")+
  theme_bw()



## Water_PC3----

m7.1 <- lmer(Water_PC3 ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m7.1) # heteroscedasticity ?
check_heteroscedasticity(m7.1)
qqnorm(resid(m7.1))
qqline(resid(m7.1))


### singularity fit
performance::check_singularity(m7.1)
performance::check_convergence(m7.1)

### multicolinearity
car::vif(m7.1) 
performance::check_collinearity(m7.1)


summary(m7.1)
car::Anova(m7.1)

### Pairwise comparisons
emmeans(m7.1, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m7.1, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = Water_PC3)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Water_PC3")+
  theme_bw()



####  remove System_type----

m7.2 <- lmer(Water_PC3 ~ #System_type +
               HI_PC1 + HI_PC2 + HI_PC3 +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m7.2) 
performance::check_collinearity(m7.2)


summary(m7.2)
car::Anova(m7.2)

set_theme(base = theme_bw(),
          axis.textsize.x = 1, axis.textsize.y = 1, axis.textcolor = "black",
          axis.title.color = "black", axis.title.size = 1.4,
          legend.pos = "top")

plot_model(m7.2, type = "pred", terms = c( "HI_PC2"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("HI_PC2", "Water_PC3")) +
  geom_point(data = Benth, 
             mapping = aes(x = HI_PC2, y = Water_PC3),  
             fill= "#0072B2", size=2, shape=21)

#### use raw HI data----

ggplot(Benth, aes(log(Highway_m), Water_PC3)) + geom_point()


m7.3 <- lmer(Water_PC3 ~ #System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
               Waste.ef + 
               Gravel_ef + 
               Agric_m_log + Highway_m_log +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m7.3) 
performance::check_collinearity(m7.3)


summary(m7.3)
car::Anova(m7.3)


### Pairwise comparisons 

# Waste.ef
emmeans(m7.3, list(pairwise ~ Waste.ef), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m7.3, list(pairwise ~ Waste.ef)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = Waste.ef, y = Water_PC3)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Waste effluent presence", y="Water_PC3")+
  theme_bw()



# Gravel_ef
emmeans(m7.3, list(pairwise ~ Gravel_ef), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m7.3, list(pairwise ~ Gravel_ef)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = Gravel_ef, y = Water_PC3)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Gravel exploitation", y="Water_PC3")+
  theme_bw()



## Macrophytes----

m8 <- lmer(Macrophyte_biomass ~ System_type +
               Water_PC1 + Water_PC2 + Water_PC3 +
               HI_PC1 + HI_PC2 + HI_PC3 +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m8) # heteroscedasticity ?
check_heteroscedasticity(m8)
qqnorm(resid(m8))
qqline(resid(m8))


m8.1 <- lmer(log(Macrophyte_biomass+1) ~ System_type +
               HI_PC1 + HI_PC2 + HI_PC3 + 
               Water_PC1 + Water_PC2 + Water_PC3 +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m8.1) # heteroscedasticity ?
check_heteroscedasticity(m8.1)
qqnorm(resid(m8.1))
qqline(resid(m8.1))


### singularity fit
performance::check_singularity(m8.1)
performance::check_convergence(m8.1)

### multicolinearity
car::vif(m8.1) 
performance::check_collinearity(m8.1)


summary(m8.1)
car::Anova(m8.1)

### Pairwise comparisons
emmeans(m8.1, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m8.1, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = Macrophyte_biomass)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Macrophyte biomass")+
  theme_bw()



####  remove System_type----

m8.2 <- lmer(log(Macrophyte_biomass+1) ~ #System_type +
               HI_PC1 + HI_PC2 + HI_PC3 +
               Water_PC1 + Water_PC2 + Water_PC3 +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m8.2) 
performance::check_collinearity(m8.2)


summary(m8.2)
car::Anova(m8.2)


Benth$Macr_log <- log(Benth$Macrophyte_biomass+1)
Benth$Water_PC2_s <- sqrt(Benth$Water_PC2+2.6)

m8.3 <- lmer(Macr_log ~ #System_type +
               HI_PC1 + HI_PC2 + HI_PC3 +
               Water_PC1 + Water_PC2 + Water_PC3 +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)

car::Anova(m8.3)

plot_model(m8.3, type = "pred", terms = c( "Water_PC2"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC2", "Macrophyte biomass")) +
  geom_point(data = Benth, 
             mapping = aes(x = Water_PC2, y = Macr_log),  
             fill= "#0072B2", size=2, shape=21)


#### use raw HI data----

ggplot(Benth, aes(log(Highway_m), Macrophyte_biomass)) + geom_point()


m8.3 <- lmer(log(Macrophyte_biomass+1) ~ 
              # System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
               Waste.ef + 
               Gravel_ef + 
               Agric_m_log + Highway_m_log +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m8.3) 
performance::check_collinearity(m8.3)


summary(m8.3)
car::Anova(m8.3)





## Fish ----

m9 <- lmer(Fish_abundance ~ System_type +
             HI_PC1 + HI_PC2 + HI_PC3 +
             Water_PC1 + Water_PC2 + Water_PC3 +
             Macrophyte_biomass +
             # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
             (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m9) # heteroscedasticity ?
check_heteroscedasticity(m9)
qqnorm(resid(m9))
qqline(resid(m9))

min(Benth$Fish_abundance)

m9.1 <- lmer(log(Fish_abundance+9) ~ System_type +
               HI_PC1 + HI_PC2 + HI_PC3 + 
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macrophyte_biomass +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)

## assumptions 
plot(m9.1) 
check_heteroscedasticity(m9.1)
qqnorm(resid(m9.1))
qqline(resid(m9.1))


### singularity fit
performance::check_singularity(m9.1)
performance::check_convergence(m9.1)

### multicolinearity
car::vif(m9.1) 
performance::check_collinearity(m9.1)


summary(m9.1)
car::Anova(m9.1)

### Pairwise comparisons
emmeans(m9.1, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m9.1, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benth, aes(x = System_type, y = Fish_abundance)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Fish abundance")+
  theme_bw()



####  remove System_type----

m9.2 <- lmer(log(Fish_abundance+9) ~ # System_type +
               HI_PC1 + HI_PC2 + HI_PC3 +
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macrophyte_biomass +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m9.2) 
performance::check_collinearity(m9.2)


summary(m9.2)
car::Anova(m9.2)

plot_model(m9.2, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Fish abundance")) +
  geom_point(data = Benth, 
             mapping = aes(x = Macrophyte_biomass, y = Fish_abundance),  
             fill= "#0072B2", size=2, shape=21)


#### use raw HI data----

ggplot(Benth, aes(log(Highway_m), Fish_abundance)) + geom_point()


Benth$Fish_log <-  log(Benth$Fish_abundance+9)


m9.3 <- lmer(Fish_log ~ 
              # System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
               Waste.ef + 
               Gravel_ef + 
               Agric_m_log + Highway_m_log +
               (1|System_id:Year:Country), data = Benth)


### multicolinearity
car::vif(m9.3) 
performance::check_collinearity(m9.3)


summary(m9.3)
car::Anova(m9.3)

plot_model(m9.3, type = "pred", terms = c( "Agric_m_log"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Agrofields distance", "Fish abundance")) +
  geom_point(data = Benth, 
             mapping = aes(x = Agric_m_log, y = Fish_log),  
             fill= "#0072B2", size=2, shape=21)


#####