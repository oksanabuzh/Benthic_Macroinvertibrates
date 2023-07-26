# Mixed effect models (G)LMM
# testing the effects of human impact and system properties on predictor variables

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

Benthic <- Benth_dat %>%
  full_join(FuncDiv, by="Sampling_site") %>% 
  mutate(System_id = as_factor(System_id), Year = as_factor(Year), 
         Country=as_factor(Country),
       #  Waste.ef=as_factor(Waste.effluent),
         Waste.ef=recode_factor(Waste.effluent,
                  "0"="absent", "1"="present"),
         Gravel_ef=as_factor(Gravel_eff),
         Agric_m_log=log(Agriculture_m),
         Highway_m_log=log(Highway_m)) 

str(Benthic)

# Human impact dependency on system type:

#Depth----
m1 <- lmer(Depth_cm ~ System_type +
  (1|System_id:Year:Country), data = Benthic)

## assumptions 
plot(m1) # heteroscedasticity  ->  transform response
qqnorm(resid(m1))
qqline(resid(m1))

m1.1 <- lmer(log(Depth_cm) ~ System_type +
             (1|System_id:Year:Country), data = Benthic)

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

ggplot(Benthic, aes(x = System_type, y = (Depth_cm))) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="System size (depth, cm)")+
  theme_bw()


#HI----

# Waste.ef
# Gravel_ef
# Agric_m_log
# Highway_m_log

## Waste.ef ----

# Sample size per system type
table(Benthic$Waste.ef)
table(Benthic$Waste.ef, Benthic$System_type)


s <- Benthic %>% 
    count(Waste.ef, System_type)  
s 

ggplot(s, aes(fill=Waste.ef, y=n, x=System_type)) + 
  geom_bar(stat="identity" ) + #, position=position_dodge()
  labs(y = "Observations, n", x = "System type", 
       fill= "Waste effluent") +
  theme_bw() +
  theme(axis.text.y=element_text(colour = "black", size=10),
        axis.text.x=element_text(colour = "black", size=9),
        axis.title=element_text(size=10),
        legend.text=element_text(colour = "black", size=7),
        legend.title=element_text(size=8),
        legend.position="right")

## Gravel_eff ----

g <- Benthic %>% 
  count(Gravel_ef, System_type)  
g 

ggplot(g, aes(fill=Gravel_ef, y=n, x=System_type)) + 
  geom_bar(stat="identity" ) + #, position=position_dodge()
  scale_fill_brewer(palette = "BuPu") +
  #  scale_color_manual(col_g)+
labs(y = "Observations, n", x = "System type", 
       fill= "Gravel exploitation") +
  theme_bw() +
  theme(axis.text.y=element_text(colour = "black", size=10),
        axis.text.x=element_text(colour = "black", size=9),
        axis.title=element_text(size=10),
        legend.text=element_text(colour = "black", size=7),
        legend.title=element_text(size=8),
        legend.position="right") 


ggplot(Benthic, aes(x = System_type, y = Gravel_eff)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Gravel Exploitation")+
  theme_bw() +
  theme(axis.text.y=element_text(colour = "black", size=10),
        axis.text.x=element_text(colour = "black", size=9),
        axis.title=element_text(size=10)) 

m2 <- lmer(Gravel_eff ~ System_type +(1|System_id:Year), data = Benthic)

car::Anova(m2)
summary(m2)


  

##Agric_m----
m3 <- lmer(Agric_m_log ~ System_type +
             (1|System_id:Year:Country), data = Benthic)

## assumptions 
plot(m3) # heteroscedasticity ?
check_heteroscedasticity(m3)
qqnorm(resid(m3))
qqline(resid(m3))


summary(m3)
car::Anova(m3)

### Pairwise comparisons
emmeans(m3, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m3, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benthic, aes(x = System_type, y = Agriculture_m)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Agriculture distance, m")+
  theme_bw()




##Highway_m_log----
m4 <- lmer(Highway_m_log ~ System_type +
             (1|System_id:Year:Country), data = Benthic)

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

ggplot(Benthic, aes(x = System_type, y = Highway_m)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Highway distance, m")+
  theme_bw()

# Water properties----

## Water_PC1----


m5 <- lmer(Water_PC1 ~ 
             System_type +
             Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), data = Benthic)
plot(m5)
check_outliers(m5)

# transform response

min(Benthic$Water_PC1)

Benthic$Water_PC1_log <- log(Benthic$Water_PC1+2.6)

m5 <- lmer(Water_PC1_log ~ 
            System_type +
             Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), data = Benthic)

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

# Plot effects:

# Highway distance

plot_model(m5, type = "pred", terms = c( "Highway_m_log"),  
           show.data=F, jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway distance", "Water PC1")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Highway_m_log, y = Water_PC1_log),  
             fill= "#0072B2", size=2, shape=21)



### Pairwise comparisons

# System_type

emmeans(m5, list(pairwise ~ System_type), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m5, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benthic, aes(x = System_type, y = Water_PC1)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Water PC1")+
  theme_bw()


# Waste.ef

emmeans(m5, list(pairwise ~ Waste.ef), adjust = "tukey", 
        Letters = letters, type="response")

cld(emmeans(m5, list(pairwise ~ Waste.ef)),  
    type="response",
    Letters = letters, adjust = "tukey")

ggplot(Benthic, aes(x = Waste.ef, y = Water_PC1)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Waste effluent", y="Water PC1")+
  theme_bw()




## Water_PC2----

m6 <- lmer(Water_PC2 ~ 
             System_type +
             Waste.ef + Gravel_ef + 
             Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), 
           data = Benthic)

## assumptions 
plot(m6) # heteroscedasticity ?
check_heteroscedasticity(m6)
check_outliers(m6)
Benthic[16, 1]


m6 <- lmer(Water_PC2 ~ 
              System_type +
             Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), 
           data = Benthic %>%  
             filter(!Sampling_site==("VI1_2016")))


plot(m6) 
check_heteroscedasticity(m6)
check_outliers(m6)

qqnorm(resid(m6))
qqline(resid(m6))


### singularity fit
performance::check_singularity(m6)
performance::check_convergence(m6)

### multicolinearity
car::vif(m6) 
performance::check_collinearity(m6)


summary(m6)
car::Anova(m6)

# remove system type

m6 <- lmer(Water_PC2 ~ 
            # System_type +
             Waste.ef + Gravel_ef + 
             Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), 
           data = Benthic %>% 
             filter(!Sampling_site==("VI1_2016")))

plot(m6) 
check_heteroscedasticity(m6)


car::Anova(m6)


# Gravel_ef
emmeans(m6, list(pairwise ~ Gravel_ef), adjust = "none", 
        Letters = letters, type="response")

cld(emmeans(m6, list(pairwise ~ Gravel_ef)),  
    type="response", Letters = letters, adjust = "none")

ggplot(Benthic %>% 
         filter(!Sampling_site==("VI1_2016")),
       aes(x = Gravel_ef, y = Water_PC2)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Gravel exploitation", y="Water PC2")+
  theme_bw()


# Waste.ef
emmeans(m6, list(pairwise ~ Waste.ef), adjust = "none", 
        Letters = letters )
cld(emmeans(m6, list(pairwise ~ Waste.ef)),  
    type="response", Letters = letters, adjust = "none")

ggplot(Benthic%>%
         filter(!Sampling_site==("VI1_2016")), 
       aes(x = Waste.ef, y = Water_PC2)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Waste effluent", y="Water PC2")+
  theme_bw()


# System_type
emmeans(m6, list(pairwise ~ System_type ), adjust = "none", 
        Letters = letters )
cld(emmeans(m6, list(pairwise ~ System_type )),  
    type="response", Letters = letters, adjust = "none")

ggplot(Benthic%>%
         filter(!Sampling_site==("VI1_2016")), 
       aes(x = System_type , y = Water_PC2)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Water PC2")+
  theme_bw()

# plot Highway

plot_model(m6, type = "pred", terms = c( "Highway_m_log"),  
           show.data=T, 
           jitter = 0, 
           title = "", dot.alpha=0.8, line.size=0.1, 
          # axis.title = c("Highway_m_log", "W PC1")
          ) +
  geom_point(data = Benthic %>% 
               filter(!Sampling_site=="III1_2016"), 
             mapping = aes(x = SR, y = FRic_log),  
             fill= "#0072B2", size=2, shape=21)




## Water_PC3----

m7 <- lmer(Water_PC3 ~ 
             System_type +
             Waste.ef + Gravel_ef + 
             Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), 
           data = Benthic)

## assumptions 
plot(m7) # heteroscedasticity ?
check_heteroscedasticity(m7)
qqnorm(resid(m7))
qqline(resid(m7))

check_outliers(m7)

Benthic[16, 1]

m7 <- lmer(Water_PC3 ~ 
             System_type +
             Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), 
             data = Benthic %>%  
               filter(!Sampling_site==("VI1_2016")))

plot(m7) # heteroscedasticity ?
check_heteroscedasticity(m7)
qqnorm(resid(m7))
qqline(resid(m7))

check_outliers(m7)

car::Anova(m7)


### Pairwise comparisons
emmeans(m7, list(pairwise ~ Waste.ef), 
        adjust = "tukey", #adjust = "none", 
        Letters = letters, type="response")

cld(emmeans(m7, list(pairwise ~ Waste.ef)),  
    type="response",
    Letters = letters, adjust = "none")


ggplot(Benthic %>%  
         filter(!Sampling_site==("VI1_2016")),
       aes(x = Waste.ef, y = Water_PC3)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="Waste effluent", y="Water PC3")+
  theme_bw()


# System_type
emmeans(m7, list(pairwise ~ System_type), 
        adjust = "tukey", #adjust = "none", 
        Letters = letters, type="response")

cld(emmeans(m7, list(pairwise ~ System_type)),  
    type="response",
    Letters = letters, adjust = "none")

ggplot(Benthic %>%  
         filter(!Sampling_site==("VI1_2016")),
       aes(x = System_type, y = Water_PC3)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System ype", y="Water PC3")+
  theme_bw()


####  remove System_type----

m7 <- lmer(Water_PC3 ~ 
            # System_type +
             Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
             (1|System_id:Year:Country), 
           data = Benthic  %>% 
             filter(!Sampling_site==("VI1_2016")))

car::Anova(m7)


## Macrophytes----

m8 <- lmer(Macrophyte_biomass ~ 
             System_type +
             Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
             Water_PC1 + Water_PC2 + Water_PC3 +
               (1|System_id:Year:Country), data = Benthic)

## assumptions 
plot(m8) # heteroscedasticity ?
check_heteroscedasticity(m8)
qqnorm(resid(m8))
qqline(resid(m8))

Benthic$Macroph_log <- log(Benthic$Macrophyte_biomass+1)

m8.1 <- lmer(Macroph_log ~ 
              # System_type +
              #  Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
               Water_PC1 + Water_PC2 + Water_PC3 +
               (1|System_id:Year:Country), data = Benthic%>% 
               filter(!Sampling_site==("VI1_2016")))

## assumptions 
plot(m8.1) # heteroscedasticity ?
check_heteroscedasticity(m8.1)
qqnorm(resid(m8.1))
qqline(resid(m8.1))
check_outliers(m8.1)


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

ggplot(Benthic%>% 
         filter(!Sampling_site==("VI1_2016")), 
       aes(x = System_type, y = Macrophyte_biomass)) +
  geom_jitter(width = 0.15, fill= "#0072B2", size=2, shape=21) +
  geom_boxplot(notch=F, alpha = 0, color = "black") +
  labs(x="System type", y="Macrophyte biomass")+
  theme_bw()


plot_model(m8.1, type = "pred", terms = c( "Water_PC1"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Water_PC1", "Macrophyte biomass")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Water_PC1, y = Macroph_log),  
             fill= "#0072B2", size=2, shape=21)


#### use raw HI data----

ggplot(Benthic, aes(log(Highway_m), Macrophyte_biomass)) + geom_point()


m8.3 <- lmer(log(Macrophyte_biomass+1) ~ 
              # System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
               Waste.ef + 
               Gravel_ef + 
               Agric_m_log + Highway_m_log +
               (1|System_id:Year:Country), data = Benthic)


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
             (1|System_id:Year:Country), data = Benthic)

## assumptions 
plot(m9) # heteroscedasticity ?
check_heteroscedasticity(m9)
qqnorm(resid(m9))
qqline(resid(m9))

min(Benthic$Fish_abundance)

m9.1 <- lmer(log(Fish_abundance+9) ~ System_type +
               HI_PC1 + HI_PC2 + HI_PC3 + 
               Water_PC1 + Water_PC2 + Water_PC3 +
               Macrophyte_biomass +
               # Waste.ef + Gravel_ef + Agriculture_m + log(Highway_m) +
               (1|System_id:Year:Country), data = Benthic)

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

ggplot(Benthic, aes(x = System_type, y = Fish_abundance)) +
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
               (1|System_id:Year:Country), data = Benthic)


### multicolinearity
car::vif(m9.2) 
performance::check_collinearity(m9.2)


summary(m9.2)
car::Anova(m9.2)

plot_model(m9.2, type = "pred", terms = c( "Macrophyte_biomass"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Macrophyte biomass", "Fish abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Macrophyte_biomass, y = Fish_abundance),  
             fill= "#0072B2", size=2, shape=21)


#### use raw HI data----

ggplot(Benthic, aes(log(Highway_m), Fish_abundance)) + geom_point()


Benthic$Fish_log <-  log(Benthic$Fish_abundance+9)


m9.3 <- lmer(Fish_log ~ 
              # System_type +
               # HI_PC1 + HI_PC2 + HI_PC3 +
               Waste.ef + 
               Gravel_ef + 
               Agric_m_log + 
                 Highway_m_log +
               (1|System_id:Year:Country), data = Benthic)


### multicolinearity
car::vif(m9.3) 
performance::check_collinearity(m9.3)


summary(m9.3)
car::Anova(m9.3)

plot_model(m9.3, type = "pred", terms = c( "Agric_m_log"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Agrofields distance", "Fish abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Agric_m_log, y = Fish_log),  
             fill= "#0072B2", size=2, shape=21)



plot_model(m9.3, type = "pred", terms = c( "Highway_m_log"),  show.data=F, 
           jitter = 0.05, 
           title = "", dot.alpha=0.8, line.size=0.1, 
           axis.title = c("Highway distance", "Fish abundance")) +
  geom_point(data = Benthic, 
             mapping = aes(x = Highway_m_log, y = Fish_log),  
             fill= "#0072B2", size=2, shape=21)
#####