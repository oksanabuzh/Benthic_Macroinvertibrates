# NMDS, PERMANOVA analyses, and envfit for functional composition

rm(list=ls(all=TRUE))

library(ggplot2)
library(tidyverse)
library(vegan)

# Data

sp_comp <- read_csv ( "Data/Com_composition.csv")
str(sp_comp)
names(sp_comp)

Benth_dat <- read_csv ("Data/Benth_PCA2.csv")

Benthic <- Benth_dat %>%
  mutate(System_id = as_factor(System_id), Year = as_factor(Year), 
         Country=as_factor(Country),
         Waste.ef=as_factor(Waste.effluent),
         Gravel_ef=as_factor(Gravel_eff),
         Agric_m_log=log(Agriculture_m),
         Highway_m_log=log(Highway_m),
         Macroph_sqrt=sqrt(Macrophyte_biomass))

str(sp_comp)
names(sp_comp)
compos <-sp_comp [,2:142] 
str(compos)
names(compos)

# NMDS analysis----

nmds1<-metaMDS(compos, distance="bray",k=2,trymax=100)
nmds1 # the stress value shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
scores(nmds1)


data <- read_csv ("Data/Benth_PCA2.csv")
str(data)


# PERMANOVA----

PERM_mod1 <- adonis2(compos ~ # System_type +
                  Agric_m_log   + Highway_m_log  + Gravel_ef +  Waste.ef,
                     data=Benthic,
                 # strata=Benthic$System_id,
                   permutations = 1000, method = "bray")

PERM_mod1

# save
write.csv(PERM_mod1 , file="Results/SpComp_PERMANOVA_Mod1.csv")

PERM_mod2 <- adonis2(compos ~ # System_type +
                        Water_PC1 + Water_PC2 + Water_PC3 +
                       Macrophyte_biomass + Fish_abundance ,
                       data=Benthic,
                    # strata=Benthic$System_id,
                     permutations = 1000, method = "bray")

PERM_mod2

# save
write.csv(PERM_mod2 , file="Results/SpComp_PERMANOVA_Mod2.csv")

## R2 coefficients----

R2_mod1 <- as.data.frame(PERM_mod1) %>%
  mutate (Variables = rownames(PERM_mod1)) %>% 
  dplyr::select(Variables, R2) %>% 
  filter(!Variables=="Residual",
         !Variables=="Total")

R2_mod2 <- as.data.frame(PERM_mod2) %>%
  mutate (Variables = rownames(PERM_mod2)) %>% 
  dplyr::select(Variables, R2) %>% 
  filter(!Variables==c("Residual"),
         !Variables==c("Total"))

SpComp_R2 <- rbind(R2_mod1, R2_mod2) %>% 
  filter(!Variables==c("System_type")) %>% 
  mutate(Variables=recode_factor(Variables, 
                                 Macrophyte_biomass="Macrophytes",
                                 Waste.ef="Waste",
                                 Gravel_ef="Gravel",
                                 Highway_m_log="Highway_m",
                                 Agric_m_log ="Agric_m",
                                 Fish_abundance ="Fish")) %>% 
  mutate(Variables=fct_relevel(Variables, c("Agric_m", "Highway_m", "Gravel", "Waste", 
                                            "Water_PC1", "Water_PC2", "Water_PC3",
                                            "Macrophytes", "Fish"))) %>% 
  arrange(Variables)

SpComp_R2
# save
write.csv(SpComp_R2 , file="Results/SpComp_R2.csv", row.names=TRUE)

# Plot---- 
### envfit----
# get coordinates
fit <- envfit(nmds1 ~  
                Agric_m_log + Highway_m_log +Waste.effluent + Gravel_eff +
                Water_PC1 + Water_PC2 + Water_PC3 +
                Macrophyte_biomass + Fish_abundance, data=Benthic, perm=1000,
              strata=Benthic$Random_effects) #



fit

var.scrs <- as.data.frame(scores(fit, display = "vectors")) #save variable intrinsic values into data-frame
var.scrs <- cbind(var.scrs, Variables = rownames(var.scrs)) #add variable names to data-frame
var.scrs <- cbind(var.scrs, pval = fit$vectors$pvals) #add p-values to data-frame so we can select variables which are significant
sig.var.scrs <- subset(var.scrs #, pval<=0.07
                       )%>%
  mutate(Variables=recode_factor(Variables, 
                                 Macrophyte_biomass="Macrophytes",
                                 Waste.effluent="Waste",
                                 Gravel_eff="Gravel",
                                 Highway_m_log="Highway_m",
                                 Agric_m_log ="Agric_m",
                                 Water_PC1 ="PC1",
                                 Water_PC2 ="PC2",
                                 Water_PC3 ="PC3",
                                 Fish_abundance ="Fish"))

# write.csv(var.scrs, file="Results/envfit_mod.csv")


# Plot NMDS----

## Plot functional composition----

#### using base plot----

x11(height=9,width=8.5)
par(mai=c(1, 1, 0.2, 0.2), mgp=c(1.8,0.5,0), cex.lab=1, cex.axis=1, pch=21, cex=1.5, lwd=3, tck=0.015)

plot(nmds1, type = "n")
#points(nmds1, display = "sites",pch = 19, cex =1.2, col="gray")
#points(nmds1, display = "species", col=col, pch = 19, cex = 0.8)
orditorp(nmds1,display="species", col="blue", air=0.0001, cex=0.5)
plot(fit, col="black")

#### using ggplot -------- 

library(ggplot2)
library(ggrepel)
library(devtools)

# 


scores(nmds1)
species.scores <- as.data.frame(scores(nmds1, "species"))  #Using the scores function from vegan to extract the variable scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of variables, from the rownames of species.scores
head(species.scores)
scores(nmds1)
sites.scores <-  as.data.frame(scores(nmds1, "sites"))
sites.scores$sites <- rownames(sites.scores)
sites.scores$System.Type <- factor(data$System_type)
sites.scores


species.scores$species




names(comp)

p1  <- ggplot(data = sites.scores, 
                aes(x = NMDS1, y = NMDS2)) +
  geom_vline(xintercept = 0, col="grey", linetype="dashed")+
  geom_hline(yintercept = 0, col="grey", linetype="dashed")+
  scale_shape_manual(values=c( 19, 1))+ 
  geom_point (data = species.scores, size = 2, pch=19, colour= "grey") + #, colour=Col, fill=Col
  geom_text_repel (data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), size = 3.6, colour= "black") + 
  theme(axis.title = element_text(size = 13, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, 
                              colour = "grey30", size=0.5), 
         axis.ticks = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"))

p1

# add elipses showing system type 

p2  <- ggplot(data = sites.scores, 
              aes(x = NMDS1, y = NMDS2)) +
  stat_ellipse(aes(colour=System.Type), alpha=.3, type='t',linewidth =0.1, geom="path")+
  stat_ellipse(aes(fill=System.Type), alpha=.1,type='t',linewidth =1, geom="polygon")+
  geom_vline(xintercept = 0, col="grey", linetype="dashed")+
  geom_hline(yintercept = 0, col="grey", linetype="dashed")+
  scale_shape_manual(values=c( 19, 1))+ 
  geom_point (data = species.scores, size = 2, pch=19, colour= "grey") + 
  geom_text_repel (data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), size = 3.6, colour= "grey") + 
  theme(axis.title = element_text(size = 13, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, 
                                                                        colour = "grey30", size=0.5), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"),
        legend.key = element_blank(), 
        legend.title = element_text(size = 11, face = "bold", colour = "black"), 
        legend.text = element_text(size = 11, colour = "black"),
        legend.position="bottom")+
  labs(col="System Type")+
  guides(fill=guide_legend("System Type"))

p2

## Add significant predictors----

coord_cont <- as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit, 
                                          rescale=TRUE, fill = 1)




# Create scores with the standard length for the posthoc plotting of arrows of the same size
# for this we would use raw vectors from envfit and not scores 
# as arrow lengths (scores) are created as NMDS1 vectors multiplied on sqrt(r2) from fit$vectors 


names(fit$vectors)
fit$vectors$arrows 

# we add standardised scores:
coord_cont <- as.data.frame(scores(fit, "vectors")) %>% 
  cbind( stand= fit$vectors$arrows)
coord_cont
# stand.NMDS1 and stand.NMDS2 are all the same lenght arrows for the posthoc plottig

coord_cont_ <- coord_cont  %>% 
  mutate(Variables = rownames(scores(fit, "vectors"))) %>% 
  mutate(Variables=recode_factor(Variables, 
                                 Macrophyte_biomass="Macrophytes",
                                 Waste.effluent="Waste",
                                 Gravel_eff="Gravel",
                                 Highway_m_log="Highway_m",
                                 Agric_m_log ="Agric_m",
                                 Water_PC1 ="PC1",
                                 Water_PC2 ="PC2",
                                 Water_PC3 ="PC3",
                                 Fish_abundance ="Fish"))  

# rescale  all arrows to fill an ordination plot, where fill =  shows proportion of plot to be filled by the arrows

coord_cont_standrd  <- coord_cont_%>% 
  mutate(stand.NMDS1=stand.NMDS1 * ordiArrowMul(fit, rescale=TRUE, fill = 0.7))%>% 
  mutate(stand.NMDS2=stand.NMDS2 * ordiArrowMul(fit, rescale=TRUE, fill = 0.7))

coord_cont_standrd




p2+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = coord_cont_, linewidth =0.7, alpha = 1, 
               colour = "black", 
               arrow = arrow(length = unit(0.03, "npc")))

p3 <- ggplot(data = sites.scores, 
             aes(x = NMDS1, y = NMDS2)) +
 stat_ellipse(aes(colour=System.Type), alpha=.3, type='t',linewidth =0.1, geom="path")+
 stat_ellipse(aes(fill=System.Type), alpha=.1,type='t',linewidth =1, geom="polygon")+
  geom_vline(xintercept = 0, col="grey", linetype="dashed")+
  geom_hline(yintercept = 0, col="grey", linetype="dashed")+
  #
  geom_segment(aes(x = 0, y = 0, xend = stand.NMDS1, yend = stand.NMDS2), 
               data = coord_cont_standrd, linewidth =0.7, alpha = 1, 
               colour = "black", 
               arrow = arrow(length = unit(0.03, "npc")))+
#
  scale_shape_manual(values=c( 19, 1))+ 
 # geom_point (data = species.scores, size = 2, pch=19, colour= "blue") + 
  geom_text_repel (data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), 
                   size = 3.2, colour= "deepskyblue4") + 
  theme(axis.title = element_text(size = 13, face = "bold", colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, 
                                                                        colour = "grey30", size=0.5), 
        axis.ticks = element_blank(), 
        axis.text = element_text(size = 13, colour = "black"),
        legend.key = element_blank(), 
        legend.title = element_text(size = 11, face = "bold", colour = "black"), 
        legend.text = element_text(size = 11, colour = "black"),
        legend.position="bottom")+
  labs(col="System Type")+
  guides(fill=guide_legend("System Type"))

p3 + 
  xlim(-1.3,1.5)+ylim(-1.1,1)

# add labels
p3 + geom_text(data = coord_cont_standrd, 
            aes(x =stand.NMDS1, y = stand.NMDS2, label = Variables), 
            colour = "black", fontface = "bold", size = 4 ,
            vjust=1, hjust=0.1 # adjust text positions
            )+
  xlim(-1.5,1.5)+ylim(-1.5,1.5)

