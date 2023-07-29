# NMDS, PERMANOVA analyses, and envfit for functional composition

rm(list=ls(all=TRUE))

library(ggplot2)
library(tidyverse)
library(vegan)

# Data

func_comp <- read_csv ( "Data/FuncComp.csv")

Benth_dat <- read_csv ("Data/Benth_PCA2.csv")

Benthic <- Benth_dat %>%
  mutate(System_id = as_factor(System_id), Year = as_factor(Year), 
         Country=as_factor(Country),
         Waste.ef=as_factor(Waste.effluent),
         Gravel_ef=as_factor(Gravel_eff),
         Agric_m_log=log(Agriculture_m),
         Highway_m_log=log(Highway_m),
         Macroph_sqrt=sqrt(Macrophyte_biomass))

str(func_comp)
names(func_comp)
compos <-func_comp [,2:36] 
str(compos)
names(compos)

# rename some traits
comp <- compos %>% 
  rename(RESP_sp=RESP_spi, 
         RLC_mult=RLC_mlt, RLC_flex=RLC_flx,
         "<0.25"=BodSz_0.25, 
         "0.25-0.5"=BodSz_0.25_0.5,
         "0.5-1"= BodSz_0.5_1,
         "1-2"=BodSz_1_2,
         "2-4"=BodSz_2_4,
         "4-8"=BodSz_4_8,
         ">8"=BodSz_8_,
         REPR_ovo=REP_ovo, REPR_fie=REP_fie, REPR_cie=REP_cie, REPR_fic=REP_fic, 
         REPR_frc=REP_frc, REPR_vec=REP_vec, REPR_tec=REP_tec, REPR_asex=REP_ase,
         FG_gra=Feed_Gra, FG_min=Feed_Min, FG_xyl=Feed_Xyl, FG_shr=Feed_Shr, 
         FG_gat=Feed_Gat, FG_aff=Feed_Aff, FG_pff=Feed_Pff, FG_pred=Feed_Pre, FG_par=Feed_Par)


names(comp)

# NMDS analysis----

nmds1<-metaMDS(comp, distance="bray",k=2,trymax=100)
nmds1 # the stress value shows how easy it was to condense multidimensional data into two dimensional space, below 0.2 is generally good
scores(nmds1)


data <- read_csv ("Data/Benth_PCA2.csv")
str(data)

Benth_PC<- data %>%
  mutate(System_id = as_factor(System_id)) %>% 
  mutate(Year = as_factor(Year)) 

names(Benth_PC)


# PERMANOVA----
##Mod1-----
PERM_mod1 <- adonis2(comp ~ 
                      System_type +
                       Waste.effluent + Gravel_eff + Agric_m_log + Highway_m_log ,
                      # strata = Benth_PC$System_id, 
                     # strata = Benth_PC$System_id:Benth_PC$Year,  
                      strata = Benth_PC$Random_effects, 
                                data=Benthic,
                          permutations = 1000, method = "bray")

PERM_mod1

# save

#library(openxlsx)
# write.xlsx(adonis1 , file="adonis1.xlsx")
write.csv(PERM_mod1 , file="Results/PERM_mod1.csv")

##Mod2-----
PERM_mod2 <- adonis2(comp ~ 
               # System_type +
               #  Waste.ef + Gravel_ef + Agric_m_log + Highway_m_log +
                 Agric_m_log + Highway_m_log +Waste.effluent + Gravel_eff +
                  Water_PC1 + Water_PC2 + Water_PC3 +
                  Macrophyte_biomass + Fish_abundance ,
                 # strata = Benth_PC$System_id, 
                     # strata = Benth_PC$System_id:Benth_PC$Year,  
                     # strata = Benth_PC$Random_effects, 
                     data=Benthic,
                     permutations = 1000, method = "bray")

PERM_mod2
# save
write.csv(PERM_mod2 , file="Results/PERM_mod2.csv")


# Plot----
### envfit----
# get coordinates
fit <- envfit(nmds1 ~  
                Agric_m_log + Highway_m_log +Waste.effluent + Gravel_eff +
                Water_PC1 + Water_PC2 + Water_PC3 +
                Macrophyte_biomass + Fish_abundance, data=Benthic, perm=1000) #



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

write.csv(var.scrs, file="Results/envfit_mod.csv")


# Plot NMDS----

## Plot functional composition----

#### using base plot----

x11(height=9,width=8.5)
par(mai=c(1, 1, 0.2, 0.2), mgp=c(1.8,0.5,0), cex.lab=1, cex.axis=1, pch=21, cex=1.5, lwd=3, tck=0.015)

plot(nmds1, type = "n")
#points(nmds1, display = "sites",pch = 19, cex =1.2, col="gray")
#points(nmds1, display = "species", col=col, pch = 19, cex = 0.8)
orditorp(nmds1,display="species", col=col,air=0.0001)
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



# Assign colours for the traits:

names(comp)

col= c("blue", "blue", "blue", "blue","blue",
         "limegreen", "limegreen",  
       "orange", "orange", "orange", "orange",
       "red", "red", "red", "red","red","red","red","red",
       "darkmagenta", "darkmagenta", "darkmagenta", "darkmagenta","darkmagenta","darkmagenta","darkmagenta",
       "darkcyan", "darkcyan", "darkcyan", "darkcyan","darkcyan","darkcyan","darkcyan","darkcyan", "darkcyan"
      )


p1  <- ggplot(data = sites.scores, 
                aes(x = NMDS1, y = NMDS2)) +
  geom_vline(xintercept = 0, col="grey", linetype="dashed")+
  geom_hline(yintercept = 0, col="grey", linetype="dashed")+
  scale_shape_manual(values=c( 19, 1))+ 
  geom_point (data = species.scores, size = 2, pch=19, colour= col) + #, colour=Col, fill=Col
  geom_text_repel (data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), size = 3.6, colour= col) + 
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
  geom_point (data = species.scores, size = 2, pch=19, colour= col) + 
  geom_text_repel (data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), size = 3.6, colour= col) + 
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

coord_cont <- as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit, rescale=F)


coord_cont_ <- coord_cont %>% 
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
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = coord_cont_, linewidth =0.7, alpha = 1, 
               colour = "black", 
               arrow = arrow(length = unit(0.03, "npc")))+
#
  scale_shape_manual(values=c( 19, 1))+ 
  geom_point (data = species.scores, size = 2, pch=19, colour= col) + 
  geom_text_repel (data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), size = 3.6, colour= col) + 
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

p3

# add labels
p3 + geom_text(data = coord_cont_, 
            aes(x = NMDS1, y = NMDS2, label = Variables), 
            colour = "black", fontface = "bold", size = 4 ,
            vjust=0, hjust=0 # adjust text positions
            )