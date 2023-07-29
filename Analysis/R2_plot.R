# Code for producing table and table chart for partial R2 for all predictors and responce variables

# packages----

library(tidyverse)
library(ggplot2)

# Data----

# read data for the R2s for all response variables
Abund_R2 <-read_csv("Results/Abund_R2.csv")%>% dplyr::select(Variables, Rsq) %>% rename(Abund=Rsq)
SR_R2 <-read_csv("Results/SR_R2.csv")%>% dplyr::select(Variables, Rsq) %>% rename(SR=Rsq)
FEve_R2 <-read_csv("Results/FEve_R2.csv")%>% dplyr::select(Variables, Rsq) %>% rename(FEve=Rsq)
FDiv_R2 <-read_csv("Results/FDiv_R2.csv")%>% dplyr::select(Variables, Rsq) %>% rename(FDiv=Rsq)
FDis_R2 <-read_csv("Results/FDis_R2.csv")%>% dplyr::select(Variables, Rsq) %>% rename(FDis=Rsq)
FuncComp_R2 <-read_csv("Results/FuncComp_R2.csv")%>% dplyr::select(Variables, R2) %>% rename(FComp=R2)

# combine the R2s for all response variables

R2_all <- Abund_R2 %>% 
  full_join(SR_R2, by="Variables") %>% 
  full_join(FEve_R2, by="Variables") %>% 
  full_join(FDiv_R2, by="Variables") %>% 
  full_join(FDis_R2, by="Variables") %>% 
  full_join(FuncComp_R2, by="Variables") %>% 
  mutate(Variables=fct_relevel(Variables, c("Agric_m", "Highway_m", "Gravel", "Waste", 
                                            "Water_PC1", "Water_PC2", "Water_PC3",
                                            "Macrophytes", "Fish"))) %>% 
  arrange(Variables) %>% 
  mutate(Drivers=c("Human impact", "Human impact", "Human impact", "Human impact", 
                   "Water properties", "Water properties", "Water properties",
                   "Macrophyte biomass", "Fish predation"))

R2_all
# save
write.csv(R2_all , file="Results/R2_all.csv", row.names=TRUE)


# Prepare data for bubble chart 
R2 <- R2_all %>%  
  group_by(Drivers) %>% 
  summarise(Abund=sum(Abund, na.rm=T), SR=sum(SR, na.rm=T),
            FEve=sum(FEve, na.rm=T), FDiv=sum(FDiv, na.rm=T),   
            FDis=sum(FDis, na.rm=T),FComp =sum(FComp, na.rm=T)) %>% 
  mutate(Drivers = fct_relevel(Drivers, c("Fish predation", "Macrophyte biomass", 
                                        "Water properties", "Human impact"))) %>% 
  arrange(Drivers) %>% 
  pivot_longer(!Drivers, names_to="Diversity", values_to = "R2") %>%
    mutate(Diversity = fct_relevel(Diversity, c("Abund","SR",
                                              "FEve","FDiv","FDis",
                                              "FComp")))
 

# Plot---- 
# Bubble chart 

ggplot(R2, aes(x = Diversity, y = Drivers, size = R2, colour = Drivers, fill = Drivers)) +
  geom_point(pch=21) +
 # scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(1, 22)) + 
   # scale_color_manual(values = c("Human impact" = "red", "Water properties" = "blue", "Fish predation" = "#FDBF6F" , "Macrophyte biomass" = "#33A02C"))+ 
   #  scale_color_brewer(palette =  "Paired") +
  labs(x = NULL, y = NULL, size=R^2~parital) + 
  theme(legend.position = "right",
        legend.key = element_rect(colour = NA, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_line(colour = "grey80"),
        axis.ticks = element_blank(),
        axis.text.y=element_text(colour = "black", size=15),
        axis.text.x=element_text(colour = "black", size=15))+  
  guides(color = "none", fill="none") # remove fill in a legend
