# Code for quantification of functional composition and functional diversity using functional traits

library(dplyr)
library(tidyverse)

citation() 
rm(list=ls(all=TRUE))

# Trait data----
traits  <- read.csv ( "Data/Traits.csv", header=T, row.names = 1)

head(traits)
str(traits)
summary(traits)

names(traits)


# check if there are NAs
which(is.na(traits))
# check if there are traits with zero appearance
colSums(traits[,-c(1:2)])

### Standardize traits ----
### Standardize the fuzzy coded traits 

traits_st <- traits %>% 
  dplyr::select(-taxa) %>% 
  mutate(RESP=RESP_teg+RESP_gil+RESP_pls+RESP_spi+RESP_ves) %>% 
  mutate(RESP_teg=RESP_teg/RESP, RESP_gil=RESP_gil/RESP, RESP_pls=RESP_pls/RESP,
         RESP_spi=RESP_spi/RESP, RESP_ves=RESP_ves/RESP) %>% 
  dplyr::select(-RESP) %>% 
  mutate(LCD=LCD_short+LCD_long) %>% 
  mutate(LCD_short=LCD_short/LCD, LCD_long=LCD_long/LCD) %>% 
  dplyr::select(-LCD)%>% 
  mutate(RLC=RLC_sev+RLC_unv+ RLC_mlt + RLC_flx) %>% 
  mutate(RLC_sev=RLC_sev/RLC, 
         RLC_unv=RLC_unv/RLC,
         RLC_mlt=RLC_mlt/RLC,
         RLC_flx=RLC_flx/RLC) %>% 
  dplyr::select(-RLC)%>% 
  mutate(REP=REP_ovo+REP_fie+REP_cie+REP_fic+REP_frc+REP_vec+REP_tec+REP_ase) %>% 
  mutate(REP_ovo=REP_ovo/REP, 
         REP_fie=REP_fie/REP,
         REP_cie=REP_cie/REP,
         REP_fic=REP_fic/REP,
         REP_frc=REP_frc/REP, 
         REP_vec=REP_vec/REP,
         REP_tec=REP_tec/REP,
         REP_ase=REP_ase/REP) %>% 
  dplyr::select(-REP) %>% 
mutate(BodSz=BodSz_0.25+ BodSz_0.25_0.5+ BodSz_0.5_1+ BodSz_1_2+ BodSz_2_4+ BodSz_4_8+ BodSz_8_)%>% 
  mutate(BodSz_0.25=BodSz_0.25/BodSz, 
         BodSz_0.25_0.5=BodSz_0.25_0.5/BodSz,
         BodSz_0.5_1=BodSz_0.5_1/BodSz,
         BodSz_1_2=BodSz_1_2/BodSz,
         BodSz_2_4=BodSz_2_4/BodSz, 
         BodSz_4_8=BodSz_4_8/BodSz,
         BodSz_8_=BodSz_8_/BodSz) %>% 
  dplyr::select(-BodSz) %>% 
  mutate(Feed=Feed_Gra+Feed_Min+Feed_Xyl+Feed_Shr+Feed_Gat+Feed_Aff+Feed_Pff+Feed_Pre+Feed_Par) %>% 
  mutate(Feed_Gra=Feed_Gra/Feed,
         Feed_Min=Feed_Min/Feed,
         Feed_Xyl=Feed_Xyl/Feed,
         Feed_Shr=Feed_Shr/Feed,
         Feed_Gat=Feed_Gat/Feed, 
         Feed_Aff=Feed_Aff/Feed,
         Feed_Pff=Feed_Pff/Feed,
         Feed_Pre=Feed_Pre/Feed,
         Feed_Par=Feed_Par/Feed) %>% 
  dplyr::select(-Feed)


    

# check if the row sums for each trait is 1
rowSums(traits_st %>% 
  dplyr::select(RESP_teg, RESP_gil, RESP_pls, RESP_spi, RESP_ves))

rowSums(traits_st %>% 
          dplyr::select(LCD_short, LCD_long))

rowSums(traits_st %>% 
          dplyr::select(RLC_sev, RLC_unv, RLC_mlt, RLC_flx))

rowSums(traits_st %>% 
          dplyr::select(REP_ovo, REP_fie, REP_cie, REP_fic, REP_frc, REP_vec,  REP_tec, REP_ase))

rowSums(traits_st %>%
          dplyr::select(BodSz_0.25, BodSz_0.25_0.5, BodSz_0.5_1, BodSz_1_2, BodSz_2_4, BodSz_4_8, BodSz_8_))

rowSums(traits_st %>%
          dplyr::select(Feed_Gra, Feed_Min, Feed_Xyl, Feed_Shr, Feed_Gat, Feed_Aff, Feed_Pff, Feed_Pre, Feed_Par))


str(traits_st)

# Community data----

Sp_comp <- read.csv ( "Data/Com_composition.csv", header=T, row.names = 1)

head(Sp_comp)
str(Sp_comp)

# check dimensions of community composition matrix and of trait matrix
dim(Sp_comp)
dim(traits_st)


# check if there are NAs
which(is.na(Sp_comp))
# check if there are species with zero appearance
rowSums(Sp_comp[,-c(1)])
colSums(Sp_comp[,-c(1)])


# prepare data for the FD analysis
trait <- as.matrix (traits_st)
Com_comp <- as.matrix(Sp_comp)

# Calculate funct. composition matrix ----

## 'functcomp' takes trait matrix and abundance matrix and for each site (pond transect) computes the community-level weighted means of each trait (waited by abundance)

library(FD)

# functcomp  - function for calculation of functional composition matrix

FuncComp <- functcomp(trait, Com_comp, 
                      CWM.type = "all", 
                      bin.num=c("BodSz_8_") #  bin.num - indicates binary traits to be treated as continuous  
                      )  
FuncComp
head(FuncComp)

str(FuncComp)

# check if the row sums for each trait is 1
rowSums(FuncComp %>% 
          dplyr::select(RESP_teg, RESP_gil, RESP_pls, RESP_spi, RESP_ves))

rowSums(FuncComp %>% 
          dplyr::select(LCD_short, LCD_long))

rowSums(FuncComp %>% 
          dplyr::select(RLC_sev, RLC_unv, RLC_mlt, RLC_flx))

rowSums(FuncComp %>% 
          dplyr::select(REP_ovo, REP_fie, REP_cie, REP_fic, REP_frc, REP_vec,  REP_tec, REP_ase))

rowSums(FuncComp %>%
          dplyr::select(BodSz_0.25, BodSz_0.25_0.5, BodSz_0.5_1, BodSz_1_2, BodSz_2_4, BodSz_4_8, BodSz_8_))

rowSums(FuncComp %>%
          dplyr::select(Feed_Gra, Feed_Min, Feed_Xyl, Feed_Shr, Feed_Gat, Feed_Aff, Feed_Pff, Feed_Pre, Feed_Par))



# save functional composition matrix
write.csv(FuncComp, file = "Data/FuncComp.csv")


# Calculate funct. diversity indices -------------

FuncD <- dbFD(trait, Com_comp, calc.CWM = F, m=6) 
FuncD

# save Functional diversity indices
write.csv(FuncD, file = "Data/FuncDiv.csv")


## Correlation among functional diversity indices----

FuncDiv <- read.csv ( "Data/FuncDiv.csv", header=T)%>% 
  dplyr::select(-qual.FRic, -sing.sp)
names(FuncDiv)

Data <- FuncDiv[,2:7]

cor <- cor(Data, method = c("pearson"))
round(cor, digits=2)

write.csv(round(corl1, 2), "correl_FunDiv.csv", row.names = T)


# Correlogram:
library(corrplot)

x11(height=9,width=8.5)
corrplot(cor, type = "upper",  
         tl.col = "black", tl.srt = 50)

