#-------------------------------------------------#
#### EPHE : Analyse de données - niveau avancé ####
# modèle sur données spatiales - version "facile"#
#-------------------------------------------------#

setwd("F:/certificat_2023/niveau3/supports/donnees_final")

library(ggplot2)
library(mapview)
library(sf)
library(patchwork)
library(tidyverse)
library(PerformanceAnalytics)
library(mgcv)
webshot::install_phantomjs() # pour les cartes
library(gratia) # une solution parmi d'autres pour les sorties graphiques de GAM
library(mgcViz) # une autre solution pour les sorties de GAM
library(spdep) # une solution parmi d'autres pour les corrélogrammes spatiaux

# voir la vignette mgcViz : https://cran.r-project.org/web/packages/mgcViz/vignettes/mgcviz.html

#---------------#
#### données ####
#---------------#

chevreuils3 = read.csv2("OFB_chevreuils2017.csv",row.names=1)

#-------------------#
#### exploration ####
#-------------------#

# répartition des données entre massifs
summary(chevreuils3)
table(chevreuils3$Nom)

# abroutissement (variable binaire)
sum(chevreuils3$Abroutissement)

# carte de l'abroutissement
chevreuils3$fabr = factor(chevreuils3$Abroutissement)
carte.abroutissement = chevreuils3%>%
  mapview(zcol = "fabr",
          xcol="X",ycol="Y",
          legend = T,
          map.types = "OpenStreetMap",
          layer.name = "Abroutissement",
          crs=2154)

carte.abroutissement

# covariables

dl= ggplot(chevreuils3)+
  aes(x=Densite_lin)+
    geom_histogram()

freq= ggplot(chevreuils3)+
  aes(x=Frequentation+1)+
  geom_histogram()

alt= ggplot(chevreuils3)+
  aes(x=Altitude)+
  geom_histogram()

pente= ggplot(chevreuils3)+
  aes(x=Pente)+
  geom_histogram()

Tirs= ggplot(chevreuils3)+
  aes(x=Tirs)+
  geom_histogram()

(dl + freq + alt) / (pente + Tirs) 

# on va log-transformer la fréquentation (trop asymétrique)
chevreuils3$lfrequentation = log(chevreuils3$Frequentation+1)

# carte des covariables

carte.freq = chevreuils3%>%
  mapview(zcol = "lfrequentation",
          xcol="X",ycol="Y",
          legend = T,
          map.types = "OpenStreetMap",
          layer.name = "log(Fréquentation",
          crs=2154)

carte.freq

# matrice de corrélation des covariables
chart.Correlation(chevreuils3[,c("Tirs","lfrequentation","Densite_lin","Altitude")]) 

# on pourra se préoccuper de la relation fréquentation - densité de chemins, mais a priori on peut séparer les effets

#-----------------#
#### le modèle ####
#-----------------#

# on commence simple
m1= glm(Abroutissement~Densite_lin + lfrequentation + Altitude + Pente + Tirs + Nom,family=binomial,data=chevreuils3)

# résidus
par(mfrow=c(2,2))
plot(m1) # ces résidus laissent imaginer qu'on explique peu de variance avec ce modèle

# pour avoir des relations plus flexibles : un GAM
m2= gam(Abroutissement~s(Densite_lin) + s(lfrequentation) + s(Altitude) + s(Pente) + s(Tirs) + Nom,family=binomial,data=chevreuils3)

plot(predict(m2),residuals(m2))
gam.check(m2)

# on peut penser que les relations sont variables par massif
chevreuils3$Nom = factor(chevreuils3$Nom)
m3= gam(Abroutissement~ 
          Nom + 
          s(Densite_lin,by = Nom) + 
          s(lfrequentation,by = Nom) + 
          s(Altitude,by = Nom) + 
          s(Pente,by = Nom) + 
          s(Tirs,by = Nom) ,family=binomial,data=chevreuils3)

plot(predict(m3),residuals(m3))
gam.check(m3)

# est-ce que l'augmentation de complexité est justifiée au regard des données
AIC(m1,m2,m3)

# reste-t-il un patron spatial?
chevreuils3$resid.gam = residuals(m3)

carte.resid = chevreuils3%>%
  mapview(zcol = "resid.gam",
          xcol="X",ycol="Y",
          legend = T,
          map.types = "OpenStreetMap",
          layer.name = "résidus du GAM",
          crs=2154)

carte.resid

# y-a-t-il une autocorrélation spatiale résiduelle? (attention, long! plus d'une dizaine de minutes)
xy.chev = chevreuils3[,c("X","Y")]
dist.chevreuil = dnearneigh(xy.chev,d1=0,d2=2500,longlat=F) # qui est voisin de qui (d2 définit une distance max. à laquelle deux points peuvent être voisins)
plot(dist.chevreuil,coords=xy.chev)

cg.chev = sp.correlogram(dist.chevreuil,var = chevreuils3$resid.gam,method="I",order=8,zero.policy=T)
plot(cg.chev) # il y a un léger patron d'autocorrélation spatiale résiduelle

# on tente de corriger un peu la variabilité spatiale
m4= gam(Abroutissement~ 
          Nom + 
          s(Densite_lin,by = Nom) + 
          s(lfrequentation,by = Nom) + 
          s(Altitude,by = Nom) + 
          s(Pente,by = Nom) + 
          s(Tirs,by = Nom) +
          s(X,Y),
        family=binomial,data=chevreuils3)

AIC(m3,m4) # on va garder le terme spatial (il y a de meilleures manières de prendre en compte des patrons spatiaux)

# NB : on a choisi d'imposer la même fonction spatiale pour tous les massifs simplement pour gagner du temps car 
# le modèle à interaction (X,Y)*massif est très long à faire tourner et ne s'ajuste pas correctement - mais dans un
# cas réel on le testerait probablement,ou on se tournerait vers des fonctions spatiales différentes

chevreuils3$resid.gam.spat = residuals(m4)

carte.resid.spat = chevreuils3%>%
  mapview(zcol = "resid.gam.spat",
          xcol="X",ycol="Y",
          legend = T,
          map.types = "OpenStreetMap",
          layer.name = "résidus du GAM",
          crs=2154)

carte.resid.spat

#---------------------------------#
#### exploration des résultats ####
#---------------------------------#

# le résumé du modèle
summary(m4)

# sorties graphiques avec gratia
appraise(m4)

draw(m4,residuals=T, select = smooths(m4)[1])
draw(m4,select = smooths(m4)[21], contour = FALSE, n = 50)

# courbes par massifs (à refaire pour chaque variable)
draw(m4,select=c("s(Densite_lin):NomHAUTES BAUGES",
                 "s(Densite_lin):NomSEMNOZ",
                 "s(Densite_lin):NomSW BAUGES")
                 ,wrap=T ,residuals=T,
     scales="free")

# différence des courbes entre massifs
draw(difference_smooths(m4,smooth="s(Densite_lin)"))

# sorties graphiques avec mgcViz

viz.m4 = getViz(m4)

plot(sm(viz.m4, 21)) +
  l_fitRaster() +
  l_fitContour() + 
  l_points()

plot(viz.m4,select=1)
