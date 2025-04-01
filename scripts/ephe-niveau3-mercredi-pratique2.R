#---------------------------------------------------#
#### EPHE : Analyse de données - niveau avancé ####
# modèle sur données spatiales - version "difficile"#
# modèles de processus ponctuels appliqués à l'écologie
#---------------------------------------------------#

setwd("C:/Users/jeany/OneDrive/Documents/EPHE_enseignement/Statistiques/2024/niveau3/supports/donnees")

library(sf)
library(spatstat)
library(terra)
library(RColorBrewer)
library(ggmap)
library(sp)
library(raster)

#--------------------------#
###  variable de réponse ###
#--------------------------#

# données opportunistes de Timon lepidus (Lézard ocellé) - coordonnées Lambert 93
rept=read.csv2("lezard_ocelle_occitanie.csv", dec = ".")

# emprise de l'analyse (= région Occitanie)
emprise=read_sf(dsn="ppm_reptiles_occitanie_emprise.shp")

#creation de la fenêtre d'analyse
emprise.sp <- as_Spatial(emprise)
fen_occit <- as.owin(emprise)


# création d'un patron de points ; on étiquette chaque point avec l'espèce correspondante
rept.ppp=ppp(rept$x, rept$y, fen_occit) 
rept.ppp=rjitter(rept.ppp, 0.01)

#-----------------#
### exploration ###
#-----------------#

# données brutes
plot(rept.ppp,main="Lézard ocellé")

# kernel de l'intensité (= espérance) du patron
plot(density(rept.ppp))

# K de Ripley (long, à éviter dès qu'il y a beaucoup de points - ici, autour d'une dizaine de minutes)
#K=Kest(rept.ppp)
# plot(K) il y a clairement un patron spatial résiduel

# comptage de points par quadrats - autre diagnostic plus simple et donc plus rapide, mais plus brut
Q=quadratcount(rept.ppp, nx = 15, ny = 15) #nx et ny donnent le nombre de quadrats en lignes et colonnes
plot(Q)
Q

#-----------------#
### covariables ###
#-----------------#

#--- la fonction as.iml.SpatRaster1 est une fonction provisoire faisant suite 
# à la disparition du package rgdal et en attendant une implémentation dans spatstat
# origine : https://stackoverflow.com/questions/77912041/convert-raster-terra-to-im-object-spatstat


as.im.SpatRaster1 <- function(X) {
  X <- X[[1]]
  rs <- terra::res(X)
  e <- as.vector(terra::ext(X))
  out <- list(
    v = as.matrix(X, wide=TRUE)[nrow(X):1, ],
    dim = dim(X)[1:2],
    xrange = e[1:2],
    yrange = e[3:4],
    xstep = rs[1],
    ystep = rs[2],
    xcol = e[1] + (1:ncol(X)) * rs[1] + 0.5 * rs[1],
    yrow = e[4] - (nrow(X):1) * rs[2] + 0.5 * rs[2],
    type = "real",
    units  = list(singular=units(X), plural=units(X), multiplier=1)
  )
  attr(out$units, "class") <- "unitname"
  attr(out, "class") <- "im"
  out
}


urb=raster("urbain_1000.tif")
urb = rast(urb)
urb=as.im.SpatRaster1(urb)
urb=(urb-mean(urb))/sd(urb) # le centrage réduction n'est pas impératif mais ici, on veut comparer les effets des différentes variables


temp=raster("temperature_moyenne.tif")
temp = rast(temp)
temp= as.im.SpatRaster1(temp)
temp=(temp-mean(temp))/sd(temp)

#------------#
### modèle ###
#------------#

# le modèle
mod.rept0=ppm(rept.ppp,~urb+temp, eps=50)

# résidus par espèce
resid.sp=residuals(mod.rept0)
plot(resid.sp,main="Lézard ocellé, résidus et localisation des données",pch=21,cex=0.5) 
# on voit une forte structure spatiale persister dans les résidus

# un autre plot de diagnostic sur les résidus qui montre bien les structures
diagnose.ppm(mod.rept0)

# corrélogramme
#qqplot.ppm(mod.rept0) # jusqu'à plusieurs heures de temps de run même pour un modèle simple

# inférence
summary(mod.rept0$internal$glmfit)
res.temp=parres(mod.rept0,covariate="temp")
plot(res.temp) # attention à l'interprétation des smoothers, tendance au surlissage

res.urb=parres(mod.rept0,covariate="urb")
plot(res.urb) 

#---------------------------------------------------------#
### estimation de l'autocorrélation spatiale résiduelle ###
#---------------------------------------------------------#

# estimation de la valeur du paramètre spatial
# !!! cette partie est essentielle pour bien identifier les paramètres de l'autocorrélation spatiale
# mais elle prend beaucoup de temps = plusieurs heures à plusieurs jours
# dans cet exemple précis ça prend environ 5 minutes....
# ici on aide un peu en faisant l'hypothèse que la structure spatiale se situe entre 0 et 1000m, avec un pas de 50m, 
# et une saturation entre 1 et 40 (voir doc du package ou Baddeley et al. pour détails)
# en tout la fonction compare 800 modèles pour voir lequel a la meilleure pseudo vraissemblance

# NB : il faut refaire tourner la fonction à chaque fois qu'on change la structure de covariables... réfléchissez bien le modèle avant!
#rs=expand.grid(r=seq(1,1001, by=50), sat=1:40)
#term.inter=profilepl(rs, 
Geyer, rept.ppp,~urb+temp)
# on obtient  r = 1001 and sat = 7


#----------------------------------------------------------#
### modèle avec correction de l'autocorrélation spatiale ###
#----------------------------------------------------------#

# eps=50 permet de rajouter des points de pseudo absence; ils n'ont aucun rôle si ce n'est d'estimer correctement
# la pseudo vraisemblance
# compter ~ 40 secondes à 1 min de run, augmente très vite si on augmente le nb d'espèces ou de variables
t1=Sys.time()
mod.rept=ppm(rept.ppp,~urb+temp, interaction = Geyer(r=1001,sat=7), eps=50)
t2=Sys.time()
t2-t1 #  temps de run (facultatif, mais utile pour les modèles longs afin de planifier l'analyse)

# résidus par espèce
resid.sp=residuals(mod.rept)
plot(resid.sp,main="Lézard ocellé, résidus et localisation des données",pch=21,cex=0.5) 
# il n'y  a plus de patron résiduel - on est même probablement surlissé.

plot(Smooth(resid.sp)) # avec une autre mise à l'échelle et un lissage par smoother : on voit à nouveau une structure résiduelle spatiale
# mais on est sur des résidus de l'ordre de 10e-08, donc pas la peine d'insister. En revanche cette structure porte une information (en particulier
# les résidus sont structurés au niveau des Pyrénées, ce qui suggère qu'il y a peut être des déterminants environnementaux un peu différents dans cette 
# zone

# inférence
summary(mod.rept$internal$glmfit) 
# l'effet de l'urbanisation est beaucoup plus faible, alors que l'effet de l'altitude reste essentiellement le même

# résidus partiels
res.temp=parres(mod.rept,covariate="temp")
res.urb=parres(mod.rept,covariate="urb")
par(mfrow=c(1,2))
plot(res.temp)
plot(res.urb)

#----------------#
### prédiction ###
#----------------#

# on veut la prédiction de la distribution de l'espèce selon le modèle
mod.rept.pred=predict(mod.rept)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))  # nombreuses autres palette, celle-ci donne des résultats assez jolis
plot(mod.rept.pred, col=h)

# la prédiction de la distribution de l'espèce selon la température uniquement --> on force l'urbanisation à être constante
coef=mod.rept$internal$glmfit$coefficients
coef_modif=coef
coef_modif["urb"]=0
mod.rept.urbctt=mod.rept
mod.rept.urbctt$internal$glmfit$coefficients=coef_modif
mod.rept.urbctt.pred=predict(mod.rept.urbctt)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))
plot(mod.rept.urbctt.pred, col=h, equal.ribbon=T) # on trouve une espèce essentiellement de plaine littorale, comme attendu

# la prédiction de la distribution de l'espèce selon l'urbanisation uniquement --> on force l'altitude à être constante
coef=mod.rept$internal$glmfit$coefficients
coef_modif=coef
coef_modif["temp"]=0
mod.rept.altictt=mod.rept
mod.rept.altictt$internal$glmfit$coefficients=coef_modif
mod.rept.altictt.pred=predict(mod.rept.altictt)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))
plot(mod.rept.altictt.pred, col=h) # on trouve une présence moindre dans les grosses agglomérations

# ces prédictions semblent relativement indépendantes de la localisation des données, ce qui suggère que le modèle ne surlisse pas trop
plot(mod.rept.pred, col=h, equal.ribbon=T)
points(rept.ppp)

# par contre, l'exploration des données et notre connaissance du processus de génération devrait nous alerter. Regardez en particulier 
# la disposition des points sur le littoral, qui forment une ligne le long de l'autoroute... Il y a aussi des regroupements de points localisés à des endroits précis
# résidus ou pas, il va falloir en tenir compte

#------------------------------#
### effort d'échantillonnage ###
#------------------------------#

# exploration de la disposition des points
coordinates(rept)=c("xcoord", "ycoord")
proj4string(rept) <- CRS("+init=epsg:2154")
rept2=spTransform(rept, CRS("+proj=longlat +datum=WGS84"))
rept.df=as.data.frame(rept2)
qmplot(x=xcoord, y=ycoord, data = rept.df) #, colour = I('red'), size = I(3), darken = .3)

rept3=data.frame(slot(rept2, "coords"), slot(rept2, "data"))
names(rept3)[1:2]=c("lon","lat")
qmplot(lon, lat, data = rept3) # on voit très bien l'alignement des points sur l'autoroute A9, les déviations de Narbonne et de Nimes,
# et un certain nombre d'autres axes routiers

# nombre de dates uniques de saisie de données
date=raster("biais_date_250.tif")
date=rast(date)
date = as.im.SpatRaster1(date) 
date=(date-mean(date))/sd(date) 
plot(date)

# nombre de données générées par des bureaux d'études
be=raster("biais_bureau_250.tif")
be<-rast(be)
be <- as.im.SpatRaster1(be)
be=(be-mean(be))/sd(be)
plot(be) # variable écrasée par quelques grosses concentrations de données

# éléments linéaires (routes et voies ferrées)
elim=raster("lineaire_1000.tif")
elim <- rast(elim)
elim=as.im.SpatRaster1(elim)
lineaire_1000=(elim-mean(elim))/sd(elim)
plot(lineaire_1000)

# on refait le modèle avec ces nouvelles variables (! pour bien faire il faudrait aussi réestimer le terme d'interaction spatial)
t1=Sys.time()
mod.rept.ech=ppm(rept.ppp,~urb+temp+date+be+lineaire_1000, interaction = Geyer(r=1001,sat=7), eps=50)
t2=Sys.time()
t2-t1 

# résidus
resid.sp.ech=residuals(mod.rept.ech)
plot(resid.sp.ech,main="Lézard vert, résidus et localisation des données",pch=21,cex=0.5) 
plot(Smooth(resid.sp.ech)) # les variables d'échantillonnage n'ont pas un gros effets sur les résidus (probable qu'ils soient dominés par le terme d'autocorrélation)

# inférence
summary(mod.rept)
summary(mod.rept.ech) # on a  réaugmenté l'impact des zones urbaines, malgré l'autocorrélation

res.temp.ech=parres(mod.rept.ech,covariate="temp")
res.urb.ech=parres(mod.rept.ech,covariate="urb")

par(mfrow=c(2,2))
plot(res.temp,main="température, sans correction",ylim=c(-4,2))
plot(res.urb,main="urbanisation, sans correction",ylim=c(-4,2))
plot(res.temp.ech,main="température, avec correction",ylim=c(-4,2))
plot(res.urb.ech,main="urbanisation, avec correction",ylim=c(-4,2))

#-----------------------------------------------------#
### prédictions à effort d'échantillonnage constant ###
# (à comparer avec le modèle sans effet échantillonnage) #
#-----------------------------------------------------#

# on veut la prédiction de la distribution de l'espèce selon le modèle, en homogénéisant l'effort pour un effort moyen
coef=mod.rept.ech$internal$glmfit$coefficients
coef_modif=coef
coef_modif["date"]=mean(date)
coef_modif["be"]=mean(be)
coef_modif["elim"]=mean(lineaire_1000)

mod.rept.effequal=mod.rept.ech
mod.rept.effequal$internal$glmfit$coefficients=coef_modif
mod.rept.effequal.pred=predict(mod.rept.effequal)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))
plot(mod.rept.effequal.pred, col=h, equal.ribbon=T,main="") 

# pour rappel, la même prédiction pour le modèle sans effort d'échantillonnage
mod.rept.pred=predict(mod.rept)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))  
x11()
plot(mod.rept.pred, col=h, equal.ribbon=T,main="") # les deux prédictions se ressemblent beaucoup, mais il y a des différences locales

#------------------------------------------------------------------#
### effet de chaque variable d'échantillonnage sur la prédiction ###
#------------------------------------------------------------------#

# effet date
coef=mod.rept.ech$internal$glmfit$coefficients
coef_modif=coef
coef_modif["date"]=mean(date)
mod.rept.effdate=mod.rept.ech
mod.rept.effdate$internal$glmfit$coefficients=coef_modif
mod.rept.effdate.pred=predict(mod.rept.effdate)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))
plot(mod.rept.effdate.pred, col=h) 

# effet bureaux d'études
coef=mod.rept.ech$internal$glmfit$coefficients
coef_modif=coef
coef_modif["be"]=mean(be)
mod.rept.effdate=mod.rept.ech
mod.rept.effdate$internal$glmfit$coefficients=coef_modif
mod.rept.effdate.pred=predict(mod.rept.effdate)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))
plot(mod.rept.effdate.pred, col=h) 

# effet éléments linéaires
coef=mod.rept.ech$internal$glmfit$coefficients
coef_modif=coef
coef_modif["lineaire_1000"]=mean(lineaire_1000)
mod.rept.effdate=mod.rept.ech
mod.rept.effdate$internal$glmfit$coefficients=coef_modif
mod.rept.effdate.pred=predict(mod.rept.effdate)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu")))
plot(mod.rept.effdate.pred, col=h) 

# c'est quand on fixe l'effet "bureaux d'études" qu'on obtient la plus grande variabilité
# dans la distribution prédite --> cette variable semble avoir le plus fort impact sur la 
# prédiction. lorsqu'on fixe l'une ou l'autre des deux autres variables, la prédiction
# devient quasiment plate sur toute l'emprise, ce qui suggère qu'elle est écrasée par
# la variation liée aux BE

#--------------------------#
### modèle multi-espèces ###
#--------------------------#

# données d'occurrences toutes espèces
dat.multsp=read.csv2("occurrences_reptiles_occitanie.csv",  dec=".", header=T, fill=T, stringsAsFactors =T)

# nombreuses espèces avec des effectifs variables (on a supprimé quelques espèces jugées trop rares)
summary(dat.multsp$species)

# processus ponctuel multiespèces : notez l'argument "marks"
pp.multsp=ppp(dat.multsp$x,dat.multsp$y,fen_occit,marks = dat.multsp$species)
pp.multsp=rjitter(pp.multsp, 0.01) # on permet un décalage très mineur des points afin d'éviter que des points superposés ne soient considérés comme doublons
pp.multsp_sp=split(pp.multsp) # on crée une table pour chaque espèce, sera utile plus tard

# le ppm (on corrige d'emblée pour les variables d'échantillonnage)
t1=Sys.time()
ppm.multsp=ppm(pp.multsp,~marks*urb+
                          marks*temp+
                        date+
                          be+
                          marks*lineaire_1000,
                              interaction = Geyer(r=980,sat=48), eps=50)
t2=Sys.time()

# résumé du modèle
summary(ppm.multsp$internal$glmfit)

# effets : se référer au résumé du modèle et reconstruire les courbes avec predict()
# les résidus partiels ne sont pas encore implémentés pour les modèles "marqués" 

# prédictions par espèces
coef=ppm.multsp$internal$glmfit$coefficients
coef_modif=coef
coef_modif["date"]=0
coef_modif["be"]=0
coef_modif["lineaire_1000"]=0
ppm.multsp.pred=ppm.multsp
ppm.multsp.pred$internal$glmfit$coefficients=coef_modif
predict_model=predict(ppm.multsp.pred)
h=colorRampPalette(rev(brewer.pal(10 , "RdYlBu"))) 
plot(predict_model, col=h) #equal.ribbon=T permet de mettre la même échelle pour toutes les espèces
par(mfrow=c(1,1)) # ici on a supprimé equal.ribbon car grosses différences entre des espèces ubiquistes vs rares

# exploration des résidus
par(mfrow=c(1,2))
res_sp<-split(residuals(ppm.multsp))
plot(res_sp$PODMUR);plot(Smooth(res_sp$PODMUR))
plot(res_sp$CHASTR);plot(Smooth(res_sp$CHASTR))
plot(res_sp$LACBIL);plot(Smooth(res_sp$LACBIL))
plot(res_sp$MALMON);plot(Smooth(res_sp$MALMON))
plot(res_sp$NATMAU);plot(Smooth(res_sp$NATMAU))
plot(res_sp$PSAALG);plot(Smooth(res_sp$PSAALG))
plot(res_sp$PSAEDW);plot(Smooth(res_sp$PSAEDW))
plot(res_sp$RHISCA);plot(Smooth(res_sp$RHISCA))
plot(res_sp$TARMAU);plot(Smooth(res_sp$TARMAU))
plot(res_sp$TIMLEP);plot(Smooth(res_sp$TIMLEP))
plot(res_sp$ZAMLON);plot(Smooth(res_sp$ZAMLON))
plot(res_sp$NATNAT);plot(Smooth(res_sp$NATNAT))
plot(res_sp$PODLIO);plot(Smooth(res_sp$PODLIO))
plot(res_sp$ANGFRA);plot(Smooth(res_sp$ANGFRA))

#------------------------------------#
### inférence sur les coefficients ###
#------------------------------------#

# comparaison multi-espèces
coefs=coef(ppm.multsp)
sp.eff=coefs[1:14] # les effets simples "espèce"
urb.eff=coefs[grep("urb",names(coefs))] # effets de l'urbanisation (!! ce sont des contrastes)
urb.eff["urb"]=0 # on va rester en contrastes pour simplifier
names(urb.eff)=levels(dat.multsp$species) # on renomme les coefficients avec les noms d'espèces
