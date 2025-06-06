#-----------------------------------------------------------------#
### Formation "Certificat en analyse de donn�es pour les �cologues"
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# version 2023
# Cas pratique : mod�les mixtes
#-----------------------------------------------------------------#

rm(list=ls())
setwd("F:/certificat_2023/niveau3/supports/donnees_final")

# packages utiles
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(cowplot)
library(lme4)
library(MuMIn)
library(visreg)

#-------------------#
#### T�tras lyre ####
#-------------------#

# donn�es
tetras=read.table("comptages_tetras_lyre.txt",header=T,sep="\t")

# cr�ation des indices de sites UN, RN, RG pour les sorties graphiques
tetras$UN=factor(paste("ID",tetras$UN,sep=""))
tetras$RN=factor(paste("ID",tetras$RN,sep=""))
tetras$RG=factor(paste("ID",tetras$RG,sep=""))

# on v�rifie la disponibilit� des donn�es par niveau de facteur
table(tetras$UN)
table(tetras$RN)
table(tetras$RG) # attention au niveau r�gions, qui n'a que deux niveaux

table(unique(tetras[,c("UN","RN")])$RN)

table(unique(tetras[,c("RN","RG")])$RG)

table(unique(tetras[,c("UN","RG")])$RG)

# mod�le 1 : effets al�atoires UN sur la pente et sur l'intercept
logpoules=log(tetras$Poules)
mod.tetras=glmer(Jeunes~NAOdjfm+I(NAOdjfm^2)+(1+NAOdjfm+I(NAOdjfm^2) |UN),family="poisson",data=tetras,offset=logpoules)

# NB : si un message "boundary (singular) fit sort, ce n'est pas inqui�tant : c'est classique dans un mod�le 
# complexe comme celui ci, et est li� � la tendance de lme4 � sur-reporter des erreurs. Ici, la raison de ce message est que vu
# l'�chelle de variation de certaines variables, on se retrouve en bordure de surface de vraisemblance. Cela n'affecte pas les 
# param�tres

# NB2 : si vous avez bien fait le travail, vous avez commenc� par ajuster un mod�le de Poisson sans l'effet UN, puis avec juste un 
# effet al�atoire UN sur l'intercept, puis avec les pentes al�atoires. Nous n'avons pas mis ces �tapes dans le script afin de ne pas
# le surcharger, mais elles sont indispensables. Ne cherchez pas � comparer ces diff�rents mod�les - il existe des outils de comparaison
# mais ce n'est pas l'objectif : ces premiers mod�les vous servent surtout � v�rifier la stabilit� des param�tres.

# r�sidus
plot(mod.tetras)

# on va chercher � comprendre le point tout en haut
tetras[which(resid(mod.tetras,type="pearson")>5),]
# il s'agit manifestement d'une UN o� la reproduction est particuli�rement �lev�e

# R�
r.squaredGLMM(mod.tetras)
# on constate la forte contribution de l'effet al�atoire

# effets
res.p=ggemmeans(mod.tetras,terms="NAOdjfm [all]")
plot(res.p,residuals=T)+ylim(0,100)

# pour repr�senter les effets al�atoires
plot_model(mod.tetras,type="eff",terms="NAOdjfm")
plot_model(mod.tetras,type="pred",terms=c("NAOdjfm [all]","UN"),pred.type="re")

# on rescale
plot(res.p,residuals=T)+ylim(0,15)

# mod�le 2
mod.tetras2=glmer(Jeunes~NAOdjfm+I(NAOdjfm^2)+(1+NAOdjfm+I(NAOdjfm^2) |RG/RN/UN),family="poisson",data=tetras,offset=logpoules)

r.squaredGLMM(mod.tetras2)
# on ne gagne pas grand chose avec l'embo�tement d'effets al�atoires : ce n'est pas tr�s �tonnant vu le peu d'information qu'ils contiennent

res.p=ggemmeans(mod.tetras2,terms="NAOdjfm [all]")
plot(res.p,residuals=T)+ylim(0,100)

plot_model(mod.tetras2,type="eff",terms="NAOdjfm")
plot_model(mod.tetras2,type="pred",terms=c("NAOdjfm [all]","UN"),pred.type="re")
# le mod�le commence � �tre vraiment complexe, on n'arrive plus � afficher les effets al�atoires embo�t�s

# autre vue de l'effet fixe NAO avec les graphiques R de base
library(effects)
ef=effect("NAOdjfm",mod.tetras)
plot(ef)

### comparaison avec le mod�le sans effets al�atoires ###
mod.tetras3=glm(Jeunes~NAOdjfm+I(NAOdjfm^2),family="poisson",data=tetras,offset=logpoules)
plot_model(mod.tetras3,type="eff",terms="NAOdjfm")

p1=plot_model(mod.tetras2,type="eff",terms="NAOdjfm")
p2=plot_model(mod.tetras3,type="eff",terms="NAOdjfm")
cowplot::plot_grid(p1,p2)

### optimum climatique ###

# mod�le mixte
y=function(x){0.38795+0.08736*x-0.08166*x*x}
xmax.mix=optimize(y, interval=range(tetras$NAOdjfm), maximum=TRUE)
xmax.mix

# glm
y=function(x){0.46064+0.04435*x-0.08280*x*x}
xmax.glm=optimize(y, interval=range(tetras$NAOdjfm), maximum=TRUE)
xmax.glm


