#-------------------------------------------------------------------------------#
### Formation "Certificat en analyse de données pour les écologues"
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# version 2023
# ce script permet la réplication des exemples du cours sur les modèles bayésiens
#-------------------------------------------------------------------------------#

rm(list=ls())

library(ggplot2)
library(ggeffects)
library(cowplot)
library(rjags)
library(PerformanceAnalytics)
library(coda)
library(bayesplot)

setwd("F:/certificat_2023/niveau3/supports/donnees_final")

#-------------------------------------#
### données baleines de B. Bouchard ###
### CEFE
#-------------------------------------#

d=read.table("baleines_Bouchard.txt",header=T,sep="\t")

#----------------------------#
### on explore les données ###
#----------------------------#

par(mfrow=c(2,2))
boxplot(d$time_z12~d$chemical)
hist(d$time_z12)
hist(log(d$time_z12+1))
dev.off()

#-----------------------------------------------#
### coder une régression linéaire en bayésien ###
#-----------------------------------------------#

## le modèle qu'on ferait en fréquentiste
modreg=lm(time_z12~chemical+log(nb_indiv)+place,data=d)

# on a transformé nb_indiv pour réduire l'échelle des valeurs sur cette variable

# exploration rapide
par(mfrow=c(2,2))
plot(modreg)
summary(modreg)
p1=plot(ggemmeans(modreg,terms="chemical"),residuals=T)
p2=plot(ggemmeans(modreg,terms="nb_indiv"),residuals=T)
p3=plot(ggemmeans(modreg,terms="place"),residuals=T)
plot_grid(p1,p2,p3)

# le même modèle en langage BUGS: script baleines_modele.R

## préparer les données: revient à mettre les données dans le format spécifié dans le script du modèle

time_z12=as.vector(d[,"time_z12"]) # on a spécifié un vecteur pour la variable de réponse
Lnb_indiv=as.vector(log(d[,"nb_indiv"])) # un autre pour le nombre d'individus, qu'on a log-transformé

chemical1=factor(d$chemical,levels=c("CTL","CLAY","DMS","KRILL")) # on ordonne les niveaux des variables facteurs de manière judicieuse (ici les deux contrôles en premier)
chemical=as.numeric(chemical1) # on doit transformer les variables facteur en variables numériques. Les valeurs correspondent aux niveaux précédents, attention à bien repérer la correspondance
corr.chemical=data.frame(levels(chemical1),sort(unique(chemical))) # juste pour conserver une trace de la correspondance entre niveaux de facteurs et leurs indices

place=as.numeric(d$place) # ici l'ordre des niveaux m'est indifférent
corr.place=data.frame(levels(d$place),sort(unique(place)))

# ne pas oublier de définir les index utilisés dans le modèle pour les boucles
ndata=length(time_z12)
nchemical=nlevels(chemical1)
nplace=nlevels(d$place)

# liste des variables à envoyer à jags, en les nommant comme elles sont nommées dans le script du modèle
xdata=list(time_z12=time_z12,chemical=chemical,place=place,Lnb_indiv=Lnb_indiv,ndata=ndata,nchemical=nchemical,nplace=nplace)

## coder les spécifications du modèles

# modèle
mod="baleines_modele.R"

# monitors
param=c("alpha","beta","gamma","delta")

# éventuellement, on peut avoir à spécifier des valeurs initiales (pas dans les exemples traités ici) 

## faire tourner le modèle sous Jags

# lancer le modèle et optimiser l'algorithme sur 1000 itérations
jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)

# burn-in 
update(jgm, n.iter=100)

# chaines de production
jgm2=coda.samples(jgm, variable.names=param, n.iter=200,thin=2)

# structure de l'objet
str(jgm2) # 3 chaines = liste de 3 éléments, chaque élément contient un data frame avec le contenu d'une chaine
dim(jgm2[[1]]) # nombre de lignes = nombre d'itérations conservées, nombre de colonnes = nombre de paramètres conservés (! ça monte vite)
head(jgm2[[1]]) 

# il peut être plus facile de séparer les chaines pour les manipuler
C1=as.data.frame(jgm2[[1]])
C2=as.data.frame(jgm2[[2]])
C3=as.data.frame(jgm2[[3]])

# on aura peut-être aussi besoin de travailler sur les 3 chaînes d'un coup
chains=rbind(C1,C2,C3)

## convergences des paramètres
plot(C1$alpha,type="l",ylab="valeur de alpha",xlab="itérations")
lines(C2$alpha,col="blue")
lines(C3$alpha,col="red") # rien qu'à regarder les chaînes brutes, on voit que le paramètre n'a pas convergé. On va donc refaire tourner le modèle sur plus d'itérations

# relancer le modèle

jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)
update(jgm, n.iter=2000)
jgm2=coda.samples(jgm, variable.names=param, n.iter=2000,thin=20)

# comme juste au dessus
C1=as.data.frame(jgm2[[1]])
C2=as.data.frame(jgm2[[2]])
C3=as.data.frame(jgm2[[3]])
chains=rbind(C1,C2,C3)

plot(C1$alpha,type="l",ylab="valeur de alpha",xlab="itérations")
lines(C2$alpha,col="blue")
lines(C3$alpha,col="red") # c'est beaucoup mieux! (il n'y avait peut-être pas besoin d'autant, mais comme ce modèle tourne vite on peut se permettre d'augmenter les chaînes)

# on regarde les autres paramètres d'intérêt
plot(C1$beta,type="l")
lines(C2$beta,col="blue")
lines(C3$beta,col="red")

plot(C1[,"delta[2]"],type="l") # bien regarder les paramètres liés aux variables facteur: ce sont souvent eux qui ont du mal à converger
lines(C2[,"delta[2]"],col="blue")
lines(C3[,"delta[2]"],col="red")

plot(C1[,"delta[3]"],type="l")
lines(C2[,"delta[3]"],col="blue")
lines(C3[,"delta[3]"],col="red")

plot(C1[,"gamma[2]"],type="l")
lines(C2[,"delta[2]"],col="blue")
lines(C3[,"delta[2]"],col="red")

# un test formel de diagnostic
gelman.plot(jgm2[,"alpha"]) # on cherche à ce que la valeur du Gelman soit en dessous de 1.1

gelman.plot(jgm2[,"beta"])
gelman.plot(jgm2[,"gamma[2]"])
gelman.plot(jgm2[,"gamma[3]"])
gelman.plot(jgm2[,"gamma[4]"])
gelman.plot(jgm2[,"delta[2]"])
gelman.plot(jgm2[,"delta[3]"])

# corrélation entre les paramètres du modèle
chart.Correlation(jgm2[,c("alpha","beta","gamma[2]","gamma[3]","gamma[4]","delta[2]","delta[3]")], histogram=TRUE, pch=19)

# on observe quelques corrélations assez fortes entre des paramètres. Il faudra chercher à les comprendre: elles créent des redondances dans le modèle
table(d[,c("place","chemical")])


# on voit par exemple qu'il n'y a pas de contrôle et de DMS en Antarctique. Peut-être que l'Antarctique est de trop dans ce modèle: il faudra essayer de le refitter sans
# l'objectif, à terme, sera d'éliminer toutes les corrélations. 

## mesurer l'ajustement

# on revient au script du modèle: prédire de nouvelles données
param=c("alpha","beta","gamma","delta","time_z12.new","residual","predicted")

# on relance le modèle en monitorant les données prédites
jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)
update(jgm, n.iter=2000)
jgm2=coda.samples(jgm, variable.names=param, n.iter=2000,thin=20)

# visualisation graphique
C1=as.data.frame(jgm2[[1]])
C2=as.data.frame(jgm2[[2]])
C3=as.data.frame(jgm2[[3]])
chains=rbind(C1,C2,C3)

# qualité de l'ajustement
Y.new=apply(chains[,236:348],2,"median")
dev.off()
plot(time_z12,Y.new,xlab="time_z12 réel",ylab="time_z12 prédit")
abline(0,1,col="red")
sum(Y.new>time_z12)/length(Y.new) # posterior probability check: une probabilité proche de 0.5 indique un bon ajustement

# et les résidus? 
predval=chains[,10:122] # valeurs prédites
predval1=apply(predval,2,"median")
resval=chains[,123:235] # résidus
resval1=apply(resval,2,"median")

par(mfrow=c(1,2))
plot(predval1,resval1,xlab="fitted",ylab="residuals")
x=loess(resval1~predval1)
y=cbind(predval1,predict(x))
y=y[order(y[,1]),]
lines(y,col="red")

# pour rappel, le plot des résidus du modèle fréquentiste
plot(modreg,1)


## interprétation du modèle: 
jgm.mcmc=as.mcmc.list(jgm2)

# estimateurs et intervalles de crédibilité
par(mfrow=c(1,2))
plot(C1$alpha,type="l",xlab="itérations",ylab="valeur de alpha")
lines(C2$alpha,col="blue")
lines(C3$alpha,col="red")
plot(density(chains$alpha),ylab="densité de probabilité",main="alpha")
# on peut ajouter le prior
lines(density(rnorm(1000,0,1000)),col="blue") # on a bien ajouté de l'info par rapport au prior!

# ce plot est aussi produit par la fonction plot.mcmc du package coda
plot(jgm2[,"alpha"]) 

# comment résumer l'info utilement?

param=c("alpha","beta","gamma","delta") # on ne garde que ce dont on a vraiment besoin et on relance le modèle
jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)
update(jgm, n.iter=2000)
jgm2=coda.samples(jgm, variable.names=param, n.iter=2000,thin=20)

# quelques outils simples pour résumer un modèle
summary(jgm2) # l'équivalent du summary d'un lm
summary(jgm2[,c("alpha","beta")]) # juste pour quelques paramètres


# et graphiquement?

# bien adapté aux effets catégoriques ou s'il y a beaucoup de coefficients à comparer entre eux
mcmc_intervals(jgm2, pars = c("alpha","gamma[2]","gamma[3]","gamma[4]"))
mcmc_intervals(jgm2, pars = c("alpha","delta[2]","delta[3]"))

# une autre représentation
mcmc_areas(jgm2,pars = c("alpha", "beta"),prob = 0.8,prob_outer = 0.99, point_est = "median")

# significativité des effets --> on oublie et on fait plus simple
sum(chains[,"alpha"]>0)/nrow(chains) # probabilité que l'intercept soit supérieur à 0
# probabilité que la pente beta soit supérieure à 0
sum(chains[,"beta"]>0)/nrow(chains) 
# probabilité que le temps de présence soit plus long en Islande qu'en Antarctique
sum(chains[,"delta[2]"]>0)/nrow(chains) 

# comparer les résultats au modèle fréquentiste
chem=factor(d$chemical,levels=c("CTL","CLAY","DMS","KRILL"))
plc=factor(d$place,levels=c("Antarctica","Iceland","Madagascar"))

modreg=lm(time_z12~log(nb_indiv)+chem+place,data=d)
regression=confint(modreg)

# plot des intervalles de confiance de la régression fréquentiste
plot(x=1:nrow(regression),coefficients(modreg),ylim=range(regression),xaxt="n",xlab="",ylab="coefficients and 95%CI",pch=21,bg="black")
axis(side=1,at=1:nrow(regression),labels=c("Intercept","log(nb indiv)","CLAY","DMS","KRILL","Iceland","Madagascar"))
arrows(x0=1:nrow(regression),x1=1:nrow(regression),y0=regression[,1],y1=regression[,2],code=3,length=0.02,angle=90)

# et des intervalles de crédibilité de la régression bayésienne
ICB=summary(jgm2)[[2]]
ICB=ICB[c("alpha","beta","gamma[2]","gamma[3]","gamma[4]","delta[2]","delta[3]"),]
points(1:nrow(regression)+0.2,ICB[,3],bg="grey",pch=21)
arrows(x0=1:nrow(regression)+0.2,x1=1:nrow(regression)+0.2,y0=ICB[,1],y1=ICB[,5],code=3,length=0.02,angle=90,col="grey")
abline(h=0,lty="dashed",col="gray30")
