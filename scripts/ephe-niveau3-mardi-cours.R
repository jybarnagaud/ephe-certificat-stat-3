#-------------------------------------------------------------------------------#
### Formation "Certificat en analyse de donn�es pour les �cologues"
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# version 2023
# ce script permet la r�plication des exemples du cours sur les mod�les bay�siens
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
### donn�es baleines de B. Bouchard ###
### CEFE
#-------------------------------------#

d=read.table("baleines_Bouchard.txt",header=T,sep="\t")

#----------------------------#
### on explore les donn�es ###
#----------------------------#

par(mfrow=c(2,2))
boxplot(d$time_z12~d$chemical)
hist(d$time_z12)
hist(log(d$time_z12+1))
dev.off()

#-----------------------------------------------#
### coder une r�gression lin�aire en bay�sien ###
#-----------------------------------------------#

## le mod�le qu'on ferait en fr�quentiste
modreg=lm(time_z12~chemical+log(nb_indiv)+place,data=d)

# on a transform� nb_indiv pour r�duire l'�chelle des valeurs sur cette variable

# exploration rapide
par(mfrow=c(2,2))
plot(modreg)
summary(modreg)
p1=plot(ggemmeans(modreg,terms="chemical"),residuals=T)
p2=plot(ggemmeans(modreg,terms="nb_indiv"),residuals=T)
p3=plot(ggemmeans(modreg,terms="place"),residuals=T)
plot_grid(p1,p2,p3)

# le m�me mod�le en langage BUGS: script baleines_modele.R

## pr�parer les donn�es: revient � mettre les donn�es dans le format sp�cifi� dans le script du mod�le

time_z12=as.vector(d[,"time_z12"]) # on a sp�cifi� un vecteur pour la variable de r�ponse
Lnb_indiv=as.vector(log(d[,"nb_indiv"])) # un autre pour le nombre d'individus, qu'on a log-transform�

chemical1=factor(d$chemical,levels=c("CTL","CLAY","DMS","KRILL")) # on ordonne les niveaux des variables facteurs de mani�re judicieuse (ici les deux contr�les en premier)
chemical=as.numeric(chemical1) # on doit transformer les variables facteur en variables num�riques. Les valeurs correspondent aux niveaux pr�c�dents, attention � bien rep�rer la correspondance
corr.chemical=data.frame(levels(chemical1),sort(unique(chemical))) # juste pour conserver une trace de la correspondance entre niveaux de facteurs et leurs indices

place=as.numeric(d$place) # ici l'ordre des niveaux m'est indiff�rent
corr.place=data.frame(levels(d$place),sort(unique(place)))

# ne pas oublier de d�finir les index utilis�s dans le mod�le pour les boucles
ndata=length(time_z12)
nchemical=nlevels(chemical1)
nplace=nlevels(d$place)

# liste des variables � envoyer � jags, en les nommant comme elles sont nomm�es dans le script du mod�le
xdata=list(time_z12=time_z12,chemical=chemical,place=place,Lnb_indiv=Lnb_indiv,ndata=ndata,nchemical=nchemical,nplace=nplace)

## coder les sp�cifications du mod�les

# mod�le
mod="baleines_modele.R"

# monitors
param=c("alpha","beta","gamma","delta")

# �ventuellement, on peut avoir � sp�cifier des valeurs initiales (pas dans les exemples trait�s ici) 

## faire tourner le mod�le sous Jags

# lancer le mod�le et optimiser l'algorithme sur 1000 it�rations
jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)

# burn-in 
update(jgm, n.iter=100)

# chaines de production
jgm2=coda.samples(jgm, variable.names=param, n.iter=200,thin=2)

# structure de l'objet
str(jgm2) # 3 chaines = liste de 3 �l�ments, chaque �l�ment contient un data frame avec le contenu d'une chaine
dim(jgm2[[1]]) # nombre de lignes = nombre d'it�rations conserv�es, nombre de colonnes = nombre de param�tres conserv�s (! �a monte vite)
head(jgm2[[1]]) 

# il peut �tre plus facile de s�parer les chaines pour les manipuler
C1=as.data.frame(jgm2[[1]])
C2=as.data.frame(jgm2[[2]])
C3=as.data.frame(jgm2[[3]])

# on aura peut-�tre aussi besoin de travailler sur les 3 cha�nes d'un coup
chains=rbind(C1,C2,C3)

## convergences des param�tres
plot(C1$alpha,type="l",ylab="valeur de alpha",xlab="it�rations")
lines(C2$alpha,col="blue")
lines(C3$alpha,col="red") # rien qu'� regarder les cha�nes brutes, on voit que le param�tre n'a pas converg�. On va donc refaire tourner le mod�le sur plus d'it�rations

# relancer le mod�le

jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)
update(jgm, n.iter=2000)
jgm2=coda.samples(jgm, variable.names=param, n.iter=2000,thin=20)

# comme juste au dessus
C1=as.data.frame(jgm2[[1]])
C2=as.data.frame(jgm2[[2]])
C3=as.data.frame(jgm2[[3]])
chains=rbind(C1,C2,C3)

plot(C1$alpha,type="l",ylab="valeur de alpha",xlab="it�rations")
lines(C2$alpha,col="blue")
lines(C3$alpha,col="red") # c'est beaucoup mieux! (il n'y avait peut-�tre pas besoin d'autant, mais comme ce mod�le tourne vite on peut se permettre d'augmenter les cha�nes)

# on regarde les autres param�tres d'int�r�t
plot(C1$beta,type="l")
lines(C2$beta,col="blue")
lines(C3$beta,col="red")

plot(C1[,"delta[2]"],type="l") # bien regarder les param�tres li�s aux variables facteur: ce sont souvent eux qui ont du mal � converger
lines(C2[,"delta[2]"],col="blue")
lines(C3[,"delta[2]"],col="red")

plot(C1[,"delta[3]"],type="l")
lines(C2[,"delta[3]"],col="blue")
lines(C3[,"delta[3]"],col="red")

plot(C1[,"gamma[2]"],type="l")
lines(C2[,"delta[2]"],col="blue")
lines(C3[,"delta[2]"],col="red")

# un test formel de diagnostic
gelman.plot(jgm2[,"alpha"]) # on cherche � ce que la valeur du Gelman soit en dessous de 1.1

gelman.plot(jgm2[,"beta"])
gelman.plot(jgm2[,"gamma[2]"])
gelman.plot(jgm2[,"gamma[3]"])
gelman.plot(jgm2[,"gamma[4]"])
gelman.plot(jgm2[,"delta[2]"])
gelman.plot(jgm2[,"delta[3]"])

# corr�lation entre les param�tres du mod�le
chart.Correlation(jgm2[,c("alpha","beta","gamma[2]","gamma[3]","gamma[4]","delta[2]","delta[3]")], histogram=TRUE, pch=19)

# on observe quelques corr�lations assez fortes entre des param�tres. Il faudra chercher � les comprendre: elles cr�ent des redondances dans le mod�le
table(d[,c("place","chemical")])


# on voit par exemple qu'il n'y a pas de contr�le et de DMS en Antarctique. Peut-�tre que l'Antarctique est de trop dans ce mod�le: il faudra essayer de le refitter sans
# l'objectif, � terme, sera d'�liminer toutes les corr�lations. 

## mesurer l'ajustement

# on revient au script du mod�le: pr�dire de nouvelles donn�es
param=c("alpha","beta","gamma","delta","time_z12.new","residual","predicted")

# on relance le mod�le en monitorant les donn�es pr�dites
jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)
update(jgm, n.iter=2000)
jgm2=coda.samples(jgm, variable.names=param, n.iter=2000,thin=20)

# visualisation graphique
C1=as.data.frame(jgm2[[1]])
C2=as.data.frame(jgm2[[2]])
C3=as.data.frame(jgm2[[3]])
chains=rbind(C1,C2,C3)

# qualit� de l'ajustement
Y.new=apply(chains[,236:348],2,"median")
dev.off()
plot(time_z12,Y.new,xlab="time_z12 r�el",ylab="time_z12 pr�dit")
abline(0,1,col="red")
sum(Y.new>time_z12)/length(Y.new) # posterior probability check: une probabilit� proche de 0.5 indique un bon ajustement

# et les r�sidus? 
predval=chains[,10:122] # valeurs pr�dites
predval1=apply(predval,2,"median")
resval=chains[,123:235] # r�sidus
resval1=apply(resval,2,"median")

par(mfrow=c(1,2))
plot(predval1,resval1,xlab="fitted",ylab="residuals")
x=loess(resval1~predval1)
y=cbind(predval1,predict(x))
y=y[order(y[,1]),]
lines(y,col="red")

# pour rappel, le plot des r�sidus du mod�le fr�quentiste
plot(modreg,1)


## interpr�tation du mod�le: 
jgm.mcmc=as.mcmc.list(jgm2)

# estimateurs et intervalles de cr�dibilit�
par(mfrow=c(1,2))
plot(C1$alpha,type="l",xlab="it�rations",ylab="valeur de alpha")
lines(C2$alpha,col="blue")
lines(C3$alpha,col="red")
plot(density(chains$alpha),ylab="densit� de probabilit�",main="alpha")
# on peut ajouter le prior
lines(density(rnorm(1000,0,1000)),col="blue") # on a bien ajout� de l'info par rapport au prior!

# ce plot est aussi produit par la fonction plot.mcmc du package coda
plot(jgm2[,"alpha"]) 

# comment r�sumer l'info utilement?

param=c("alpha","beta","gamma","delta") # on ne garde que ce dont on a vraiment besoin et on relance le mod�le
jgm=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)
update(jgm, n.iter=2000)
jgm2=coda.samples(jgm, variable.names=param, n.iter=2000,thin=20)

# quelques outils simples pour r�sumer un mod�le
summary(jgm2) # l'�quivalent du summary d'un lm
summary(jgm2[,c("alpha","beta")]) # juste pour quelques param�tres


# et graphiquement?

# bien adapt� aux effets cat�goriques ou s'il y a beaucoup de coefficients � comparer entre eux
mcmc_intervals(jgm2, pars = c("alpha","gamma[2]","gamma[3]","gamma[4]"))
mcmc_intervals(jgm2, pars = c("alpha","delta[2]","delta[3]"))

# une autre repr�sentation
mcmc_areas(jgm2,pars = c("alpha", "beta"),prob = 0.8,prob_outer = 0.99, point_est = "median")

# significativit� des effets --> on oublie et on fait plus simple
sum(chains[,"alpha"]>0)/nrow(chains) # probabilit� que l'intercept soit sup�rieur � 0
# probabilit� que la pente beta soit sup�rieure � 0
sum(chains[,"beta"]>0)/nrow(chains) 
# probabilit� que le temps de pr�sence soit plus long en Islande qu'en Antarctique
sum(chains[,"delta[2]"]>0)/nrow(chains) 

# comparer les r�sultats au mod�le fr�quentiste
chem=factor(d$chemical,levels=c("CTL","CLAY","DMS","KRILL"))
plc=factor(d$place,levels=c("Antarctica","Iceland","Madagascar"))

modreg=lm(time_z12~log(nb_indiv)+chem+place,data=d)
regression=confint(modreg)

# plot des intervalles de confiance de la r�gression fr�quentiste
plot(x=1:nrow(regression),coefficients(modreg),ylim=range(regression),xaxt="n",xlab="",ylab="coefficients and 95%CI",pch=21,bg="black")
axis(side=1,at=1:nrow(regression),labels=c("Intercept","log(nb indiv)","CLAY","DMS","KRILL","Iceland","Madagascar"))
arrows(x0=1:nrow(regression),x1=1:nrow(regression),y0=regression[,1],y1=regression[,2],code=3,length=0.02,angle=90)

# et des intervalles de cr�dibilit� de la r�gression bay�sienne
ICB=summary(jgm2)[[2]]
ICB=ICB[c("alpha","beta","gamma[2]","gamma[3]","gamma[4]","delta[2]","delta[3]"),]
points(1:nrow(regression)+0.2,ICB[,3],bg="grey",pch=21)
arrows(x0=1:nrow(regression)+0.2,x1=1:nrow(regression)+0.2,y0=ICB[,1],y1=ICB[,5],code=3,length=0.02,angle=90,col="grey")
abline(h=0,lty="dashed",col="gray30")
