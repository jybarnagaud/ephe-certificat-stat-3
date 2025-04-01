#-----------------------------------------------------------------#
### Formation "Certificat en analyse de données pour les écologues"
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# version 2023
# Cas pratique : modèles en cadre bayésien
#-----------------------------------------------------------------#

rm(list=ls())

library(ggplot2)
library(rjags) # interfaçage de JAGS avec R
library(coda) # utilitaires pour lire des chaînes MCMC
library(bayesplot) # d'autres utilitaires (pour les amateurs de ggplot)
library(MCMCvis) # encore d'autres utilitaires

#---------------------------------------#
#### modèle écologique sur le Pinson ####
#---------------------------------------#

## Préparation des données

# jeu de données: abondance de pinsons dans les forêts de Nouvelle Zélande: 3 réplicats par plot, variables année, jour de l'année, altitude, proportion de forêt native
pinson = read.csv2("donnees/nouvelle_zelande_pinson.csv")

# exploration des variables explicatives : altitude et proportion de forêt native
ggplot(pinson)+
  aes(x = altitude,y = PROPNATFOR)+
    geom_point()

# préparer la variable de réponse : comptage maxi sur 3*5 minutes
obs=as.matrix(pinson[,3:5])
obs.naif = apply(pinson[,3:5],1,"max")

# préparer les variables explicatives
PROPNATFOR=as.vector(scale(pinson$PROPNATFOR)) # on centre-réduit pour limiter les risques de mauvaises convergences (également applicable aux modèles fréquentistes)
altitude=as.vector(scale(pinson$altitude))
year=as.vector(as.numeric(factor(pinson$year))) # on a passé les années en facteur. 1=2005, 2=2006
julian_date=as.vector(scale(pinson$julian_date))

# les indices dont on aura besoin
npoints=nrow(obs)
nrepl=ncol(obs)
npoints.naif = length(obs.naif)
nyears=length(unique(year))

# liste des données
nzdata.naif=list(N=obs.naif,PROPNATFOR=PROPNATFOR,altitude=altitude,year=year,npoints=npoints,nyears=nyears)

# choisir les paramètres à suivre
monitor.naif=c("alpha","beta","delta","esp","r.naif","E","E.rep")

# lancer le modèle (run d'adaptation)
jnz.naif=jags.model(file="donnees/nouvelle_zelande_modele_ecol.R", data=nzdata.naif, n.chains=3, n.adapt=1000)

# burn in 
update(jnz.naif, n.iter=5000)

# codas (= chaines de Markov)
jnz.naif.coda=coda.samples(jnz.naif, variable.names=monitor.naif, n.iter=5000,thin=50)

# explorer les chaines
color_scheme_set("viridis")
mcmc_trace(jnz.naif.coda,pars=c("alpha[1]","alpha[2]","beta","delta"))

# explorer les paramètres de plusieurs manières
mcmc_areas(jnz.naif.coda,pars=c("alpha[1]","alpha[2]","beta","delta"))
mcmc_intervals(jnz.naif.coda,pars=c("alpha[1]","alpha[2]","beta","delta"))
mcmc_dens_overlay(jnz.naif.coda, pars = c("alpha[1]","alpha[2]"))

# explorer les corrélations entre paramètres (! on peut s'inquiéter de la relation
# entre altitude et forêt native!)
mcmc_pairs(jnz.naif.coda, pars=c("alpha[1]","alpha[2]","beta","delta"),
           off_diag_args = list(size = 1.5))

# explorer les résidus
resid.nz = MCMCchains(
  jnz.naif.coda,
  params = c("r.naif"))

med.resid.nz = apply(resid.nz,2,median)

esp.nz = MCMCchains(
  jnz.naif.coda,
  params = c("esp"))

med.esp.nz = apply(esp.nz,2,median)

# calcul des résidus de pearson
r.pearson = med.resid.nz/sqrt(med.esp.nz)
plot(med.esp.nz,r.pearson) # on retrouve le diagnostic classique de résidus

# ajustement
E.rep.nz = MCMCchains(
  jnz.naif.coda,
  params = c("E.rep"))

E.nz = MCMCchains(
  jnz.naif.coda,
  params = c("E"))

E.med = apply(E.nz,2,median)
E.rep.med = apply(E.rep.nz,2,median)
plot(E.med,E.rep.med) # cet ajustement n'est vraiment pas bon. 

#--------------------------------------#
#### modèle n-mixture sur le pinson ####
#--------------------------------------#

# liste des données
nzdata=list(obs=obs,PROPNATFOR=PROPNATFOR,altitude=altitude,year=year,julian_date=julian_date,npoints=npoints,nrepl=nrepl,nyears=nyears)

# choisir les paramètres à suivre
monitor=c("pmean","fit","fit.new","N","beta.det","beta","delta","mualpha","muyear","det.alpha")

# lancer le modèle
jnz=jags.model(file="nouvelle_zelande_modele.R", data=nzdata, n.chains=3, n.adapt=1000)

# ce message d'erreur suggère qu'il va falloir initialiser quelque chose - en l'occurrence, N
# on choisit d'initialiser N avec la valeur maximum observée par point
N=apply(obs,1,max) 

inis=list(N=N)

jnz=jags.model(file="donnees/nouvelle_zelande_modele.R", data=nzdata, inits=inis,n.chains=3, n.adapt=1000)

update(jnz, n.iter=5000) # ça va prendre un peu de temps!! bien réfléchir à son modèle... 
jnz2=coda.samples(jnz, variable.names=monitor, n.iter=5000,thin=50)

# mesure de convergence 
library(coda)
gelman.plot(jnz2)

# on va aussi regarder les chaines --> on pourrait faire un peu plus stable
color_scheme_set("viridis")
mcmc_trace(jnz2,pars=c("beta","delta","beta.det[1]",
                                "beta.det[2]","mualpha[1]",
                                "mualpha[2]","muyear","pmean"))

# corrélations entre les paramètres
mcmc_pairs(jnz2, pars=c("beta","delta","beta.det[1]",
                         "beta.det[2]","mualpha[1]",
                         "mualpha[2]","muyear","pmean"),
           off_diag_args = list(size = 1.5))

# mesure d'ajustement
fit.new.nz = MCMCchains(
  jnz2,
  params = c("fit.new"))

fit.nz = MCMCchains(
  jnz2,
  params = c("fit"))

plot(fit.nz,fit.new.nz)
abline(0,1,col="blue") # bon ajustement, qu'on peut quantifier: 
pval=sum(fit.nz>fit.new.nz)/length(fit.nz) # posterior predictive check pour quantifier l'ajustement - on veut être autour de 0.5
pval

# on veut quand même voir si les estimations sont réalistes
Npred.nz = MCMCchains(
  jnz2,
  params = "N")

med.N=apply(Npred.nz,2,median)
hist(med.N) # on a de 0 à 12 pinsons avec surtout des 0 - 6, réaliste

plot(obs.naif,med.N) # bonne relation entre le max sur 3 comptages et le N estimé

# quelle est la probabilité de détection
pdet.nz = MCMCchains(
  jnz2,
  params = "pmean")

hist(pdet.nz)  # probabilité de détection moyenne sur les points
quantile(pdet.nz,p=c(0.025,0.5,0.975)) # intervalle de crédibilité

# à chaque passage, on détecte un peu plus d'un pinson sur deux. 

# étude des résultats 
summary(jnz2)

# et graphiquement

# on commence par regarder les variables de la couche de process (effets de la forêt et de l'altitude)
mcmc_intervals(jnz2, pars = c("beta","delta"))

# ensuite, on va voir celles de la couche d'obs (intercept et interaction date*année)
mcmc_intervals(jnz2, pars = c("beta.det[1]","beta.det[2]","det.alpha"))

# on va maintenant vérifier les effets aléatoires année sur la couche de process
mcmc_intervals(jnz2, pars = c("mualpha[1]","mualpha[2]"))

#-------------------#
#### Tétras lyre ####
#-------------------#

# Dans cet autre exemple, on recode le modèle mixte sur les tétras lyre (cours du lundi)
tetras=read.table("comptages_tetras_lyre.txt",header=T,sep="\t")

# !!! Attention à l'ordre des lignes pour les colonnes liées aux années et aux UN
tetras=tetras[order(tetras$Annee),]

# matrice Je
Je=reshape(tetras[,c("Jeunes","Annee","UN")],v.names="Jeunes",idvar="UN",timevar="Annee",direction="wide")
rownames(Je)=Je[,1]
Je=Je[,-1]
colnames(Je)=substring(colnames(Je),first=8,nchar(colnames(Je)))
Je=as.matrix(Je)

# matrice Po
Po=reshape(tetras[,c("Poules","Annee","UN")],v.names="Poules",idvar="UN",timevar="Annee",direction="wide")
rownames(Po)=Po[,1]
Po=Po[,-1]
colnames(Po)=substring(colnames(Po),first=8,nchar(colnames(Po)))
Po=as.matrix(Po)
poul.moy=apply(log(Po),1,"mean",na.rm=T)
poul.var=apply(log(Po),1,"var",na.rm=T)
poul.var[is.na(poul.var)]=rep(var(as.numeric(Po),na.rm=T),sum(is.na(poul.var))) # pour les sites avec une seule donnée, on applique la variance de la population

# années
nan=ncol(Je)

# sites
nsites=nrow(Je)


# NAO
NAO0=unique(tetras[,c("Annee","NAOdjfm")])
NAO=as.numeric(NAO0[,2])

# données
xdata=list(Je=Je,Po=Po,NAO=NAO,nan=nan,nsites=nsites)

#-------------------#
### run du modèle ###
#-------------------#

# modèle
mod="tetras_mjags.R"

# paramètres à suivre
param=c("mualpha","mubeta","mubeta2","sigalpha","sigbeta","sigbeta2","Je.rep","Po.rep")

# adaptation
jags.tetras=jags.model(file=mod, data=xdata, n.chains=3, n.adapt=1000)

# burn in
update(jags.tetras, n.iter=1000)

# chaines de production
jags.tetras2=coda.samples(jags.tetras, variable.names=param, n.iter=2000,thin=20)

#------------------#
### convergences ###
#------------------#

# hyperparamètres (= effets aléatoires)
plot(jags.tetras2[,"mualpha"])
plot(jags.tetras2[,"mubeta"])
plot(jags.tetras2[,"mubeta2"])

# la convergence des hyperparamètres nous suffit, mais en toute rigueur il serait préférable de sonder quelques 
# paramètres locaux afin de vérifier qu'il n'y a pas de grosses aberrations

# on sonde quelques estimations de Po (uniquement sur des valeurs estimées, i.e. des NA dans la matrice originale)
po.na=which(is.na(Po),arr.ind=T)
xsamp=sample(po.na,10)
xdim=paste("Po[",po.na[xsamp,1],",",po.na[xsamp,2],"]",sep="")
plot(jags.tetras2[,xdim[8]])
# ça ne va pas : les chaines ne sont pas cohérente (probablement surfittage). On va tenter un modèle où on contraint le prior

#--------------------#
### nouveau modèle ###
#--------------------#

xdata=list(Je=Je,Po=Po,NAO=NAO,nan=nan,nsites=nsites,poul.moy=poul.moy,poul.var=poul.var)
param=c("mualpha","mubeta","mubeta2","sigalpha","sigbeta","sigbeta2","Je.rep","Po")

mod2="tetras_mjags2.R"

# adaptation
jags.tetras3=jags.model(file=mod2, data=xdata, n.chains=3, n.adapt=1000)

# burn in
update(jags.tetras3, n.iter=1000)

# chaines de production
jags.tetras4=coda.samples(jags.tetras3, variable.names=param, n.iter=2000,thin=20)

#-----------------------------------#
### on revérifie les convergences ###
#-----------------------------------#

# hyperparamètres (= effets aléatoires)
plot(jags.tetras4[,"mualpha"])
plot(jags.tetras4[,"mubeta"])
plot(jags.tetras4[,"mubeta2"])

# on sonde quelques estimations de Po (uniquement sur des valeurs estimées, i.e. des NA dans la matrice originale)
po.na=which(is.na(Po),arr.ind=T)
xsamp=sample(po.na,10)
xdim=paste("Po[",po.na[xsamp,1],",",po.na[xsamp,2],"]",sep="")
par(mfrow=c(5,2))
for(i in 1:10){
  plot(jags.tetras4[,xdim[i]])[1]
}
# c'est un peu mieux avec ce prior informatif, mais loin d'être parfait... 

#----------------#
### ajustement ###
#----------------#

# on monitore un nouveau jeu de paramètres (Po.rep au lieu de Po)
param2=c("mualpha","mubeta","mubeta2","sigalpha","sigbeta","sigbeta2","Je.rep","Po.rep")
jags.tetras5=coda.samples(jags.tetras3, variable.names=param2, n.iter=2000,thin=20)

# on passe les paramètres estimés en data.frame pour faciliter la manipulation
C1=as.data.frame(jags.tetras5[[1]])
C2=as.data.frame(jags.tetras5[[2]])
C3=as.data.frame(jags.tetras5[[3]])
chaines=rbind(C1,C2,C3)
je.ind=grep("Je",colnames(chaines))
je.rep=apply(chaines[,je.ind],2,"mean")
je.num=as.vector(Je)
je.rep[is.na(je.num)]=NA
plot(je.num,je.rep,xlab="Nombres de jeunes, valeurs observées",ylab="Nombres de jeunes, valeurs estimées")
abline(0,1,col="red")
sum(je.rep>je.num,na.rm=T)/length(je.rep)
# le modèle est légèrement biaisé vers des sous-estimations, mais ce n'est pas dramatique pour autant
# c'est en revanche un point à prendre en compte pour améliorer le modèle

po.ind=grep("Po",colnames(chaines))
po.rep=apply(chaines[,po.ind],2,"mean")
po.num=as.vector(Po)
po.rep[is.na(po.num)]=NA
plot(po.num,po.rep,xlab="Nombres de Poules, valeurs observées",ylab="Nombres de Poules, valeurs estimées")
abline(0,1,col="red")
sum(po.rep>po.num,na.rm=T)/length(po.rep)
# on surfitte clairement le nombre de poules. Ce n'est pas très grave pour notre modèle : on veut justement des 
# nombres de poules au plus près possible de ce qu'ils sont censés être, il n'y a aucune inférence dessus, donc 
# pour une fois c'est même préférable. Par contre, notre modèle ne serait pas utilisable pour prédire des nombres de
# poules ailleurs, sa transférabilité est donc faible. Si nous voulions un modèle plus transférable, nous essaierions de 
# formuler les variations d'abondances de poules de manière plus fine - par exemple avec des covariables

#---------------#
### inférence ###
#---------------#

# et graphiquement?
library(bayesplot)
pred1=mcmc_intervals(jags.tetras4, pars = c("mualpha","mubeta","mubeta2"))
pred1 # permet un aperçu graphique des paramètres mais pas terrible pour les effets continus

## on trace la courbe à la main

# courbe pour une itération
lambda.pred=chaines[1,"mualpha"]+chaines[1,"mubeta"]*NAO+chaines[1,"mubeta2"]*NAO*NAO
z=cbind(NAO,lambda.pred)
z=z[order(z[,1]),]
plot(z[,1],z[,2],type="l")

# courbes pour toutes les itérations (attention : on sort une courbe de nombres de jeunes par poules)
lambda.pred=data.frame()
for(i in 1:nrow(chaines)){
  lpred=exp(chaines[i,"mualpha"]+chaines[i,"mubeta"]*NAO+chaines[i,"mubeta2"]*NAO*NAO)
  lambda.pred=rbind(lambda.pred,lpred)
}

# pour chaque colonne (= chaque valeur de NAO), on prend le lambda prédit médian et son intervalle de crédibilité
icred=apply(lambda.pred,2,FUN="quantile",p=c(0.025,0.5,0.975))

# on trace la courbe médiane et son IC
z=cbind(NAO,t(icred))
z=z[order(z[,1]),]
plot(z[,1],z[,3],xlab="NAO",ylab="Nombre de jeunes par poule marginal",type="l",ylim=range(z[,-1]))
lines(z[,1],z[,2],lty="dashed")
lines(z[,1],z[,4],lty="dashed")
# on peut obtenir une courbe plus lisse en simulant des valeurs de NAO (uniquement dans la gamme de NAO de la variable observée)

# comparer avec la courbe du modèle fréquentiste correspondant
library(lme4)
ofpoules=log(tetras$Poules)
mix.tetras=glmer(Jeunes~NAOdjfm+I(NAOdjfm^2)+(1|UN),data=tetras,family="poisson",offset=ofpoules)
library(sjPlot)
plot_model(mix.tetras,type="pred",terms="NAOdjfm [all]")

